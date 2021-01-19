import pandas as pd
import numpy as np
import os
from tqdm import tqdm
import tempfile
import logging
from glob import glob
from pyarrow import ArrowIOError
from pyarrow.lib import ArrowInvalid
import time
import random
from pandas.api.types import is_numeric_dtype
from polyfun_utils import Logger, check_package_versions, set_snpid_index, configure_logger
from scipy.optimize import nnls

def splash_screen():
    print('*********************************************************************')
    print('* PolyPred (POLYgenic risk PREDiction)')
    print('* Version 1.0.0')
    print('* (C) 2020-2021 Omer Weissbrod')
    print('*********************************************************************')
    print()
    


def create_plink_range_file(df_betas, temp_dir, num_jk=200):
    #make sure that the df_betas is ordered
    is_chr_ordered = np.all(np.diff(df_betas['CHR']) >= 0)
    is_bp_ordered = True
    for chr_num, df_chr in df_betas.groupby('CHR'):        
        is_bp_ordered = is_bp_ordered & (np.all(np.diff(df_chr['BP']) >= 0))
        if not is_bp_ordered: break
    if not is_bp_ordered or not is_chr_ordered:
        df_betas.sort_values(['CHR', 'BP'], inplace=True)
        
    #bound num_jk by num_snps if we have very few SNPs
    num_jk = np.minimum(num_jk, df_betas.shape[0])
    
    #create df_ranges
    ranges_file = os.path.join(temp_dir, 'ranges.txt')
    ranges_list = [{'block':'block%d'%(block), 'lb':str(block), 'ub':str(block+1)} for block in range(1,num_jk+1)]
    df_ranges = pd.DataFrame(ranges_list)
    df_ranges = df_ranges[['block', 'lb', 'ub']]
    df_ranges.to_csv(ranges_file, sep='\t', header=False, index=False)
    
    #create df_scores
    scores_file = os.path.join(temp_dir, 'snp_scores.txt')
    separators = np.floor(np.linspace(0, df_betas.shape[0], num_jk+1)).astype(int)
    df_betas['score'] = 0
    is_in_range = np.zeros(df_betas.shape[0], dtype=np.bool)
    for i in range(len(separators)-1):
        is_in_range[separators[i] : separators[i+1]] = True
        df_betas.loc[is_in_range, 'score'] = i+1.5
        is_in_range[:] = False
    assert np.all(df_betas['score'] > 0)
        
    return ranges_file


def compute_prs_for_file(args,
                     plink_file,
                     df_betas,
                     temp_dir,
                     ranges_file=None,
                     keep_file=None
                     ):
                     
    #read the bim file
    plink_file_prefix = plink_file[:plink_file.rfind('.')]
    df_bim = pd.read_csv(plink_file_prefix+'.bim',
                           header=None,
                           names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'],
                           delim_whitespace=True)
    df_bim = set_snpid_index(df_bim)
    
    #keep only relevant SNPs
    df_betas = df_betas.merge(df_bim[['SNP']],
                              left_index=True,
                              right_index=True,
                              suffixes=('_betas', '_bim')
                              )
    if df_betas.shape[0] == 0:
        raise ValueError('No betas found for SNPs in plink file %s'%(plink_file_prefix))
        
    #create temp files
    betas_file = os.path.join(temp_dir, next(tempfile._get_candidate_names()))
    outfile_prs_temp = os.path.join(temp_dir, next(tempfile._get_candidate_names()))

    #save the betas to a file
    df_betas[['SNP_bim', 'A1', 'BETA']].to_csv(betas_file, header=False, index=False, sep='\t', float_format='%0.8e')
        
    #Run plink
    plink_exe = args.plink_exe if (args.plink_exe is not None) else args.plink2_exe
    plink_cmd = '%s --allow-no-sex --extract %s --out %s --memory %d --threads %d'%(
        plink_exe,
        betas_file,
        outfile_prs_temp,
        args.memory*1024,
        args.threads
    )
    if plink_file.endswith('.pgen'):
        plink_cmd += ' --bpfile %s --score %s cols=scoresums'%(plink_file_prefix, betas_file)
    elif plink_file.endswith('.bed'):
        plink_cmd += ' --bfile %s --score %s sum'%(plink_file_prefix, betas_file)
    else:
        raise ValueError('neither --bed nor --pgen specified')
    if ranges_file is not None:
        scores_file = os.path.join(temp_dir, next(tempfile._get_candidate_names()))
        df_betas[['SNP_bim', 'score']].drop_duplicates('SNP_bim').to_csv(scores_file, sep='\t', header=False, index=False)
        plink_cmd += ' --q-score-range %s %s'%(ranges_file, scores_file)
    if keep_file is not None:
        plink_cmd += ' --keep %s'%(keep_file)
    os.system(plink_cmd)
    
    #read results
    if ranges_file is None:
        if plink_file.endswith('.bed'):
            df_prs = pd.read_csv(outfile_prs_temp+'.profile', delim_whitespace=True)
        elif plink_file.endswith('.pgen'):
            df_prs = pd.read_csv(outfile_prs_temp+'.sscore', delim_whitespace=True)
            df_prs.rename(columns={'#IID':'IID', 'SCORE1_SUM':'SCORESUM'}, inplace=True)
            df_prs['FID'] = df_prs['IID']
        else:
            raise ValueError('neither --bed nor --pgen specified')
        df_prs.set_index('IID', inplace=True, drop=True)
        if np.any(df_prs.index.duplicated()):
            raise ValueError('duplicated iids found in %s'%(plink_file_prefix))
    else:
        df_prs = None
        if plink_file.endswith('.bed'): jk_files = glob(outfile_prs_temp+'.*.profile')
        elif plink_file.endswith('.pgen'): jk_files = glob(outfile_prs_temp+'.*.sscore')
        else: raise ValueError('neither --bed nor --pgen specified')
        for jk_file in jk_files:
            df_jk = pd.read_csv(jk_file, delim_whitespace=True)
            df_jk.rename(columns={'#IID':'IID', 'SCORE1_SUM':'SCORESUM'}, inplace=True)
            df_jk['FID'] = df_jk['IID']
            df_jk.set_index('IID', inplace=True, drop=True)
            if np.any(df_jk.index.duplicated()):
                raise ValueError('duplicated iids found in %s'%(plink_file_prefix))
            jk_file_basename = os.path.basename(jk_file)
            block_name = jk_file_basename.split('.')[-2]
            block_num = int(block_name[5:])
            scoresum_colname = 'SCORESUM.jk%d'%(block_num)
            df_jk.rename(columns={'SCORESUM':scoresum_colname}, inplace=True)
            if df_prs is None:
                df_prs = df_jk
                df_prs['SCORESUM'] = 0
            else:
                assert np.all(df_jk.index == df_prs.index)
                df_prs[scoresum_colname] = df_jk[scoresum_colname]
            df_prs['SCORESUM'] += df_jk[scoresum_colname]

    
    if df_prs is None:
        raise ValueError('The following plink command failed:\n%s'%(plink_cmd))
    return df_prs
    
    
def load_betas_files(betas_file, verbose=True):

    if verbose:
        logging.info('Loading betas file %s...'%(betas_file))
        t0 = time.time()
    try:
        df_betas = pd.read_parquet(betas_file)
        if len(df_betas.index.names) > 1:
            df_betas.reset_index(inplace=True)
    except (ArrowIOError, ArrowInvalid):
        if betas_file.endswith('.parquet'):
            raise IOError('corrupt parquet file: %s'%(betas_file))
        df_betas = pd.read_csv(betas_file, delim_whitespace=True)
    if verbose:
        logging.info('done in %0.2f seconds'%(time.time() - t0))

    #rename columns if needed
    df_betas.rename(columns={'sid':'SNP', 'nt1':'A1', 'nt2':'A2', 'BETA_MEAN':'BETA', 'ldpred_inf_beta':'BETA', 'chrom':'CHR', 'Chrom':'CHR', 'pos':'BP'}, inplace=True, errors='ignore')
    if not is_numeric_dtype(df_betas['CHR']):
        if df_betas['CHR'].str.startswith('chrom_').all():
            df_betas['CHR'] = df_betas['CHR'].str[6:].astype(np.int)
        else:
            raise ValueError('unknown CHR format')
    df_betas.rename(columns={'BETA_joint':'BETA', 'ALLELE1':'A1', 'ALLELE0':'A2', 'beta_mean':'BETA', 'MAF_BOLT':'A1Frq', 'Name':'SNP', 'A1Effect':'BETA', 'Name':'SNP', 'Chrom':'CHR', 'Position':'BP', 'beta':'BETA'}, inplace=True, errors='ignore')
    
    #create index
    df_betas = set_snpid_index(df_betas)
    
    #subset SNPs according to extract file
    if args.extract is not None:
        df_extract = pd.read_csv(args.extract, header=None, squeeze=True)
        df_betas = df_betas.loc[df_betas['SNP'].isin(df_extract)]
        if df_betas.shape[0]==0:
            raise ValueError('No SNPs remained after applying --extract')
        if verbose:
            logging.info('#SNPs after --extract: %s'%(df_betas.shape[0]))

    return df_betas

def computs_prs_all_files(args, betas_file, disable_jackknife=False, keep_file=None):

    #Check that all input files are valid
    found_pgen = False
    found_bed = False
    for file in args.files:
        if file.endswith('.bed'):
            if found_pgen: raise ValueError('Cannot combined .bed and .pgen files')
            if args.plink_exe is None: raise ValueError('You must provide --plink-exe if you use .bed files')
            found_bed = True
        elif file.endswith('.pgen'):
            if found_bed: raise ValueError('Cannot combined .bed and .pgen files')
            if args.plink2_exe is None: raise ValueError('You must provide --plink2-exe if you use .pgen files')
            found_pgen = True
        else:
            raise ValueError('File suffix must be either .bed or .pgen')

    #create a temp dir
    temp_dir = tempfile.mkdtemp()
    
    #load df_betas
    df_betas = load_betas_files(betas_file)

    #create range files for prs
    if args.num_jk == 0 or disable_jackknife:
        ranges_file = None
    else:
        ranges_file = create_plink_range_file(df_betas, temp_dir, num_jk=args.num_jk)

    #compute per-chromosome PRS
    df_prs_sum = None
    tqdm_files = tqdm(args.files)
    for plink_file in tqdm_files:
        if not os.path.exists(plink_file):
            raise IOError('%s doesn\'t exist'%(plink_file))
        tqdm_files.set_description(plink_file)
        df_prs_file = compute_prs_for_file(args, 
                     plink_file,
                     df_betas,
                     temp_dir,
                     ranges_file=ranges_file,
                     keep_file=keep_file
                     )
                     
        
        #add up the file-specific PRSs
        assert df_prs_file.shape[0]>0
        if df_prs_sum is None:
            df_prs_sum = df_prs_file
        else:
            assert np.all(df_prs_sum.index == df_prs_file.index)
            for c in df_prs_file.columns:
                if not c.startswith('SCORESUM'): continue
                if not c in df_prs_sum.columns: df_prs_sum[c] = 0
                df_prs_sum[c] += df_prs_file[c]
    
    #compute jackknife-block PRS
    if args.num_jk > 0 and not disable_jackknife:
        jk_columns = df_prs_sum.columns[df_prs_sum.columns.str.startswith('SCORESUM.jk')]
        assert np.allclose(df_prs_sum[jk_columns].sum(axis=1), df_prs_sum['SCORESUM'], atol=1e-3)
        for c in jk_columns:
            df_prs_sum[c] = df_prs_sum['SCORESUM'] - df_prs_sum[c]
    
    return df_prs_sum
    
    
def nonneg_lstsq(X, y):

    assert np.all(X.std(axis=0)>0)
    y_mean = y.mean()
    X_mean = X.mean(axis=0)
    y_c = y-y_mean
    X_c = X-X_mean
    coef, _ = nnls(X_c, y_c)
    y_hat = X.dot(coef)
    intercept = np.mean(y-y_hat)
    return coef, intercept
    
    
def estimate_mixing_weights(args):

    #read phenotypes
    df_pheno = pd.read_csv(args.pheno, names=['FID', 'IID', 'PHENO'], index_col='IID', delim_whitespace=True)
    
    #make sure that we didn't include a header line
    try:
        float(df_pheno['PHENO'].iloc[0])
    except:
        df_pheno = df_pheno.iloc[1:]
        df_pheno['PHENO'] = df_pheno['PHENO'].astype(np.float)
    if np.any(df_pheno.index.duplicated()):
        raise ValueError('duplicate ids found in %s'%(args.pheno))

    #compute a PRS for each beta file
    beta_files = args.betas.split(',')
    df_prs_sum_list = []
    for betas_file in beta_files:
        df_prs_sum = computs_prs_all_files(args, betas_file, disable_jackknife=True, keep_file=args.pheno)
        df_prs_sum_list.append(df_prs_sum[['SCORESUM']])
    for df_prs_sum in df_prs_sum_list:
        assert np.all(df_prs_sum.index == df_prs_sum_list[0].index)
    df_prs_sum_all = pd.concat(df_prs_sum_list, axis=1)
    
    #sync df_pheno and df_prs_sum_all
    index_shared = df_prs_sum_all.index.intersection(df_pheno.index)
    assert len(index_shared)>0
    if len(index_shared) < df_prs_sum_all.shape[0]:
        df_prs_sum_all = df_prs_sum_all.loc[index_shared]
    if df_pheno.shape[0] != df_prs_sum_all.shape[0] or np.any(df_prs_sum_all.index != df_pheno.index):
        df_pheno = df_pheno.loc[df_prs_sum_all.index]
    
    #compute mixing weights
    mix_weights, intercept = nonneg_lstsq(df_prs_sum_all.values, df_pheno['PHENO'].values)
    
    #create and print df_coef, and save it to disk
    df_coef = pd.Series(mix_weights, index=beta_files)
    df_coef.loc['intercept'] = intercept
    mix_weights_file = args.output_prefix+'.mixweights'
    df_coef.to_frame(name='mix_weight').to_csv(mix_weights_file, sep='\t')
    logging.info('Writing mixing weights to %s'%(mix_weights_file))
    
    #compute weighted betas
    df_betas_weighted = None
    for betas_file, mix_weight in zip(beta_files, mix_weights):
        df_betas = load_betas_files(betas_file)
        df_betas = df_betas[['SNP', 'CHR', 'BP', 'A1', 'A2', 'BETA']]
        df_betas['BETA'] *= mix_weight
        if df_betas_weighted is None:
            df_betas_weighted = df_betas
            continue
        
        index_shared = df_betas.index.intersection(df_betas_weighted.index)
        df_betas['BETA2'] = df_betas['BETA']
        df_new = df_betas_weighted.loc[index_shared].merge(df_betas.loc[index_shared, ['BETA2']], left_index=True, right_index=True)
        df_new['BETA'] += df_new['BETA2']
        del df_new['BETA2']
        del df_betas['BETA2']
        df_list = [df_new, df_betas.loc[~df_betas.index.isin(index_shared)], df_betas_weighted.loc[~df_betas_weighted.index.isin(index_shared)]]
        df_betas_weighted = pd.concat(df_list, axis=0)
    df_betas_weighted.sort_values(['CHR', 'BP', 'A1'], inplace=True)
    
    #save output to file
    df_betas_weighted[['SNP', 'CHR', 'BP', 'A1', 'A2', 'BETA']].to_csv(args.output_prefix+'.betas', sep='\t', index=False, float_format='%0.6e')
    logging.info('Saving weighted betas to %s'%(args.output_prefix+'.betas'))
    

def compute_prs(args):
    
    if args.betas.count(',') > 0:
        raise ValueError('--predict can only be used with a single betas file')
    df_prs_sum = computs_prs_all_files(args, args.betas, disable_jackknife=False, keep_file=args.keep)
    df_prs_sum.reset_index(inplace=True, drop=False)
    df_prs_sum.columns = df_prs_sum.columns.str.replace('SCORESUM', 'PRS')
    df_prs_sum_main = df_prs_sum[['FID', 'IID', 'PRS']]
    df_prs_sum_jk = df_prs_sum[['FID', 'IID'] + [c for c in df_prs_sum.columns if c.startswith('PRS.')]]
    
    df_prs_sum_main.to_csv(args.output_prefix+'.prs', sep='\t', index=False, float_format='%0.5f')
    if df_prs_sum_jk.shape[1]>1:
        df_prs_sum_jk.to_csv(args.output_prefix+'.prs_jk', sep='\t', index=False, float_format='%0.5f')
    logging.info('Saving PRS to %s'%(args.output_prefix+'.prs'))   


def check_args(args):
    if int(args.predict) + int(args.combine_betas) != 1:
        raise ValueError('you must specify either --predict or --combine-betas (but not both)')
    if args.plink_exe is None and args.plink2_exe is None:
        raise ValueError('you must specify either --plink-exe or --plink2-exe')
    if args.combine_betas:
        if args.keep is not None:
            raise ValueError('you cannot provide both --combine-betas and --keep')
        if args.pheno is None:
            raise ValueError('you must provide --pheno if you specify --combine-betas')
        if args.betas.count(',')==0:
            raise ValueError('you must provide multiple files in --betas if you specify --combine-betas')
    if args.num_jk<0:
        raise ValueError('--num-jk must be >=0')
    if args.pheno is not None and args.predict:
        raise ValueError('--pheno can only be used with --combine-betas')
        
    if len(list(args.files)) == 0:
        raise ValueError('no input files specified')
    
if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('--betas', required=True, help='files with SNP effect sizes (comma separated). A1 is the effect allele.')
    parser.add_argument('--output-prefix', required=True, help='Prefix of output file')
    
    parser.add_argument('--combine-betas', default=False, action='store_true', help='If specified, PolyPred will estimate mixing weights')
    parser.add_argument('--predict', default=False, action='store_true', help='If specified, PolyPred will compute PRS')
    parser.add_argument('--pheno', default=None, help='Phenotype file (required for estimating mixing weights)')
    parser.add_argument('--plink-exe', default=None, help='Plink 1.9 executable (required for .bed files)')
    parser.add_argument('--plink2-exe', default=None, help='Plink 2 executable (required for .pgen files)')
    
    parser.add_argument('--extract', default=None, help='A text file with rsids of SNPs to use (one per line)')
    parser.add_argument('--keep', default=None, help='A text file with ids of individuals to use (two columns per line, each containing FID,IID)')
    parser.add_argument('--num-jk', type=int, default=200, help='number of genomic jackknife blocks')
    
    parser.add_argument('--memory', type=int, default=2, help='Maximum memory usage (in GB)')
    parser.add_argument('--threads', type=int, default=1, help='Number of CPU threads')
    
    parser.add_argument('files', nargs='*', help='plink files (either bed files or .pgen files)')
    
    #check package versions
    check_package_versions()
    
    #show splash screen
    splash_screen()

    #extract args
    args = parser.parse_args()
    
    #check that the output directory exists
    if len(os.path.dirname(args.output_prefix))>0 and not os.path.exists(os.path.dirname(args.output_prefix)):
        raise ValueError('output directory %s doesn\'t exist'%(os.path.dirname(args.output_prefix)))
    
    #configure logger
    configure_logger(args.output_prefix)
    
    #check arguments
    check_args(args)
    
    #estimate mixing weights if needed
    if args.combine_betas:
        estimate_mixing_weights(args)

    #compute PRS if needed
    if args.predict:
        compute_prs(args)
    
    print()
    

    
