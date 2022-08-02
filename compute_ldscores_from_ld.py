import numpy as np
import pandas as pd
import os
import sys
import time
import logging
from tqdm import tqdm
import urllib.request
from pyarrow import ArrowIOError
from pyarrow.lib import ArrowInvalid
import tempfile
import scipy.sparse as sparse
from pandas.api.types import is_numeric_dtype
from polyfun_utils import configure_logger, set_snpid_index, TqdmUpTo, LONG_RANGE_LD_REGIONS
from ldstore.bcor import bcor


UKBB_LD_URL = 'https://data.broadinstitute.org/alkesgroup/UKBB_LD'
REGION_LENGTH = 3000000
UKB_N=337545
META_COLUMNS = ['SNP', 'CHR', 'BP', 'A1', 'A2']

    
def read_annot(annot_file):
    try:
        df_annot = pd.read_parquet(annot_file)
    except (ArrowIOError, ArrowInvalid):
        df_annot = pd.read_table(annot_file, sep='\s+')
    
    assert 'CHR' in df_annot.columns
    assert 'SNP' in df_annot.columns
    assert 'BP' in df_annot.columns
    assert 'A1' in df_annot.columns
    assert 'A2' in df_annot.columns
    
    for c in df_annot.columns:
        if c in META_COLUMNS: continue
        if not is_numeric_dtype(df_annot[c]):
            raise ValueError('Annotation %s does not have numeric values'%(c))
        
    
    df_annot = set_snpid_index(df_annot)
    
    return df_annot
    
    
def compute_R2_unbiased(R, n):
    np.power(R, 2, out=R)
    R *= ((n-1)/(n-2))
    R -= 1/(n-2)
    return R
            
    
    
def load_ld_npz(ld_dir, ld_prefix):
    
    #load the SNPs metadata
    gz_file = os.path.join(ld_dir, '%s.gz'%(ld_prefix))
    try:
        df_ld_snps = pd.read_table(gz_file, sep='\s+')
    except (ArrowIOError, ArrowInvalid):
        raise IOError('Corrupt file downloaded')
    df_ld_snps.rename(columns={'rsid':'SNP', 'chromosome':'CHR', 'position':'BP', 'allele1':'A1', 'allele2':'A2'}, inplace=True, errors='ignore')
    assert 'SNP' in df_ld_snps.columns
    assert 'CHR' in df_ld_snps.columns
    assert 'BP' in df_ld_snps.columns
    assert 'A1' in df_ld_snps.columns
    assert 'A2' in df_ld_snps.columns
    df_ld_snps = set_snpid_index(df_ld_snps)
        
    #load the LD matrix
    npz_file = os.path.join(ld_dir, '%s.npz'%(ld_prefix))
    logging.info('Loading LD from file %s'%(npz_file))
    t0 = time.time()
    try: 
        R = sparse.load_npz(npz_file).toarray()
        R += R.T
    except ValueError:
        raise IOError('Corrupt file downloaded')
    logging.info('Done in %0.2f seconds'%(time.time() - t0))
        
    #create df_R and return it
    df_R = pd.DataFrame(R, index=df_ld_snps.index, columns=df_ld_snps.index)
    return df_R, df_ld_snps
         

def get_bcor_meta(bcor_obj):
    df_ld_snps = bcor_obj.getMeta()
    df_ld_snps.rename(columns={'rsid':'SNP', 'position':'BP', 'chromosome':'CHR', 'allele1':'A1', 'allele2':'A2'}, inplace=True, errors='raise')
    df_ld_snps['CHR'] = df_ld_snps['CHR'].astype(np.int64)
    df_ld_snps['BP'] = df_ld_snps['BP'].astype(np.int64)
    df_ld_snps = set_snpid_index(df_ld_snps)
    return df_ld_snps
                
                
def load_ld_bcor(bcor_file):
    if not os.path.exists(bcor_file):
        raise IOError('%s not found'%(bcor_file))
    logging.info('Loading LD file %s'%(bcor_file))
    t0 = time.time()
    bcor_obj = bcor(bcor_file)
    df_ld_snps = get_bcor_meta(bcor_obj)
    ld_arr = bcor_obj.readCorr([])
    logging.info('Done in %0.2f seconds'%(time.time() - t0))
    df_R = pd.DataFrame(ld_arr, index=df_ld_snps.index, columns=df_ld_snps.index)
    return df_R, df_ld_snps
    
    
def load_ld(ld_file):
    if not ld_file.endswith('.bcor'):
        raise NotImplementedError('only .bcor files are currenty supported')
    df_R, df_ld_snps = load_ld_bcor(ld_file)
    chr_num = df_ld_snps['CHR'].iloc[0]
    region_start = df_ld_snps['BP'].min()
    region_end = df_ld_snps['BP'].max()
    
    return df_R, chr_num, region_start, region_end
        
    
def download_ukb_ld_file(chr_num, region_start, overwrite=False, ld_dir=None, no_cache=False):
    region_end = region_start + REGION_LENGTH
    ld_prefix = 'chr%d_%d_%d'%(chr_num, region_start, region_end)
    
    #if the files already exist, simply return them
    gz_file = os.path.join(ld_dir, '%s.gz'%(ld_prefix))
    npz_file = os.path.join(ld_dir, '%s.npz'%(ld_prefix))    
    if not overwrite and os.path.exists(gz_file) and os.path.exists(npz_file):
        try:
            df_R, _ = load_ld_npz(ld_dir, ld_prefix)
            return df_R
        except IOError:
            pass
    
    ### if we got here, we need to download the LD files ###
        
    #download the region files
    for suffix in ['npz', 'gz']:
    
        suffix_file = os.path.join(ld_dir, '%s.%s'%(ld_prefix, suffix))
        url = '%s/%s.%s'%(UKBB_LD_URL, ld_prefix, suffix)
        
        #special handling for long-range LD regions
        try:
            urllib.request.urlopen(url)
        except urllib.request.HTTPError:
            url += '2'
            urllib.request.urlopen(url)        
        
        with TqdmUpTo(unit='B', unit_scale=True, unit_divisor=1024, miniters=1, desc='downloading %s'%(url)) as t:
            urllib.request.urlretrieve(url, filename=suffix_file, reporthook=t.update_to)
            
            
    #load the LD matrix to memory
    df_R, _ = load_ld_npz(ld_dir, ld_prefix)
    
    #delete the downloaded files if requested
    if no_cache:
        for suffix in ['npz', 'gz']:
            suffix_file = os.path.join(ld_dir, '%s.%s'%(ld_prefix, suffix))
            os.remove(suffix_file)
    
    return df_R
    
    
    
    

def compute_ldscores_region(df_R_region, df_annot_region, n, is_binary, chr_num, region_start, region_end):
        
    #sync df_R_region and df_annot_region
    index_intersect = df_R_region.index.intersection(df_annot_region.index)
    if len(index_intersect)==0:
        raise ValueError('no SNPs in chromosome %d BP %d-%d had LD info'%(chr_num, region_start, region_end))
    if len(index_intersect) < df_R_region.shape[0]:
        logging.warning('Only %d/%d SNPs in chromosome %d BP %d-%d have annotations info. This may severely down-bias the LD-scores'%(len(index_intersect), df_R_region.shape[0], chr_num, region_start, region_end))
        is_keep = df_R_region.index.isin(index_intersect)
        df_R_region = df_R_region.loc[is_keep, is_keep]
    if df_annot_region.shape[0] != df_R_region.shape[0] or np.any(df_annot_region.index != df_R_region.index):
        df_annot_region = df_annot_region.loc[df_R_region.index]
    assert np.all(df_annot_region.index == df_R_region.index)
    
    #compute R2 (the unbiased estimator of R squared)
    R2_region = compute_R2_unbiased(df_R_region.values, n=n)        
    
    #compute LD scores
    df_annot_raw_region = df_annot_region.drop(columns=META_COLUMNS, errors='raise')
    annot_region = df_annot_raw_region.values
    if is_binary:
        annot_sparse = sparse.csc_matrix(annot_region)
        ld_scores_region = (annot_sparse.T.dot(R2_region)).T
    else:
        ld_scores_region = R2_region.dot(annot_region)
    df_ldscores_region = pd.DataFrame(ld_scores_region, index=df_annot_raw_region.index, columns=df_annot_raw_region.columns)
    df_ldscores_region = pd.concat((df_annot_region[META_COLUMNS], df_ldscores_region), axis=1)
    df_ldscores_region['distance_from_center'] = np.abs(df_ldscores_region['BP'] - ((region_start+region_end)//2))
    
    return df_ldscores_region


def compute_ldscores_chr(df_annot_chr, ld_dir=None, use_ukb=False, n=None, ld_files=None, no_cache=False):

    #create a temp output dir if required
    if use_ukb:
        if ld_dir is None: ld_dir = tempfile.mkdtemp()
        if not os.path.exists(ld_dir): raise IOError('LD directory %s doesn\'t exist'%(ld_dir))

    #infer chromosome number
    assert len(df_annot_chr['CHR'].unique()) == 1
    chr_num = df_annot_chr['CHR'].unique()[0]
    
    # #remove long-range LD regions
    # for r in LONG_RANGE_LD_REGIONS:
        # if r['chr'] != chr_num: continue
        # is_in_r = df_annot_chr['BP'].between(r['start'], r['end'])
        # if not np.any(is_in_r): continue
        # logging.warning('Removing %d SNPs from long-range LD region on chromosome %d BP %d-%d'%(is_in_r.sum(), r['chr'], r['start'], r['end']))
        # df_annot_chr = df_annot_chr.loc[~is_in_r]
    # if df_annot_chr.shape[0]==0:
        # raise ValueError('No SNPs found in chromosome %d (after removing long-range LD regions)'%(chr_num))
    
    #sort the SNPs by BP if needed
    if not np.all(np.diff(df_annot_chr['BP'])>=0):
        df_annot_chr = df_annot_chr.sort_values('BP', ascending=True)
    assert np.all(np.diff(df_annot_chr['BP'])>=0)
    
    #check if the data is binary
    df_annot_chr_raw = df_annot_chr.drop(columns=META_COLUMNS, errors='raise')
    if np.all(df_annot_chr_raw.dtypes == bool):
        is_binary = True
    elif np.all([len(np.unique(df_annot_chr_raw[c]))<=2 for c in df_annot_chr_raw.columns]):
        is_binary = True
    else:
        is_binary = False    
    
    #iterate over regions
    df_ldscores_regions_list = []
    
    #iterate over regions - UKB
    if use_ukb:
        for region_start in tqdm(range(1, df_annot_chr['BP'].max()+1, REGION_LENGTH)):
        
            #extract annotations for this region only
            region_end = region_start+REGION_LENGTH
            df_annot_region = df_annot_chr.query('%d <= BP <= %d'%(region_start, region_end))
            if df_annot_region.shape[0]==0: continue
            
            #skip over HLA region
            if chr_num==6 and region_start in [28000001, 29000001, 30000001]: continue
            
            #download the LD data
            df_R_region = download_ukb_ld_file(chr_num, region_start, ld_dir=ld_dir, no_cache=no_cache)
            
            df_ldscores_region = compute_ldscores_region(df_R_region, df_annot_region, n=UKB_N, is_binary=is_binary,
                                                         chr_num=chr_num, region_start=region_start, region_end=region_end)
            df_ldscores_regions_list.append(df_ldscores_region)

    #iterate over LD files
    else:
        for ld_file in ld_files:
            df_R_region, chr_num_ld, region_start, region_end = load_ld(ld_file)
            assert chr_num_ld == chr_num
            df_annot_region = df_annot_chr.query('%d <= BP <= %d'%(region_start, region_end))
            if df_annot_region.shape[0]==0: continue
            df_ldscores_region = compute_ldscores_region(df_R_region, df_annot_region, n=n, is_binary=is_binary,
                                                         chr_num=chr_num, region_start=region_start, region_end=region_end)
            df_ldscores_regions_list.append(df_ldscores_region)
        
    #keep the best ld-score for each SNP (the one closest to its region center)
    df_ldscores_chr = pd.concat(df_ldscores_regions_list, axis=0)
    df_ldscores_chr.sort_values('distance_from_center', ascending=True, inplace=True)
    df_ldscores_chr = df_ldscores_chr.loc[~df_ldscores_chr.index.duplicated(keep='first')]
    del df_ldscores_chr['distance_from_center']
    df_ldscores_chr.sort_values(['CHR', 'BP'], ascending=True, inplace=True)
    
    return df_ldscores_chr
    
    
    
def compute_ldscores_main(args, df_annot):
    
    #iterate over chromosomes
    df_ldscores_chr_list = []
    for chr_num, df_annot_chr in df_annot.groupby('CHR'):
        #df_ldscores_chr = compute_ldscores_chr(args, df_annot_chr)
        df_ldscores_chr = compute_ldscores_chr(df_annot_chr, ld_dir=args.ld_dir, use_ukb=args.ukb, n=args.n, ld_files=args.files, no_cache=args.no_cache)
        df_ldscores_chr_list.append(df_ldscores_chr)
        
    df_ldscores = pd.concat((df_ldscores_chr_list), axis=0)
    return df_ldscores
    
    
def main(args):

    #read annotations
    df_annot = read_annot(args.annot)
    
    #comptue LD-scores
    df_ldscores = compute_ldscores_main(args, df_annot)
    
    #save LD-scores to output file
    if args.gz_out:
        df_ldscores.to_csv(args.out, sep='\t', compression='gzip', index=False, float_format='%0.3f')
    else:
        df_ldscores.to_parquet(args.out, index=False)


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()    
    parser.add_argument('--annot', required=True, help='annotations file')
    parser.add_argument('--out', required=True, help='output file')
    parser.add_argument('--ukb', default=False, action='store_true', help='if specified, the script will download and use UK Biobank LD files')
    parser.add_argument('--n', type=int, help='Sample size used to compute ')
    parser.add_argument('--gz-out', default=False, action='store_true', help='if specified, the output file will be a gzipped text file instead of a parquet file')
    parser.add_argument('--ld-dir', default=None, help='the path of an LD files directory. If not provided, LD files will be downloaded to a temporary directory')
    parser.add_argument('--no-cache', default=False, action='store_true', help='If this flag is specified, the LD files will be removed from the ld-dir after downloading them to save disk space')
    parser.add_argument('files', nargs='*', help='bcor files')
    args = parser.parse_args()
    
    configure_logger(args.out)
    
    #check input arguments
    if args.ld_dir is not None and not os.path.exists(args.ld_dir):
        raise ValueError('Specified LD directory %s doesn\'t exist'%(args.ld_dir))
    if len(os.path.dirname(args.out))>0 and not os.path.exists(os.path.dirname(args.out)):
        raise ValueError('output directory %s doesn\'t exist'%(os.path.dirname(args.out)))
    if args.ukb:
        if args.n is not None:
            raise ValueError('--n cannot be specified together with --ukb')
        if len(args.files) > 0:
            raise ValueError('input LD files cannot be provided when using --ukb')
    else:
        if args.n is None:
            raise ValueError('--n must be specified when not using --ukb')
        if args.ld_dir is not None:
            raise ValueError('--ld-dir cannot be specified when not using --ukb')
        if len(args.files)==0:
            raise ValueError('no input LD files specified. Did you mean to provide the flag --ukb?')
    
    main(args)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

