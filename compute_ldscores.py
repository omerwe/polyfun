import numpy as np
import pandas as pd
import os
import sys
import time
from ldsc_polyfun import ldscore, parse
import logging
from pandas.api.types import is_numeric_dtype
from polyfun_utils import configure_logger, set_snpid_index, SNP_COLUMNS
from pyarrow import ArrowIOError



def __filter__(fname, noun, verb, merge_obj):
    merged_list = None
    if fname:
        f = lambda x,n: x.format(noun=noun, verb=verb, fname=fname, num=n)
        x = parse.FilterFile(fname)
        c = 'Read list of {num} {noun} to {verb} from {fname}'
        print(f(c, len(x.IDList)))
        merged_list = merge_obj.loj(x.IDList)
        len_merged_list = len(merged_list)
        if len_merged_list > 0:
            c = 'After merging, {num} {noun} remain'
            print(f(c, len_merged_list))
        else:
            error_msg = 'No {noun} retained for analysis'
            raise ValueError(f(error_msg, 0))

        return merged_list


        
def compute_ldscores(args):

    #read bim/snp
    array_snps = parse.PlinkBIMFile(args.bfile+'.bim')
    df_bim = array_snps.df
    if len(df_bim['CHR'].unique()) > 1:
        raise ValueError('plink file includes multiple chromosomes. Please specify a plink file with a single chromosome')
    df_bim = set_snpid_index(df_bim)
    
    #read annotations
    keep_snps = None
    if args.annot is not None:
    
        try:
            df_annot = pd.read_parquet(args.annot)
        except ArrowIOError:
            df_annot = pd.read_table(args.annot, delim_whitespace=True)
        
        #Remove annotations of SNPs that are not in the .bim file
        df_annot = set_snpid_index(df_annot)
        df_annot = df_annot.loc[df_annot.index.isin(df_bim.index)]

        #make sure that all SNPs have annotations
        if np.any(~df_bim.index.isin(df_annot.index)):
            error_msg = 'Not all SNPs have annotation values'
            if args.allow_missing:
                is_good_snp = df_bim.index.isin(df_annot.index)
                if not np.any(is_good_snp):
                    raise ValueError('No SNPs have annotations')
                keep_snps = np.where(is_good_snp)[0]
                logging.warning(error_msg)
                logging.warning('Keeping only %d/%d SNPs that have annotations'%(is_good_snp.sum(), len(is_good_snp)))
            else:
                raise ValueError(error_msg + '. If you wish to omit the missing SNPs, please use the flag --allow-missing')

        #make sure that all of the annotations are numeric
        for c in df_annot.columns:
            if c in SNP_COLUMNS: continue
            if not is_numeric_dtype(df_annot[c]):
                raise ValueError('Annotation %s does not have numeric values'%(c))

    #find #individuals in bfile
    fam_file = args.bfile+'.fam'
    df_fam = pd.read_table(fam_file, header=None, usecols=[5], delim_whitespace=True)
    n = df_fam.shape[0]

    #find keep_indivs    
    if args.keep is None:
        keep_indivs= None
    else:
        array_indivs = parse.PlinkFAMFile(args.bfile+'.fam')
        keep_indivs = __filter__(args.keep, 'individuals', 'include', array_indivs)
        logging.info('after applying --keep, %d individuals remain'%(len(keep_indivs)))
    
    #read plink file    
    bed_file = args.bfile+'.bed'
    geno_array = ldscore.PlinkBEDFile(bed_file, n, array_snps, keep_snps=keep_snps,
        keep_indivs=keep_indivs, mafMin=None)
    
    #remove omitted SNPs from df_bim
    if len(geno_array.kept_snps) < df_bim.shape[0]:
        assert np.all(np.array(geno_array.kept_snps) == np.sort(np.array(geno_array.kept_snps)))
        assert geno_array.kept_snps[-1] < df_bim.shape[0]
        df_bim = df_bim.iloc[geno_array.kept_snps]
        
    #rearrange annotations to match the order of SNPs in the plink file
    if args.annot is not None:
        assert df_annot.shape[0] >= df_bim.shape[0]
        if (df_annot.shape[0] > df_bim.shape[0]) or np.any(df_annot.index != df_bim.index):
            assert np.all(df_bim.index.isin(df_annot.index))
            df_annot = df_annot.loc[df_bim.index]
                    
    # determine block widths
    num_wind_args = np.array((args.ld_wind_snps, args.ld_wind_kb, args.ld_wind_cm), dtype=bool)
    if np.sum(num_wind_args) != 1:
        raise ValueError('Must specify exactly one --ld-wind option')
    if args.ld_wind_snps:
        max_dist = args.ld_wind_snps
        coords = np.array(list(range(geno_array.m)))
    elif args.ld_wind_kb:
        max_dist = args.ld_wind_kb*1000
        coords = np.array(df_bim['BP'])
        if len(np.unique(coords)) == 1:
            raise ValueError('bim file has no basepair data --- please use a different ld-wind option')
    elif args.ld_wind_cm:
        max_dist = args.ld_wind_cm
        coords = np.array(df_bim['CM'])
        if len(np.unique(coords)) == 1:
            raise ValueError('bim file has no CM data --- please use a different ld-wind option')
        
    #compute LD-scores
    block_left = ldscore.getBlockLefts(coords, max_dist)
    if block_left[len(block_left)-1] == 0:
        error_msg = 'Only a single block selected - this is probably a mistake'
        raise ValueError(error_msg)
    t0 = time.time()
    geno_array._currentSNP = 0
    annot_values = (None if args.annot is None else df_annot.drop(columns=SNP_COLUMNS).values)
    ldscores = geno_array.ldScoreVarBlocks(block_left, args.chunk_size, annot=annot_values)
    
    #create an ldscores df
    if args.annot is None:
        df_ldscores = pd.DataFrame(ldscores, columns=['base'])
    else:
        df_ldscores = pd.DataFrame(ldscores, columns=df_annot.drop(columns=SNP_COLUMNS).columns)
        
    #add SNP identifier columns
    for c in SNP_COLUMNS:
        df_ldscores[c] = df_bim[c].values
    df_ldscores = df_ldscores[SNP_COLUMNS + list(df_ldscores.drop(columns=SNP_COLUMNS).columns)]
    
    return df_ldscores
        


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()    
    parser.add_argument('--bfile', required=True, help='Plink file')
    parser.add_argument('--out', required=True, help='output file')
    parser.add_argument('--annot', default=None, help='annotations file (if not specified, compute the global LD-score)')
    parser.add_argument('--ld-wind-cm', type=float, default=None, help='window size to be used for estimating LD-scores in units of centiMorgans (cM).')
    parser.add_argument('--ld-wind-kb', type=int, default=None, help='window size to be used for estimating LD-scores in units of Kb.')
    parser.add_argument('--ld-wind-snps', type=int, default=None, help='window size to be used for estimating LD-scores in units of SNPs.')
    parser.add_argument('--chunk-size',  type=int, default=50, help='chunk size for LD-scores calculation')
    #parser.add_argument('--extract',  default=None, help='File with rsids of SNP to use')
    parser.add_argument('--keep',  default=None, help='File with ids of individuals to use')
    parser.add_argument('--allow-missing', default=False, action='store_true', help='If specified, The script will not terminate if some SNPs with sumstats are not found in the annotations files')
    args = parser.parse_args()
    
    configure_logger(args.out)
    
    if args.ld_wind_cm is None and args.ld_wind_kb is None and args.ld_wind_snps is None:
        args.ld_wind_cm = 1.0
        logging.warning('no ld-wind argument specified.  We will use --ld-cm 1.0')
        
    df_ldscores = compute_ldscores(args)
    df_ldscores.to_parquet(args.out)
    
    
    
    
    
    