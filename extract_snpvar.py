import pandas as pd
import numpy as np
import os
import logging
import time
from pyarrow import ArrowIOError
from pyarrow.lib import ArrowInvalid
from polyfun_utils import check_package_versions, configure_logger, set_snpid_index, SNP_COLUMNS



if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--sumstats', required=True, help='Whitespace-delimited file with SNPs to extract. Must include columns CHR,BP,A1,A2')
    parser.add_argument('--out', required=True, help='Prefix of the name of the output file')
    parser.add_argument('--allow-missing', default=False, action='store_true', help='If specified, the script will not terminate if some SNPs are not found in the meta file')
    parser.add_argument('--q', type=int, default=100, help='The maximum ratio between the largest and smallest prior causal probabilities')
    args = parser.parse_args()
    
    #check package versions
    check_package_versions()
    
    #configure logger
    configure_logger(args.out)
    
    #read sumtats file
    logging.info('Loading sumstats files...')
    t0 = time.time()
    try:
        df_snps = pd.read_parquet(args.sumstats)
    except (ArrowIOError, ArrowInvalid):
        df_snps = pd.read_table(args.sumstats, sep='\s+')
    if 'A1' not in df_snps.columns:
        raise ValueError('missing column A1')
    if 'A2' not in df_snps.columns:
        raise ValueError('missing column A2')
    if 'CHR' not in df_snps.columns:
        raise ValueError('missing column CHR')
    if 'BP' not in df_snps.columns:
        raise ValueError('missing column BP')
            
    #set index
    df_snps = set_snpid_index(df_snps)
    logging.info('Done in %0.2f seconds'%(time.time() - t0))
    
    #make sure there aren't any duplicated SNPs
    if np.any(df_snps.index.duplicated()):
        raise ValueError('duplicate SNPs found in output - please make sure there aren\'t any duplicate SNPs in your sumstats file')
    
            
    #read df_meta
    logging.info('Loading meta-analyzed per-SNP-h2 files...')
    t0 = time.time()
    script_dir = os.path.dirname(os.path.realpath(__file__))
    df_meta1 = pd.read_parquet(os.path.join(script_dir, 'snpvar_meta.chr1_7.parquet'))
    df_meta2 = pd.read_parquet(os.path.join(script_dir, 'snpvar_meta.chr8_22.parquet'))
    df_meta = pd.concat([df_meta1, df_meta2], axis=0)
    del df_meta1
    del df_meta2
    df_meta = set_snpid_index(df_meta)
    logging.info('Done in %0.2f seconds'%(time.time() - t0))
    
    #truncate the ratio between the largest and smallest per-SNP h^2    
    min_snpvar = df_meta['snpvar_bin'].max() / args.q
    snpvar_sum = df_meta['snpvar_bin'].sum()
    df_meta.loc[df_meta['snpvar_bin'] < min_snpvar, 'snpvar_bin'] = min_snpvar
    snpvar_sum2 = df_meta['snpvar_bin'].sum()
    df_meta['snpvar_bin'] *= snpvar_sum / snpvar_sum2
    snpvar_sum3 = df_meta['snpvar_bin'].sum()
    assert np.isclose(snpvar_sum3, snpvar_sum)    
    
    #Make sure there aren't any shared columns (except for SNP columns)
    for c in df_snps.columns:
        if c not in SNP_COLUMNS and c in df_meta.columns:
            raise ValueError('Column %s appears in both the sumstats files and in the meta-analysis files'%(c))
    
    #merge the dfs
    logging.info('Merging sumstats with per-SNP h2 data...')
    t0 = time.time()
    df_meta = df_meta.loc[df_meta.index.isin(df_snps.index)]
    df_snps = df_snps.rename(columns = {'A1': 'A_eff'})
    df = df_meta.merge(df_snps.drop(columns=SNP_COLUMNS, errors='ignore'), left_index=True, right_index=True)
    
    
    #flip Z-sign if A1 of prior not match A1 of sumstats
    if 'Z' in df.columns:
        is_flipped = df['A2'] == df['A_eff']
        if is_flipped.sum() > 0:
            df.loc[is_flipped, 'Z'] *= -1
            logging.info('Flipping the Z-sign of %d SNPs that A1 in sumstats = A2 in the per-SNP h2 data'%(is_flipped.sum()))
            
    #drop unneeded columns
    df = df.drop(columns = 'A_eff')
        
    #report processing time
    logging.info('Done in %0.2f seconds'%(time.time() - t0))
        
    #If we didn't find everything, write a list of missing SNPs to an output file
    if df.shape[0] < df_snps.shape[0]:
        df_miss = df_snps.loc[~df_snps.index.isin(df_meta.index)]
        df_miss.to_csv(args.out+'.miss.gz', sep='\t', index=False, compression='gzip')
        error_message = 'Not all SNPs in the SNPs file were found in the meta file. Wrote a list of missing SNPs to %s'%(args.out+'.miss.gz')
        if args.allow_missing:
            logging.warning(error_message)
        else:
            raise ValueError(error_message)
        
    #normalize the prior-causal probabilities
    df['SNPVAR'] = df['snpvar_bin'] #/ df['snpvar_bin'].sum()
    #assert df['SNPVAR'].sum() == 1
    del df['snpvar_bin']
    
    #make sure there aren't any duplicated SNPs
    if np.any(df.index.duplicated()):
        raise ValueError('duplicate SNPs found in output - please make sure there aren\'t any duplicate SNPs in your sumstats file')
        
    #write output to file
    logging.info('Writing output file to %s'%(args.out))
    if args.out.endswith('.parquet'):
        df.to_parquet(args.out, index=False)
    else:
        df.to_csv(args.out, sep='\t', index=False, float_format='%0.4e')
    
        
    
    
    
