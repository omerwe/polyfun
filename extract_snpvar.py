import pandas as pd
import numpy as np
import os
import logging
from polyfun import configure_logger


def check_package_versions():
    from pkg_resources import parse_version
    if parse_version(pd.__version__) < parse_version('0.24.0'):
        raise ValueError('your pandas version is too old --- please update pandas')
    

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--snps', required=True, help='Whitespace-delimited file with SNPs to extract. Must include columns A1,A2 and either (1) SNP or (2) both CHR and BP')
    parser.add_argument('--out', required=True, help='Prefix of the name of the output file')
    parser.add_argument('--allow-missing', default=False, action='store_true', help='If specified, the script will not terminate if some SNPs are not found in the meta file')
    parser.add_argument('--q', type=int, default=100, help='The maximum ratio between the largest and smallest prior causal probabilities')
    args = parser.parse_args()
    
    #check package versions
    check_package_versions()
    
    #configure logger
    configure_logger(args.out)
    
    #read snps file
    df_snps = pd.read_table(args.snps, delim_whitespace=True)    
    if 'A1' not in df_snps.columns:
        raise ValueError('missing column A1')
    if 'A2' not in df_snps.columns:
        raise ValueError('missing column A2')
    if 'SNP' not in df_snps.columns:
        if 'CHR' not in df_snps.columns:
            raise ValueError('You must provide either a SNP or a CHR column')
        if 'BP' not in df_snps.columns:
            raise ValueError('You must provide either a SNP or a BP column')
            
    #read df_meta        
    script_dir = os.path.dirname(os.path.realpath(__file__))
    df_meta1 = pd.read_parquet(os.path.join(script_dir, 'snpvar_meta.chr1_7.parquet'))
    df_meta2 = pd.read_parquet(os.path.join(script_dir, 'snpvar_meta.chr8_22.parquet'))
    df_meta = pd.concat([df_meta1, df_meta2], axis=0)
    del df_meta1
    del df_meta2
    
    #truncate the ratio between the largest and smallest per-SNP h^2    
    min_snpvar = df_meta['snpvar_bin'].max() / args.q
    snpvar_sum = df_meta['snpvar_bin'].sum()
    df_meta.loc[df_meta['snpvar_bin'] < min_snpvar, 'snpvar_bin'] = min_snpvar
    snpvar_sum2 = df_meta['snpvar_bin'].sum()
    df_meta['snpvar_bin'] *= snpvar_sum / snpvar_sum2
    snpvar_sum3 = df_meta['snpvar_bin'].sum()
    assert np.isclose(snpvar_sum3, snpvar_sum)    
    
    #subset df_meta to include a small superset of df_snps
    if 'SNP' in df_snps.columns:
        df_meta = df_meta.loc[df_meta['SNP'].isin(df_snps['SNP'])]
    else:
        df_meta = df_meta.loc[df_meta['BP'].isin(df_snps['BP'])]
    
    #duplicate df_meta to include every SNP twice, with alternating A1/A2 alleles (to handle allele flips)
    df_meta2 = df_meta.copy()
    df_meta2['A1'] = df_meta['A2']
    df_meta2['A2'] = df_meta['A1']
    df_meta = pd.concat([df_meta, df_meta2], axis=0)    
    
    #merge the dfs
    if 'SNP' in df_snps.columns:
        df = df_meta.merge(df_snps, on=['SNP', 'A1', 'A2'], how='inner')
    else:
        df = df_meta.merge(df_snps, on=['CHR', 'BP', 'A1', 'A2'], how='inner')
        
    #If we didn't find everything, write a list of missing SNPs to an output file
    if df.shape[0] < df_snps.shape[0]:
        if 'SNP' in df_snps.columns:
            df_snps.index = df_snps['SNP'] + '.' \
                          + df_snps['A1'] + '.' \
                          + df_snps['A2']
            df_meta.index = df_meta['SNP'] + '.' \
                          + df_meta['A1'] + '.' \
                          + df_meta['A2']
        else:
            df_snps.index = df_snps['CHR'].astype(str) + '.' \
                          + df_snps['BP'].astype(str) + '.' \
                          + df_snps['A1'] + '.' \
                          + df_snps['A2']
            df_meta.index = df_meta['CHR'].astype(str) + '.' \
                          + df_meta['BP'].astype(str) + '.' \
                          + df_meta['A1'] + '.' \
                          + df_meta['A2']
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
        
    #write output to file
    logging.info('Writing output file to %s'%(args.out))
    df.to_csv(args.out, sep='\t', index=False, float_format='%0.4e')
    
        
    
    
    
