import pandas as pd
import numpy as np
import os
import logging
import time
from pyarrow import ArrowIOError
from pyarrow.lib import ArrowInvalid
from polyfun_utils import check_package_versions, configure_logger, set_snpid_index



if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--pips', required=True, help='Whitespace-delimited file with SNPs ids and PPs. Must include columns CHR,BP,A1,A2,PIP')
    parser.add_argument('--annot', required=True, help='Annotations file')
    parser.add_argument('--out', required=True, help='Prefix of the name of the output file')
    parser.add_argument('--allow-missing', default=False, action='store_true', help='If specified, the script will not terminate if some SNPs are not found in the annotation file')
    parser.add_argument('--pip-cutoff', type=float, default=0.95, help='Only extract annotations for SNPs with PIP greater than this cutoff (default is 0.95)')
    args = parser.parse_args()
    
    #check package versions
    check_package_versions()
    
    #configure logger
    configure_logger(args.out)
    
    if not 0<=args.pip_cutoff<=1:
        raise ValueError('--pip-cutoff must be between 0 and 1')
    
    #read pips file
    try:
        df_snps = pd.read_parquet(args.pips)
    except (ArrowIOError, ArrowInvalid):
        df_snps = pd.read_table(args.pips, delim_whitespace=True)
    if 'A1' not in df_snps.columns:
        raise ValueError('missing column A1')
    if 'A2' not in df_snps.columns:
        raise ValueError('missing column A2')
    if 'CHR' not in df_snps.columns:
        raise ValueError('missing column CHR')
    if 'BP' not in df_snps.columns:
        raise ValueError('missing column BP')
    if 'PIP' not in df_snps.columns:
        raise ValueError('missing column PIP')
            
    #set index
    df_snps = set_snpid_index(df_snps)
    
    #restrict to SNPs with a large PIP
    df_snps = df_snps.query('PIP>=%s'%(args.pip_cutoff))
    if df_snps.shape[0]==0:
        raise ValueError('No SNPs with PIP>=%s found'%(args.pip_cutoff))
                
    #read df_annot
    logging.info('Loading annotations file...')
    t0 = time.time()
    try:
        df_annot = pd.read_parquet(args.annot)
    except (ArrowIOError, ArrowInvalid):
        df_annot = pd.read_table(args.annot, delim_whitespace=True)
    df_annot = set_snpid_index(df_annot)
    logging.info('Done in %0.2f seconds'%(time.time() - t0))
    
    #make sure that all SNPs have annotations
    is_annot_missing = ~(df_snps.index.isin(df_annot.index))
    if np.any(is_annot_missing):
        err_msg = 'The annotations of %d SNPs were not found in the annotation files'%(is_annot_missing.sum())
        if args.allow_missing:
            logging.warning(err_msg)
            df_snps = df_snps.loc[~is_annot_missing]
        else:
            raise ValueError(err_msg + '. To proceed without these SNPs, use the flag --allow-missing')
            
    #extract the annotations of the relevant SNPs
    df_annot = df_snps[['PIP']].merge(df_annot, left_index=True, right_index=True)

    #write output to file
    logging.info('Writing annotations to %s'%(args.out))
    if args.out.endswith('.parquet'):
        df_annot.to_parquet(args.out, index=False)
    else:
        df_annot.to_csv(args.out, sep='\t', index=False, float_format='%0.4e')
    
        
    
    
    
