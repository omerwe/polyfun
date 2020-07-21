import numpy as np; np.set_printoptions(precision=4, linewidth=200)
import pandas as pd; pd.set_option('display.width', 200)
import os
import logging
from tqdm import tqdm
from polyfun import configure_logger, check_package_versions
from polyfun_utils import set_snpid_index
from pyarrow import ArrowIOError
from pyarrow.lib import ArrowInvalid
from polyfun_utils import DEFAULT_REGIONS_FILE


def main(args):
    
    #read sumstats file
    try:
        df_sumstats = pd.read_parquet(args.sumstats)
    except (ArrowIOError, ArrowInvalid):
        df_sumstats = pd.read_table(args.sumstats, delim_whitespace=True)
        
    #read regions file
    df_regions = pd.read_table(args.regions_file)
    df_regions = df_regions.loc[df_regions.apply(lambda r: np.any((df_sumstats['CHR']==r['CHR']) & (df_sumstats['BP'].between(r['START'], r['END']))), axis=1)]
    
    #aggregate outputs
    df_sumstats_list = []
    logging.info('Aggregating results...')
    for _, r in tqdm(df_regions.iterrows()):
        chr_num, start, end, url_prefix = r['CHR'], r['START'], r['END'], r['URL_PREFIX']
        output_file_r = '%s.chr%s.%s_%s.gz'%(args.out_prefix, chr_num, start, end)
        if not os.path.exists(output_file_r):
            err_msg = 'output file for chromosome %d bp %d-%d doesn\'t exist'%(chr_num, start, end)
            if args.allow_missing_jobs:
                logging.warning(err_msg)
                continue
            else:
                raise IOError(err_msg)
        df_sumstats_r = pd.read_table(output_file_r)
        
        #mark distance from center
        middle = (start+end)//2
        df_sumstats_r['DISTANCE_FROM_CENTER'] = np.abs(df_sumstats_r['BP'] - middle)
        df_sumstats_list.append(df_sumstats_r)
    
    #keep only the most central result for each SNP
    df_sumstats = pd.concat(df_sumstats_list, axis=0)
    df_sumstats.sort_values('DISTANCE_FROM_CENTER', inplace=True, ascending=True)
    df_sumstats = set_snpid_index(df_sumstats, allow_duplicates=True)
    df_sumstats = df_sumstats.loc[~df_sumstats.index.duplicated(keep='first')]
    del df_sumstats['DISTANCE_FROM_CENTER']
    df_sumstats.sort_values(['CHR', 'BP'], inplace=True, ascending=True)
    
    #write output file
    df_sumstats.to_csv(args.out, sep='\t', index=False)
    logging.info('Wrote aggregated results to %s'%(args.out))
        


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    #general parameters
    parser.add_argument('--sumstats', required=True, help='Name of sumstats file')
    parser.add_argument('--out-prefix', required=True, help='prefix of output files')
    parser.add_argument('--out', required=True, help='name of the aggregated output files')
    parser.add_argument('--allow-missing-jobs', default=False, action='store_true', help='whether to allow missing jobs')
    parser.add_argument('--regions-file', default=DEFAULT_REGIONS_FILE, help='name of file of regions and their URLs')
    
    #check package versions
    check_package_versions()

    #extract args
    args = parser.parse_args()
    
    #check that the output directory exists
    if len(os.path.dirname(args.out))>0 and not os.path.exists(os.path.dirname(args.out)):
        raise ValueError('output directory %s doesn\'t exist'%(os.path.dirname(args.out)))    

    #configure logger
    configure_logger(args.out_prefix)

    #invoke main function
    main(args)
    