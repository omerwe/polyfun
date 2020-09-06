import numpy as np; np.set_printoptions(precision=4, linewidth=200)
import pandas as pd; pd.set_option('display.width', 200)
import os
import logging
import scipy.stats as stats
from polyfun import configure_logger, check_package_versions
from pyarrow import ArrowIOError
from pyarrow.lib import ArrowInvalid
from polyfun_utils import DEFAULT_REGIONS_FILE


FINEMAPPER_SCRIPT = 'finemapper.py'


def create_finemapper_cmd(args, chr_num, start, end, url_prefix):

    output_file = '%s.chr%s.%s_%s.gz'%(args.out_prefix, chr_num, start, end)
    cmd = '%s %s --chr %s --start %s --end %s --out %s'%(args.python3, FINEMAPPER_SCRIPT, chr_num, start, end, output_file)
    if args.max_num_causal>1 and args.geno is None:
        cmd += ' --ld %s'%(url_prefix)
    
    #add command line arguments
    for key, value in vars(args).items():
        if key in ['python3', 'regions_file', 'out_prefix', 'jobs_file', 'chr', 'pvalue_cutoff']: continue
        key = key.replace('_', '-')
        if type(value)==bool:
            if value:
                cmd += ' --%s'%(key)                
        elif value is not None:
            cmd += ' --%s %s'%(key, value)
        
    return cmd
    

def main(args):
    
    #read sumstats file
    try:
        df_sumstats = pd.read_parquet(args.sumstats)
    except (ArrowIOError, ArrowInvalid):
        df_sumstats = pd.read_table(args.sumstats, sep='\s+')
        
    #compute p-values if needed
    if args.pvalue_cutoff is not None:
        df_sumstats['P'] = stats.chi2(1).sf(df_sumstats['Z']**2)

    #read regions file
    df_regions = pd.read_table(args.regions_file)
    if args.chr is not None:
        df_regions = df_regions.query('CHR==%d'%(args.chr))
        if df_regions.shape[0]==0: raise ValueError('no SNPs found in chromosome %d'%(args.chr))
    df_regions = df_regions.loc[df_regions.apply(lambda r: np.any((df_sumstats['CHR']==r['CHR']) & (df_sumstats['BP'].between(r['START'], r['END']))), axis=1)]
    
    
    
    #create jobs
    with open(args.jobs_file, 'w') as f:
        for _, r in df_regions.iterrows():
            chr_num, start, end, url_prefix = r['CHR'], r['START'], r['END'], r['URL_PREFIX']
            
            #apply p-value filter if needed
            if args.pvalue_cutoff is not None:
                df_sumstats_r = df_sumstats.query('CHR==%d & %d <= BP <= %d'%(chr_num, start, end))
                if np.all(df_sumstats_r['P'] > args.pvalue_cutoff): continue
            
            #create and write the fine-mapping command for this region
            cmd = create_finemapper_cmd(args, chr_num, start, end, url_prefix)
            f.write(cmd + '\n')
    
    logging.info('Wrote fine-mapping commands to %s'%(args.jobs_file))
        


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    #general parameters
    parser.add_argument('--method', required=True, help='Fine-mapping method (currently susie and finemap are supported)')
    parser.add_argument('--sumstats', required=True, help='Name of sumstats file')
    parser.add_argument('--n', required=True, type=int, help='Sample size')
    parser.add_argument('--geno', default=None, help='Genotypes file (plink or bgen format)')
    parser.add_argument('--chr', default=None, type=int, help='Target chromosome (if not provided, all chromosomes will be considered)')
    
    #LDstore related parameters
    parser.add_argument('--finemap-exe', default=None, help='Path to FINEMAP v1.4 executable file')
    parser.add_argument('--memory', type=int, default=1, help='Maximum amount of memory in GB to allocate to LDStore')
    parser.add_argument('--threads', type=int, default=None, help='The number of CPU cores LDstore will use (if not specified, LDstore will use the max number of CPU cores available')
    
    parser.add_argument('--max-num-causal', required=True, type=int, help='Number of causal SNPs')
    parser.add_argument('--non-funct', action='store_true', default=False, help='Perform non-functionally informed fine-mapping')
    parser.add_argument('--hess', action='store_true', default=False, help='If specified, estimate causal effect variance via HESS')
    parser.add_argument('--verbose', action='store_true', default=False, help='If specified, show verbose output')
    parser.add_argument('--allow-missing', default=False, action='store_true', help='If specified, SNPs with sumstats that are not \
                            found in the LD panel will be omitted. This is not recommended, because the omitted SNPs may be causal,\
                            which could lead to false positive results')
    
    parser.add_argument('--regions-file', default=DEFAULT_REGIONS_FILE, help='name of file of regions and their URLs')
    parser.add_argument('--python3', default='python3', help='python3 executable')
    parser.add_argument('--out-prefix', required=True, help='prefix of the output files')
    parser.add_argument('--jobs-file', required=True, help='name of file with fine-mapping commands')
    parser.add_argument('--pvalue-cutoff', type=float, default=None, help='only consider regions that have at least one SNP with a p-value greater than this cutoff')
    
    #check package versions
    check_package_versions()

    #extract args
    args = parser.parse_args()
    
    #check that the output directory exists
    if len(os.path.dirname(args.out_prefix))>0 and not os.path.exists(os.path.dirname(args.out_prefix)):
        raise ValueError('output directory %s doesn\'t exist'%(os.path.dirname(args.out_prefix)))    
    if len(os.path.dirname(args.jobs_file))>0 and not os.path.exists(os.path.dirname(args.jobs_file)):
        raise ValueError('output directory %s doesn\'t exist'%(os.path.dirname(args.jobs_file)))    

    #configure logger
    configure_logger(args.out_prefix)

    #invoke main function
    main(args)
    