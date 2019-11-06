import numpy as np
import pandas as pd
from finemapper import SUSIE_Wrapper
from polyfun import configure_logger
import logging

def splash_screen():
    print('*********************************************************************')
    print('* Fine-mapping Wrapper')
    print('* Version 1.0.0')
    print('* (C) 2019 Omer Weissbrod')
    print('*********************************************************************')
    print()
    

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()

    #general parameters
    parser.add_argument('--method', required=True, help='Fine-mapping method (currently only susie is supported)')
    parser.add_argument('--geno', required=True, help='Genotypes file (plink or bgen format)')
    parser.add_argument('--sumstats', required=True, help='Name of sumstats file')
    parser.add_argument('--chr', required=True, type=int, help='Target chromosome')
    parser.add_argument('--start', required=True, type=int, help='First base-pair in the region to finemap')
    parser.add_argument('--end', required=True, type=int, help='Last base-pair in the region to finemap')
    parser.add_argument('--n', required=True, type=int, help='Sample size')
    
    #LDstore related parameters
    parser.add_argument('--ldstore', required=True, help='Path to an LDstore executable file')
    parser.add_argument('--threads', type=int, default=None, help='The number of CPU cores LDstore will use (if not specified, LDstore will use the max number of CPU cores available')
    parser.add_argument('--cache-dir', default=None, help='If specified, this is a path of a directory that will cache LD matrices that have already been computed')
    
    parser.add_argument('--max-num-causal', required=True, type=int, help='Number of causal SNPs')
    parser.add_argument('--non-funct', action='store_true', default=False, help='Perform non-functionally informed fine-mapping')
    parser.add_argument('--hess', action='store_true', default=False, help='If specified, estimate causal effect variance via HESS')
    parser.add_argument('--verbose', action='store_true', default=False, help='If specified, show verbose output')
    
    parser.add_argument('--sample-file', default=None, help='BGEN files must be used together with a sample file')
    parser.add_argument('--incl-samples', default=None, help='A single-column text file specifying the ids of individuals to exclude from fine-mapping')
    parser.add_argument('--out', required=True, help='name of the output file')
    
    #show splash screen
    splash_screen()
    
    #extract args
    args = parser.parse_args()

    #configure logger
    configure_logger(args.out)
    
    
    if args.method == 'susie':
        #import ipdb; ipdb.set_trace()
        finemap_obj = SUSIE_Wrapper(args.geno, args.sumstats, n=args.n, chr_num=args.chr, sample_file=args.sample_file, incl_samples=args.incl_samples, ldstore_exe=args.ldstore, n_threads=args.threads, cache_dir=args.cache_dir)
        df_finemap = finemap_obj.finemap(args.start, args.end, num_causal_snps=args.max_num_causal, use_prior_causal_prob=not args.non_funct, prior_var=None, residual_var=None, hess=args.hess, verbose=args.verbose)
    elif args.method == 'finemap':
        raise ValueError('FINEMAP is not yet supported')
    else:
        raise ValueError('unknown method specified in --method')
    logging.info('Writing fine-mapping results to %s'%(args.out))
    df_finemap.sort_values('PIP', ascending=False, inplace=True)
    df_finemap.to_csv(args.out, sep='\t', index=False, float_format='%0.5e')

