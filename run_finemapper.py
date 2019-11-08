import numpy as np
import pandas as pd
from finemapper import SUSIE_Wrapper
from polyfun import configure_logger, check_package_versions
import logging
import os
import scipy.sparse as sparse


def splash_screen():
    print('*********************************************************************')
    print('* Fine-mapping Wrapper')
    print('* Version 1.0.0')
    print('* (C) 2019 Omer Weissbrod')
    print('*********************************************************************')
    print()
    
    
    
def load_ld(ld_prefix):
        
    #load SNPs info
    snps_filename_parquet = ld_prefix+'.parquet'
    snps_filename_gz = ld_prefix+'.gz'
    if os.path.exists(snps_filename_parquet):
        df_ld_snps = pd.read_parquet(snps_filename_parquet)
    elif os.path.exists(snps_filename_gz):
        df_ld_snps = pd.read_table(snps_filename_gz, delim_whitespace=True)
    else:
        raise ValueError('couldn\'t find SNPs file %s or %s'%(snps_filename_parquet, snps_filename_gz))
        
    #load LD matrix
    R_filename = ld_prefix+'.npz'
    if not os.path.exists(R_filename):
        raise IOError('%s not found'%(R_filename))
    R = sparse.load_npz(R_filename).toarray()
    R = R+R.T
    assert np.allclose(np.diag(R), 1.0)
    logging.info('Loaded an LD matrix for %d SNPs from %s'%(R.shape[0], R_filename))
    
    #sanity checks
    assert R.shape[0] == R.shape[1]
    if R.shape[0] != df_ld_snps.shape[0]:
        raise ValueError('LD matrix has a different number of SNPs than the SNPs file')
    
    return R, df_ld_snps
    
    

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()

    #general parameters
    parser.add_argument('--method', required=True, help='Fine-mapping method (currently only susie is supported)')
    parser.add_argument('--sumstats', required=True, help='Name of sumstats file')
    parser.add_argument('--chr', required=True, type=int, help='Target chromosome')
    parser.add_argument('--start', required=True, type=int, help='First base-pair in the region to finemap')
    parser.add_argument('--end', required=True, type=int, help='Last base-pair in the region to finemap')
    parser.add_argument('--n', required=True, type=int, help='Sample size')
    parser.add_argument('--geno', default=None, help='Genotypes file (plink or bgen format)')
    parser.add_argument('--ld', default=None, help='prefix of LD matrix file names')
    
    #LDstore related parameters
    parser.add_argument('--ldstore', default=None, help='Path to an LDstore executable file')
    parser.add_argument('--threads', type=int, default=None, help='The number of CPU cores LDstore will use (if not specified, LDstore will use the max number of CPU cores available')
    parser.add_argument('--cache-dir', default=None, help='If specified, this is a path of a directory that will cache LD matrices that have already been computed')
    
    parser.add_argument('--max-num-causal', required=True, type=int, help='Number of causal SNPs')
    parser.add_argument('--non-funct', action='store_true', default=False, help='Perform non-functionally informed fine-mapping')
    parser.add_argument('--hess', action='store_true', default=False, help='If specified, estimate causal effect variance via HESS')
    parser.add_argument('--verbose', action='store_true', default=False, help='If specified, show verbose output')
    
    parser.add_argument('--sample-file', default=None, help='BGEN files must be used together with a sample file')
    parser.add_argument('--incl-samples', default=None, help='A single-column text file specifying the ids of individuals to exclude from fine-mapping')
    parser.add_argument('--out', required=True, help='name of the output file')
    
    #check package versions
    check_package_versions()
    
    #show splash screen
    splash_screen()
    
    #extract args
    args = parser.parse_args()

    #configure logger
    configure_logger(args.out)
    
    #check params
    if args.geno is None:
        if args.ld is None:
            raise ValueError('must specify either --geno or --ld')
        if args.ldstore is not None:
            raise ValueError('cannot specify both --ld and --ldstore')
    if args.geno is not None:
        if args.ld is not None:
            raise ValueError('cannot specify both --geno and --ld')
        if args.ldstore is None:
            raise ValueError('You must specify --ldstore together with --geno')
    
        
    #load LD matrix if requested
    if args.ld is None:
        ld, df_ld_snps = None, None
    else:
        ld, df_ld_snps = load_ld(args.ld)
    
    #start fine-mapping
    if args.method == 'susie':
        finemap_obj = SUSIE_Wrapper(genotypes_file=args.geno, sumstats_file=args.sumstats, n=args.n, chr_num=args.chr, 
                                    sample_file=args.sample_file, incl_samples=args.incl_samples,
                                    ldstore_exe=args.ldstore, n_threads=args.threads,
                                    cache_dir=args.cache_dir)
        df_finemap = finemap_obj.finemap(locus_start=args.start, locus_end=args.end, num_causal_snps=args.max_num_causal,
                     use_prior_causal_prob=not args.non_funct, prior_var=None, residual_var=None, hess=args.hess,
                     verbose=args.verbose, ld=ld, df_ld_snps=df_ld_snps)
    elif args.method == 'finemap':
        raise ValueError('FINEMAP is not yet supported')
    else:
        raise ValueError('unknown method specified in --method')
    logging.info('Writing fine-mapping results to %s'%(args.out))
    df_finemap.sort_values('PIP', ascending=False, inplace=True)
    df_finemap.to_csv(args.out, sep='\t', index=False, float_format='%0.5e')

