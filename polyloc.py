import numpy as np
import pandas as pd
import os
import sys
import time
import logging
import scipy.stats as stats
from polyfun import PolyFun
from polyfun_utils import configure_logger, set_snpid_index, get_file_name
from polyfun_utils import SNP_COLUMNS
from pyarrow import ArrowIOError
from pyarrow.lib import ArrowInvalid


def splash_screen():
    print('*********************************************************************')
    print('* PolyLoc (POLYgenic LOCalization of complex trait heritability')
    print('* Version 1.0.0')
    print('* (C) 2019-2022 Omer Weissbrod')
    print('*********************************************************************')
    print()
    
    
def check_args(args):
    
    #verify that the requested computations are valid
    mode_params = np.array([args.compute_partitions, args.compute_ldscores, args.compute_polyloc])
    if np.sum(mode_params)==0:
        raise ValueError('must specify at least one of --compute-partitions, --compute-ldscores, --compute-polyloc')
    if args.compute_partitions and args.compute_polyloc and not args.compute_ldscores:
        raise ValueError('cannot use both --compute-partitions and --compute_polyloc without also specifying --compute-ldscores')
    if args.chr is not None:
        if args.compute_partitions or args.compute_polyloc:
            raise ValueError('--chr can only be specified when using only --compute-ldscores')
    if args.bfile_chr is not None:
        if not args.compute_ldscores and not args.compute_partitions:
            raise ValueError('--bfile-chr can only be specified when using --compute-partitions or --compute-ldscores')
    if args.ld_ukb:
        if not args.compute_ldscores and not args.compute_partitions:
            raise ValueError('--ld-ukb can only be specified when using --compute-partitions or --compute-ldscores')
    if args.compute_ldscores and args.compute_polyloc and not args.compute_partitions:
        raise ValueError('cannot use both --compute-ldscores and --compute_polyloc without also specifying --compute-partitions')    
        
    if args.posterior is not None and not args.compute_partitions:
        raise ValueError('--posterior can only be specified together with --compute-partitions')        
    if args.sumstats is not None and not args.compute_polyloc:
        raise ValueError('--sumstats can only be specified together with --compute-polyloc')
    if args.ld_dir is not None and not args.ld_ukb:
        raise ValueError('You cannot specify --ld-dir without also specifying --ld-ukb')
        
    
    #verify partitioning parameters
    if args.skip_Ckmedian and (args.num_bins is None or args.num_bins<=0):
        raise ValueError('You must specify --num-bins when using --skip-Ckmedian')        

    #partitionining-related parameters
    if args.compute_partitions:
        if args.bfile_chr is None:
            #raise ValueError('You must specify --bfile-chr when you specify --compute-partitions')
            logging.warning('You did not specify --bfile-chr with --compute-partitions. PolyLoc will only use SNPs provided in the posterior files')
        if args.posterior is None:
            raise ValueError('--posterior must be specified when using --compute-partitions')
            
    #verify LD-score related parameters
    if args.compute_ldscores:
        if not not args.ld_ukb and args.bfile_chr is None:
            raise ValueError('You must specify either --ld-ukb or --bfile-chr when using --compute-ldscores')
        if not args.ld_ukb and (args.ld_wind_cm is None and args.ld_wind_kb is None and args.ld_wind_snps is None):
            args.ld_wind_cm = 1.0
            logging.warning('no ld-wind argument specified.  PolyLoc will use --ld-cm 1.0')

    if not args.compute_ldscores:
        if not (args.ld_wind_cm is None and args.ld_wind_kb is None and args.ld_wind_snps is None):
            raise ValueError('--ld-wind parameters can only be specified together with --compute-ldscores')
        if args.keep is not None:
            raise ValueError('--keep can only be specified together with --compute-ldscores')
        if args.chr is not None:
            raise ValueError('--chr can only be specified together with --compute-ldscores')

    if args.compute_polyloc:
        if args.sumstats is None:
            raise ValueError('--sumstats must be specified when using --compute-polyloc')    
        if args.w_ld_chr is None:
            raise ValueError('--w-ld-chr must be specified when using --compute-polyloc')    
            
    return args

def check_files(args):

    if args.compute_partitions:
        if not os.path.exists(args.posterior):
            raise IOError('%s not found'%(args.posterior))
        
    #check that required input files exist
    if args.compute_ldscores or args.compute_partitions:
        if args.chr is None: chr_range = range(1,23)            
        else: chr_range = range(args.chr, args.chr+1)
        
        for chr_num in chr_range:
            if args.bfile_chr is not None:
                get_file_name(args, 'bim', chr_num, verify_exists=True)
            if args.compute_ldscores and not args.ld_ukb:
                if args.bfile_chr is None:
                    raise ValueError('--bfile-chr not provided')
                get_file_name(args, 'fam', chr_num, verify_exists=True)
                get_file_name(args, 'bed', chr_num, verify_exists=True)
            if not args.compute_partitions:
                get_file_name(args, 'bins', chr_num, verify_exists=True)
                
    if args.compute_polyloc:    
        for chr_num in range(1,23):
            get_file_name(args, 'w-ld', chr_num, verify_exists=True)
            if not args.compute_partitions:
                get_file_name(args, 'bins', chr_num, verify_exists=True)
            
        

    

class PolyLoc(PolyFun):
    def __init__(self):
        pass
        
    def load_posterior_betas(self, args):
        try:
            df_posterior = pd.read_parquet(args.posterior)
        except (ArrowIOError, ArrowInvalid):
            df_posterior = pd.read_table(args.posterior, sep='\s+')
            
        #preprocess columns
        df_posterior.columns = df_posterior.columns.str.upper()
        
        #make sure that all required columns are found
        has_missing_col = False
        for column in SNP_COLUMNS + ['BETA_MEAN', 'BETA_SD']:
            if column not in df_posterior.columns:
                logging.error('%s has a missing column: %s'%(args.posterior, column))
                has_missing_col = True
        if has_missing_col:
            raise ValueError('%s has missing columns'%(args.posterior))
            
        df_posterior['SNPVAR'] = df_posterior['BETA_MEAN']**2 + df_posterior['BETA_SD']**2
        self.df_snpvar = df_posterior
        
        
    def polyloc_partitions(self, args):
    
        self.load_posterior_betas(args)    
        self.partition_snps_to_bins(args, use_ridge=False)
        self.df_bins = set_snpid_index(self.df_bins)
        
        #add another partition for all SNPs not in the posterior file
        if args.bfile_chr is not None:
            df_bim_list = []
            for chr_num in range(1,23):
                df_bim_chr = pd.read_table(args.bfile_chr+'%d.bim'%(chr_num), sep='\s+', names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'], header=None)
                df_bim_list.append(df_bim_chr)
            df_bim = pd.concat(df_bim_list, axis=0)
            df_bim = set_snpid_index(df_bim)        
            
            #make sure that all variants in the posterior file are also in the plink files
            if np.any(~self.df_bins.index.isin(df_bim.index)):
                raise ValueError('Found variants in posterior file that are not found in the plink files')
                
            #add a new bin for SNPs that are not found in the posterior file (if there are any)
            if df_bim.shape[0] > self.df_bins.shape[0]:
                new_snps = df_bim.index[~df_bim.index.isin(self.df_bins.index)]
                df_bins_new = df_bim.loc[new_snps, SNP_COLUMNS].copy()
                for colname in self.df_bins.drop(columns=SNP_COLUMNS).columns:
                    df_bins_new[colname] = False
                new_colname = 'snpvar_bin%d'%(df_bins_new.shape[1] - len(SNP_COLUMNS)+1)
                self.df_bins[new_colname] = False
                df_bins_new[new_colname] = True
                self.df_bins = pd.concat([self.df_bins, df_bins_new], axis=0)
        
        #save the bins to disk
        self.save_bins_to_disk(args)
        
        #save the bin sizes to disk
        df_binsize = pd.DataFrame(index=np.arange(1,self.df_bins.shape[1] - len(SNP_COLUMNS)+1))
        df_binsize.index.name='BIN'
        df_binsize['BIN_SIZE'] = [self.df_bins[c].sum() for c in self.df_bins.drop(columns=SNP_COLUMNS).columns]   #saves memory
        df_binsize.to_csv(args.output_prefix+'.binsize', sep='\t', index=True)
        
        
    def compute_per_bin_h2(self, prop_h2, prop_h2_jk, df_binsize):
        
        #compute required quantities
        sum_prop_h2 = np.cumsum(prop_h2)
        prop_h2_stderr = np.std(prop_h2_jk, axis=1, ddof=0) * np.sqrt(prop_h2_jk.shape[1]-1)
        sum_prop_h2_jk = np.cumsum(prop_h2_jk, axis=0)
        sum_prop_h2_stderr = np.std(sum_prop_h2_jk, axis=1, ddof=0) * np.sqrt(prop_h2_jk.shape[1]-1)
        
        #create df_bin_h2
        df_bin_h2 = df_binsize.copy()
        df_bin_h2['SUM_BINSIZE'] = df_bin_h2['BIN_SIZE'].cumsum()
        df_bin_h2['%H2'] = prop_h2
        df_bin_h2['%H2_STDERR'] = prop_h2_stderr
        df_bin_h2['SUM_%H2'] = sum_prop_h2
        df_bin_h2['SUM_%H2_STDERR'] = sum_prop_h2_stderr
        
        return df_bin_h2
    
    
    def compute_Mp(self, p, cumsum_prop_h2, cumnum_binsize):

        if p==1: p=0.99999999
        assert np.all(np.any(cumsum_prop_h2 >= p, axis=0))
        num_jk = cumsum_prop_h2.shape[1]
        
        last_bin_index = np.argmax(cumsum_prop_h2 >= p, axis=0)
        num_snps_bin1 = np.zeros(num_jk, dtype=np.int64)
        h2_bin1 = np.zeros(num_jk)
        num_snps_bin1[last_bin_index != 0] = cumnum_binsize[last_bin_index[last_bin_index != 0] - 1]
        h2_bin1[last_bin_index != 0] = cumsum_prop_h2[last_bin_index[last_bin_index != 0] - 1, np.arange(num_jk)[last_bin_index != 0]]
            
        num_snps_bin2 = cumnum_binsize[last_bin_index]        
        h2_bin2 = cumsum_prop_h2[last_bin_index, np.arange(num_jk)]
        slope = (num_snps_bin2-num_snps_bin1).astype(np.float64) / (h2_bin2-h2_bin1)
        assert not np.any(np.isnan(slope))
        Mp = np.ceil(num_snps_bin1 + slope * (p - h2_bin1)).astype(np.int64)        
        
        return Mp
            
    
    def compute_Mp_df(self, prop_h2, prop_h2_jk, cumnum_binsize, outlier_coef=25):    
    
        cumsum_prop_h2 = np.cumsum(prop_h2)
        cumsum_prop_h2_jk = np.cumsum(prop_h2_jk, axis=0)
        dicts_list = []
        for p in np.arange(1,101):
            Mp = self.compute_Mp(p/100., np.reshape(cumsum_prop_h2, (cumsum_prop_h2.shape[0], 1)), cumnum_binsize)[0]
            Mp_jk = self.compute_Mp(p/100., cumsum_prop_h2_jk, cumnum_binsize)
            is_outlier_jk = np.abs(Mp_jk - np.median(Mp_jk)) > stats.iqr(Mp_jk) * outlier_coef
            Mp_jk = Mp_jk[~is_outlier_jk]
            Mp_stderr = np.std(Mp_jk, ddof=0) * np.sqrt(len(Mp_jk)-1)
            dicts_list.append({'p':p, 'Mp':Mp, 'Mp_STDERR':Mp_stderr})
        df_Mp = pd.DataFrame(dicts_list)
        return df_Mp
        
    def compute_polyloc(self, args):
    
        #run S-LDSC and compute taus
        self.run_ldsc(args, use_ridge=False, nn=True, evenodd_split=False, keep_large=True, n_blocks=200)
        hsqhat = self.hsqhat
        jknife = hsqhat.jknife
        taus = jknife.est[0, :hsqhat.n_annot] / hsqhat.Nbar
        taus_jk = jknife.delete_values[:, :hsqhat.n_annot] / hsqhat.Nbar
        
        #load bin sizes
        df_binsize = pd.read_table(args.output_prefix+'.binsize', sep='\t')
        
        #compute prop_h2 for the main analysis and for each jackknife block
        prop_h2 = taus * df_binsize['BIN_SIZE'].values
        prop_h2 /= prop_h2.sum()
        prop_h2_jk = taus_jk.T * df_binsize['BIN_SIZE'].values[:, np.newaxis]
        prop_h2_jk = prop_h2_jk / prop_h2_jk.sum(axis=0)
        
        #compuute per-bin h^2 and stderrs
        df_bin_h2 = self.compute_per_bin_h2(prop_h2, prop_h2_jk, df_binsize)
        
        #compute df_Mp
        cumnum_binsize = df_binsize['BIN_SIZE'].cumsum().values
        df_Mp = self.compute_Mp_df(prop_h2, prop_h2_jk, cumnum_binsize)
        
        #write df_bin_h2 to output file
        outfile = args.output_prefix+'.bin_h2'
        df_bin_h2.to_csv(outfile, sep='\t', float_format='%0.5f', index=False)
        logging.info('Wrote per-bin heritability to %s'%(outfile))
        
        #write df_Mp to output file
        outfile = args.output_prefix+'.Mp'
        df_Mp.to_csv(outfile, sep='\t', float_format='%0.5f', index=False)
        logging.info('Wrote Mp estimates to %s'%(outfile))
        
        

    def polyloc_main(self, args):
    
        #compute snp variances using L2-regularized S-LDSC with an odd/even chromosome split
        if args.compute_partitions:
            self.polyloc_partitions(args)
            
        #compute LD-scores of SNP partitions
        if args.compute_ldscores:
            self.compute_ld_scores(args)
        
        #compute polygenic localization
        if args.compute_polyloc:
            self.compute_polyloc(args)
            
            
            
    
        
    
        


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()

    #partitioning-related parameters
    parser.add_argument('--num-bins', type=int, default=None, help='Number of bins to partition SNPs into. If not specified, PolyLoc will automatically select this number based on a BIC criterion')
    parser.add_argument('--skip-Ckmedian', default=False, action='store_true', help='If specified, use a regular K-means algorithm instead of the R Ckmeans.1d.dp package')
    
    #mode related parameters
    parser.add_argument('--compute-partitions', default=False, action='store_true', help='If specified, PolyLoc will compute per-SNP h2 using L2-regularized S-LDSC')
    parser.add_argument('--compute-ldscores', default=False, action='store_true', help='If specified, PolyLoc will compute LD-scores of SNP bins')
    parser.add_argument('--compute-polyloc', default=False, action='store_true', help='If specified, PolyLoc will perform polygenic localization of SNP heritability')
    
    #ld-score related parameters
    parser.add_argument('--chr', type=int, default=None, help='Chromosome number (only applicable when only specifying --ldscores). If not set, PolyLoc will compute LD-scores for all chromosomes')
    #parser.add_argument('--npz-prefix', default=None, help='Prefix of npz files that encode LD matrices (used to compute LD-scores)')
    parser.add_argument('--ld-wind-cm', type=float, default=None, help='window size to be used for estimating LD-scores in units of centiMorgans (cM).')
    parser.add_argument('--ld-wind-kb', type=int, default=None, help='window size to be used for estimating LD-scores in units of Kb.')
    parser.add_argument('--ld-wind-snps', type=int, default=None, help='window size to be used for estimating LD-scores in units of SNPs.')
    parser.add_argument('--chunk-size',  type=int, default=50, help='chunk size for LD-scores calculation')
    parser.add_argument('--keep',  default=None, help='File with ids of individuals to use when computing LD-scores')
    
    #data input/output parameters
    parser.add_argument('--sumstats', help='Input summary statistics file')
    parser.add_argument('--posterior', help='Input file with posterior means and variances of causal effect sizes')
    parser.add_argument('--w-ld-chr', help='Suffix of LD-score weights files (as in ldsc)')
    parser.add_argument('--bfile-chr', default=None, help='Prefix of plink files (used to compute LD-scores)')
    parser.add_argument('--output-prefix', required=True, help='Prefix of all PolyLoc out file names')    
    parser.add_argument('--ld-ukb', default=False, action='store_true', help='If specified, PolyLoc will use UKB LD matrices to compute LD-scores')
    parser.add_argument('--ld-dir', default=None, help='The path of a directory with UKB LD files (if not specified PolyLoc will create a temporary directory)')
    
    
    
    #LDSC parameters
    parser.add_argument('--nnls-exact', default=False, action='store_true', help='If specified, S-LDSC will estimate non-negative taus using an exact instead of an approximate solver (this will be slower but slightly more accurate)')
    
    #show splash screen
    splash_screen()

    #extract args
    args = parser.parse_args()
    
    #check that the output directory exists
    if len(os.path.dirname(args.output_prefix))>0 and not os.path.exists(os.path.dirname(args.output_prefix)):
        raise ValueError('output directory %s doesn\'t exist'%(os.path.dirname(args.output_prefix)))    
    
    #configure logger
    configure_logger(args.output_prefix)
        
    #check and fix args
    args = check_args(args)
    check_files(args)
    args.anno = None
    
    #create and run PolyLoc object
    polyloc_obj = PolyLoc()
    polyloc_obj.polyloc_main(args)
    
    print()
    
