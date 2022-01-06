import numpy as np
import pandas as pd
import os
import time
from ldsc_polyfun import jackknife, regressions, sumstats, ldscore, parse
import logging
from copy import deepcopy
from tqdm import tqdm
from polyfun_utils import Logger, check_package_versions, set_snpid_index, configure_logger, get_file_name
from polyfun_utils import SNP_COLUMNS
from pyarrow import ArrowIOError
from pyarrow.lib import ArrowInvalid
from compute_ldscores_from_ld import compute_ldscores_chr
import tempfile


MAX_CHI2=80


def __filter__(fname, noun, verb, merge_obj):
    merged_list = None
    if fname:
        f = lambda x,n: x.format(noun=noun, verb=verb, fname=fname, num=n)
        x = parse.FilterFile(fname)
        c = 'Read list of {num} {noun} to {verb} from {fname}'
        logging.info(f(c, len(x.IDList)))
        merged_list = merge_obj.loj(x.IDList)
        len_merged_list = len(merged_list)
        if len_merged_list > 0:
            c = 'After merging, {num} {noun} remain'
            logging.info(f(c, len_merged_list))
        else:
            error_msg = 'No {noun} retained for analysis'
            raise ValueError(f(error_msg, 0))

        return merged_list


def splash_screen():
    print('*********************************************************************')
    print('* PolyFun (POLYgenic FUNctionally-informed fine-mapping)')
    print('* Version 1.0.0')
    print('* (C) 2019-2021 Omer Weissbrod')
    print('*********************************************************************')
    print()
    
    
def check_args(args):    
    
    #verify that the requested computations are valid
    mode_params = np.array([args.compute_h2_L2, args.compute_ldscores, args.compute_h2_bins])
    if np.sum(mode_params)==0:
        raise ValueError('must specify at least one of --compute-h2-L2, --compute-ldscores, --compute-h2-bins')
    if args.compute_h2_L2 and args.compute_h2_bins and not args.compute_ldscores:
        raise ValueError('cannot use both --compute-h2_L2 and --compute_h2_bins without also specifying --compute-ldscores')
    if args.chr is not None:
        if args.compute_h2_L2 or args.compute_h2_bins:
            raise ValueError('--chr can only be specified when using only --compute-ldscores')
    if args.bfile_chr is not None:
        if not args.compute_ldscores:
            raise ValueError('--bfile-chr can only be specified when using --compute-ldscores')
    if args.ld_ukb:
        if not args.compute_ldscores:
            raise ValueError('--ld-ukb can only be specified when using --compute-ldscores')
    if args.no_partitions:
        if not args.compute_h2_L2:
            raise ValueError('cannot specify --no-partitions without specifying --compute-h2-L2')
        if args.compute_ldscores:
            raise ValueError('cannot specify both --no-partitions and --compute-ldscores')    
        if args.compute_h2_bins:
            raise ValueError('cannot specify both --no-partitions and --compute-h2-bins')
    if args.compute_ldscores and args.compute_h2_bins and not args.compute_h2_L2:
        raise ValueError('cannot use both --compute-ldscores and --compute_h2_bins without also specifying --compute-h2-L2')    
    
    #verify partitioning parameters
    if args.skip_Ckmedian and (args.num_bins is None or args.num_bins<=0):
        raise ValueError('You must specify --num-bins when using --skip-Ckmedian')        
    
    #verify LD-score related parameters
    if args.ld_dir is not None and not args.ld_ukb:
        raise ValueError('You cannot specify --ld-dir without also specifying --ld-ukb')
    if args.bfile_chr is not None and args.ld_ukb:
        raise ValueError('You can specify only one of --bfile-chr and --ld-ukb')

    if args.compute_ldscores:
        if args.bfile_chr is None and not args.ld_ukb:
            raise ValueError('You must specify either --bfile-chr or --ld-ukb when you specify --compute-ldscores')    
        if not args.ld_ukb and (args.ld_wind_cm is None and args.ld_wind_kb is None and args.ld_wind_snps is None):
            args.ld_wind_cm = 1.0
            logging.warning('no ld-wind argument specified.  PolyFun will use --ld-cm 1.0')
            
    if not args.compute_ldscores:
        if not (args.ld_wind_cm is None and args.ld_wind_kb is None and args.ld_wind_snps is None):
            raise ValueError('--ld-wind parameters can only be specified together with --compute-ldscores')
        if args.keep is not None:
            raise ValueError('--keep can only be specified together with --compute-ldscores')
        if args.chr is not None:
            raise ValueError('--chr can only be specified together with --compute-ldscores')
        
            
    
    if args.compute_h2_L2:
        if args.sumstats is None:
            raise ValueError('--sumstats must be specified when using --compute-h2-L2')
        if args.ref_ld_chr is None:
            raise ValueError('--ref-ld-chr must be specified when using --compute-h2-L2')
        if args.w_ld_chr is None:
            raise ValueError('--w-ld-chr must be specified when using --compute-h2-L2')
            
    if args.compute_h2_bins:
        if args.sumstats is None:
            raise ValueError('--sumstats must be specified when using --compute-h2-bins')
        if args.w_ld_chr is None:
            raise ValueError('--w-ld-chr must be specified when using --compute-h2-bins')
        if args.ref_ld_chr is not None and not args.compute_ldscores:
            raise ValueError('--ref-ld-chr should not be specified when using --compute-h2-bins, unless you also use --compute-ldscores')
            
            
    return args

def check_files(args):
        
    #check that required input files exist
    if args.compute_h2_L2:
        if not os.path.exists(args.sumstats):
            raise IOError('Cannot find sumstats file %s'%(args.sumstats))
        for chr_num in range(1,23):
            get_file_name(args, 'ref-ld', chr_num, verify_exists=True, allow_multiple=True)
            get_file_name(args, 'w-ld', chr_num, verify_exists=True)
            get_file_name(args, 'annot', chr_num, verify_exists=True, allow_multiple=True)
            
    if args.compute_ldscores:
        if args.chr is None: chr_range = range(1,23)            
        else: chr_range = range(args.chr, args.chr+1)
        
        for chr_num in chr_range:
            if args.bfile_chr is not None:
                get_file_name(args, 'bim', chr_num, verify_exists=True)
                get_file_name(args, 'fam', chr_num, verify_exists=True)
                get_file_name(args, 'bed', chr_num, verify_exists=True)
            if not args.compute_h2_L2:
                get_file_name(args, 'snpvar_ridge', chr_num, verify_exists=True)
                get_file_name(args, 'bins', chr_num, verify_exists=True)
                
    if args.compute_h2_bins and not args.compute_ldscores:
        for chr_num in range(1,23):
            get_file_name(args, 'w-ld', chr_num, verify_exists=True)
            if not args.compute_h2_L2:
                get_file_name(args, 'bins', chr_num, verify_exists=True)
            

        

    
    
    
class PolyFun:
    def __init__(self):
        pass
        
    def run_ldsc(self, args, use_ridge, nn, keep_large, evenodd_split, n_blocks=2):

        #prepare LDSC objects
        log = Logger()
        args.h2 = args.sumstats
        args.ref_ld = None
        args.w_ld = None
        args.n_blocks = n_blocks
        args.M = None
        args.not_M_5_50 = True
        
        #if not ridge, the we'll use the LD-scores of our bins
        if not use_ridge:
            args = deepcopy(args)
            args.ref_ld_chr = args.output_prefix+'.'
        
        #read input data
        if use_ridge or not args.compute_ldscores or True:
            M_annot, w_ld_cname, ref_ld_cnames, df_sumstats, _ = sumstats._read_ld_sumstats(args, log, args.h2)
        else:
            #TODO: Don't reload files if we don't have to...
            M_annot = self.M
            w_ld_cname = 'w_ld'
            ref_ld_cnames = self.df_bins.columns
            try:
                df_sumstats = pd.read_parquet(args.sumstats)            
            except (ArrowIOError, ArrowInvalid):
                df_sumstats = pd.read_table(args.sumstats, sep='\s+')            
            ###merge everything together...
            
        #prepare LD-scores for S-LDSC run
        ref_ld = np.array(df_sumstats[ref_ld_cnames], dtype=np.float32)
        sumstats._check_ld_condnum(args, log, ref_ld_cnames)
        if df_sumstats.shape[0] < 200000:
            logging.warning('number of SNPs is smaller than 200k; this is almost always bad.')
        n_snp = len(df_sumstats)
        n_blocks = np.minimum(n_snp, args.n_blocks)
        n_annot = len(ref_ld_cnames)
        if n_annot<=1:
            raise ValueError('Only one annotation found')
        chisq_max = max(0.001*df_sumstats['N'].max(), MAX_CHI2)

        #prepare chi2 statistics
        s = lambda x: np.array(x).reshape((n_snp, 1))
        chisq = s(df_sumstats.Z**2).astype(np.float32)
        ii = np.ravel(chisq < chisq_max)
        df_sumstats = df_sumstats.loc[ii, :]
        if np.any(~ii):
            logging.info('Removed {M} SNPs with chi^2 > {C} ({N} SNPs remain)'.format(
                    C=chisq_max, N=np.sum(ii), M=n_snp-np.sum(ii)))
        n_snp = np.sum(ii)  # lambdas are late-binding, so this works
        ref_ld = np.array(df_sumstats[ref_ld_cnames], dtype=np.float32)
        chisq = chisq[ii].reshape((n_snp, 1))

        #Run S-LDSC
        self.ref_ld_cnames = [c for c in ref_ld_cnames.str[:-2] if c not in SNP_COLUMNS]
        hsqhat = regressions.Hsq(chisq, 
            ref_ld,
            s(df_sumstats[w_ld_cname]),
            s(df_sumstats.N),
            M_annot, n_blocks=n_blocks, intercept=None,
            twostep=None, old_weights=True,
            chr_num=df_sumstats['CHR'],
            loco=use_ridge, ridge_lambda=None,
            standardize_ridge=True,
            approx_ridge=True,
            num_chr_sets=2,
            evenodd_split=evenodd_split,
            nn=nn,
            keep_large=keep_large,
            nnls_exact=args.nnls_exact
            )
            
        #save the results object
        if use_ridge:
            self.hsqhat_ridge = hsqhat
        else:
            self.hsqhat = hsqhat
            

    
    def load_annotations_file(self, args, chr_num, use_ridge):
        
        #load annotations file for this chromosome
        if use_ridge:
            annot_filenames = get_file_name(args, 'annot', chr_num, allow_multiple=True)
        else:
            annot_filenames = [get_file_name(args, 'bins', chr_num)]
        
        #load annotation file(s)
        df_annot_chr_list = []
        for annot_filename in annot_filenames:
            try:
                df_annot_chr = pd.read_parquet(annot_filename)
            except (ArrowIOError, ArrowInvalid):
                df_annot_chr = pd.read_table(annot_filename)
            df_annot_chr_list.append(df_annot_chr)
        if len(df_annot_chr_list)==1:
            df_annot_chr = df_annot_chr_list[0]
        else:
            for df in df_annot_chr_list[1:]:
                for snp_col in SNP_COLUMNS:
                    if (df.shape[0] != df_annot_chr_list[0].shape[0]) or (np.any(df[snp_col] != df_annot_chr_list[0][snp_col])):
                        raise ValueError('Different annotation files of chromosome %d must be perfectly aligned'%(chr_num))
                df.drop(columns=['CM'], inplace=True, errors='ignore')
                df.drop(columns=SNP_COLUMNS, inplace=True, errors='raise')
            df_annot_chr = pd.concat(df_annot_chr_list, axis=1)
        
        #make sure all required columns were found
        df_annot_chr.drop(columns=['CM'], inplace=True, errors='ignore')
        found_missing_col = False
        for colname in SNP_COLUMNS:
            if colname not in df_annot_chr.columns:
                logging.error('%s has a missing column: %s'%(annot_filename, colname))
                found_missing_col = True
        if found_missing_col:
            raise ValueError('Missing columns found in %s'%(annot_filename))
            
        #subset annotations if requested
        if args.anno is not None:
            anno_to_use = args.anno.split(',')
            assert np.all(np.isin(anno_to_use, df_annot_chr.columns))
            df_annot_chr = df_annot_chr[SNP_COLUMNS + anno_to_use]
            
        #if we have more annotations that ref-ld, it might mean that some annotations were removed, so remove them from here as well
        if not np.all(np.isin(self.ref_ld_cnames, df_annot_chr.columns)):
            raise ValueError('Annotation names in annotations file do not match the one in the LD-scores file')
        if len(self.ref_ld_cnames) < len(df_annot_chr.columns) - len(SNP_COLUMNS):            
            df_annot_chr = df_annot_chr[SNP_COLUMNS + self.ref_ld_cnames]

        #make sure that we get the same columns as the ones in the LD-score files
        if not np.all([c for c in df_annot_chr.columns if c not in SNP_COLUMNS ]== self.ref_ld_cnames):
            raise ValueError('Annotation names in annotations file do not match the one in the LD-scores file')            

        return df_annot_chr


    def compute_snpvar_chr(self, args, chr_num, use_ridge):        

        #load annotations file from disk
        df_annot_chr = self.load_annotations_file(args, chr_num, use_ridge)
            
        #extract taus from a jknife object
        if use_ridge:
            hsqhat = self.hsqhat_ridge
            jknife = hsqhat.jknife_ridge
            
            #make sure that the chromosome exists in one set
            found_chrom = np.any([chr_num in chr_set for chr_set in jknife.chromosome_sets])
            if not found_chrom:
                raise ValueError('not all chromosomes have a taus estimate - please make sure that the intersection of SNPs with sumstats and with annotations data spans all 22 human chromosomes')

            #find the relevant set number
            set_num=None
            for chr_set_i, chr_set in enumerate(jknife.chromosome_sets):
                if chr_num not in chr_set:
                    assert set_num is None
                    set_num = chr_set_i
            if set_num is None:
                raise ValueError('Could not find Ridge predictions for chromosome %d'%(chr_num))
            
            #compute and return snpvar
            taus = jknife.est_loco_ridge[set_num][:hsqhat.n_annot] / hsqhat.Nbar

        else:
            hsqhat = self.hsqhat
            jknife = hsqhat.jknife
            if len(jknife.est_loco) != 22:
                raise ValueError('not all chromosomes have a taus estimate - please make sure that the intersection of SNPs with sumstats and with annotations data spans all 22 human chromosomes')
            taus = jknife.est_loco[chr_num-1][:hsqhat.n_annot] / hsqhat.Nbar
            
        #save the taus to disk
        taus_output_file = get_file_name(args, ('taus_ridge' if use_ridge else 'taus_nn'), chr_num, verify_exists=False)
        df_taus = pd.Series(taus, index=df_annot_chr.drop(columns=SNP_COLUMNS, errors='raise').columns)
        df_taus.index.name = 'ANNOTATION'
        df_taus.name = 'ANNOTATION_COEFFICIENT'
        df_taus.to_csv(taus_output_file, header=True, index=True, sep='\t')
            
        #compute and return the snp variances
        df_snpvar_chr = df_annot_chr.drop(columns=SNP_COLUMNS, errors='raise').dot(taus)
        df_snpvar_chr = df_snpvar_chr.to_frame(name='SNPVAR')
        df_snpvar_chr = pd.concat((df_annot_chr[SNP_COLUMNS], df_snpvar_chr), axis=1)
        return df_snpvar_chr
        

    def compute_snpvar(self, args, use_ridge):
        logging.info('Computing per-SNP h^2 for each chromosome...')
        
        #iterate over chromosomes
        df_snpvar_chr_list = []
        for chr_num in tqdm(range(1,23)):
            df_snpvar_chr = self.compute_snpvar_chr(args, chr_num, use_ridge=use_ridge)            
            df_snpvar_chr_list.append(df_snpvar_chr)
        df_snpvar = pd.concat(df_snpvar_chr_list, axis=0)
        df_snpvar.reset_index(inplace=True, drop=True)
        
        #save snpvar to a class member
        if use_ridge:
            self.df_snpvar_ridge = df_snpvar
        else:
            self.df_snpvar = df_snpvar
        
        
        
    def create_df_bins(self, bin_sizes, df_snpvar, df_snpvar_sorted=None, min_bin_size=10):
    
        #sort df_snpvar if needed
        if df_snpvar_sorted is None:
            df_snpvar_sorted = df_snpvar['SNPVAR'].sort_values()
        assert bin_sizes.sum() == df_snpvar_sorted.shape[0]
        
        
        #rearrange bins to prevent very small bins
        bin_i = len(bin_sizes)-1
        while True:
        
            #if the current bin is large enough, proceed to the previous one
            if bin_sizes[bin_i] >= min_bin_size:
                bin_i -= 1
                if bin_i==0: break
                continue
            
            #Compare the effects of the weakest bin in the current bin, and the strongest bin in the previous bin
            bin_start_ind = bin_sizes[:bin_i].sum()
            weakest_bin_snp = df_snpvar_sorted.iloc[::-1].iloc[bin_start_ind]
            strongest_lastbin_snp = df_snpvar_sorted.iloc[::-1].iloc[bin_start_ind-1]
                
            num_snps_to_transfer = np.minimum(min_bin_size-bin_sizes[bin_i], bin_sizes[bin_i-1])
            bin_sizes[bin_i] += num_snps_to_transfer
            bin_sizes[bin_i-1] -= num_snps_to_transfer
            
            #if we emptied the previous bin, delete it
            if bin_sizes[bin_i-1]==0:
                bin_sizes = np.concatenate((bin_sizes[:bin_i-1], bin_sizes[bin_i:]))
                bin_i -= 1
            
            #if the current bin is large enough, move to the previous one
            if bin_sizes[bin_i] >= min_bin_size:
                bin_i -= 1
                
            if bin_i==0: break
                         
        
        #create df_bins
        ind=0
        df_bins = pd.DataFrame(index=df_snpvar_sorted.index)        
        for bin_i, bin_size in enumerate(bin_sizes):
            snpvar_bin = np.zeros(df_bins.shape[0], dtype=np.bool)
            snpvar_bin[ind : ind+bin_size] = True
            df_bins['snpvar_bin%d'%(len(bin_sizes) - bin_i)] = snpvar_bin
            ind += bin_size
        assert np.all(df_bins.sum(axis=0) == bin_sizes)
        df_bins = df_bins.iloc[:, ::-1]
        assert df_bins.shape[0] == df_snpvar.shape[0]
        assert np.all(df_bins.sum(axis=1)==1)

        #reorder df_bins
        df_bins = df_bins.loc[df_snpvar.index]
        df_bins = pd.concat((df_snpvar[SNP_COLUMNS], df_bins), axis=1)
        assert np.all(df_bins.index == df_snpvar.index)
        
        return df_bins
            
        
        
    def partition_snps_Ckmedian(self, args, use_ridge):
        logging.info('Clustering SNPs into bins using the R Ckmeans.1d.dp package')
        
        #try loading the Ckmeans.1d.dp package
        try:
            import rpy2
            import rpy2.robjects.numpy2ri as numpy2ri
            try:
                from importlib import reload
                reload(rpy2.robjects.numpy2ri)
            except:
                pass
            import rpy2.robjects as ro
            ro.conversion.py2ri = numpy2ri
            numpy2ri.activate()
            from rpy2.robjects.packages import importr
            importr('Ckmeans.1d.dp')
            median_seg_func = ro.r('Ckmedian.1d.dp')
            mean_seg_func = ro.r('Ckmeans.1d.dp')            
        except:
            logging.error('Could not load the R package Ckmeans.1d.dp. Either install it or rerun PolyFun with --skip-Ckmedian')
            logging.error('')
            raise    

        #access the right class member
        if use_ridge:
            df_snpvar = self.df_snpvar_ridge
        else:
            df_snpvar = self.df_snpvar

        #sort df_snpvar
        df_snpvar_sorted = df_snpvar['SNPVAR'].sort_values()
        
        #perform the segmentation
        if args.num_bins is None or args.num_bins<=0:
            logging.info('Determining the optimal number of bins (if this is slow, consider using --num-bins 20 (or some other number))')
            seg_obj = median_seg_func(df_snpvar_sorted.values, k=np.array([5,30]))
        else:
            seg_obj = median_seg_func(df_snpvar_sorted.values, k=args.num_bins)
        bin_sizes = np.array(seg_obj.rx2('size')).astype(np.int)
        num_bins = len(bin_sizes)
        logging.info('Ckmedian.1d.dp partitioned SNPs into %d bins'%(num_bins))        

        #define df_bins
        df_bins = self.create_df_bins(bin_sizes, df_snpvar, df_snpvar_sorted=df_snpvar_sorted)
        return df_bins
        
        
        
        
    def partition_snps_Kmeans(self, args, use_ridge):
        logging.info('Clustering SNPs into bins using K-means clustering with %d bins'%(args.num_bins))

        #make sure that we can run K-means clustering
        assert args.num_bins is not None and args.num_bins>0
        try:
            from sklearn.cluster import KMeans
        except ImportError:          
            raise ImportError('sklearn not properly installed. Please reinstall it')
            
        #access the right class member
        if use_ridge: df_snpvar = self.df_snpvar_ridge            
        else: df_snpvar = self.df_snpvar
            
            
        #perform K-means clustering
        kmeans_obj = KMeans(n_clusters=args.num_bins)
        kmeans_obj.fit(df_snpvar[['SNPVAR']])
        assert kmeans_obj.cluster_centers_.shape[0] == args.num_bins

        #Make sure that clusters are contiguous
        bins_order = np.argsort(kmeans_obj.cluster_centers_[:,0])
        for bin_i, cluster_label in enumerate(bins_order[:-1]):
            next_cluster_label = bins_order[bin_i+1]
            assert df_snpvar.loc[kmeans_obj.labels_==cluster_label, 'SNPVAR'].max() <= df_snpvar.loc[kmeans_obj.labels_==next_cluster_label, 'SNPVAR'].min()
        
        #define bin_sizes
        bin_sizes = np.bincount(kmeans_obj.labels_)[bins_order]
            
        #define df_bins
        df_bins = self.create_df_bins(bin_sizes, df_snpvar, df_snpvar_sorted=None)
        return df_bins
        

    def partition_snps_to_bins(self, args, use_ridge):
    
        #if skip_ckmedian was specified, run regular K-means
        if args.skip_Ckmedian:
            self.df_bins = self.partition_snps_Kmeans(args, use_ridge=use_ridge)
        else:
            self.df_bins = self.partition_snps_Ckmedian(args, use_ridge=use_ridge)
        
        
    def save_bins_to_disk(self, args):
        logging.info('Saving SNP-bins to disk')
        for chr_num in tqdm(range(1,23)):

            #save bins file to disk
            df_bins_chr = self.df_bins.query('CHR==%d'%(chr_num))
            bins_chr_file = get_file_name(args, 'bins', chr_num, verify_exists=False)
            df_bins_chr.to_parquet(bins_chr_file, index=False)
            
            #save M files to disk
            M_chr_file = get_file_name(args, 'M', chr_num, verify_exists=False)
            M_chr = df_bins_chr.drop(columns=SNP_COLUMNS).sum(axis=0).values
            np.savetxt(M_chr_file, M_chr.reshape((1, M_chr.shape[0])), fmt='%i')
            
            
        
    
    def save_snpvar_to_disk(self, args, use_ridge, constrain_range):
        if constrain_range:
            logging.info('Saving constrained SNP variances to disk')
        else:
            logging.info('Saving SNP variances to disk')
        
        #determine which df_snpvar to use
        if use_ridge: df_snpvar = self.df_snpvar_ridge            
        else: df_snpvar = self.df_snpvar

        #constrain the ratio between the largest and smallest snp-var
        if constrain_range:
            df_snpvar = df_snpvar.copy()
            h2_total = df_snpvar['SNPVAR'].sum()
            min_snpvar = df_snpvar['SNPVAR'].max() / args.q
            df_snpvar.loc[df_snpvar['SNPVAR'] < min_snpvar, 'SNPVAR'] = min_snpvar
            df_snpvar['SNPVAR'] *= h2_total / df_snpvar['SNPVAR'].sum()
            assert np.isclose(df_snpvar['SNPVAR'].sum(), h2_total)
            
        #merge snpvar with sumstats
        try:
            df_sumstats = pd.read_parquet(args.sumstats)
        except (ArrowIOError, ArrowInvalid):
            df_sumstats = pd.read_table(args.sumstats, sep='\s+')
        df_sumstats.drop(columns=['SNP'], errors='ignore', inplace=True)
        for col in ['CHR', 'BP', 'A1', 'A2']:
            if col not in df_sumstats.columns:
                raise ValueError('sumstats file has a missing column: %s'%(col))
        df_snpvar = set_snpid_index(df_snpvar, copy=True)
        df_sumstats = set_snpid_index(df_sumstats)
        svpvar_cols = df_snpvar.columns.copy()
        df_snpvar.drop(columns=['CHR', 'BP', 'A1', 'A2'], inplace=True)
        df_snpvar = df_snpvar.merge(df_sumstats, left_index=True, right_index=True)
        df_snpvar = df_snpvar[list(svpvar_cols) + [c for c in df_sumstats.columns if c not in list(svpvar_cols)]]
        if df_snpvar.shape[0] < df_sumstats.shape[0]:
            error_message = 'not all SNPs in the sumstats file and/or in the LD reference files are also in the annotations file'
            if args.allow_missing:
                logging.warning(error_message + '. Keeping %d/%d SNPs'%(df_snpvar.shape[0], df_sumstats.shape[0]))
            else:
                raise ValueError(error_message + '. If you wish to omit the missing SNPs, please use the flag --allow-missing')

        #iterate over chromosomes 
        for chr_num in tqdm(range(1,23)):
        
            #define output file name
            output_fname = 'snpvar'
            if use_ridge: output_fname += '_ridge'
            if constrain_range: output_fname += '_constrained'
            snpvar_chr_file = get_file_name(args, output_fname, chr_num, verify_exists=False)
                
            #save snpvar to file
            df_snpvar_chr = df_snpvar.query('CHR==%d'%(chr_num))
            df_snpvar_chr.to_csv(snpvar_chr_file, index=False, sep='\t', compression='gzip', float_format='%0.4e')
        
        
        
    def polyfun_h2_L2(self, args):
        #run Ridge regression
        self.run_ldsc(args, use_ridge=True, nn=False, evenodd_split=False, keep_large=False)

        #compute per-SNP h^2 based on L2-regularized S-LDSC coefficients
        self.compute_snpvar(args, use_ridge=True)
        
        #save L2-regularized S-LDSC per-SNP h^2 to disk
        self.save_snpvar_to_disk(args, use_ridge=True, constrain_range=True)
        self.save_snpvar_to_disk(args, use_ridge=True, constrain_range=False)

        #partitions SNPs into bins and save to disk, unless explicitly requested not to
        if not args.no_partitions:
            self.partition_snps_to_bins(args, use_ridge=True)
            self.save_bins_to_disk(args)
    
        
    def load_bins_chr(self, args, chr_num):
        bins_file = get_file_name(args, 'bins', chr_num)
        df_bins_chr = pd.read_parquet(bins_file)
        return df_bins_chr
        
        
    def compute_ld_scores(self, args):    
        #define the range of chromosomes to iterate over
        if args.chr is None:
            chr_range = range(1,23)
        else:
            chr_range = range(args.chr, args.chr+1)
        
        #iterate over chromosomes and compute LD-scores
        ###df_ldscores_chr_list = []
        for chr_num in tqdm(chr_range, disable=len(chr_range)==1):
        
            #load or extract the bins for the current chromosome
            try:
                df_bins_chr = self.df_bins.query('CHR==%d'%(chr_num))
            except AttributeError:
                df_bins_chr = self.load_bins_chr(args, chr_num)
                
            #compute LD-scores for this chromosome
            if args.ld_ukb:
                if args.ld_dir is None: ld_dir = tempfile.mkdtemp()
                else: ld_dir = args.ld_dir
                df_bins_chr = set_snpid_index(df_bins_chr)
                df_ldscores_chr = compute_ldscores_chr(df_bins_chr, ld_dir=ld_dir, use_ukb=True)
            elif args.bfile_chr is not None:
                df_ldscores_chr = self.compute_ldscores_plink_chr(args, chr_num, df_bins_chr)
            else:
                raise ValueError('no LDscore computation method specified')
                
            #save the LD-scores to disk
            ldscores_output_file = get_file_name(args, 'ldscores', chr_num, verify_exists=False)
            df_ldscores_chr.to_parquet(ldscores_output_file, index=False)
            
            # #add the ldscores to the LDscores list
            # df_ldscores_chr_list.append(df_ldscores_chr)
        
        # #concatenate all the LD-score dfs
        # if len(df_ldscores_chr_list)==1:
            # self.df_bin_ldscores = df_ldscores_chr_list[0]
        # else:
            # self.df_bin_ldscores = pd.concat(df_ldscores_chr_list, axis=0)
            # self.df_bin_ldscores.reset_index(inplace=True, drop=True)
        
                
                    
                

        
    def compute_ldscores_plink_chr(self, args, chr_num, df_bins_chr):
    
        # read bim/snp
        bim_file = get_file_name(args, 'bim', chr_num)
        array_snps = parse.PlinkBIMFile(bim_file)
        df_bim = array_snps.df
        df_bim = set_snpid_index(df_bim)
        
        #Remove annotations of SNPs that are not in the .bim file
        df_bins_chr = set_snpid_index(df_bins_chr)
        df_bins_chr = df_bins_chr.loc[df_bins_chr.index.isin(df_bim.index)]

        #make sure that all SNPs have a bin
        keep_snps = None
        if np.any(~df_bim.index.isin(df_bins_chr.index)):
            error_msg = 'Not all SNPs were assigned a bin (meaning some SNPS in the summary statistics and/or in the LD reference files are not in the annotation files)'
            if args.allow_missing:
                is_good_snp = df_bim.index.isin(df_bins_chr.index)
                if not np.any(is_good_snp):
                    raise ValueError('No SNPs in chromosome %d have annotations'%(chr_num))
                keep_snps = np.where(is_good_snp)[0]
                logging.warning(error_msg)
                logging.warning('Keeping only %d/%d SNPs in chromosome %d that have annotations'%(df_bim.shape[0], len(is_good_snp), chr_num))
            else:
                raise ValueError(error_msg + '. If you wish to omit the missing SNPs, please use the flag --allow-missing')

        #find #individuals in bfile
        fam_file = get_file_name(args, 'fam', chr_num)
        df_fam = pd.read_table(fam_file, header=None, usecols=[5], sep='\s+')
        n = df_fam.shape[0]
        
        #find keep_indivs
        if args.keep is None:
            keep_indivs= None
        else:
            array_indivs = parse.PlinkFAMFile(args.bfile+'.fam')
            keep_indivs = __filter__(args.keep, 'individuals', 'include', array_indivs)
            logging.info('after applying --keep, %d individuals remain'%(len(keep_indivs)))

        #read plink file
        logging.info('Loading SNP file...')
        bed_file = get_file_name(args, 'bed', chr_num)
        geno_array = ldscore.PlinkBEDFile(bed_file, n, array_snps, keep_snps=keep_snps,
            keep_indivs=keep_indivs, mafMin=None)
        
        #remove omitted SNPs from df_bim
        if len(geno_array.kept_snps) != df_bim.shape[0]:
            assert np.all(np.array(geno_array.kept_snps) == np.sort(np.array(geno_array.kept_snps)))
            assert geno_array.kept_snps[-1] < df_bim.shape[0]
            df_bim = df_bim.iloc[geno_array.kept_snps]
            
        #rearrange annotations to match the order of SNPs in the plink file
        assert df_bins_chr.shape[0] >= df_bim.shape[0]
        if (df_bins_chr.shape[0] > df_bim.shape[0]) or np.any(df_bins_chr.index != df_bim.index):
            assert np.all(df_bim.index.isin(df_bins_chr.index))
            df_bins_chr = df_bins_chr.loc[df_bim.index]
            
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
        logging.info('Computing LD scores for chromosome %d'%(chr_num))
        ldscores = geno_array.ldScoreVarBlocks(block_left, args.chunk_size, annot=df_bins_chr.drop(columns=SNP_COLUMNS).values)
        
        #create an ldscores df
        df_ldscores = pd.DataFrame(ldscores, index=df_bins_chr.index, columns=df_bins_chr.drop(columns=SNP_COLUMNS).columns)
        df_ldscores = pd.concat((df_bins_chr[SNP_COLUMNS], df_ldscores), axis=1)
        return df_ldscores
        

    
    def compute_h2_bins(self, args, constrain_range):
        #run S-LDSC 
        self.run_ldsc(args, use_ridge=False, nn=True, evenodd_split=True, keep_large=False)

        #compute per-SNP h^2
        self.compute_snpvar(args, use_ridge=False)
        
        #save per-SNP h^2 to disk
        self.save_snpvar_to_disk(args, use_ridge=False, constrain_range=constrain_range)

    
        

                
    def polyfun_main(self, args):
    
        #compute snp variances using L2-regularized S-LDSC with an odd/even chromosome split
        if args.compute_h2_L2:
            self.polyfun_h2_L2(args)
            
        #compute LD-scores of SNP partitions
        if args.compute_ldscores:
            self.compute_ld_scores(args)
        
        if args.compute_h2_bins:
            self.compute_h2_bins(args, constrain_range=True)
        
    


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()

    #partitioning-related parameters
    parser.add_argument('--num-bins', type=int, default=None, help='Number of bins to partition SNPs into. If not specified, PolyFun will automatically select this number based on a BIC criterion')
    parser.add_argument('--anno', default=None, help='Comma-delimited list of annotations to use (if not specified, will use all annotations)')
    parser.add_argument('--skip-Ckmedian', default=False, action='store_true', help='If specified, use a regular K-means algorithm instead of the R Ckmeans.1d.dp package')
    
    #mode related parameters
    parser.add_argument('--compute-ldscores', default=False, action='store_true', help='If specified, PolyFun will compute LD-scores of SNP bins')
    parser.add_argument('--compute-h2-L2', default=False, action='store_true', help='If specified, PolyFun will compute per-SNP h2 using L2-regularized S-LDSC')
    parser.add_argument('--compute-h2-bins', default=False, action='store_true', help='If specified, PolyFun will robustly compute per-SNP h2 based on SNP bins')
    parser.add_argument('--no-partitions', default=False, action='store_true', help='If specified, PolyFun will not partition SNPs into bins (do this only if you want to use the per-SNP h^2 from L2-regularized S-LDSC as prior causal probabilities, which is typically not recommended)')
    
    #ld-score related parameters
    parser.add_argument('--chr', type=int, default=None, help='Chromosome number (only applicable when only specifying --ldscores). If not set, PolyFun will compute LD-scores for all chromosomes')
    parser.add_argument('--ld-wind-cm', type=float, default=None, help='window size to be used for estimating LD-scores in units of centiMorgans (cM).')
    parser.add_argument('--ld-wind-kb', type=int, default=None, help='window size to be used for estimating LD-scores in units of Kb.')
    parser.add_argument('--ld-wind-snps', type=int, default=None, help='window size to be used for estimating LD-scores in units of SNPs.')
    parser.add_argument('--chunk-size',  type=int, default=50, help='chunk size for LD-scores calculation')
    parser.add_argument('--keep',  default=None, help='File with ids of individuals to use when computing LD-scores')
    
    #per-SNP h2 related parameters
    parser.add_argument('--q', type=float, default=100, help='The maximum ratio between the largest and smallest truncated per-SNP heritabilites')

    #data input/output parameters
    parser.add_argument('--sumstats', help='Input summary statistics file')
    parser.add_argument('--ref-ld-chr', help='Suffix of LD-score files (as in ldsc)')
    parser.add_argument('--w-ld-chr', help='Suffix of LD-score weights files (as in ldsc)')
    parser.add_argument('--bfile-chr', default=None, help='Prefix of plink files (used to compute LD-scores)')
    parser.add_argument('--ld-ukb', default=False, action='store_true', help='If specified, PolyFun will use UKB LD matrices to compute LD-scores')
    parser.add_argument('--ld-dir', default=None, help='The path of a directory with UKB LD files (if not specified PolyFun will create a temporary directory)')
    parser.add_argument('--output-prefix', required=True, help='Prefix of all PolyFun output file names')    
    parser.add_argument('--allow-missing', default=False, action='store_true', help='If specified, PolyFun will not terminate if some SNPs with sumstats are not found in the annotations files')
    
    #LDSC parameters
    parser.add_argument('--nnls-exact', default=False, action='store_true', help='If specified, S-LDSC will estimate non-negative taus using an exact instead of an approximate solver (this will be slower but slightly more accurate)')
    
    
    #check package versions
    check_package_versions()
    
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
    
    #create and run PolyFun object
    polyfun_obj = PolyFun()
    polyfun_obj.polyfun_main(args)
    
    print()
    
