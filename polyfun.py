import numpy as np
import pandas as pd
import os
import sys
import time
from ldsc_polyfun import jackknife, regressions, sumstats, ldscore, parse
import logging
from copy import deepcopy
from tqdm import tqdm

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
    print('* (C) 2019 Omer Weissbrod')
    print('*********************************************************************')
    print()
    
    
def check_args(args):    
    mode_params = np.array([args.compute_h2_L2, args.compute_ldscores, args.compute_h2_bins])
    
    #verify that the requested computations are valid
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
    if args.compute_ldscores:
        if args.bfile_chr is None:
            raise ValueError('You must specify --bfile-chr when you specify --compute-ldscores')    
        if args.ld_wind_cm is None and args.ld_wind_kb is None and args.ld_wind_snps is None:
            args.ld_wind_cm = 1.0
            logging.warning('no ld-wind argument specified.  PolyFun will use --ld-cm 1.0')
            
    
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
            
            
    return args

def check_files(args):
        
    #check that required input files exist
    if args.compute_h2_L2:
        if not os.path.exists(args.sumstats):
            raise IOError('Cannot find sumstats file %s'%(args.sumstats))
        for chr_num in range(1,23):
            get_file_name(args, 'ref-ld', chr_num, verify_exists=True)
            get_file_name(args, 'w-ld', chr_num, verify_exists=True)
            get_file_name(args, 'annot', chr_num, verify_exists=True)
            
    if args.compute_ldscores:
        if args.chr is None: chr_range = range(1,23)            
        else: chr_range = range(args.chr, args.chr+1)
        
        if args.bfile_chr is not None:
            for chr_num in chr_range:
                get_file_name(args, 'bim', chr_num, verify_exists=True)
                get_file_name(args, 'fam', chr_num, verify_exists=True)
                get_file_name(args, 'bed', chr_num, verify_exists=True)
        
        if not args.compute_h2_L2:
            for chr_num in chr_range:
                get_file_name(args, 'snpvar_ridge', chr_num, verify_exists=True)
                
    if args.compute_h2_bins and not args.compute_ldscores:
        for chr_num in range(1,23):
            get_file_name(args, 'bins', chr_num, verify_exists=True)
            
        

    

''' Logger class (for compatability with LDSC code)'''
class Logger(object):
    def __init__(self):
        pass
    def log(self, msg):
        logging.info(msg)
        
        
class TqdmHandler(logging.StreamHandler):
    def __init__(self):
        logging.StreamHandler.__init__(self)

    def emit(self, record):
        msg = self.format(record)
        tqdm.write(msg)


def configure_logger(out_prefix):

    logFormatter = logging.Formatter("[%(levelname)s]  %(message)s")
    logger = logging.getLogger()
    logger.setLevel(logging.NOTSET)
    
    consoleHandler = TqdmHandler()
    consoleHandler.setFormatter(logFormatter)
    logger.addHandler(consoleHandler)
    
    fileHandler = logging.FileHandler(out_prefix+'.log')
    fileHandler.setFormatter(logFormatter)
    logger.addHandler(fileHandler)
    
        
def get_file_name(args, file_type, chr_num, verify_exists=True):
    if file_type == 'ldscores':
        file_name = args.output_prefix + '.%d.l2.ldscore.parquet'%(chr_num)
    elif file_type == 'snpvar_ridge':
        file_name = args.output_prefix + '.%d.snpvar_ridge.gz'%(chr_num)
    elif file_type == 'snpvar_ridge_constrained':
        file_name = args.output_prefix + '.%d.snpvar_ridge_constrained.gz'%(chr_num)        
    elif file_type == 'snpvar_constrained':
        file_name = args.output_prefix + '.%d.snpvar_constrained.gz'%(chr_num)        
    elif file_type == 'bins':
        file_name = args.output_prefix + '.%d.bins.parquet'%(chr_num)
    elif file_type == 'M':
        file_name = args.output_prefix + '.%d.l2.M'%(chr_num)
        
    elif file_type == 'annot':
        assert verify_exists
        file_name = args.ref_ld_chr + '%d.annot.gz'%(chr_num)
        if not os.path.exists(file_name):
            file_name = args.ref_ld_chr + '%d.annot.parquet'%(chr_num)
        
    elif file_type == 'ref-ld':
        assert verify_exists
        file_name = args.ref_ld_chr + '%d.l2.ldscore.gz'%(chr_num)
        if not os.path.exists(file_name):
            file_name = args.ref_ld_chr + '%d.l2.ldscore.parquet'%(chr_num)
        
    elif file_type == 'w-ld':
        assert verify_exists
        file_name = args.w_ld_chr + '%d.l2.ldscore.gz'%(chr_num)
        if not os.path.exists(file_name):
            file_name = args.w_ld_chr + '%d.l2.ldscore.parquet'%(chr_num)
    
    
    elif file_type == 'bim':
        file_name = args.bfile_chr + '%d.bim'%(chr_num)
    elif file_type == 'fam':
        file_name = args.bfile_chr + '%d.fam'%(chr_num)
    elif file_type == 'bed':
        file_name = args.bfile_chr + '%d.bed'%(chr_num)
    else:
        raise ValueError('unknown file type')
        
    if verify_exists:
        if not os.path.exists(file_name):
            raise IOError('%s file not found: %s'%(file_type, file_name))
            
    return file_name
    
    
    
class PolyFun:
    def __init__(self):
        pass
        
    def run_ldsc(self, args, use_ridge):

        #prepare LDSC objects
        log = Logger()
        args.h2 = args.sumstats
        args.ref_ld = None
        args.w_ld = None
        args.n_blocks = 2
        args.M = None
        args.not_M_5_50 = True
        
        #if not ridge, the LD-scores are our bins
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
            df_sumstats = pd.read_parquet(args.sumstats)
            ###merge everything together...
            
        #prepare LD-scores for S-LDSC run
        ref_ld = np.array(df_sumstats[ref_ld_cnames], dtype=np.float32)
        sumstats._check_ld_condnum(args, log, ref_ld_cnames)
        if df_sumstats.shape[0] < 200000:
            logging.warning('number of SNPs less than 200k; this is almost always bad.')
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
        self.ref_ld_cnames = ref_ld_cnames
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
            evenodd_split=(not use_ridge),
            nn=(not use_ridge),
            keep_large=(not use_ridge)
            )
            
        #save the results object
        if use_ridge:
            self.hsqhat_ridge = hsqhat
        else:
            self.hsqhat = hsqhat
            

    
    def load_annotations_file(self, args, chr_num, use_ridge):
        
        #load annotations file for this chromosome
        if use_ridge:
            annot_filename = get_file_name(args, 'annot', chr_num)
        else:
            annot_filename = get_file_name(args, 'bins', chr_num)
        
        #load parquet file
        if annot_filename.endswith('.parquet'):
            df_annot_chr = pd.read_parquet(annot_filename)
        else:
            df_annot_chr = pd.read_table(annot_filename)
        
        #make sure all required columns were found
        df_annot_chr.drop(columns=['CM'], inplace=True, errors='ignore')
        if not ('A1' in df_annot_chr.index.names):
            found_missing_col = False
            for colname in ['CHR', 'SNP', 'BP', 'A1', 'A2']:
                if colname not in df_annot_chr.columns:
                    logging.error('%s has a missing column: %s'%(annot_filename, colname))
                    found_missing_col = True
            if found_missing_col:
                raise ValueError('Missing columns found in %s'%(annot_filename))
            df_annot_chr.set_index(['CHR', 'BP', 'SNP', 'A1', 'A2'], inplace=True, drop=True)
            
        #subset annotations if requested
        if args.anno is not None:
            anno_to_use = args.anno.split(',')
            assert np.all(np.isin(anno_to_use, df_annot_chr.columns))
            df_annot_chr = df_annot_chr[anno_to_use]
            
        #if we have more annotations that ref-ld, it might mean that some annotations were removed, so remove them from here as well
        if not np.all(self.ref_ld_cnames.str[:-2].isin(df_annot_chr.columns)):
            raise ValueError('Annotation names in annotations file do not match the one in the LD-scores file')
        if len(self.ref_ld_cnames) < len(df_annot_chr.columns):
            assert np.all(np.isin(self.ref_ld_cnames.str[:-2], df_annot_chr.columns))
            df_annot_chr = df_annot_chr[self.ref_ld_cnames.str[:-2]]

        #make sure that we get the same columns as the ones in the LD-score files
        if not np.all(df_annot_chr.columns == self.ref_ld_cnames.str[:-2]):
            raise ValueError('Annotation names in annotations file do not match the one in the LD-scores file')            

        return df_annot_chr


    def compute_snpvar_chr(self, args, chr_num, use_ridge):        

        #load annotations file from disk
        df_annot_chr = self.load_annotations_file(args, chr_num, use_ridge)
            
        #extract taus from a jknife object
        if use_ridge:
            hsqhat = self.hsqhat_ridge
            jknife = hsqhat.jknife_ridge

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
            taus = jknife.est_loco[chr_num-1][:hsqhat.n_annot] / hsqhat.Nbar
            
        #compute and return the snp variances
        df_snpvar_chr = df_annot_chr.dot(taus)
        return df_snpvar_chr
        

    def compute_snpvar(self, args, use_ridge):
        logging.info('Computing per-SNP h^2 for each chromosome...')
        
        #iterate over chromosomes
        df_snpvar_chr_list = []
        for chr_num in tqdm(range(1,23)):
            df_snpvar_chr = self.compute_snpvar_chr(args, chr_num, use_ridge=use_ridge)
            df_snpvar_chr = df_snpvar_chr.to_frame(name='snpvar')
            df_snpvar_chr_list.append(df_snpvar_chr)
        df_snpvar = pd.concat(df_snpvar_chr_list, axis=0)
        
        #save snpvar to a class member
        if use_ridge:
            self.df_snpvar_ridge = df_snpvar
        else:
            self.df_snpvar = df_snpvar
        
        
        
    def partition_snpvar_Ckmedian(self, args, use_ridge):
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
        df_snpvar_sorted = df_snpvar['snpvar'].sort_values()
        df_bins = pd.DataFrame(index=df_snpvar_sorted.index)
        
        #perform the segmentation
        if args.num_bins is None or args.num_bins<=0:
            logging.info('Determining the optimal number of bins (if this is slow, consider using --num-bins 20 (or some other number))')
            seg_obj = median_seg_func(df_snpvar_sorted.values, k=np.array([5,30]))
        else:
            seg_obj = median_seg_func(df_snpvar_sorted.values, k=args.num_bins)
        sizes = np.array(seg_obj.rx2('size')).astype(np.int)
        assert sizes.sum() == df_snpvar_sorted.shape[0]

        #create df_bins
        num_bins = len(sizes)
        logging.info('Ckmedian.1d.dp partitioned SNPs into %d bins'%(num_bins))
        ind=0
        for s_i, s in enumerate(sizes):
            snpvar_bin = np.zeros(df_bins.shape[0], dtype=np.bool)
            snpvar_bin[ind : ind+s] = True
            df_bins['snpvar_bin%d'%(num_bins-s_i)] = snpvar_bin
            ind += s
        assert df_bins.shape[0] == df_snpvar.shape[0]
        assert np.all(df_bins.sum(axis=1)==1)

        #reorder df_bins (multi-index based sorting is so slow...)
        df_bins.reset_index(inplace=True)
        df_bins.index = df_bins['CHR'].astype(str) + '.' \
                      + df_bins['BP'].astype(str) + '.' \
                      + df_bins['A1'] + '.' \
                      + df_bins['A2']
        snpvar_index = df_snpvar.index.get_level_values('CHR').astype(str) + '.' \
                      + df_snpvar.index.get_level_values('BP').astype(str) + '.' \
                      + df_snpvar.index.get_level_values('A1') + '.' \
                      + df_snpvar.index.get_level_values('A2')
        df_bins = df_bins.loc[snpvar_index]
        df_bins.set_index(['CHR', 'BP', 'SNP', 'A1', 'A2'], inplace=True)
        assert np.all(df_bins.index == df_snpvar.index)
        
        return df_bins
        
        
        
    def partition_snpvar_Kmeans(self, args, use_ridge):
        logging.info('Clustering SNPs into bins using K-means clustering with %d bins'%(args.num_bins))

        #make sure that we can run K-means clustering
        assert args.num_bins is not None and args.num_bins>0
        try:
            from sklearn.cluster import KMeans
        except ImportError:          
            raise ImportError('sklearn not properly installed. Please reinstall it')
            
        #access the right class member
        if use_ridge:
            df_snpvar = self.df_snpvar_ridge
        else:
            df_snpvar = self.df_snpvar
            
        #perform K-means clustering
        kmeans_obj = KMeans(n_clusters=args.num_bins)
        kmeans_obj.fit(df_snpvar[['snpvar']])
        assert kmeans_obj.cluster_centers_.shape[0] == args.num_bins

        #Create df_bins
        bins_order = np.argsort(kmeans_obj.cluster_centers_[:,0])[::-1]
        df_bins = pd.DataFrame(index=df_snpvar.index)
        for bin_i, cluster_label in enumerate(bins_order):        
            df_bins['snpvar_bin%d'%(bin_i+1)] = (kmeans_obj.labels_==cluster_label)
            assert df_bins['snpvar_bin%d'%(bin_i+1)].any()
        assert np.all(df_bins.sum(axis=1)==1)
        
        return df_bins
        


    def partition_snps_to_bins(self, args, use_ridge):
    
        #if skip_ckmedian was specified, run regular K-means
        if args.skip_Ckmedian:
            self.df_bins = self.partition_snpvar_Kmeans(args, use_ridge=use_ridge)
        else:
            self.df_bins = self.partition_snpvar_Ckmedian(args, use_ridge=use_ridge)
        
        
    def save_bins_to_disk(self, args):
        logging.info('Saving SNP-bins to disk')
        for chr_num in tqdm(range(1,23)):

            #save bins file to disk
            df_bins_chr = self.df_bins.query('CHR==%d'%(chr_num))            
            bins_chr_file = get_file_name(args, 'bins', chr_num, verify_exists=False)
            df_bins_chr.reset_index().to_parquet(bins_chr_file, index=False)
            
            #save M file to disk
            M_chr_file = get_file_name(args, 'M', chr_num, verify_exists=False)
            M_chr = df_bins_chr.sum(axis=0).values
            np.savetxt(M_chr_file, M_chr.reshape((1, M_chr.shape[0])), fmt='%i')
            
            
        
    
    def save_snpvar_to_disk(self, args, use_ridge, constrain_range):
        if constrain_range:
            logging.info('Saving constrained SNP variances to disk')
        else:
            logging.info('Saving SNP variances to disk')
        
        #determine which df_snpvar to use
        if use_ridge:
            df_snpvar = self.df_snpvar_ridge
        else:
            df_snpvar = self.df_snpvar
            assert constrain_range

        #constrain the ratio between the largest and smallest snp-var
        if constrain_range:
            df_snpvar = df_snpvar.copy()
            h2_total = df_snpvar['snpvar'].sum()
            min_snpvar = df_snpvar['snpvar'].max() / args.q
            df_snpvar.loc[df_snpvar['snpvar'] < min_snpvar, 'snpvar'] = min_snpvar
            df_snpvar['snpvar'] *= h2_total / df_snpvar['snpvar'].sum()
            assert np.isclose(df_snpvar['snpvar'].sum(), h2_total)
            
        #merge snpvar with sumstats
        df_sumstats = pd.read_parquet(args.sumstats)
        for col in ['SNP', 'A1', 'A2']:
            if col not in df_sumstats.columns:
                raise ValueError('sumstats file has a missing column: %s'%(col))
        sumstats_index = df_sumstats['SNP'] + df_sumstats['A1'] + df_sumstats['A2']
        df_snpvar_index = df_snpvar.index.get_level_values('SNP') + df_snpvar.index.get_level_values('A1') + df_snpvar.index.get_level_values('A2')
        if np.any(~sumstats_index.isin(df_snpvar_index)):
            raise ValueError('not all SNPs in the sumstats file are also in the annotations file')
        df_snpvar = df_snpvar.reset_index()
        df_snpvar = df_snpvar.merge(df_sumstats, on=['SNP', 'A1', 'A2'])
        assert df_snpvar.shape[0] == df_sumstats.shape[0]

        #iterate over chromosomes 
        for chr_num in tqdm(range(1,23)):
        
            #define output file name
            if use_ridge:
                if constrain_range:
                    snpvar_chr_file = get_file_name(args, 'snpvar_ridge_constrained', chr_num, verify_exists=False)
                else:
                    snpvar_chr_file = get_file_name(args, 'snpvar_ridge', chr_num, verify_exists=False)
            else:
                snpvar_chr_file = get_file_name(args, 'snpvar_constrained', chr_num, verify_exists=False)
                
            #save snpvar to file
            df_snpvar_chr = df_snpvar.query('CHR==%d'%(chr_num))
            df_snpvar_chr.to_csv(snpvar_chr_file, index=False, sep='\t', compression='gzip', float_format='%0.4e')
        
        
        
    def polyfun_h2_L2(self, args):
        #run Ridge regression
        self.run_ldsc(args, use_ridge=True)

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
        df_ldscores_chr_list = []
        for chr_num in tqdm(chr_range, disable=len(chr_range)==1):
        
            #load or extract the bins for the current chromosome
            if args.compute_h2_L2:
                df_bins_chr = self.df_bins.reset_index().query('CHR==%d'%(chr_num))
            else:
                df_bins_chr = self.load_bins_chr(args, chr_num)
                
            #compute LD-scores for this chromosome
            if args.bfile_chr is not None:
                df_ldscores_chr = self.compute_ldscores_plink_chr(args, chr_num, df_bins_chr)
            # elif args.npz_prefix is not None:
                # raise NotImplementedError('--npz-prefix is not yet supported')
            else:
                raise ValueError('no LDscore computation method specified')
                
            #save the LD-scores to disk
            ldscores_output_file = get_file_name(args, 'ldscores', chr_num, verify_exists=False)
            df_ldscores_chr.reset_index().to_parquet(ldscores_output_file)
            
            #add the ldscores to the LDscores list
            df_ldscores_chr_list.append(df_ldscores_chr)
        
        #concatenate all the LD-score dfs
        if len(df_ldscores_chr_list)==1:
            self.df_bin_ldscores = df_ldscores_chr_list[0]
        else:
            self.df_bin_ldscores = pd.concat(df_ldscores_chr_list, axis=0)
        
                
                    
                

        
    def compute_ldscores_plink_chr(self, args, chr_num, df_bins_chr):
    
        # read bim/snp
        bim_file = get_file_name(args, 'bim', chr_num)
        array_snps = parse.PlinkBIMFile(bim_file)
        df_bim = array_snps.df
        
        #heuristically reduce df_bins_chr to a small superset of the relevant SNPs        
        df_bins_chr = df_bins_chr.loc[df_bins_chr['SNP'].isin(df_bim['SNP'])]
        
        #duplicate df_bins_chr to make sure we have no flipped alleles that will cause a mess
        bins_index1 =     df_bins_chr['CHR'].astype(str) + '.' \
                        + df_bins_chr['BP'].astype(str) + '.' \
                        + df_bins_chr['A1'] + '.' \
                        + df_bins_chr['A2']
        df_bins_chr2 = df_bins_chr.copy()
        df_bins_chr2['A2'] = df_bins_chr['A1'].copy()
        df_bins_chr2['A1'] = df_bins_chr['A2'].copy()
        df_bins_chr.index  = df_bins_chr['CHR'].astype(str) + '.' + df_bins_chr['BP'].astype(str) + '.' + df_bins_chr['A1'] + '.' + df_bins_chr['A2']
        df_bins_chr2.index = df_bins_chr2['CHR'].astype(str) + '.' + df_bins_chr2['BP'].astype(str) + '.' + df_bins_chr2['A1'] + '.' + df_bins_chr2['A2']        
        df_bins_chr = pd.concat([df_bins_chr, df_bins_chr2], axis=0)
        
        #make sure that all SNPs have a bin
        df_bim.index = df_bim['CHR'].astype(str) + '.' + df_bim['BP'].astype(str) + '.' + df_bim['A1'] + '.' + df_bim['A2']
        if np.any(~df_bim.index.isin(df_bins_chr.index)):
            raise ValueError('Not all SNPs were assigned a bin (meaning some SNPS are not in the annotation files)')
            
        #rearrange df_bins_chr to match the order of SNPs in the plink file
        if (df_bins_chr.shape[0] > df_bim.shape[0]) or np.any(df_bins_chr.index != df_bim.index):
            assert np.all(df_bim.index.isin(df_bins_chr.index))
            df_bins_chr = df_bins_chr.loc[df_bim.index]
        assert np.all(df_bins_chr['BP'] == df_bim['BP'].values)
        assert np.all(df_bins_chr['A1'] == df_bim['A1'].values)
        assert np.all(df_bins_chr['A2'] == df_bim['A2'].values)
        df_bins_chr.set_index(['CHR', 'BP', 'SNP', 'A1', 'A2'], drop=True, inplace=True)

        #find #individuals in bfile
        fam_file = get_file_name(args, 'fam', chr_num)
        df_fam = pd.read_table(fam_file, header=None, usecols=[5], delim_whitespace=True)
        n = df_fam.shape[0]
        
        #read plink file    
        logging.info('Loading SNP file...')
        bed_file = get_file_name(args, 'bed', chr_num)
        geno_array = ldscore.PlinkBEDFile(bed_file, n, array_snps, keep_snps=None,
            keep_indivs=None, mafMin=None)

            
        # determine block widths
        num_wind_args = np.array((args.ld_wind_snps, args.ld_wind_kb, args.ld_wind_cm), dtype=bool)
        if np.sum(num_wind_args) != 1:
            raise ValueError('Must specify exactly one --ld-wind option')
        if args.ld_wind_snps:
            max_dist = args.ld_wind_snps
            coords = np.array(list(range(geno_array.m)))
        elif args.ld_wind_kb:
            max_dist = args.ld_wind_kb*1000
            coords = np.array(df_bim['BP'])[geno_array.kept_snps]
            if len(np.unique(coords)) == 1:
                raise ValueError('bim file has no basepair data --- please use a different ld-wind option')
        elif args.ld_wind_cm:
            max_dist = args.ld_wind_cm
            coords = np.array(df_bim['CM'])[geno_array.kept_snps]
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
        ldscores = geno_array.ldScoreVarBlocks(block_left, args.chunk_size, annot=df_bins_chr.values)
        
        #create an ldscores df
        df_ldscores = pd.DataFrame(ldscores, index=df_bins_chr.index, columns=df_bins_chr.columns)
        return df_ldscores
        

    
    def compute_h2_bins(self, args):
        #run S-LDSC 
        self.run_ldsc(args, use_ridge=False)

        #compute per-SNP h^2 based on L2-regularized S-LDSC coefficients
        self.compute_snpvar(args, use_ridge=False)
        
        #save L2-regularized S-LDSC per-SNP h^2 to disk
        self.save_snpvar_to_disk(args, use_ridge=False, constrain_range=True)

    
        

                
    def polyfun_main(self, args):
    
        #compute snp variances using L2-regularized S-LDSC with an odd/even chromosome split
        if args.compute_h2_L2:
            self.polyfun_h2_L2(args)
            
        #compute LD-scores of SNP partitions
        if args.compute_ldscores:
            self.compute_ld_scores(args)
        
        if args.compute_h2_bins:
            self.compute_h2_bins(args)
        
    


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
    #parser.add_argument('--npz-prefix', default=None, help='Prefix of npz files that encode LD matrices (used to compute LD-scores)')
    parser.add_argument('--ld-wind-cm', type=float, default=None, help='window size to be used for estimating LD-scores in units of centiMorgans (cM).')
    parser.add_argument('--ld-wind-kb', type=int, default=None, help='window size to be used for estimating LD-scores in units of Kb.')
    parser.add_argument('--ld-wind-snps', type=int, default=None, help='window size to be used for estimating LD-scores in units of SNPs.')
    parser.add_argument('--chunk-size',  type=int, default=50, help='chunk size for LD-scores calculation')
    
    #per-SNP h2 related parameters
    parser.add_argument('--q', type=float, default=100, help='The maximum ratio between the largest and smallest estimated per-SNP heritability')

    #data input/output parameters
    parser.add_argument('--sumstats', help='Input summary statistics file')
    parser.add_argument('--ref-ld-chr', help='Suffix of LD-score files (as in ldsc)')
    parser.add_argument('--w-ld-chr', help='Suffix of LD-score weights files (as in ldsc)')
    parser.add_argument('--bfile-chr', default=None, help='Prefix of plink files (used to compute LD-scores)')
    parser.add_argument('--output-prefix', required=True, help='Prefix of all PolyFun files')    
    
    #show splash screen
    splash_screen()

    #extract args
    args = parser.parse_args()
    
    #check that the output directory exists
    if os.path.isabs(args.output_prefix) and not os.path.exists(os.path.dirname(args.output_prefix)):
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
    
