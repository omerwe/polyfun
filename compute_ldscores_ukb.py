import numpy as np
import pandas as pd
import os
import sys
import time
import logging
from tqdm import tqdm
import urllib.request
from pyarrow import ArrowIOError
import tempfile
import scipy.sparse as sparse

UKBB_LD_URL = 'https://data.broadinstitute.org/alkesgroup/UKBB_LD'
REGION_LENGTH = 3000000
UKB_N=337545
META_COLUMNS = ['SNP', 'CHR', 'BP', 'A1', 'A2']


class TqdmHandler(logging.StreamHandler):
    def __init__(self):
        logging.StreamHandler.__init__(self)

    def emit(self, record):
        msg = self.format(record)
        tqdm.write(msg)
        
        
class TqdmUpTo(tqdm):
    """
        taken from: https://github.com/tqdm/tqdm/blob/master/examples/tqdm_wget.py
    """

    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None: self.total = tsize            
        self.update(b * bsize - self.n)


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
    
    
def set_snpid_index(df, copy=False):
    if copy:
        df = df.copy()
    df['A1_first'] = (df['A1'] < df['A2']) | (df['A1'].str.len()>1) | (df['A2'].str.len()>1)
    df['A1s'] = df['A2'].copy()
    df.loc[df['A1_first'], 'A1s'] = df.loc[df['A1_first'], 'A1'].copy()
    df['A2s'] = df['A1'].copy()
    df.loc[df['A1_first'], 'A2s'] = df.loc[df['A1_first'], 'A2'].copy()
    df.index = df['CHR'].astype(str) + '.' + df['BP'].astype(str) + '.' + df['A1s'] + '.' + df['A2s']
    df.index.name = 'snpid'
    df.drop(columns=['A1_first', 'A1s', 'A2s'], inplace=True)
    return df
    
    
def read_annot(annot_file):
    try:
        df_annot = pd.read_parquet(annot_file)
    except ArrowIOError:
        df_annot = pd.read_table(annot_file, delim_whitespace=True)
    
    assert 'CHR' in df_annot.columns
    assert 'SNP' in df_annot.columns
    assert 'BP' in df_annot.columns
    assert 'A1' in df_annot.columns
    assert 'A2' in df_annot.columns
    
    df_annot = set_snpid_index(df_annot)
    
    return df_annot
    
    
def compute_R2_unbiased(R, n):
    np.power(R, 2, out=R)
    R *= ((n-1)/(n-2))
    R -= 1/(n-2)
    return R
            
    
    
def load_ld_matrix(ld_dir, ld_prefix):
    
    #load the SNPs metadata
    gz_file = os.path.join(ld_dir, '%s.gz'%(ld_prefix))
    try:
        df_ld_snps = pd.read_table(gz_file)
    except ArrowIOError:
        raise IOError('Corrupt file downloaded')
    df_ld_snps.rename(columns={'rsid':'SNP', 'chromosome':'CHR', 'position':'BP', 'allele1':'A1', 'allele2':'A2'}, inplace=True, errors='raise')
    df_ld_snps = set_snpid_index(df_ld_snps)
        
    #load the LD matrix
    npz_file = os.path.join(ld_dir, '%s.npz'%(ld_prefix))
    try: 
        R = sparse.load_npz(npz_file).toarray()
        R += R.T
    except ValueError:
        raise IOError('Corrupt file downloaded')
        
    #create df_R and return it
    df_R = pd.DataFrame(R, index=df_ld_snps.index, columns=df_ld_snps.index)    
    
    
    return df_R
                
    
    
def download_ukb_ld_file(chr_num, region_start, overwrite=False, ld_dir=None):
    region_end = region_start + REGION_LENGTH
    ld_prefix = 'chr%d_%d_%d'%(chr_num, region_start, region_end)
    
    #create a temp output dir if required
    if ld_dir is None:
        ld_dir = tempfile.mkdtemp()
    
    #if the files already exist, simply return them
    gz_file = os.path.join(ld_dir, '%s.gz'%(ld_prefix))
    npz_file = os.path.join(ld_dir, '%s.npz'%(ld_prefix))    
    if not overwrite and os.path.exists(gz_file) and os.path.exists(npz_file):
        try:
            df_R = load_ld_matrix(ld_dir, ld_prefix)
            return df_R
        except IOError:
            pass
        
    #download the region files
    for suffix in ['npz', 'gz']:
        suffix_file = os.path.join(ld_dir, '%s.%s'%(ld_prefix, suffix))
        url = '%s/%s.%s'%(UKBB_LD_URL, ld_prefix, suffix)
        with TqdmUpTo(unit='B', unit_scale=True, unit_divisor=1024, miniters=1, desc='downloading %s'%(url)) as t:                  
            urllib.request.urlretrieve(url, filename=suffix_file, reporthook=t.update_to)
            
    #load the LD matrix to memory
    df_R = load_ld_matrix(ld_dir, ld_prefix)    
    return df_R
    
    


def compute_ldscores_chr(df_annot_chr, ld_dir):
    assert len(df_annot_chr['CHR'].unique()) == 1
    chr_num = df_annot_chr['CHR'].unique()[0]
    
    #sort the SNPs by BP if needed
    if not np.all(np.diff(df_annot_chr['BP'])>=0):
        df_annot_chr = df_annot_chr.sort_values('BP', ascending=True)
    assert np.all(np.diff(df_annot_chr['BP'])>=0)
    
    #check if the data is binary
    df_annot_chr_raw = df_annot_chr.drop(columns=META_COLUMNS, errors='raise')
    if np.all(df_annot_chr_raw.dtypes == np.bool):
        is_binary = True
    elif np.all([len(np.unique(df_annot_chr_raw[c]))<=2 for c in df_annot_chr_raw.columns]):
        is_binary = True
    else:
        is_binary = False    
    
    #iterate over regions
    df_ldscores_regions_list = []
    for region_start in tqdm(range(1, df_annot_chr['BP'].max()+1, REGION_LENGTH)):
    
        #extract annotations for this region only
        region_end = region_start+REGION_LENGTH
        df_annot_region = df_annot_chr.query('%d <= BP <= %d'%(region_start, region_end))
        if df_annot_region.shape[0]==0: continue
        
        #download the LD data
        df_R_region = download_ukb_ld_file(chr_num, region_start, ld_dir=ld_dir)
        
        #sync df_R_region and df_annot_region
        index_intersect = df_R_region.index.intersection(df_annot_region.index)
        if len(index_intersect)==0:
            raise ValueError('no SNPs in chromosome %d BP %d-%d had LD info'%(chr_num, region_start, region_end))
        if len(index_intersect) < df_R_region.shape[0]:
            is_keep = df_R_region.index.isin(index_intersect)
            df_R_region = df_R_region.loc[is_keep, is_keep]
        if df_annot_region.shape[0] != df_R_region.shape[0] or np.any(df_annot_region.index != df_R_region.index):
            df_annot_region = df_annot_region.loc[df_R_region.index]
        assert np.all(df_annot_region.index == df_R_region.index)
        
        #compute R2 (the unbiased estimator of R squared)
        R2_region = compute_R2_unbiased(df_R_region.values, n=UKB_N)        
        
        #compute LD scores
        df_annot_raw_region = df_annot_region.drop(columns=META_COLUMNS, errors='raise')
        annot_region = df_annot_raw_region.values
        if is_binary:
            annot_sparse = sparse.csc_matrix(annot_region)
            ld_scores_region = (annot_sparse.T.dot(R2_region)).T
        else:
            ld_scores_region = R2_region.dot(annot_region)
        df_ldscores_region = pd.DataFrame(ld_scores_region, index=df_annot_raw_region.index, columns=df_annot_raw_region.columns)
        df_ldscores_region = pd.concat((df_annot_region[META_COLUMNS], df_ldscores_region), axis=1)
        df_ldscores_region['distance_from_center'] = np.abs(df_ldscores_region['BP'] - ((region_start+region_end)//2))
        df_ldscores_regions_list.append(df_ldscores_region)
        
    #keep the best ld-score for each SNP (the one closest to its region center)
    df_ldscores_chr = pd.concat(df_ldscores_regions_list, axis=0)
    df_ldscores_chr.sort_values('distance_from_center', ascending=True, inplace=True)
    df_ldscores_chr = df_ldscores_chr.loc[~df_ldscores_chr.index.duplicated(keep='first')]
    del df_ldscores_chr['distance_from_center']
    df_ldscores_chr.sort_values(['CHR', 'BP'], ascending=True, inplace=True)
    
    return df_ldscores_chr
    
    
    
def compute_ldscores_main(df_annot, ld_dir=None):
    
    #iterate over chromosomes
    df_ldscores_chr_list = []
    for chr_num, df_annot_chr in df_annot.groupby('CHR'):
        df_ldscores_chr = compute_ldscores_chr(df_annot_chr, ld_dir=ld_dir)
        df_ldscores_chr_list.append(df_ldscores_chr)
        
    df_ldscores = pd.concat((df_ldscores_chr_list), axis=0)
    return df_ldscores
    
    

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()    
    parser.add_argument('--annot', required=True, help='annotations file')
    parser.add_argument('--out', required=True, help='output file')
    parser.add_argument('--gz-out', default=False, action='store_true', help='if specified, the output file will be a gzipped text file instead of a parquet file')
    parser.add_argument('--cache', default=None, help='the path of an LD cache directory. If not provided, LD files will be downloaded to a temporary directory')
    args = parser.parse_args()
    
    configure_logger(args.out)
    
    #check input arguments
    if args.cache is not None and not os.path.exists(args.cache):
        raise ValueError('cache directory %s doesn\'t exist'%(args.cache))
    if len(os.path.dirname(args.out))>0 and not os.path.exists(os.path.dirname(args.out)):
        raise ValueError('output directory %s doesn\'t exist'%(os.path.dirname(args.out)))
    
    #read annotations
    df_annot = read_annot(args.annot)
    
    #comptue LD-scores
    df_ldscores = compute_ldscores_main(df_annot, ld_dir=args.cache)
    
    #save LD-scores to output file
    if args.gz_out:
        df_ldscores.to_csv(args.out, sep='\t', compression='gzip', index=False, float_format='%0.3f')
    else:
        df_ldscores.to_parquet(args.out, index=False)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

