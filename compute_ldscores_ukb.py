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
from pandas.api.types import is_numeric_dtype
from polyfun_utils import configure_logger, set_snpid_index


UKBB_LD_URL = 'https://data.broadinstitute.org/alkesgroup/UKBB_LD'
REGION_LENGTH = 3000000
UKB_N=337545
META_COLUMNS = ['SNP', 'CHR', 'BP', 'A1', 'A2']



LONG_RANGE_LD_REGIONS = []
LONG_RANGE_LD_REGIONS.append({'chr':6, 'start':25500000, 'end':33500000})
LONG_RANGE_LD_REGIONS.append({'chr':8, 'start':8000000, 'end':12000000})
LONG_RANGE_LD_REGIONS.append({'chr':11, 'start':46000000, 'end':57000000})

        
class TqdmUpTo(tqdm):
    """
        taken from: https://github.com/tqdm/tqdm/blob/master/examples/tqdm_wget.py
    """

    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None: self.total = tsize            
        self.update(b * bsize - self.n)

    
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
    
    for c in df_annot.columns:
        if c in META_COLUMNS: continue
        if not is_numeric_dtype(df_annot[c]):
            raise ValueError('Annotation %s does not have numeric values'%(c))
        
    
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
        df_ld_snps = pd.read_table(gz_file, delim_whitespace=True)
    except ArrowIOError:
        raise IOError('Corrupt file downloaded')
    df_ld_snps.rename(columns={'rsid':'SNP', 'chromosome':'CHR', 'position':'BP', 'allele1':'A1', 'allele2':'A2'}, inplace=True, errors='ignore')
    assert 'SNP' in df_ld_snps.columns
    assert 'CHR' in df_ld_snps.columns
    assert 'BP' in df_ld_snps.columns
    assert 'A1' in df_ld_snps.columns
    assert 'A2' in df_ld_snps.columns
    df_ld_snps = set_snpid_index(df_ld_snps)
        
    #load the LD matrix
    npz_file = os.path.join(ld_dir, '%s.npz'%(ld_prefix))
    logging.info('Loading LD from file %s'%(npz_file))
    t0 = time.time()
    try: 
        R = sparse.load_npz(npz_file).toarray()
        R += R.T
    except ValueError:
        raise IOError('Corrupt file downloaded')
    logging.info('Done in %0.2f seconds'%(time.time() - t0))
        
    #create df_R and return it
    df_R = pd.DataFrame(R, index=df_ld_snps.index, columns=df_ld_snps.index)    
    
    
    return df_R
                
    
    
def download_ukb_ld_file(chr_num, region_start, overwrite=False, ld_dir=None, no_cache=False):
    region_end = region_start + REGION_LENGTH
    ld_prefix = 'chr%d_%d_%d'%(chr_num, region_start, region_end)
    
    #if the files already exist, simply return them
    gz_file = os.path.join(ld_dir, '%s.gz'%(ld_prefix))
    npz_file = os.path.join(ld_dir, '%s.npz'%(ld_prefix))    
    if not overwrite and os.path.exists(gz_file) and os.path.exists(npz_file):
        try:
            df_R = load_ld_matrix(ld_dir, ld_prefix)
            return df_R
        except IOError:
            pass
            
    
    
    ### if we got here, we need to download the LD files ###
        
    #download the region files
    for suffix in ['npz', 'gz']:
        suffix_file = os.path.join(ld_dir, '%s.%s'%(ld_prefix, suffix))
        url = '%s/%s.%s'%(UKBB_LD_URL, ld_prefix, suffix)
        with TqdmUpTo(unit='B', unit_scale=True, unit_divisor=1024, miniters=1, desc='downloading %s'%(url)) as t:                  
            urllib.request.urlretrieve(url, filename=suffix_file, reporthook=t.update_to)
            
    #load the LD matrix to memory
    df_R = load_ld_matrix(ld_dir, ld_prefix)
    
    #delete the downloaded files if requested
    if no_cache:
        for suffix in ['npz', 'gz']:
            suffix_file = os.path.join(ld_dir, '%s.%s'%(ld_prefix, suffix))
            os.remove(suffix_file)
    
    return df_R
    
    


def compute_ldscores_chr(df_annot_chr, ld_dir, no_cache=False):

    #create a temp output dir if required
    if ld_dir is None:
        ld_dir = tempfile.mkdtemp()

    if not os.path.exists(ld_dir):
        raise IOError('LD directory %s doesn\'t exist'%(ld_dir))

    assert len(df_annot_chr['CHR'].unique()) == 1
    chr_num = df_annot_chr['CHR'].unique()[0]
    
    #remove long-range LD regions
    for r in LONG_RANGE_LD_REGIONS:
        if r['chr'] != chr_num: continue
        is_in_r = df_annot_chr['BP'].between(r['start'], r['end'])
        if not np.any(is_in_r): continue
        logging.warning('Removing %d SNPs from long-range LD region on chromosome %d BP %d-%d'%(is_in_r.sum(), r['chr'], r['start'], r['end']))
        df_annot_chr = df_annot_chr.loc[~is_in_r]
    if df_annot_chr.shape[0]==0:
        raise ValueError('No SNPs found in chromosome %d (ater removing long-range LD regions)'%(chr_num))
    
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
        df_R_region = download_ukb_ld_file(chr_num, region_start, ld_dir=ld_dir, no_cache=no_cache)
        
        #sync df_R_region and df_annot_region
        index_intersect = df_R_region.index.intersection(df_annot_region.index)
        if len(index_intersect)==0:
            raise ValueError('no SNPs in chromosome %d BP %d-%d had LD info'%(chr_num, region_start, region_end))
        if len(index_intersect) < df_R_region.shape[0]:
            logging.warning('Only %d/%d SNPs in chromosome %d BP %d-%d have annotations info. This may severely down-bias the LD-scores'%(len(index_intersect), df_R_region.shape[0], chr_num, region_start, region_end))
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
    
    
    
def compute_ldscores_main(df_annot, ld_dir=None, no_cache=False):
    
    #iterate over chromosomes
    df_ldscores_chr_list = []
    for chr_num, df_annot_chr in df_annot.groupby('CHR'):
        df_ldscores_chr = compute_ldscores_chr(df_annot_chr, ld_dir=ld_dir, no_cache=no_cache)
        df_ldscores_chr_list.append(df_ldscores_chr)
        
    df_ldscores = pd.concat((df_ldscores_chr_list), axis=0)
    return df_ldscores
    
    
def main(args):

    #read annotations
    df_annot = read_annot(args.annot)
    
    #comptue LD-scores
    df_ldscores = compute_ldscores_main(df_annot, ld_dir=args.ld_dir, no_cache=args.no_cache)
    
    #save LD-scores to output file
    if args.gz_out:
        df_ldscores.to_csv(args.out, sep='\t', compression='gzip', index=False, float_format='%0.3f')
    else:
        df_ldscores.to_parquet(args.out, index=False)


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()    
    parser.add_argument('--annot', required=True, help='annotations file')
    parser.add_argument('--out', required=True, help='output file')
    parser.add_argument('--gz-out', default=False, action='store_true', help='if specified, the output file will be a gzipped text file instead of a parquet file')
    parser.add_argument('--ld-dir', default=None, help='the path of an LD files directory. If not provided, LD files will be downloaded to a temporary directory')
    parser.add_argument('--no-cache', default=False, action='store_true', help='If this flag is specified, the LD files will be removed from the ld-dir after downloading them to save disk space')
    args = parser.parse_args()
    
    configure_logger(args.out)
    
    #check input arguments
    if args.ld_dir is not None and not os.path.exists(args.ld_dir):
        raise ValueError('Specified LD directory %s doesn\'t exist'%(args.ld_dir))
    if len(os.path.dirname(args.out))>0 and not os.path.exists(os.path.dirname(args.out)):
        raise ValueError('output directory %s doesn\'t exist'%(os.path.dirname(args.out)))
    
    main(args)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

