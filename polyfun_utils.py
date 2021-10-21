import numpy as np
import pandas as pd
import os
import sys
import logging
from tqdm import tqdm


SNP_COLUMNS = ['CHR', 'SNP', 'BP', 'A1', 'A2']
LONG_RANGE_LD_REGIONS = []
LONG_RANGE_LD_REGIONS.append({'chr':6, 'start':25500000, 'end':33500000})
LONG_RANGE_LD_REGIONS.append({'chr':8, 'start':8000000, 'end':12000000})
LONG_RANGE_LD_REGIONS.append({'chr':11, 'start':46000000, 'end':57000000})
DEFAULT_REGIONS_FILE = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'ukb_regions.tsv.gz')



class TqdmUpTo(tqdm):
    """
        taken from: https://github.com/tqdm/tqdm/blob/master/examples/tqdm_wget.py
    """

    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None: self.total = tsize            
        self.update(b * bsize - self.n)



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
        

def check_package_versions():
    from pkg_resources import parse_version
    if parse_version(pd.__version__) < parse_version('0.25.0'):
        raise ValueError('your pandas version is too old --- please update pandas')
        
    try:
        import pandas_plink
    except (ImportError, ModuleNotFoundError):
        raise ValueError('\n\nPlease install the python package pandas_plink (using either "pip install pandas-plink" or "conda install -c conda-forge pandas-plink")\n\n')
    
    
def set_snpid_index(df, copy=False, allow_duplicates=False, allow_swapped_indel_alleles=False):
    if copy:
        df = df.copy()
    is_indel = (df['A1'].str.len()>1) | (df['A2'].str.len()>1)
    alleles_are_alphabetical = df['A1'] < df['A2']
    if allow_swapped_indel_alleles:
        df['A1_first'] = alleles_are_alphabetical
    else:
        df['A1_first'] = alleles_are_alphabetical | is_indel
    df['A1s'] = df['A2'].copy()
    df.loc[df['A1_first'], 'A1s'] = df.loc[df['A1_first'], 'A1'].copy()
    df['A2s'] = df['A1'].copy()
    df.loc[df['A1_first'], 'A2s'] = df.loc[df['A1_first'], 'A2'].copy()
    df.index = df['CHR'].astype(int).astype(str) + '.' + df['BP'].astype(str) + '.' + df['A1s'] + '.' + df['A2s']
    df.index.name = 'snpid'
    df.drop(columns=['A1_first', 'A1s', 'A2s'], inplace=True)
    
    #check for duplicate SNPs
    if not allow_duplicates:
        is_duplicate_snp = df.index.duplicated()
        if np.any(is_duplicate_snp):
            df_dup_snps = df.loc[is_duplicate_snp]
            snp_colums = [c for c in ['SNP', 'CHR', 'BP', 'A1', 'A2'] if c in df.columns]
            df_dup_snps = df_dup_snps.loc[~df_dup_snps.index.duplicated(), snp_colums]
            error_msg = 'Duplicate SNPs were found in the input data:\n%s'%(df_dup_snps)
            raise ValueError(error_msg)
    return df


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



def get_file_name(args, file_type, chr_num, verify_exists=True, allow_multiple=False):
    if file_type == 'ldscores':
        file_name = args.output_prefix + '.%d.l2.ldscore.parquet'%(chr_num)
    elif file_type == 'snpvar_ridge':
        file_name = args.output_prefix + '.%d.snpvar_ridge.gz'%(chr_num)
    elif file_type == 'taus_ridge':
        file_name = args.output_prefix + '.annot_coeff_ridge.%d.txt'%(chr_num)
    elif file_type == 'taus_nn':
        file_name = args.output_prefix + '.annot_coeff_nn.%d.txt'%(chr_num)
    elif file_type == 'snpvar_ridge_constrained':
        file_name = args.output_prefix + '.%d.snpvar_ridge_constrained.gz'%(chr_num)        
    elif file_type == 'snpvar_constrained':
        file_name = args.output_prefix + '.%d.snpvar_constrained.gz'%(chr_num)        
    elif file_type == 'snpvar':
        file_name = args.output_prefix + '.%d.snpvar.gz'%(chr_num)        
    elif file_type == 'bins':
        file_name = args.output_prefix + '.%d.bins.parquet'%(chr_num)
    elif file_type == 'M':
        file_name = args.output_prefix + '.%d.l2.M'%(chr_num)
        
    elif file_type == 'annot':
        assert verify_exists
        assert allow_multiple
        file_name = []
        for ref_ld_chr in args.ref_ld_chr.split(','):
            file_name_part = ref_ld_chr + '%d.annot.gz'%(chr_num)
            if not os.path.exists(file_name_part):
                file_name_part = ref_ld_chr + '%d.annot.parquet'%(chr_num)
            file_name.append(file_name_part)
        
    elif file_type == 'ref-ld':
        assert verify_exists
        assert allow_multiple
        file_name = []
        for ref_ld_chr in args.ref_ld_chr.split(','):
            file_name_part = ref_ld_chr + '%d.l2.ldscore.gz'%(chr_num)
            if not os.path.exists(file_name_part):
                file_name_part = ref_ld_chr + '%d.l2.ldscore.parquet'%(chr_num)
            file_name.append(file_name_part)
        
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
        if allow_multiple:
            for fname in file_name:
                if not os.path.exists(fname):
                    raise IOError('%s file not found: %s'%(file_type, fname))
        else:
            if not os.path.exists(file_name):
                raise IOError('%s file not found: %s'%(file_type, file_name))
            
    return file_name
