import numpy as np
import pandas as pd
import os
import sys
import logging
from tqdm import tqdm


SNP_COLUMNS = ['CHR', 'SNP', 'BP', 'A1', 'A2']




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
    
    
def set_snpid_index(df, copy=False):
    if copy:
        df = df.copy()
    is_indel = (df['A1'].str.len()>1) | (df['A2'].str.len()>1)
    df['A1_first'] = (df['A1'] < df['A2']) | is_indel
    df['A1s'] = df['A2'].copy()
    df.loc[df['A1_first'], 'A1s'] = df.loc[df['A1_first'], 'A1'].copy()
    df['A2s'] = df['A1'].copy()
    df.loc[df['A1_first'], 'A2s'] = df.loc[df['A1_first'], 'A2'].copy()
    df.index = df['CHR'].astype(str) + '.' + df['BP'].astype(str) + '.' + df['A1s'] + '.' + df['A2s']
    df.index.name = 'snpid'
    df.drop(columns=['A1_first', 'A1s', 'A2s'], inplace=True)
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
