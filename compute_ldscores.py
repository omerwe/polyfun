import numpy as np
import pandas as pd
import os
import sys
import time
from ldsc_polyfun import ldscore, parse
import logging
from tqdm import tqdm


def __filter__(fname, noun, verb, merge_obj):
    merged_list = None
    if fname:
        f = lambda x,n: x.format(noun=noun, verb=verb, fname=fname, num=n)
        x = parse.FilterFile(fname)
        c = 'Read list of {num} {noun} to {verb} from {fname}'
        print(f(c, len(x.IDList)))
        merged_list = merge_obj.loj(x.IDList)
        len_merged_list = len(merged_list)
        if len_merged_list > 0:
            c = 'After merging, {num} {noun} remain'
            print(f(c, len_merged_list))
        else:
            error_msg = 'No {noun} retained for analysis'
            raise ValueError(f(error_msg, 0))

        return merged_list



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
    
        
        
def compute_ldscores(args):

    # read bim/snp
    array_snps = parse.PlinkBIMFile(args.bfile+'.bim')
    df_bim = array_snps.df
    
    #read annotations 
    if args.annot is not None:
        if args.annot.endswith('.parquet'):
            df_annot = pd.read_parquet(args.annot)
        else:
            df_annot = pd.read_table(args.annot, delim_whitespace=True)
    
        #heuristically reduce df_annot to a small superset of the relevant SNPs        
        df_annot = df_annot.loc[df_annot['SNP'].isin(df_bim['SNP'])]
        
        #duplicate df_annot to make sure we have no flipped alleles that will cause a mess
        bins_index1 =     df_annot['CHR'].astype(str) + '.' \
                        + df_annot['BP'].astype(str) + '.' \
                        + df_annot['A1'] + '.' \
                        + df_annot['A2']
        df_annot2 = df_annot.copy()
        df_annot2['A2'] = df_annot['A1'].copy()
        df_annot2['A1'] = df_annot['A2'].copy()
        df_annot.index  = df_annot['CHR'].astype(str) + '.' + df_annot['BP'].astype(str) + '.' + df_annot['A1'] + '.' + df_annot['A2']
        df_annot2.index = df_annot2['CHR'].astype(str) + '.' + df_annot2['BP'].astype(str) + '.' + df_annot2['A1'] + '.' + df_annot2['A2']        
        df_annot = pd.concat([df_annot, df_annot2], axis=0)
        del df_annot2
        
        #make sure that all SNPs have a bin
        df_bim.index = df_bim['CHR'].astype(str) + '.' + df_bim['BP'].astype(str) + '.' + df_bim['A1'] + '.' + df_bim['A2']
        if np.any(~df_bim.index.isin(df_annot.index)):
            raise ValueError('Not all SNPs have annotation values')
            
        #rearrange df_annot to match the order of SNPs in the plink file
        if (df_annot.shape[0] > df_bim.shape[0]) or np.any(df_annot.index != df_bim.index):
            assert np.all(df_bim.index.isin(df_annot.index))
            df_annot = df_annot.loc[df_bim.index]
        assert np.all(df_annot['BP'] == df_bim['BP'].values)
        assert np.all(df_annot['A1'] == df_bim['A1'].values)
        assert np.all(df_annot['A2'] == df_bim['A2'].values)
        df_annot.set_index(['CHR', 'BP', 'SNP', 'A1', 'A2'], drop=True, inplace=True)

    #find #individuals in bfile
    fam_file = args.bfile+'.fam'
    df_fam = pd.read_table(fam_file, header=None, usecols=[5], delim_whitespace=True)
    n = df_fam.shape[0]

    #find keep_indivs    
    if args.keep:
        array_indivs = parse.PlinkFAMFile(args.bfile+'.fam')
        keep_indivs = __filter__(args.keep, 'individuals', 'include', array_indivs)
        logging.info('after applying --keep, %d individuals remain'%(len(keep_indivs)))
    else:
        keep_indivs= None
    
    #read plink file    
    bed_file = args.bfile+'.bed'
    geno_array = ldscore.PlinkBEDFile(bed_file, n, array_snps, keep_snps=None,
        keep_indivs=keep_indivs, mafMin=None)
    
    #if we have missing SNPs, take only relevant rows from df_annotation
    if args.annot is not None and len(geno_array.kept_snps) < df_annot.shape[0]:
        assert np.all(np.array(geno_array.kept_snps) == np.sort(np.array(geno_array.kept_snps)))
        df_annot = df_annot.iloc[geno_array.kept_snps]
        df_bim = df_bim.iloc[geno_array.kept_snps]
        
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
    annot_values = (None if args.annot is None else df_annot.values)
    ldscores = geno_array.ldScoreVarBlocks(block_left, args.chunk_size, annot=annot_values)
    
    #create an ldscores df
    if args.annot is None:
        df_ldscores = pd.DataFrame(ldscores, columns=['base'])
    else:
        df_ldscores = pd.DataFrame(ldscores, columns=df_annot.columns)
        
    #add SNP identifier columns
    df_ldscores['SNP'] = df_bim['SNP'].values
    df_ldscores['CHR'] = df_bim['CHR'].values
    df_ldscores['BP'] = df_bim['BP'].values
    df_ldscores['A1'] = df_bim['A1'].values
    df_ldscores['A2'] = df_bim['A2'].values
    df_ldscores = df_ldscores[['CHR', 'SNP', 'BP', 'A1', 'A2'] + [c for c in df_ldscores.columns if c not in ['CHR', 'SNP', 'BP', 'A1', 'A2']]]
    
    return df_ldscores
        


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()    
    parser.add_argument('--bfile', required=True, help='Plink file')
    parser.add_argument('--out', required=True, help='output file')
    parser.add_argument('--annot', default=None, help='annotations file (if not specified, compute the global LD-score)')
    parser.add_argument('--ld-wind-cm', type=float, default=None, help='window size to be used for estimating LD-scores in units of centiMorgans (cM).')
    parser.add_argument('--ld-wind-kb', type=int, default=None, help='window size to be used for estimating LD-scores in units of Kb.')
    parser.add_argument('--ld-wind-snps', type=int, default=None, help='window size to be used for estimating LD-scores in units of SNPs.')
    parser.add_argument('--chunk-size',  type=int, default=50, help='chunk size for LD-scores calculation')
    #parser.add_argument('--extract',  default=None, help='File with rsids of SNP to use')
    parser.add_argument('--keep',  default=None, help='File with ids of individuals to use')
    args = parser.parse_args()
    
    configure_logger(args.out)
    
    if args.ld_wind_cm is None and args.ld_wind_kb is None and args.ld_wind_snps is None:
        args.ld_wind_cm = 1.0
        logging.warning('no ld-wind argument specified.  We will use --ld-cm 1.0')
        
    df_ldscores = compute_ldscores(args)
    df_ldscores.to_parquet(args.out)
    
    
    
    
    
    