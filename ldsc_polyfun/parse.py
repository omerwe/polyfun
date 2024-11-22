'''
(c) 2014 Brendan Bulik-Sullivan and Hilary Finucane

This module contains functions for parsing various ldsc-defined file formats.

'''


import numpy as np
import pandas as pd
import os
from tqdm import tqdm
import logging

def series_eq(x, y):
    '''Compare series, return False if lengths not equal.'''
    return len(x) == len(y) and (x == y).all()


def read_csv(fh, **kwargs):
    if fh.endswith('.hdf'):
        df = pd.read_hdf(fh)
        if 'usecols' in kwargs.keys():
            df = df[kwargs['usecols']]
    elif fh.endswith('.parquet'):
        df = pd.read_parquet(fh)
        if 'usecols' in kwargs.keys():
            df = df[kwargs['usecols']]
    else:        
        df = pd.read_csv(fh, delim_whitespace=True, na_values='.', **kwargs)
    
    return df
    
def set_snpid_index(df):
    df['A1_first'] = (df['A1'] < df['A2']) | (df['A1'].str.len()>1) | (df['A2'].str.len()>1)
    df['A1s'] = df['A2'].copy()
    df.loc[df['A1_first'], 'A1s'] = df.loc[df['A1_first'], 'A1'].copy()
    df['A2s'] = df['A1'].copy()
    df.loc[df['A1_first'], 'A2s'] = df.loc[df['A1_first'], 'A2'].copy()
    s_chr = df['CHR'].map(lambda c: int(c) if str(c)[0] in ['0','1','2','3','4','5,','6','7','8','9'] else c).astype(str)
    s_bp = df['BP'].astype(int).astype(str)
    df.index = s_chr + '.' + s_bp + '.' + df['A1s'] + '.' + df['A2s']
    df.index.name = 'snpid'
    df.drop(columns=['A1_first', 'A1s', 'A2s'], inplace=True)
    return df



def sub_chr(s, chr):
    '''Substitute chr for @, else append chr to the end of str.'''
    if '@' not in s:
        s += '@'

    return s.replace('@', str(chr))


def which_compression(fh):
    '''Given a file prefix, figure out what sort of compression to use.'''
    #import ipdb; ipdb.set_trace()
    if os.access(fh + '.parquet', 4):
        suffix = '.parquet'
        compression = 'parquet'        
    elif os.access(fh + '.hdf', 4):
        suffix = '.hdf'
        compression = 'hdf'
    elif os.access(fh + '.bz2', 4):
        suffix = '.bz2'
        compression = 'bz2'
    elif os.access(fh + '.gz', 4):
        suffix = '.gz'
        compression = 'gzip'
    elif os.access(fh, 4):
        suffix = ''
        compression = None
    else:
        raise IOError('Could not open {F}[./gz/bz2/.hdf/.parquet]'.format(F=fh))

    return suffix, compression


def get_compression(fh):
    '''Which sort of compression should we use with read_csv?'''
    if fh.endswith('gz'):
        compression = 'gzip'
    elif fh.endswith('bz2'):
        compression = 'bz2'
    else:
        compression = None

    return compression


def read_cts(fh, match_snps):
    '''Reads files for --cts-bin.'''
    compression = get_compression(fh)
    cts = read_csv(fh, compression=compression, header=None, names=['SNP', 'ANNOT'])
    if not series_eq(cts.SNP, match_snps):
        raise ValueError('--cts-bin and the .bim file must have identical SNP columns.')

    return cts.ANNOT.values


def sumstats(fh, alleles=True, dropna=True):
    '''Parses .sumstats files. See docs/file_formats_sumstats.txt.'''
    dtype_dict = {'SNP': str,   'Z': float, 'N': float, 'A1': str, 'A2': str}
    compression = get_compression(fh)
    usecols = ['SNP', 'CHR', 'BP', 'Z', 'N']
    #if alleles:
    usecols += ['A1', 'A2']

    try:
        x = read_csv(fh, usecols=usecols, dtype=dtype_dict, compression=compression)
    except (AttributeError, ValueError) as e:
        raise ValueError('Improperly formatted sumstats file: ' + str(e.args))

    if dropna:
        x = x.dropna(how='any')
        
    x = set_snpid_index(x)
    x.drop(columns=['CHR', 'BP'], inplace=True)


    return x


def ldscore_fromlist(flist, num=None):
    '''Sideways concatenation of a list of LD Score files.'''
    ldscore_array = []
    for fh_i, fh in enumerate(flist):
        y = ldscore(fh, num)
        if len(ldscore_array)>0:
        
            #make sure that all files contain the same SNPs in the same order
            if ((not series_eq(y.index, ldscore_array[0].index) or not series_eq(y.SNP, ldscore_array[0].SNP))):
                all_baseline_snps_found = np.all(ldscore_array[0].index.isin(y.index))
                if not all_baseline_snps_found:
                    raise ValueError('Some SNPs in the first set of annotations are not found in the one of the other sets of annotations')
                extra_snps_found = np.any(~(y.index.isin(ldscore_array[0].index)))
                if extra_snps_found:
                    logging.warning('some SNPs in one of the sets of annotations are not found in the first set of annotations. We will ignore these SNPs')
                    
                #reorder the SNPs to make sure that they're in the corret order
                y = y.loc[ldscore_array[0].index]
                assert series_eq(y.index, ldscore_array[0].index) and series_eq(y.SNP, ldscore_array[0].SNP)
            
            # keep SNP and CHR column from only the first file
            y = y.drop(columns=['SNP', 'CHR'], axis=1)

        new_col_dict = {c: c + '_' + str(fh_i) for c in y.columns if c not in ['SNP', 'CHR']}
        y.rename(columns=new_col_dict, inplace=True)
        ldscore_array.append(y)

    if len(ldscore_array)==1:
        ldscores_all = ldscore_array[0]
    else:
        #ldscores_all = pd.concat(ldscore_array, axis=1)
        ldscores_all = pd.concat(ldscore_array, axis=1)
    return ldscores_all


def l2_parser(fh, compression):
    '''Parse LD Score files'''
    x = read_csv(fh, header=0, compression=compression)
    if 'SNP' not in x.columns:
        x.reset_index(inplace=True)
        assert 'SNP' in x.columns
        assert 'CHR' in x.columns
        assert 'BP' in x.columns
    if 'A1' in x.columns and 'A2' in x.columns:
        x = set_snpid_index(x)
    else:
        ###x.set_index('SNP', inplace=True, drop=False)
        logging.warning('%s doesn\'t have A1,A2 columns'%(fh))
    x.drop(columns=['MAF', 'CM', 'A1', 'A2'], errors='ignore', inplace=True)
    return x


def annot_parser(fh, compression, frqfile_full=None, compression_frq=None, anno=None):
    '''Parse annot files'''
    df_annot = read_csv(fh, header=0, compression=compression)
    df_annot.drop(columns=['SNP', 'BP', 'CM', 'CHR', 'A1', 'A2'], inplace=True, errors='ignore')
    df_annot = df_annot.astype(float)
    if anno is not None:
        df_annot = df_annot.loc[:, [c for c in df_annot.columns if (c=='SNP' or c in anno)]]
    if frqfile_full is not None:
        df_frq = frq_parser(frqfile_full, compression_frq)
        df_annot = df_annot[(.95 > df_frq.FRQ) & (df_frq.FRQ > 0.05)]
    return df_annot


def frq_parser(fh, compression):
    '''Parse frequency files.'''
    df = read_csv(fh, header=0, compression=compression)
    if 'MAF' in df.columns:
        df.rename(columns={'MAF': 'FRQ'}, inplace=True)
    return df[['SNP', 'FRQ']]


def ldscore(fh, num=None):
    '''Parse .l2.ldscore files, split across num chromosomes. See docs/file_formats_ld.txt.'''
    suffix = '.l2.ldscore'
    if num is not None:  # num files, e.g., one per chromosome
        first_fh = sub_chr(fh, 1) + suffix
        s, compression = which_compression(first_fh)
        ###compression = None
        chr_ld = []
        for i in tqdm(range(1, num+1)):
            chr_ld.append(l2_parser(sub_chr(fh, i) + suffix + s, compression))
        x = pd.concat(chr_ld, axis=0)  # automatically sorted by chromosome
        del chr_ld
    else:  # just one file
        s, compression = which_compression(fh + suffix)
        x = l2_parser(fh + suffix + s, compression)

    #sort array if needed
    is_sorted = True
    for c in x['CHR'].unique():
        is_sorted = np.all(np.diff(x.loc[x['CHR']==c, 'BP']) >= 0)
        if not is_sorted: break            
    if not is_sorted:
        x.sort_values(by=['CHR', 'BP'], inplace=True) # SEs will be wrong unless sorted
        
        
    x.drop(columns=['BP'], inplace=True)
    
    if x.index.name == 'snpid':
        is_duplicate_snp = x.index.duplicated()
        if np.any(is_duplicate_snp):
            index_dup_snps = x.index[is_duplicate_snp]
            index_dup_snps = index_dup_snps[~index_dup_snps.duplicated()]
            error_msg = 'Duplicate SNPs were found in the input data:\n%s'%(index_dup_snps)
            raise ValueError(error_msg)
    else:
        if np.any(x['SNP'].duplicated()):
            x.drop_duplicates(subset='SNP', inplace=True)    
    return x
    


def M(fh, num=None, N=2, common=False):
    '''Parses .l{N}.M files, split across num chromosomes. See docs/file_formats_ld.txt.'''
    parsefunc = lambda y: [float(z) for z in open(y, 'r').readline().split()]
    suffix = '.l' + str(N) + '.M'
    if common:
        suffix += '_5_50'

    if num is not None:
        x = np.sum([parsefunc(sub_chr(fh, i) + suffix) for i in range(1, num + 1)], axis=0)
    else:
        x = parsefunc(fh + suffix)

    return np.array(x).reshape((1, len(x)))


def M_fromlist(flist, num=None, N=2, common=False):
    '''Read a list of .M* files and concatenate sideways.'''
    return np.hstack([M(fh, num, N, common) for fh in flist])


def annot(fh_list, num=None, frqfile=None, anno=None):
    '''
    Parses .annot files and returns an overlap matrix. See docs/file_formats_ld.txt.
    If num is not None, parses .annot files split across [num] chromosomes (e.g., the
    output of parallelizing ldsc.py --l2 across chromosomes).

    '''
    annot_suffix = ['.annot' for fh in fh_list]
    annot_compression = []
    if num is not None:  # 22 files, one for each chromosome
        for i, fh in enumerate(fh_list):
            first_fh = sub_chr(fh, 1) + annot_suffix[i]
            annot_s, annot_comp_single = which_compression(first_fh)
            annot_suffix[i] += annot_s
            annot_compression.append(annot_comp_single)

        if frqfile is not None:
            frq_suffix = '.frq'
            first_frqfile = sub_chr(frqfile, 1) + frq_suffix
            frq_s, frq_compression = which_compression(first_frqfile)
            frq_suffix += frq_s

        y = []
        M_tot = 0
        logging.info('reading .annot files to compute the overlap matrix...')
        for chr in tqdm(range(1, num + 1)):
            if frqfile is not None:
                df_annot_chr_list = [annot_parser(sub_chr(fh, chr) + annot_suffix[i], annot_compression[i],
                                                  sub_chr(frqfile, chr) + frq_suffix, frq_compression, anno=anno)
                                     for i, fh in enumerate(fh_list)]
            else:
                df_annot_chr_list = [annot_parser(sub_chr(fh, chr) + annot_suffix[i], annot_compression[i], anno=anno)
                                     for i, fh in enumerate(fh_list)]

            if anno is not None:
                list_list_c = [list(df.columns) for df in df_annot_chr_list]
                list_c = [c for c_list in list_list_c for c in c_list]
                for a in anno:
                    assert a in list_c, 'Annotation %s was not found in the annotations file'%(a)
                    
            annot_matrix_chr_list = [np.matrix(df_annot_chr) for df_annot_chr in df_annot_chr_list]
            if len(annot_matrix_chr_list)==1:
                annot_matrix_chr = annot_matrix_chr_list[0]
            else:
                annot_matrix_chr = np.hstack(annot_matrix_chr_list)
                
            y.append(np.dot(annot_matrix_chr.T, annot_matrix_chr))
            M_tot += len(df_annot_chr_list[0])

        x = sum(y)
    else:  # just one file
        for i, fh in enumerate(fh_list):
            annot_s, annot_comp_single = which_compression(fh + annot_suffix[i])
            annot_suffix[i] += annot_s
            annot_compression.append(annot_comp_single)

        if frqfile is not None:
            frq_suffix = '.frq'
            frq_s, frq_compression = which_compression(frqfile + frq_suffix)
            frq_suffix += frq_s

            df_annot_list = [annot_parser(fh + annot_suffix[i], annot_compression[i],
                                          frqfile + frq_suffix, frq_compression, anno=anno) for i, fh in enumerate(fh_list)]

        else:
            df_annot_list = [annot_parser(fh + annot_suffix[i], annot_compression[i], anno=anno)
                             for i, fh in enumerate(fh_list)]

        if anno is not None:
            list_list_c = [list(df.columns) for df in df_annot_list]
            list_c = [c for c_list in list_list_c for c in c_list]
            for a in anno:
                assert a in list_c, 'Annotation %s was not found in the annotations file'%(a)
                
        annot_matrix_list = [np.matrix(y) for y in df_annot_list]
        annot_matrix = np.hstack(annot_matrix_list)
        x = np.dot(annot_matrix.T, annot_matrix)
        M_tot = len(df_annot_list[0])

    return x, M_tot


def __ID_List_Factory__(colnames, keepcol, fname_end, header=None, usecols=None):

    class IDContainer(object):

        def __init__(self, fname):
            self.__usecols__ = usecols
            self.__colnames__ = colnames
            self.__keepcol__ = keepcol
            self.__fname_end__ = fname_end
            self.__header__ = header
            self.__read__(fname)
            self.n = len(self.df)

        def __read__(self, fname):
            end = self.__fname_end__
            if end and not fname.endswith(end):
                raise ValueError('{f} filename must end in {f}'.format(f=end))

            comp = get_compression(fname)
            self.df = pd.read_csv(fname, header=self.__header__, usecols=self.__usecols__,
                                  delim_whitespace=True, compression=comp)

            if self.__colnames__:
                self.df.columns = self.__colnames__

            if self.__keepcol__ is not None:
                self.IDList = self.df.iloc[:, [self.__keepcol__]].astype('object')

        def loj(self, externalDf):
            '''Returns indices of those elements of self.IDList that appear in exernalDf.'''
            r = externalDf.columns[0]
            l = self.IDList.columns[0]
            merge_df = externalDf.iloc[:, [0]]
            merge_df['keep'] = True
            z = pd.merge(self.IDList, merge_df, how='left', left_on=l, right_on=r,
                         sort=False)
            ii = z['keep'] == True
            return np.nonzero(ii)[0]

    return IDContainer


PlinkBIMFile = __ID_List_Factory__(['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'], 1, '.bim', usecols=[0, 1, 2, 3, 4, 5])
PlinkFAMFile = __ID_List_Factory__(['IID'], 0, '.fam', usecols=[1])
FilterFile = __ID_List_Factory__(['ID'], 0, None, usecols=[0])
AnnotFile = __ID_List_Factory__(None, 2, None, header=0, usecols=None)
ThinAnnotFile = __ID_List_Factory__(None, None, None, header=0, usecols=None)
