import numpy as np; np.set_printoptions(precision=4, linewidth=200)
import pandas as pd; pd.set_option('display.width', 200)
import os
import time
import scipy.stats as stats
import logging
import gzip
from tqdm import tqdm
import tempfile
import shutil
import glob
import subprocess
from importlib import reload
from polyfun_utils import set_snpid_index, TqdmUpTo
from pyarrow import ArrowIOError
from pyarrow.lib import ArrowInvalid
from ldstore.bcor import bcor
import scipy.sparse as sparse
from pandas_plink import read_plink
from sklearn.impute import SimpleImputer
from polyfun import configure_logger, check_package_versions
import urllib.request
from urllib.parse import urlparse
from packaging.version import Version

from pathlib import Path

def is_finemap_tools_installed():
    try:
        import finemap_tools
        from finemap_tools.reader import is_tabix_installed
        if  not is_tabix_installed():
            logging.warning('finemap_tools is installed but tabix is not installed. Please install tabix to use the finemap_tools package to read tabix-indexed files')
            return False
        
        return True
    except:
        return False


def splash_screen():
    print('*********************************************************************')
    print('* Fine-mapping Wrapper')
    print('* Version 1.0.0')
    print('* (C) 2019-2022 Omer Weissbrod')
    print('*********************************************************************')
    print()


def uri_validator(x):
    '''
    code taken from: https://stackoverflow.com/questions/7160737/python-how-to-validate-a-url-in-python-malformed-or-not
    '''
    try:
        result = urlparse(x)
        return all([result.scheme, result.netloc, result.path])
    except:
        return False


def load_ld_npz(ld_prefix):

    logging.info('Loading LD file %s'%(ld_prefix))
    t0 = time.time()

    # load SNPs info
    snps_filename_parquet = ld_prefix+'.parquet'
    snps_filename_gz = ld_prefix+'.gz'
    if os.path.exists(snps_filename_parquet):
        df_ld_snps = pd.read_parquet(snps_filename_parquet)
    elif os.path.exists(snps_filename_gz):
        df_ld_snps = pd.read_table(snps_filename_gz, sep="\s+")

        # df_ld_snps.rename(columns={'allele1':'A1', 'allele2':'A2', 'position':'BP', 'chromosome':'CHR', 'rsid':'SNP'}, inplace=True, errors='ignore')
    else:
        raise ValueError('couldn\'t find SNPs file %s or %s'%(snps_filename_parquet, snps_filename_gz))

    # load LD matrix
    R_filename = ld_prefix+'.npz'
    if not os.path.exists(R_filename):
        raise IOError('%s not found'%(R_filename))
    ld_arr = sparse.load_npz(R_filename).toarray()
    ld_arr = ld_arr+ld_arr.T
    assert np.allclose(np.diag(ld_arr), 1.0)
    # assert np.all(~np.isnan(ld_arr))

    # sanity checks
    assert ld_arr.shape[0] == ld_arr.shape[1]
    if ld_arr.shape[0] != df_ld_snps.shape[0]:
        raise ValueError('LD matrix has a different number of SNPs than the SNPs file')

    logging.info('Done in %0.2f seconds'%(time.time() - t0))
    return ld_arr, df_ld_snps


def get_bcor_meta(bcor_obj):
    df_ld_snps = bcor_obj.getMeta()
    df_ld_snps.rename(columns={'rsid':'SNP', 'position':'BP', 'chromosome':'CHR', 'allele1':'A1', 'allele2':'A2'}, inplace=True, errors='raise')
    ###df_ld_snps['CHR'] = df_ld_snps['CHR'].astype(np.int64)
    df_ld_snps['BP'] = df_ld_snps['BP'].astype(np.int64)
    return df_ld_snps


def ldstore_txt_to_npz(ldstore_txt_file, z_file, output_file):
    pass


def load_ld_bcor(ld_prefix):
    bcor_file = ld_prefix+'.bcor'
    if not os.path.exists(bcor_file):
        raise IOError('%s not found'%(bcor_file))
    logging.info('Loading LD file %s'%(bcor_file))
    t0 = time.time()
    bcor_obj = bcor(bcor_file)
    df_ld_snps = get_bcor_meta(bcor_obj)
    ld_arr = bcor_obj.readCorr([])
    # assert np.all(~np.isnan(ld_arr))
    logging.info('Done in %0.2f seconds'%(time.time() - t0))
    return ld_arr, df_ld_snps


def read_ld_from_file(ld_file):
    # if ld_file is a prefix, make it into a full file name
    if not ld_file.endswith('.bcor') and not ld_file.endswith('.npz'):
        if os.path.exists(ld_file+'.npz'):
            ld_file = ld_file + '.npz'
        elif os.path.exists(ld_file+'.bcor'):
            ld_file = ld_file + '.bcor'
        else:
            raise IOError('No suitable LD file found')

    # read the LD file
    if ld_file.endswith('.bcor'):
        ld_arr, df_ld_snps = load_ld_bcor(ld_file[:-5])  # TODO: modify
    elif ld_file.endswith('.npz'):
        ld_arr, df_ld_snps = load_ld_npz(ld_file[:-4])  # TODO:modify
    else:
        raise ValueError('unknown LD format')
    is_na_ld = np.all(~np.isnan(ld_arr)) # only keep this and I suppose this check is no need, only thing could do is to avoid the nan in the ld_arr, and drop them all.
    logging.warning(f"there are {np.isnan(ld_arr).sum(axis=1).shape[0]} nan in R matrix in the ld_arr")

    return ld_arr, df_ld_snps


def download_ld_file(url_prefix):
    temp_dir = tempfile.mkdtemp()
    filename_prefix = os.path.join(temp_dir, 'ld')
    for suffix in ['npz', 'gz']:
        url = url_prefix + '.' + suffix
        suffix_file = filename_prefix + '.' + suffix
        with TqdmUpTo(unit='B', unit_scale=True, unit_divisor=1024, miniters=1, desc='downloading %s'%(url)) as t:                  
            try:
                urllib.request.urlretrieve(url, filename=suffix_file, reporthook=t.update_to)
            except urllib.error.HTTPError as e:
                if e.code == 404:
                    raise ValueError('URL %s wasn\'t found'%(url))
                else:
                    raise

    return filename_prefix


def run_executable(cmd, description, good_returncode=0, measure_time=True, check_errors=True, show_output=False, show_command=False):
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    logging.info('Running %s...'%(description))
    if show_command:
        logging.info('Command: %s'%(' '.join(cmd)))
    t0 = time.time()
    stdout = []
    if show_output:
        for line in proc.stdout:
            if len(line.strip()) > 0:
                line_str = line.strip().decode("utf-8")
                stdout.append(line_str)
                print(line_str)
        print()  #  Why: to add a newline after the last line of the output
        stdout = '\n'.join(stdout)
        _, stderr = proc.communicate()
    else:
        stdout, stderr = proc.communicate()
        if stdout is not None:
            stdout = stdout.decode('ascii')
            if len(stdout)==0: stdout=None
    if stderr is not None:
        stderr = stderr.decode('ascii')
        if len(stderr)==0: stderr=None        

    # if (stderr is not None or proc.returncode != good_returncode):
    if proc.returncode != good_returncode:
        if stderr is not None:            
            logging.error('stderr:\n%s'%(stderr))
        if stdout is not None and not show_output:            
            logging.error('stdout:\n%s'%(stdout))
        raise RuntimeError('%s error'%(description))
    if measure_time:
        logging.info('done in %0.2f seconds'%(time.time() - t0))

    if check_errors and stdout is not None:        
        for l in stdout.split('\n'):
            if 'error' in l.lower():
                logging.error(l)
                raise RuntimeError('%s reported an error'%(description))
    if check_errors and stderr is not None:
        for l in stderr.split("\n"):
            if "error" in l.lower():
                logging.error(l)
                raise RuntimeError("%s reported an error" % (description))

    return stdout, stderr


def save_ld_to_npz(ld_arr, df_ld_snps, npz_file):

    assert npz_file.endswith('.npz')
    logging.info('Saving LD file %s'%(npz_file))
    t0 = time.time()

    # save meta file
    meta_file = npz_file[:-4] + '.gz'
    df_ld_snps.to_csv(meta_file, sep='\t', index=False)

    # save .npz file
    R = np.tril(ld_arr).astype(np.float64)
    np.fill_diagonal(R, np.diag(R)/2.0)    
    R = sparse.coo_matrix(R)
    sparse.save_npz(npz_file, R, compressed=True)
    logging.info('Done in %0.2f seconds'%(time.time() - t0))


def load_sumstats(sumstats_file, chr_num, allow_swapped_indel_alleles=False):
    # read sumstats and filter to target chromosome only
    logging.info("Loading sumstats file...")
    t0 = time.time()


    if (
        sumstats_file.endswith(".gz")
        and Path(sumstats_file + ".tbi").exists()
        and is_finemap_tools_installed()
    ):

        from finemap_tools.reader import tabix_reader

        df_sumstats = tabix_reader(sumstats_file, region=f"{chr_num}")
        if all([isinstance(i, str) for i in df_sumstats.iloc[0].values]):
            logging.info(
                f"this tabix with header info in simply load and the first row is {df_sumstats.iloc[:1].values}"
            )

    else:
        try:
            df_sumstats = pd.read_parquet(sumstats_file)
        except (ArrowIOError, ArrowInvalid):
            df_sumstats = pd.read_table(sumstats_file, sep='\s+')
        if not np.any(df_sumstats['CHR'] == chr_num):
            raise IOError('sumstats file does not include any SNPs in chromosome %s'%(chr_num))
        if np.any(df_sumstats['CHR'] != chr_num):
            df_sumstats = df_sumstats.query('CHR==%s'%(chr_num)).copy()

    df_sumstats = set_snpid_index(
        df_sumstats, allow_swapped_indel_alleles=allow_swapped_indel_alleles
    )

    if "P" not in df_sumstats.columns:
        df_sumstats["P"] = stats.chi2(1).sf(df_sumstats["Z"] ** 2)
    logging.info(
        "Loaded sumstats for %d SNPs in %0.2f seconds"
        % (df_sumstats.shape[0], time.time() - t0)
    )
    ## filter pipline provided by finemap_tools (if not installed will pass)

    if is_finemap_tools_installed:
        from finemap_tools.utils import add_ID
        from finemap_tools.snpfilter import filter_pipline

        df_sumstats["added_id"] = add_ID(df_sumstats, ["CHR", "BP", "A1", "A2"])

        logging.info(
            f"filtering SNP by finemap_tools with {df_sumstats.shape[0]} SNP at begining........"
        )
        df_sumstats = filter_pipline(sumstats=df_sumstats, id_col="added_id")

        logging.info(f"after filtering, left {df_sumstats.shape[0]} SNP")
        df_sumstats = df_sumstats.drop(columns=["added_id"])
    return df_sumstats


class Fine_Mapping(object):
    def __init__(
        self,
        genotypes_file,
        sumstats_file,
        n,
        chr_num,
        ldstore_exe,
        sample_file=None,
        incl_samples=None,
        cache_dir=None,
        cache_format=None,
        n_threads=None,
        memory=None,
        allow_swapped_indel_alleles=False,
    ):

        # check that data is valid
        if genotypes_file is not None:
            if genotypes_file.endswith(".bgen"):
                if sample_file is None:
                    raise IOError("sample-file must be provided with a bgen file")

        df_sumstats = load_sumstats(
            sumstats_file,
            chr_num,
            allow_swapped_indel_alleles=allow_swapped_indel_alleles,
        )

        # # read sumstats and filter to target chromosome only
        # logging.info('Loading sumstats file...')
        # t0 = time.time()
        # try:
        #     df_sumstats = pd.read_parquet(sumstats_file)
        # except (ArrowIOError, ArrowInvalid):
        #     df_sumstats = pd.read_table(sumstats_file, sep='\s+')
        # if not np.any(df_sumstats['CHR'] == chr_num):
        #     raise IOError('sumstats file does not include any SNPs in chromosome %s'%(chr_num))
        # if np.any(df_sumstats['CHR'] != chr_num):
        #     df_sumstats = df_sumstats.query('CHR==%s'%(chr_num)).copy()
        # df_sumstats = set_snpid_index(df_sumstats, allow_swapped_indel_alleles=allow_swapped_indel_alleles)
        # if 'P' not in df_sumstats.columns:
        #     df_sumstats['P'] = stats.chi2(1).sf(df_sumstats['Z']**2)
        # logging.info('Loaded sumstats for %d SNPs in %0.2f seconds'%(df_sumstats.shape[0], time.time()-t0))

        # ## filter pipline provided by finemap_tools (if not installed will pass)
        # try:
        #     from finemap_tools.utils import add_ID
        #     from finemap_tools.snpfilter import filter_pipline

        #     df_z["added_id"] = add_ID(df_z, ["CHR", "BP", "A1", "A2"])

        #     logging.info(
        #         f"filtering SNP by finemap_tools with {df_z.shape[0]} SNP at begining........"
        #     )
        #     df_z = filter_pipline(df_z, id_col="added_id")

        #     logging.info(f"after filtering, left {df_z.shape[0]} SNP")
        #     df_z = df_z.drop(columns=["added_id"])
        # except:
        #     logging.info("finemap_tools not installed, will not filter SNPs")
        #     pass

        # save class members
        self.genotypes_file = genotypes_file
        self.n = n
        self.sample_file = sample_file
        self.df_sumstats = df_sumstats
        self.incl_samples = incl_samples
        self.ldstore_exe = ldstore_exe
        self.cache_dir = cache_dir
        self.cache_format = cache_format
        self.n_threads = n_threads
        self.chr = chr_num
        self.memory = memory
        self.allow_swapped_indel_alleles = allow_swapped_indel_alleles

    def sync_ld_sumstats(self, ld_arr, df_ld_snps, allow_missing=False):
        df_ld_snps = set_snpid_index(df_ld_snps, allow_swapped_indel_alleles=self.allow_swapped_indel_alleles)

        if ld_arr is None:
            df_ld = pd.DataFrame(np.zeros(len(df_ld_snps.index), dtype=np.int64), index=df_ld_snps.index, columns=['dummy'])
        else:
            assert ld_arr.shape[0] == df_ld_snps.shape[0]
            assert ld_arr.shape[0] == ld_arr.shape[1]
            df_ld = pd.DataFrame(ld_arr, index=df_ld_snps.index, columns=df_ld_snps.index)
        # TODO: rm some LD with NaNs

        have_na_snp = df_ld[df_ld.isna().sum() >= 1].index.tolist()
        logging.info(
            f"remove {len(have_na_snp)} SNPs with NaNs in LD, this should be think carefully."
        )
        df_ld = df_ld.drop(index=have_na_snp, columns=have_na_snp)

        # make sure that all SNPs in the sumstats file are in the LD file
        snps_in_ld_file = self.df_sumstats_locus.index.isin(df_ld.index)
        if not np.all(snps_in_ld_file):
            # Could the missing SNPs be due to mismatched indel alleles?
            df_sumstats_missing = self.df_sumstats_locus.loc[~snps_in_ld_file]
            num_missing_is_indel = np.sum((df_sumstats_missing['A1'].str.len()>1) | (df_sumstats_missing['A2'].str.len()>1))
            if allow_missing:
                num_missing = np.sum(~snps_in_ld_file)
                logging.warning('%d variants with sumstats were not found in the LD file and will be omitted (please note that this may lead to false positives if the omitted SNPs are causal!)'%(num_missing))            
                if num_missing_is_indel > 0 and not self.allow_swapped_indel_alleles:
                    logging.warning('%d of the missing variants were indels. Check that the allele order (A1/A2) matches between the sumstats and the LD file. Also see the flag --allow-swapped-indel-alleles'%(num_missing_is_indel))
                self.df_sumstats_locus = self.df_sumstats_locus.loc[snps_in_ld_file]
                assert np.all(self.df_sumstats_locus.index.isin(df_ld.index))
            else:
                error_msg = ('not all SNPs in the sumstats file were found in the LD matrix!'
                            ' You could drop the missing SNPs with the flag --allow-missing, but please note that'
                            ' these omitted SNPs may be causal, in which case you may get false positive results...'
                            ' If there should be no missing SNPs (e.g. you are using insample LD), see the flag --allow-swapped-indel-alleles')
                raise ValueError(error_msg)

        # make sure that df_sumstats_locus is not empty
        assert self.df_sumstats_locus.shape[0] > 0, 'no SNPs found in df_sumstats_locus. Please double-check that the SNPs in the LD file and in the sumstats file have the exact same positions'

        # filter LD to only SNPs found in the sumstats file
        assert not np.any(self.df_sumstats_locus.index.duplicated())
        if df_ld.shape[0] != self.df_sumstats_locus.shape[0] or np.any(df_ld.index != self.df_sumstats_locus.index):
            if ld_arr is None:
                df_ld = df_ld.loc[self.df_sumstats_locus.index]
            else:
                df_ld = df_ld.loc[self.df_sumstats_locus.index, self.df_sumstats_locus.index]
            df_ld_snps = df_ld_snps.loc[df_ld.index]

        # do a final verification that we're synced
        assert np.all(df_ld.index == self.df_sumstats_locus.index)
        assert np.all(df_ld_snps.index == self.df_sumstats_locus.index)

        # add leading zero to sumstats CHR column if needed
        if np.any(df_ld_snps['CHR'].astype(str).str.startswith('0')):
            self.df_sumstats_locus = self.df_sumstats_locus.copy()
            self.df_sumstats_locus['CHR'] = self.df_sumstats_locus['CHR'].astype(str)
            is_1digit = self.df_sumstats_locus['CHR'].str.len()==1
            self.df_sumstats_locus.loc[is_1digit, 'CHR'] = '0' + self.df_sumstats_locus.loc[is_1digit, 'CHR']

        # update self.df_ld
        self.df_ld = df_ld
        self.df_ld_snps = df_ld_snps
        assert self.df_ld.notnull().all().all()

    def find_cached_ld_file(self, locus_start, locus_end, need_bcor=False):

        # if there's no cache dir, return None
        if self.cache_dir is None:
            return None

        if self.incl_samples is None:
            fname_pattern = '%s.%d'%(os.path.basename(self.genotypes_file), self.chr)
        else:
            fname_pattern = '%s.%s.%d'%(os.path.basename(self.genotypes_file), os.path.basename(self.incl_samples), self.chr)

        # search for suitable LD files
        bcor_files = glob.glob(os.path.join(self.cache_dir, fname_pattern+'*.bcor'))
        npz_files = glob.glob(os.path.join(self.cache_dir, fname_pattern+'*.npz'))
        if need_bcor:
            ld_files = bcor_files + npz_files
        else:
            ld_files = npz_files + bcor_files

        for ld_file in ld_files:
            if os.stat(ld_file).st_size==0:
                os.remove(ld_file)
                continue
            ld_basename = os.path.basename(ld_file)
            bp1 = int(ld_basename.split('.')[-3])
            bp2 = int(ld_basename.split('.')[-2])
            assert bp1 < bp2
            if (bp1 > locus_start) or (bp2 < locus_end): continue

            # get the list of SNPs in the LD file
            if ld_file.endswith('.npz'):
                meta_file = ld_file[:-4] + '.gz'
                if not os.path.exists(meta_file): continue
                df_ld_snps = pd.read_table(
                    meta_file,
                    sep="\s+",
                    usecols=["A1", "A2", "BP", "CHR", "SNP"],
                )
            elif ld_file.endswith('.bcor'):
                bcor_obj = bcor(ld_file)
                df_ld_snps = bcor_obj.getMeta()
                del bcor_obj
                df_ld_snps.rename(columns={'rsid':'SNP', 'position':'BP', 'chromosome':'CHR', 'allele1':'A1', 'allele2':'A2'}, inplace=True, errors='raise')
                ###df_ld_snps['CHR'] = df_ld_snps['CHR'].astype(np.int64)
                df_ld_snps['BP'] = df_ld_snps['BP'].astype(np.int64)
            else:
                raise IOError("unknown file extension")

            df_ld_snps = set_snpid_index(df_ld_snps, allow_swapped_indel_alleles=self.allow_swapped_indel_alleles)

            # make sure that the LD file includes data for all the SNPs in the locus
            if not np.all(self.df_sumstats_locus.index.isin(df_ld_snps.index)):
                logging.warning('The available cached LD file was ignored because it does not contain data for all the SNPs in the locus')
                continue

            # if we got here than we found a suitable d file
            logging.info('Found a cached LD file containing all SNPs with sumstats in chromosome %d BP %d-%d: %s'%(self.chr, locus_start, locus_end, ld_file))
            return ld_file

    def get_ld_output_file_prefix(self, locus_start, locus_end, output_dir=None):
        if self.cache_dir is None:
            if output_dir is None:
                output_dir = tempfile.mkdtemp()    
            output_prefix = os.path.join(output_dir, 'ld')
        else:
            if self.incl_samples is None:
                output_prefix = os.path.join(self.cache_dir, '%s.%d.%d.%d'%(os.path.basename(self.genotypes_file), self.chr, locus_start, locus_end))
            else:
                output_prefix = os.path.join(self.cache_dir, '%s.%s.%d.%d.%d'%(os.path.basename(self.genotypes_file), os.path.basename(self.incl_samples), self.chr, locus_start, locus_end))

        return output_prefix

    def compute_ld_bgen(self, locus_start, locus_end, verbose=False):

        # create df_z
        df_z = self.df_sumstats_locus[['SNP', 'CHR', 'BP', 'A1', 'A2']].copy()

        try:
            import bgen
        except (ImportError, ModuleNotFoundError):
            raise ValueError('\n\nPlease install the bgen package (using "pip install bgen")')
        from bgen.reader import BgenFile
        bfile = BgenFile(self.genotypes_file)
        bgen_chromosomes = bfile.chroms()
        if bgen_chromosomes[0].startswith('0'):
            df_z['CHR'] = '0' + df_z['CHR'].astype(str)

        # sync the order of the alleles between the sumstats and the bgen file
        list_bgen = []
        # rsids = bfile.rsids()
        # small change reduces the time for bgen processing
        # the previous implementation would iterate through all the SNPs in the bgen file
        # this implementation loops over just the snps in the locus
        rsids = bfile.fetch(self.chr, locus_start, locus_end)
        for snp_i, rsid in enumerate(rsids):
            #             if rsid not in df_z['SNP'].values: continue
            #             snp_alleles = bfile[snp_i].alleles
            #             snp_chrom = bfile[snp_i].chrom
            #             snp_pos = bfile[snp_i].pos
            if rsid.rsid not in df_z["SNP"].values:
                continue  # NOTE: so this indeed take the intersection of the two sets
            snp_alleles = rsid.alleles
            snp_chrom = rsid.chrom
            snp_pos = rsid.pos
            rsid = rsid.rsid  # NOTE: this is the change
            assert len(snp_alleles) == 2, 'cannot handle SNPs with more than two alleles'
            df_snp = df_z.query('SNP == "%s"' % (rsid))

            assert df_snp.shape[0]==1
            a1, a2 = df_snp['A1'].iloc[0], df_snp['A2'].iloc[0]
            snp_a1, snp_a2 = snp_alleles[0], snp_alleles[1]
            if set([a1,a2]) != set([snp_a1, snp_a2]):
                raise ValueError('The alleles for SNP %s are different in the sumstats and in the bgen file:\n \
                                 bgen:     A1=%s  A2=%s\n \
                                 sumstats: A1=%s  A2=%s \
                                '%(rsid, snp_alleles[0], snp_alleles[1], a1, a2))
            d = {'SNP':rsid, 'CHR':snp_chrom, 'BP':snp_pos, 'A1':snp_a1, 'A2':snp_a2}
            list_bgen.append(d)
        df_bgen = pd.DataFrame(list_bgen)
        df_bgen = set_snpid_index(df_bgen, allow_swapped_indel_alleles=self.allow_swapped_indel_alleles)
        df_z = set_snpid_index(df_z, allow_swapped_indel_alleles=self.allow_swapped_indel_alleles)
        df_z = df_z[[]].merge(df_bgen, left_index=True, right_index=True)
        df_z = df_z[['SNP', 'CHR', 'BP', 'A1', 'A2']]

        # rename columns
        df_z.rename(columns={'SNP':'rsid', 'CHR':'chromosome', 'BP':'position', 'A1':'allele1', 'A2':'allele2'}, inplace=True)

        # Create LDstore input files
        temp_dir = tempfile.mkdtemp()    
        incl_file = os.path.join(
            temp_dir, "incl.incl"
        )  # NOTE: remove this line if not needed
        master_file = os.path.join(temp_dir, 'master.master')
        z_file = os.path.join(temp_dir, 'chr%s.%s_%s.z'%(self.chr, locus_start, locus_end))
        dose_file = os.path.join(temp_dir, 'dosages.bdose')
        df_z.to_csv(z_file, sep=' ', index=False)

        # find number of samples
        if (
            self.incl_samples is None
        ):  # NOTE: looks prety  bad as this is only the reason of wheather there is a header in it
            num_samples = pd.read_table(self.sample_file).shape[0]-1
        else:
            num_samples = pd.read_table(self.incl_samples, header=None).shape[0]

        # get output file name
        bcor_file = os.path.join(
            temp_dir, "chr%s.%s_%s.bcor" % (self.chr, locus_start, locus_end)
        )
        ## parse for cache_format
        if self.cache_format == "npz":  #
            ld_file = os.path.join(
                temp_dir, "chr%s.%s_%s.ld" % (self.chr, locus_start, locus_end)
            )

        elif self.cache_format == "bcor":
            bcor_file = (
                self.get_ld_output_file_prefix(locus_start, locus_end, temp_dir)
                + ".bcor"
            )
        else:
            raise ValueError(f"unknown cache format {self.cache_format}")

        # Create LDstore master file
        df_master = pd.DataFrame(columns=['z','bgen','bgi','bcor','dose','sample','n_samples'])
        df_master['z'] = [z_file]
        df_master['bgen'] = [self.genotypes_file]
        df_master['bgi'] = [self.genotypes_file+'.bgi']
        df_master['bcor'] = [bcor_file]
        df_master["ld"] = [ld_file]
        df_master['bdose'] = [dose_file]
        df_master['sample'] = [self.sample_file]
        df_master['n_samples'] = num_samples
        if self.incl_samples is not None:
            df_master['incl'] = self.incl_samples
        df_master.to_csv(master_file, sep=';', header=True, index=False)    

        # run LDstore
        ldstore_cmd = [
            self.ldstore_exe,
            "--in-files",
            master_file,
            # "--write-bcor",
            # "--write-text",
            "--write-bdose",
            "--bdose-version",
            "1.0",
        ]  # TODO: maybe for checking big files or for bdose 1.1

        if self.cache_format == "npz":
            ldstore_cmd += ["--write-text"]
        elif self.cache_format == "bcor":
            ldstore_cmd += ["--write-bcor"]

        if self.memory is not None:
            ldstore_cmd += ['--memory', str(self.memory)]
        if self.n_threads is not None:
            ldstore_cmd += ['--n-threads', str(self.n_threads)]
        run_executable(ldstore_cmd, 'LDStore', measure_time=True, show_output=verbose, show_command=verbose)

        if self.cache_format == "bcor":
            if not os.path.exists(bcor_file):
                raise IOError("Could not find output BCOR file")
            return bcor_file
        elif (
            self.cache_format == "npz"
        ):  # load txt file by np and return ld_arr and df_z
            if not os.path.exists(ld_file):
                raise IOError("Could not find output LD file")

            ld_arr = np.loadtxt(ld_file)
            assert ld_arr.shape[0] == ld_arr.shape[1] == df_z.shape[0]
            df_ld_snps = df_z.rename(
                columns={
                    "allele1": "A1",
                    "allele2": "A2",
                    "position": "BP",
                    "chromosome": "CHR",
                    "rsid": "SNP",
                },
                inplace=False,
                errors="ignore",
            )
            # save_ld_to_npz(ld_arr, df_z, npz_file) # NOTE: after this function will auto save to npz file, so only need to return ld_arr, df_z
            # return
            return ld_arr, df_ld_snps

    def read_plink_genotypes(self, bed):
        X = bed.compute().astype(np.float64)
        if np.any(np.isnan(X)):
            imp = SimpleImputer(missing_values=np.nan, strategy='mean', copy=False)
            imp.fit(X)
            X = imp.transform(X)
        X -= X.mean(axis=0)
        assert not np.any(np.isnan(X))
        is_polymorphic = X.std(axis=0)>0
        X[:, is_polymorphic] /= X[:, is_polymorphic].std(axis=0)
        return X

    def compute_ld_plink_pgen(self, locus_start, locus_end, verbose):
        from finemap_tools.plink import plink2_cal_LD
        from finemap_tools.reader.plink import read_pvar
        logging.info(
            f"Computing LD from pgen fileset {self.genotypes_file} chromosome {self.chr} region {locus_start}-{locus_end}"
        )
        t0 = time.time()


        ld_snp_df = read_pvar(Path(self.genotypes_file).with_suffix(".pvar")).rename(
            columns={
                "#CHROM": "CHR",
                "ID": "SNP",
                "POS": "BP",
                "REF": "A2",  # REF as A2 is to make sure the A1 is the minor allele and consistent with the sumstats later
                "ALT": "A1",
            }
        )
        ld_snp_df = set_snpid_index(
            ld_snp_df, allow_swapped_indel_alleles=self.allow_swapped_indel_alleles
        )
        logging.info(f"read {ld_snp_df.shape[0]} SNPs from pvar file")

        # used sumstats to include the snp in analysis
        df_z = self.df_sumstats_locus[["SNP", "CHR", "BP", "A1", "A2"]].copy()
        df_z = set_snpid_index(
            df_z, allow_swapped_indel_alleles=self.allow_swapped_indel_alleles
        )
        df_z = df_z[[]].merge(ld_snp_df, left_index=True, right_index=True)
        ld_snp_df = df_z[["SNP", "CHR", "BP", "A1", "A2"]]
        del df_z

        logging.info(
            f"only keep {ld_snp_df.shape[0]} SNPs in the sumstats file to cal LD to avoid the long time of cal LD"
        )

        # get_ld
        with tempfile.TemporaryDirectory() as tmpdirname:
            print("created temporary directory", tmpdirname)
            ld_df = plink2_cal_LD(
                pgen=self.genotypes_file,
                snplist=ld_snp_df["SNP"].tolist(),
                outSuffix=tmpdirname + "/test",
                thread=self.n_threads,
                memory=self.memory * 1024,
            )

        assert ld_df.shape[0] == ld_snp_df.shape[0]
        assert (
            ld_df.index.tolist() == ld_snp_df["SNP"].tolist() == ld_df.columns.tolist()
        )

        return ld_df.values, ld_snp_df

    def compute_ld_plink(self, locus_start, locus_end, verbose):
        logging.info('Computing LD from plink fileset %s chromosome %s region %s-%s'%(self.genotypes_file, self.chr, locus_start, locus_end))
        t0 = time.time()

        # read the plink file
        df_bim, df_fam, bed = read_plink(self.genotypes_file)
        df_bim.rename(columns={'snp':'SNP', 'pos':'BP', 'chrom':'CHR', 'a0':'A2', 'a1':'A1'}, inplace=True)
        df_bim['A1'] = df_bim['A1'].astype('str')
        df_bim['A2'] = df_bim['A2'].astype('str')
        df_bim['CHR'] = df_bim['CHR'].astype(np.int64)
        del df_bim['i']
        del df_bim['cm']
        bed = bed.T

        # zoom in on target locus
        is_snp_in_region = (df_bim['BP'].between(locus_start, locus_end)) & (df_bim['CHR']==self.chr)
        df_bim = df_bim.loc[is_snp_in_region]
        df_ld_snps = df_bim
        bed = bed[:, is_snp_in_region.values]
        assert bed.shape[1]>0, 'No SNPs found in the target region'

        # compute chunk size, using the formula MEM = bed.shape[0] * chunk_size * 4 / 2**30
        if self.memory is None:
            mem_limit = 1
        else:
            mem_limit = self.memory
        chunk_size = np.int64((np.float64(mem_limit) * 0.8) / bed.shape[0] / 4 * (2**30))
        if chunk_size==0: chunk_size=1
        if chunk_size > bed.shape[1]: chunk_size = bed.shape[1]
        num_chunks = np.int64(np.ceil(bed.shape[1] / chunk_size))
        if num_chunks>1:
            assert chunk_size * (num_chunks-2) < bed.shape[1]-1
        if chunk_size * (num_chunks-1) >= bed.shape[1]:
            num_chunks-=1

        # compute LD in chunks
        logging.info('Found %d SNPs in target region. Computing LD in %d chunks...'%(bed.shape[1], num_chunks))
        ld_arr = np.empty((bed.shape[1], bed.shape[1]), dtype=np.float64)
        for chunk_i in tqdm(range(num_chunks)):
            chunk_i_start = chunk_i*chunk_size
            chunk_i_end = np.minimum(chunk_i_start+chunk_size, bed.shape[1])
            X_i = self.read_plink_genotypes(bed[:, chunk_i_start:chunk_i_end])
            ld_arr[chunk_i_start:chunk_i_end, chunk_i_start:chunk_i_end] = X_i.T.dot(X_i) / X_i.shape[0]
            for chunk_j in range(chunk_i+1, num_chunks):
                chunk_j_start = chunk_j*chunk_size
                chunk_j_end = np.minimum(chunk_j_start+chunk_size, bed.shape[1])
                X_j = self.read_plink_genotypes(bed[:, chunk_j_start:chunk_j_end])
                ld_arr[chunk_i_start:chunk_i_end, chunk_j_start:chunk_j_end] = X_i.T.dot(X_j) / X_i.shape[0]
                ld_arr[chunk_j_start:chunk_j_end, chunk_i_start:chunk_i_end] = ld_arr[chunk_i_start:chunk_i_end, chunk_j_start:chunk_j_end].T
        ld_arr = np.nan_to_num(ld_arr, copy=False)
        ld_diag = np.diag(ld_arr).copy()
        if np.any(np.isclose(ld_diag, 0.0)):
            ld_diag[np.isclose(ld_diag, 0.0)] = 1.0
            np.fill_diagonal(ld_arr, ld_diag)

        logging.info('Done in %0.2f seconds'%(time.time() - t0))
        return ld_arr, df_ld_snps

    def set_locus(self, locus_start, locus_end):

        # update self.df_sumstats_locus
        self.df_sumstats_locus = self.df_sumstats.query('%d <= BP <= %d'%(locus_start, locus_end))
        num_snps = self.df_sumstats_locus.shape[0]
        if num_snps < 2:
            raise ValueError('%d SNP(s) found in sumstats file in the BP range %d-%d'%(num_snps, locus_start, locus_end))

    def get_ld_data(self, locus_start, locus_end, need_bcor=False, verbose=False):

        ld_arr, df_ld_snps, ld_file = None, None, None

        # check if we already have a suitable LD file in the cache dir
        ld_file = self.find_cached_ld_file(locus_start, locus_end, need_bcor=need_bcor)

        # compute LD if we couldn't find a suitable LD file
        if ld_file is None:
            if self.genotypes_file.endswith(".bgen"):  # this won't return None
                if not os.path.exists(self.genotypes_file):
                    raise IOError('%s doesn\'t exist'%(self.genotypes_file))
                if self.cache_format == "bcor":
                    ld_file = self.compute_ld_bgen(
                        locus_start, locus_end, verbose=verbose
                    )
                elif self.cache_format == "npz":
                    ld_arr, df_ld_snps = self.compute_ld_bgen(
                        locus_start, locus_end, verbose=verbose
                    )
                # ld_file = self.compute_ld_bgen(locus_start, locus_end, verbose=verbose)
            elif os.path.exists(self.genotypes_file+'.bed'):
                ld_arr, df_ld_snps = self.compute_ld_plink(locus_start, locus_end, verbose=verbose)
            elif os.path.exists(self.genotypes_file + ".pgen"):
                ld_arr, df_ld_snps = self.compute_ld_plink_pgen(
                    locus_start, locus_end, verbose=verbose
                )
            else:
                raise ValueError('no suitable file found for: %s'%(self.genotypes_file))

        # arrange the LD data
        assert ld_file is None or (ld_arr is None and df_ld_snps is None)

        # if there is no LD file, return the LD data directly

        if ld_file is None:  #  NOTE: this is not possiblely be None as the code
            # cache output if possible
            if self.cache_dir is not None:
                npz_file = self.get_ld_output_file_prefix(locus_start, locus_end) + '.npz'
                save_ld_to_npz(ld_arr, df_ld_snps, npz_file)
            return ld_arr, df_ld_snps
        # if we have an LD file, return it if it's a bcor and we want a bcor, or return the LD data directly otherwise
        else:  # NOTE: this code is quiet strange
            if ld_file.endswith('.bcor'):

                if need_bcor:
                    return ld_file
                else:
                    ld_arr, df_ld_snps = read_ld_from_file(ld_file)

                    # cache output if possible
                    if self.cache_dir is not None and ld_file.endswith('.bcor'):
                        npz_file = self.get_ld_output_file_prefix(locus_start, locus_end) + '.npz'
                        save_ld_to_npz(ld_arr, df_ld_snps, npz_file)
                    return ld_arr, df_ld_snps

            else:
                ld_arr, df_ld_snps = read_ld_from_file(ld_file)
                return ld_arr, df_ld_snps

    def finemap(self):
        raise NotImplementedError()

    def estimate_h2_hess(self, prop_keep=0.005, R_cutoff=0.99, pvalue_bound=None):
        '''
            prop_keep:  Proportion of SNPs to use in the estimation (only the ones with the smallest p-values)
            R_cutoff: Exclude one of each pair of SNPs with with magnitude of correlation greater than this value
            pvalue_bound: An upper bound on the p-value cutoff (i.e., SNPs with P greater than this cutoff will never be used in the estimation)

            The modified HESS equation implemented below is

            $$ \frac{n \alpha R^{-1} \alpha - m}{n} = \alpha R^{-1} \alpha - \frac{m}{n} $$

            where $\alpha$ is a vector of marginal effect size estimates for $m$ standardized SNPs,
            $R$ is a matrix of summary LD information, and $n$ is the sample size.

            This is a biased estimator (denominator of $n$) with a smaller estimation variance compared
            with the unbiased estimator (denominator of $n - m$) used in the original HESS publication
            (Shi et al., 2014; https://doi.org/10.1016/j.ajhg.2016.05.013).
        '''

        # keep only potential causal SNPs
        pvalue_cutoff = self.df_sumstats_locus['P'].quantile(prop_keep)
        if pvalue_cutoff==0:
            pvalue_cutoff = np.min(self.df_sumstats_locus['P'].loc[lambda p:p>0])
        if pvalue_bound is not None and pvalue_cutoff>pvalue_bound:
            pvalue_cutoff = pvalue_bound
        is_potential_csnp = self.df_sumstats_locus['P'].values<pvalue_cutoff
        if np.any(is_potential_csnp):
            R_pot_csnp = self.df_ld.loc[is_potential_csnp, is_potential_csnp].values
        else:
            return 0

        # take a maximally independent subset
        np.fill_diagonal(R_pot_csnp,0)
        import networkx as nx
        G = nx.from_numpy_array(np.abs(R_pot_csnp)>R_cutoff)
        np.fill_diagonal(R_pot_csnp,1)
        inds = np.sort(nx.maximal_independent_set(G))

        # estimate h2 using HESS
        R_subset = R_pot_csnp[np.ix_(inds, inds)]
        alpha_subset = self.df_sumstats_locus.loc[is_potential_csnp, 'Z'].iloc[inds].values / np.sqrt(self.n)
        h2_hess = alpha_subset.dot(np.linalg.solve(R_subset, alpha_subset)) - R_subset.shape[0]/self.n

        return h2_hess

    def estimate_h2_hess_wrapper(self, prop_keep=0.005, R_cutoff=0.99, min_h2=None, num_samples=100):
        '''
            prop_keep:  Proprtion of SNPs to use in the estimation (only the ones with the smallest p-values)
            R_cutoff: Exclude one of each pair of SNPs with with magnitude of correlation greater than this value
            min_h2: Exclude SNPs that tag less than this amount of heritability
            num_samples: Number of random samples of indepdendent SNPs to draw        
        '''

        if min_h2 is None:
            pvalue_bound = None
        else:
            assert min_h2 > 0 and min_h2 < 1, \
                'The minimum proportion of heritability to exclude SNPs from HESS estimation must be between 0 and 1'
            pvalue_bound = stats.chi2(1).sf(min_h2 * self.n)

        assert num_samples > 0, 'Number of random samples must be a positive integer'

        h2_hess_list = [self.estimate_h2_hess(prop_keep=prop_keep, R_cutoff=R_cutoff, pvalue_bound=pvalue_bound) \
                        for try_num in range(num_samples)]
        h2_hess = np.mean(h2_hess_list)
        return h2_hess        


class SUSIE_Wrapper(Fine_Mapping):

    def __init__(
        self,
        genotypes_file,
        sumstats_file,
        n,
        chr_num,
        ldstore_exe,
        sample_file=None,
        incl_samples=None,
        cache_dir=None,
        cache_format=None,
        n_threads=None,
        memory=None,
        allow_swapped_indel_alleles=False,
    ):

        super(SUSIE_Wrapper, self).__init__(
            genotypes_file,
            sumstats_file,
            n,
            chr_num,
            ldstore_exe=ldstore_exe,
            sample_file=sample_file,
            incl_samples=incl_samples,
            cache_dir=cache_dir,
            cache_format=cache_format,
            n_threads=n_threads,
            memory=memory,
            allow_swapped_indel_alleles=allow_swapped_indel_alleles,
        )

        # load SuSiE R package
        import rpy2
        import rpy2.robjects.numpy2ri as numpy2ri
        import rpy2.robjects as ro
        ro.conversion.py2ri = numpy2ri
        numpy2ri.activate()
        from rpy2.robjects.packages import importr
        self.susieR = importr('susieR')
        self.R_null = ro.rinterface.NULL
        # self.RNULLType = rpy2.rinterface.RNULLType

    def finemap(self, locus_start, locus_end, num_causal_snps, use_prior_causal_prob=True, prior_var=None, residual_var=None, residual_var_init=None, hess_resvar=False, hess=False, hess_iter=100, hess_min_h2=None, susie_max_iter=100, verbose=False, ld_file=None, debug_dir=None, allow_missing=False, susie_outfile=None, finemap_dir=None):

        # check params
        if use_prior_causal_prob and 'SNPVAR' not in self.df_sumstats.columns:
            raise ValueError('SNPVAR column not found in sumstats file')
        if hess_resvar:
            assert hess, 'hess_resvar cannot be specified if hess is FALSE'

        # set locus
        self.set_locus(locus_start, locus_end)

        # download LD file if it's a url
        if uri_validator(ld_file):
            ld_file = download_ld_file(ld_file)
            delete_ld_files_on_exit = True
        else:
            delete_ld_files_on_exit = False

        # Load LD data into memory if num_causal_snps>1
        if num_causal_snps==1:
            if hess:
                raise ValueError('Cannot use HESS-based variance estimator when assuming a single causal SNP per locus')
            self.df_ld = pd.DataFrame(np.eye(self.df_sumstats_locus.shape[0]), index=self.df_sumstats_locus.index, columns=self.df_sumstats_locus)
            self.df_ld_snps = self.df_sumstats_locus
        else:
            if ld_file is None:
                ld_arr, df_ld_snps = self.get_ld_data(locus_start, locus_end, need_bcor=False, verbose=verbose)
            else:
                ld_arr, df_ld_snps = read_ld_from_file(ld_file)

            # assert np.all(~np.isnan(ld_arr)) # NOTE: this may often happned with some ld is NaN, I supposed to rm them instead of raise errors; And code in self.sync_ld_sumstats will do this, so i comments this code here

            self.sync_ld_sumstats(ld_arr, df_ld_snps, allow_missing=allow_missing)
            del ld_arr
            del df_ld_snps

        # define prior causal probabilities
        if use_prior_causal_prob:
            prior_weights = self.df_sumstats_locus['SNPVAR'].copy().values
            prior_weights /= prior_weights.sum()
            assert np.isclose(prior_weights.sum(), 1)

        # flip effect sizes if needed
        assert np.all(self.df_ld_snps['BP'] == self.df_sumstats_locus['BP'])
        is_flipped = self.df_ld_snps['A1'] == self.df_sumstats_locus['A2']
        is_not_flipped = self.df_ld_snps['A1'] == self.df_sumstats_locus['A1']
        assert np.all(is_flipped | is_not_flipped)
        bhat = self.df_sumstats_locus['Z'].values.copy()
        if np.any(is_flipped):
            bhat[is_flipped.values] *= -1
            logging.info('Flipping the effect-sign of %d SNPs that are flipped compared to the LD panel'%(is_flipped.sum()))

        # Use HESS to estimate causal effect sizes
        if hess:
            if prior_var is not None:
                raise ValueError('cannot specify both hess and a custom prior_var')
            if self.n < 20000:
                logging.warning('HESS method is intended for studies with large sample sizes (i.e. >20K)')
            if hess_min_h2 is None:
                logging.warning('For best results, you should consider setting --hess-min-h2 to exclude SNPs with low heritability from the HESS estimation. You will need to experiment with your data to find a suitable heritability threshold. To start, try --hess-min-h2 1e-4')
            else:
                logging.info('Excluding SNPs with heritability less than %0.4e from the HESS estimation'%(hess_min_h2))
            h2_hess = self.estimate_h2_hess_wrapper(min_h2=hess_min_h2, num_samples=hess_iter)
            logging.info('Average local SNP heritability estimated by modified HESS over %d iterations: %0.4e'%(hess_iter, h2_hess))
            if h2_hess > 10:
                logging.warning('The HESS estimator is unconstrained, and the estimate is an order of magnitude greater than the expected max of 1. Use with caution')
            prior_var = h2_hess / num_causal_snps
            if prior_var <= 0:
                raise ValueError('HESS estimates that the locus causally explains zero heritability')
            if prior_var >= 1:
                raise ValueError('HESS-estimated prior-var >1. The HESS estimator cannot be used in this locus.')
            logging.info('HESS estimated causal effect size variance: %0.4e'%(prior_var))

            if hess_resvar:
                residual_var = 1 - h2_hess
                logging.info('Residual variance using the HESS estimate: %0.4e'%(residual_var))
                assert residual_var>=0

        # rpy2 bug fix
        import rpy2.robjects.numpy2ri as numpy2ri
        reload(numpy2ri)
        numpy2ri.activate()

        # run SuSiE
        t0 = time.time()
        m = self.df_sumstats_locus.shape[0]        
        logging.info('Starting %s SuSiE fine-mapping for chromosome %d BP %d-%d (%d SNPs)'%(
            ('functionally-informed' if use_prior_causal_prob else 'non-functionally informed'),
            self.chr,
            locus_start,
            locus_end,
            self.df_ld.shape[0]
            ))

        # save variables to debug dir if needed
        if debug_dir is not None:
            os.makedirs(debug_dir, exist_ok=True)
            logging.info('Saving debug info to: %s'%(debug_dir))
            self.df_sumstats_locus.to_csv(os.path.join(debug_dir, 'df_sumstats_locus.txt'), index=False, sep='\t')
            np.savetxt(os.path.join(debug_dir, 'bhat.txt'), bhat)
            # np.savez_compressed(os.path.join(debug_dir, 'R.npz'), R=self.df_ld.values)
            np.savetxt(os.path.join(debug_dir, 'n.txt'), [self.n])
            np.savetxt(os.path.join(debug_dir, 'L.txt'), [num_causal_snps])
            np.savetxt(os.path.join(debug_dir, 'residual_var.txt'), [np.nan] if (residual_var is None) else [residual_var])
            np.savetxt(os.path.join(debug_dir, 'prior_var.txt'), [np.nan] if (prior_var is None) else [prior_var])
            np.savetxt(os.path.join(debug_dir, 'prior_weights.txt'), prior_weights if use_prior_causal_prob else [np.nan])

            # create a zipped debug file
            import zipfile
            debug_files = glob.glob(os.path.join(debug_dir, '*.txt'))
            zip_file = os.path.join(debug_dir, 'debug.zip')
            zf = zipfile.ZipFile(zip_file, mode='w')
            for debug_file in debug_files:
                zf.write(
                    debug_file,
                    os.path.basename(debug_file),
                    compress_type=zipfile.ZIP_DEFLATED,
                )

        assert self.df_ld.notnull().all().all()
        if residual_var is not None: residual_var_init = residual_var

        if hasattr(self.susieR, 'susie_suff_stat'):
            logging.info('Using susieR::susie_suff_stat()')
            susie_obj = self.susieR.susie_suff_stat(
                    bhat=bhat.reshape((m,1)),
                    shat=np.ones((m,1)),
                    R=self.df_ld.values,
                    n=self.n,
                    L=num_causal_snps,
                    scaled_prior_variance=(0.0001 if (prior_var is None) else prior_var),
                    estimate_prior_variance=(prior_var is None),
                    residual_variance=(self.R_null if (residual_var_init is None) else residual_var_init),
                    estimate_residual_variance=(residual_var is None),
                    max_iter=susie_max_iter,
                    verbose=verbose,
                    prior_weights=(prior_weights.reshape((m,1)) if use_prior_causal_prob else self.R_null)
                )
        elif hasattr(self.susieR, 'susie_bhat'):
            logging.info('Using susieR::susie_bhat()')
            susie_obj = self.susieR.susie_bhat(
                    bhat=bhat.reshape((m,1)),
                    shat=np.ones((m,1)),
                    R=self.df_ld.values,
                    n=self.n,
                    L=num_causal_snps,
                    scaled_prior_variance=(0.0001 if (prior_var is None) else prior_var),
                    estimate_prior_variance=(prior_var is None),
                    residual_variance=(self.R_null if (residual_var is None) else residual_var),
                    estimate_residual_variance=(residual_var is None),
                    max_iter=susie_max_iter,
                    verbose=verbose,
                    prior_weights=(prior_weights.reshape((m,1)) if use_prior_causal_prob else self.R_null)
                )
        else:
            raise NotImplementedError('Only susie_suff_stat() and susie_bhat() are supported. Check your version of susieR')
        susie_time = time.time()-t0        
        logging.info('Done in %0.2f seconds'%(susie_time))

        # extract pip and beta_mean
        pip = np.array(self.susieR.susie_get_pip(susie_obj))
        beta_mean = np.array(self.susieR.coef_susie(susie_obj)[1:])
        assert np.allclose(beta_mean, np.sum(np.array(susie_obj.rx2('mu')) * np.array(susie_obj.rx2('alpha')), axis=0) / np.array(susie_obj.rx2('X_column_scale_factors')))

        # compute the posterior mean of beta^2
        s_alpha = np.array(susie_obj.rx2('alpha'))
        s_mu = np.array(susie_obj.rx2('mu'))
        s_mu2 = np.array(susie_obj.rx2('mu2'))
        s_X_column_scale_factors = np.array(susie_obj.rx2('X_column_scale_factors'))
        beta_var = np.sum(s_alpha*s_mu2 - (s_alpha*s_mu)**2, axis=0) / (s_X_column_scale_factors**2)
        assert np.all(beta_var>=0)

        # create output df
        df_susie = self.df_sumstats_locus.copy()
        df_susie['PIP'] = pip
        df_susie['BETA_MEAN'] = beta_mean
        # flip back the finemap BETA, as the alleles are in original order
        df_susie.loc[is_flipped, 'BETA_MEAN'] *= (-1)
        df_susie['BETA_SD'] = np.sqrt(beta_var)

        # add distance from center
        start = df_susie['BP'].min()
        end = df_susie['BP'].max()
        middle = (start+end)//2
        df_susie['DISTANCE_FROM_CENTER'] = np.abs(df_susie['BP'] - middle)        

        # mark causal sets
        import rpy2
        logging.info('Using rpy2 version %s'%(rpy2.__version__))
        if Version(rpy2.__version__) >= Version('3.5.9'):
            snames = (susie_obj.names).tolist()
            self.susie_dict = {key: np.array(susie_obj.rx2(key), dtype=object) for key in snames}
        else:
            self.susie_dict = {key:np.array(susie_obj.rx2(key), dtype=object) for key in list(susie_obj.names)}
        df_susie['CREDIBLE_SET'] = 0
        susie_sets = self.susie_dict['sets'][0]
        # if type(susie_sets) != self.RNULLType:
        try:
            for set_i, susie_set in enumerate(susie_sets):
                is_in_set = np.zeros(df_susie.shape[0], dtype=bool)
                is_in_set[np.array(susie_set)-1] = True
                is_in_set[df_susie['CREDIBLE_SET']>0] = False
                df_susie.loc[is_in_set, 'CREDIBLE_SET'] = set_i+1
        except TypeError:
            pass

        # save SuSiE object if requested
        if susie_outfile is not None:
            from rpy2.robjects.packages import importr
            R_base = importr('base', robject_translations = {'print.me': 'print_dot_me', 'print_me': 'print_uscore_me'})
            R_base.saveRDS(susie_obj, file=susie_outfile)
            logging.info('Saved SuSiE object to RDS file: %s'%(susie_outfile))

        # delete the LD file if needed
        if delete_ld_files_on_exit:
            ld_file_dir = os.path.dirname(ld_file)
            if os.path.exists(ld_file_dir): shutil.rmtree(ld_file_dir)

        return df_susie


class FINEMAP_Wrapper(Fine_Mapping):

    def __init__(
        self,
        genotypes_file,
        sumstats_file,
        n,
        chr_num,
        finemap_exe,
        ldstore_exe,
        sample_file=None,
        incl_samples=None,
        cache_dir=None,
        cache_format=None,
        n_threads=None,
        memory=None,
        allow_swapped_indel_alleles=False,
    ):

        super(FINEMAP_Wrapper, self).__init__(
            genotypes_file,
            sumstats_file,
            n,
            chr_num,
            ldstore_exe=ldstore_exe,
            sample_file=sample_file,
            incl_samples=incl_samples,
            cache_dir=cache_dir,
            cache_format=cache_format,
            n_threads=n_threads,
            memory=memory,
            allow_swapped_indel_alleles=allow_swapped_indel_alleles,
        )
        self.finemap_exe = finemap_exe

    def finemap(self, locus_start, locus_end, num_causal_snps, use_prior_causal_prob=True, prior_var=None, residual_var=None, hess=False, hess_iter=100, hess_min_h2=None, susie_max_iter=100, verbose=False, ld_file=None, debug_dir=None, allow_missing=False, susie_outfile=None, residual_var_init=None, hess_resvar=False, finemap_dir=None):

        # check params
        if use_prior_causal_prob and 'SNPVAR' not in self.df_sumstats.columns:
            raise ValueError('SNPVAR column not found in sumstats file')
        if hess:
            raise ValueError('FINEMAP cannot be used with a HESS-based variance estimator')
        if residual_var is not None:
            raise ValueError('cannot specify residual_var for FINEMAP')
        if debug_dir is not None:
            raise NotImplementedError('FINEMAP object does not support --debug-dir')
        if hess_resvar:
            raise NotImplementedError('FINEMAP object does not support --susie-resvar-hess')
        if residual_var_init is not None:
            raise NotImplementedError('FINEMAP object does not support --susie-resvar-init')
        # if allow_missing:
        # raise ValueError('FINEMAP object does not support --allow-missing')

        # download LD file if it's a url
        if uri_validator(ld_file):
            ld_file = download_ld_file(ld_file)            

        # create prefix of output files
        if finemap_dir is None:
            finemap_dir = tempfile.mkdtemp()
        else:
            os.makedirs(finemap_dir, exist_ok=True)
            logging.info('Saving FINEMAP files to directory: %s'%(finemap_dir))
        assert os.path.isdir(finemap_dir)
        finemap_output_prefix = os.path.join(finemap_dir, 'finemap')

        # set locus
        self.set_locus(locus_start, locus_end)

        # find or create a suitable ld_file
        if num_causal_snps==1:
            if ld_file is not None:
                raise ValueError('cannot specify an ld file when assuming a single causal SNP per locus')
            ld_file = finemap_output_prefix+'.ld'
            np.savetxt(ld_file, np.eye(self.df_sumstats_locus.shape[0], dtype=np.int64), fmt='%s')
        else:
            if ld_file is None:
                ld_data = self.get_ld_data(locus_start, locus_end, need_bcor=True, verbose=verbose)
                if isinstance(ld_data, str):
                    ld_file = ld_data
                    assert ld_file.endswith('.bcor')
                    assert os.path.exists(ld_file)
                elif isinstance(ld_data, tuple):
                    assert len(ld_data)==2
                    ld_arr, df_ld_snps = ld_data[0], ld_data[1]
                    self.sync_ld_sumstats(ld_arr, df_ld_snps, allow_missing=allow_missing)
                    del ld_arr, df_ld_snps
                    ld_file = finemap_output_prefix + '.ld'
                    np.savetxt(ld_file, self.df_ld.values, fmt='%0.5f')
            elif ld_file.endswith('.bcor'):
                pass
            elif ld_file.endswith('.npz') or os.path.exists(ld_file+'.npz'):
                ld_arr, df_ld_snps = read_ld_from_file(ld_file)
                self.sync_ld_sumstats(ld_arr, df_ld_snps, allow_missing=allow_missing)
                del ld_arr, df_ld_snps
                ld_file = finemap_output_prefix + ".ld"
                np.savetxt(ld_file, self.df_ld.values, fmt="%0.5f")
            else:
                raise ValueError("unknown LD file format for file: %s" % (ld_file))

        # define file names
        master_file = finemap_output_prefix+'.master'
        snp_filename = finemap_output_prefix+'.snp'
        config_filename = finemap_output_prefix+'.config'
        cred_filename = finemap_output_prefix+'.cred'
        log_filename = finemap_output_prefix+'.log'        
        z_filename = finemap_output_prefix+'.z'

        # flip some of the alleles
        if num_causal_snps == 1:
            is_flipped = np.zeros(self.df_sumstats_locus.shape[0], dtype=bool)
        else:
            if ld_file.endswith('.bcor'):
                bcor_obj = bcor(ld_file)
                df_ld_snps = get_bcor_meta(bcor_obj)
                self.sync_ld_sumstats(None, df_ld_snps, allow_missing=allow_missing)
                del df_ld_snps
            assert np.all(self.df_ld_snps['BP'] == self.df_sumstats_locus['BP'])
            is_flipped = self.df_ld_snps['A1'] == self.df_sumstats_locus['A2']
            is_not_flipped = self.df_ld_snps['A1'] == self.df_sumstats_locus['A1']
            assert np.all(is_flipped | is_not_flipped)
            if np.any(is_flipped):
                logging.info('Flipping the effect-sign of %d SNPs that are flipped compared to the LD panel'%(is_flipped.sum()))

        # create df_z and save it to disk
        df_z = self.df_sumstats_locus[['SNP', 'CHR', 'BP', 'A1', 'A2', 'Z']].copy()
        df_z.loc[is_flipped, 'A1'] = self.df_sumstats_locus.loc[is_flipped, 'A2']
        df_z.loc[is_flipped, 'A2'] = self.df_sumstats_locus.loc[is_flipped, 'A1']
        df_z.loc[is_flipped, 'Z'] *= (-1)

        df_z.rename(columns={'SNP':'rsid', 'CHR':'chromosome', 'BP':'position', 'A1':'allele1', 'A2':'allele2', 'Z':'beta'}, inplace=True, errors='ignore')
        df_z['se'] = 1
        df_z['maf'] = 0.05
        df_z = df_z[['rsid', 'chromosome', 'position', 'allele1', 'allele2', 'maf', 'beta', 'se']]
        if use_prior_causal_prob:
            df_z['prob'] = self.df_sumstats_locus['SNPVAR'] / self.df_sumstats_locus['SNPVAR'].sum()
        df_z.to_csv(z_filename, header=True, index=False, sep=' ', float_format='%0.5f')

        # create the master file
        df_master = pd.DataFrame()
        df_master['z'] = [z_filename]
        df_master['snp'] = [snp_filename]
        df_master['config'] = [config_filename]
        df_master['cred'] = [cred_filename]
        df_master['log'] = [log_filename]
        df_master['n_samples'] = [self.n]
        if ld_file.endswith('.bcor'):
            df_master['bcor'] = [ld_file]
        elif ld_file.endswith('.ld'):
            df_master['ld'] = [ld_file]
        else:
            raise ValueError('Illegal LD file format')
        df_master.to_csv(master_file, sep=';', header=True, index=False)

        # prepare the FINEMAP command
        finemap_cmd = [self.finemap_exe]        
        finemap_cmd += ['--in-files', master_file, '--sss']
        finemap_cmd += ['--force-n-samples']
        finemap_cmd += ['--log', log_filename]
        finemap_cmd += ['--n-causal-snps', str(num_causal_snps)]
        finemap_cmd += ['--std-effects']
        ###finemap_cmd += ['--flip-beta']
        if self.n_threads is not None:
            finemap_cmd += ['--n-threads', str(self.n_threads)]
        if prior_var is not None:
            finemap_cmd += ['--prior-std', str(np.sqrt(prior_var))]
        if use_prior_causal_prob:
            finemap_cmd += ['--prior-snps']

        # run FINEMAP
        t0 = time.time()
        m = self.df_sumstats_locus.shape[0]        
        logging.info('Starting %s FINEMAP fine-mapping for chromosome %d BP %d-%d (%d SNPs)'%(
            ('functionally-informed' if use_prior_causal_prob else 'non-functionally informed'),
            self.chr,
            locus_start,
            locus_end,
            self.df_sumstats_locus.shape[0]
            ))
        run_executable(finemap_cmd, 'FINEMAP', measure_time=True, show_output=verbose, show_command=verbose)
        if not os.path.exists(log_filename+'_sss'):
            raise IOError('FINEMAP output files not found')

        # load log file
        found_post_csnps = False
        with open(log_filename+'_sss') as f:
            for line in f:
                if line.startswith('- Post-expected # of causal SNPs'):
                    post_mean_num_csnps = float(line[line.index(': ')+2:-1])
                    found_post_csnps = True
                    break
        if not found_post_csnps:
            raise IOError('corrupt log file found: %s'%(log_filename+'_sss'))

        # load results
        df_finemap = pd.read_table(snp_filename, sep=' ', usecols=['rsid', 'chromosome', 'position', 'allele1', 'allele2', 'prob', 'mean', 'sd'])
        df_finemap.rename(columns={'rsid':'SNP', 'position':'BP', 'chromosome':'CHR', 'prob':'PIP', 'mean':'BETA_MEAN', 'sd':'BETA_SD', 'allele1':'A1', 'allele2':'A2'}, inplace=True, errors='raise')

        # read log10bf
        log10bf = None
        with open(log_filename+'_sss', 'r') as f:
            for line in f:
                if line.startswith('- Log10-BF'):
                    log10bf = float(line.split()[-1])
                    break
        if log10bf is None:
            raise ValueError('FINEMP did not report Log10-BF')

        # add distance from center
        start = df_finemap['BP'].min()
        end = df_finemap['BP'].max()
        middle = (start+end)//2
        df_finemap['DISTANCE_FROM_CENTER'] = np.abs(df_finemap['BP'] - middle)

        # add causal set info
        df_finemap['CREDIBLE_SET'] = 0
        cred_file = None
        for m in range(num_causal_snps, 0, -1):
            if os.path.exists(cred_filename+str(m)):
                cred_file = cred_filename+str(m)
                break
        if cred_file is None:
            raise IOError('cred file not found')
        df_cred = pd.read_table(cred_file, sep=' ', usecols=(lambda c: c.startswith('cred')), comment='#')
        df_finemap.set_index('SNP', inplace=True, drop=False)
        for c_i, c in enumerate(df_cred.columns):
            df_finemap.loc[df_cred[c].dropna(), 'CREDIBLE_SET'] = c_i+1
        df_finemap.reset_index(inplace=True, drop=True)

        finemap_time = time.time()-t0
        logging.info('Done in %0.2f seconds'%(finemap_time))

        return df_finemap


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()

    # general parameters
    parser.add_argument('--method', required=True, help='Fine-mapping method (currently susie and finemap are supported)')
    parser.add_argument('--sumstats', required=True, help='Name of sumstats file')
    parser.add_argument('--chr', required=True, type=int, help='Target chromosome')
    parser.add_argument('--start', required=True, type=int, help='First base-pair in the region to finemap')
    parser.add_argument('--end', required=True, type=int, help='Last base-pair in the region to finemap')
    parser.add_argument('--n', required=True, type=int, help='Sample size')
    parser.add_argument(
        "--geno",
        default=None,
        help="Genotypes file (plink1 format, plink2 format(pgen) or bgen format)",
    )
    parser.add_argument('--ld', default=None, help='prefix or fill name of an LD matrix file')
    parser.add_argument('--out', required=True, help='name of the output file')
    parser.add_argument('--verbose', action='store_true', default=False, help='If specified, show verbose output')
    parser.add_argument('--debug-dir', default=None, help='If specified, this is a path of a directory that will include files for debugging problems')
    parser.add_argument('--sample-file', default=None, help='BGEN files must be used together with a sample file')
    parser.add_argument('--incl-samples', default=None, help='A single-column text file specifying the ids of individuals to include in fine-mapping')

    # fine-mapping parameters
    parser.add_argument('--max-num-causal', required=True, type=int, help='Number of causal SNPs')
    parser.add_argument('--non-funct', action='store_true', default=False, help='Perform non-functionally informed fine-mapping')
    parser.add_argument('--allow-missing', default=False, action='store_true', help='If specified, SNPs with sumstats that are not \
                            found in the LD panel will be omitted. This is not recommended, because the omitted SNPs may be causal,\
                            which could lead to false positive results')
    parser.add_argument('--allow-swapped-indel-alleles', default=False, action='store_true',
                        help='If specified, indels whose alleles are swapped between the sumstats and LD matrix \
                            are kept for fine-mapping. The default behavior considers indels at the same position \
                            with swapped alleles to be different variants, and thus removes them. Use with caution. \
                            This is intended for use only when you are confident that the indels are identical, \
                            e.g. when using insample LD')
    parser.add_argument('--no-sort-pip', default=False, action='store_true',
                        help='Do **not** sort results by PIP. Recommended for use with --susie-outfile')

    # LDstore related parameters
    parser.add_argument('--ldstore2', default=None, help='Path to an LDstore 2.0 executable file')
    parser.add_argument('--finemap-exe', default=None, help='Path to FINEMAP v1.4 executable file')
    parser.add_argument('--memory', type=int, default=1, help='Maximum amount of memory in GB to allocate to LDStore')
    parser.add_argument('--threads', type=int, default=None, help='The number of CPU cores LDstore will use (if not specified, LDstore will use the max number of CPU cores available')
    parser.add_argument('--cache-dir', default=None, help='If specified, this is a path of a directory that will cache LD matrices that have already been computed')
    parser.add_argument(
        "--cache-format",
        default=None,
        help="Format of the LDstore cache files (default: bcor); npz will with a metafile containing snpname named like .npz.gz should have columns: A1 (effect allele) A2 BP, SNP, chromosome",
        choices=["bcor", "npz", None],
    )
    # FINEMAP-specific parameters
    parser.add_argument('--finemap-dir', default=None, help='If specified, the FINEMAP files will be saved to this directory')

    # SuSiE-specific parameters
    parser.add_argument('--susie-outfile', default=None, help='If specified, the SuSiE object will be saved to an output file (also see --no-sort-pip for help merging with main output file)')
    parser.add_argument('--susie-resvar', default=None, type=float, help='If specified, SuSiE will use this value of the residual variance')
    parser.add_argument('--susie-resvar-init', default=None, type=float, help='If specified, SuSiE will use this initial value of the residual variance')
    parser.add_argument('--susie-resvar-hess', default=False, action='store_true', help='If specified, SuSiE will specify the residual variance using the HESS estimate')
    parser.add_argument('--susie-max-iter', default=100, type=int, help='SuSiE argument max_iter which controls the max number of IBSS iterations to perform (default: 100)')
    parser.add_argument('--hess', action='store_true', default=False, help='If specified, estimate causal effect variance via HESS')
    parser.add_argument('--hess-iter', type=int, default=100, help='Average HESS over this number of iterations (default: 100)')
    parser.add_argument('--hess-min-h2', type=float, default=None, help='When estimating causal effect variance via HESS, exclude SNPs that tag less than this amount of heritability (default: None)')

    # twas fine-mapping parameters
    # comming soon

    # check package versions
    check_package_versions()

    # show splash screen
    splash_screen()

    # extract args
    args = parser.parse_args()

    # check that the output directory exists
    if len(os.path.dirname(args.out))>0 and not os.path.exists(os.path.dirname(args.out)):
        raise ValueError('output directory %s doesn\'t exist'%(os.path.dirname(args.out)))    
    if args.susie_outfile is not None and len(os.path.dirname(args.susie_outfile))>0 and not os.path.exists(os.path.dirname(args.susie_outfile)):
        raise ValueError('output directory %s doesn\'t exist'%(os.path.dirname(args.susie_outfile)))    

    # configure logger
    configure_logger(args.out)

    # check params
    if args.max_num_causal==1:
        if args.geno is not None or args.ld is not None:
            raise ValueError('When max_num_causal=1 fine-mapping please omit the flags --geno and --ld (we cannot use LD information in that setting)')
    else:
        if args.geno is None:
            if args.ld is None:
                raise ValueError('must specify either --geno or --ld')
            if args.ldstore2 is not None:
                raise ValueError('cannot specify both --ld and --ldstore2')
        if args.geno is not None:
            if args.ld is not None:
                raise ValueError('cannot specify both --geno and --ld')
            if args.geno.endswith('.bgen') and args.ldstore2 is None:
                raise ValueError('You must specify --ldstore2 when --geno that points to a bgen file')

    if args.susie_outfile is not None and not args.no_sort_pip:
        logging.warning('--susie-outfile was set but not --no-sort-pip. This will make it difficult to assign SNP names to the SuSiE R object')

    # Create a fine-mapping class member
    if args.method == 'susie':
        if args.finemap_dir is not None:
            raise ValueError('--finemap-dir cannot be specified with susie method')
        finemap_obj = SUSIE_Wrapper(
            genotypes_file=args.geno,
            sumstats_file=args.sumstats,
            n=args.n,
            chr_num=args.chr,
            sample_file=args.sample_file,
            incl_samples=args.incl_samples,
            ldstore_exe=args.ldstore2,
            n_threads=args.threads,
            cache_dir=args.cache_dir,
            cache_format=args.cache_format,
            memory=args.memory,
            allow_swapped_indel_alleles=args.allow_swapped_indel_alleles,
        )
    elif args.method == 'finemap':
        if args.susie_outfile is not None:
            raise ValueError('--susie-outfile cannot be specified with finemap method')
        if args.finemap_exe is None:
            raise ValueError('need to specify --finemap-exe')
        if args.hess:
            raise ValueError('FINEMAP cannot be used with --hess')
        finemap_obj = FINEMAP_Wrapper(
            genotypes_file=args.geno,
            sumstats_file=args.sumstats,
            n=args.n,
            chr_num=args.chr,
            sample_file=args.sample_file,
            incl_samples=args.incl_samples,
            ldstore_exe=args.ldstore2,
            finemap_exe=args.finemap_exe,
            n_threads=args.threads,
            cache_dir=args.cache_dir,
            cache_format=args.cache_format,
            memory=args.memory,
            allow_swapped_indel_alleles=args.allow_swapped_indel_alleles,
        )
    else:
        raise ValueError('unknown method specified in --method')

    # run fine-mapping
    df_finemap = finemap_obj.finemap(locus_start=args.start, locus_end=args.end, num_causal_snps=args.max_num_causal,
                 use_prior_causal_prob=not args.non_funct, prior_var=None,
                 hess=args.hess, hess_iter=args.hess_iter, hess_min_h2=args.hess_min_h2,
                 verbose=args.verbose, ld_file=args.ld, debug_dir=args.debug_dir, allow_missing=args.allow_missing,
                 susie_outfile=args.susie_outfile, finemap_dir=args.finemap_dir,
                 residual_var=args.susie_resvar, residual_var_init=args.susie_resvar_init, hess_resvar=args.susie_resvar_hess,
                 susie_max_iter=args.susie_max_iter)
    logging.info('Writing fine-mapping results to %s'%(args.out))
    if not args.no_sort_pip:
        df_finemap.sort_values('PIP', ascending=False, inplace=True)
    if args.out.endswith('.parquet'):
        df_finemap.to_parquet(args.out, index=False)
    else:
        df_finemap.to_csv(args.out, sep="\t", index=False, float_format="%0.5e")
