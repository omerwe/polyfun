import pandas as pd
import numpy as np
import os
import time
import scipy.stats as stats
import logging
from pandas.api.types import is_integer_dtype
from polyfun_utils import check_package_versions, configure_logger


def compute_Neff(df_sumstats, n, chi2_cutoff=30):
    df_sumstats_chi2_large = df_sumstats.query('CHISQ_BOLT_LMM > %s'%(chi2_cutoff))
    if df_sumstats_chi2_large.shape[0]==0:
        return n
    Neff = int(np.median(df_sumstats_chi2_large['CHISQ_BOLT_LMM'] / df_sumstats_chi2_large['CHISQ_LINREG']) * n)
    return Neff
    

def find_df_column(df, strings_to_find, allow_missing=False):
    
    if isinstance(strings_to_find, str):
        strings_to_find = [strings_to_find]
        
    is_relevant_col = np.zeros(df.shape[1], dtype=np.bool)
    for str_to_find in strings_to_find:
        is_relevant_col = is_relevant_col | (df.columns.str.upper() == str_to_find.upper())
    if np.sum(is_relevant_col)==0:
        if allow_missing:
            return ''
        else:
            raise ValueError('No matching column found among: %s'%str(strings_to_find))
    elif np.sum(is_relevant_col)>1:
        raise ValueError('Too many matching columns found among: %s'%str(strings_to_find))
    else:
        return df.columns[is_relevant_col][0]
        
        
def rename_df_columns(df_sumstats):
    chr_column = find_df_column(df_sumstats, ['CHR', 'CHROMOSOME', 'CHROM'])
    bp_column = find_df_column(df_sumstats, ['BP', 'POS', 'POSITION', 'COORDINATE', 'BASEPAIR'])
    snp_column = find_df_column(df_sumstats, ['SNP', 'RSID', 'RS', 'NAME'])
    a1freq_col = find_df_column(df_sumstats, ['A1FREQ', 'freq', 'MAF', 'FRQ'], allow_missing=True)
    info_col = find_df_column(df_sumstats, 'INFO', allow_missing=True)
    beta_col = find_df_column(df_sumstats, ['BETA', 'EFF', 'EFFECT', 'EFFECT_SIZE', 'OR'], allow_missing=True)
    se_col = find_df_column(df_sumstats, ['SE'], allow_missing=True)
    pvalue_col = find_df_column(df_sumstats, ['P_BOLT_LMM', 'P', 'PVALUE', 'P-VALUE', 'P_value', 'PVAL'], allow_missing=True)
    z_col = find_df_column(df_sumstats, ['Z', 'ZSCORE', 'Z_SCORE'], allow_missing=True)    
    n_col = find_df_column(df_sumstats, ['N', 'sample_size'], allow_missing=True)    
    ncase_col = find_df_column(df_sumstats, ['N_cases', 'Ncase'], allow_missing=True)    
    ncontrol_col = find_df_column(df_sumstats, ['N_controls', 'Ncontrol'], allow_missing=True)    
    try:
        allele1_col = find_df_column(df_sumstats, ['ALLELE1', 'A1'])
        allele0_col = find_df_column(df_sumstats, ['ALLELE0', 'A0'])
    except ValueError:
        allele1_col = find_df_column(df_sumstats, ['ALLELE1', 'A1'])
        allele0_col = find_df_column(df_sumstats, ['ALLELE2', 'A2'])
    
    return df_sumstats.rename(columns={snp_column:'SNP', allele1_col:'A1', 
              allele0_col:'A2', a1freq_col:'MAF', bp_column:'BP', 
             chr_column:'CHR', info_col:'INFO', beta_col:'BETA', 
             se_col:'SE', pvalue_col:'P', z_col:'Z', n_col:'N',
             ncase_col:'N_CASES', ncontrol_col:'N_CONTROLS'}, errors='ignore')



def compute_z(df_sumstats):

    #make sure that we have required fields
    if 'BETA' not in df_sumstats.columns:
        raise ValueError('Beta column not found in sumstats file (required to compute Z-scores)')
    if 'P' not in df_sumstats.columns:
        raise ValueError('P-value column not found in sumstats file (required to compute Z-scores)')

    #compute z-scores
    df_sumstats['Z'] = stats.norm(0,1).isf(df_sumstats['P'] / 2.0) * np.sign(df_sumstats['BETA'])
    
    #Use LDpred-funct trick to estimate Z for SNPs with P=0
    is_zero_pval = np.isinf(df_sumstats['Z'])    
    if np.any(is_zero_pval):
    
        #make sure that we have required fields
        if 'MAF' not in df_sumstats.columns:
            raise ValueError('MAF column not found in sumstats file (required to compute Z-scores for SNPs with P=0)')
        
        #estimate sigma2pheno
        df_sumstats_nonzero = df_sumstats.loc[~is_zero_pval]
        df_snp_var_nonzero = 2 * df_sumstats_nonzero['MAF'] * (1-df_sumstats_nonzero['MAF'])
        z_prop = df_sumstats_nonzero['BETA'] * np.sqrt(df_snp_var_nonzero)
        assert np.corrcoef(z_prop.values, df_sumstats_nonzero['Z'].values)[0,1] > 0.6
        sqrt_sigma2pheno = np.median(df_sumstats_nonzero['Z'].values / z_prop)
        assert not np.isnan(sqrt_sigma2pheno)
    
        #compute Z for SNPs with P=0
        df_sumstats_iszero = df_sumstats.loc[is_zero_pval]
        df_snp_var_zero = 2 * df_sumstats_iszero['MAF'] * (1-df_sumstats_iszero['MAF'])
        df_sumstats.loc[is_zero_pval, 'Z'] = df_sumstats_iszero['BETA'] * np.sqrt(df_snp_var_zero) * sqrt_sigma2pheno
        assert df_sumstats.loc[is_zero_pval, 'Z'].notnull().all()
    
    return df_sumstats
        
            


def filter_sumstats(df_sumstats, min_info_score=None, min_maf=None, remove_strand_ambig=False, keep_hla=False):

    logging.info('%d SNPs are in the sumstats file'%(df_sumstats.shape[0]))
    is_good_snp = np.ones(df_sumstats.shape[0], dtype=np.bool)
    
    #remove 'bad' BOLT-LMM SNPs
    if 'CHISQ_BOLT_LMM' in df_sumstats.columns:
        is_good_chi2_snp = df_sumstats['CHISQ_BOLT_LMM']>0
        is_good_snp = is_good_snp & is_good_chi2_snp
        if np.any(~is_good_chi2_snp):
            logging.info('Removing %d SNPs with BOLT CHI2=0'%(np.sum(~is_good_chi2_snp)))

    #Filter SNPs based on INFO score    
    if min_info_score is not None and min_info_score>0:        
        if 'INFO' not in df_sumstats.columns:
            logging.warning('Could not find INFO column. Please set --min-info 0 to omit this warning.')
        else:
            is_good_info_snp = (df_sumstats['INFO'] >= min_info_score)
            is_good_snp = is_good_snp & is_good_info_snp
            if np.any(~is_good_info_snp):
                logging.info('Removing %d SNPs with INFO<%0.2f'%(np.sum(~is_good_info_snp), min_info_score))
    
    #filter SNPs based on MAF
    if min_maf is not None and min_maf>0:        
        if 'MAF' not in df_sumstats.columns:
            logging.warning('Could not find MAF column. Please set --min-maf 0 to omit this warning.')
        else:
            is_good_maf_snp = (df_sumstats['MAF'].between(min_maf, 1-min_maf))
            is_good_snp = is_good_snp & is_good_maf_snp
            if np.any(~is_good_maf_snp):
                logging.info('Removing %d SNPs with MAF<%s'%(np.sum(~is_good_maf_snp), min_maf))
        

    #find strand ambiguous summary statistics
    if remove_strand_ambig:
        is_strand_ambig = np.zeros(df_sumstats.shape[0], dtype=np.bool)
        for ambig_pairs in [('A', 'T'), ('T', 'A'), ('C', 'G'), ('G', 'C')]:                    
            is_strand_ambig = is_strand_ambig | ((df_sumstats['A2']==ambig_pairs[0]) & (df_sumstats['A1']==ambig_pairs[1]))         
        is_good_snp = is_good_snp & (~is_strand_ambig)            
        if np.any(is_strand_ambig):
            logging.info('Removing %d SNPs with strand ambiguity'%(is_strand_ambig.sum()))
            
            
    #remove HLA SNPs
    if not keep_hla:
        is_hla = (df_sumstats['CHR']==6) & (df_sumstats['BP'].between(28000000, 34000000))
        is_good_snp = is_good_snp & (~is_hla)
        if np.any(is_hla):
            logging.info('Removing %d HLA SNPs'%(is_hla.sum()))            
            

    #finally do the actual filtering
    if np.any(~is_good_snp):
        if not np.any(is_good_snp):
            raise ValueError('No SNPs remained after all filtering stages')
        df_sumstats = df_sumstats.loc[is_good_snp]        
        logging.info('%d SNPs with sumstats remained after all filtering stages'%(df_sumstats.shape[0]))
        
    return df_sumstats
    
    

def compute_casecontrol_neff(df_sumstats):
    logging.info('Computing the effective sample size for case-control data...')    
    Neff = (4.0 / (1.0/df_sumstats['N_CASES'] + 1.0/df_sumstats['N_CONTROLS'])).astype(np.int)
    return Neff


def sanity_checks(df_sumstats):
    if np.any(df_sumstats['CHR'].isnull()):
        raise ValueError('Some SNPs have a null chromosome value')
    if np.any(df_sumstats['BP'].isnull()):
        raise ValueError('Some SNPs have a null base-pair value')
    if not is_integer_dtype(df_sumstats['CHR']):
        raise ValueError('Some chromosome values are not integers. Please double-check your input')
    if not is_integer_dtype(df_sumstats['BP']):
        raise ValueError('Some base-pair values are not integers. Please double-check your input')
        
    #check for duplicates
    df_snp = df_sumstats['CHR'].astype('str') + '.' + \
             df_sumstats['BP'].astype('str') + '.' + \
             df_sumstats['A1'].astype('str') + '.' + \
             df_sumstats['A2'].astype('str')
    if np.any(df_snp.duplicated()):
        raise ValueError('The input file includes duplicate SNPs')
    
    #compute Z    
    
    
def convert_odds_ratio_to_log(df_sumstats):
    if 'OR' in df_sumstats.columns and np.all(df_sumstats['OR']>0):
        df_sumstats['OR'] = np.log(df_sumstats['OR'])
        logging.info('Converting OR column to log-odds')
    return df_sumstats
        
    
    
if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()    
    parser.add_argument('--sumstats', required=True, help='Input summary statistics file')
    parser.add_argument('--out', required=True, help='Name of output file')
    parser.add_argument('--n', type=int, default=None, help='Sample size. If not specified, will try to infer this from the input file')
    parser.add_argument('--min-info', type=float, default=0.6, help='Minimum INFO score (set to zero to avoid INFO-based filtering)')
    parser.add_argument('--min-maf', type=float, default=0.001, help='Minimum MAF (set to zero to avoid MAF-based filtering)')
    parser.add_argument('--remove-strand-ambig', default=False, action='store_true', help='If specified, strand-ambigous SNPs will be removed')
    parser.add_argument('--chi2-cutoff', type=float, default=30, help ='Chi2 cutoff for effective sample size computations')
    parser.add_argument('--keep-hla', default=False, action='store_true', help='If specified, Keep SNPs in the HLA region')
    parser.add_argument('--no-neff', default=False, action='store_true', help='If specified, use the true rather than the effective sample size in BOLT-LMM runs')
    args = parser.parse_args()
    
    #check package versions
    check_package_versions()    
    
    #configure the logger
    configure_logger(args.out)
    
    #read sumstats file
    logging.info('Reading sumstats file...')
    t0 = time.time()
    df_sumstats = pd.read_table(args.sumstats, delim_whitespace=True)
    logging.info('Done in %0.2f seconds'%(time.time()-t0))
    
    #convert odds-ratio to log-odds ratio if needed
    df_sumstats = convert_odds_ratio_to_log(df_sumstats)
    
    #rename df_sumstats columns
    df_sumstats = rename_df_columns(df_sumstats)

    #filter sumstats
    df_sumstats = filter_sumstats(df_sumstats, min_info_score=args.min_info, min_maf=args.min_maf, remove_strand_ambig=args.remove_strand_ambig, keep_hla=args.keep_hla)
    
    #do some sanity checks
    sanity_checks(df_sumstats)

    #compute Neff    
    if 'CHISQ_BOLT_LMM' in df_sumstats.columns and not args.no_neff:
        if args.n is None:
            raise ValueError('--n must be specified with BOLT input files')
        Neff = compute_Neff(df_sumstats, args.n, args.chi2_cutoff)
        logging.info('Effective sample size is %s'%(Neff))
        df_sumstats['N']  = Neff
    elif args.n is not None:
        if 'N' in df_sumstats.columns:
            raise ValueError('cannot both specify --n and have an N column in the sumstats file')
        if 'N_CASES' in df_sumstats.columns or 'N_CONTROLS' in df_sumstats.columns:
            raise ValueError('cannot both specify --n and have an N_cases/N_controls column in the sumstats file')
        df_sumstats['N']  = args.n
    elif 'N' in df_sumstats.columns:
        if 'N_CASES' in df_sumstats.columns or 'N_CONTROLS' in df_sumstats.columns:
            raise ValueError('cannot both have an N column and N_cases/N_controls columns in the sumstats file')
        pass
    elif 'N_CASES' in df_sumstats.columns and 'N_CONTROLS' in df_sumstats.columns:
        df_sumstats['N']  = compute_casecontrol_neff(df_sumstats)
    else:
        raise ValueError('must specify sample size, via either (1) --n flag; (2) N column in the sumstats file; or (3) Two columns N_cases, N_controls in the sumstats file')
    
    # #create SNP string
    # df_sumstats['SNP'] = df_sumstats['SNP'].astype('str') + '.' + \
                         # df_sumstats['CHR'].astype('str') + '.' + \
                         # df_sumstats['BP'].astype('str') + '.' + \
                         # df_sumstats['A1'].astype('str') + '.' + \
                         # df_sumstats['A2'].astype('str')
    
    #compute Z
    if 'Z' in df_sumstats.columns:
        pass
    elif 'CHISQ_BOLT_LMM' in df_sumstats.columns:
        df_sumstats['Z'] = np.sqrt(df_sumstats['CHISQ_BOLT_LMM']) * np.sign(df_sumstats['BETA'])
    elif 'BETA' in df_sumstats.columns and 'SE' in df_sumstats.columns:
        if np.any(df_sumstats['SE']==0):
            raise ValueError('Found SNPs with BETA stderr==0')
        df_sumstats['Z'] = df_sumstats['BETA'] / df_sumstats['SE']
    elif 'P' in df_sumstats.columns:
        df_sumstats = compute_z(df_sumstats)
    else:
        raise ValueError('Sumstats file must include a p-value, Z-score or chi2 column to compute Z-scores')
    
    #write output
    logging.info('Saving munged sumstats of %d SNPs to %s'%(df_sumstats.shape[0], args.out))
    df_sumstats[['SNP', 'CHR', 'BP', 'A1', 'A2', 'Z', 'N']].to_parquet(args.out)
    logging.info('Done')
    
    
