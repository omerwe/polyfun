'''
(c) 2014 Brendan Bulik-Sullivan and Hilary Finucane

This module deals with getting all the data needed for LD Score regression from files
into memory and checking that the input makes sense. There is no math here. LD Score
regression is implemented in the regressions module.

'''

import numpy as np
import pandas as pd
from scipy import stats
import itertools as it
from . import parse as ps
from . import regressions as reg
import sys
import traceback
import copy
import os


_N_CHR = 22
# complementary bases
COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
# bases
BASES = list(COMPLEMENT.keys())
# true iff strand ambiguous
STRAND_AMBIGUOUS = {''.join(x): x[0] == COMPLEMENT[x[1]]
                    for x in it.product(BASES, BASES)
                    if x[0] != x[1]}
# SNPS we want to keep (pairs of alleles)
VALID_SNPS = {x for x in [''.join(y) for y in it.product(BASES, BASES)]
              if x[0] != x[1] and not STRAND_AMBIGUOUS[x]}
# T iff SNP 1 has the same alleles as SNP 2 (allowing for strand or ref allele flip).
MATCH_ALLELES = {x for x in [''.join(y) for y in it.product(VALID_SNPS, VALID_SNPS)]
                 # strand and ref match
                 if ((x[0] == x[2]) and (x[1] == x[3])) or
                 # ref match, strand flip
                 ((x[0] == COMPLEMENT[x[2]]) and (x[1] == COMPLEMENT[x[3]])) or
                 # ref flip, strand match
                 ((x[0] == x[3]) and (x[1] == x[2])) or
                 ((x[0] == COMPLEMENT[x[3]]) and (x[1] == COMPLEMENT[x[2]]))}  # strand and ref flip
# T iff SNP 1 has the same alleles as SNP 2 w/ ref allele flip.
FLIP_ALLELES = {''.join(x):
                ((x[0] == x[3]) and (x[1] == x[2])) or  # strand match
                # strand flip
                ((x[0] == COMPLEMENT[x[3]]) and (x[1] == COMPLEMENT[x[2]]))
                for x in MATCH_ALLELES}


def _splitp(fstr):
    flist = fstr.split(',')
    flist = [os.path.expanduser(os.path.expandvars(x)) for x in flist]
    return flist


def _select_and_log(x, ii, log, msg):
    '''Fiter down to rows that are True in ii. Log # of SNPs removed.'''
    new_len = ii.sum()
    if new_len == 0:
        raise ValueError(msg.format(N=0))
    else:
        x = x[ii]
        log.log(msg.format(N=new_len))
    return x


def smart_merge(x, y):
    '''Check if SNP columns are equal. If so, save time by using concat instead of merge.'''
    assert len(set(x.columns).intersection(set(y.drop(columns=['SNP']).columns))) == 0
    if len(x) == len(y) and (x.index == y.index).all() and (x.SNP == y.SNP).all():
        #x = x.reset_index(drop=True)
        #y = y.reset_index(drop=True).drop(columns=['SNP'])
        #out = pd.concat([x, y], axis=1)        
        out = pd.concat([x, y.drop(columns=['SNP'])], axis=1)
    else:
        if x.index.name == 'snpid' and y.index.name == 'snpid':
            out = pd.merge(x, y.drop(columns=['SNP']), how='inner', left_index=True, right_index=True)
        else:
            out = pd.merge(x, y, how='inner', on='SNP')
    return out


def _read_ref_ld(args, log):
    '''Read reference LD Scores.'''
    ref_ld = _read_chr_split_files(args.ref_ld_chr, args.ref_ld, log,
                                   'reference panel LD Score', ps.ldscore_fromlist)
    log.log(
        'Read reference panel LD Scores for {N} SNPs.'.format(N=len(ref_ld)))
    return ref_ld


def _read_annot(args, log):
    '''Read annot matrix.'''
    if (args.anno is not None): annotations = args.anno.split(',')
    else: annotations = None
    try:
        if args.ref_ld is not None:
            overlap_matrix, M_tot = _read_chr_split_files(args.ref_ld_chr, args.ref_ld, log,
                                                          'annot matrix', ps.annot, frqfile=args.frqfile, anno=annotations)
        elif args.ref_ld_chr is not None:
            overlap_matrix, M_tot = _read_chr_split_files(args.ref_ld_chr, args.ref_ld, log,
                                                      'annot matrix', ps.annot, frqfile=args.frqfile_chr, anno=annotations)
    except Exception:
        log.log('Error parsing .annot file.')
        raise

    return overlap_matrix, M_tot


def _read_M(args, log, n_annot):
    '''Read M (--M, --M-file, etc).'''
    if args.M:
        try:
            M_annot = [float(x) for x in _splitp(args.M)]
        except ValueError as e:
            raise ValueError('Could not cast --M to float: ' + str(e.args))
    else:
        if args.ref_ld:
            M_annot = ps.M_fromlist(
                _splitp(args.ref_ld), common=(not args.not_M_5_50))
        elif args.ref_ld_chr:
            M_annot = ps.M_fromlist(
                _splitp(args.ref_ld_chr), _N_CHR, common=(not args.not_M_5_50))

    try:
        M_annot = np.array(M_annot).reshape((1, n_annot))
    except ValueError as e:
        raise ValueError(
            '# terms in --M must match # of LD Scores in --ref-ld.\n' + str(e.args))

    return M_annot


def _read_w_ld(args, log):
    '''Read regression SNP LD.'''
    if (args.w_ld and ',' in args.w_ld) or (args.w_ld_chr and ',' in args.w_ld_chr):
        raise ValueError(
            '--w-ld must point to a single fileset (no commas allowed).')
    w_ld = _read_chr_split_files(args.w_ld_chr, args.w_ld, log,
                                 'regression weight LD Score', ps.ldscore_fromlist)
    w_ld.drop(['CHR'], axis=1, inplace=True)
    if len(w_ld.columns) != 2:
        raise ValueError('--w-ld may only have one LD Score column.')
    w_ld.columns = ['SNP', 'LD_weights']  # prevent colname conflicts w/ ref ld
    log.log(
        'Read regression weight LD Scores for {N} SNPs.'.format(N=len(w_ld)))
    return w_ld


def _read_chr_split_files(chr_arg, not_chr_arg, log, noun, parsefunc, **kwargs):
    '''Read files split across 22 chromosomes (annot, ref_ld, w_ld).'''
    try:
        if not_chr_arg:
            log.log('Reading {N} from {F} ...'.format(F=not_chr_arg, N=noun))
            out = parsefunc(_splitp(not_chr_arg), **kwargs)
        elif chr_arg:
            f = ps.sub_chr(chr_arg, '[1-22]')
            log.log('Reading {N} from {F} ...'.format(F=f, N=noun))
            out = parsefunc(_splitp(chr_arg), _N_CHR, **kwargs)
    except ValueError as e:
        log.log('Error parsing {N}.'.format(N=noun))
        raise e

    return out


def _read_sumstats(args, log, fh, alleles=True, dropna=False):
    '''Parse summary statistics.'''
    log.log('Reading summary statistics from {S} ...'.format(S=fh))
    sumstats = ps.sumstats(fh, alleles=alleles, dropna=dropna)
    log_msg = 'Read summary statistics for {N} SNPs.'
    log.log(log_msg.format(N=len(sumstats)))
    if np.any(sumstats.index.duplicated()):
        m = len(sumstats)
        sumstats = sumstats.loc[~sumstats.index.duplicated()]
        log.log('Dropped {M} duplicated SNPs.'.format(M=m - len(sumstats)))

    return sumstats


def _check_ld_condnum(args, log, ref_ld):
    '''Check condition number of LD Score matrix.'''
    if len(ref_ld.shape) >= 2:
        cond_num = int(np.linalg.cond(ref_ld))
        if cond_num > 100000:
            if args.invert_anyway:
                warn = "WARNING: LD Score matrix condition number is {C}. "
                warn += "Inverting anyway because the --invert-anyway flag is set."
                log.log(warn.format(C=cond_num))
            else:
                warn = "WARNING: LD Score matrix condition number is {C}. "
                warn += "Remove collinear LD Scores. "
                raise ValueError(warn.format(C=cond_num))


def _check_variance(log, M_annot, ref_ld):
    '''Remove zero-variance LD Scores.'''
    ###ii = ref_ld.iloc[:, 2:].var(axis=0) == 0  # NB there is a SNP and CHR column here
    ii = np.array([(ref_ld[c].var() == 0) for c in ref_ld.columns[2:]]) #This command uses way way less memory
    if ii.all():
        raise ValueError('All LD Scores have zero variance.')
    elif ii.any():
        log.log('Removing partitioned LD Scores with zero variance: %s'%(','.join(ref_ld.columns[2:][ii])))
        ii_snp = np.array([True, True] + list(~ii))
        ii_m = np.array(~ii)
        ref_ld = ref_ld.loc[:, ii_snp]
        M_annot = M_annot[:, ii_m]

    return M_annot, ref_ld, ii


def _warn_length(log, sumstats):
    if len(sumstats) < 200000:
        log.log(
            'WARNING: number of SNPs less than 200k; this is almost always bad.')


def _print_cov(ldscore_reg, ofh, log):
    '''Prints covariance matrix of slopes.'''
    log.log(
        'Printing covariance matrix of the estimates to {F}.'.format(F=ofh))
    np.savetxt(ofh, ldscore_reg.coef_cov)


def _print_delete_values(ldscore_reg, ofh, log):
    '''Prints block jackknife delete-k values'''
    log.log('Printing block jackknife delete values to {F}.'.format(F=ofh))
    np.savetxt(ofh, ldscore_reg.tot_delete_values)

def _print_part_delete_values(ldscore_reg, ofh, log):
    '''Prints partitioned block jackknife delete-k values'''
    log.log('Printing partitioned block jackknife delete values to {F}.'.format(F=ofh))
    np.savetxt(ofh, ldscore_reg.part_delete_values)


def _merge_and_log(ld, sumstats, noun, log):
    '''Wrap smart merge with log messages about # of SNPs.'''
    sumstats = smart_merge(ld, sumstats)
    msg = 'After merging with {F}, {N} SNPs remain.'
    if len(sumstats) == 0:
        msg += ' Please make sure that your annotation files include the SNPs in your sumstats files (please see the PolyFun wiki for details on downloading functional annotations)'
        raise ValueError(msg.format(N=len(sumstats), F=noun))
    else:
        log.log(msg.format(N=len(sumstats), F=noun))

    return sumstats


def _read_ld_sumstats(args, log, fh, alleles=True, dropna=True):
    sumstats = _read_sumstats(args, log, fh, alleles=alleles, dropna=dropna)
    ref_ld = _read_ref_ld(args, log)
    n_annot = len(ref_ld.columns) - 2 #Changed to -2 because we also have chromosome column now
    M_annot = _read_M(args, log, n_annot)
    
    #keep only requested annotations if --anno was specified
    if args.anno is not None:
        cols_to_keep = np.zeros(len(ref_ld.columns), dtype=np.bool)        
        annotations = args.anno.split(',')
        if np.any(~np.isin(annotations, ref_ld.columns.str[:-2])):
            raise ValueError('Not all annotations specified with --anno are found in the LD scores file')        
        cols_to_keep = (ref_ld.columns.str[:-2].isin(annotations)) | (ref_ld.columns.str[:-4].isin(annotations)) | (ref_ld.columns.isin(['CHR', 'SNP']))
        assert np.sum(cols_to_keep) == len(annotations)+2
        cols_nochrsnp = ref_ld.drop(columns=['CHR', 'SNP']).columns
        M_cols_to_keep = (cols_nochrsnp.str[:-2].isin(annotations)) | (cols_nochrsnp.str[:-4].isin(annotations))
        assert np.sum(M_cols_to_keep) == len(annotations)
        ref_ld = ref_ld.loc[:, cols_to_keep]
        M_annot = M_annot[:, M_cols_to_keep]
        log.log('Keeping only annotations specified with --anno')
    
    M_annot, ref_ld, novar_cols = _check_variance(log, M_annot, ref_ld)
    w_ld = _read_w_ld(args, log)
    sumstats = _merge_and_log(ref_ld, sumstats, 'reference panel LD', log)
    sumstats = _merge_and_log(sumstats, w_ld, 'regression SNP LD', log)
    w_ld_cname = sumstats.columns[-1]
    ref_ld_cnames = ref_ld.drop(columns=['CHR', 'SNP']).columns
    return M_annot, w_ld_cname, ref_ld_cnames, sumstats, novar_cols

def estimate_h2(args, log):
    '''Estimate h2 and partitioned h2.'''
    args = copy.deepcopy(args)
    if args.samp_prev is not None and args.pop_prev is not None:
        args.samp_prev, args.pop_prev = list(map(
            float, [args.samp_prev, args.pop_prev]))
    if args.intercept_h2 is not None:
        args.intercept_h2 = float(args.intercept_h2)
    if args.no_intercept:
        args.intercept_h2 = 1
    M_annot, w_ld_cname, ref_ld_cnames, sumstats, novar_cols = _read_ld_sumstats(
        args, log, args.h2)
    ref_ld = np.array(sumstats[ref_ld_cnames], dtype=np.float32)
    if not args.skip_cond_check:
        _check_ld_condnum(args, log, ref_ld_cnames)
    _warn_length(log, sumstats)
    n_snp = len(sumstats)
    n_blocks = min(n_snp, args.n_blocks)
    n_annot = len(ref_ld_cnames)
    chisq_max = args.chisq_max
    old_weights = False
    if n_annot == 1:
        if args.two_step is None and args.intercept_h2 is None:
            args.two_step = 30
    else:
        old_weights = True
        if args.chisq_max is None:
            chisq_max = max(0.001*sumstats.N.max(), args.max_chi2)

    s = lambda x: np.array(x).reshape((n_snp, 1))
    chisq = s(sumstats.Z**2).astype(np.float32)
    if chisq_max is not None and not args.keep_large:
        ii = np.ravel(chisq < chisq_max)
        sumstats = sumstats.loc[ii, :]
        log.log('Removed {M} SNPs with chi^2 > {C} ({N} SNPs remain)'.format(
                C=chisq_max, N=np.sum(ii), M=n_snp-np.sum(ii)))
        n_snp = np.sum(ii)  # lambdas are late-binding, so this works
        ref_ld = np.array(sumstats[ref_ld_cnames], dtype=np.float32)
        chisq = chisq[ii].reshape((n_snp, 1))

    if args.two_step is not None:
        log.log('Using two-step estimator with cutoff at {M}.'.format(M=args.two_step))

    hsqhat = reg.Hsq(chisq, ref_ld, s(sumstats[w_ld_cname]), s(sumstats.N),
                    M_annot, n_blocks=n_blocks, intercept=args.intercept_h2,
                    twostep=args.two_step, old_weights=old_weights,
                    chr_num=sumstats['CHR'],
                    loco=args.loco, ridge_lambda=args.ridge_lambda,
                    standardize_ridge=not args.no_standardize_ridge,
                    approx_ridge=not args.reestimate_lambdas,
                    skip_ridge_jackknife=not args.ridge_jackknife,
                    num_chr_sets = args.num_chr_sets,
                    evenodd_split=args.evenodd_split,
                    nn=args.nn,
                    keep_large=args.keep_large,
                    nnls_exact=args.nnls_exact
                    )

    if args.print_cov:
        _print_cov(hsqhat, args.out + '.cov', log)
    if args.print_delete_vals:
        _print_delete_values(hsqhat, args.out + '.delete', log)
        _print_part_delete_values(hsqhat, args.out + '.part_delete', log)
        
        
    #save ridge-regression lambdas if possible
    if args.loco and args.ridge_jackknife and args.reestimate_lambdas:
        np.savetxt(args.out+'.out_of_chrom_r2.txt', [hsqhat.jknife_ridge.r2_best_lambda])
        df = pd.Series(hsqhat.jknife_ridge.best_r2_jk_noblock)
        df.to_csv(args.out+'.out_of_chrom_r2_jk.txt', index=False, header=False)

    log.log(hsqhat.summary(ref_ld_cnames, P=args.samp_prev, K=args.pop_prev, overlap = args.overlap_annot))
    if args.overlap_annot:
        overlap_matrix, M_tot = _read_annot(args, log)

        # overlap_matrix = overlap_matrix[np.array(~novar_cols), np.array(~novar_cols)]#np.logical_not
        df_results = hsqhat._overlap_output(ref_ld_cnames, overlap_matrix, M_annot, M_tot, args.print_coefficients)
        df_results.to_csv(args.out+'.results', sep="\t", index=False, na_rep='NA', float_format='%0.4e')
        log.log('Results printed to '+args.out+'.results')

    return hsqhat


def estimate_rg(args, log):
    '''Estimate rg between trait 1 and a list of other traits.'''
    args = copy.deepcopy(args)
    rg_paths, rg_files = _parse_rg(args.rg)
    n_pheno = len(rg_paths)
    f = lambda x: _split_or_none(x, n_pheno)
    args.intercept_h2, args.intercept_gencov, args.samp_prev, args.pop_prev = list(map(f,
        (args.intercept_h2, args.intercept_gencov, args.samp_prev, args.pop_prev)))
    list(map(lambda x: _check_arg_len(x, n_pheno), ((args.intercept_h2, '--intercept-h2'),
                                               (args.intercept_gencov, '--intercept-gencov'),
                                               (args.samp_prev, '--samp-prev'),
                                               (args.pop_prev, '--pop-prev'))))
    if args.no_intercept:
        args.intercept_h2 = [1 for _ in range(n_pheno)]
        args.intercept_gencov = [0 for _ in range(n_pheno)]
    p1 = rg_paths[0]
    out_prefix = args.out + rg_files[0]
    M_annot, w_ld_cname, ref_ld_cnames, sumstats, _ = _read_ld_sumstats(args, log, p1,
                                                                        alleles=True, dropna=True)
    RG = []
    n_annot = M_annot.shape[1]
    if n_annot == 1 and args.two_step is None and args.intercept_h2 is None:
        args.two_step = 30
    if args.two_step is not None:
        log.log('Using two-step estimator with cutoff at {M}.'.format(M=args.two_step))

    for i, p2 in enumerate(rg_paths[1:n_pheno]):
        log.log(
            'Computing rg for phenotype {I}/{N}'.format(I=i + 2, N=len(rg_paths)))
        try:
            loop = _read_other_sumstats(args, log, p2, sumstats, ref_ld_cnames)
            rghat = _rg(loop, args, log, M_annot, ref_ld_cnames, w_ld_cname, i)
            RG.append(rghat)
            _print_gencor(args, log, rghat, ref_ld_cnames, i, rg_paths, i == 0)
            out_prefix_loop = out_prefix + '_' + rg_files[i + 1]
            if args.print_cov:
                _print_rg_cov(rghat, out_prefix_loop, log)
            if args.print_delete_vals:
                _print_rg_delete_values(rghat, out_prefix_loop, log)

        except Exception:  # keep going if phenotype 50/100 causes an error
            msg = 'ERROR computing rg for phenotype {I}/{N}, from file {F}.'
            log.log(msg.format(I=i + 2, N=len(rg_paths), F=rg_paths[i + 1]))
            ex_type, ex, tb = sys.exc_info()
            log.log(traceback.format_exc(ex) + '\n')
            if len(RG) <= i:  # if exception raised before appending to RG
                RG.append(None)

    log.log('\nSummary of Genetic Correlation Results\n' +
            _get_rg_table(rg_paths, RG, args))
    return RG


def _read_other_sumstats(args, log, p2, sumstats, ref_ld_cnames):
    loop = _read_sumstats(args, log, p2, alleles=True, dropna=False)
    loop = _merge_sumstats_sumstats(args, sumstats, loop, log)
    loop = loop.dropna(how='any')
    alleles = loop.A1 + loop.A2 + loop.A1x + loop.A2x
    if not args.no_check_alleles:
        loop = _select_and_log(loop, _filter_alleles(alleles), log,
                               '{N} SNPs with valid alleles.')

    loop['Z2'] = _align_alleles(loop.Z2, alleles)    
    loop = loop.drop(['A1', 'A1x', 'A2', 'A2x'], axis=1)
    _check_ld_condnum(args, log, loop[ref_ld_cnames])
    _warn_length(log, loop)
    return loop


def _get_rg_table(rg_paths, RG, args):
    '''Print a table of genetic correlations.'''
    t = lambda attr: lambda obj: getattr(obj, attr, 'NA')
    x = pd.DataFrame()
    x['p1'] = [rg_paths[0] for i in range(1, len(rg_paths))]
    x['p2'] = rg_paths[1:len(rg_paths)]
    x['rg'] = list(map(t('rg_ratio'), RG))
    x['se'] = list(map(t('rg_se'), RG))
    x['z'] = list(map(t('z'), RG))
    x['p'] = list(map(t('p'), RG))
    if args.samp_prev is not None and args.pop_prev is not None and\
            all((i is not None for i in args.samp_prev)) and all((i is not None for it in args.pop_prev)):
        c = reg.h2_obs_to_liab(1, args.samp_prev[1], args.pop_prev[1])
        x['h2_liab'] = [c * x for x in list(map(t('tot'), list(map(t('hsq2'), RG))))]
        x['h2_liab_se'] = [c * x for x in list(map(t('tot_se'), list(map(t('hsq2'), RG))))]
    else:
        x['h2_obs'] = list(map(t('tot'), list(map(t('hsq2'), RG))))
        x['h2_obs_se'] = list(map(t('tot_se'), list(map(t('hsq2'), RG))))

    x['h2_int'] = list(map(t('intercept'), list(map(t('hsq2'), RG))))
    x['h2_int_se'] = list(map(t('intercept_se'), list(map(t('hsq2'), RG))))
    x['gcov_int'] = list(map(t('intercept'), list(map(t('gencov'), RG))))
    x['gcov_int_se'] = list(map(t('intercept_se'), list(map(t('gencov'), RG))))
    return x.to_string(header=True, index=False) + '\n'


def _print_gencor(args, log, rghat, ref_ld_cnames, i, rg_paths, print_hsq1):
    l = lambda x: x + ''.join(['-' for i in range(len(x.replace('\n', '')))])
    P = [args.samp_prev[0], args.samp_prev[i + 1]]
    K = [args.pop_prev[0], args.pop_prev[i + 1]]
    if args.samp_prev is None and args.pop_prev is None:
        args.samp_prev = [None, None]
        args.pop_prev = [None, None]
    if print_hsq1:
        log.log(l('\nHeritability of phenotype 1\n'))
        log.log(rghat.hsq1.summary(ref_ld_cnames, P=P[0], K=K[0]))

    log.log(
        l('\nHeritability of phenotype {I}/{N}\n'.format(I=i + 2, N=len(rg_paths))))
    log.log(rghat.hsq2.summary(ref_ld_cnames, P=P[1], K=K[1]))
    log.log(l('\nGenetic Covariance\n'))
    log.log(rghat.gencov.summary(ref_ld_cnames, P=P, K=K))
    log.log(l('\nGenetic Correlation\n'))
    log.log(rghat.summary() + '\n')


def _merge_sumstats_sumstats(args, sumstats1, sumstats2, log):
    '''Merge two sets of summary statistics.'''
    sumstats1.rename(columns={'N': 'N1', 'Z': 'Z1'}, inplace=True)
    sumstats2.rename(
        columns={'A1': 'A1x', 'A2': 'A2x', 'N': 'N2', 'Z': 'Z2'}, inplace=True)
    x = _merge_and_log(sumstats1, sumstats2, 'summary statistics', log)
    return x


def _filter_alleles(alleles):
    '''Remove bad variants (mismatched alleles, non-SNPs, strand ambiguous).'''
    ii = alleles.apply(lambda y: y in MATCH_ALLELES)
    return ii


def _align_alleles(z, alleles):
    '''Align Z1 and Z2 to same choice of ref allele (allowing for strand flip).'''
    try:
        z *= (-1) ** alleles.apply(lambda y: FLIP_ALLELES[y])
    except KeyError as e:
        msg = 'Incompatible alleles in .sumstats files: %s. ' % e.args
        msg += 'Did you forget to use --merge-alleles with munge_sumstats.py?'
        raise KeyError(msg)
    return z


def _rg(sumstats, args, log, M_annot, ref_ld_cnames, w_ld_cname, i):
    '''Run the regressions.'''
    n_snp = len(sumstats)
    s = lambda x: np.array(x).reshape((n_snp, 1))
    if args.chisq_max is not None:
        ii = sumstats.Z1**2*sumstats.Z2**2 < args.chisq_max**2
        n_snp = np.sum(ii)  # lambdas are late binding, so this works
        sumstats = sumstats[ii]
    n_blocks = min(args.n_blocks, n_snp)
    ref_ld = sumstats.as_matrix(columns=ref_ld_cnames)
    intercepts = [args.intercept_h2[0], args.intercept_h2[
        i + 1], args.intercept_gencov[i + 1]]
    rghat = reg.RG(s(sumstats.Z1), s(sumstats.Z2),
                   ref_ld, s(sumstats[w_ld_cname]), s(
                       sumstats.N1), s(sumstats.N2), M_annot,
                   intercept_hsq1=intercepts[0], intercept_hsq2=intercepts[1],
                   intercept_gencov=intercepts[2], n_blocks=n_blocks, twostep=args.two_step)

    return rghat


def _parse_rg(rg):
    '''Parse args.rg.'''
    rg_paths = _splitp(rg)
    rg_files = [x.split('/')[-1] for x in rg_paths]
    if len(rg_paths) < 2:
        raise ValueError(
            'Must specify at least two phenotypes for rg estimation.')

    return rg_paths, rg_files


def _print_rg_delete_values(rg, fh, log):
    '''Print block jackknife delete values.'''
    _print_delete_values(rg.hsq1, fh + '.hsq1.delete', log)
    _print_delete_values(rg.hsq2, fh + '.hsq2.delete', log)
    _print_delete_values(rg.gencov, fh + '.gencov.delete', log)


def _print_rg_cov(rghat, fh, log):
    '''Print covariance matrix of estimates.'''
    _print_cov(rghat.hsq1, fh + '.hsq1.cov', log)
    _print_cov(rghat.hsq2, fh + '.hsq2.cov', log)
    _print_cov(rghat.gencov, fh + '.gencov.cov', log)


def _split_or_none(x, n):
    if x is not None:
        y = list(map(float, x.replace('N', '-').split(',')))
    else:
        y = [None for _ in range(n)]
    return y


def _check_arg_len(x, n):
    x, m = x
    if len(x) != n:
        raise ValueError(
            '{M} must have the same number of arguments as --rg/--h2.'.format(M=m))
