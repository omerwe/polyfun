'''
(c) 2014 Brendan Bulik-Sullivan and Hilary Finucane

Fast block jackknives.

Everything in this module deals with 2D numpy arrays. 1D data are represented as arrays
with dimension (N, 1) or (1, N), to avoid bugs arising from numpy treating (N, ) as
a fundamentally different shape from (N, 1). The convention in this module is for the
first dimension to represent # of data points (or # of blocks in a block jackknife, since
a block is like a datapoint), and for the second dimension to represent the dimensionality
of the data.

'''


import numpy as np
from scipy.optimize import nnls
np.seterr(divide='raise', invalid='raise')
from tqdm import tqdm
from sklearn.linear_model import Lasso
import logging
import warnings
warnings.filterwarnings('ignore', message='Coordinate descent with alpha=0 may lead to unexpected results and is discouraged.')
warnings.filterwarnings('ignore', message='Objective did not converge. You might want to increase the number of iterations. Fitting data with very small alpha may cause precision problems.')

from sklearn.metrics import r2_score

def _check_shape(x, y):
    '''Check that x and y have the correct shapes (for regression jackknives).'''
    if len(x.shape) != 2 or len(y.shape) != 2:
        raise ValueError('x and y must be 2D arrays.')
    if x.shape[0] != y.shape[0]:
        raise ValueError(
            'Number of datapoints in x != number of datapoints in y.')
    if y.shape[1] != 1:
        raise ValueError('y must have shape (n_snp, 1)')
    n, p = x.shape
    if p > n:
        raise ValueError('More dimensions than datapoints.')

    return (n, p)


def _check_shape_block(xty_block_values, xtx_block_values):
    '''Check that xty_block_values and xtx_block_values have correct shapes.'''
    if xtx_block_values.shape[0:2] != xty_block_values.shape:
        raise ValueError(
            'Shape of xty_block_values must equal shape of first two dimensions of xty_block_values.')
    if len(xtx_block_values.shape) < 3:
        raise ValueError('xtx_block_values must be a 3D array.')
    if xtx_block_values.shape[1] != xtx_block_values.shape[2]:
        raise ValueError(
            'Last two axes of xtx_block_values must have same dimension.')

    return xtx_block_values.shape[0:2]


class Jackknife(object):

    '''
    Base class for jackknife objects. Input involves x,y, so this base class is tailored
    for statistics computed from independent and dependent variables (e.g., regressions).
    The __delete_vals_to_pseudovalues__ and __jknife__ methods will still be useful for other
    sorts of statistics, but the __init__ method will need to be overriden.

    Parameters
    ----------
    x : np.matrix with shape (n, p)
        Independent variable.
    y : np.matrix with shape (n, 1)
        Dependent variable.
    n_blocks : int
        Number of jackknife blocks
    *args, **kwargs :
        Arguments for inheriting jackknives.

    Attributes
    ----------
    n_blocks : int
        Number of jackknife blocks
    p : int
        Dimensionality of the independent varianble
    N : int
        Number of datapoints (equal to x.shape[0])

    Methods
    -------
    jknife(pseudovalues):
        Computes jackknife estimate and variance from the jackknife pseudovalues.
    delete_vals_to_pseudovalues(delete_vals, est):
        Converts delete values and the whole-data estimate to pseudovalues.
    get_separators():
        Returns (approximately) evenly-spaced jackknife block boundaries.
    '''

    def __init__(self, x, y, n_blocks=None, separators=None):
        self.N, self.p = _check_shape(x, y)
        if separators is not None:
            if max(separators) != self.N:
                raise ValueError(
                    'Max(separators) must be equal to number of data points.')
            if min(separators) != 0:
                raise ValueError('Max(separators) must be equal to 0.')
            self.separators = sorted(separators)
            self.n_blocks = len(separators) - 1
        elif n_blocks is not None:
            self.n_blocks = n_blocks
            self.separators = self.get_separators(self.N, self.n_blocks)
        else:
            raise ValueError('Must specify either n_blocks are separators.')

        if self.n_blocks > self.N:
            raise ValueError('More blocks than data points.')

    @classmethod
    def jknife(cls, pseudovalues):
        '''
        Converts pseudovalues to jackknife estimate and variance.

        Parameters
        ----------
        pseudovalues : np.matrix pf floats with shape (n_blocks, p)

        Returns
        -------
        jknife_est : np.matrix with shape (1, p)
            Jackknifed estimate.
        jknife_var : np.matrix with shape (1, p)
            Variance of jackknifed estimate.
        jknife_se : np.matrix with shape (1, p)
            Standard error of jackknifed estimate, equal to sqrt(jknife_var).
        jknife_cov : np.matrix with shape (p, p)
            Covariance matrix of jackknifed estimate.

        '''
        
        n_blocks = pseudovalues.shape[0]
        jknife_cov = np.atleast_2d(np.cov(pseudovalues.T, ddof=1) / n_blocks)
        jknife_var = np.atleast_2d(np.diag(jknife_cov))
        jknife_se = np.atleast_2d(np.sqrt(jknife_var))
        jknife_est = np.atleast_2d(np.mean(pseudovalues, axis=0))
        return (jknife_est, jknife_var, jknife_se, jknife_cov)

    @classmethod
    def delete_values_to_pseudovalues(cls, delete_values, est):
        '''
        Converts whole-data estimate and delete values to pseudovalues.

        Parameters
        ----------
        delete_values : np.matrix with shape (n_blocks, p)
            Delete values.
        est : np.matrix with shape (1, p):
            Whole-data estimate.

        Returns
        -------
        pseudovalues : np.matrix with shape (n_blocks, p)
            Psuedovalues.

        Raises
        ------
        ValueError :
            If est.shape != (1, delete_values.shape[1])

        '''
        
        n_blocks, p = delete_values.shape
        if est.shape != (1, p):
            raise ValueError(
                'Different number of parameters in delete_values than in est.')

        return n_blocks * est - (n_blocks - 1) * delete_values

    @classmethod
    def get_separators(cls, N, n_blocks):
        '''Define evenly-spaced block boundaries.'''
        return np.floor(np.linspace(0, N, n_blocks + 1)).astype(int)



        
class LstsqJackknifeSlow(Jackknife):

    '''
    Slow linear-regression block jackknife. This class computes delete values directly,
    rather than forming delete values from block values. Useful for testing and for
    non-negative least squares (which as far as I am aware does not admit a fast block
    jackknife algorithm).

    Inherits from Jackknife class.

    Parameters
    ----------
    x : np.matrix with shape (n, p)
        Independent variable.
    y : np.matrix with shape (n, 1)
        Dependent variable.
    n_blocks : int
        Number of jackknife blocks
    nn: bool
        Non-negative least-squares?

    Attributes
    ----------
    est : np.matrix with shape (1, p)
        FWLS estimate.
    jknife_est : np.matrix with shape (1, p)
        Jackknifed estimate.
    jknife_var : np.matrix with shape (1, p)
        Variance of jackknifed estimate.
    jknife_se : np.matrix with shape (1, p)
        Standard error of jackknifed estimate, equal to sqrt(jknife_var).
    jknife_cov : np.matrix with shape (p, p)
        Covariance matrix of jackknifed estimate.
    delete_vals : np.matrix with shape (n_blocks, p)
        Jackknife delete values.
    '''

    @classmethod
    def delete_values(cls, x, y, func, s):
        '''
        Compute delete values by deleting one block at a time.

        Parameters
        ----------
        x : np.matrix with shape (n, p)
            Independent variable.
        y : np.matrix with shape (n, 1)
            Dependent variable.
        func : function (n, p) , (n, 1) --> (1, p)
            Function of x and y to be jackknived.
        s : list of ints
            Block separators.

        Returns
        -------
        delete_values : np.matrix with shape (n_blocks, p)
            Delete block values (with n_blocks blocks defined by parameter s).

        Raises
        ------
        ValueError :
            If x.shape[0] does not equal y.shape[0] or x and y are not 2D.

        '''
        _check_shape(x, y)
        d = []
        logging.info('Starting non-negative jackknife...')
        for i in tqdm(range(len(s) - 1)):
            jk_est = func(np.vstack([x[0:s[i], ...], x[s[i + 1]:, ...]]), np.vstack([y[0:s[i], ...], y[s[i + 1]:, ...]]))
            d.append(jk_est)
        
        return np.concatenate(d, axis=0)
        
    
    def __init__(self, x, y, is_large_chi2, n_blocks=None, nn=False, separators=None, chr_num=None, evenodd_split=False):
    
        Jackknife.__init__(self, x, y, n_blocks, separators)
    
        #estimate taus
        if nn:  # non-negative least squares
            #func = lambda x, y: np.atleast_2d(nnls(x, np.array(y).T[0])[0])
            xtx = x.T.dot(x)
            lasso = Lasso(alpha=1e-100, fit_intercept=False, normalize=False, precompute=xtx, positive=True, max_iter=10000)
            self.est = lasso.fit(x,y[:,0]).coef_.reshape((1, x.shape[1]))
        else:
            func = lambda x, y: np.atleast_2d(
                np.linalg.lstsq(x, np.array(y).T[0])[0])
            self.est = func(x, y)        
        
        #move large_chi2 SNPs to the end of x and y (don't include them in the separator definition, so that they'll never get removed during jackknife)
        if np.any(is_large_chi2):
            x_large = x[is_large_chi2]
            y_large = y[is_large_chi2]
            x = x[~is_large_chi2]
            y = y[~is_large_chi2]
            Jackknife.__init__(self, x, y, n_blocks, separators)
            x = np.concatenate((x,x_large), axis=0)
            y = np.concatenate((y,y_large), axis=0)
        
        #jackknife
        if nn:
            d = []
            s = self.separators
            for i in tqdm(range(len(s) - 1), disable=False):
                x_noblock = np.delete(x, slice(s[i], s[i+1]), axis=0)
                y_noblock = np.delete(y, slice(s[i], s[i+1]), axis=0)
                x_block = x[s[i] : s[i+1]]
                xtx_noblock = xtx - x_block.T.dot(x_block)
                lasso_noblock = Lasso(alpha=1e-100, fit_intercept=False, normalize=False, precompute=xtx_noblock, positive=True, max_iter=10000)
                jk_est = lasso_noblock.fit(x_noblock, y_noblock[:,0]).coef_.reshape((1, x.shape[1]))
                ###z = nnls(x_noblock, y_noblock[:,0])[0]
                ###assert np.allclose(z, jk_est[0])
                d.append(jk_est)    
            self.delete_values = np.concatenate(d, axis=0)
        else:
            self.delete_values = self.delete_values(x, y, func, self.separators)
        
        
        self.pseudovalues = self.delete_values_to_pseudovalues(
            self.delete_values, self.est)
        (self.jknife_est, self.jknife_var, self.jknife_se, self.jknife_cov) =\
            self.jknife(self.pseudovalues)
    
        if evenodd_split:
            assert y.shape[1]==1
            assert chr_num is not None            
            assert len(np.unique(chr_num)) > 1                        
            self.chr_list = np.sort(np.unique(chr_num))
            self.est_loco = np.empty((len(self.chr_list), x.shape[1]), dtype=np.float32)
            for chr_i, left_out_chr in enumerate(tqdm(self.chr_list)):
                is_loco = ((chr_num%2)==(left_out_chr%2)) & (chr_num != left_out_chr)
                x_loco = x[is_loco]
                y_loco = y[is_loco]
                self.est_loco[chr_i, :] = nnls(x_loco, y_loco[:,0])[0]

        
        
class LstsqJackknifeFast(Jackknife):
    
    def __init__(self, x, y, is_large_chi2, n_blocks=None, separators=None, chr_num=None, evenodd_split=False):
    
        #compute jackknife estimates using all SNPs
        Jackknife.__init__(self, x, y, n_blocks, separators)
        xty, xtx = self.block_values(x, y, self.separators)
        self.est = self.block_values_to_est(xty, xtx)
        
        #compute xtx_tot and xty_tot
        xty_tot = np.sum(xty, axis=0)
        xtx_tot = np.sum(xtx, axis=0)
    
        #exclude large-chi2 SNPs from xtx and xty for the jackknife
        if np.any(is_large_chi2):
            x = x[~is_large_chi2]
            y = y[~is_large_chi2]
            Jackknife.__init__(self, x, y, n_blocks, separators)
            xty, xtx = self.block_values(x, y, self.separators)
        
        self.delete_values = self.block_values_to_delete_values(xty, xtx, xtx_tot=xtx_tot, xty_tot=xty_tot)
        self.pseudovalues = self.delete_values_to_pseudovalues(
            self.delete_values, self.est)
        (self.jknife_est, self.jknife_var, self.jknife_se, self.jknife_cov) =\
            self.jknife(self.pseudovalues)
            
            
        if evenodd_split:
            assert y.shape[1]==1
            assert chr_num is not None
            assert len(np.unique(chr_num)) > 1            
            
            x_even = x[chr_num %2 == 0]
            y_even = y[chr_num %2 == 0]
            XTX_even = x_even.T.dot(x_even)
            XTy_even = y_even[:,0].dot(x_even)
            del x_even, y_even
            
            x_odd = x[chr_num %2 != 0]
            y_odd = y[chr_num %2 != 0]
            XTX_odd = x_odd.T.dot(x_odd)
            XTy_odd = y_odd[:,0].dot(x_odd)
            del x_odd, y_odd            
            
            assert np.allclose(XTy_even + XTy_odd, y[:,0].dot(x))        
            assert np.allclose(XTX_even + XTX_odd, x.T.dot(x))
            
            self.chr_list = np.sort(np.unique(chr_num))
            #self.est_chr  = np.empty((len(self.chr_list), x.shape[1]), dtype=np.float32)
            self.est_loco = np.empty((len(self.chr_list), x.shape[1]), dtype=np.float32)            
            for chr_i, left_out_chr in enumerate(self.chr_list):
                x_chr = x[chr_num == left_out_chr]
                y_chr = y[chr_num == left_out_chr, 0]
                XTX_chr = x_chr.T.dot(x_chr)
                XTy_chr = y_chr.dot(x_chr)
                
                if left_out_chr %2 == 0:
                    XTX_loco = XTX_even - XTX_chr
                    XTy_loco = XTy_even - XTy_chr
                else:
                    XTX_loco = XTX_odd - XTX_chr
                    XTy_loco = XTy_odd - XTy_chr
                
                self.est_loco[chr_i, :] = np.linalg.solve(XTX_loco, XTy_loco)
                #self.est_chr[chr_i, :] = np.linalg.solve(XTX_chr, XTy_chr)
                    
            
            
    @classmethod
    def block_values_to_delete_values(cls, xty_block_values, xtx_block_values, xtx_tot, xty_tot):
        n_blocks, p = _check_shape_block(xty_block_values, xtx_block_values)
        delete_values = np.zeros((n_blocks, p))
        for j in range(n_blocks):
            delete_xty = xty_tot - xty_block_values[j]
            delete_xtx = xtx_tot - xtx_block_values[j]
            delete_values[j, ...] = np.linalg.solve(
                delete_xtx, delete_xty).reshape((1, p))

        return delete_values
        
        
    @classmethod
    def block_values(cls, x, y, s):
        '''
        Compute block values.

        Parameters
        ----------
        x : np.matrix with shape (n, p)
            Independent variable.
        y : np.matrix with shape (n, 1)
            Dependent variable.
        n_blocks : int
            Number of jackknife blocks
        s : list of ints
            Block separators.

        Returns
        -------
        xty_block_values : np.matrix with shape (n_blocks, p)
            Block values of X^T Y.
        xtx_block_values : 3d np array with shape (n_blocks, p, p)
            Block values of X^T X.

        Raises
        ------
        ValueError :
            If x.shape[0] does not equal y.shape[0] or x and y are not 2D.

        '''
        n, p = _check_shape(x, y)
        n_blocks = len(s) - 1
        xtx_block_values = np.zeros((n_blocks, p, p))
        xty_block_values = np.zeros((n_blocks, p))
        for i in range(n_blocks):
            xty_block_values[i, ...] = np.dot(
                x[s[i]:s[i + 1], ...].T, y[s[i]:s[i + 1], ...]).reshape((1, p))
            xtx_block_values[i, ...] = np.dot(
                x[s[i]:s[i + 1], ...].T, x[s[i]:s[i + 1], ...])

        return (xty_block_values, xtx_block_values)        
        
        
    @classmethod
    def block_values_to_est(cls, xty_block_values, xtx_block_values):
        '''
        Converts block values to the whole-data linear regression estimate.

        Parameters
        ----------
        xty_block_values : np.matrix with shape (n_blocks, p)
            Block values of X^T Y.
        xtx_block_values : 3D np.array with shape (n_blocks, p, p)
            Block values of X^T X

        Returns
        -------
        est : np.matrix with shape (1, p)
            Whole data estimate.

        Raises
        ------
        LinAlgError :
            If design matrix is singular.
        ValueError :
            If the last two dimensions of xtx_block_values are not equal or if the first two
        dimensions of xtx_block_values do not equal the shape of xty_block_values.

        '''
        n_blocks, p = _check_shape_block(xty_block_values, xtx_block_values)
        xty = np.sum(xty_block_values, axis=0)
        xtx = np.sum(xtx_block_values, axis=0)
        return np.linalg.solve(xtx, xty).reshape((1, p))
        
            


class RatioJackknife(Jackknife):

    '''
    Block jackknife ratio estimate.

    Jackknife.

    Parameters
    ----------
    est : float or np.array with shape (1, p)
        Whole data ratio estimate
    numer_delete_values : np.matrix with shape (n_blocks, p)
        Delete values for the numerator.
    denom_delete_values: np.matrix with shape (n_blocks, p)
        Delete values for the denominator.

    Methods
    -------
    delete_vals_to_pseudovalues(est, denom, num):
        Converts denominator/ numerator delete values and the whole-data estimate to
        pseudovalues.

    Raises
    ------
    FloatingPointError :
        If any entry of denom_delete_values is zero.

    Note that it is possible for the denominator to cross zero (i.e., be both positive
    and negative) and still have a finite ratio estimate and SE, for example if the
    numerator is fixed to 0 and the denominator is either -1 or 1. If the denominator
    is noisily close to zero, then it is unlikely that the denominator will yield zero
    exactly (and therefore yield an inf or nan), but delete values will be of the form
    (numerator / close to zero) and -(numerator / close to zero), i.e., (big) and -(big),
    and so the jackknife will (correctly) yield huge SE.

    '''

    def __init__(self, est, numer_delete_values, denom_delete_values):
        if numer_delete_values.shape != denom_delete_values.shape:
            raise ValueError(
                'numer_delete_values.shape != denom_delete_values.shape.')
        if len(numer_delete_values.shape) != 2:
            raise ValueError('Delete values must be matrices.')
        if len(est.shape) != 2 or est.shape[0] != 1 or est.shape[1] != numer_delete_values.shape[1]:
            raise ValueError(
                'Shape of est does not match shape of delete values.')

        self.n_blocks = numer_delete_values.shape[0]
        self.est = est
        self.pseudovalues = self.delete_values_to_pseudovalues(self.est,
                                                               denom_delete_values, numer_delete_values)
        (self.jknife_est, self.jknife_var, self.jknife_se, self.jknife_cov) =\
            self.jknife(self.pseudovalues)

    @classmethod
    def delete_values_to_pseudovalues(cls, est, denom, numer):
        '''
        Converts delete values to pseudovalues.

        Parameters
        ----------
        est : np.matrix with shape (1, p)
            Whole-data ratio estimate.
        denom : np.matrix with shape (n_blocks, p)
            Denominator delete values.
        numer : np.matrix with shape (n_blocks, p)
            Numerator delete values.

        Returns
        -------
        pseudovalues :
            Ratio Jackknife Pseudovalues.

        Raises
        ------
        ValueError :
            If numer.shape != denom.shape.

        '''
        n_blocks, p = denom.shape
        pseudovalues = np.zeros((n_blocks, p))
        for j in range(0, n_blocks):
            pseudovalues[j, ...] = n_blocks * est - \
                (n_blocks - 1) * numer[j, ...] / denom[j, ...]

        return pseudovalues




class Jackknife_Ridge(Jackknife):

    def __init__(self, x, y, n_blocks=None, separators=None, chr_num=None, verbose=True,
        num_lambdas=100, approx_ridge=False, 
        ridge_lambda=None, use_1se=False, has_intercept=False, standardize=True,
        skip_ridge_jackknife=True, num_chr_sets=2):
        
        #sanity checks
        assert chr_num is not None
        assert len(np.unique(chr_num)) > 1
        
        #init stuff
        Jackknife.__init__(self, x, y, n_blocks=n_blocks, separators=separators)
        self.use_1se = use_1se
        self.verbose=verbose        
        self.has_intercept = has_intercept
        
        ###define chromosome sets
        assert num_chr_sets>1
        
        if num_chr_sets == 2:
            #Use the good old fashioned odd/even chromosome split
            chromosomes = np.sort(np.unique(chr_num))
            self.chromosome_sets = []
            self.chromosome_sets.append(chromosomes[chromosomes%2==0])
            self.chromosome_sets.append(chromosomes[chromosomes%2!=0])
        elif num_chr_sets == 22:
            self.chromosome_sets = [np.array([c]) for c in range(1,23)]
        else:        
            chr_sizes = np.bincount(chr_num)[1:]
            assert num_chr_sets<=len(chr_sizes)
            chr_assignments = self._divide_chromosomes_to_sets(chr_sizes, num_chr_sets)
            self.chromosome_sets = []
            for set_i in range(num_chr_sets):
                self.chromosome_sets.append(np.where(chr_assignments==set_i)[0]+1)

        #make sure we work with numpy arrays, not dataframes
        try: x=x.values
        except: pass
        try: y=y.values
        except: pass
        try: constraints=constraints.values
        except: pass
        try: chr_num=chr_num.values
        except: pass
        
        #make y look like a vector
        assert y.shape[1]==1
        y = y[:,0]
      
        #standardize x
        if standardize:
            x_l2 = np.sqrt(np.einsum('ij,ij->j', x, x))
            x /= x_l2
        else:
            x_l2 = None
        
        #Create a set of ridge lambdas to evaluate
        XTX_all = x.T.dot(x)
        XTy_all = y.dot(x)
        mean_diag = np.mean(np.diag(XTX_all))
        self.ridge_lambdas = np.logspace(np.log10(mean_diag*1e-8), np.log10(mean_diag*1e2), num=num_lambdas)
        
        #find best lambda (using off-chromosome estimation) and estimate taus
        if ridge_lambda is not None:
            assert self.approx_ridge            
            best_lambda = ridge_lambda
        else:
            best_lambda = self._find_best_lambda(x, y, XTX_all, XTy_all, chr_num)
        self.est = np.atleast_2d(self._est_ridge(XTX_all, XTy_all, best_lambda))
        if standardize:
            self.est /= x_l2
        
        
        #LOCO (leave one chromosome out) computations        
        self.est_chr_lstsq, self.est_chr_ridge, self.est_loco_lstsq, self.est_loco_ridge = \
            self._est_taus_loco(x, y, XTX_all, XTy_all, chr_num, best_lambda, standardize, x_l2)

        #run jackknife
        if not skip_ridge_jackknife:
            self.delete_values = np.empty((len(self.separators)-1, self.est.shape[1]), dtype=np.float32)
            self.est_chr_lstsq_jk_list = []
            self.est_chr_ridge_jk_list = []
            self.est_loco_lstsq_jk_list = []
            self.est_loco_ridge_jk_list = []
            
            logging.info('Running ridge jackknife...')
            for block_i in tqdm(range(len(self.separators) - 1)):
                        
                #prepare data structures
                x_block = x[self.separators[block_i]:self.separators[block_i+1], ...]
                y_block = y[self.separators[block_i]:self.separators[block_i+1], ...]
                XTX_noblock = XTX_all - x_block.T.dot(x_block)
                XTy_noblock = XTy_all - y_block.dot(x_block)
                slice_block = slice(self.separators[block_i], self.separators[block_i+1])
                x_noblock = np.delete(x, slice_block, axis=0)
                y_noblock = np.delete(y, slice_block, axis=0)
                chr_noblock = np.delete(chr_num, slice_block, axis=0)

                #find best lambda for this jackknife block
                if approx_ridge:
                    best_lambda_noblock = best_lambda
                else:
                    best_lambda_noblock = self._find_best_lambda(x_noblock, y_noblock, XTX_noblock, XTy_noblock, chr_noblock)
                
                #main jackknife estimation
                est_block = self._est_ridge(XTX_noblock, XTy_noblock, best_lambda_noblock)            
                self.delete_values[block_i, ...] = est_block
                
                #jackknife LOCO computation
                est_chr_lstsq, est_chr_ridge, est_loco_lstsq, est_loco_ridge = \
                    self._est_taus_loco(x_noblock, y_noblock, XTX_noblock, XTy_noblock, 
                                        chr_noblock, best_lambda_noblock, standardize, x_l2)
                self.est_chr_lstsq_jk_list.append(est_chr_lstsq)
                self.est_chr_ridge_jk_list.append(est_chr_ridge)
                self.est_loco_lstsq_jk_list.append(est_loco_lstsq)
                self.est_loco_ridge_jk_list.append(est_loco_ridge)
            if standardize: self.delete_values /= x_l2

            #compute jackknife pseudo-values
            self.pseudovalues = self.delete_values_to_pseudovalues(self.delete_values, self.est)
            (self.jknife_est, self.jknife_var, self.jknife_se, self.jknife_cov) = self.jknife(self.pseudovalues)
        
        #restore original x
        if standardize: x *= x_l2
        
        
        
    def _divide_chromosomes_to_sets(self, chr_sizes, num_sets):
        chr_order = np.argsort(chr_sizes)[::-1]     #np.arange(len(chr_sizes))
        chr_assignments = np.zeros(22, dtype=np.int) - 1
        chr_assignments[chr_order[:num_sets]] = np.arange(num_sets)
        set_sizes = chr_sizes[chr_order[:num_sets]].copy()    
        for c_i in chr_order[num_sets : len(chr_sizes)]:
            smallest_set = np.argmin(set_sizes)
            chr_assignments[c_i] = smallest_set
            set_sizes[smallest_set] += chr_sizes[c_i]
        assert set_sizes.sum() == chr_sizes.sum()    
        return chr_assignments
                
            

    def _est_taus_loco(self, x, y, XTX, XTy, chr_num, ridge_lambda, standardize, x_l2=None, reestimate_lambda=False):        
        chromosomes = np.sort(np.unique(chr_num))
        est_set_lstsq  = np.empty((len(self.chromosome_sets), x.shape[1]), dtype=np.float32)
        est_noset_lstsq = np.empty((len(self.chromosome_sets), x.shape[1]), dtype=np.float32)
        est_set_ridge  = np.empty((len(self.chromosome_sets), x.shape[1]), dtype=np.float32)
        est_noset_ridge = np.empty((len(self.chromosome_sets), x.shape[1]), dtype=np.float32)
        tqdm_chr_sets = tqdm(self.chromosome_sets)
        logging.info('Estimating annotation coefficients for each chromosomes set')
        for set_i, chromosome_set in enumerate(tqdm_chr_sets):            
            is_in_set = np.isin(chr_num, chromosome_set)
            if not np.any(is_in_set): continue
            x_set = x[is_in_set]
            y_set = y[is_in_set]
            XTX_set = x_set.T.dot(x_set)
            XTy_set = y_set.dot(x_set)
            XTX_noset = XTX - XTX_set
            XTy_noset = XTy - XTy_set                
            
            if (not reestimate_lambda) or (len(chromosomes) <= 2):
                best_lambda_noset = ridge_lambda
                best_lambda_set = ridge_lambda
            else:
                x_loco = x[~is_in_set]
                y_loco = y[~is_in_set]
                chr_loco = chr_num[~is_in_set]                
                best_lambda_noset = self._find_best_lambda(x_loco, y_loco, XTX_noset, XTy_noset, chr_loco)
                if len(chromosome_set) == 1:
                    best_lambda_set = ridge_lambda
                else:
                    best_lambda_set = self._find_best_lambda(x_set, y_set, XTX_set, XTy_set, chr_num[is_in_set])
            est_set_lstsq[set_i, :]   = self._est_ridge(XTX_set, XTy_set, ridge_lambda=0)
            est_set_ridge[set_i, :]   = self._est_ridge(XTX_set, XTy_set, best_lambda_set)
            est_noset_lstsq[set_i, :] = self._est_ridge(XTX_noset, XTy_noset, ridge_lambda=0)
            est_noset_ridge[set_i, :] = self._est_ridge(XTX_noset, XTy_noset, best_lambda_noset)
            
        if standardize:
            est_set_lstsq /= x_l2
            est_set_ridge /= x_l2            
            est_noset_lstsq /= x_l2
            est_noset_ridge /= x_l2
                
        return est_set_lstsq, est_set_ridge, est_noset_lstsq, est_noset_ridge
        
        
    def _find_best_lambda(self, x, y, XTX, XTy, chr_num):
        chromosomes = np.sort(np.unique(chr_num))
        assert len(chromosomes) > 1
        num_lambdas = len(self.ridge_lambdas)
        y_pred_lambdas = np.empty((chr_num.shape[0], num_lambdas), dtype=np.float32)
        if self.verbose:
            y_pred_lambdas_lstsq = np.empty(chr_num.shape[0], dtype=np.float32)
        logging.info('iterating over chromosomes to compute XTX, XTy...')
        for chr_i, left_out_chr in enumerate(tqdm(chromosomes)):        
            is_chr = (chr_num == left_out_chr)
            chr_inds = np.where(is_chr)[0]
            assert np.all(chr_inds == np.arange(chr_inds[0], chr_inds[-1]+1))
            chr_slice = slice(chr_inds[0], chr_inds[-1]+1)
            x_chr = x[chr_slice]
            y_chr = y[chr_slice]            
            XTX_loco = XTX - x_chr.T.dot(x_chr)
            XTy_loco = XTy - y_chr.dot(x_chr)
            y_pred_lambdas[chr_slice, :] = self._predict_lambdas(XTX_loco, XTy_loco, x_chr)
            
            if self.verbose:
                tau_lstsq_loco = self._est_ridge(XTX_loco, XTy_loco, 0)
                y_pred_lambdas_lstsq[chr_slice] = x_chr.dot(tau_lstsq_loco)
            
        #Assign an r2 score to each lambda
        score_lambdas = np.empty(num_lambdas, dtype=np.float32)
        logging.info('Evaluating Ridge lambdas...')
        for r_i in tqdm(range(num_lambdas)):
            score_lambdas[r_i] = r2_score(y, y_pred_lambdas[:,r_i])        
        
        #choose lambda based on the 1SE rule?
        if self.use_1se:
            score_folds = np.empty(len(chromosomes), dtype=np.float32)
            for chr_i, left_out_chr in enumerate(chromosomes):
                is_chr = (chr_num == left_out_chr)
                score_folds[chr_i] = r2_score(y[is_chr], y_pred_lambdas[is_chr, best_lambda_index])
            scores_std = np.std(score_folds)
            best_score = score_lambdas[best_lambda_index]
            assert np.isclose(best_score, np.max(score_lambdas))
            best_lambda_index = np.where(score_lambdas > best_score - scores_std)[0][-1]
        else:
            best_lambda_index = np.argmax(score_lambdas)
        
        best_lambda = self.ridge_lambdas[best_lambda_index]        
        if self.verbose:
            score_lstsq = r2_score(y, y_pred_lambdas_lstsq)
            logging.info('Selected ridge lambda: %0.4e (%d/%d)  score: %0.4e  score lstsq: %0.4e'%(best_lambda, 
                best_lambda_index+1, num_lambdas, score_lambdas[best_lambda_index], score_lstsq))
            
        return best_lambda
        
        
    def _predict_lambdas(self, XTX_train, XTy_train, X_validation):
        tau_est_ridge = np.empty((XTX_train.shape[0], len(self.ridge_lambdas)), dtype=np.float32)
        for r_i, r in enumerate(self.ridge_lambdas):
            tau_est_ridge[:, r_i] = self._est_ridge(XTX_train, XTy_train, r)
        y_pred = X_validation.dot(tau_est_ridge)
        return y_pred
        
    
    def _est_ridge(self, XTX, XTy, ridge_lambda):
        I = np.eye(XTX.shape[0]) * ridge_lambda
        if self.has_intercept: I[-1,-1]=0
        return np.linalg.solve(XTX+I, XTy)



