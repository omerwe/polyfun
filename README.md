# PolyFun  and  PolyLoc
**PolyFun** (POLYgenic FUNctionally-informed fine-mapping)
<br>
**PolyLoc** (POLYgenic LOCalization of complex trait heritability)

This page contains the code of the methods **PolyFun** for functionally-informed fine-mapping and **PolyLoc** for polygenic localization of complex trait heritability, described in [Weissbrod et al. 2019 bioRxiv](https://www.biorxiv.org/content/10.1101/807792v2).
<br><br>
**PolyFun** estimates prior causal probabilities for SNPs, which can then be used by fine-mapping methods like [SuSiE](https://github.com/stephenslab/susieR) or [FINEMAP](http://www.christianbenner.com/). Unlike previous methods for functionally-informed fine-mapping, **PolyFun** can aggregate polygenic data from across the entire genome and hundreds of functional annotations.
<br><br>
**PolyLoc** generalizes fine-mapping by constructing minimal sets of SNPs that causally explain a given proportion (e.g. 50%) of SNP heritability.

We also provide a script called **finemapper** that facilitates fine-mapping with methods like SuSiE, saving many of the preprocessing steps often required to perform fine-mapping (e.g. handling allelic flips between the summary statistics and reference genotypes).

<br>
<br>

# Manual
We provide a detailed manual of PolyFun, PolyLoc and finemapper in the [Wiki page](https://github.com/omerwe/polyfun/wiki).



<br><br>
# Installation
PolyFun and PolyLoc are designed for Python 3, and require the following freely available Python packages:
* [numpy](http://www.numpy.org/) and [scipy](http://www.scipy.org/)
* [scikit-learn](http://scikit-learn.org/stable/)
* [pandas](https://pandas.pydata.org/getpandas.html) (version >=0.24.0)
* [tqdm](https://github.com/tqdm/tqdm)
* [pyarrow](https://arrow.apache.org/docs/python/install.html)
* [networkx](https://github.com/networkx/networkx) (only required for HESS-based estimation of effect size variance)

It is recommended (but not required) to also install the following:
* [rpy2](https://rpy2.bitbucket.io/)  (a Python package)
* [R version 3.5.1 or higher](https://www.r-project.org/)
* [Ckmeans.1d.dp](https://cran.r-project.org/web/packages/Ckmeans.1d.dp/index.html) (a package for R, that will be invoked from python via the rpy2 package).

If rpy2 or Ckmeans.1d.dp are not installed, PolyFun and PolyLoc will fallback to suboptimal clustering via scikit-learn.

The `finemapper` script also requires the  [LDstore](http://www.christianbenner.com/) program, which should be installed on your system. 


We recommend running PolyFun/PolyLoc via the [Anaconda Python distribution](https://www.anaconda.com/download/). In Anaconda, you can install all the Python packages with the command "conda install \<package_name\>". Alternatively, the Python packages can be installed with the command "pip install --user \<package_name\>".

Once all the prerequisite packages are installed, you can install PolyFun/PolyLoc on a git-enabled machine by typing:
```
git clone https://github.com/omerwe/polyfun
```

<br><br>

# Testing the installation
We recommend testing PolyFun by invoking the script:
```
test_polyfun.py --ldstore <ldstore_executable>
```
where `ldstore_executable` is the path to the [LDstore](http://www.christianbenner.com/) executable on your system. If the script completes without an error, everything is fine. If you omit the `--ldstore` flag, the `finemapper` functionality will not be tested.




<br><br>
# Contact
For questions and comments, please open a Github issue (preferred) or contact Omer Weissbrod at oweissbrod[at]hsph.harvard.edu



