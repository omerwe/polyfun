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
We provide a detailed manual of PolyFun, PolyLoc and finemapper in the [Wiki page](https://github.com/omerwe/polyfun/wiki). If you run into any issues, please check [the FAQ](https://github.com/omerwe/polyfun/wiki/5.-FAQ) first.



<br><br>
# Installation

We provide several installation options.

## Install option 1: Create an Anaconda environment
The easiest way to install polyfun is by creating a dedicated environment through the [Anaconda Python distribution](https://www.anaconda.com/download). To do this, please install Anaconda on your machine and then type the following commands:
```
git clone https://github.com/omerwe/polyfun
cd polyfun
conda env create -f polyfun.yml
conda activate polyfun
```
This will install all the dependencies except for [SuSiE](https://github.com/stephenslab/susieR) 
and [LDstore](http://www.christianbenner.com).
You can use PolyFun without these packages to compute prior causal probabilities, but you won't be able to apply the actual fine-mapping. Please see installation instructions for these two packages below.

After the installation, you can always invoke the PolyFun environment with the command `conda activate polyfun`.
We recommend that you frequently make sure you have the latest version of polyfun installed by going to the polyfun directory and typing `git pull`.


## Install option 2: Manually install packages
PolyFun and PolyLoc are designed for Python >=3.6 and require the following freely available Python packages:
* [numpy](http://www.numpy.org/) and [scipy](http://www.scipy.org/)
* [scikit-learn](http://scikit-learn.org/stable/)
* [pandas](https://pandas.pydata.org/getpandas.html) (version >=0.25.0)
* [tqdm](https://github.com/tqdm/tqdm)
* [pyarrow](https://arrow.apache.org/docs/python/install.html)
* [bitarray](https://github.com/ilanschnell/bitarray)
* [networkx](https://github.com/networkx/networkx) (only required for HESS-based estimation of effect size variance)

It is recommended (but not required) to also install the following:
* [rpy2](https://rpy2.bitbucket.io/)  (a Python package)
* [R version 3.5.1 or higher](https://www.r-project.org/)
* [Ckmeans.1d.dp](https://cran.r-project.org/web/packages/Ckmeans.1d.dp/index.html) (a package for R, that will be invoked from python via the rpy2 package).

If rpy2 or Ckmeans.1d.dp are not installed, PolyFun and PolyLoc will fallback to suboptimal clustering via scikit-learn.

The `finemapper` script also requires the R package [susieR](https://github.com/stephenslab/susieR) and the program [LDstore](http://www.christianbenner.com/) program, which should be installed on your system. Please see installation instructions for these packages below.

We recommend running PolyFun/PolyLoc via the [Anaconda Python distribution](https://www.anaconda.com/download/). In Anaconda, you can install all the Python packages with the command "conda install \<package_name\>". Alternatively, the Python packages can be installed with the command "pip install --user \<package_name\>".

Once all the prerequisite packages are installed, you can install PolyFun/PolyLoc on a git-enabled machine by typing:
```
git clone https://github.com/omerwe/polyfun
```
We recommend that you frequently make sure you have the latest version of polyfun installed by going to the polyfun directory and typing `git pull`.



## Installing SuSiE
To install SuSiE, please start an R shell (usually by typing `R`) and then type: <br>
```
devtools::install_github("stephenslab/susieR@0.8.0",build_vignettes=FALSE)
```
If this doesn't work, please refer to the [SuSiE website](https://github.com/stephenslab/susieR) for more information.

## Installing LDstore
To install LDstore, please type one of the following two commands:
<br>
If you use Linux:
```
wget http://www.christianbenner.com/ldstore_v1.1_x86_64.tgz
tar xzvf ldstore_v1.1_x86_64.tgz
```
If you use Mac OS X :
```
wget http://www.christianbenner.com/ldstore_v1.1_MacOSX.tgz
tar xzvf ldstore_v1.1_MacOSX.tgz
```


<br><br>

## Testing the installation
We recommend testing PolyFun by invoking the script:
```
python test_polyfun.py --python3 <python3_exe> --ldstore <ldstore_executable>
```
where `python3_exe` (optional) is the command you type to start a python3 session (default is `python`) and `ldstore_executable` is the path to the [LDstore](http://www.christianbenner.com/) executable on your system. If the script completes without an error, everything is fine. If you omit the `--ldstore` flag, the `finemapper` functionality will not be tested. If you see any errors, please consult the [FAQ](https://github.com/omerwe/polyfun/wiki/5.-FAQ).




<br><br>
# Contact
For questions and comments, please open a Github issue (preferred) or contact Omer Weissbrod at oweissbrod[at]hsph.harvard.edu



