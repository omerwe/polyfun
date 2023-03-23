# PolyFun, PolyLoc and PolyPred
**PolyFun** (POLYgenic FUNctionally-informed fine-mapping)
<br>
**PolyLoc** (POLYgenic LOCalization of complex trait heritability)
<br>
**PolyPred** (POLYgenic Prediction of complex traits)

This page contains the code of the methods **PolyFun** for functionally-informed fine-mapping, **PolyLoc** for polygenic localization of complex trait heritability, and **PolyPred** for complex trait prediction. **PolyFun** and **PolyLoc** are described in [Weissbrod et al. 2020 Nat Genet](https://www.nature.com/articles/s41588-020-00735-5). **PolyPred** is described in [Weissbrod*, Kanai*, Shi* et al. 2022 Nat Genet](https://www.nature.com/articles/s41588-022-01036-9).
<br><br>
**PolyFun** estimates prior causal probabilities for SNPs, which can then be used by fine-mapping methods like [SuSiE](https://github.com/stephenslab/susieR) or [FINEMAP](http://www.christianbenner.com/). Unlike previous methods for functionally-informed fine-mapping, **PolyFun** can aggregate polygenic data from across the entire genome and hundreds of functional annotations.
<br><br>
**PolyLoc** generalizes fine-mapping by constructing minimal sets of SNPs that causally explain a given proportion (e.g. 50%) of SNP heritability.
<br><br>
**PolyPred** exploits fine-mapping to improve cross-population polygenic risk scores, by predicting using causal effect estimates intead of tagging effect estimates.

We also provide a script called **finemapper** that facilitates fine-mapping with methods like SuSiE, saving many of the preprocessing steps often required to perform fine-mapping (e.g. handling allelic flips between the summary statistics and reference genotypes).

The files in the `ldstore` directory are an adaptation of the [ldstore package](https://pypi.org/project/ldstore/) (written by [Christian Benner](http://www.christianbenner.com)) to Python 3.

<br>
<br>


# Manual
We provide a detailed manual of PolyFun, PolyLoc and PolyPred in the [Wiki page](https://github.com/omerwe/polyfun/wiki). If you run into any issues, please check [the FAQ](https://github.com/omerwe/polyfun/wiki/7.-FAQ) first.



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
 ⚡ **Note**: You can speed up the installation by ~100x by [installing mamba](https://github.com/mamba-org/mamba) and then replacing mamba with conda in the commands above.


This will install all the dependencies except for [FINEMAP](http://www.christianbenner.com) and [LDstore](http://www.christianbenner.com)
This will allow you to perform fine-mapping using SuSiE, but not using FINEMAP. Please see installation instructions for FINEMAP and LDstore below.

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
* [pandas-plink](https://github.com/limix/pandas-plink)
* [r-susier](https://anaconda.org/conda-forge/r-susier)

It is recommended (but not required) to also install the following:
* [rpy2](https://rpy2.bitbucket.io/)  (a Python package)
* [R version 3.5.1 or higher](https://www.r-project.org/)
* [Ckmeans.1d.dp](https://cran.r-project.org/web/packages/Ckmeans.1d.dp/index.html) (a package for R, that will be invoked from python via the rpy2 package).

If rpy2 or Ckmeans.1d.dp are not installed, PolyFun and PolyLoc will fallback to suboptimal clustering via scikit-learn.

If you'd like to use FINEMAP instead of SuSiE for fine-mappping, you will also require:
1. [FINEMAP v1.4.1](http://www.christianbenner.com).
2. (optional) The program [LDstore 2.0](http://www.christianbenner.com) for computing LD directly from .bgen files (imputed genotypes)

Please see installation instructions for these packages below.

We recommend running PolyFun/PolyLoc via the [Anaconda Python distribution](https://www.anaconda.com/download/). In Anaconda, you can install all the Python packages with the command "conda install \<package_name\>". Alternatively, the Python packages can be installed with the command "pip install --user \<package_name\>".

Once all the prerequisite packages are installed, you can install PolyFun/PolyLoc on a git-enabled machine by typing:
```
git clone https://github.com/omerwe/polyfun
```
We recommend that you frequently make sure you have the latest version of polyfun installed by going to the polyfun directory and typing `git pull`.



## Installing FINEMAP v1.4.1
To install FINEMAP v1.4.1, please type one of the following two commands:
<br>
If you use Linux:
```
wget http://christianbenner.com/finemap_v1.4.1_x86_64.tgz
tar xvf finemap_v1.4.1_x86_64.tgz
```
If you use Mac OS X :
```
wget http://christianbenner.com/finemap_v1.4.1_MacOSX.tgz
tar xvf finemap_v1.4.1_MacOSX.tgz
```

## Installing LDstore 2.0
To install LDstore, please type one of the following two commands:
<br>
If you use Linux:
```
wget http://www.christianbenner.com/ldstore_v2.0_x86_64.tgz
tar xzvf ldstore_v2.0_x86_64.tgz
```
If you use Mac OS X :
```
wget http://www.christianbenner.com/ldstore_v2.0_MacOSX.tgz
tar xzvf ldstore_v2.0_MacOSX.tgz
```


<br><br>

## Testing the installation
We recommend testing PolyFun by invoking the script:
```
python test_polyfun.py --python3 <python3_exe>
```
where `python3_exe` (optional) is the command you type to start a python3 session (default is `python`). If the script completes without an error, everything is fine. If you see any errors, please consult the [FAQ](https://github.com/omerwe/polyfun/wiki/7.-FAQ).

To test FINEMAP integration, provide the path to the executable when invoking the script (if not provided, the script will not test FINEMAP integration):
```
python test_polyfun.py --python3 <python3_exe> --finemap-exe <finemap_exe>
```



<br><br>
# Contact
For questions and comments, please open a Github issue.



