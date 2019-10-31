# PolyFun
PolyFun (POLYgenic FUNctionally-informed fine-mapping)

PolyFun is a framework for functionally-informed fine-mapping [(see preprint)](https://www.biorxiv.org/content/10.1101/807792v2). PolyFun estimates prior causal probability for SNPs, which can then be used by fine-mapping methods like [SuSiE](https://github.com/stephenslab/susieR) or [FINEMAP](http://www.christianbenner.com/). Unlike previous methods, PolyFun can aggregate polygenic data from across the entire genome and hundreds of functional annotations.

<br><br>
# Installation
PolyFun is designed for Python 3, and requires the following freely available Python packages:
* [numpy](http://www.numpy.org/) and [scipy](http://www.scipy.org/)
* [scikit-learn](http://scikit-learn.org/stable/)
* [pandas](https://pandas.pydata.org/getpandas.html) (version >=0.24.0)
* [tqdm](https://github.com/tqdm/tqdm)
* [pyarrow](https://arrow.apache.org/docs/python/install.html)

It is recommended (but not required) to also install the following:
* [rpy2](https://rpy2.bitbucket.io/)  (a Python package)
* [R version 3.5.1 or higher](https://www.r-project.org/)
* [Ckmeans.1d.dp](https://cran.r-project.org/web/packages/Ckmeans.1d.dp/index.html) (a package for R, that will be invoked from python via the rpy2 package).
If these are not installed, PolyFun will fallback to sub-optimal clustering via scikit-learn.

We recommend running PolyFun via the [Anaconda Python distribution](https://www.anaconda.com/download/). In Anaconda, you can install all the packages with the command "conda install \<package_name\>". Alternatively, the Python packages can be installed with the command "pip install --user \<package_name\>".

Once all the prerequisite packages are installed, you can install PolyFun on a git-enabled machine by typing:
```
git clone https://github.com/omerwe/polyfun
```



<br><br>
# Usage overview
There are three ways to run PolyFun:
1. **Using precomputed prior causal probabilities based on a meta-analysis of 15 [UK Biobank](https://www.ukbiobank.ac.uk) traits**. This is the simplest approach, but the prior causal probabilities may not be optimal for some traits.
2. **Computing prior causal probabilities via the [baseline-LF model annotations](https://www.nature.com/articles/s41588-018-0231-8)**. This is a relatively simple approach, but the prior causal probabilities may not be robust to modeling misspecification.
3. **Computing prior causal probabilities non-parametrically**. This is the most robust approach, but it is computationally intensive and requires access to individual-level genotypic data from a large reference panel (optimally >3,000 population-matched individuals).

Below are instructions on how to use each of these approaches.

<br><br>

# Using precomputed prior causal probabilities based on a meta-analysis of 15 UK Bionank traits
Here, all you need to do is provide a file with SNP identifiers. PolyFun will extract the prior causal probabilities for this SNP. To do this, use the following command:
```
python extract_snpvar.py --snps <snps_file> --out <output_prefix>
```
The snps_file should be a whitespace-delimited file (that can be gzipped) with a header line and at least one of the following two combinations of columns:
1. SNP, A1, A2 - SNP name, reference allele, alternative allele
2. CHR, BP, A1, A2 - chromosome, basepair position (in hg19 coordinates), reference allele, alternative allele

Here is a toy example you can try:
```
python extract_snpvar.py --snps snps_to_finemap.txt.gz --out snps_with_priors
zcat snps_with_priors.snpvar.gz | head
```
The output should be:
```
CHR  BP        SNP                    A1        A2  prior_causal_prob
1    10000006  rs186077422            G         A   1.750133e-05
1    10000179  1:10000179_AAAAAAAC_A  AAAAAAAC  A   1.750133e-05
1    10000400  rs1237370              T         A   1.750133e-05
1    10000476  rs182770070            A         T   1.750133e-05
1    10000553  rs574892739            T         G   1.750133e-05
1    10000732  rs563811805            T         C   1.750133e-05
1    10000804  rs114880362            T         C   1.750133e-05
1    10001239  rs68058227             G         T   1.750133e-05
1    10001401  rs60132751             C         T   1.750133e-05
```
