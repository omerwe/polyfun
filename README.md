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


We recommend running PolyFun/PolyLoc via the [Anaconda Python distribution](https://www.anaconda.com/download/). In Anaconda, you can install all the Python packages with the command "conda install \<package_name\>". Alternatively, the Python packages can be installed with the command "pip install --user \<package_name\>".

Once all the prerequisite packages are installed, you can install PolyFun/PolyLoc on a git-enabled machine by typing:
```
git clone https://github.com/omerwe/polyfun
```



<br><br>
# Overview of PolyFun
There are three ways to run PolyFun:
1. **Using precomputed prior causal probabilities of 19 million imputed [UK Biobank](https://www.ukbiobank.ac.uk) SNPs with MAF>0.1%, based on a meta-analysis of 15 UK Biobank traits**. This is the simplest approach, but it may not include all your SNPs of interest (especially when analyzing non-European populations) and the prior causal probabilities may not be optimal for some traits.
2. **Computing prior causal probabilities via an L2-regularized extension of [stratified LD-score regression (S-LDSC)](https://www.nature.com/articles/ng.3404)**. This is a relatively simple approach, but the prior causal probabilities may not be robust to modeling misspecification.
3. **Computing prior causal probabilities non-parametrically**. This is the most robust approach, but it is computationally intensive and requires access to individual-level genotypic data from a large reference panel (optimally >10,000 population-matched individuals).

Below are instructions and examples on how to use each approach. We recommend that you run these examples to become familiar with PolyFun. The examples use small datasets and run very quickly (typically <1 minute)
<br>
### A note on file formats
PolyFun uses input files that are very similar to [the input files of S-LDSC](https://github.com/bulik/ldsc/wiki/LD-File-Formats). The main differences are:
1. The .annot files **must** contain two additional columns called A1,A2 which encode the identifies of the reference and alternative allele
2. The .l2.ldscore files **may** contain the additional columns A1,A2. We strongly encourage including these columns.
3. Polyfun supports files in [.parquet format](https://parquet.apache.org) in addition to .gzip/.bzip2 formats. Parquet files can be loaded substantially faster than alternative formats, at the cost of slightly larger file sizes.

<br><br>

## PolyFun approach 1: Using precomputed prior causal probabilities based on a meta-analysis of 15 UK Biobank traits
Here, all you need to do is provide a file with SNP identifiers. PolyFun will extract the per-SNP heritabilities of these SNPs. To do this, use the following command:
```
python extract_snpvar.py --snps <snps_file> --out <output_file>
```
The snps_file should be a whitespace-delimited file (that can be gzipped) with a header line and at least one of the following two combinations of columns:
1. SNP, A1, A2 - SNP name, reference allele, alternative allele
2. CHR, BP, A1, A2 - chromosome, base pair position (in hg19 coordinates), reference allele, alternative allele

Here is a toy example you can try:
```
mkdir -p output
python extract_snpvar.py --snps example_data/snps_to_finemap.txt.gz --out output/snps_with_var.gz
zcat output/snps_with_var.gz | head
```
The top lines of the output should be:
```
CHR  BP        SNP                    A1        A2  SNPVAR
1    10000006  rs186077422            G         A   4.0733e-09
1    10000179  1:10000179_AAAAAAAC_A  AAAAAAAC  A   4.0733e-09
1    10000400  rs1237370              T         A   4.0733e-09
1    10000476  rs182770070            A         T   4.0733e-09
1    10000553  rs574892739            T         G   4.0733e-09
1    10000732  rs563811805            T         C   4.0733e-09
1    10000804  rs114880362            T         C   4.0733e-09
1    10001239  rs68058227             G         T   4.0733e-09
1    10001401  rs60132751             C         T   4.0733e-09

```
The column `SNPVAR` contains the per-SNP heritabilities, which are proportional to prior causal probabilities. These per-SNP heritabilities can be used directly as prior causal probabilities in fine-mapping (see below for details).

<br><br>

## PolyFun approach 2: Computing prior causal probabilities via an L2-regularized extension of S-LDSC
This is done in two stages:

#### 1. Create a munged summary statistics file in a PolyFun-friendly [parquet](https://parquet.apache.org) format.
To do this, use the script `munge_polyfun_sumstats.py`, which takes an input summary statistics file and creates a munged output file. The script tries to be flexible and accommodate multiple file formats and column names. It generally requires only a sample size parameter (n) and a whitespace-delimited input file with SNP rsids, chromosome and base pair info, and either a p-value, an effect size estimate and its standard error, a Z-score or a p-value.

Here is a usage example:
```
python munge_polyfun_sumstats.py \
  --sumstats example_data/boltlmm_sumstats.gz \
  --n 327209 \
  --out example_data/sumstats_munged.parquet \
  --min-info 0.6 \
  --min-maf 0.001
```
This takes the input BOLT-LMM file `example_data/boltlmm_sumstats.gz` and converts it to the parquet file `example_data/sumstats_munged.parquet`, excluding SNPs with INFO score<0.6, with MAF<0.001 or in the MHC region. It will additionally compute the [BOLT-LMM effective sample size](https://www.nature.com/articles/s41588-018-0144-6). You can see other possible arguments with the command `python munge_polyfun_sumstats.py --help`. You can see the output file by opening the parquet file through python with the command `df = pd.read_parquet('example_data/sumstats_munged.parquet')`

#### 2. Run PolyFun with L2-regularized S-LDSC
In this stage PolyFun will estimate per-SNP heritabilities for SNPs on odd (resp. even) chromosomes by applying L2-regularized S-LDSC to even (resp. odd) chromosomes. To do this, run the script `polyfun.py`. This script handles all possible uses of PolyFun, but here we'll only compute prior causal probabilities with L2-extended S-LDSC, using a subset of the [baseline-LF model annotations](https://www.nature.com/articles/s41588-018-0231-8). Here is an example command, that uses 8 annotations from the baseline-LF model:
```
mkdir -p output

python polyfun.py \
    --compute-h2-L2 \
    --no-partitions \
    --output-prefix output/testrun \
    --sumstats example_data/sumstats.parquet \
    --ref-ld-chr example_data/annotations. \
    --w-ld-chr example_data/weights.
```
This will create 2 output files for each chromosome: `output/testrun.<CHR>.snpvar_ridge.gz` and `output/testrun.<CHR>.snpvar_ridge_constrained.gz`. The first contains estimated per-SNP heritabilities for all SNPs (which can be used for downstream analysis with PolyFun; see below), and the second contains truncated per-SNP heritabilities, which can be used directly as prior causal probabilities in fine-mapping. For example, here is the output for the top 10 SNPs in chromosome 1: (seen with `zcat output/testrun.22.snpvar_ridge_constrained.gz | head`)
```
CHR  SNP          BP        A1  A2  SNPVAR      Z            N
22   rs139069276  16866502  G   A   7.1906e-09  -1.2188e-02  383290
22   rs34747326   16870173  A   G   7.1906e-09  -1.8948e+00  383290
22   rs4010550    16900134  G   A   1.4878e-08  1.3657e+00   383290
22   rs5994099    16905044  G   A   7.1906e-09  7.2224e-01   383290
22   rs59750689   16936369  T   C   7.1906e-09  9.4481e-02   383290
22   rs148021587  16939232  C   T   1.4878e-08  8.0397e-01   383290
22   rs3954635    17024983  A   C   1.4878e-08  -3.4335e-01  383290
22   rs371202053  17034394  G   A   1.4878e-08  -7.8644e-01  383290
22   rs200071370  17037779  C   T   1.4878e-08  3.5049e-01   383290
```
The column called 'SNPVAR' contains truncated per-SNP heritabilities, which can be used directly as prior causal probabilities in fine-mapping (see below).

The parameters we provided are the following:
1. `--compute-h2-L2` - this tells PolyFun to compute per-SNP heritabilities via an L2-regularized S-LDSC
2. `--no-partitions` - this tells PolyFun to **not** partition SNPs into bins based on their estimated per-SNP heritabilities. This makes the computations slightly faster. You should only provide this flag if you are only interested in L2-regularized estimation of per-SNP heritabilities.
3. `--output-prefix output/testrun` - this specifies the prefix of all the PolyFun output files.
4. `--sumstats` - this specifies an input summary statistics file (created via the `munge_polyfun_sumstats.py` script).
5. `--ref-ld-chr` - this is the prefix of the LD-score and annotation files that S-LDSC uses. These are similar to the standard [S-LDSC  input files](https://github.com/bulik/ldsc/wiki/Partitioned-Heritability) with an important addition: The annotation files **must** include columns called A1,A2 for reference and alternative alleles (because unfortunatley SNP rsid is not a unique SNP identifier). Additionally, it is strongly recommdended that the LD-score files also include columns called A1,A2, to prevent excluding multiple SNPs with the same rsid from the estimation stage. PolyFun will accept files with either .gz or .parquet extension (parquet is faster)
6. `--w-ld-chr` - this is the prefix of the [LD-score weight files](https://github.com/bulik/ldsc/wiki/Partitioned-Heritability). These files should be similar to standard LD-score files, but with only one 'annotation' that includes weights. These weights should be equal to the (non-stratified) LD-scores of the SNPs, when computed using plink files that include only the set of SNPs used for fitting S-LDSC. As before, it is strongly recommdended that these files include A1,A2 columns.

We strongly encourage that you look at the input files provided in the `example_data` directory to get a sense of their structure.

<br><br>
## PolyFun approach 3: Computing prior causal probabilities non-parametrically
This is done in four stages:

#### 1. Create a munged summary statistics file in a PolyFun-friendly [parquet](https://parquet.apache.org) format.
This is done exactly as in step 1 of Approach 2 (see above).

#### 2. Run PolyFun with L2-regularized S-LDSC
In this stage PolyFun will estimate per-SNP heritabilities for SNPs on odd (resp. even) chromosomes by applying L2-regularized S-LDSC to even (resp. odd) chromosomes, and will then partition the SNPs into bins. This is done similarly to step 2 of Approach 2 (see above), except that you should **remove** the `--no-partitions` flag. There are several additional flags that you can specify:
1. `--skip-Ckmedian` - This tells PolyFun to partition SNPs into bins using scikits-learn instead of Ckmedian. This is a suboptimal clustering, so Ckmedian is preferable. You should only specify this argument if rpy2 and/or Ckmeans.1d.dp are not installed on your machine or you can't get them to run.
2. `-num-bins <K>` - this specifies the number of bins to partition SNPs into. By default PolyFun will try to estimate this number. You should specify this number if either (a) you specified `--skip-Ckmedian` (because scikits-learn cannot estimate the number of bins) or (b) the estimation is too slow.

Here is a usage example:
```
mkdir -p output

python polyfun.py \
    --compute-h2-L2 \
    --output-prefix output/testrun \
    --sumstats example_data/sumstats.parquet \
    --ref-ld-chr example_data/annotations. \
    --w-ld-chr example_data/weights.
```


#### 3. Compute LD-scores for each SNP bin
In this stage PolyFun will compute LD-scores for each SNP bin. This is the most computationally intensive part of PolyFun. PolyFun can compute LD-scores for all chromosomes in the same run, or for only one chromosome at a time. The second option is recommended if you can run multiple jobs on a cluster. Here is a usage example that computes LD-scores for only chromosome 1:
```
python polyfun.py \
    --compute-ldscores \
    --output-prefix output/testrun \
    --bfile-chr example_data/reference. \
    --chr 1
```
###### Please note the following:
1. You must specify the same `--output-prefix` argument that you provided in stage 2, because PolyFun requires intermediate files that were created in stage 2.
2. If you remove the flag `--chr `, PolyFun will iterate over all chromosomes and compute LD-scores for all of them
3. This stage requires individual-level genotypic data from a large reference panel that is population-matched to your study. Ideally this data should come from your study directly. In this example we used a small subset of SNPs of European-ancestry individuals from the [1000 genomes project](https://www.internationalgenome.org).
4. There are various parameters that you can use to control the LD-score computations, analogue to the respective parameters in the [ldsc package](https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial). Please type `python polyfun.py --help` to see all available parameters. The parameter `--keep <keep file>` can be especially useful if you have a very large reference panel and would like to speed-up the computations by using only a subset of individuals.
5. You can run stages 2 and 3 together by invoking `polyfun.py` with both the flags `--compute-h2-L2` and `--compute-ldscores`.

#### 4. Re-estimate per-SNP heritabilities via S-LDSC 
In this stage PolyFun will re-estimate per-SNP heritabilities robustly  by estimating the heritability causally explained by each bin, using only the chromosomes that were not used to partition SNPs into bins. Technically, PolyFun  will estimate the heritability causally explained by SNPs in each bin that are in an even (resp. odd) target chromosome, by applying S-LDSC with non-negativity constraints to even (resp. odd) chromosomes, excluding the target chromosome, and dividing the bin heritability by the number of SNPs in the bin. You can only run this stage after computing LD-scores for all chromosomes. Here is a usage example:
```
python polyfun.py \
    --compute-h2-bins \
    --output-prefix output/testrun \
    --sumstats example_data/sumstats.parquet \
    --w-ld-chr example_data/weights.
```
This script will output files with re-estimated per-SNP heritabilities that can be used directly for fine-mapping. Here is the output for chromosome 1 (seen via `zcat output/testrun.22.snpvar_constrained.gz | head`):
```
CHR  SNP          BP        A1  A2  SNPVAR      Z            N
22   rs139069276  16866502  G   A   4.5173e-06  -1.2188e-02  383290
22   rs34747326   16870173  A   G   3.5684e-06  -1.8948e+00  383290
22   rs4010550    16900134  G   A   4.9595e-06  1.3657e+00   383290
22   rs5994099    16905044  G   A   3.5684e-06  7.2224e-01   383290
22   rs59750689   16936369  T   C   4.5173e-06  9.4481e-02   383290
22   rs148021587  16939232  C   T   4.9595e-06  8.0397e-01   383290
22   rs3954635    17024983  A   C   4.9595e-06  -3.4335e-01  383290
22   rs371202053  17034394  G   A   4.9595e-06  -7.8644e-01  383290
22   rs200071370  17037779  C   T   4.9595e-06  3.5049e-01   383290
```
The `SNPVAR` column contains per-SNP heritabilities. These can be used directly as prior causal probabilities in fine-mapping (see below).


<br><br>
# Using prior causal probabilities in fine-mapping
Below we explain how to use the estimated prior causal probabilities with SuSiE and FINEMAP. We recommend using the script **run_finemapper.py**, which saves many of the preprocessing steps requires to perform fine-mapping (e.g. handling allelic flips between the summary statistics and reference genotypes). Alternatively, you can run SuSiE or FINEMAP directly with the prior causal probabilities computed by PolyFun, as described below.

## Using prior causal probabilities using the run_finemapper script
The script `run_finemapper` takes an input a file with summary statistics and a file with genotypes from a reference panel, and performs functionally-informed fine-mapping using methods like SuSiE or FINEMAP. It works seamlessly with PolyFun by taking input files created by `munge_polyfun_sumstats.py` or by PolyFun itself. `run_finemapper` computes an LD matrix using [LDstore](http://www.christianbenner.com), which must be installed on your system. `run_finemapper` can cache LD matrices on disk, which can save substantial time and effort when re-analyzing the same data multiple times with different configurations (which always happens).

To run `run_finemapper` with SuSiE, you need to install [rpy2](https://rpy2.bitbucket.io/) and [the SuSiE package](https://github.com/stephenslab/susieR) on your system. To run it with FINEMAP, you need to install [the FINEMAP software](http://www.christianbenner.com) on your system. 

We will first show a use example and then describe all the command line arguments:
```
mkdir -p LD_cache
mkdir -o output

python run_finemapper.py \
    --geno example_data/chr1 \
    --sumstats example_data/chr1.finemap_sumstats.txt.gz \
    --n 383290 \
    --chr 1 \
    --start 46000000 \
    --end 49000000 \
    --method susie \
    --max-num-causal 5 \
    --cache-dir LD_cache \
    --out output/finemap.1.46000000.49000000.gz \
    --ldstore <PATH_TO_LDSTORE_EXECUTABLE>
```
This command takes an input plink file (`example_data/chr1.bed`) and an input summary statistics file (corresponding to an analysis of n=383,290 individuals). The script performs fine-mapping in chromosome 1, in the locus spanning basepair positions 46000000-49000000, using SuSiE. It prints the output to the file `output/finemap.1.46000000.49000000.gz`, and saves the computed LD matrix in the directory `LD_cache`. The argument `--max-num-causal 5` tells SuSiE to assume that there are exactly 5 causal SNPs in the locus (the argument name is general, but for SuSiE it specifies an exact rather than a max number). Here are the first few lines of the output (seen with `zcat output/finemap.1.35000000.55000000.gz | head`):
```
CHR  SNP         BP        A1  A2  SNPVAR       Z             N       P            PIP          BETA_MEAN     BETA_SD      CREDIBLE_SET
1    rs2088102   46032974  T   C   1.70060e-06  1.25500e+01   383290  3.97510e-36  1.00000e+00  1.91941e-02   1.60994e-03  1
1    rs7528714   47966058  G   A   1.18040e-06  5.14320e+00   383290  2.70098e-07  1.00000e+00  1.84424e-02   1.60920e-03  4
1    rs870357    47912628  T   C   1.70060e-06  -3.85430e+00  383290  1.16061e-04  1.00000e+00  -1.61992e-02  1.60781e-03  3
1    rs7528075   47870271  G   A   1.18040e-06  4.40160e+00   383290  1.07456e-05  9.90661e-01  6.54552e-03   1.68805e-03  5
1    rs212968    48734666  G   A   1.70060e-06  -3.01130e+00  383290  2.60132e-03  9.51385e-01  -5.53938e-03  1.97850e-03  2
1    rs4511165   48293181  G   A   1.70060e-06  -1.18940e+00  383290  2.34282e-01  3.09750e-02  -1.24825e-04  7.52260e-04  0
1    rs3766196   47284526  C   A   6.93040e-06  -5.92360e-02  383290  9.52764e-01  8.78914e-03  1.51326e-05   2.18600e-04  0
1    rs12567716  48197570  T   C   1.18040e-06  2.14810e+00   383290  3.17058e-02  7.28725e-03  2.39232e-05   3.09911e-04  0
1    rs4927234   48384796  G   A   1.70060e-06  -7.59600e-01  383290  4.47494e-01  1.48081e-03  -1.52787e-06  7.26440e-05  0
```
Columns 1-9 describe the input summary statistics (and are based on data from the input files). The rows are sorted according to PIP in descending order. Columns 10-14 contain the following fields:
1. **PIP** - posterior causal probability
2. **BETA_MEAN** - posterior mean of causal effect size
3. **BETA_SD** - posterior standard deviation of causal effect size
4. **CREDIBLE_SET** - the index of the first (typically smallest) credible set that the SNP belongs to (0 means none).

We now describe the command-lime arguments of `run_finemapper` in detail:
1. **--geno** - The name of a .bgen file, or the *prefix* of the name of a plink file (without the suffix .bed). `run_finemapper` will compute an LD matrix using the genotypes in this file. **Warning: this file should ideally contain the same individuals used to generate summary statistics, or at least very closely matched individuals. Using an external reference panel in fine-mapping is strongly discouraged and can lead to severe false-positive results** (see [Benner et al. 2017 AJHG](https://www.cell.com/ajhg/fulltext/S0002-9297(17)30334-8), [Ulirsch et al. 2019 Nat Genet](https://www.nature.com/articles/s41588-019-0362-6) for an investigation of this issue).
2. **--sumstats** - The name of a summary statistics file, which must include the columns `SNP`, `CHR`, `BP`, `A1`, `A2`, `Z` (z-score). This file can also include a column called `SNPVAR` that specifies prior per-SNP heritability. If it exists (and unless requested otherwise), `run_finemapper` will use this column to perform functionally-informed fine-mapping. We recommend using the output files of PolyFun as input sumstats files for `run_finemapper`.
3. **--n** - the sample size used to generate summary statistics. In case the summary statistics were computed with [BOLT-LMM](https://data.broadinstitute.org/alkesgroup/BOLT-LMM), we recommend specifying the [effective sample size](https://www.nature.com/articles/s41588-018-0144-6) (this quantity is automatically computed by `munge_polyfun_sumstats.py`).
4. **--chr** - the target chromosome to fine-map.
5. **--start**, **--end** - the start and end positions of the target locus to finemap (base pair coordinates).
6. **--method** - the fine-mapping method. `run_finemapper` currently supports only `--method susie`, but FINEMAP support is coming soon.
7. **--max-num-causal** - the max number of causal SNPs that can be modeled (for FINEMAP) or the exact number (for SuSiE).
8. **--cache-dir** - a directory that will cache LD matrices for future reuse. If not specified, LD matrices will be saved to a temp directory and deleted after the script terminates.
9. **--out** - the name of the output file.
10. **--ldstore** - the path of the [LDstore](http://www.christianbenner.com) executable on your system.
11. **--non-funct** - if specified, `run_finemapper` will perform non-functionally informed fine-mapping. In this case it does not require that the sumstats file has a column called `SNPVAR`.
12. **--verbose** - if specified, `run_finemapper` and the tools it runs will print more detailed output messages.
13. **--hess** - if specified, the prior causal effect size variance will be determined using a modified [HESS](https://www.sciencedirect.com/science/article/pii/S0002929716301483) procedure, as described in the PolyFun paper. Otherwise, the causal effect size variance will be estimated by SuSiE and/or FINEMAP.
14. **--sample-file** - if you provide a bgen file for `--geno`, you must also provide a SNPTEST2 sample file, like for other tools that use bgen data. This is a simple text file without a header and a single column that contains individual ids. If `--geno` is a plink file, you do not need to provide this argument. Please see the [LDstore website](http://www.christianbenner.com) for more information and examples.
15. **--incl-samples** - an optional text file without a header and with a single column that includes the ids of a subset of individuals to use in the LD matrix computation. This can be provided for both plink and bgen files. If not provided, all individuals will be used. Please see the [LDstore website](http://www.christianbenner.com) for more information and examples.
16. **--threads** - the number of CPU threads that LDstore will use to compute LD matrices (if not specified, use the max number of available CPU cores).

## Using prior causal probabilities in SuSiE directly
All you have to do is provide SuSiE the flag **prior_weights** with per-SNP heritability estimates from PolyFun (i.e., the contents of the column `SNPVAR`).

## Using prior causal probabilities in FINEMAP directly
This functionality is not implemented yet - please check back soon...

<br><br>
# Using and creating functional annotations
You can either download existing functional annotation files, or create your own:

### Downloading existing functional annotation files
We provide [functional annotations for ~19 million UK Biobank imputed SNPs with MAF>0.1%, based on the baseline-LF 2.2.UKB annotations](https://data.broadinstitute.org/alkesgroup/LDSCORE/baselineLF_v2.2.UKB.polyfun.tar.gz) **(WARNING: this is a large download, requiring 30GB)**. This is a broad set of coding, conserved, regulatory and LD-related annotations, based on [Gazal et al. 2018 Nat Genet](https://www.nature.com/articles/s41588-018-0231-8) and
described in Supplementary Table 1 of [Weissbrod et al. 2019 bioRxiv](https://www.biorxiv.org/content/10.1101/807792v2).

### Creating your own annotations
You can easily create your own annotations. The only requirement is to create 22 files (one for each chromosome), each containing columns for CHR, BP, SNP, A1, A2 and arbitrary other columns representing your annotations. These files can be either .parquet or .gz files (we recommend using .parquet files).
<br>
To see an example file, type the following commands from within python:
```
import pandas as pd
df = pd.read_parquet('example_data/annotations.22.annot.parquet')
print(df.head())
```
The output should be:
```
           SNP  CHR        BP A1 A2  ...  Conserved_LindbladToh_common  Conserved_LindbladToh_lowfreq  Repressed_Hoffman_common  Repressed_Hoffman_lowfreq  base
0  rs139069276   22  16866502  G  A  ...                             0                              0                         1                          0     1
1   rs34747326   22  16870173  A  G  ...                             0                              0                         0                          1     1
2    rs4010550   22  16900134  G  A  ...                             0                              0                         0                          0     1
3    rs5994099   22  16905044  G  A  ...                             0                              0                         0                          1     1
4   rs59750689   22  16936369  T  C  ...                             0                              0                         1                          0     1
```


After creating these files, you should compute LD-scores for these annotations (one file for each chromosome). You can do this using the script `compute_ldscores.py`. Here is a use example for chromosome 1:
```
mkdir -p output

python compute_ldscores.py \
  --bfile example_data/reference.1 \
  --annot example_data/annotations.1.l2.ldscore.parquet \
  --out output/ldscores_example.parquet
```
Here, `--bfile` is the prefix of a plink .bed file of a reference panel with chromosome 1 SNPs, `--annot` is the name of an annotations file, and `--out` is the name of an output file.  The script also accepts a `--keep <keep file>` parameter to use a subset of individuals for faster computation. This script accepts annotations in either .parquet or .gz format (parquet is much faster). Please note that you can also use S-LDSC to compute LD-scores. However, S-LDSC requires python 2 and does not use the columns A1, A2 in the LD-score and annotation files.


<br><br>
# Overview of PolyLoc
PolyLoc takes an input file with posterior means and standard deviations of causal effect sizes (estimated by `run_finemapper`). PolyLoc uses this file to partition SNPs into bins of similar posterior per-SNP heritability, and estimates the heritability causally explained by each bin. PolyLoc consists of three stages:
1. Partition SNPs into bins of similar posterior per-SNP heritability
2. Compute LD-scores for the SNP bins
3. Estimate the heritability casaully explained by each bin. **This stage requires summary statistics based on different data than the data used to run PolyFun** (see [Weissbrod et al. 2019 bioRxiv for details](https://www.biorxiv.org/content/10.1101/807792v2)).
<br>
PolyLoc and PolyFun have similar input files and they share many command-line arguments. You can see all the PolyLoc options by typing `python polyloc.py --help`. We now describe each of the stages of PolyLoc in detail.

#### PolyLoc stage 1: Partition SNPs into bins of similar posterior per-SNP heritability
This stage requires a file with posterior causal effect sizes of SNPs (ideally all genome-wide SNPs, or at least all the ones in genome-wide significant loci), as created by `run_finemapper.py`. Here is an example file (seen via `zcat example_data/posterior_betas.gz | head`):
```
CHR  SNP         BP       A1  A2  Z         N       BETA_MEAN  BETA_SD
1    rs3748597   888659   T   C   0.85233   383290  0.00021    0.00005
1    rs6661956   2465912  C   T   0.98513   383290  0.00161    0.00006
1    rs2376821   2974852  T   C   -0.73402  383290  -0.00186   0.00111
1    rs10797386  3168280  A   G   -0.59635  383290  -0.00057   0.00006
1    rs7548372   4943523  C   A   2.60293   383290  0.00643    0.00014
1    rs6678838   5003147  G   T   1.10294   383290  0.00135    0.00008
1    rs6675995   5226752  G   A   -1.22154  383290  -0.00027   0.00034
1    rs2012852   5541205  C   A   0.30166   383290  0.00022    0.00027
1    rs10864271  7135268  A   G   2.03717   383290  0.00471    0.00031
```
The column `BETA_MEAN` contains the posterior means of causal effect sizes (as estimated by PolyFun), and the column `BETA_SD` contains their posterior standard deviation. The other required columns are `CHR`, `SNP`, `BP`, `A1`, `A2`.

Here is an example command that uses this file:
```
mkdir -p output

python polyloc.py \
    --compute-partitions \
    --output-prefix output/polyloc_test \
    --posterior example_data/posterior_betas.gz \
    --bfile-chr example_data/reference. 
```
The parameters we provided are the following:
1. `--compute-partitions` - this tells PolyLoc to partition SNPs into bins of similar posterior per-SNP heritability
2. `--output-prefix` - this specifies the prefix of all PolyLoc output file names
3. `--posterior` - this specifies the input file with posterior means and standard deviations of causal effect sizes. This file should ideally include information for all SNPs. Every SNP not included in this file will be treated as if its posterior causal effect size is approximately zero, potentially leading to suboptimal polygenic localization if it's an important SNP. As in PolyFun, this file can either be a text file (possibly gzipped) or a parquet file, which allows faster loading.
4. `--bfile-chr` - this is the prefix of plink files that PolyLoc will use to assign SNPs not reported in the `--posterior` file into bins. As in PolyFun, there must be one plink file for each chromosome.

<br>

Additional optional parameters are:
1. `--skip-Ckmedian` - This tells PolyLoc to partition SNPs into bins using scikits-learn instead of Ckmedian. This is a suboptimal clustering, so Ckmedian is preferable. You should only specify this argument if rpy2 and/or Ckmeans.1d.dp are not installed on your machine or you can't get them to run.
2. `-num-bins <K>` - this specifies the number of bins to partition SNPs into. By default PolyLoc will try to estimate this number. You should specify this number if either (a) you specified `--skip-Ckmedian` (because scikits-learn cannot estimate the number of bins) or (b) the estimation is too slow.

###### Other notes and comments:
1. If you use `run_finemapper` to fine-map multiple loci, you can concatenate the output files of these loci (making sure there aren't duplicate SNPs that appear twice) and provide the concatenated file as an input posterior file.

## PolyLoc stage 2: Compute LD-scores for the SNP bins
This stage is similar to the LD-score computation stage in PolyFun. In this stage PolyLoc will compute LD-scores for each SNP bin. This is the most computationally intensive part of PolyLoc. PolyLoc can compute LD-scores for all chromosomes in the same run, or for only one chromosome at a time. The second option is recommended if you can run multiple jobs on a cluster. Here is an example that computes LD-scores for only chromosome 1:
```
python polyloc.py \
    --output-prefix output/polyloc_test \
    --compute-ldscores \
    --bfile-chr example_data/reference. \
    --chr 1
```
###### Please note the following:
1. PolyLoc accepts the same parameters as PolyFun for LD-scores computations.
2. You must specify the same `--output-prefix` argument that you provided in stage 1, because PolyLoc requires intermediate files that were created in stage 1.
3. If you remove the flag `--chr`, PolyLoc will iterate over all chromosomes and compute LD-scores for all of them, which may take a long time.
4. This stage requires individual-level genotypic data from a large reference panel that is population-matched to your study. Ideally this data should come from your study directly. In this example we used a small subset of SNPs of European-ancestry individuals from the [1000 genomes project](https://www.internationalgenome.org).
5. There are various parameters that you can use to control the LD-score computations, analogue to the respective parameters in the [ldsc package](https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial). Please type `python polyloc.py --help` to see all available parameters. The parameter `--keep <keep file>` can be especially useful if you have a very large reference panel and would like to speed-up the computations by using only a subset of individuals.
6. You can run stages 1 and 2 together by invoking `polyloc.py` with both the flags `--compute-partitions` and `--compute-ldscores`.


## PolyLoc stage 3: Estimate the heritability casaully explained by each bin
This stage requires LD-score weights files and a summary statistics file **that is different from the one used to compute posterior causal effect sizes (to prevent biased estimates due to winner's curse)**. Here is an example command:
```
python polyloc.py \
    --output-prefix output/polyloc_test \
    --compute-polyloc \
    --w-ld-chr example_data/weights. \
    --sumstats example_data/sumstats2.parquet
```
The output of this command is a polygenic localization table. Here is the output of this example, which you can see by typing `cat output/polyloc_test.polyloc`:
```
BIN  BIN_SIZE  %H2      SUM_%H2
1    10        0.17250  0.17250
2    32        0.13510  0.30760
3    70        0.12648  0.43408
4    110       0.13666  0.57074
5    193       0.15129  0.72203
6    255       0.09629  0.81832
7    315       0.07411  0.89243
8    410       0.05337  0.94580
9    492       0.05057  0.99637
10   639       0.00363  1.00000
11   949       0.00000  1.00000
```
The output shows that PolyLoc partitioned SNPs into 11 bins of similar posterior per-SNP heritability. The first bin includes 10 SNPs that jointly explain 17.25% of the total SNP heritability, the second bin includes 32 SNPs that jointly explain 13.5% of the total SNP heritability, and so on. The identities of the SNPs in each bin are the SNPs in the posterior effect sizes file, ranked according to their posterior per-SNP heritability estimates (i.e., the sum of their squared posterior mean and their squared posterior standard deviation). That is, the SNPs in bin 1 are the 10 SNPs with the largest posterior per-SNP heritability, the SNPs in bin 2 are the next top ranked 32 SNPs, and so on.

As in PolyFun, you can run jointly multiple stages of PolyLoc by using several mode parameters (e.g. `python polyloc.py ----compute-partitions --compute-ldscores --compute-polyloc`).

## General notes about PolyLoc
1. PolyLoc estimates per-bin heritability with respect to all SNPs in the plink files. If you would like to exclude to e.g. only common SNPs (as done in the PolyFun paper), please create plink files with only common SNPs and subset your posterior file accordingly (no need to subset the sumstats file).


<br><br>
# FAQ
**Q**: Should I create a base annotations that includes only the number 1 for all SNPs?
<br>
**A**: Typically yes. However, in some cases the LD-scores for this annotation may be linearly dependent on the LD-scores of your other annotations, in which case you don't need to create this annotation. This can happen if (1) the vector of ones [1.0 1.0 ... 1.0] is linearly dependent on your other annotations (which holds for the baseline-LF annotations); and (2) The LD-score that you compute for each SNP is based on (almost) exactly the same set of SNPs as your set of summary statistics SNPs. Hence, we did not include a base annotation in our version of the baseline-LF annotations.
<br>
<br>
**Q**: Can I add extra annotations on top of the baseline-LF annotations, without creating huge new files from scratch?
<br>
**A**: Yes. The flag `--ref-ld-chr` accepts a comma-separated list of file name prefixes, just like standard S-LDSC. For example, you can create a set of annotation files called my_annot.1.annot.parquet, ... my_annot.22.annot.parquet, and then invoke polyfun as follows:
```
python polyfun.py \
    --compute-h2-L2 \
    --output-prefix output/testrun \
    --sumstats example_data/sumstats.parquet \
    --ref-ld-chr example_data/annotations.,my_annot. \
    --w-ld-chr example_data/weights.
```

<br><br>
# Contact
For questions and comments, please open a Github issue (preferred) or contact Omer Weissbrod at oweissbrod[at]hsph.harvard.edu



