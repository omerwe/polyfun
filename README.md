# PolyFun
PolyFun (POLYgenic FUNctionally-informed fine-mapping)

This page contains the code of the method **PolyFun** for functionally-informed fine-mapping, described in [Weissbrod et al. 2019 bioRxiv](https://www.biorxiv.org/content/10.1101/807792v2). PolyFun estimates prior causal probabilities for SNPs, which can then be used by fine-mapping methods like [SuSiE](https://github.com/stephenslab/susieR) or [FINEMAP](http://www.christianbenner.com/). Unlike previous methods for functionally-informed fine-mapping, PolyFun can aggregate polygenic data from across the entire genome and hundreds of functional annotations.

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

If rpy2 or Ckmeans.1d.dp are not installed, PolyFun will fallback to suboptimal clustering via scikit-learn.


We recommend running PolyFun via the [Anaconda Python distribution](https://www.anaconda.com/download/). In Anaconda, you can install all the Python packages with the command "conda install \<package_name\>". Alternatively, the Python packages can be installed with the command "pip install --user \<package_name\>".

Once all the prerequisite packages are installed, you can install PolyFun on a git-enabled machine by typing:
```
git clone https://github.com/omerwe/polyfun
```



<br><br>
# Usage overview
There are three ways to run PolyFun:
1. **Using precomputed prior causal probabilities of 19 million imputed [UK Biobank](https://www.ukbiobank.ac.uk) SNPs with MAF>0.1%, based on a meta-analysis of 15 UK Biobank traits**. This is the simplest approach, but it may not include all your SNPs of interest (especially when analyzing non-European populations) and the prior causal probabilities may not be optimal for some traits.
2. **Computing prior causal probabilities via an L2-regularized extension of [stratified LD-score regression (S-LDSC)](https://www.nature.com/articles/ng.3404)**. This is a relatively simple approach, but the prior causal probabilities may not be robust to modeling misspecification.
3. **Computing prior causal probabilities non-parametrically**. This is the most robust approach, but it is computationally intensive and requires access to individual-level genotypic data from a large reference panel (optimally >10,000 population-matched individuals).

Below are instructions and examples on how to use each approach. We recommend that you run these examples to become familiar with PolyFun. The examples are based on small datasets and run very quickly (typically <1 minute)
<br>
### A note on file formats
PolyFun uses input files that are very similar to [the input files of S-LDSC](https://github.com/bulik/ldsc/wiki/LD-File-Formats). The main differences are:
1. The .annot files **must** contain two addditinal columns called A1,A2 which encode the identifies of the reference and alternative allele
2. The .l2.ldscore files **may** contain the additional columns A1,A2. We strongly encourage including these columns.
3. Polyfun supports files in [.parquet format](https://parquet.apache.org) in addition to .gzip/.bzip2 formats. Parquet files can be loaded substantially faster than alternative formats, at the cost of slightly larger file sizes.

<br><br>

# Approach 1: Using precomputed prior causal probabilities based on a meta-analysis of 15 UK Biobank traits
Here, all you need to do is provide a file with SNP identifiers. PolyFun will extract the prior causal probabilities for these SNPs. To do this, use the following command:
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
The top lines of the output should be:
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


<br><br>

# Approach 2: Computing prior causal probabilities via an L2-regularized extension of S-LDSC
This is done in two stages:

#### 1. Create a munged summary statistics file in a PolyFun-friendly [parquet](https://parquet.apache.org) format.
To do this, use the script `munge_polyfun_sumstats.py`, which takes an input summary statistics file and creates a munged output file. The script tries to be flexible and accomodate multiple file formats and column names. It generally requires only a sample size parameter (n) and a whitespace-delimited input file with SNP rsids, chromosome and basepair info, and either a p-value, an effect size estimate and its standard error, a Z-score or a p-value.

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
In this stage PolyFun will estimate per-SNP heritabilities for SNPs on odd (resp. even) chromosomes by applying L2-regularized S-LDSC to even (resp. odd) chromosomes. To do this, run the scripy `polyfun.py`. This script handles all possible uses of PolyFun, but here we'll only compute prior causal probabilities with L2-extended S-LDSC, using a subset of the [baseline-LF model annotations](https://www.nature.com/articles/s41588-018-0231-8). Here is an example command, that uses 8 annotations from the baseline-LF model:
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
This will create 2 output files for each chromosome: `output/testrun.<CHR>.snpvar_ridge.gz` and `output/testrun.<CHR>.snpvar_ridge_constrained.gz`. The first contains estimated per-SNP heritabilities for all SNPs (which can be used for downstream analysis with PolyFun; see below), and the second contains truncated per-SNP heritabilities, which can be used directly as prior causal probabilities in fine-mapping. For example, here is the output for the top 10 SNPs in chromosome 1: (seen with `zcat output/testrun.1.snpvar_ridge_constrained.gz | head`)
```
CHR  BP      SNP                              A1                    A2  snpvar      Z            N
1    737125  rs151055642                      T                     A   1.3502e-08  4.5924e-01   383290
1    741833  rs148581628                      C                     T   1.3502e-08  -9.4801e-01  383290
1    745642  1:745642_AC_A                    AC                    A   6.5501e-09  5.6848e-02   383290
1    772437  rs138499329                      C                     T   1.3502e-08  8.6788e-01   383290
1    797281  rs76631953                       G                     C   6.5501e-09  -1.2842e+00  383290
1    814300  rs80022136                       T                     A   6.5501e-09  -9.9172e-01  383290
1    821862  1:821862_CACAGCAGCTGTGCTGTGTT_C  CACAGCAGCTGTGCTGTGTT  C   1.3502e-08  3.9271e-02   383290
1    845273  rs117039017                      G                     A   1.3502e-08  5.9879e-01   383290
1    846398  rs58781670                       G                     A   1.3502e-08  2.9464e+00   383290
```
The column called 'snpvar' contains truncated per-SNP heritabilities, which can be used directly as prior causal probabilities in fine-mapping (see below).

The parameters we provided are the following:
1. `--compute-h2-L2` - this tells PolyFun to compute per-SNP heritabilities via an L2-regularized S-LDSC
2. `--no-partitions` - this tells PolyFun to **not** partition SNPs into bins based on their estimated per-SNP heritabilities. You should only provide this flag if you are only interested in L2-regularized estimation of per-SNP heritabilities.
3. `--output-prefix output/testrun` - this specifies the prefix of all the PolyFun output files.
4. `--sumstats` - this specifies an input summary statistics file (created via the `munge_polyfun_sumstats.py` script).
5. `--ref-ld-chr` - this is the prefix of the LD-score and annotation files that S-LDSC uses. These are similar to the standard [S-LDSC  input files](https://github.com/bulik/ldsc/wiki/Partitioned-Heritability) with an important addition: The annotation files **must** include columns called A1,A2 for reference and alternative alleles (because unfortunatley SNP rsid is not a unique SNP identifier). Additionally, it is strongly recommdended that the LD-score files also include columns called A1,A2, to prevent excluding multiple SNPs with the same rsid from the estimation stage. PolyFun will accept files with either .gz or .parquet extension (parquet is faster)
6. `--w-ld-chr` - this is the prefix of the [LD-score weight files](https://github.com/bulik/ldsc/wiki/Partitioned-Heritability), which are generally equal to the regular LD-scores of the SNPs, but restricted to only the set of SNPs used for fitting S-LDSC. As before, it is strongly recommdended that these files include A1,A2 columns.

We strongly encourage that you look at the input files provided in the `example_data` directory to get a sense of their structure.

<br><br>
# Approach 3: Computing prior causal probabilities non-parametrically
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
4. There are various parameters that you can use to control the LD-score computations, analogue to the respective parameters in the [ldsc package](https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial). Please type `python polyfun.py --help` to see all available parameters.
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
This script will output files with re-estimated per-SNP heritabilities that can be used directly for fine-mapping. Here is the output for chromosome 1 (seen via `zcat output/testrun.1.snpvar_constrained.gz | head`):
```
CHR  BP      SNP                              A1                    A2  snpvar      Z            N
1    737125  rs151055642                      T                     A   5.7732e-06  4.5924e-01   383290
1    741833  rs148581628                      C                     T   5.7732e-06  -9.4801e-01  383290
1    745642  1:745642_AC_A                    AC                    A   1.5774e-06  5.6848e-02   383290
1    772437  rs138499329                      C                     T   5.7732e-06  8.6788e-01   383290
1    797281  rs76631953                       G                     C   4.0586e-06  -1.2842e+00  383290
1    814300  rs80022136                       T                     A   4.0586e-06  -9.9172e-01  383290
1    821862  1:821862_CACAGCAGCTGTGCTGTGTT_C  CACAGCAGCTGTGCTGTGTT  C   5.7732e-06  3.9271e-02   383290
1    845273  rs117039017                      G                     A   5.7732e-06  5.9879e-01   383290
1    846398  rs58781670                       G                     A   5.7732e-06  2.9464e+00   383290
```
The snpvar column contains per-SNP heritabilities. These can be used directly as prior causal probabilities in fine-mapping (see below).

<br><br>
# Using prior causal probabilities in fine-mapping
Below we explain how to use the estimated prior causal probabilities with SuSiE and FINEMAP

### Using prior causal probabilities in SuSiE
All you have to do is provide SuSiE the flag **prior_weights** with per-SNP heritability estimates from PolyFun (i.e., the contents of the column `snpvar`).

### Using prior causal probabilities in FINEMAP
This functionality is not implemented yet - please check back soon...

<br><br>
# Using and creating functional annotations
You can either download existing functional annotation files, or create your own:

### Downloading existing functional annotation files
We provide functional annotations for ~19 million UK Biobank imputed SNPs with MAF>0.1%, based on the baseline-LF 2.2.UKB annotations. This is a broad set of coding, conserved, regulatory and LD-related annotations. You can download these annotations and their LD-scores [here](https://data.broadinstitute.org/alkesgroup/LDSCORE/baselineLF_v2.2.UKB.polyfun.tar.gz) (WARNING: this is a large download, requiring 30GB).

### Creating your own annotations
You can easily create your own annotations. The only requirement is to create 22 files (one for each chromosome), each containing columns for CHR, BP, SNP, A1, A2 and arbitrary other columns. These fies can be either .parquet files or gzipped txt files. After creating these files, you should compute LD-scores in each chromosome. You can do this using the script `compute_ldscores.py`. Here is a usage example:
```
mkdir -p output

python compute_ldscores.py \
  --bfile example_data/reference.1 \
  --annot example_data/annotations.1.l2.ldscore.parquet \
  --out output/ldscores_example.parquet
```
This script accepts annotations in either .parquet or plain gzipped text file (parquet is much faster). Please note that you can also use S-LDSC to compute LD-scores. However, S-LDSC does not use the columns A1, A2 in the LD-score and annotation files. Please read the section called "Creating an annot file" in the [S-LDSC wiki](https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial) for more information and instructions.

<br><br>
# FAQ
**Q**: How can I create my own annotations?<br>
**A**: Please read the section called "Creating an annot file" in the [S-LDSC wiki](https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial) for instructions. Note that PolyFun requires adding columns called A1,A2 to uniquely identify SNPs. PolyFun accepts annotation files in either .gzipped text or .parquet format (.parquet is much faster). Please see the `example_data` directory for examples of annotation files. Please note that you will also need to compute LD-scores for your annotations...

<br><br>
# Contact
For questions and comments, please contact Omer Weissbrod at oweissbrod[at]hsph.harvard.edu



