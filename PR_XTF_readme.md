## Pull Request Title: Enhancements and New Features for Genetic Data Processing
Code only modified for `finemapper.py`
### Summary of Changes

This PR introduces several updates and new features to improve the handling and processing of genetic data. Key changes include:

1. **Support for BGEN to  cal LD and save to .npz Conversion to accelate the loading speed compared to .bcor:**
    
    - Implemented functionality for converting BGEN files directly to LD matrices stored in .npz format, aligning with polyFun standards. Use `--geno genofile --ldstore2 $(which ldstore) --cache-dir ./ --cache-format npz` to save in npz format by default.
2. **NPZ File Reading Capability:**
    
    - Added capability to read npz files using `--ld your_npz_prefix`.
3. **PGEN File Support:**
	- This may be slower than ldstore2 did!
    - Introduced support for using pgen files as input via `--geno pgen_file_prefix`.
    - Integrated `finemap_tools` for invoking Plink2 (version must be later than PLINK v2.00a6LM 64-bit Intel, dated 2 Mar 2024) to compute LD, using the command template: `plink2 --r2-unphased square`.
    - **Note:** The `--geno` option matches files using prefixes, with bed files having higher priority over pgen to avoid conflicts when both file types are present.
5. **Improvements in LD Matrix Handling:**
    
    - Removed the assertion that LD matrices must not contain any NaN values during reading. Instead, i modified `sync_ld_sumstats` function to exclude SNPs with NA values.
    - Addressed cases where an SNP's LD calculation might result in NaN due to unfiltered genotype data or extreme scenarios. Users are advised to consider stricter QC or automatically exclude SNPs with NA in any LD calculation based on the number of SNPs dropped, as indicated by the output.
5. **Enhancements in Summary Statistics (sumstats) Loading:**
    
    - For sumstats files that are bgz compressed and have an associated .tbi file, implemented reading via the `tabix` command-line tool. This approach is particularly efficient for genome-wide sumstats, allowing direct retrieval of data by chromosome, significantly reducing loading times.
    - In scenarios where `finemap_tools` is unavailable, the original logic of reading the entire file will be followed.
    - Sumstats files organized according to polyFun requirements (columns: SNP, A1, A2, BP, CHR, with SNP in the first column) can be processed using `tabix -s 2 -b 3 -e 3 -c S sumstats_with_bgz_compressed.bgz`.
6. **Integration of `finemap_tools` Package:**
    
    - Included `finemap_tools` for filtering bialleic and ambiguous alleles during sumstats reading.
7. **Code Formatting Updates:**
    
    - Applied code formatting improvements using Black.

### Future Developments

Further development and updates will continue in my own repository and will not be submitted as pull requests to this project.