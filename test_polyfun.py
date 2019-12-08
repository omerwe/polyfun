import numpy as np; np.set_printoptions(precision=3, linewidth=200)
import pandas as pd; pd.set_option('display.width', 200) 
import sys
import tempfile
import subprocess
import os
import glob
from pandas.api.types import is_string_dtype
from pandas.api.types import is_numeric_dtype

        
def compare_dfs(dir1, dir2, filename):
    file1 = os.path.join(dir1, filename)
    file2 = os.path.join(dir2, filename)
    if file1.endswith('.parquet'): df1 = pd.read_parquet(file1)
    else: df1 = pd.read_table(file1, delim_whitespace=True)
    if file2.endswith('.parquet'): df2 = pd.read_parquet(file2)
    else: df2 = pd.read_table(file2, delim_whitespace=True)
    assert np.all(df1.shape == df2.shape), 'found dimension mismatch between %s and %s'%(file1, file2)
    assert np.all(df1.columns == df2.columns), 'found mismatch between %s and %s'%(file1, file2)
    for c in df1.columns:
        if is_numeric_dtype(df1[c]):
            assert np.allclose(df1[c], df2[c]), 'found mismatch between %s and %s in column %s'%(file1, file2, c)
        else:
            assert np.all(df1[c] == df2[c]), 'found mismatch between %s and %s in column %s'%(file1, file2, c)

    
    
def test_munge_sumstats(tmpdir):
    script_path = os.path.dirname(os.path.realpath(__file__))
    example_dir = os.path.join(script_path, 'example_data')
    gold_dir = os.path.join(script_path, 'gold')
    script_exe = os.path.join(script_path, 'munge_polyfun_sumstats.py')
    input_file = os.path.join(example_dir, 'boltlmm_sumstats.gz')
    outfile = 'sumstats_munged.parquet'
    output_file = os.path.join(tmpdir, outfile)
    gold_file = os.path.join(gold_dir, outfile)
    
    os.system('python %s --sumstats %s --out %s --n 327209 --min-info 0.6 --min-maf 0.001'%(script_exe, input_file, output_file))
    compare_dfs(tmpdir, gold_dir, outfile)
    
    
def test_polyfun(tmpdir):
    script_path = os.path.dirname(os.path.realpath(__file__))
    example_dir = os.path.join(script_path, 'example_data')
    gold_dir = os.path.join(script_path, 'gold')
    script_exe = os.path.join(script_path, 'polyfun.py')
    sumstats_file = os.path.join(example_dir, 'sumstats.parquet')    
    output_prefix = os.path.join(tmpdir, 'testrun')    
    ref_ld_prefix = os.path.join(example_dir, 'annotations.')  
    w_ld_prefix = os.path.join(example_dir, 'weights.')  
    plink_prefix =  os.path.join(example_dir, 'reference.')
    
    #polyfun stage 2
    os.system('python %s --compute-h2-L2 --output-prefix %s --sumstats %s --ref-ld-chr %s --w-ld-chr %s'%
             (script_exe, output_prefix, sumstats_file, ref_ld_prefix, w_ld_prefix))    
    for outfile in ['testrun.22.bins.parquet', 'testrun.22.snpvar_ridge_constrained.gz', 'testrun.22.snpvar_ridge.gz', 'testrun.22.l2.M']:
        compare_dfs(tmpdir, gold_dir, outfile)
        
    #polyfun stage 3
    os.system('python %s --compute-ldscores --output-prefix %s --bfile-chr %s'%
             (script_exe, output_prefix, plink_prefix))    
    compare_dfs(tmpdir, gold_dir, 'testrun.22.l2.ldscore.parquet')
        
    #polyfun stage 4
    os.system('python %s --compute-h2-bins --output-prefix %s --sumstats %s --w-ld-chr %s'%
             (script_exe, output_prefix, sumstats_file, w_ld_prefix))    
    compare_dfs(tmpdir, gold_dir, 'testrun.22.snpvar_constrained.gz')
    
    
    
def test_polyloc(tmpdir):
    script_path = os.path.dirname(os.path.realpath(__file__))
    example_dir = os.path.join(script_path, 'example_data')
    gold_dir = os.path.join(script_path, 'gold')
    script_exe = os.path.join(script_path, 'polyloc.py')
    sumstats_file = os.path.join(example_dir, 'sumstats2.parquet')
    output_prefix = os.path.join(tmpdir, 'polyloc_test')
    w_ld_prefix = os.path.join(example_dir, 'weights.')
    plink_prefix =  os.path.join(example_dir, 'reference.')
    posterior_file = os.path.join(example_dir, 'posterior_betas.gz')
    
    #polyloc stage 1
    os.system('python %s --compute-partitions --output-prefix %s --posterior %s --bfile-chr %s'%
             (script_exe, output_prefix, posterior_file, plink_prefix))    
    for outfile in ['polyloc_test.22.bins.parquet', 'polyloc_test.22.l2.M']:
        compare_dfs(tmpdir, gold_dir, outfile)
    
    #polyloc stage 2
    os.system('python %s --compute-ldscores --output-prefix %s --bfile-chr %s'%
             (script_exe, output_prefix, plink_prefix))        
    compare_dfs(tmpdir, gold_dir, 'polyloc_test.22.l2.ldscore.parquet')
    
    #polyloc stage 3
    os.system('python %s --compute-polyloc --output-prefix %s --w-ld-chr %s --sumstats %s'%
             (script_exe, output_prefix, w_ld_prefix, sumstats_file))
    for outfile in ['polyloc_test.bin_h2', 'polyloc_test.Mp']:
        compare_dfs(tmpdir, gold_dir, outfile)
        
        

   
    
        
def test_extract_snpvar(tmpdir):
    script_path = os.path.dirname(os.path.realpath(__file__))
    example_dir = os.path.join(script_path, 'example_data')
    gold_dir = os.path.join(script_path, 'gold')
    script_exe = os.path.join(script_path, 'extract_snpvar.py')
    input_file = os.path.join(example_dir, 'snps_to_finemap.txt.gz')
    outfile = 'snps_with_var.gz'
    output_file = os.path.join(tmpdir, outfile)
    gold_file = os.path.join(gold_dir, outfile)    
    
    os.system('python %s --snps %s --out %s'%(script_exe, input_file, output_file))
    compare_dfs(tmpdir, gold_dir, outfile)


def test_finemapper(tmpdir, ldstore_exe):
    script_path = os.path.dirname(os.path.realpath(__file__))
    example_dir = os.path.join(script_path, 'example_data')
    gold_dir = os.path.join(script_path, 'gold')
    script_exe = os.path.join(script_path, 'run_finemapper.py')
    sumstats_file = os.path.join(example_dir, 'chr1.finemap_sumstats.txt.gz')
    plink_file = os.path.join(example_dir, 'chr1')
    outfile = 'finemap.1.46000001.49000001.gz'
    output_file = os.path.join(tmpdir, outfile)
    
    os.system('python %s \
       --geno %s \
       --sumstats %s \
       --n 383290 \
       --chr 1 \
       --start 46000001 \
       --end 49000001 \
       --method susie \
       --max-num-causal 5 \
        --out %s \
        --ldstore %s \
       '
       %(script_exe, plink_file, sumstats_file, output_file, ldstore_exe))
   
    compare_dfs(tmpdir, gold_dir, outfile)    
   
   
    

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('--ldstore', default=None, help='Path to ldstore executable')
    args = parser.parse_args()
    temp_dir = tempfile.mkdtemp()
    
    if args.ldstore is None:
        print('Skipping finemapper tests because --ldstore was not specified')
    else:
        test_finemapper(temp_dir, args.ldstore)
        
    test_polyloc(temp_dir)
    test_polyfun(temp_dir)
    test_extract_snpvar(temp_dir)
    test_munge_sumstats(temp_dir)
    
    
    