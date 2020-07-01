import numpy as np; np.set_printoptions(precision=3, linewidth=200)
import pandas as pd; pd.set_option('display.width', 200) 
import sys
import tempfile
import subprocess
import os
import glob
from pandas.api.types import is_string_dtype
from pandas.api.types import is_numeric_dtype

        
def compare_dfs(dir1, dir2, filename, sort_column=None):
    file1 = os.path.join(dir1, filename)
    file2 = os.path.join(dir2, filename)
    if not os.path.exists(file1):
        raise IOError('%s not found'%(file1))
    if not os.path.exists(file2):
        raise IOError('%s not found'%(file2))
    if file1.endswith('.parquet'): df1 = pd.read_parquet(file1)
    else: df1 = pd.read_table(file1, delim_whitespace=True)
    if file2.endswith('.parquet'): df2 = pd.read_parquet(file2)
    else: df2 = pd.read_table(file2, delim_whitespace=True)
    assert np.all(df1.shape == df2.shape), 'found dimension mismatch between %s and %s'%(file1, file2)
    assert np.all(df1.columns == df2.columns), 'found mismatch between %s and %s'%(file1, file2)
    if sort_column is not None:
        df1 = df1.sort_values(sort_column).reset_index(drop=True)
        df2 = df2.sort_values(sort_column).reset_index(drop=True)
    for c in df1.columns:
        if c=='CREDIBLE_SET':
            continue
        if is_numeric_dtype(df1[c]):
            if c=='Mp_STDERR': atol=1e-1
            elif c=='BETA_MEAN': atol=1e-3
            elif c=='BETA_SD': atol=1e-3
            elif c=='PIP': atol=1e-2
            else: atol=1e-4
            if c=='PIP':
                assert np.corrcoef(df1[c], df2[c])[0,1]>0.97, 'found major mismatch between %s and %s in column %s'%(file1, file2, c)
                pip_cutoff = 0.8
                is_large_pip = df2[c] > pip_cutoff
                assert np.allclose(df1.loc[is_large_pip,c], df2.loc[is_large_pip,c], atol=atol), 'found mismatch between %s and %s in column %s'%(file1, file2, c)
            else:
                assert np.allclose(df1[c], df2[c], atol=atol), 'found mismatch between %s and %s in column %s'%(file1, file2, c)
        else:
            assert np.all(df1[c] == df2[c]), 'found mismatch between %s and %s in column %s'%(file1, file2, c)
            
     
          
def is_susie_installed():
    try:
        import rpy2
        import rpy2.robjects.numpy2ri as numpy2ri
        import rpy2.robjects as ro
        ro.conversion.py2ri = numpy2ri
        numpy2ri.activate()
        from rpy2.robjects.packages import importr
        susieR = importr('susieR')
    except:
        return False
    return True

           
def is_Ckmeans_installed():        
    try:
        import rpy2
        import rpy2.robjects.numpy2ri as numpy2ri
        try:
            from importlib import reload
            reload(rpy2.robjects.numpy2ri)
        except:
            pass
        import rpy2.robjects as ro
        ro.conversion.py2ri = numpy2ri
        numpy2ri.activate()
        from rpy2.robjects.packages import importr
        importr('Ckmeans.1d.dp')
        median_seg_func = ro.r('Ckmedian.1d.dp')
        mean_seg_func = ro.r('Ckmeans.1d.dp')            
    except:
        return False
    return True

    
    
def test_munge_sumstats(tmpdir, python3_exe):
    script_path = os.path.dirname(os.path.realpath(__file__))
    example_dir = os.path.join(script_path, 'example_data')
    gold_dir = os.path.join(script_path, 'gold')
    script_exe = os.path.join(script_path, 'munge_polyfun_sumstats.py')
    input_file = os.path.join(example_dir, 'boltlmm_sumstats.gz')
    outfile = 'sumstats_munged.parquet'
    output_file = os.path.join(tmpdir, outfile)
    gold_file = os.path.join(gold_dir, outfile)
    
    retval = os.system('%s %s --sumstats %s --out %s --n 327209 --min-info 0.6 --min-maf 0.001'%(python3_exe, script_exe, input_file, output_file))
    if retval != 0:
        raise ValueError('munge_sumstats command failed')
    
    compare_dfs(tmpdir, gold_dir, outfile)
    
    
def test_polyfun(tmpdir, python3_exe):

    if not is_Ckmeans_installed():
        print('Skipping polyfun tests becuase CkMeans is not installed')
        return


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
    retval = os.system('%s %s --compute-h2-L2 --output-prefix %s --sumstats %s --ref-ld-chr %s --w-ld-chr %s --nnls-exact'%
             (python3_exe, script_exe, output_prefix, sumstats_file, ref_ld_prefix, w_ld_prefix))
    if retval != 0:
        raise ValueError('PolyFun command failed')
    for outfile in ['testrun.22.bins.parquet', 'testrun.22.snpvar_ridge_constrained.gz', 'testrun.22.snpvar_ridge.gz', 'testrun.22.l2.M']:
        compare_dfs(tmpdir, gold_dir, outfile)
        
    #polyfun stage 3
    retval = os.system('%s %s --compute-ldscores --output-prefix %s --bfile-chr %s --nnls-exact'%
             (python3_exe, script_exe, output_prefix, plink_prefix))    
    if retval != 0:
        raise ValueError('PolyFun command failed')
             
    compare_dfs(tmpdir, gold_dir, 'testrun.22.l2.ldscore.parquet')
        
    #polyfun stage 4
    retval = os.system('%s %s --compute-h2-bins --output-prefix %s --sumstats %s --w-ld-chr %s --nnls-exact'%
             (python3_exe, script_exe, output_prefix, sumstats_file, w_ld_prefix))
    if retval != 0:
        raise ValueError('PolyFun command failed')
    compare_dfs(tmpdir, gold_dir, 'testrun.22.snpvar_constrained.gz')
    
    
    
def test_polyloc(tmpdir, python3_exe):

    if not is_Ckmeans_installed():
        print('Skipping polyloc tests becuase CkMeans is not installed')
        return

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
    polyloc_cmd = ('%s %s --compute-partitions --output-prefix %s --posterior %s --bfile-chr %s --nnls-exact'%
                    (python3_exe, script_exe, output_prefix, posterior_file, plink_prefix))
    retval = os.system(polyloc_cmd)
    if retval != 0:
        raise ValueError('PolyLoc command failed')
    for outfile in ['polyloc_test.22.bins.parquet', 'polyloc_test.22.l2.M']:
        compare_dfs(tmpdir, gold_dir, outfile)
        
    
    #polyloc stage 2
    polyloc_cmd = ('%s %s --compute-ldscores --output-prefix %s --bfile-chr %s --nnls-exact'%
                    (python3_exe, script_exe, output_prefix, plink_prefix))
    retval = os.system(polyloc_cmd)
    if retval != 0:
        raise ValueError('PolyLoc command failed')
    compare_dfs(tmpdir, gold_dir, 'polyloc_test.22.l2.ldscore.parquet')
    
    #polyloc stage 3
    polyloc_cmd = ('%s %s --compute-polyloc --output-prefix %s --w-ld-chr %s --sumstats %s --nnls-exact'%
                    (python3_exe, script_exe, output_prefix, w_ld_prefix, sumstats_file))
    retval = os.system(polyloc_cmd)
    if retval != 0:
        raise ValueError('PolyLoc command failed')
    for outfile in ['polyloc_test.bin_h2', 'polyloc_test.Mp']:
        compare_dfs(tmpdir, gold_dir, outfile)
        
        

   
    
        
def test_extract_snpvar(tmpdir, python3_exe):
    script_path = os.path.dirname(os.path.realpath(__file__))
    example_dir = os.path.join(script_path, 'example_data')
    gold_dir = os.path.join(script_path, 'gold')
    script_exe = os.path.join(script_path, 'extract_snpvar.py')
    input_file = os.path.join(example_dir, 'snps_to_finemap.txt.gz')
    outfile = 'snps_with_var.gz'
    output_file = os.path.join(tmpdir, outfile)
    gold_file = os.path.join(gold_dir, outfile)    
    
    retval = os.system('%s %s --sumstats %s --out %s'%(python3_exe, script_exe, input_file, output_file))
    if retval != 0:
        raise ValueError('extract_snpvar command failed')
    compare_dfs(tmpdir, gold_dir, outfile)


def test_finemapper(tmpdir, python3_exe):

    if not is_susie_installed():
        print('Skipping fine-mapping test because either rpy2 or SuSiE are not properly installed. Please try reinstalling rpy2 and/or SuSiE')
        return

    script_path = os.path.dirname(os.path.realpath(__file__))
    example_dir = os.path.join(script_path, 'example_data')
    gold_dir = os.path.join(script_path, 'gold')
    script_exe = os.path.join(script_path, 'run_finemapper.py')
    sumstats_file = os.path.join(example_dir, 'chr1.finemap_sumstats.txt.gz')
    plink_file = os.path.join(example_dir, 'chr1')
    outfile = 'finemap.1.46000001.49000001.gz'
    output_file = os.path.join(tmpdir, outfile)
    
    finemapper_cmd = '%s %s \
       --geno %s \
       --sumstats %s \
       --n 383290 \
       --chr 1 \
       --start 46000001 \
       --end 49000001 \
       --method susie \
       --max-num-causal 5 \
       --out %s \
       --threads 1 \
       --verbose \
       ' \
       %(python3_exe, script_exe, plink_file, sumstats_file, output_file)
    
    #print(finemapper_cmd)
    retval = os.system(finemapper_cmd)
    if retval != 0:
        raise ValueError('finemapper command failed')
    compare_dfs(tmpdir, gold_dir, outfile, sort_column='SNP')
   
   
    

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('--python3', default='python', help='python 3 executable')
    args = parser.parse_args()
    temp_dir = tempfile.mkdtemp()
    
    test_finemapper(temp_dir, args.python3)
        
    test_polyloc(temp_dir, args.python3)
    test_polyfun(temp_dir, args.python3)
    test_extract_snpvar(temp_dir, args.python3)
    test_munge_sumstats(temp_dir, args.python3)
    
    print()
    print()
    print()
    print('---------------------------------------')
    print('All tests completed successfully, hooray!!!')
    print('---------------------------------------')
    print()
    
    
    