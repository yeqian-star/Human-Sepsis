#cellphonedbv4分析method2
import pandas as pd
import sys
import os
import gzip
import argparse
pd.set_option('display.max_columns', 100)
#传参
parser = argparse.ArgumentParser(description="manual to this script")
parser.add_argument('--outpath',type=str, default=None)
parser.add_argument('--cpdbpath',type=str, default=None)
parser.add_argument('--celllist',type=str, default=None)
parser.add_argument('--countsfile',type=str, default=None)
parser.add_argument('--threshold',type=float, default=None)
parser.add_argument('--precision',type=int, default=None)
parser.add_argument('--iterations',type=int, default=None)
parser.add_argument('--thread',type=int, default=None)
#parser.add_argument('--outpath',type=str, default=None)
args = parser.parse_args()
#输入文件路径
#cellphone数据库路径
cpdb_file_path = args.cpdbpath
#单细胞分析分群信息
meta_file_path = args.celllist
#单细胞分析count表达量的h5ad
counts_file_path = args.countsfile
#微环境信息，暂时不输入
microenvs_file_path = None
#输出路径设置
out_path = args.outpath
#计算结果精度
precision =args.precision
#表达细胞通讯基因的细胞占比阈值
threshold = args.threshold
#计算p值迭代次数
iterations =args.iterations
thread = args.thread
os.chdir(out_path)
#if os.path.exists('SeuratForCPdb.h5ad'):
	#counts_file_path = 'SeuratForCPdb.h5ad'

#if os.path.exists('counts.txt'):
#	bgzfile_path = counts_file_path
#	with gzip.open(bgzfile_path, 'rb') as f:
#		uncompressed_data = f.read()
#	output_file_path = '/tmp/uncompressed_counts_file.txt'
#	with open(output_file_path, 'wb') as f:
#		f.write(uncompressed_data)
#	counts_file_path = 'counts.txt'



from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
#运行cellphone statistical_analysis
deconvoluted, means, pvalues, significant_means = cpdb_statistical_analysis_method.call(
    cpdb_file_path = cpdb_file_path,                 # mandatory: CellPhoneDB database zip file.
    meta_file_path = meta_file_path,                 # mandatory: tsv file defining barcodes to cell label.
    counts_file_path = counts_file_path,             # mandatory: normalized count matrix.
    counts_data = 'ensembl',                     # defines the gene annotation in counts matrix.
    microenvs_file_path = microenvs_file_path,       # optional (default: None): defines cells per microenvironment.
    iterations = iterations,                               # denotes the number of shufflings performed in the analysis.
    threshold = threshold,                                 # defines the min % of cells expressing a gene for this to be employed in the analysis.
    threads = thread,                                     # number of threads to use in the analysis.
    debug_seed = 42,                                 # debug randome seed. To disable >=0.
    result_precision = precision,                            # Sets the rounding for the mean values in significan_means.
    pvalue = 0.05,                                   # P-value threshold to employ for significance.
    subsampling = False,                             # To enable subsampling the data (geometri sketching).
    subsampling_log = False,                         # (mandatory) enable subsampling log1p for non log-transformed data inputs.
    subsampling_num_pc = 100,                        # Number of componets to subsample via geometric skectching (dafault: 100).
    subsampling_num_cells = 1000,                    # Number of cells to subsample (integer) (default: 1/3 of the dataset).
    separator = '_',                                 # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
    debug = False,                                   # Saves all intermediate tables employed during the analysis in pkl format.
    output_path = out_path,                          # Path to save results.
    output_suffix = 'Novelbio'                       # Replaces the timestamp in the output files by a user defined string in the  (default: None).
    )