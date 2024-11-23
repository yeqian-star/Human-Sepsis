# 20230505 决定把结果写在容器根目录./下 但是问题在于 没办法到summary中压缩总的cc/ci文件夹
# 所以要在每个prefix的分容器中 只保留prefix命名的文件夹 再到xml中新建cc/ci文件夹保存
# 最终输出结果的时候 几个容器中的cc/ci路径会被合并 这样就可以再到summary中压缩cc/ci总文件夹
# 结果中  cc/ci.zip解压完是分prefix的压缩包
# 20230504 由于task图片太多 需要把除了all文件夹以外的都压缩起来 最终效果是
#  最终效果是 cc.zip + cc/prefix/prefix.png/Cluster/AllCluster 
# 就需要把结果先写在tmp里 然后整体压缩并cp出来 但是dat文件和csv还是要写在output路径 所以要添加一个输出路径参数 
#20230106  pdf打包 png展示  文件太多就不输出png了   文件太多文件夹会打不开 
#20221216 task目标是 展示需要根据cluster自定义的图 并且加参数是否画分样本的图 根据不同的逻辑分文件夹输出

#20240708 添加一个小需求
# 由于重分析需要重新aggr,如果新aggr和单细胞最初的aggr的输入的样本顺序不相同,那么细胞名-后面的数字就会不同;
# 这个顺序调整由技术部自己负责,本次只需要添加:如果单细胞和tcr的细胞交集少于30% 则task不报错但返回蓝色;
# 如果交集小于10%,则task不报错但返回紫色;返回蓝色是223,返回紫色是233

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import numpy as np
import pandas as pd
import scanpy as sc
import scirpy as ir
from matplotlib import pyplot as plt, cm as mpl_cm
import cycler
import argparse
import os
import pickle
import logging
import skbio

from pandas.api.types import CategoricalDtype
from typing import cast

import sys

sc.set_figure_params(figsize=(10, 10))
sc.settings.verbosity = 2  
#打印使用到的软件包
sc.logging.print_header()

parser = argparse.ArgumentParser(description='manual to this script')

parser.add_argument('--path', type=str, default=None)#容器当前目录 ./路径 本质和tmp没区别 结果都要cp出来
parser.add_argument('--outpath', type=str, default=None)#输出路径 用于保存dat文件

parser.add_argument('--prefix', type=str, default=None)#prefix
parser.add_argument('--seuseth5', type=str, default=None)#包含cluster信息的单细胞adata
parser.add_argument('--adatair', type=str, default = None)#norm生成的adatair

parser.add_argument('--min_nodes2', type=int, default = None)
parser.add_argument('--size_power2', type=float, default = None)
parser.add_argument('--min_cells2', type=int, default = None)
parser.add_argument('--base_size2', type=int, default = None)#可能为空值
parser.add_argument('--edges_width2', type=float, default = None)
parser.add_argument('--label_fontsize2', type=int, default = None)
parser.add_argument('--label_alpha2', type=float, default = None)

parser.add_argument('--clusterorder', type=str, default=None)#单列绘图时cluster顺序表
parser.add_argument('--sampleorder', type=str, default = None)#单列绘图时sample顺序表

args = parser.parse_args()

path1 = args.path
outpath1 = args.outpath

prefix = args.prefix
seuseth5 = args.seuseth5
adatair = args.adatair

min_nodes2 = args.min_nodes2
size_power2 = args.size_power2
min_cells2 = args.min_cells2
base_size2 = args.base_size2
edges_width2 = args.edges_width2
label_fontsize2 = args.label_fontsize2
label_alpha2 = args.label_alpha2

clusterorder = args.clusterorder
sampleorder = args.sampleorder

##################################################################新建文件夹

# prefix = "CD4T"
# adatair = "/media/nbc1/lijialun/scirpy/zhenbo/TCRNrom/cc/adatair_cc.dat"
# path1 = "/media/nbc1/lijialun/scirpy/zhenbo/TCRCluster/cc/CD4T"  
# seuseth5 = "/media/nbc1/lijialun/scirpy/zhenbo/inputfile/CD4T.seuratsub.h5ad"

#20230112 文件夹结构 

# -prefix
# --prefix.png
# ---Cluster
# ----AllCluster
# -----QC_Graph
# -----clonotype_cluster_network
# -----clonal_expansion
# ------umapinfo
# ------clonal_expansion_stackedbarplot
# ------clonotype_diversity
# ------clonetpye_size_abundance
# ----Cluster1
# -----QC_Graph
# -----clonotype_cluster_network
# -----clonal_expansion
# ------clonal_expansion_stackedbarplot
# ------clonetpye_size_abundance
# ---Sample
# ---AllSample#类比Cluster结构
# -----QC_Graph
# -----clonotype_cluster_network
# -----clonal_expansion
# ------umapinfo
# ------clonal_expansion_stackedbarplot
# ------clonotype_diversity
# ------clonetpye_size_abundance
# ---Sample1#类比Cluster1结构
# -----QC_Graph
# -----clonotype_cluster_network
# -----clonal_expansion
# ------clonal_expansion_stackedbarplot
# ------clonetpye_size_abundance
# --prefix.pdf#类比png结构 

#新建总文件夹 
www = "/"
path = path1

#20230505 文件夹结构改动 不再带上最外层的cc文件夹
# os.mkdir(path1+www+"cc")
# path1 = path1+www+"cc"+www

os.mkdir(path1+www+prefix)
link = path1+www+prefix+www

#设计一个以rds_prefix.png命名的只存放png展示的文件夹 
link1 = link+ prefix+".png"+www
os.mkdir(link1)

#设计一个以rds_prefix.pdf命名的只存放png展示的文件夹 
link2 = link+ prefix+".pdf"+www
os.mkdir(link2)

####存放cluster信息的总文件夹 
#rds_prefix.png/Cluster
link11 =  link1+"Cluster"+www 
os.mkdir(link11)

# rds_prefix.pdf/Cluster
link21 =  link2+"Cluster"+www
os.mkdir(link21)

#存放Sample信息的总文件夹
#rds_prefix.png/Sample
link12 =  link1+"Sample"+www
os.mkdir(link12)

#rds_prefix.pdf/Sample
link22 =  link2+"Sample"+www
os.mkdir(link22)

####绘制全体Cluster的图
#rds_prefix.png/Cluster/AllCluster
link111 =  link11+"AllCluster"+www 
os.mkdir(link111)

#rds_prefix.pdf/Cluster/AllCluster
link211 =  link21+"AllCluster"+www 
os.mkdir(link211)

####绘制全体Sample的图

#rds_prefix.png/Sample/AllSample
link122 =  link12+"AllSample"+www
os.mkdir(link122)

#rds_prefix.pdf/Sample/AllSample
link222 =  link22+"AllSample"+www
os.mkdir(link222)

####绘制全体Cluster的QC图

#rds_prefix.png/Cluster/AllCluster/QC_Graph
link1111 =  link111+"QC_Graph"+www 
os.mkdir(link1111)

#rds_prefix.pdf/Cluster/AllCluster/QC_Graph
link2111 =  link211+"QC_Graph"+www 
os.mkdir(link2111)

####绘制全体Cluster的clonotype_cluster_network图 

#rds_prefix.png/Cluster/AllCluster/clonotype_cluster_network
link1112 =  link111+"clonotype_cluster_network"+www 
os.mkdir(link1112)

#rds_prefix.pdf/Cluster/AllCluster/clonotype_cluster_network
link2112 =  link211+"clonotype_cluster_network"+www 
os.mkdir(link2112)

####绘制全体Cluster的clonotype_expansion图

#rds_prefix.png/Cluster/AllCluster/clonotype_expansion
link1113 =  link111+"clonotype_expansion"+www 
os.mkdir(link1113)

#rds_prefix.pdf/Cluster/AllCluster/clonotype_expansion
link2113 =  link211+"clonotype_expansion"+www 
os.mkdir(link2113)

###全体Cluster的clonotype_expansion中的不同类图分别存储

#rds_prefix.png/Cluster/AllCluster/clonotype_expansion/umapinfo
link11131 =  link1113+"umapinfo"+www 
os.mkdir(link11131)

#rds_prefix.pdf/Cluster/AllCluster/clonotype_expansion/umapinfo
link21131 =  link2113+"umapinfo"+www 
os.mkdir(link21131)

#rds_prefix.png/Cluster/AllCluster/clonotype_expansion/clonal_expansion_stackedbarplot
link11132 =  link1113+"clonal_expansion_stackedbarplot"+www 
os.mkdir(link11132)

#rds_prefix.pdf/Cluster/AllCluster/clonotype_expansion/clonal_expansion_stackedbarplot
link21132 =  link2113+"clonal_expansion_stackedbarplot"+www 
os.mkdir(link21132)

#rds_prefix.png/Cluster/AllCluster/clonotype_expansion/clonotype_diversity
link11133 =  link1113+"clonotype_diversity"+www 
os.mkdir(link11133)

#rds_prefix.pdf/Cluster/AllCluster/clonotype_expansion/clonotype_diversity
link21133 =  link2113+"clonotype_diversity"+www 
os.mkdir(link21133)

#rds_prefix.png/Cluster/AllCluster/clonotype_expansion/clonetpye_size_abundance
link11134 =  link1113+"clonetpye_size_abundance"+www 
os.mkdir(link11134)

#rds_prefix.pdf/Cluster/AllCluster/clonotype_expansion/clonetpye_size_abundance
link21134 =  link2113+"clonetpye_size_abundance"+www 
os.mkdir(link21134)

####绘制全体Sample的QC图

#rds_prefix.png/Sample/AllSample/QC_Graph
link1221 =  link122+"QC_Graph"+www
os.mkdir(link1221)

#rds_prefix.pdf/Sample/AllSample/QC_Graph
link2221 =  link222+"QC_Graph"+www
os.mkdir(link2221)

####绘制全体Sample的clonotype_cluster_network

#rds_prefix.png/Sample/AllSample/clonotype_cluster_network
link1222 =  link122+"clonotype_cluster_network"+www
os.mkdir(link1222)

#rds_prefix.pdf/Sample/AllSample/clonotype_cluster_network
link2222 =  link222+"clonotype_cluster_network"+www
os.mkdir(link2222)

####绘制全体Cluster的clonotype_expansion图

#rds_prefix.png/Sample/AllSample/clonotype_expansion
link1223 =  link122+"clonotype_expansion"+www
os.mkdir(link1223)

#rds_prefix.pdf/Sample/AllSample/clonotype_expansion
link2223 =  link222+"clonotype_expansion"+www
os.mkdir(link2223)

###全体Sample的clonotype_expansion中的不同类图分别存储

#rds_prefix.png/Sample/AllSample/clonotype_expansion/umapinfo
link12231 =  link1223+"umapinfo"+www 
os.mkdir(link12231)

#rds_prefix.pdf/Sample/AllSample/clonotype_expansion/umapinfo
link22231 =  link2223+"umapinfo"+www 
os.mkdir(link22231)

#rds_prefix.png/Sample/AllSample/clonotype_expansion/clonal_expansion_stackedbarplot
link12232 =  link1223+"clonal_expansion_stackedbarplot"+www 
os.mkdir(link12232)

#rds_prefix.pdf/Sample/AllSample/clonotype_expansion/clonal_expansion_stackedbarplot
link22232 =  link2223+"clonal_expansion_stackedbarplot"+www 
os.mkdir(link22232)

#rds_prefix.png/Sample/AllSample/clonotype_expansion/clonotype_diversity
link12233 =  link1223+"clonotype_diversity"+www 
os.mkdir(link12233)

#rds_prefix.pdf/Sample/AllSample/clonotype_expansion/clonotype_diversity
link22233 =  link2223+"clonotype_diversity"+www 
os.mkdir(link22233)

#rds_prefix.png/Sample/AllSample/clonotype_expansion/clonetpye_size_abundance
link12234 =  link1223+"clonetpye_size_abundance"+www 
os.mkdir(link12234)

#rds_prefix.pdf/Sample/AllSample/clonotype_expansion/clonetpye_size_abundance
link22234 =  link2223+"clonetpye_size_abundance"+www 
os.mkdir(link22234)

#################################################################加载数据 

#20230111 这里合并ir和单细胞数据后 要把ir里面的uns转存到新的anndata中 然后测试一下重算布局最大的clonotype有没有变化 
# 这样可以知道有没有按细胞名检索uns中的clonotype_id 如果这里的uns成功继承 那么全套task可以沿用最初norm里的id

adatair = pickle.load(open(adatair,"rb"))

#读取单细胞数据
adata = sc.read(seuseth5)
#以单细胞数据的细胞和基因为准 将ir信息加载到adata中
ir.pp.merge_with_ir(adata,adatair)
adata.uns = adatair.uns

#这里合并虽然以单细胞为准 但是并不是所有细胞都有ir信息 
# 20230207 一个比较严重的问题 单细胞数据集和ir数据里的细胞真的全部有重合吗？
# 如果不是mutliaggr下来的可能会不一样 
# 考虑在这添加一步取交集 把adata里的细胞全部替换成出现在ir数据里的 
aaa = set(list(adata.obs.index))
bbb = set(list(adatair.obs.index))
commons = list(aaa.intersection(bbb))
adata = adata[commons,:]

print("commoncellsnum of GEX and IR :" + str(len(commons)))

# 20240708 在这里统计overlap的细胞数占TCR的总细胞数的占比
returnvalue = 0
if len(commons)/len(list(bbb)) < 0.3:
  returnvalue = 223
if len(commons)/len(list(bbb)) < 0.1:
  returnvalue = 233

#函数返回adata 
# 有可能出现 ir里的样本名和rds里的样本名不同 那么以rds里的样本名为准
# 在rds里把orig.ident改为sample

#20230310 前面技术部测试发现 size是鉴定clonetype这一步函数算好的 从norm传到cluster继承的还是norm中的数目 
# 但是和单细胞合并后会少很多细胞 由于我们不重新鉴定clonetype这一步就和真实结果有出入
# 所以这里要扒一下代码重新算一遍
#这个函数比较有用 这个步骤有点类似于table一下再匹配着填回去  有机会要研究一下
clonotype_cluster_series = adata.obs['cc_aa_alignment']
clonotype_cluster_size_series = clonotype_cluster_series.groupby(clonotype_cluster_series).transform("count")
adata.obs['cc_aa_alignment_size'] = clonotype_cluster_size_series.astype("float64")

###################################################按照输入文件修改cluster和sample顺序
#20230426需求 要传入一列cluster一列color的表 调整横坐标顺序和图例颜色
#20230417需求 画图需要自定义cluster/sample顺序 输入两个新文件 要发布新版本

if clusterorder is not None:
  # 读取改名文件 默认为单列无行名 
  # 20230418修改 平台生成文件默认带表头
  df11 = pd.read_csv(clusterorder,sep="\t",header=0)
  haha = df11.iloc[:,0].tolist()
  # 防止出现数字cluster这种 
  haha = [str(num) for num in haha]
  #修改level顺序的方式
  adata.obs['seurat_clusters'] = adata.obs['seurat_clusters'].astype(str)
  cat_type = CategoricalDtype(categories=haha, ordered=False)
  adata.obs['seurat_clusters'] = adata.obs['seurat_clusters'].astype(cat_type)
  if df11.shape[1] > 1:
    # 如果包含第二列 那么就是颜色 
    coco = df11.iloc[:,1].tolist()
    adata.uns['seurat_clusters_colors'] = coco

########################改sample顺序

if sampleorder is not None:
  # 读取改名文件 默认为单列无行名 
  # 20230418修改 平台生成文件默认带表头
  df1 = pd.read_csv(sampleorder,sep="\t",header=0)
  haha1 = df1.iloc[:,0].tolist()
  # 防止出现数字cluster这种 
  haha1 = [str(num) for num in haha1]
  #修改level顺序的方式
  adata.obs['sample'] = adata.obs['sample'].astype(str)
  cat_type = CategoricalDtype(categories=haha1, ordered=False)
  adata.obs['sample'] = adata.obs['sample'].astype(cat_type)
  if df1.shape[1] > 1:
    # 如果包含第二列 那么就是颜色 
    coco = df1.iloc[:,1].tolist()
    adata.uns['sample_colors'] = coco



##################################################################重新按照cluster绘制qc分布  分样本也要画

######################################为了方便之后绘图 这里会给整个数据按照样本/celltype拆分 

#20221230发现 这里每个subadata 在后续的计算中都是独立的 所以最后需要把每个sub的obs矩阵全部写出来
Sample = adata.obs['sample'].unique()
Cluster = adata.obs['seurat_clusters'].unique()
#这里新增一个规则 如果只有一种样本 就不画分样本的图
print("Cluster num :"+str(len(Cluster)))
print("Sample num :"+str(len(Sample)))

a = {}

if len(Cluster) > 1:
  for i in Cluster:
    subadata = adata[adata.obs['seurat_clusters'] == i]
    a[i] = subadata
    ####绘制Cluster_i的图
    #rds_prefix.png/Cluster/Cluster_i
    link11_i =  link11+"Cluster_"+i+www 
    os.mkdir(link11_i)
    
    #rds_prefix.pdf/Cluster/Cluster_i
    link21_i =  link21+"Cluster_"+i+www 
    os.mkdir(link21_i)
    
    ####绘制Cluster_i的QC图
    
    #rds_prefix.png/Cluster/Cluster_i/QC_Graph
    link11_i1 =  link11_i+"QC_Graph"+www 
    os.mkdir(link11_i1)
    
    #rds_prefix.pdf/Cluster/Cluster_i/QC_Graph
    link21_i1 =  link21_i+"QC_Graph"+www 
    os.mkdir(link21_i1)
    
    ####绘制Cluster_i的clonotype_cluster_network图 
    
    #rds_prefix.png/Cluster/Cluster_i/clonotype_cluster_network
    link11_i2 =  link11_i+"clonotype_cluster_network"+www 
    os.mkdir(link11_i2)
    
    #rds_prefix.pdf/Cluster/Cluster_i/clonotype_cluster_network
    link21_i2 =  link21_i+"clonotype_cluster_network"+www 
    os.mkdir(link21_i2)
    
    ####绘制Cluster_i的clonotype_expansion图
    
    #rds_prefix.png/Cluster/Cluster_i/clonotype_expansion
    link11_i3 =  link11_i+"clonotype_expansion"+www 
    os.mkdir(link11_i3)
    
    #rds_prefix.pdf/Cluster/Cluster_i/clonotype_expansion
    link21_i3 =  link21_i+"clonotype_expansion"+www 
    os.mkdir(link21_i3)
    
    ###Cluster_i的clonotype_expansion中的不同类图分别存储
    
    #rds_prefix.png/Cluster/Cluster_i/clonotype_expansion/umapinfo
    link11_i32 =  link11_i3+"umapinfo"+www 
    os.mkdir(link11_i32)
    
    #rds_prefix.pdf/Cluster/Cluster_i/clonotype_expansion/umapinfo
    link21_i32 =  link21_i3+"umapinfo"+www 
    os.mkdir(link21_i32)
    
    #rds_prefix.png/Cluster/Cluster_i/clonotype_expansion/clonal_expansion_stackedbarplot
    link11_i32 =  link11_i3+"clonal_expansion_stackedbarplot"+www 
    os.mkdir(link11_i32)
    
    #rds_prefix.pdf/Cluster/Cluster_i/clonotype_expansion/clonal_expansion_stackedbarplot
    link21_i32 =  link21_i3+"clonal_expansion_stackedbarplot"+www 
    os.mkdir(link21_i32)
    
    #rds_prefix.png/Cluster/Cluster_i/clonotype_expansion/clonetpye_size_abundance
    link11_i34 =  link11_i3+"clonetpye_size_abundance"+www 
    os.mkdir(link11_i34)
    
    #rds_prefix.pdf/Cluster/Cluster_i/clonotype_expansion/clonetpye_size_abundance
    link21_i34 =  link21_i3+"clonetpye_size_abundance"+www 
    os.mkdir(link21_i34)


b = {}
if len(Sample) > 1:
  for j in Sample:
    subadata = adata[adata.obs['sample'] == j]
    b[j] = subadata
    ####绘制Sample_j的图
    #rds_prefix.png/Sample/Sample_j
    link12_j =  link12+"Sample_"+j+www 
    os.mkdir(link12_j)
    
    #rds_prefix.pdf/Sample/Sample_j
    link22_j =  link22+"Sample_"+j+www 
    os.mkdir(link22_j)
    
    ####绘制Sample_i的QC图
    
    #rds_prefix.png/Sample/Sample_j/QC_Graph
    link12_j1 =  link12_j+"QC_Graph"+www 
    os.mkdir(link12_j1)
    
    #rds_prefix.pdf/Sample/Sample_j/QC_Graph
    link22_j1 =  link22_j+"QC_Graph"+www 
    os.mkdir(link22_j1)
    
    ####绘制Sample_j的clonotype_cluster_network图 
    
    #rds_prefix.png/Sample/Sample_j/clonotype_cluster_network
    link12_j2 =  link12_j+"clonotype_cluster_network"+www 
    os.mkdir(link12_j2)
    
    #rds_prefix.pdf/Sample/Sample_j/clonotype_cluster_network
    link22_j2 =  link22_j+"clonotype_cluster_network"+www 
    os.mkdir(link22_j2)
    
    ####绘制Sample_j的clonotype_expansion图
    
    #rds_prefix.png/Sample/Sample_j/clonotype_expansion
    link12_j3 =  link12_j+"clonotype_expansion"+www 
    os.mkdir(link12_j3)
    
    #rds_prefix.pdf/Sample/Sample_j/clonotype_expansion
    link22_j3 =  link22_j+"clonotype_expansion"+www 
    os.mkdir(link22_j3)
    
    ###Sample_j的clonotype_expansion中的不同类图分别存储
    
    #rds_prefix.png/Sample/Sample_j/clonotype_expansion/umapinfo
    link12_j31 =  link12_j3+"umapinfo"+www 
    os.mkdir(link12_j31)
    
    #rds_prefix.pdf/Sample/Sample_j/clonotype_expansion/umapinfo
    link22_j31 =  link22_j3+"umapinfo"+www 
    os.mkdir(link22_j31)
    
    #rds_prefix.png/Sample/Sample_j/clonotype_expansion/clonal_expansion_stackedbarplot
    link12_j32 =  link12_j3+"clonal_expansion_stackedbarplot"+www 
    os.mkdir(link12_j32)
    
    #rds_prefix.pdf/Sample/Sample_j/clonotype_expansion/clonal_expansion_stackedbarplot
    link22_j32 =  link22_j3+"clonal_expansion_stackedbarplot"+www 
    os.mkdir(link22_j32)
    
    #rds_prefix.png/Sample/Sample_j/clonotype_expansion/clonetpye_size_abundance
    link12_j34 =  link12_j3+"clonetpye_size_abundance"+www 
    os.mkdir(link12_j34)
    
    #rds_prefix.pdf/Sample/Sample_j/clonotype_expansion/clonetpye_size_abundance
    link22_j34 =  link22_j3+"clonetpye_size_abundance"+www 
    os.mkdir(link22_j34)

if len(Cluster) > 1:
  for i in Cluster:
    print("Cluster:_"+i+"_cellnum_"+str(len(a[i].obs.index)))

if len(Sample) > 1:
  for j in Sample:
    print("Sample_"+j+"_cellnum_"+str(len(b[j].obs.index)))


#######################################################绘制QC图

#rds_prefix.png/Cluster/AllCluster/QC_Graph
#link1111 
#rds_prefix.pdf/Cluster/AllCluster/QC_Graph
#link2111 

#rds_prefix.png/Cluster/Cluster_i/QC_Graph
#link11_i1
#rds_prefix.pdf/Cluster/Cluster_i/QC_Graph
#link21_i1

ax = ir.pl.group_abundance(adata, groupby="receptor_subtype", target_col="seurat_clusters",fig_kws = {'figsize': (15, 12), 'dpi': 300})

handles, labels = plt.gca().get_legend_handles_labels()
list1 = list(adata.obs['seurat_clusters'].cat.categories)
indices = [labels.index(x) for x in list1]
ncol1=(1 if len(list1) <= 21 else 2 if len(list1) <= 30 else 3)
place = 1.5
if ncol1 == 2:
  place = 1.8

if ncol1 == 3:
  place = 2.0

plt.legend([handles[i] for i in indices], [labels[i] for i in indices],loc = 'center right',frameon = False,title = "celltype",bbox_to_anchor  =(place,0.5),ncol=ncol1)

aa =  link1111 + "ALLCluster_qc_receptor_subtype.png"
aa1 = link2111 + "ALLCluster_qc_receptor_subtype.pdf"   
plt.savefig(aa,bbox_inches='tight')
plt.savefig(aa1,bbox_inches='tight')

ax = ir.pl.group_abundance(adata, groupby="receptor_type", target_col="seurat_clusters",fig_kws = {'figsize': (15, 12), 'dpi': 300})

handles, labels = plt.gca().get_legend_handles_labels()
list1 = list(adata.obs['seurat_clusters'].cat.categories)
indices = [labels.index(x) for x in list1]
ncol1=(1 if len(list1) <= 21 else 2 if len(list1) <= 30 else 3)
place = 1.5
if ncol1 == 2:
  place = 1.8

if ncol1 == 3:
  place = 2.0

plt.legend([handles[i] for i in indices], [labels[i] for i in indices],loc = 'center right',frameon = False,title = "celltype",bbox_to_anchor  =(place,0.5),ncol=ncol1)

aa =  link1111 + "ALLCluster_qc_receptor_type.png"
aa1 = link2111 + "ALLCluster_qc_receptor_type.pdf"   
plt.savefig(aa,bbox_inches='tight')
plt.savefig(aa1,bbox_inches='tight')

ax = ir.pl.group_abundance(adata, groupby="chain_pairing",target_col="seurat_clusters",fig_kws = {'figsize': (15,12), 'dpi': 300})

handles, labels = plt.gca().get_legend_handles_labels()
list1 = list(adata.obs['seurat_clusters'].cat.categories)
indices = [labels.index(x) for x in list1]
ncol1=(1 if len(list1) <= 21 else 2 if len(list1) <= 30 else 3)
place = 1.5
if ncol1 == 2:
  place = 1.8

if ncol1 == 3:
  place = 2.0

plt.legend([handles[i] for i in indices], [labels[i] for i in indices],loc = 'center right',frameon = False,title = "celltype",bbox_to_anchor  =(place,0.5),ncol=ncol1)

aa =  link1111 + "ALLCluster_qc_chain_pairing.png"
aa1 = link2111 + "ALLCluster_qc_chain_pairing.pdf"   
plt.savefig(aa,bbox_inches='tight')
plt.savefig(aa1,bbox_inches='tight')

if len(Cluster) > 1:
  for i in Cluster:
    folder =  link11+"Cluster_"+i+"/QC_Graph/"
    folder1 = link21+"Cluster_"+i+"/QC_Graph/"
    ax = ir.pl.group_abundance(a[i], groupby="receptor_subtype", target_col="sample",fig_kws = {'figsize': (15,12), 'dpi': 300})
    
    handles, labels = plt.gca().get_legend_handles_labels()
    list1 = list(a[i].obs['sample'].cat.categories)
    indices = [labels.index(x) for x in list1]
    ncol1=(1 if len(list1) <= 21 else 2 if len(list1) <= 30 else 3)
    place = 1.5
    if ncol1 == 2:
      place = 1.8
    
    if ncol1 == 3:
      place = 2.0
    
    plt.legend([handles[i] for i in indices], [labels[i] for i in indices],loc = 'center right',frameon = False,title = "sample",bbox_to_anchor  =(place,0.5),ncol=ncol1)
    
    aa =  folder  + "Cluster_" + i + "_qc_receptor_subtype_groupby_sample.png"
    aa1 = folder1 + "Cluster_" + i + "_qc_receptor_subtype_groupby_sample.pdf"   
    plt.savefig(aa,bbox_inches='tight')
    plt.savefig(aa1,bbox_inches='tight')
    
    ax = ir.pl.group_abundance(a[i], groupby="receptor_type", target_col="sample",fig_kws = {'figsize': (15,12), 'dpi': 300})
    
    handles, labels = plt.gca().get_legend_handles_labels()
    list1 = list(a[i].obs['sample'].cat.categories)
    indices = [labels.index(x) for x in list1]
    ncol1=(1 if len(list1) <= 21 else 2 if len(list1) <= 30 else 3)
    place = 1.5
    if ncol1 == 2:
      place = 1.8
    
    if ncol1 == 3:
      place = 2.0
    
    plt.legend([handles[i] for i in indices], [labels[i] for i in indices],loc = 'center right',frameon = False,title = "sample",bbox_to_anchor  =(place,0.5),ncol=ncol1)
    
    aa =  folder  + "Cluster_" + i + "_qc_receptor_type_groupby_sample.png"
    aa1 = folder1 + "Cluster_" + i + "_qc_receptor_type_groupby_sample.pdf"   
    plt.savefig(aa,bbox_inches='tight')
    plt.savefig(aa1,bbox_inches='tight')
    
    ax = ir.pl.group_abundance(a[i], groupby="chain_pairing",target_col="sample",fig_kws = {'figsize': (15,12), 'dpi': 300})
    
    handles, labels = plt.gca().get_legend_handles_labels()
    list1 = list(a[i].obs['sample'].cat.categories)
    indices = [labels.index(x) for x in list1]
    ncol1=(1 if len(list1) <= 21 else 2 if len(list1) <= 30 else 3)
    place = 1.5
    if ncol1 == 2:
      place = 1.8
    
    if ncol1 == 3:
      place = 2.0
    
    plt.legend([handles[i] for i in indices], [labels[i] for i in indices],loc = 'center right',frameon = False,title = "sample",bbox_to_anchor  =(place,0.5),ncol=ncol1)
    
    aa =  folder  + "Cluster_" + i + "_qc_chain_pairing_groupby_sample.png"
    aa1 = folder1 + "Cluster_" + i + "_qc_chain_pairing_groupby_sample.pdf"   
    plt.savefig(aa,bbox_inches='tight')
    plt.savefig(aa1,bbox_inches='tight')

########################################绘制样本图

#rds_prefix.png/Sample/AllSample/QC_Graph
#link1221
#rds_prefix.pdf/Sample/AllSample/QC_Graph
#link2221

ax = ir.pl.group_abundance(adata, groupby="receptor_subtype", target_col="sample",fig_kws = {'figsize': (15,12), 'dpi': 300})

handles, labels = plt.gca().get_legend_handles_labels()
# 这里的handles和labels参数分别表示图例的外观和标签 
# 我们这里使用的方式是先让plt创建好图例 再调整图例的顺序 
# list2.sort(key=lambda x: list1.index(x))
list1 = list(adata.obs['sample'].cat.categories)
indices = [labels.index(x) for x in list1]

ncol1=(1 if len(list1) <= 21 else 2 if len(list1) <= 30 else 3)
place = 1.5
if ncol1 == 2:
  place = 1.8

if ncol1 == 3:
  place = 2.0

plt.legend([handles[i] for i in indices], [labels[i] for i in indices],loc = 'center right',frameon = False,title = "sample",bbox_to_anchor  =(place,0.5),ncol=ncol1)

aa =  link1221 + "ALLSample_qc_receptor_subtype.png"
aa1 = link2221 + "ALLSample_qc_receptor_subtype.pdf"   
plt.savefig(aa,bbox_inches='tight')
plt.savefig(aa1,bbox_inches='tight')
  
ax = ir.pl.group_abundance(adata, groupby="receptor_type", target_col="sample",fig_kws = {'figsize': (15,12), 'dpi': 300})

handles, labels = plt.gca().get_legend_handles_labels()
list1 = list(adata.obs['sample'].cat.categories)
indices = [labels.index(x) for x in list1]

ncol1=(1 if len(list1) <= 21 else 2 if len(list1) <= 30 else 3)
place = 1.5
if ncol1 == 2:
  place = 1.8

if ncol1 == 3:
  place = 2.0

plt.legend([handles[i] for i in indices], [labels[i] for i in indices],loc = 'center right',frameon = False,title = "sample",bbox_to_anchor  =(place,0.5),ncol=ncol1)

aa = link1221 + "ALLSample_qc_receptor_type.png"
aa1 = link2221 + "ALLSample_qc_receptor_type.pdf"   
plt.savefig(aa,bbox_inches='tight')
plt.savefig(aa1,bbox_inches='tight')
  
ax = ir.pl.group_abundance(adata, groupby="chain_pairing",target_col="sample",fig_kws = {'figsize': (15,12), 'dpi': 300})

handles, labels = plt.gca().get_legend_handles_labels()
list1 = list(adata.obs['sample'].cat.categories)
indices = [labels.index(x) for x in list1]

ncol1=(1 if len(list1) <= 21 else 2 if len(list1) <= 30 else 3)
place = 1.5
if ncol1 == 2:
  place = 1.8

if ncol1 == 3:
  place = 2.0

plt.legend([handles[i] for i in indices], [labels[i] for i in indices],loc = 'center right',frameon = False,title = "sample",bbox_to_anchor  =(place,0.5),ncol=ncol1)

aa =  link1221 + "ALLSample_qc_chain_pairing.png"
aa1 = link2221 + "ALLSample_qc_chain_pairing.pdf"   
plt.savefig(aa,bbox_inches='tight')
plt.savefig(aa1,bbox_inches='tight')
  
if len(Sample) > 1:
  for j in Sample:
    folder =  link12+"Sample_"+j+"/QC_Graph/"
    folder1 = link22+"Sample_"+j+"/QC_Graph/"
    ax = ir.pl.group_abundance(b[j], groupby="receptor_subtype", target_col="seurat_clusters",fig_kws = {'figsize': (15,12), 'dpi': 300})
    
    handles, labels = plt.gca().get_legend_handles_labels()
    list1 = list(b[j].obs['seurat_clusters'].cat.categories)
    indices = [labels.index(x) for x in list1]
    
    ncol1=(1 if len(list1) <= 21 else 2 if len(list1) <= 30 else 3)
    place = 1.5
    if ncol1 == 2:
      place = 1.8
    
    if ncol1 == 3:
      place = 2.0
    
    plt.legend([handles[i] for i in indices], [labels[i] for i in indices],loc = 'center right',frameon = False,title = "celltype",bbox_to_anchor  =(place,0.5),ncol=ncol1)
    
    aa =   folder + "Sample_" + j + "_qc_receptor_subtype_groupby_celltype.png"
    aa1 = folder1 + "Sample_" + j + "_qc_receptor_subtype_groupby_celltype.pdf"   
    plt.savefig(aa,bbox_inches='tight')
    plt.savefig(aa1,bbox_inches='tight')
    
    ax = ir.pl.group_abundance(b[j], groupby="receptor_type", target_col="seurat_clusters",fig_kws = {'figsize': (15,12), 'dpi': 300})
    
    handles, labels = plt.gca().get_legend_handles_labels()
    list1 = list(b[j].obs['seurat_clusters'].cat.categories)
    indices = [labels.index(x) for x in list1]
    
    ncol1=(1 if len(list1) <= 21 else 2 if len(list1) <= 30 else 3)
    place = 1.5
    if ncol1 == 2:
      place = 1.8
    
    if ncol1 == 3:
      place = 2.0
    
    plt.legend([handles[i] for i in indices], [labels[i] for i in indices],loc = 'center right',frameon = False,title = "celltype",bbox_to_anchor  =(place,0.5),ncol=ncol1)
    
    aa =   folder + "Sample_" + j + "_receptor_type_groupby_celltype.png"
    aa1 = folder1 + "Sample_" + j + "_receptor_type_groupby_celltype.pdf"   
    plt.savefig(aa,bbox_inches='tight')
    plt.savefig(aa1,bbox_inches='tight')
    
    ax = ir.pl.group_abundance(b[j], groupby="chain_pairing",target_col="seurat_clusters",fig_kws = {'figsize': (15,12), 'dpi': 300})
    
    handles, labels = plt.gca().get_legend_handles_labels()
    list1 = list(b[j].obs['seurat_clusters'].cat.categories)
    indices = [labels.index(x) for x in list1]
    
    ncol1=(1 if len(list1) <= 21 else 2 if len(list1) <= 30 else 3)
    place = 1.5
    if ncol1 == 2:
      place = 1.8
    
    if ncol1 == 3:
      place = 2.0
    
    plt.legend([handles[i] for i in indices], [labels[i] for i in indices],loc = 'center right',frameon = False,title = "celltype",bbox_to_anchor  =(place,0.5),ncol=ncol1)
    
    aa =   folder + "Sample_" + j + "_qc_chain_pairing_groupby_celltype.png"
    aa1 = folder1 + "Sample_" + j + "_qc_chain_pairing_groupby_celltype.pdf"   
    plt.savefig(aa,bbox_inches='tight')
    plt.savefig(aa1,bbox_inches='tight')


################################################################绘制clonotype_cluster布局图 

# min_cells2 = 12
# min_nodes2 = 3
# base_size2 = None
# size_power2 = 0.6
# label_fontsize2 = 9
# label_alpha2 = 0.6
# edges_width2 = 0.5

###############################################################clonotype_cluster的布局图

####################################首先计算allcluster/allsample

######################################################建立clonotype_cluster关系图

#20230112 检查发现 确实能够在图中继承clonotype_id 但是uns的结果输出dict不止有一个cloneid对应细胞名的矩阵 还有clonotype之间的距离矩阵
# 这个距离矩阵是按照删减之前的细胞数计算的 所以继承之前的uns虽然画图没问题 但是点的位置布局和之前完全相同 比例也是固定的
# 缺失细胞导致消失的clonotype会在图中原位置扣掉 如果细胞删减的太少 图会很空旷 但是好处是以后所有涉及到clonotype的图 完全能够统一起来

#rds_prefix.png/Cluster/AllCluster/clonotype_cluster_network
#link1112
#rds_prefix.pdf/Cluster/AllCluster/clonotype_cluster_network
#link2112
#rds_prefix.png/Sample/AllSample/clonotype_cluster_network
#link1222
#rds_prefix.pdf/Sample/AllSample/clonotype_cluster_network
#link2222

#明确一点 调用鉴定cc生成的 cc_aa_alignment pl固定调用uns中的clonotype_network

ir.tl.clonotype_network(
  adata, 
  min_cells=min_cells2, 
  sequence='aa', 
  metric='alignment',
  min_nodes = min_nodes2,
  size_power = size_power2,
  base_size = base_size2,
  clonotype_key  = "cc_aa_alignment",
  inplace = True
  )

#########################################################绘制clonotype_cluster布局图

ir.pl.clonotype_network(
  adata, 
  color="seurat_clusters",
  label_fontsize=label_fontsize2,
  panel_size=(12, 10),
  label_alpha = label_alpha2,
  base_size = base_size2,
  edges_width = edges_width2
)

aa =  link1112 + "AllCluster_clonotype_cluster_network.png"
aa1 = link2112 + "AllCluster_clonotype_cluster_network.pdf"   
plt.savefig(aa,bbox_inches='tight')
plt.savefig(aa1,bbox_inches='tight')

ir.pl.clonotype_network(adata,color="sample",label_fontsize=label_fontsize2, panel_size=(12, 10),label_alpha = label_alpha2,base_size = base_size2,edges_width = edges_width2)
  
aa =  link1222 + "AllSample_pl_clonotype_cluster_network.png"
aa1 = link2222 + "AllSample_pl_clonotype_cluster_network.pdf"   
plt.savefig(aa,bbox_inches='tight')
plt.savefig(aa1,bbox_inches='tight')


###################################分样本/cluster绘图
#20230103 出现bug 有些细胞数太少的在mincell参数设为2的时候tlnetwork报错 即没有两个及以上的细胞表达同一种cd3序列
#在这写了一个trycatch 注意python的trycatch只能判断一条命令 而且要事先清楚错误类型 比如我这里就是ValueError 
#logging模块的作用是记录报错信息但是不停止代码运行
#20220114 这里分样本绘图就染上cluster的颜色 分cluster同理
if len(Cluster) > 1:
  for i in Cluster:
    folder =  link11+"Cluster_"+i+"/clonotype_cluster_network/"
    folder1 = link21+"Cluster_"+i+"/clonotype_cluster_network/"
    try:
      print("try ir.tl.clonotype_cluster_network_Cluster_" + i)
      ir.tl.clonotype_network(
        a[i],
        min_cells=min_cells2, 
        sequence='aa', 
        metric='alignment',
        min_nodes = min_nodes2,
        size_power = size_power2,
        base_size = base_size2,
        clonotype_key  = "cc_aa_alignment",
        inplace = True
      )
      ir.pl.clonotype_network(
        a[i], 
        color="sample",
        label_fontsize=label_fontsize2,
        panel_size=(12, 10),
        label_alpha = label_alpha2,
        base_size = base_size2,
        edges_width = edges_width2
      )
      aa =   folder + "Cluster_" + i + "_clonotype_cluster_network.png"
      aa1 = folder1 + "Cluster_" + i + "_clonotype_cluster_network.pdf"   
      plt.savefig(aa,bbox_inches='tight')
      plt.savefig(aa1,bbox_inches='tight')
    except ValueError as e:
      logging.exception(e)
      print('except:', e)
      print("except ir.tl.clonotype_cluster_network_Cluster_" + i)
      continue


# 20230115新报错 再测试过程中我发现有一些图不在tl报错 但是在pl报错 try except需要把tl和pl全部卸写进来
# ValueError: zero-size array to reduction operation 报错信息含义应该是 tl之后pl没办法通过布局矩阵计算出细胞坐标 
# 之所以逻辑产生变化 因为第一版会给每个subadata重算一次clonotype 这样tl就会直接算不了 
# V2版本不重算 距离矩阵是直接继承的总表 tl用总表运算 就是在原位保留细胞位置  但是pl时候会因为分的太开或者细胞太少了画不了图 细胞坐标也是空的 

if len(Sample) > 1:
  for j in Sample:
    folder =  link12+"Sample_"+j+"/clonotype_cluster_network/"
    folder1 = link22+"Sample_"+j+"/clonotype_cluster_network/"
    try:
      print("try ir.tl.clonotype_cluster_network_Sample_" + j)
      ir.tl.clonotype_network(
        b[j],
        min_cells=min_cells2, 
        sequence='aa', 
        metric='alignment',
        min_nodes = min_nodes2,
        size_power = size_power2,
        base_size = base_size2,
        clonotype_key  = "cc_aa_alignment",
        inplace = True
      )
      ir.pl.clonotype_network(
        b[j], 
        color="seurat_clusters",
        label_fontsize=label_fontsize2,
        panel_size=(12, 10),
        label_alpha = label_alpha2,
        base_size = base_size2,
        edges_width = edges_width2
      ) 
      aa =   folder + "Sample_" + j + "_clonotype_cluster_network.png"
      aa1 = folder1 + "Sample_" + j + "_clonotype_cluster_network.pdf"   
      plt.savefig(aa,bbox_inches='tight')
      plt.savefig(aa1,bbox_inches='tight')
    except ValueError as e:
      logging.exception(e)
      print('except:', e)
      print("except ir.tl.clonotype_network_Sample_" + j)
      continue

###############################################################对于cc的独特功能：
# 限定V基因相同后 cc中有些clonotype_cluster会被分成几个新的cc

# 这里再加一个需求 如果总细胞大于一个阈值 那么就不加这个功能了 
# 下个版本再加一个判断参数 如果确定不做直接pass
dosamev  = True
maxcell = 100000

if dosamev:
  if len(adata.obs.index) < maxcell:
    if "cc_aa_alignment_same_v" in adata.obs.columns:
      ct_different_v = adata.obs.groupby("cc_aa_alignment").apply(
          lambda x: x["cc_aa_alignment_same_v"].nunique() > 1
      )
      ct_different_v = ct_different_v[ct_different_v].index.values.tolist()
      ct_different_v
      
      ct_different = adata.obs.loc[
          adata.obs["cc_aa_alignment"].isin(ct_different_v),
          [
              "cc_aa_alignment",
              "cc_aa_alignment_same_v",
              "IR_VJ_1_v_call",
              "IR_VDJ_1_v_call",
          ],
      ].sort_values("cc_aa_alignment").drop_duplicates().reset_index(drop=True)
      
      #此处这个表分别写到allcluster和allsample中
      
      aa = link1112 + "clonotypes_cluster_difference.txt"
      ct_different.to_csv(aa,sep='\t',index = False)
      
      aa = link2112 + "clonotypes_cluster_difference.txt"
      ct_different.to_csv(aa,sep='\t',index = False)
      
      aa = link1222 + "clonotypes_cluster_difference.txt"
      ct_different.to_csv(aa,sep='\t',index = False)
      
      aa = link2222 + "clonotypes_cluster_difference.txt"
      ct_different.to_csv(aa,sep='\t',index = False)

#对于分样本/分cluster来说 
# 一般来说这里不会超限 那么就判定一下dosamev
if dosamev:
  if len(Sample) > 1:
    for j in Sample:
      if "cc_aa_alignment_same_v" in b[j].obs.columns:
        folder =  link12+"Sample_"+j+"/clonotype_cluster_network/"
        folder1 = link22+"Sample_"+j+"/clonotype_cluster_network/"
        ct_different_v = b[j].obs.groupby("cc_aa_alignment").apply(
        lambda x: x["cc_aa_alignment_same_v"].nunique() > 1
        )
        ct_different_v = ct_different_v[ct_different_v].index.values.tolist()   
        ct_different = b[j].obs.loc[
        b[j].obs["cc_aa_alignment"].isin(ct_different_v),
        [
            "cc_aa_alignment",
            "cc_aa_alignment_same_v",
            "IR_VJ_1_v_call",
            "IR_VDJ_1_v_call",
        ],
        ].sort_values("cc_aa_alignment").drop_duplicates().reset_index(drop=True)
        aa = folder + "Sample_" + j + "clonotypes_cluster_difference.txt"
        aa1 = folder1 + "Sample_" + j + "clonotypes_cluster_difference.txt"
        ct_different.to_csv(aa,sep='\t',index = False)
        ct_different.to_csv(aa1,sep='\t',index = False)
  
  if len(Cluster) > 1:
    for i in Cluster:
      if "cc_aa_alignment_same_v" in a[i].obs.columns:
        folder =  link11+"Cluster_"+i+"/clonotype_cluster_network/"
        folder1 = link21+"Cluster_"+i+"/clonotype_cluster_network/"
        ct_different_v = a[i].obs.groupby("cc_aa_alignment").apply(
        lambda x: x["cc_aa_alignment_same_v"].nunique() > 1
        )
        ct_different_v = ct_different_v[ct_different_v].index.values.tolist()   
        ct_different = a[i].obs.loc[
        a[i].obs["cc_aa_alignment"].isin(ct_different_v),
        [
            "cc_aa_alignment",
            "cc_aa_alignment_same_v",
            "IR_VJ_1_v_call",
            "IR_VDJ_1_v_call",
        ],
        ].sort_values("cc_aa_alignment").drop_duplicates().reset_index(drop=True)
        aa  =  folder + "Cluster_" + i + "clonotypes_cluster_difference.txt"
        aa1 = folder1 + "Cluster_" + i + "clonotypes_cluster_difference.txt"
        ct_different.to_csv(aa,sep='\t',index = False)
        ct_different.to_csv(aa1,sep='\t',index = False)


##########################################################################umap图展示基本信息
# clonal_expansion 指拓展类别 一般指单细胞克隆型 双细胞克隆型 大于两个的克隆型
# clonotype_size 表示clone中包含的细胞数 

ir.tl.clonal_expansion(adata,target_col = "cc_aa_alignment")

sc.set_figure_params(figsize = (15,15))
#clonal_expansion 名字是固定输出 
sc.pl.umap(adata, color="clonal_expansion",size=60)
aa = link11131 + "UMAP_clonal_expansion.png"
aa1 = link21131 + "UMAP_clonal_expansion.pdf"
plt.savefig(aa,bbox_inches='tight')
plt.savefig(aa1,bbox_inches='tight')

sc.set_figure_params(figsize = (15,15))
sc.pl.umap(adata, color="cc_aa_alignment_size",size=60)
aa = link11131 + "UMAP_clone_id_size.png"
aa1 = link21131 + "UMAP_clone_id_size.pdf"
plt.savefig(aa,bbox_inches='tight')
plt.savefig(aa1,bbox_inches='tight')

sc.set_figure_params(figsize = (15,15))
sc.pl.umap(adata, color="seurat_clusters",size=60)
aa = link11131 + "UMAP_seurat_clusters.png"
aa1 = link21131 + "UMAP_seurat_clusters.pdf"
plt.savefig(aa,bbox_inches='tight')
plt.savefig(aa1,bbox_inches='tight')

sc.set_figure_params(figsize = (15,15))
sc.pl.umap(adata, color="sample",size=60)
aa = link11131 + "UMAP_sample.png"
aa1 = link21131 + "UMAP_sample.pdf"
plt.savefig(aa,bbox_inches='tight')
plt.savefig(aa1,bbox_inches='tight')

#20230113 发现一个bug 这个umap四张图是不分样本和cluster的 结果我这里是做了区分的 
# 那么在allsample和allcluster中也要画这四张图 方便提取结果 
# 设计是 按样本拆分的subadata 就只画sample染色 按sample拆分的同理 

#rds_prefix.png/Sample/AllSample/clonotype_expansion/umapinfo
#link12231
#rds_prefix.pdf/Sample/AllSample/clonotype_expansion/umapinfo
#link22231

sc.set_figure_params(figsize = (15,15))
sc.pl.umap(adata, color="clonal_expansion",size=60)
aa = link12231 + "UMAP_clonal_expansion.png"
aa1 = link22231 + "UMAP_clonal_expansion.pdf"
plt.savefig(aa,bbox_inches='tight')
plt.savefig(aa1,bbox_inches='tight')

sc.set_figure_params(figsize = (15,15))
sc.pl.umap(adata, color="cc_aa_alignment_size",size=60)
aa = link12231 + "UMAP_clone_id_size.png"
aa1 = link22231 + "UMAP_clone_id_size.pdf"
plt.savefig(aa,bbox_inches='tight')
plt.savefig(aa1,bbox_inches='tight')

sc.set_figure_params(figsize = (15,15))
sc.pl.umap(adata, color="seurat_clusters",size=60)
aa = link12231 + "UMAP_seurat_clusters.png"
aa1 = link22231 + "UMAP_seurat_clusters.pdf"
plt.savefig(aa,bbox_inches='tight')
plt.savefig(aa1,bbox_inches='tight')

sc.set_figure_params(figsize = (15,15))
sc.pl.umap(adata, color="sample",size=60)
aa = link12231 + "UMAP_sample.png"
aa1 = link22231 + "UMAP_sample.pdf"
plt.savefig(aa,bbox_inches='tight')
plt.savefig(aa1,bbox_inches='tight')

#######################################分样本/cluster绘制umap图展示信息 
#20230113 发现bug 这里分的subdata  因为没跑过ir.tl.clonal_expansion(adata) 导致没有记录clonotype_size这些信息  所以这里要重新拆一遍subdata

a = {}
if len(Cluster) > 1:
  for i in Cluster:
    folder =  link11+"Cluster_"+i+"/clonotype_expansion/umapinfo/"
    folder1 = link21+"Cluster_"+i+"/clonotype_expansion/umapinfo/"
    subadata = adata[adata.obs['seurat_clusters'] == i]
    a[i] = subadata
    sc.set_figure_params(figsize = (15,15))
    sc.pl.umap(a[i], color="clonal_expansion",size=60)
    aa =  folder + "Cluster_"+ i + "_UMAP_clonal_expansion.png"
    aa1 = folder1 + "Cluster_"+ i + "_UMAP_clonal_expansion.pdf"
    plt.savefig(aa,bbox_inches='tight')
    plt.savefig(aa1,bbox_inches='tight')
    
    sc.set_figure_params(figsize = (15,15))
    sc.pl.umap(a[i], color="cc_aa_alignment_size",size=60)
    aa =   folder + "Cluster_"+ i + "_UMAP_clone_id_size.png"
    aa1 = folder1 + "Cluster_"+ i + "_UMAP_clone_id_size.pdf"
    plt.savefig(aa,bbox_inches='tight')
    plt.savefig(aa1,bbox_inches='tight')
    
    sc.set_figure_params(figsize = (15,15))
    sc.pl.umap(a[i], color="sample",size=60)
    aa =   folder + "Cluster_"+ i + "_UMAP_sample.png"
    aa1 = folder1 + "Cluster_"+ i + "_UMAP_sample.pdf"
    plt.savefig(aa,bbox_inches='tight')
    plt.savefig(aa1,bbox_inches='tight')

b = {}

if len(Sample) > 1:
  for j in Sample:
    folder =  link12+"Sample_"+j+"/clonotype_expansion/umapinfo/"
    folder1 = link22+"Sample_"+j+"/clonotype_expansion/umapinfo/"
    subadata = adata[adata.obs['sample'] == j]
    b[j] = subadata
    sc.set_figure_params(figsize = (15,15))
    sc.pl.umap(b[j], color="clonal_expansion",size=60)
    aa =  folder + "Sample_"+ j + "_UMAP_clonal_expansion.png"
    aa1 = folder1 + "Sample_"+ j + "_UMAP_clonal_expansion.pdf"
    plt.savefig(aa,bbox_inches='tight')
    plt.savefig(aa1,bbox_inches='tight')
    
    sc.set_figure_params(figsize = (15,15))
    sc.pl.umap(b[j], color="cc_aa_alignment_size",size=60)
    aa =   folder + "Sample_"+ j + "_UMAP_clone_id_size.png"
    aa1 = folder1 + "Sample_"+ j + "_UMAP_clone_id_size.pdf"
    plt.savefig(aa,bbox_inches='tight')
    plt.savefig(aa1,bbox_inches='tight')
    
    sc.set_figure_params(figsize = (15,15))
    sc.pl.umap(b[j], color="seurat_clusters",size=60)
    aa =   folder + "Sample_"+ j + "_UMAP_seurat_clusters.png"
    aa1 = folder1 + "Sample_"+ j + "_UMAP_seurat_clusters.pdf"
    plt.savefig(aa,bbox_inches='tight')
    plt.savefig(aa1,bbox_inches='tight')

############################################################用条形图分cluster展示拓展类别(即当前clonotype有多个copy)的比例 

#rds_prefix.png/Cluster/AllCluster/clonotype_expansion/clonal_expansion_stackedbarplot
#link11132
#rds_prefix.pdf/Cluster/AllCluster/clonotype_expansion/clonal_expansion_stackedbarplot
#link21132
#rds_prefix.png/Sample/AllSample/clonotype_expansion/clonal_expansion_stackedbarplot
#link12232
#rds_prefix.pdf/Sample/AllSample/clonotype_expansion/clonal_expansion_stackedbarplot
#link22232

if len(adata.obs['seurat_clusters'].cat.categories) > 21:
  kwargsdict = {'fig_kws':{'figsize': (18, 15), 'dpi': 300},'style_kws':{'title':'Clonal_expansion Distribution groupby Celltype','legend_title': "clonal expansion",'ylab':"Cellnum of Celltypes"}}
else:
  kwargsdict = {'fig_kws':{'figsize': (12, 10), 'dpi': 200},'style_kws':{'title':'Clonal_expansion Distribution groupby Celltype','legend_title': "clonal expansion",'ylab':"Cellnum of Celltypes"}}

ir.pl.clonal_expansion(adata, groupby="seurat_clusters", clip_at=7, normalize=False,target_col = "cc_aa_alignment",**kwargsdict)
aa = link11132 + "AllCluster_clonal_expansion_stackedbarplot.png"
plt.savefig(aa,bbox_inches='tight')
aa1 = link21132 + "AllCluster_clonal_expansion_stackedbarplot.pdf"
plt.savefig(aa1,bbox_inches='tight')

###

if len(adata.obs['sample'].cat.categories) > 21:
  kwargsdict = {'fig_kws':{'figsize': (18, 15), 'dpi': 300},'style_kws':{'title':'Clonal_expansion Distribution groupby Sample','legend_title': "clonal expansion",'ylab':"Cellnum of Samples"}}
else:
  kwargsdict = {'fig_kws':{'figsize': (12, 10), 'dpi': 200},'style_kws':{'title':'Clonal_expansion Distribution groupby Sample','legend_title': "clonal expansion",'ylab':"Cellnum of Samples"}}

ir.pl.clonal_expansion(adata, groupby="sample", clip_at=7, normalize=False,target_col = "cc_aa_alignment",**kwargsdict)
aa = link12232 + "AllSample_clonal_expansion_stackedbarplot.png"
plt.savefig(aa,bbox_inches='tight')
aa1 = link22232 + "AllSample_clonal_expansion_stackedbarplot.pdf"
plt.savefig(aa1,bbox_inches='tight')

###########################################################归一化

if len(adata.obs['seurat_clusters'].cat.categories) > 21:
  kwargsdict = {'fig_kws':{'figsize': (18, 15), 'dpi': 300},'style_kws':{'title':'Normlized Clonal_expansion Fraction groupby Celltype','legend_title': "clonal expansion",'ylab':"Fraction of Expanded Clonotype_clusters"}}
else:
  kwargsdict = {'fig_kws':{'figsize': (12, 10), 'dpi': 200},'style_kws':{'title':'Normlized Clonal_expansion Fraction groupby Celltype','legend_title': "clonal expansion",'ylab':"Fraction of Expanded Clonotype_clusters"}}

ir.pl.clonal_expansion(adata, clip_at=7, groupby="seurat_clusters",target_col = "cc_aa_alignment",**kwargsdict)

aa = link11132 + "AllCluster_clonotype_expansion_stackedbarplot_normalized.png"
plt.savefig(aa,bbox_inches='tight')
aa1 = link21132 + "AllCluster_clonotype_expansion_stackedbarplot_normalized.pdf"
plt.savefig(aa1,bbox_inches='tight')
#克隆扩增是某个效应T细胞克隆进行阳性选择的标志 

if len(adata.obs['sample'].cat.categories) > 21:
  kwargsdict = {'fig_kws':{'figsize': (18, 15), 'dpi': 300},'style_kws':{'title':'Normlized Clonal_expansion Fraction groupby Sample','legend_title': "clonal expansion",'ylab':"Fraction of Expanded Clonotype_clusters"}}
else:
  kwargsdict = {'fig_kws':{'figsize': (12, 10), 'dpi': 200},'style_kws':{'title':'Normlized Clonal_expansion Fraction groupby Sample','legend_title': "clonal expansion",'ylab':"Fraction of Expanded Clonotype_clusters"}}

ir.pl.clonal_expansion(adata, clip_at=7,groupby="sample",target_col = "cc_aa_alignment",**kwargsdict)
aa = link12232 + "AllSample_clonotype_expansion_stackedbarplot_normalized.png"
plt.savefig(aa,bbox_inches='tight')
aa1 = link22232 + "AllSample_clonotype_expansion_stackedbarplot_normalized.pdf"
plt.savefig(aa1,bbox_inches='tight')

######################################################分样本/cluster绘制

if len(Cluster) > 1:
  for i in Cluster:
    folder =  link11+"Cluster_"+i+"/clonotype_expansion/clonal_expansion_stackedbarplot/"
    folder1 = link21+"Cluster_"+i+"/clonotype_expansion/clonal_expansion_stackedbarplot/"
    
    if len(a[i].obs['sample'].cat.categories) > 21:
      kwargsdict = {'fig_kws':{'figsize': (18, 15), 'dpi': 300},'style_kws':{'title':'Celltype '+ i + ': Clonal_expansion Distribution groupby Sample','legend_title': "clonal expansion",'ylab':"Cellnum of Samples"}}
    else:
      kwargsdict = {'fig_kws':{'figsize': (12, 10), 'dpi': 200},'style_kws':{'title':'Celltype '+ i + ': Clonal_expansion Distribution groupby Sample','legend_title': "clonal expansion",'ylab':"Cellnum of Samples"}}
    
    ir.pl.clonal_expansion(a[i], groupby="sample", clip_at=7,normalize=False,target_col = "cc_aa_alignment",**kwargsdict)
    aa =   folder + "Cluster_" + i + "_clonal_expansion_stackedbarplot_groupby_sample.png"
    aa1 = folder1 + "Cluster_" + i + "_clonal_expansion_stackedbarplot_groupby_sample.pdf"   
    plt.savefig(aa,bbox_inches='tight')
    plt.savefig(aa1,bbox_inches='tight')
    
    if len(a[i].obs['sample'].cat.categories) > 21:
      kwargsdict = {'fig_kws':{'figsize': (18, 15), 'dpi': 300},'style_kws':{'title':'Celltype '+ i + ': Normlized Clonal_expansion Distribution groupby Sample','legend_title': "clonal expansion",'ylab':"Fraction of Expanded Clonotype_clusters"}}
    else:
      kwargsdict = {'fig_kws':{'figsize': (12, 10), 'dpi': 200},'style_kws':{'title':'Celltype '+ i + ': Normlized Clonal_expansion Distribution groupby Sample','legend_title': "clonal expansion",'ylab':"Fraction of Expanded Clonotype_clusters"}}
    
    ir.pl.clonal_expansion(a[i], groupby="sample", clip_at=7,normalize=True,target_col = "cc_aa_alignment",**kwargsdict)
    aa =   folder + "Cluster_" + i + "_clonotype_expansion_stackedbarplot_normalized_groupby_sample.png"
    aa1 = folder1 + "Cluster_" + i + "_clonotype_expansion_stackedbarplot_normalized_groupby_sample.pdf"   
    plt.savefig(aa,bbox_inches='tight')
    plt.savefig(aa1,bbox_inches='tight')

if len(Sample) > 1:
  for j in Sample:
    folder =  link12+"Sample_"+j+"/clonotype_expansion/clonal_expansion_stackedbarplot/"
    folder1 = link22+"Sample_"+j+"/clonotype_expansion/clonal_expansion_stackedbarplot/"
    
    if len(b[j].obs['seurat_clusters'].cat.categories) > 21:
      kwargsdict = {'fig_kws':{'figsize': (18, 15), 'dpi': 300},'style_kws':{'title':'Sample '+ j +': Clonal_expansion Distribution groupby Celltype','legend_title': "clonal expansion",'ylab':"Cellnum of Celltypes"}}
    else:
      kwargsdict = {'fig_kws':{'figsize': (12, 10), 'dpi': 200},'style_kws':{'title':'Sample '+ j +': Clonal_expansion Distribution groupby Celltype','legend_title': "clonal expansion",'ylab':"Cellnum of Celltypes"}}
    
    ir.pl.clonal_expansion(b[j], clip_at=7,groupby="seurat_clusters",normalize=False,target_col = "cc_aa_alignment")
    aa =   folder + "Sample_" + j + "_clonal_expansion_stackedbarplot_groupby_celltype.png"
    aa1 = folder1 + "Sample_" + j + "_clonal_expansion_stackedbarplot_groupby_celltype.pdf"   
    plt.savefig(aa,bbox_inches='tight')
    plt.savefig(aa1,bbox_inches='tight')
    
    if len(b[j].obs['seurat_clusters'].cat.categories) > 21:
      kwargsdict = {'fig_kws':{'figsize': (18, 15), 'dpi': 300},'style_kws':{'title':'Sample '+ j +': Normlized Clonal_expansion Distribution groupby Celltype','legend_title': "clonal expansion",'ylab':"Fraction of Expanded Clonotypes"}}
    else:
      kwargsdict = {'fig_kws':{'figsize': (12, 10), 'dpi': 200},'style_kws':{'title':'Sample '+ j +': Normlized Clonal_expansion Distribution groupby Celltype','legend_title': "clonal expansion",'ylab':"Fraction of Expanded Clonotypes"}}
    
    ir.pl.clonal_expansion(b[j], clip_at=7,groupby="seurat_clusters",normalize=True,target_col = "cc_aa_alignment")
    aa =   folder + "Sample_" + j + "_clonotype_expansion_stackedbarplot_normalized_groupby_celltype.png"
    aa1 = folder1 + "Sample_" + j + "_clonotype_expansion_stackedbarplot_normalized_groupby_celltype.pdf"   
    plt.savefig(aa,bbox_inches='tight')
    plt.savefig(aa1,bbox_inches='tight')


##################################################################

###############################################################计算每组细胞的clonotype多样性  

# 此处多样性越低的细胞分组 代表其免疫细胞克隆越多 这里可以只看cluster  
# normalized_shannon_entropy 这里的shannon_entropy是一种作为多样性度量的指标  按细胞分组归一化

#rds_prefix.png/Cluster/AllCluster/clonotype_expansion/clonotype_diversity
#link11133
#rds_prefix.pdf/Cluster/AllCluster/clonotype_expansion/clonotype_diversity
#link21133
#rds_prefix.png/Sample/AllSample/clonotype_expansion/clonotype_diversity
#link12233
#rds_prefix.pdf/Sample/AllSample/clonotype_expansion/clonotype_diversity
#link22233

##20230821发现 默认的香农系数的函数没有按照level顺序排列横坐标
# 扒一下代码重写一下函数 

############################################################# 绘制香农系数图

def _is_na2(x):
    return pd.isnull(x) or x in ("NaN", "nan", "None", "N/A", "")

_is_na = np.vectorize(_is_na2, otypes=[bool])

def _shannon_entropy(counts: np.ndarray):
    freqs = counts / np.sum(counts)
    np.testing.assert_almost_equal(np.sum(freqs), 1)
    
    if len(freqs) == 1:
        # the formula below is not defined for n==1
        return 0
    else:
        return -np.sum((freqs * np.log(freqs)) / np.log(len(freqs)))

def shannon_entropy(adata,groupby,target_col):
    ir_obs = adata.obs.loc[~_is_na(adata.obs[target_col]), :]
    clono_counts = ir_obs.groupby([groupby, target_col], observed=True).size().reset_index(name="count")
    # sample cc_aa_alignment  count 三列表 
    diversity = {}
    for k in list(adata.obs[groupby].cat.categories):
        tmp_counts = cast(
            np.ndarray,
            cast(pd.Series, clono_counts.loc[clono_counts[groupby] == k, "count"]).values,
        )
        diversity[k] = _shannon_entropy(tmp_counts)
    
    ## 把shannon entropy 添加到adata中 按照每个细胞所属的样本
    key_added = f"{'shannon_entropy_of_'}{target_col}{'_groupby_'}{groupby}"
    adata.obs[key_added] = adata.obs[groupby].map(diversity) 
    data = pd.DataFrame().from_dict(diversity, orient="index")
    
    #自己生成画图代码
    default_style_kws = {
            "title": "Alpha diversity of {} by {}".format(target_col,groupby),
            "ylab": f"{'Shannon_entropy of '}{target_col}{' groupby '}{groupby}"
    }  
    ax = ir.pl.base.bar(data, style_kws=default_style_kws,fig_kws = {'figsize': (15, 12), 'dpi': 240})
    # 删除图例 因为这里的图例由于没有染色 显示的是一个 0 很突兀 这里移除掉
    ax.get_legend().remove()
    # 刷新图形并移除图例的空白区域 
    plt.gcf().canvas.draw()
    
    return ax

ax = shannon_entropy(adata,"seurat_clusters","cc_aa_alignment")
aa = link11133 + "AllCluster_shannon_entropy_diversity.png"
plt.savefig(aa,bbox_inches='tight')
aa1 = link21133 + "AllCluster_shannon_entropy_diversity.pdf"
plt.savefig(aa1,bbox_inches='tight')

ax = shannon_entropy(adata,"sample","cc_aa_alignment")
aa = link12233 + "AllSample_shannon_entropy_diversity.png"
plt.savefig(aa,bbox_inches='tight')
aa1 = link22233 + "AllSample_shannon_entropy_diversity.pdf"
plt.savefig(aa1,bbox_inches='tight')

########################################################### 绘制giniindex图

# 20230404发现 在按照客户需求添加gini_index作为clonosize多样性指标时 自带的基尼系数由于输入数据没有从小到大排序
# 导致计算的值有负数 所以按照基尼系数定义写一个小函数

#计算基尼系数的函数
# def gini_coefficient(x):
#   x_sorted = np.sort(x)
#   cum = np.cumsum(x_sorted)
#   trapes = sum((cum[i] + cum[i+1]) / 2 for i in range(len(cum)-1))
#   triangle = len(x_sorted) * cum[-1] / 2
#   gini_index = (triangle - trapes) / triangle
#   return gini_index  

# 此处做出修改 改为使用skbio提供的gini函数 
# 自己写的这个 应该为0的时候只是无限接近0 计算的数值也有出入 
# 计算gini就用该用曲线下面积

def get_gini(adata,groupby,target_col):
  #这是生成clonosize列的过程
  # 抽出去除na的adata的obs_DataFrame
  ir_obs = adata.obs.loc[~_is_na(adata.obs[target_col]), :]
  
  clono_counts = (
          ir_obs.groupby([groupby, target_col], observed=True)
          .size()
          .reset_index(name="count")
  )
  #最终是抽出seurat_clusters clone_id  count三列 按照clone_id非na 然后对clone_id进行计数 
  
  diversity  = dict()
  for k in sorted(ir_obs[groupby].unique()):
    # k为cluster名字 抽出每个cluster下的每个clonotype的count列。
    # 结果是一列数字 这列数字用来计算每个cluster的giniindex 
    tmp_counts = clono_counts.loc[clono_counts[groupby] == k, "count"].values
    diversity[k] = skbio.diversity.alpha.gini_index(tmp_counts, method="rectangles")
  
  key_added = f"{'gini_index_of_'}{target_col}{'_groupby_'}{groupby}"
  adata.obs[key_added] = adata.obs[groupby].map(diversity) 
  data = pd.DataFrame().from_dict(diversity, orient="index")
  data = data.reindex(list(adata.obs[groupby].cat.categories))
  #自己生成画图代码
  default_style_kws = {
          "title": "Alpha diversity of {} by {}".format(target_col,groupby),
          "ylab": f"{'Gini_index of '}{target_col}{' groupby '}{groupby}"
  }
  ax = ir.pl.base.bar(data, style_kws=default_style_kws,fig_kws = {'figsize': (15, 12), 'dpi': 240})
  # 删除图例 因为这里的图例由于没有染色 显示的是一个 0 很突兀 这里移除掉
  ax.get_legend().remove()
  # 刷新图形并移除图例的空白区域 
  plt.gcf().canvas.draw()
  
  return ax

ax = get_gini(adata,"seurat_clusters","cc_aa_alignment")
aa = link11133 + "AllCluster_gini_index_diversity.png"
plt.savefig(aa,bbox_inches='tight')
aa1 = link21133 + "AllCluster_gini_index_diversity.pdf"
plt.savefig(aa1,bbox_inches='tight')

ax = get_gini(adata,"sample","cc_aa_alignment")
aa = link12233 + "AllSample_gini_index_diversity.png"
plt.savefig(aa,bbox_inches='tight')
aa1 = link22233 + "AllSample_gini_index_diversity.pdf"
plt.savefig(aa1,bbox_inches='tight')


#######################################################clonotype_size分析/克隆型丰度分析

#rds_prefix.png/Cluster/AllCluster/clonotype_expansion/clonetpye_size_abundance
#link11134
#rds_prefix.pdf/Cluster/AllCluster/clonotype_expansion/clonetpye_size_abundance
#link21134
#rds_prefix.png/Sample/AllSample/clonotype_expansion/clonetpye_size_abundance
#link12234
#rds_prefix.pdf/Sample/AllSample/clonotype_expansion/clonetpye_size_abundance
#link22234

##############################################################展示细胞最多的n种clonotype 和细胞分类在其中的分布 

ir.pl.group_abundance(adata, groupby="cc_aa_alignment", target_col="seurat_clusters", max_cols=15,fig_kws = {'figsize': (15, 12), 'dpi': 300})#只展示前15个

handles, labels = plt.gca().get_legend_handles_labels()
list1 = list(adata.obs['seurat_clusters'].cat.categories)
indices = [labels.index(x) for x in list1]
ncol1=(1 if len(list1) <= 21 else 2 if len(list1) <= 30 else 3)
place = 1.5
if ncol1 == 2:
  place = 1.8

if ncol1 == 3:
  place = 2.0

plt.legend([handles[i] for i in indices], [labels[i] for i in indices],loc = 'center right',frameon = False,title = "celltype",bbox_to_anchor  =(place,0.5),ncol=ncol1) 

aa = link11134 + "AllCluster_Top15_clonotype_size.png"
plt.savefig(aa,bbox_inches='tight')
aa1 = link21134 + "AllCluster_Top15_clonotype_size.pdf"
plt.savefig(aa1,bbox_inches='tight')

ir.pl.group_abundance(adata, groupby="cc_aa_alignment", target_col="sample", max_cols=15,fig_kws = {'figsize': (15, 12), 'dpi': 300})

handles, labels = plt.gca().get_legend_handles_labels()
list1 = list(adata.obs['sample'].cat.categories)
indices = [labels.index(x) for x in list1]
ncol1=(1 if len(list1) <= 21 else 2 if len(list1) <= 30 else 3)
place = 1.5
if ncol1 == 2:
  place = 1.8

if ncol1 == 3:
  place = 2.0

plt.legend([handles[i] for i in indices], [labels[i] for i in indices],loc = 'center right',frameon = False,title = "sample",bbox_to_anchor  =(place,0.5),ncol=ncol1) 

aa =  link12234 + "AllSample_Top15_clonotype_size.png"
plt.savefig(aa,bbox_inches='tight')
aa1 = link22234 + "AllSample_Top15_clonotype_size.pdf"
plt.savefig(aa1,bbox_inches='tight')

##############################################################分样本/cluster展示细胞最多的clonotype

if len(Cluster) > 1:
  for i in Cluster:  
    folder =  link11+"Cluster_"+i+"/clonotype_expansion/clonetpye_size_abundance/"
    folder1 = link21+"Cluster_"+i+"/clonotype_expansion/clonetpye_size_abundance/"
    ir.pl.group_abundance(a[i], groupby="cc_aa_alignment", target_col="sample", max_cols=15,fig_kws = {'figsize': (15, 12), 'dpi': 300})
    
    handles, labels = plt.gca().get_legend_handles_labels()
    list1 = list(a[i].obs['sample'].cat.categories)
    indices = [labels.index(x) for x in list1]
    
    ncol1=(1 if len(list1) <= 21 else 2 if len(list1) <= 30 else 3)
    place = 1.5
    if ncol1 == 2:
      place = 1.8
    
    if ncol1 == 3:
      place = 2.0
    
    plt.legend([handles[i] for i in indices], [labels[i] for i in indices],loc = 'center right',frameon = False,title = "sample",bbox_to_anchor  =(place,0.5),ncol=ncol1)
    
    aa =  folder + "Cluster_" + i + "_Top15_clonotype_size.png"
    aa1 = folder1 + "Cluster_" + i + "_Top15_clonotype_size.pdf"   
    plt.savefig(aa,bbox_inches='tight')
    plt.savefig(aa1,bbox_inches='tight')

if len(Sample) > 1:
  for j in Sample:
    folder =  link12+"Sample_"+j+"/clonotype_expansion/clonetpye_size_abundance/"
    folder1 = link22+"Sample_"+j+"/clonotype_expansion/clonetpye_size_abundance/"
    ir.pl.group_abundance(b[j], groupby="cc_aa_alignment", target_col="seurat_clusters", max_cols=15,fig_kws = {'figsize': (15, 12), 'dpi': 300})
    
    handles, labels = plt.gca().get_legend_handles_labels()
    list1 = list(b[j].obs['seurat_clusters'].cat.categories)
    indices = [labels.index(x) for x in list1]
    
    ncol1=(1 if len(list1) <= 21 else 2 if len(list1) <= 30 else 3)
    place = 1.5
    if ncol1 == 2:
      place = 1.8
    
    if ncol1 == 3:
      place = 2.0
    
    plt.legend([handles[i] for i in indices], [labels[i] for i in indices],loc = 'center right',frameon = False,title = "celltype",bbox_to_anchor = (place,0.5),ncol=ncol1)
    
    aa =   folder + "Sample_" + j + "_Top15_clonotype_size.png"
    aa1 = folder1 + "Sample_" + j + "_Top15_clonotype_size.pdf"   
    plt.savefig(aa,bbox_inches='tight')
    plt.savefig(aa1,bbox_inches='tight')

##############################################################按样本标化占比

# It might be beneficial to normalize the counts to the number of cells per sample to mitigate biases due to different sample sizes: 原文是这个 理解不了
# 不理解的是从数值标化为占比 这个占比是什么的占比 
# 最新的解释是 因为有的样本细胞数多 有的细胞数少 就会有librarysize 为了避免这个偏差 让每个细胞占比代表的含义更真实，这里是相当于给每个细胞一个相对于它本身sample的标化后的值，如：
# cell1 1/sample1num
# cell2 1/sample2num 
# 再用这个值加和成每个clonotype的总数 再在总数上按照细胞类型染色
# 这里举个栗子就是 cell1 + cell2 + cell3（假设cell3属于sample2）属于celltype1 那么他们代表的占比总数就是 (1/sample1num + 1/sample2num + 1/sample2num) 
# 所以在counts和freq表中 每列条形图的高度不同 celltype的占比也不同 
# 由于图中的数值都是0.1~0.5这种 推测是百分比。

ir.pl.group_abundance(
    adata, groupby="cc_aa_alignment", target_col="seurat_clusters", max_cols=15, normalize="sample",fig_kws = {'figsize': (15, 12), 'dpi': 300}
)

handles, labels = plt.gca().get_legend_handles_labels()
list1 = list(adata.obs['seurat_clusters'].cat.categories)
indices = [labels.index(x) for x in list1]
ncol1=(1 if len(list1) <= 21 else 2 if len(list1) <= 30 else 3)
place = 1.5
if ncol1 == 2:
  place = 1.8

if ncol1 == 3:
  place = 2.0

plt.legend([handles[i] for i in indices], [labels[i] for i in indices],loc = 'center right',frameon = False,title = "celltype",bbox_to_anchor  =(place,0.5),ncol=ncol1) 

aa = link11134 + "AllCluster_fraction_of_Top15clone_size_normalizedby_Sample.png"
plt.savefig(aa,bbox_inches='tight')
aa1 = link21134 + "AllCluster_fraction_of_Top15clone_size_normalizedby_Sample.pdf"
plt.savefig(aa1,bbox_inches='tight')

ir.pl.group_abundance(
    adata, groupby="cc_aa_alignment", target_col="sample", max_cols=15, normalize="seurat_clusters",figsize=(12, 6)
)

handles, labels = plt.gca().get_legend_handles_labels()
list1 = list(adata.obs['sample'].cat.categories)
indices = [labels.index(x) for x in list1]
ncol1=(1 if len(list1) <= 21 else 2 if len(list1) <= 30 else 3)
place = 1.5
if ncol1 == 2:
  place = 1.8

if ncol1 == 3:
  place = 2.0

plt.legend([handles[i] for i in indices], [labels[i] for i in indices],loc = 'center right',frameon = False,title = "sample",bbox_to_anchor  =(place,0.5),ncol=ncol1) 

aa = link12234 + "AllSample_fraction_of_Top15clone_size_normalizedby_Cluster.png"
plt.savefig(aa,bbox_inches='tight')
aa1 = link22234 + "AllSample_fraction_of_Top15clone_size_normalizedby_Cluster.pdf"
plt.savefig(aa1,bbox_inches='tight')

# 研究表明 有些clone是独属于某些组织的 有些是共同拥有的 但是患者/样本之间很少共享clonotype 这就解答了我前面一个疑问 最开始画clonotype图的时候完全没有占比的饼图 现在知道样本之间的clonotype基本是独立的

############################################################写出adata

#此处由于uns里还是记录了距离矩阵 所以不能用h5保存出来
obs = adata.obs
obs.index.name = 'cellnames'
aa = outpath1 + www + prefix + ".TCRcluster_cc_result.txt"
obs.to_csv(aa,sep = "\t")
aa = outpath1 + www + prefix + ".TCRcluster_cc_adata.dat"
pickle.dump(adata,open(aa,"wb"))
# loaded_model = pickle.load(open(aa,"rb"))

# 把带prefix的dat文件放在根目录(./) 

### 返回 0/223/233 注意这样会终止脚本运行
# python默认将脚本的最后一行作为返回值 也可以直接放一个223/233在这
sys.exit(returnvalue)
