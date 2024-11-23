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

sc.set_figure_params(figsize=(10, 10))
sc.settings.verbosity = 2  
#打印使用到的软件包
sc.logging.print_header()

# 20230110 本次新版主要功能：
# 首先 10X的filter_contig文件已经自带了clonotype 就是根据CDR3序列是否相同计算的   但是和nt+identity结果几乎相同 
# 另外 这个task只鉴定clonotype_cluster 简称cc

parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--prefix', type=str, default=None)#prefix
parser.add_argument('--path', type=str, default=None)#输出路径
parser.add_argument('--inputfilename', type=str, default = None)#aggr task里的filtered_contig_annotations.csv文件  

parser.add_argument('--receptor_arms2', type=str, default = None)
parser.add_argument('--dual_ir2', type=str, default = None)
parser.add_argument('--same_v_gene2', type=str, default = None)
parser.add_argument('--within_group2', type=str, default = None)
parser.add_argument('--min_nodes2', type=int, default = None)
parser.add_argument('--size_power2', type=float, default = None)
parser.add_argument('--min_cells2', type=int, default = None)
parser.add_argument('--base_size2', type=int, default = None)#可能为空值
parser.add_argument('--edges_width2', type=float, default = None)
parser.add_argument('--label_fontsize2', type=int, default = None)
parser.add_argument('--label_alpha2', type=float, default = None)
parser.add_argument('--cutoff2', type=int, default = None)

parser.add_argument('--receptor_type_cut', type=str, default = None)
parser.add_argument('--receptor_subtype_cut', type=str, default = None)
parser.add_argument('--chain_pairing_cut', type=str, default = None)

parser.add_argument('--threads', type=int, default = None)
parser.add_argument('--chunksize2', type=int, default = None)

args = parser.parse_args()
prefix = args.prefix
path1 = args.path
inputfilename = args.inputfilename

receptor_arms2 = args.receptor_arms2
dual_ir2 = args.dual_ir2
same_v_gene2 = args.same_v_gene2 
within_group2 = args.within_group2
min_nodes2 = args.min_nodes2
size_power2 = args.size_power2
min_cells2 = args.min_cells2
base_size2 = args.base_size2
edges_width2 = args.edges_width2
label_fontsize2 = args.label_fontsize2
label_alpha2 = args.label_alpha2
cutoff2 = args.cutoff2

receptor_type_cut = args.receptor_type_cut
receptor_subtype_cut = args.receptor_subtype_cut
chain_pairing_cut = args.chain_pairing_cut

threads = args.threads
chunksize2 = args.chunksize2

if same_v_gene2 == 'true':
  same_v_gene2 = True
if same_v_gene2 == 'false':
  same_v_gene2 = False

############################################QC参数组 

group = receptor_type_cut.split(',')
group1 = []
for i in group:
  t = i.replace('_', ' ')
  group1.append(t)
receptor_type_cut = group1

group = receptor_subtype_cut.split(',')
group1 = []
for i in group:
  t = i.replace('_', ' ')
  group1.append(t)
receptor_subtype_cut = group1

group = chain_pairing_cut.split(',')
group1 = []
for i in group:
  t = i.replace('_', ' ')
  group1.append(t)
chain_pairing_cut = group1

# chain_pairing_cut = ['ambiguous','multichain','no IR']
# receptor_subtype_cut = ['TRA+TRB','TRG+TRD','multichain','ambiguous','no IR']
# receptor_type_cut = ['multichain','ambiguous','TCR','no IR']

#########################################################加载ir信息 
#注意 这里读取文件 没有带入raw_clonotype_id 就是10X根据CDR3序列算的根本没记录在里面 

# inputfilename = "/media/nbc1/lijialun/scirpy/zhenbo/inputfile/ALL.csv"
# path1 = "/media/nbc1/lijialun/scirpy/zhenbo/TCRNrom/cc"  
# prefix = 'tcr'
# receptor_arms2 = 'all'
# dual_ir2 = "any"
# within_group2 = 'receptor_type'
# min_cells2 = 12
# min_nodes2 = 3
# base_size2 = None
# size_power2 = 0.6
# label_fontsize2 = 9
# label_alpha2 = 0.6
# edges_width2 = 0.5
# cutoff2 =  10

####################################################

adata = ir.io.read_10x_vdj(inputfilename)

#从输入文件中抽样本名 给整个dataframe去重复 直接添加到obs
irfile = pd.read_csv(inputfilename,sep = ',')
if "origin" in irfile.columns:
  irfile = irfile.drop_duplicates("barcode")
  irfile = irfile[['barcode',"origin"]]
  irfile.index = irfile['barcode']
  irfile.index.name = None
  irfile = irfile.reindex(list(adata.obs.index))
  adata.obs['orig.ident'] = irfile['origin'].astype('category')
else:
  adata.obs['orig.ident'] = prefix

#################################################################新建结果文件夹

#20230130 发现 实际使用要把dat文件写在根目录 但是我三个task传的都是prefix命名的文件夹 所以要改成传根目录进来 然后还要传prefix名字 在这里新建一个prefix名字命名的文件夹 

www = "/"
path = path1

os.mkdir(path1+www+"cc")
path1 = path1+www+"cc"+www

os.mkdir(path1+www+prefix)
path1 = path1+www+prefix+www


#设计一个单独保存pdf的文件夹
www = "/"
pdfgraph = "/ALLPDF"
os.mkdir(path1+pdfgraph)
pdfgraph = path1 + pdfgraph + www


#qc过程中展示ir基本信息的图
www = "/"
qcgraph = "/QC_Graph"
os.mkdir(path1+qcgraph)
qcgraph_png = path1 + qcgraph + www

os.mkdir(pdfgraph+"/QC_Graph")
qcgraph_pdf = pdfgraph + qcgraph + www

os.mkdir(qcgraph_pdf+"/QC_raw")
os.mkdir(qcgraph_pdf+"/QC_filter")
os.mkdir(qcgraph_png+"/QC_raw")
os.mkdir(qcgraph_png+"/QC_filter")

link1 = qcgraph_png + "/QC_raw/"
link2 = qcgraph_pdf + "/QC_raw/"
link3 = qcgraph_png + "/QC_filter/"
link4 = qcgraph_pdf + "/QC_filter/"

#鉴定clonotype/clonotype_cluster的图
dcgraph = "/definite_clonotype_cluster"
os.mkdir(path1+dcgraph)
dcgraph_png = path1 + dcgraph + www

os.mkdir(pdfgraph+dcgraph)
dcgraph_pdf = pdfgraph + dcgraph + www

################################################################## qc

ir.tl.chain_qc(adata)

################################################################质控图

ax = ir.pl.group_abundance(adata, groupby="receptor_subtype", target_col="orig.ident",fig_kws = {'figsize': (12, 9), 'dpi': 300})
aa =  link1 + "QC_receptor_subtype_Sample.png"
aa1 = link2 + "QC_receptor_subtype_Sample.pdf"   
plt.savefig(aa,bbox_inches='tight')
plt.savefig(aa1,bbox_inches='tight')

ax = ir.pl.group_abundance(adata, groupby="receptor_type", target_col="orig.ident",fig_kws = {'figsize': (12, 9), 'dpi': 300})
aa =  link1 + "QC_receptor_type_Sample.png"
aa1 = link2 + "QC_receptor_type_Sample.pdf"   
plt.savefig(aa,bbox_inches='tight')
plt.savefig(aa1,bbox_inches='tight')

ax = ir.pl.group_abundance(adata, groupby="chain_pairing",target_col="orig.ident",fig_kws = {'figsize': (12, 9), 'dpi': 300})

aa =  link1 + "QC_chain_pairing_Sample.png"
aa1 = link2 + "QC_chain_pairing_Sample.pdf"   
plt.savefig(aa,bbox_inches='tight')
plt.savefig(aa1,bbox_inches='tight')

#包含超过一对t细胞受体的细胞占比：
print("Fraction of cells with more than one pair of BCRs: {:.2f}".format(np.sum(adata.obs["chain_pairing"].isin(["extra VJ", "extra VDJ", "two full chains"]))/ adata.n_obs))

#过滤
# 20230130 bcr和tcr过滤不太一样 10X和BD的过滤策略也不太相同 
# 10X的BCR 在receptor_subtype中出现了只有IGH的选项 而TCR只有配对的选项 只有BCR把只有重链的情况单独作为一个过滤条目 
# 所以决定把每个过滤条目都列出来 做下拉框筛选 

# receptor_type可以是以下之一
# TCR（包含 TRA/TRB/TRG/TRD 链的任意组合但不含 IGH/IGK/IGL 链的所有细胞）
# BCR（包含 IGH/IGK/IGL 链的任意组合但不含 TCR 链的所有细胞）
# ambiguous（所有包含 BCR 和 TCR 链的细胞）
# multichain（所有具有两个以上 VJ 或两个以上 VDJ 链的细胞）
# no IR（所有未检测到免疫受体的细胞）

# receptor_subtype可以是以下之一
# TRA+TRB（所有只有 TRA 和/或 TRB 链的细胞）
# TRG+TRD（所有只有 TRG 和/或 TRD 链的细胞）
# IGH（所有只有 IGH 链但没有 IGL 或 IGK 的细胞）
# IGH+IGL（所有只有 IGH 和 IGL 链的细胞）
# IGH+IGK（所有只有 IGH 和 IGK 链的细胞）
# multichain（所有具有两个以上 VJ 或两个以上 VDJ 链的细胞）
# ambiguous（所有非上述类型的细胞，例如 TRA+TRD、TRA+IGH 或 IGH+IGK 作为主要受体，IGH+IGL 作为次要受体）
# no IR（所有未检测到免疫受体的细胞）

# chain_pairing可以是以下之一
# single pair（所有具有完全匹配的 VJ 和 VDJ 链的单元格）
# orphan VJ（所有只有一个 VJ 链的单元格）
# orphan VDJ（所有只有一个 VDJ 链的细胞）
# extra VJ（所有具有一对匹配的 VJ 和 VDJ 链以及一个额外的 VJ 链的细胞）
# extra VDJ（类比）
# two full chains（所有具有两对匹配的 VJ 和 VDJ 链的细胞）
# ambiguous（所有具有不匹配链的细胞，即已被归类为ambiguous receptor_subtype 的细胞）
# multichain（所有具有两个以上 VJ 或两个以上 VDJ 链的细胞）
# no IR（所有没有免疫受体链的链）

#20230130 发现 在chain_pairing中出现的orphan VJ比在receptor_subtype出现的ambiguous多 
#最后发现 IGH+IGL/IGH+IGK并不表示配对链 而是包含只有IGL/IGK单链的细胞+IGL/IGK与IGH配对的细胞 单独只有IGH单链的
#至于TCR 并没有单独把只有VDJ链的算作单独的过滤条目 所以差异不明显

#20230130 定义参数: 
# receptor_type_cut 表示删除receptor_type下的某些项目 多选下拉框 逗号隔开 
# TCR 仅包含TCR链组合且不包含BCR链
# BCR 仅包含BCR链组合且不包含TCR链
# ambiguous 同时包含TCR BCR链
# multichain 具有两个以上VJ或两个以上VDJ链
# no IR 未检测到免疫受体的细胞

# receptor_subtype_cut 表示删除receptor_subtype下的某些项目 多选下拉框 逗号隔开 
# TRA+TRB 包含TRA/TRB/TRA+TRB链
# TRG+TRD 包含TRG/TRD/TRG+TRD链
# IGH 仅包含TGH链
# IGH+IGL 包含IGL/TGL+IGH链
# IGH+IGK 包含IGK/IGK+IGH链
# multichain具有两个以上VJ或两个以上VDJ链
# ambiguous 非上述组合的其他细胞(具体来说TRA+TRD,不同类的TCR受体错配；TRA+IGH,BCR和TCR受体错配；同时包含IGL+IGh与IGK+IGH并不会被算作双受体对 证据是ambiguous和two full chains内的细胞没有交集)
# no IR 未检测到免疫受体的细胞

# chain_pairing_cut 表示删除chain_pairing下的某些项目 多选下拉框 逗号隔开 
# single pair 具有一对完全匹配的VJ和VDJ链
# orphan VJ 只有一个VJ链
# orphan VDJ 只有一个VDJ链
# extra VJ 具有一对匹配的VJ和VDJ链和一个额外VJ链
# extra VDJ 具有一对匹配的VJ和VDJ链和一个额外VDJ链
# two full chains 具有两对匹配的VJ和VDJ链
# ambiguous 被归类为receptor_subtype_ambiguous的细胞
# multichain 具有两个以上VJ或两个以上VDJ链
# no IR 未检测到免疫受体的细胞

#20230131发现 参数传递没法传带空格的字符串  所以参数在task界面传进来要用下划线代替空格 所以在脚本里要写替换 把_替换为空格 

#过滤
#20230111 几个问题 
# 首先 BD两万细胞如果删除只有VJ或者只有VDJ的细胞 直接变成一万了 
# 其次 BD有多种奇怪的组合 例如a链匹配d链 b链匹配G链 等等 有三十多种 这些不常见的组合有些文章说确实存在 
# 就和双受体类似 学术界没有共识
# 这个输入文件表里目前只标了TRA+TRB 后续会把多种链组合拆分掉重新画图 
# 所以对于BD数据 一定不能干掉单链细胞  
# 保留单链细胞会导致少了一大部分clonotype 因为距离矩阵是VJ和VDJ分开算 鉴定也是 这样只测到一条链的细胞会被归类到某些clonotype中 
# 好处是和单转联合的时候overlap较多 但是风险较大 准确率不够
# 对于10X的数据 6500细胞只删掉500左右的单链 具体这里要不要删除还是要等技术部调研 目前是固定不删除 
# 另外 BD的测序也是包含三条链或者四条链的数据 但是每种只保留了一条 
adata = adata[~adata.obs["receptor_type"].isin(receptor_type_cut), :].copy()
adata = adata[~adata.obs["receptor_subtype"].isin(receptor_subtype_cut), :].copy()
adata = adata[~adata.obs["chain_pairing"].isin(chain_pairing_cut), :].copy()


#正常来说 10X的下机数据不删除单链细胞 不删除ambiguous 就什么不会被过滤掉 

#过滤后重新绘制展示图
ax = ir.pl.group_abundance(adata, groupby="chain_pairing",target_col="orig.ident",fig_kws = {'figsize': (12, 9), 'dpi': 300})
aa =  link3 + "QC_filter_chain_pairing_Sample.png"
aa1 = link4 + "QC_filter_chain_pairing_Sample.pdf"   
plt.savefig(aa,bbox_inches='tight')
plt.savefig(aa1,bbox_inches='tight')

ax = ir.pl.group_abundance(adata, groupby="receptor_type",target_col="orig.ident",fig_kws = {'figsize': (12, 9), 'dpi': 300})
aa =  link3 + "QC_filter_receptor_type_Sample.png"
aa1 = link4 + "QC_filter_receptor_type_Sample.pdf"   
plt.savefig(aa,bbox_inches='tight')
plt.savefig(aa1,bbox_inches='tight')

ax = ir.pl.group_abundance(adata, groupby="receptor_subtype",target_col="orig.ident",fig_kws = {'figsize': (12, 9), 'dpi': 300})
aa =  link3 + "QC_filter_receptor_subtype_Sample.png"
aa1 = link4 + "QC_filter_receptor_subtype_Sample.pdf"   
plt.savefig(aa,bbox_inches='tight')
plt.savefig(aa1,bbox_inches='tight')

#########################################################################鉴定clonotype

# 没有暴露的参数：
# 1.metric 表示对于距离的度量算法 ci默认为identity  cc默认为alignment
# identity  1表示相同的序列，否则为0 该指标意味着cutoff为0  0意味着cutoff参数无意义 这里建议每种度量都使用其默认的cutoff 
# alignment 原理是给序列对打分 这个对齐距离有个算法 具体需要查看替换矩阵(蛋白质之间的替换能够在打分矩阵中找到对应的分数) 对齐距离是实际对齐分数与最大值之间的差异。
# 最大值是两条链自己和自己比对的分数中的最小值  Min(S1vs1 S2vs2) - S1vs2
# 替换矩阵名为BLOSUM62矩阵 每个蛋白之间的替换都有打分 这是非常常用的序列比对表
#    A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
# A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
# R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
# N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
# D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
# C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
# Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
# E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
# G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
# H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
# I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
# L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
# K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
# M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
# F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
# P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
# S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 
# T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
# W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
# Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
# V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
# B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
# Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
# X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 
# * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 
# 暂时没有明确的判断标准表明距离超过多少就不再识别用一抗原 不过不建议超过15 此处默认为10 
# 2.cutoff 参数表示任何distance>cutoff值的都会被干掉 但是identity度量方式导致这个参数无意义 在这个情况下默认忽略掉cutoff并且cutoff为0 对于其他情况直接设置为None会按照每种情况的默认阈值设置。alignment算法要看情况调整
# 3.sequence ci默认为nt cc默认为aa 即使用氨基酸序列还是蛋白质序列来鉴定clonotype 


#20230111  BD这组数据细胞太多  卡mincell=3全部重叠在一起 试了一下mincell卡4  好了一点但是还是很多 另外发现有很多1个细胞的点 需要卡一下min_nodes
# 这里之所以显得很混乱是因为 mincell卡的是cc总的细胞数 但是minnodes卡的是每个点 这里很多只有一个细胞的受体配置跟大的点之间有连线 导致图很混乱 

# 鉴定clonotype_cluster默认参数：
# within_group参数 如果设为redeptor_type 含义应该是只能有TCR或者BCR 设为 receptor_subtype这种可能就是clonotype不能共存于ab链和dk链之间 卡的太细了
# 20230206 研究发现 在鉴定cc时：
# 鉴定cc的步骤首先是
# 1.根据CDR3序列区分类型(包含相同CDR3序列的细胞会被折叠起来，算作同一种独特的受体配置)
# 2.分别计算VJ和VDJ的距离矩阵。在这个过程中引入aligment算法，cutoff表示算法的阈值，该算法会根据蛋白质的替换给出打分，默认为10，实测发现10能够很好地区分不同长度的氨基酸序列
# 20230207 建议10X为15 因为结果很好 BD为10 BD单链过多 避免产生误判
# 此处距离大于cutoff的匹配结果会被删除。
# 3.根据两个距离矩阵和receptor_arms2，dual_ir2两个参数匹配各个独特受体配置。
# 而receptor_arms2表示在判定cc匹配结果时，是只考虑VJ/VDJ还是全部都要考虑(all any)
# 但是all和any含义有所不同 根据上文得知 VJ和VDJ会给出两种距离(暂时不考虑多对受体，多受体由dual_ir控制)
# 此处要引入一个假设：两条受体间距离小于某个值时，可以识别同种抗原；那么any代表只有一个受体匹配就可以，all表示VJ和VDJ都要匹配。
# 所以any是用两部分距离的最小值来匹配受体对，all是用两部分距离的最大值来匹配。

# 至于参数dual_ir，则表示只考虑最丰富的一对受体还是考虑两对(all any)

# 此处和cc出现了歧义  dual_ir在ci中设置为primary_only 表示只看第一对链 因为第一对链表示更可信
# cc图上每个点代表独特的受体配置 即所有链相同的才会被视为同一个点 而ci鉴定只看第一对链 这样导致同个ci在cc图可能分开成两个点 
#至于ci图为什么没有带连线的点 因为只考虑第一对受体 区分受体配置的时候就不会因为多链的情况分成多个点。

# 此处all表示两对受体都要小于距离矩阵的截止值，any表示任意一対满足要求即可。
# 注意 此处同样适用于单链受体和一対多一条链的类型。
# 此处两个默认参数较为合理，receptor_arms = 'all'表示把VJ和VDJ的结果结合起来鉴定cc；
# dual_ir = "any"表示双受体之一匹配即可被归类(因为实际生物学过程中没人知道使用哪种受体识别了抗原，且2+1类型也不一定代表多余一条单链，也可能是单链的匹配链没测到)
# 4.上面匹配的结果会定义成一个点图 根据图中点的距离查找连接的点，连接起来的点被认为是一个cc。

# BD由于单链受体过多，会聚类成很大群体的cc(400细胞左右)，并且cutoff需要严格控制；10X细胞数最多的一个cc也才20几个细胞，包含的ci种类最多也就四个。
# 后续绘制sequence图的时候要注意区分这里

# same_v_gene 是否强制包含完全相同的v基因 10X结果带有cdr1和cdr2 是否要囊括进去 
# within_group 可以选择obs里任意一列 默认为receptor_type 目前是选择QC三个参数其中之一 
# label_fontsize 每个克隆型标签的大小 建议默认值为9
# label_alpha 标签透明度 默认0.6
# base_size 表示每一个细胞代表的基础大小 默认的话会使用tl里输入的值（在tl里这个参数是计算时候自动赋值的）一般不更改
# size_power 细胞数目代表的点大小的权重 一般更改这个调整点重叠问题
# edges_width 连线的宽度 当前设置只有cc有这参数 默认为0.4 

# receptor_arms2 = 'all'
# dual_ir2 = "any"
# within_group2 = 'receptor_type'
# min_cells2 = 5
# min_nodes2 = 3
# base_size2 = None
# size_power2 = 0.9
# label_fontsize2 = 9
# label_alpha2 = 0.6
# edges_width2 = 0.5
# cutoff2 =  10

#20230111  BD这组数据细胞太多  卡mincell=3全部重叠在一起 试了一下mincell卡4  好了一点但是还是很多 另外发现有很多1个细胞的点 需要卡一下min_nodes
# 这里之所以显得很混乱是因为 mincell卡的是cc总的细胞数 但是minnodes卡的是每个点 这里很多只有一个细胞的受体配置跟大的点之间有连线 导致图很混乱 

########################################################计算distance矩阵

ir.pp.ir_dist(adata,metric = 'alignment',cutoff = cutoff2,sequence = 'aa',n_jobs = threads)

#这里全部是默认参数 那么存储位置为 uns中的'ir_dist_aa_alignment'

#################################################################鉴定cc

#这里是调用ir_dist_aa_alignment 存储为obs和uns中的 cc_aa_alignment
#20230228发现 此处注意 n_jobs参数可以控制这个函数使用多少个cpu
# 默认是使用全部 在私有云测试时占用了全部的cpu去跑 直接卡死
# chunksize参数表示一个线程处理几个任务 默认是2000 在细胞数小于2*chunksize的时候 只用一个cpu
ir.tl.define_clonotype_clusters(
  adata, 
  sequence='aa',
  metric='alignment', 
  receptor_arms=receptor_arms2, 
  dual_ir=dual_ir2,
  inplace = True,
  same_v_gene = False,
  within_group = within_group2,
  distance_key = "ir_dist_aa_alignment",
  n_jobs = threads,
  chunksize = chunksize2,
)


ir.tl.clonotype_network(
  adata, 
  min_cells=min_cells2, 
  sequence='aa', 
  metric='alignment',
  min_nodes = min_nodes2,
  size_power = size_power2,
  base_size = base_size2,
  inplace = True,
  clonotype_key  = "cc_aa_alignment"
  )

#########################################################绘制clonotype_cluster布局图

ir.pl.clonotype_network(
  adata, 
  color="orig.ident",
  label_fontsize=label_fontsize2,
   panel_size=(12, 10),
   label_alpha = label_alpha2,
   base_size = base_size2,
   edges_width = edges_width2
   )

aa =  dcgraph_png + "clonotype_cluster_network_Sample.png"
aa1 = dcgraph_pdf + "clonotype_cluster_network_Sample.pdf"   
plt.savefig(aa,bbox_inches='tight')
plt.savefig(aa1,bbox_inches='tight')

#############################################################设定V基因相同后 观察cc变化

ir.tl.define_clonotype_clusters(
  adata, 
  sequence="aa",
  metric="alignment", 
  receptor_arms=receptor_arms2, 
  dual_ir=dual_ir2,
  inplace = True,
  same_v_gene = True,
  within_group = within_group2,
  key_added="cc_aa_alignment_same_v"
)

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

aa = dcgraph_png + "clonotypes_cluster_difference.txt"
ct_different.to_csv(aa,sep='\t',index = False)

############################################################写出adatair 保存obs文件

obs = adata.obs

obs.index.name = 'cellnames'
aa = path + www + prefix + ".scirpyBCRNorm_cc_result.txt"
obs.to_csv(aa,sep = "\t")

aa = path + www + prefix + ".BCR_adatair_cc.dat"
pickle.dump(adata,open(aa,"wb"))
# loaded_model = pickle.load(open(aa,"rb"))

