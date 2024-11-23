#20221216 task目标 在不引入单细胞数据的情况下 进行qc 鉴定clonetype和cc 另外 对于样本名 csv文件中有一列写了这个信息 把它读取出来 放到obs["orig.ident"]里
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

sc.set_figure_params(figsize=(10, 10))
sc.settings.verbosity = 2  
#打印使用到的软件包
sc.logging.print_header()

parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--path', type=str, default=None)#输出路径
parser.add_argument('--inputfilename', type=str, default = None)#aggr task里的filtered_contig_annotations.csv文件  
parser.add_argument('--metric1', type=str, default = None)
parser.add_argument('--metric2', type=str, default = None)
parser.add_argument('--sequence1', type=str, default = None)
parser.add_argument('--sequence2', type=str, default = None)
parser.add_argument('--cutoff1', type=int, default = None)#可能为空值 尝试了一下 不输入任何值 到python里就是None
parser.add_argument('--cutoff2', type=int, default = None)
parser.add_argument('--receptor_arms1', type=str, default = None)
parser.add_argument('--receptor_arms2', type=str, default = None)
parser.add_argument('--dual_ir1', type=str, default = None)
parser.add_argument('--dual_ir2', type=str, default = None)
parser.add_argument('--same_v_gene1', type=str, default = None)
parser.add_argument('--same_v_gene2', type=str, default = None)
parser.add_argument('--within_group1', type=str, default = None)
parser.add_argument('--within_group2', type=str, default = None)
parser.add_argument('--min_nodes1', type=int, default = None)
parser.add_argument('--min_nodes2', type=int, default = None)
parser.add_argument('--size_power1', type=float, default = None)
parser.add_argument('--size_power2', type=float, default = None)
parser.add_argument('--min_cells1', type=int, default = None)
parser.add_argument('--min_cells2', type=int, default = None)
parser.add_argument('--base_size1', type=int, default = None)#可能为空值
parser.add_argument('--base_size2', type=int, default = None)
parser.add_argument('--edges_width1', type=float, default = None)
parser.add_argument('--edges_width2', type=float, default = None)
parser.add_argument('--label_fontsize1', type=int, default = None)
parser.add_argument('--label_fontsize2', type=int, default = None)
parser.add_argument('--label_alpha1', type=float, default = None)
parser.add_argument('--label_alpha2', type=float, default = None)


args = parser.parse_args()
path1 = args.path
inputfilename = args.inputfilename
metric1 = args.metric1
metric2 = args.metric2
sequence1 = args.sequence1
sequence2 = args.sequence2
cutoff1 = args.cutoff1
cutoff2 = args.cutoff2
receptor_arms1 = args.receptor_arms1
receptor_arms2 = args.receptor_arms2
dual_ir1 = args.dual_ir1
dual_ir2 = args.dual_ir2
same_v_gene1 = args.same_v_gene1 
same_v_gene2 = args.same_v_gene2
within_group1 = args.within_group1
within_group2 = args.within_group2
min_nodes1 = args.min_nodes1
min_nodes2 = args.min_nodes2
size_power1 = args.size_power1
size_power2 = args.size_power2
min_cells1 = args.min_cells1
min_cells2 = args.min_cells2
base_size1 = args.base_size1
base_size2 = args.base_size2
edges_width1 = args.edges_width1
edges_width2 = args.edges_width2
label_fontsize1 = args.label_fontsize1
label_fontsize2 = args.label_fontsize2
label_alpha1 = args.label_alpha1
label_alpha2 = args.label_alpha2

if same_v_gene1 == 'true':
  same_v_gene1 = True
if same_v_gene1 == 'false':
  same_v_gene1 = False

if same_v_gene2 == 'true':
  same_v_gene2 = True
if same_v_gene2 == 'false':
  same_v_gene2 = False

#################################################################加载数据 

# inputfilename = "/media/nbc1/lijialun/scirpy/inputfile/tcr_filtered_contig_annotations.csv"
# path1 = "/media/nbc1/lijialun/scirpy/TCRNorm"  

#加载ir信息  来源自cellranger T/B分开做
adata = ir.io.read_10x_vdj(inputfilename)

#从输入文件中抽样本名 给整个dataframe去重复 直接添加到obs
irfile = pd.read_csv(inputfilename,sep = ',')
irfile = irfile.drop_duplicates("barcode")
irfile = irfile[['barcode',"origin"]]
irfile.index = irfile['barcode']
adata.obs['orig.ident'] = irfile['origin']

#################################################################新建结果文件夹

#qc过程中展示ir基本信息的图
www = "/"
qcgraph = "/QC_Graph"
os.mkdir(path1+qcgraph)
qcgraph = path1 + qcgraph + www

#鉴定clonetype/clonetype_cluster的图
dcgraph = "/definite_clonetype"
os.mkdir(path1+dcgraph)
dcgraph = path1 + dcgraph + www

################################################################对ir信息做质控

ir.tl.chain_qc(adata)

################################################################质控图

ax = ir.pl.group_abundance(adata, groupby="receptor_subtype", target_col="orig.ident",fig_kws = {'figsize': (12, 9), 'dpi': 300})
aa = qcgraph + "qc_group_abundance_receptor_subtype_Sample.png"
aa1 = qcgraph + "qc_group_abundance_receptor_subtype_Sample.pdf"   
plt.savefig(aa,bbox_inches='tight')
plt.savefig(aa1,bbox_inches='tight')

ax = ir.pl.group_abundance(adata, groupby="receptor_type", target_col="orig.ident",fig_kws = {'figsize': (12, 9), 'dpi': 300})
aa = qcgraph + "qc_group_abundance_receptor_type_Sample.png"
aa1 = qcgraph + "qc_group_abundance_receptor_type_Sample.pdf"   
plt.savefig(aa,bbox_inches='tight')
plt.savefig(aa1,bbox_inches='tight')

ax = ir.pl.group_abundance(adata, groupby="chain_pairing",target_col="orig.ident",fig_kws = {'figsize': (12, 9), 'dpi': 300})

aa =  qcgraph + "qc_group_abundance_chain_pairing_Sample.png"
aa1 = qcgraph + "qc_group_abundance_chain_pairing_Sample.pdf"   
plt.savefig(aa,bbox_inches='tight')
plt.savefig(aa1,bbox_inches='tight')

#包含超过一对t细胞受体的细胞占比：
print("Fraction of cells with more than one pair of TCRs: {:.2f}".format(np.sum(adata.obs["chain_pairing"].isin(["extra VJ", "extra VDJ", "two full chains"]))/ adata.n_obs))

#在此处对这些细胞进行过滤 
adata = adata[adata.obs["receptor_type"] != "multichain", :].copy()
adata = adata[adata.obs["receptor_type"] != "ambiguous", :].copy()
adata = adata[adata.obs["receptor_type"] != "BCR", :].copy()

#过滤 只含有vj和只含有vdj的被过滤掉了

adata = adata[~adata.obs["chain_pairing"].isin(["orphan VDJ", "orphan VJ"]), :].copy()
#示例数据中 删掉了600左右的细胞

#过滤后重新绘制展示图
ax = ir.pl.group_abundance(adata, groupby="chain_pairing",target_col="orig.ident",fig_kws = {'figsize': (12, 9), 'dpi': 300})
aa =  qcgraph + "qc_group_abundance_re_chain_pairing_Sample.png"
aa1 = qcgraph + "qc_group_abundance_re_chain_pairing_Sample.pdf"   
plt.savefig(aa,bbox_inches='tight')
plt.savefig(aa1,bbox_inches='tight')

#############################################################################################鉴定clonetype

# #严格鉴定 的默认参数：
# metric1 = 'identity'
# cutoff1 = None
# sequence1 = 'nt'
# #define_clonotypes参数
# receptor_arms1 = 'all'
# dual_ir1 = "primary_only"
# within_group1 = 'receptor_type'
# same_v_gene1 = False
# #clonotype_network
# min_cells1 = 2
# min_nodes1 = 1
# base_size1 = None
# size_power1 = 0.9
# label_fontsize1 = 9
# label_alpha1 = 0.6
# edges_width1 = 0.5

########################################################计算distance矩阵

ir.pp.ir_dist(adata,metric = metric1,cutoff = cutoff1,sequence = sequence1)

distance_key1 = "ir_dist_"+ sequence1 + "_" + metric1

#########################################################鉴定clonetype

# 根据CDR3核酸序列同一性严格定义克隆型
ir.tl.define_clonotypes(
  adata, 
  receptor_arms=receptor_arms1, 
  dual_ir=dual_ir1,
  inplace = True,
  same_v_gene = same_v_gene1,
  within_group = within_group1,
  distance_key = distance_key1
)
#20221214检查发现 clonetype_cluster函数中存在一个叫distance_key的参数  ir_dist计算的距离矩阵存在adata.uns中 并且是默认按照规则命名的 
#而clonetype中如果不调用这个key参数 还是会按照最基础的ir_dist_nt_identity来获取

######################################################建立clonetype关系图

ir.tl.clonotype_network(
  adata, 
  sequence = sequence1,
  metric = metric1,
  min_cells=min_cells1,
  min_nodes = min_nodes1,
  base_size = base_size1,
  size_power = size_power1,
  clonotype_key = "clone_id",
  inplace = True
)
#20221214检查发现 此处有个问题 前面计算clonetype默认存储名称为clone_id 
# 但是tl_network只有在 nt + iderntity的时候才会调用clone_id 其余时间调用cc_{sequence}_{metric} 
#改进为 前面函数调整过metric和sequence后 依然存储为clone_id  此处也写死调用clone_id 

####################################################################绘制克隆型分布图

ir.pl.clonotype_network(
  adata, 
  color="orig.ident",
  label_fontsize=label_fontsize1, 
  panel_size=(12, 10),
  label_alpha = label_alpha1,
  base_size = base_size1,
  edges_width = edges_width1
  )

aa =  dcgraph + "define_clonotypes_pl_clonotype_network_Sample.png"
aa1 = dcgraph + "define_clonotypes_pl_clonotype_network_Sample.pdf"   
plt.savefig(aa,bbox_inches='tight')
plt.savefig(aa1,bbox_inches='tight')

##########################################################################鉴定clonetype_cluster(以下简称cc)

###############################################宽松鉴定的默认参数

# metric2 = 'alignment'
# sequence2 = 'aa'
# #注意 这里布局和sequence写死 原因在于首先alignment和aa配套 另外后续有需要用特定obs列名获取cc信息的 
# cutoff2 = 15
# #define_clonotypes参数
# receptor_arms2 = 'all'
# dual_ir2 = "any"
# within_group2 = 'receptor_type'
# min_cells2 = 3
# min_nodes2 = 1
# base_size2 = None
# size_power2 = 0.9
# label_fontsize2 = 9
# label_alpha2 = 0.6
# edges_width2 = 0.5

###############################################计算距离矩阵

ir.pp.ir_dist(adata,metric=metric2,sequence=sequence2,cutoff=15)

###############################################鉴定cc

#此处不需要声明调用obs里面哪一列 这里默认是ir_aa_alignment

ir.tl.define_clonotype_clusters(
  adata, 
  sequence=sequence2,
  metric=metric2, 
  receptor_arms=receptor_arms2, 
  dual_ir=dual_ir2,
  inplace = True,
  same_v_gene = False,
  within_group = within_group2
)

# 原理性内容：BCR受体包含2x轻链 由V+J基因编码/2x重链 由V+D+J编码 
# TCR受体有两套 95%是α和β 外周血中的T细胞（1-5%）和在粘膜部位发现的T细胞表达的γδ链（γδ TCR）组成  其中α和γ是V+J β和δ是V+D+J
# 所以TCR和BCR流程差不多 但是可能obs里信息会多一点 外加qc时候要做一定的改写

ir.tl.clonotype_network(
  adata, 
  min_cells=min_cells2, 
  sequence=sequence2, 
  metric=metric2,
  min_nodes = min_nodes2,
  size_power = size_power2,
  inplace = True
  )

#########################################################绘制clonetype_cluster布局图

ir.pl.clonotype_network(
  adata, 
  color="orig.ident",
  label_fontsize=label_fontsize2,
   panel_size=(12, 10),
   label_alpha = label_alpha2,
   base_size = base_size2,
   edges_width = edges_width2
   )

aa =  dcgraph + "define_clonotypes_cluster_pl_clonotype_network_Sample.png"
aa1 = dcgraph + "define_clonotypes_cluster_pl_clonotype_network_Sample.pdf"   
plt.savefig(aa,bbox_inches='tight')
plt.savefig(aa1,bbox_inches='tight')

#############################################################设定V基因相同后 观察cc变化

ir.tl.define_clonotype_clusters(
  adata, 
  sequence=sequence2,
  metric=metric2, 
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

a = dcgraph + "clonotypes_cluster_difference.txt"
ct_different.to_csv(a,sep='\t',index = False)

############################################################写出adatair
#此处由于uns里还是记录了距离矩阵 所以不能用h5保存出来

obs = adata.obs

obs.index.name = 'cellnames'
aa = path1 + "/scirpy_ir_result.txt"
obs.to_csv(aa,sep = "\t")

#找到一个新模块 叫pickle 专门保存奇怪的机器学习模型 刚好可以存这个 特定的规则 保存再读取不会破坏anndata结构

aa = path1 + "/adatair.dat"
pickle.dump(adata,open(aa,"wb"))
# loaded_model = pickle.load(open(aa,"rb"))


#docker rmi -f 192.168.0.171:5000/bioscirpy.root:20221212

# docker pull 192.168.0.171:5000/bioscirpy.root:20221212
# docker pull 192.168.0.171:5000/bioscirpy:20221212









