library(sceasy)
library(getopt)
library(dplyr)
library(Seurat)

#202301011 重新设计V2版本  
# 1.此处首先处理带有Cluster信息的rds 样本信息以rds记录的为准  
# 2.根据输入的summarycell文件 筛选数据 重新分配celltype  不输入就用数字cluster
# 3.根据参数 挑选一部分cluster/celltype做分析 

options(stringsAsFactors = FALSE)
options(future.globals.maxSize= 30*1024*1024*1024)
command = matrix(c('inputrds','a',1,"character",#输入的跑过cluster的rds 后续转成adata到python中 
                   'summarycell','b',2,"character", #带有cluster/celltype信息的表 如果rds里没有seurat_cluster信息 且不传这个表 直接报错
                   'path','c',1,"character",
                   'selectcluster','d',2,"character")#中间文件夹 暂时定为tmp   
                 ,byrow=TRUE, ncol=4)

args=getopt(command)

inputrds1 = args$inputrds;#输入文件 norm过得rds
summarycell1 = args$summarycell;
path1 = args$path;
selectcluster = args$selectcluster
#path1 = "/media/nbc1/lijialun/scirpy/xuna/zhale"
setwd(path1)

# inputrds1 = "/media/nbc1/lijialun/scirpy/xuna/zhale/Manual.seuset.rds"
# summarycell1 = "/media/nbc1/lijialun/scirpy/xuna/zhale/B_cells_type_OnlyMBC.txt"
# selectcluster = "MBC"
seuset<-readRDS(inputrds1)

# 有可能出现 ir里的样本名和rds里的样本名不同 那么以rds里的样本名为准 
# 所以需要在rds里把orig.ident改成一个其他名字
seuset$sample<-seuset$orig.ident
#20230214发现 seurat——clusters只会存数字类型的cluster 所以预先转存一下
seuset$seurat_clusters<-seuset@active.ident

if(!is.null(summarycell1)){
        summarycell = read.table(summarycell1,sep = "\t",header = T,check.names = F)
        needcell = intersect(summarycell[,1],colnames(seuset))
        rownames(summarycell)<-summarycell[,1]
        seuset<-subset(seuset,cell = needcell)
        summarycell<-summarycell[needcell,]
        summarycell = summarycell[,2]
        names(summarycell)<-needcell
        summarycell = as.factor(summarycell)
        seuset$seurat_clusters = summarycell
}

if(!("seurat_clusters" %in% names(seuset@meta.data))){
    stop("no cluster info,check rds/summarycell")
    }

if(!is.null(selectcluster)){
    selectcluster = unlist(strsplit(selectcluster, "[,]"))
	cluster = unique(seuset$seurat_clusters)
	select = cluster %in% selectcluster
	if(!(TRUE %in% select)){
		stop("Please check your selected clusters")
	}
    if(TRUE %in% select){
        selectcluster = cluster[select]
        sum = seuset$seurat_clusters
        index = sum %in% selectcluster
        cellname = names(sum[index])
        seuset<-subset(seuset,cell = cellname)
    }
}

#这里ident也要改 不然这个seurat_clusters的level里还是包括全部的celltype
#20230214发现 还是有这个问题 如果把active.ident赋值到seurat_clusters中 那么本身就是factor格式 
#需要先转为character再转到factor重置一下
seuset$seurat_clusters<-as.factor(as.character(seuset$seurat_clusters))

#但是有个问题 不知道哪里处理出问题了 没有转存seurat_clusters出来
# 问题找到了 只有一个cluster的根本不允许转存 
#Dropping single category variables:seurat_clusters
#drop_single_values = F 这个参数防止丢掉只有一个cluster的情况 
sceasy::convertFormat(seuset, from="seurat", to="anndata",outFile='seuratsub.h5ad',drop_single_values = F)





