library(sceasy)
library(getopt)
library(dplyr)
library(Seurat)

#脚本目的：判断输入文件是否带有cluster信息+载入细胞分类表筛选数据+选定某些cluster (此处如果输入细胞分类表,选定cluster的参数以表中为准)

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
#path = "/media/nbc1/lijialun/scirpy/TCRCluster"
setwd(path1)

# inputrds1 = "/media/nbc1/lijialun/scirpy/inputfile/Allsample_GraphClust.seuset.rds"
# summarycell1 = "/media/nbc1/lijialun/scirpy/inputfile/xifensummrycell.txt"
# seurat_clusters = 
seuset<-readRDS(inputrds1)

# 有可能出现 ir里的样本名和rds里的样本名不同 那么以rds里的样本名为准 
# 所以需要在rds里把orig.ident改成一个其他名字
seuset$sample<-seuset$orig.ident

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
seuset$seurat_clusters<-as.factor(seuset$seurat_clusters)

sceasy::convertFormat(seuset, from="seurat", to="anndata",outFile='seuratsub.h5ad')









