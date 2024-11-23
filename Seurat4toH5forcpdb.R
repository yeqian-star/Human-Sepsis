#利用SeuratDisk包将SeuratObject保存为h5ad文件,并处理Celltype信息
library(Seurat)
library(SeuratDisk)
library(getopt)

command = matrix(c('rdsfile','r',1,"character",
        'outpath','o',1,"character",
		'celltypefile','c',2,"character",
		'dbgenefile','d',1,"character")
		,byrow=TRUE, ncol=4)
args = getopt(command)
options(stringsAsFactors = FALSE)
rdsfile = args$rdsfile;
outpath = args$outpath;
celltypefile = args$celltypefile;
dbgenefile = args$dbgenefile;
#prefix = args$prefix;

setwd(outpath)
if(endsWith(rdsfile,"rds")){
	seuset <- readRDS(rdsfile)
	#不将scaledata赋为空的话，annadata会将h5ad读成scaledata的结果
	seuset@assays$RNA@scale.data = new(Class="matrix")

if(!is.null(celltypefile))
{
	newcluster <- read.table(celltypefile,sep="\t",header=F)
	if(TRUE %in% (newcluster[,1] %in% colnames(seuset))){
		realcluster = newcluster[is.finite(match(newcluster[,1],colnames(seuset))),]
		seuset = subset(seuset,cells=as.character(realcluster[,1]))	
		seuset = SetIdent(object = seuset, cells= as.character(realcluster[,1]), as.character(realcluster[,2]))
	}else if(TRUE %in% (newcluster[,1] %in% Idents(seuset))){
                realcluster = newcluster[is.finite(match(newcluster[,1],Idents(seuset))),]
                cellIn = Idents(seuset) %in% realcluster[,1]
                seuset = seuset[,cellIn]
                cluster = as.character(realcluster[,2])
                names(cluster) = realcluster[,1]
                cluster = cluster[as.character(Idents(seuset))]
                seuset = SetIdent(object = seuset, cells= colnames(seuset), cluster)
        }else{
                stop("Please check your celltypeFile")
        }

}

seuset$seurat_clusters = seuset@active.ident
seuset@meta.data$orig.ident <- as.character(seuset@meta.data$orig.ident)
seuset@meta.data$seurat_clusters <- as.character(seuset@meta.data$seurat_clusters)

SaveH5Seurat(seuset, filename = paste0("SeuratForCPdb.h5Seurat"))
Convert(paste0("SeuratForCPdb.h5Seurat"),dest="h5ad")
}
#将celltype第二列中数字部分转成字符串，防止cellphone自动读入时出现类型报错
newcluster <- read.table(celltypefile,sep="\t",header=F)
if("0"%in%(newcluster[,2])){
	newcluster[,2] <- paste0("Cluster",newcluster[,2])
}
write.table(newcluster,paste0('CellClusterListFile.txt'),sep="\t",quote=F,row.names=F,col.names=F)
if(endsWith(rdsfile,"bgz")){
	counts <- read.table(paste0(outpath,"counts.txt"),sep="\t",header=T)
	title <- read.table(paste0(outpath,"counts.txt"),nrows=1,header=F)
	dbgene <- read.table(dbgenefile,sep=",",header=T)
	targetgene <- intersect(dbgene[,4],counts[,1])
	rownames(counts)<- counts[,1]
	cutcounts <- counts[targetgene,]
	colnames(cutcounts)<-title
	write.table(cutcounts,paste0(outpath,"cutcounts.txt"),sep="\t",row.names=F,col.names=T,quote=F)
}
