library(monocle)
library(Seurat)
library(getopt)
library(plyr)
library(reshape2)
library(dplyr)
library(magrittr)
library(tibble)
args <- commandArgs(T);

command = matrix(c('rds','r',1,"character",
    'outpath','o',1,"character",
    'clusters','c',2,"character",
	'celllist','C',2,"character",
    'gene','g',2,"character",
	'randomType','T',2,"character",
	'ReductionType','t',2,"character",
	'CellsPercent','P',2,"numeric",
	'MinCellsPerCluster','M',2,"numeric",
	'RandomCells','R',0,"logical",
    'method','m',1,"character")
    ,byrow=TRUE, ncol=4)
args=getopt(command)

rdsFile = args$rds;
outpath = args$outpath;
geneFile = args$gene;
reductionmethod = args$method;
cluster = args$clusters;
celllist = args$celllist;
if(is.null(args$CellsPercent)){args$CellsPercent=0}
if(is.null(args$MinCellsPerCluster)){args$MinCellsPerCluster=0}
if(is.null(args$RandomCells)){args$RandomCells=FALSE}
if(is.null(args$randomType)){args$randomType="random"}
if(is.null(args$ReductionType)){args$ReductionType="tsne"}

ReductionType=args$ReductionType
randomType = args$randomType
cellsPercent <<- args$CellsPercent
minCellsPerCluster <<- args$MinCellsPerCluster

if(cellsPercent==0&minCellsPerCluster==0){args$RandomCells=FALSE}
RandomCells = args$RandomCells


ImportSeurat4CDS = function (otherCDS){
        requireNamespace("Seurat")
	if("SCT" %in% names(seuset@assays)){	
        data <- GetAssayData(otherCDS, slot = "counts", assay="SCT")
	}else{
		data <- GetAssayData(otherCDS, slot = "counts", assay="RNA")
	}
    if (class(data) == "data.frame"){
        data <- as(as.matrix(data), "sparseMatrix")
    }
    pd <- tryCatch({
        pd <- new("AnnotatedDataFrame", data = otherCDS@meta.data)
        pd
    }, error = function(e) {
        pData <- data.frame(cell_id = colnames(data), row.names = colnames(data))
        pd <- new("AnnotatedDataFrame", data = pData)
        message("This Seurat object doesn't provide any meta data")
        pd
    })
    if (length(setdiff(colnames(data), rownames(pd))) > 0) {
        data <- data[, rownames(pd)]
    }
    fData <- data.frame(gene_short_name = row.names(data),
    row.names = row.names(data))
    fd <- new("AnnotatedDataFrame", data = fData)
    lowerDetectionLimit <- 0
    if (all(data == floor(data))) {
        expressionFamily <- negbinomial.size()
    }else if(any(data < 0)){
        expressionFamily <- uninormal()
    }else{
        expressionFamily <- tobit()
    }
    valid_data <- data[, row.names(pd)]
    monocle_cds <- newCellDataSet(data, phenoData = pd, featureData = fd,
    lowerDetectionLimit = lowerDetectionLimit, expressionFamily = expressionFamily)
    if ("Monocle" %in% names(otherCDS@misc)) {
        otherCDS@misc$Monocle@auxClusteringData$seurat <- NULL
        otherCDS@misc$Monocle@auxClusteringData$scran <- NULL
        monocle_cds <- otherCDS@misc$Monocle
        mist_list <- otherCDS
    }else{
        mist_list <- otherCDS
    }

    if ("var.genes" %in% slotNames(otherCDS)) {
        var.genes <- setOrderingFilter(monocle_cds, otherCDS@var.genes)
    }
    monocle_cds@auxClusteringData$seurat <- mist_list
	return(monocle_cds)
}


setNum = function(number){
	selectNum = round(number*cellsPercent)
    if(selectNum>=minCellsPerCluster){
        finalNum = selectNum
    }else{
        finalNum = min(number,minCellsPerCluster)
    }
	return(finalNum)
}

randomSample = function(cluster,size=NULL){
	samples = c()
	id = unique(cluster)
	for(i in id){
		clusterInfo = cluster[cluster==i]
		if(is.null(size)){
			finalNum = setNum(length(clusterInfo))
		}else{
			finalNum = size
		}
		insample = sample(names(clusterInfo),finalNum,replace = FALSE)
		samples = c(samples,insample)
	}
	return(samples)
}

setwd(outpath)
seuset = readRDS(rdsFile)

if(substring(seuset@version,1,1)<=3){
	stop("Please check out the version of rdsFile.")
}#检查seurat版本

if(!is.null(celllist)){
	newcluster = read.table(celllist,sep="\t",header=F)
	realcluster = newcluster[is.finite(match(newcluster[,1],colnames(seuset))),]
	seuset = subset(seuset,cells=as.character(realcluster[,1]))
        seuset = SetIdent(object = seuset, cells= as.character(realcluster[,1]), as.character(realcluster[,2]))

}else if(!is.null(cluster)){
	clusters = unlist(strsplit(cluster, "[,]"))
	celllist = names(seuset@active.ident[is.finite(match(seuset@active.ident,clusters))])
	seuset = subset(seuset,cells=celllist)
}
identType = unique(seuset@active.ident)

if(RandomCells){
	if(randomType=="random"){
		newcelllist = randomSample(seuset@active.ident)
	}
	if(randomType=="kmean"){
		newcelllist = c()
		for(i in identType){
			if(ReductionType=="tsne"){
				reductions = seuset@reductions$tsne@cell.embeddings[names(seuset@active.ident[seuset@active.ident==i]),]
			}else if(ReductionType=="umap"){
				reductions = seuset@reductions$umap@cell.embeddings[names(seuset@active.ident[seuset@active.ident==i]),]
			}
			finalNum = setNum(nrow(reductions))
			if(finalNum>=nrow(reductions)){
				newcelllist = c(newcelllist,rownames(reductions))
			}else{
				clusterid = kmeans(reductions, centers = finalNum)$clust
				newcelllist = c(newcelllist,randomSample(clusterid,1))
			}
		}
	}
	if(randomType=="hclust"){
		newcelllist = c()
		for(i in identType){
			if(ReductionType=="tsne"){
                reductions = seuset@reductions$tsne@cell.embeddings[names(seuset@active.ident[seuset@active.ident==i]),]
        	}else if(ReductionType=="umap"){
            	reductions = seuset@reductions$umap@cell.embeddings[names(seuset@active.ident[seuset@active.ident==i]),]
        	}
        	finalNum = setNum(nrow(reductions))
        	if(finalNum>=nrow(reductions)){
                newcelllist = c(newcelllist,rownames(reductions))
            }else{
				dist.r = dist(reductions)
				hc.r=hclust(dist.r)
				clusterid = cutree(hc.r,finalNum)
				newcelllist = c(newcelllist,randomSample(clusterid,1))
			}
		}
	}
	seuset = subset(seuset,cells=as.character(newcelllist))
	write.table(data.frame(Cell = colnames(seuset),Cluster=seuset@active.ident),file=paste0("RandomCells.txt"),sep="\t",row.names=F,quot=F)
}

seuset <- AddMetaData(object = seuset, metadata = seuset@active.ident, col.name = "FinalCluster")
seuCDSAll = ImportSeurat4CDS(seuset)

if(!is.null(geneFile)){
	genelist = unique(read.table(geneFile,sep="\t",header=T)[,1])
}else{
	genelist = VariableFeatures(seuset)
}

seuCDSAll <- setOrderingFilter(seuCDSAll, genelist)
seuCDSAll <- estimateSizeFactors(seuCDSAll)
seuCDS <- reduceDimension(seuCDSAll, max_components = 2, reduction_method = reductionmethod)
seuCDS <- orderCells(seuCDS)
table = table(seuCDS$FinalCluster,seuCDS$State)
print(table)
table = dcast(data.frame(table),Var1~Var2)
colnames(table)[1]="Cluster"
write.table(table,file=paste0("Pseudotime.Statistics.txt"),sep="\t",row.names=F,quot=F)
phenoData = data.frame(CellName=rownames(seuCDS@phenoData@data),seuCDS@phenoData@data)
write.table(phenoData,file=paste0("Pseudotime.Summary_Cell.txt"),sep="\t",row.names=F,quot=F)

seuCDS$FinalCluster = as.character(seuCDS$FinalCluster)
uniquecluster = unique(seuCDS$FinalCluster)
if(!NA %in% as.numeric(uniquecluster)){
	uniquecluster = as.numeric(uniquecluster)
}
sort = sort(uniquecluster)
seuCDS$Cluster = factor(seuCDS$FinalCluster,levels=sort)

gg =  plot_cell_trajectory(seuCDS, color_by = "State")+theme(axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
ggsave(paste0("Pseudotime.State.png"))
ggsave(paste0("Pseudotime.State.pdf"))
gg =  plot_cell_trajectory(seuCDS, color_by = "Cluster")+theme(axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
ggsave(paste0("Pseudotime.Cluster.png"))
ggsave(paste0("Pseudotime.Cluster.pdf"))
gg =  plot_cell_trajectory(seuCDS, color_by = "Pseudotime")+theme(axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
ggsave(paste0("Pseudotime.Pseudotime.png"))
ggsave(paste0("Pseudotime.Pseudotime.pdf"))
gg =  plot_cell_trajectory(seuCDS, color_by = "Cluster") + facet_wrap(~Cluster, nrow = floor(sqrt(length(unique(seuCDS$Cluster)))))+theme(axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
ggsave(paste0("Pseudotime.ClusterEach.png"))
ggsave(paste0("Pseudotime.ClusterEach.pdf"))
gg =  plot_cell_trajectory(seuCDS, color_by = "State") + facet_wrap(~State, nrow = floor(sqrt(length(unique(seuCDS$State)))))+theme(axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
ggsave(paste0("Pseudotime.StateEach.png"))
ggsave(paste0("Pseudotime.StateEach.pdf"))

gg = plot_complex_cell_trajectory(seuCDS, color_by = 'Cluster', show_branch_points = T, cell_size = 0.5, cell_link_size = 0.3)+theme(axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
ggsave(paste0("Pseudotime.tree.png"))
ggsave(paste0("Pseudotime.tree.pdf"))

if(length(unique(seuCDS$orig.ident))>1){
seuCDS$Sample = as.factor(seuCDS$orig.ident)
gg =  plot_cell_trajectory(seuCDS, color_by = "Sample")+theme(axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
ggsave(paste0("Pseudotime.Sample.png"))
ggsave(paste0("Pseudotime.Sample.pdf"))
gg =  plot_cell_trajectory(seuCDS, color_by = "Cluster") + facet_wrap(~Sample, nrow = floor(sqrt(length(unique(seuCDS$Sample)))))+theme(axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
ggsave(paste0("Pseudotime.SampleEach.png"))
ggsave(paste0("Pseudotime.SampleEach.pdf"))
gg = plot_complex_cell_trajectory(seuCDS, color_by = 'Cluster', show_branch_points = T, cell_size = 0.5, cell_link_size = 0.3) + facet_wrap(~Sample, nrow = floor(sqrt(length(unique(seuCDS$Sample))))) + scale_size(range = c(0.2, 0.2))+theme(axis.text.x = element_text(angle = 30, hjust = 1))+theme(axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))
ggsave(paste0("SampleEach.tree.png"))
ggsave(paste0("SampleEach.tree.pdf"))
}

dir.create(paste0("../Pseudotime.coords"))

sample_state <- pData(seuCDS)$State;
lib_info_with_pseudo <- pData(seuCDS);
if (seuCDS@dim_reduce_type == "ICA") {
    reduced_dim_coords <- reducedDimS(seuCDS)
}else if(seuCDS@dim_reduce_type %in% c("simplePPT", "DDRTree")) {
    reduced_dim_coords <- reducedDimK(seuCDS)
}
x = 1;
y = 2;
ica_space_df <- Matrix::t(reduced_dim_coords) %>% as.data.frame() %>%
        select_(prin_graph_dim_1 = x, prin_graph_dim_2 = y) %>%
        mutate(sample_name = rownames(.), sample_state = rownames(.))
dp_mst <- minSpanningTree(seuCDS)
edge_df <- dp_mst %>% igraph::as_data_frame() %>% select_(source = "from",
        target = "to") %>% left_join(ica_space_df %>% select_(source = "sample_name",
        source_prin_graph_dim_1 = "prin_graph_dim_1", source_prin_graph_dim_2 = "prin_graph_dim_2"),
        by = "source") %>% left_join(ica_space_df %>% select_(target = "sample_name",
        target_prin_graph_dim_1 = "prin_graph_dim_1", target_prin_graph_dim_2 = "prin_graph_dim_2"),
        by = "target")
data_df <- t(monocle::reducedDimS(seuCDS)) %>% as.data.frame() %>%
        select_(Component_1 = x, Component_2 = y) %>% rownames_to_column("CellName")

mst_branch_nodes <- seuCDS@auxOrderingData[[seuCDS@dim_reduce_type]]$branch_points
if (seuCDS@dim_reduce_type == "DDRTree") {
branch_point_df <- ica_space_df %>% slice(match(mst_branch_nodes,
            sample_name)) %>% mutate(branch_point_idx = seq_len(n()))
write.table(branch_point_df[,c(5,1,2)],paste0("../Pseudotime.coords/label.2d.txt"),quote=F,row.names=F,sep="\t")
}
write.table(edge_df,paste0("../Pseudotime.coords/line.2d.txt"),quote=F,row.names=F,sep="\t")
write.table(data_df,paste0("../Pseudotime.coords/coord.2d.txt"),quote=F,row.names=F,sep="\t")
saveRDS(seuCDS, file = paste0("../Pseudotime.monocle.rds"))

seuCDS <- reduceDimension(seuCDSAll, max_components = 3, reduction_method = reductionmethod)
seuCDS <- orderCells(seuCDS)

sample_state <- pData(seuCDS)$State;
lib_info_with_pseudo <- pData(seuCDS);
if (seuCDS@dim_reduce_type == "ICA") {
        reduced_dim_coords <- reducedDimS(seuCDS)
    }else if (seuCDS@dim_reduce_type %in% c("simplePPT", "DDRTree")) {
        reduced_dim_coords <- reducedDimK(seuCDS)
    }
x = 1;
y = 2;
z = 3;
ica_space_df <- Matrix::t(reduced_dim_coords) %>% as.data.frame() %>%
        select_(prin_graph_dim_1 = x, prin_graph_dim_2 = y, prin_graph_dim_3 = z) %>%
        mutate(sample_name = rownames(.), sample_state = rownames(.))
dp_mst <- minSpanningTree(seuCDS)
edge_df <- dp_mst %>% igraph::as_data_frame() %>% select_(source = "from",
        target = "to") %>% left_join(ica_space_df %>% select_(source = "sample_name",
        source_prin_graph_dim_1 = "prin_graph_dim_1", source_prin_graph_dim_2 = "prin_graph_dim_2", source_prin_graph_dim_3 = "prin_graph_dim_3"),
        by = "source") %>% left_join(ica_space_df %>% select_(target = "sample_name",
        target_prin_graph_dim_1 = "prin_graph_dim_1", target_prin_graph_dim_2 = "prin_graph_dim_2", target_prin_graph_dim_3 = "prin_graph_dim_3"),
        by = "target")
data_df <- t(monocle::reducedDimS(seuCDS)) %>% as.data.frame() %>%
        select_(Component_1 = x, Component_2 = y, Component_3 = z) %>% rownames_to_column("CellName")

mst_branch_nodes <- seuCDS@auxOrderingData[[seuCDS@dim_reduce_type]]$branch_points
if (seuCDS@dim_reduce_type == "DDRTree") {
branch_point_df <- ica_space_df %>% slice(match(mst_branch_nodes,
            sample_name)) %>% mutate(branch_point_idx = seq_len(n()))
write.table(branch_point_df[,c(6,1,2,3)],paste0("../Pseudotime.coords/label.3d.txt"),quote=F,row.names=F,sep="\t")
}
write.table(edge_df,paste0("../Pseudotime.coords/line.3d.txt"),quote=F,row.names=F,sep="\t")
write.table(data_df,paste0("../Pseudotime.coords/coord.3d.txt"),quote=F,row.names=F,sep="\t")



