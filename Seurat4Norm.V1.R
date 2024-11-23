
library(SingleCellExperiment)
library(scater)
library(plyr)
library(reshape2)
library(Seurat)
library(dplyr)
library(getopt)
library(stringr)
library(cowplot)
library(ggplot2)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize= 30*1024*1024*1024)
command = matrix(c('counts','r',1,"character",
                'prefix','s',1,"character",
		        'celllist','L',2,"character",
                'outpath','o',1,"character",
                'mingene','g',2,"numeric",
                'mincell','c',2,"numeric",
                'maxgene','G',2,"numeric",
                'maxMT','m',2,"numeric",
                'maxPT','b',2,"numeric",
				'HouseKeepingPercent','K',2,"numeric",
				'maxDAPercent','a',2,"numeric",
                'removeGene','M',2,"character",
                'scalemodule','S',1,"character",
				'ExtraAnalysis','A',2,"character",
				'NormType','B',1,"character",
                'regressOut','R',2,"character",
				'databasepath','e',2,"character",
                'tmppath','T',1,"character",
				'parameter','D',1,"character",
                'runFast','F',0,"logical",
				'writeMatrix','W',0,"logical",
				'samplelist','E',2,"character",
				'pcadim','P',2,"character",
				'nneighbors','N',2,"numeric",
				'computePCs','p',2,"numeric",
				'AsNorm','n',0,"logical",
				'varGeneNum','V',2,"numeric")
                ,byrow=TRUE, ncol=4)

args=getopt(command)
if(is.null(args$removeGene)){args$removeGene = ""}
if(is.null(args$ExtraAnalysis)){args$ExtraAnalysis = ""}
countsFile = args$counts;
prefix = args$prefix;
outpath = args$outpath;
geneExpMinCell = args$mincell;
cellMinGene = args$mingene;
cellMaxGene = args$maxgene;
MTlist = args$MTList;
MaxMTPercent = args$maxMT;
MaxPTPercent = args$maxPT;
HouseKeepingPercent = args$HouseKeepingPercent;
maxDAPercent = args$maxDAPercent;
ScaleModel = args$scalemodule;
removeGene = unlist(strsplit(args$removeGene,"[,]"))
ExtraAnalysis = unlist(strsplit(args$ExtraAnalysis,"[,]"))
RemoveMTGene = args$removeMT;
tmppath = args$tmppath;
if(is.null(args$runFast)){args$runFast = FALSE}
if(is.null(args$AsNorm)){args$AsNorm = FALSE}
if(is.null(args$writeMatrix)){args$writeMatrix = FALSE}
RunFast = args$runFast;
pcadim = args$pcadim
writeMatrix = args$writeMatrix
computePCs = args$computePCs
varGeneNum = args$varGeneNum
nneighbors = args$nneighbors
NormType = args$NormType
regressOut = args$regressOut
parameter = args$parameter
databasepath = args$databasepath
celllist = args$celllist

getPCDim = function(x){
        pcDim = c()
        x = unlist(strsplit(x,"[,]"))
        for(i in x){
                info = unlist(strsplit(i,"[:]"))
                if(length(info)!=1){
                        info = seq(info[1],info[2])
                }
                pcDim = c(pcDim,info)
        }
	pcDim = as.matrix(pcDim)
        print(paste("usePCs:",paste(pcDim,collapse=",")))
        return(pcDim)
}

metadataplot = function (object, features, dims = c(1, 2),pt.size = 1,cols=c("lightgrey", "blue"),reduction = "tsne"){
        dims <- paste0(Key(object = object[[reduction]]), dims)
        data <- FetchData(object = object, vars = c(dims,features))
        gg = ggplot(data=data,mapping = aes_string(x = dims[1],y=dims[2],color=features))+geom_point(size=pt.size)+scale_color_gradientn(colours = cols)+labs(color = NULL)+ theme_cowplot()+labs(title = features)+ theme(plot.title = element_text(hjust = 0.5,size = 30))+theme(axis.title.x=element_text(size=25),axis.title.y=element_text(size=25),axis.text.x=element_text(size=15),axis.text.y=element_text(size=15))

        return(gg)
}

#默认参数设置
if(parameter=="Default"){
	if(NormType=="scale"){
		pcadim = "1:10"
		computePCs = 50
		nneighbors = 30
		varGeneNum = 2000		
	}else if(NormType=="sctransform"){
		pcadim = "1:30"
                computePCs = 50
                nneighbors = 30
                varGeneNum = 3000
	}
}

print("Start")
setwd(outpath)

if(dir.exists(countsFile)){
        allcounts = Read10X(data.dir=countsFile)
        if(class(allcounts)=="list"){
        	allcounts = allcounts[[1]]
        }
}else{
		if(endsWith(countsFile,"csv")||endsWith(countsFile,"csv.gz")){
            allcounts = read.csv(countsFile,header=T,row.names=1)
        	title = read.csv(countsFile, he=F, nrow = 1,colClasses='character',row.names=1)
	        colnames(allcounts) = title
        }else if(endsWith(countsFile,"h5")){
            allcounts = Read10X_h5(countsFile)
			if(class(allcounts)=="list"){
         		allcounts = allcounts[[1]]
	        }
        }else{
            allcounts = read.table(countsFile,sep="\t",header=T,row.names=1)
			title = read.table(countsFile, he=F, nrow = 1,colClasses='character',sep="\t",row.names=1)
	        colnames(allcounts) = title
        }
}

realOrigIdents = NULL
if(!is.null(celllist)){
    celllist = read.table(args$celllist,sep="\t",header=F,comment.char = "")
	rownames(celllist) = celllist[,1]
	realcells = intersect(colnames(allcounts),celllist[,1])
	allcounts = allcounts[,realcells]
	if(ncol(celllist)>1){
		realcells = celllist[realcells,]
		realOrigIdents = realcells[,2]
		names(realOrigIdents) = realcells[,1]
    }
}
dir.create("1.Normalization")
dir.create("2.PCAAnalysis")
dir.create("3.TSNEAnalysis")
dir.create("4.UMAPAnalysis")
dir.create(paste0(tmppath,"/",prefix),recursive = TRUE)


allcounts = allcounts[,colSums(allcounts)>0,drop=F]

if("MT" %in% ExtraAnalysis){
	mito.genes = read.table(paste0(databasepath,"/MT.txt"), sep = "\t", header = FALSE)[,1]
	mito.genes = mito.genes[is.finite(match(mito.genes,rownames(allcounts)))]
	percent.mito <- Matrix::colSums(allcounts[mito.genes, ])/Matrix::colSums(allcounts)
	print("MT Gene Percent.")
	print(quantile(percent.mito,probs = seq(0, 1, 0.05)))
}

if("PT" %in% ExtraAnalysis){
	pt.genes = read.table(paste0(databasepath,"/PT.txt"), sep = "\t", header = FALSE)[,1]
	pt.genes = pt.genes[is.finite(match(pt.genes,rownames(allcounts)))]
	percent.pt <- Matrix::colSums(allcounts[pt.genes, ])/Matrix::colSums(allcounts)
	print("PT Gene Percent.")
	print(quantile(percent.pt,probs = seq(0, 1, 0.05)))
}

if("HouseKeeping" %in% ExtraAnalysis){
    HouseKeepinggenes = read.table(paste0(databasepath,"/HouseKeeping.txt"), sep = "\t", header = FALSE)[,1]
	num = length(HouseKeepinggenes)
    HouseKeepinggenes = HouseKeepinggenes[is.finite(match(HouseKeepinggenes,rownames(allcounts)))]
    percent.HouseKeeping = apply(allcounts[HouseKeepinggenes, ],2,FUN=function(x){return(sum(x > 0) / num)})
	print("HouseKeeping Gene Percent.")
	print(quantile(percent.HouseKeeping,probs = seq(0, 1, 0.05)))
	forplot = data.frame(percent.HouseKeeping)
	gg = ggplot(data=forplot,mapping=aes(x=percent.HouseKeeping))+geom_density()+theme_cowplot()
	ggsave("1.Normalization/HouseKeepingGenePercent.pdf",gg,width = 8,height = 6)
	ggsave("1.Normalization/HouseKeepingGenePercent.png",gg,width = 8,height = 6)
}

if("dagene" %in% ExtraAnalysis){
        dagenes = read.table(paste0(databasepath,"/Dissociation-Associated.txt"), sep = "\t", header = FALSE)[,1]
        dagenes = dagenes[is.finite(match(dagenes,rownames(allcounts)))]
        percent.dagene <- Matrix::colSums(allcounts[dagenes, ])/Matrix::colSums(allcounts)
	print("Dissociation Associated Gene Percent.")
        print(quantile(percent.dagene,probs = seq(0, 1, 0.05)))
}
#构建seuratobject
seuset <- CreateSeuratObject(counts = allcounts,project = prefix,min.cells = geneExpMinCell, min.features = 1)
if("MTGene" %in% removeGene){
    print("Remove MTGene")
	seuset = subset(seuset,features = setdiff(row.names(seuset),mito.genes))
}
if("PTGene" %in% removeGene){
    print("Remove PTGene")
	seuset = subset(seuset,features = setdiff(row.names(seuset),pt.genes))
}
if("HouseKeeping" %in% removeGene){
    print("Remove HouseKeeping Gene")
	seuset = subset(seuset,features = setdiff(row.names(seuset),HouseKeepinggenes))
}
if("DAgene" %in% removeGene){
    print("Remove DissociationAssociated Gene")
	seuset = subset(seuset,features = setdiff(row.names(seuset),dagenes))
}
print("nFeature_RNA")
quantile(seuset@meta.data$nFeature_RNA,probs = seq(0, 1, 0.05))
if("MT" %in% ExtraAnalysis){
	seuset <- AddMetaData(object = seuset, metadata = percent.mito, col.name = "percent.mito")
	mtInfo = seuset@meta.data[,c("nCount_RNA","percent.mito")]
	mtInfo$percent.mito = round(100*mtInfo$percent.mito)/100
	cellNum = aggregate(mtInfo$nCount_RNA,by=data.frame(mtInfo$percent.mito),FUN=length)
	rownames(cellNum) = cellNum[,1]
	filter = cellNum[,2]>=(0.005*nrow(mtInfo))
	filter = rownames(cellNum)[filter]

	result = aggregate(mtInfo$nCount_RNA,by=data.frame(mtInfo$percent.mito),FUN=function(x){quantile(x)["75%"]})
	rownames(result) = result[,1]
	result = result[filter,,drop=F]

	index = order(result$x,decreasing = T)
	colnames(result) = c("percent.mito","UpperQuantileValue")
	percent = result[index[1],,drop=F][1,1]

	result = result[index,,drop=F]
	print(result)
	if(MaxMTPercent<percent){
 	    stop("Please check mito.percent")
	}else if(MaxMTPercent<percent+0.05){
        returnValue = 223
	}else{
	    returnValue = 0
	}
}else{
	seuset <- AddMetaData(object = seuset, metadata = 0, col.name = "percent.mito")
}
if("PT" %in% ExtraAnalysis){
	seuset <- AddMetaData(object = seuset, metadata = percent.pt, col.name = "percent.pt")
	ptInfo = seuset@meta.data[,c("nCount_RNA","percent.pt")]
	ptInfo$percent.pt = round(100*ptInfo$percent.pt)/100
	cellNum = aggregate(ptInfo$nCount_RNA,by=data.frame(ptInfo$percent.pt),FUN=length)
	rownames(cellNum) = cellNum[,1]
	filter = cellNum[,2]>=(0.005*nrow(ptInfo))
	filter = rownames(cellNum)[filter]

	result = aggregate(ptInfo$nCount_RNA,by=data.frame(ptInfo$percent.pt),FUN=function(x){quantile(x)["75%"]})
	rownames(result) = result[,1]
	result = result[filter,,drop=F]

	index = order(result$x,decreasing = T)
	colnames(result) = c("percent.pt","UpperQuantileValue")
	percent = result[index[1],,drop=F][1,1]

	result = result[index,,drop=F]
	print(result)
	if(MaxPTPercent<percent){
 	    stop("Please check pt.percent")
	}else if(MaxPTPercent<percent+0.05){
        returnValue = 223
	}else{
	    returnValue = 0
	}
}else{
	seuset <- AddMetaData(object = seuset, metadata = 0, col.name = "percent.pt")
}

if(!is.null(args$samplelist)){
    samplelist = args$samplelist;
    sampleInfo = read.table(samplelist, sep = "\t", header = FALSE,row.names=1)
    samples = sapply(strsplit(colnames(seuset), "-"), function(v) return(v[2]))
    if(length(unique(samples))==1){
        samples[is.na(samples)]<-sampleInfo[1,1]
        sampleInfo = samples
    }else{
        sampleInfo = sampleInfo[samples,]
    }
    names(sampleInfo) = colnames(seuset)
    seuset <- AddMetaData(object = seuset, metadata = sampleInfo, col.name = "orig.ident")
}

if(!is.null(realOrigIdents)){
    seuset <- AddMetaData(object = seuset, metadata = realOrigIdents, col.name = "orig.ident")
}

nSample = length(unique(seuset@meta.data$orig.ident))

if("HouseKeeping" %in% ExtraAnalysis){
    seuset <- AddMetaData(object = seuset, metadata = percent.HouseKeeping, col.name = "percent.HouseKeeping")
	gg = VlnPlot(object = seuset, features = "percent.HouseKeeping", group.by="orig.ident", pt.size = 0.5,ncol = 1) + NoLegend()
	ggsave("1.Normalization/percent.HouseKeeping.png",gg,width = 3+nSample/3,height = 7)
	ggsave("1.Normalization/percent.HouseKeeping.pdf",gg,width = 3+nSample/3,height = 7)
}
if("dagene" %in% ExtraAnalysis){
    seuset <- AddMetaData(object = seuset, metadata = percent.dagene, col.name = "percent.DissociationAssociated")
	gg = VlnPlot(object = seuset, features = "percent.DissociationAssociated", group.by="orig.ident", pt.size = 0.5, ncol = 1) + NoLegend()
        ggsave("1.Normalization/percent.DissociationAssociated.png",gg,width = 3+nSample/3,height = 7)
        ggsave("1.Normalization/percent.DissociationAssociated.pdf",gg,width = 3+nSample/3,height = 7)
}
write.table(data.frame(CellName=rownames(seuset@meta.data),seuset@meta.data),file="1.Normalization/CellInfosRaw.txt",sep="\t",quote=F,row.names=F)
write.table(data.frame(CellName=rownames(seuset@meta.data),seuset@meta.data),file="CellInfosRaw.txt",sep="\t",quote=F,row.names=F)

seuset = SetIdent(object = seuset, cells = colnames(seuset), value = seuset@meta.data$orig.ident)
gg = VlnPlot(object = seuset, features = c("nFeature_RNA", "nCount_RNA"), group.by="orig.ident", pt.size = 0.5, ncol = 2) + NoLegend()
ggsave("1.Normalization/Counts_Distribution.png",gg,width = 6+nSample/2,height = 7)
ggsave("1.Normalization/Counts_Distribution.pdf",gg,width = 6+nSample/2,height = 7)

gg = VlnPlot(object = seuset, features = "percent.mito", group.by="orig.ident", ncol = 1,pt.size = 0.5) + NoLegend()
ggsave("1.Normalization/percent.mito.png",gg,width = 3+nSample/3,height = 7)
ggsave("1.Normalization/percent.mito.pdf",gg,width = 3+nSample/3,height = 7)

plot1 <- FeatureScatter(object = seuset, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(object = seuset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot3 <- FeatureScatter(object = seuset, feature1 = "nFeature_RNA", feature2 = "percent.mito")

gg = CombinePlots(plots = list(plot1, plot2),ncol=2)
ggsave("1.Normalization/Counts_correlation.png",gg,width = 15,height = 7)
ggsave("1.Normalization/Counts_correlation.pdf",gg,width = 15,height = 7)

if("PT" %in% ExtraAnalysis){
	gg = VlnPlot(object = seuset, features = "percent.pt", group.by="orig.ident", ncol = 1,pt.size = 0.5) + NoLegend()
	ggsave("1.Normalization/percent.pt.png",gg,width = 3+nSample/3,height = 7)
	ggsave("1.Normalization/percent.pt.pdf",gg,width = 3+nSample/3,height = 7)
}

table = table(seuset@meta.data$orig.ident)
print(table)

print(paste0("nGeneFilterRegion:",cellMinGene,"-",cellMaxGene))
seuset <- subset(x = seuset, subset = nFeature_RNA > cellMinGene & nFeature_RNA < cellMaxGene)
if("MT" %in% ExtraAnalysis){
	print(paste0("MTPercentFilterRegion:<",MaxMTPercent))
	seuset <- subset(x = seuset, subset = percent.mito < MaxMTPercent)
}
if("PT" %in% ExtraAnalysis){
	print(paste0("PTPercentFilterRegion:<",MaxPTPercent))
	seuset <- subset(x = seuset, subset = percent.pt < MaxPTPercent)
}
if("HouseKeeping" %in% ExtraAnalysis){
        print(paste0("HouseKeepingPercentFilterRegion:>",HouseKeepingPercent))
	seuset <- subset(x = seuset, subset = percent.HouseKeeping > HouseKeepingPercent)
}
if("dagene" %in% ExtraAnalysis){
        print(paste0("Dissociation-AssociatedPercentFilterRegion:<",maxDAPercent))
	seuset <- subset(x = seuset, subset = percent.DissociationAssociated < maxDAPercent)
}
table = table(seuset@meta.data$orig.ident)
print(table)

if(!args$AsNorm){
	seuset <- NormalizeData(object = seuset, normalization.method = "LogNormalize", scale.factor = 10000)
}
seuset <- FindVariableFeatures(seuset, selection.method = "vst", nfeatures = varGeneNum)

if("CellCycle" %in% ExtraAnalysis){
        ccRDS = paste0(databasepath,"/CCGenes.rds");
        ccgene = readRDS(ccRDS)
        seuset = CellCycleScoring(seuset, s.features =ccgene$s.genes, g2m.features = ccgene$g2m.genes,  set.ident = FALSE)
        seuset$CC.Difference <- seuset$S.Score - seuset$G2M.Score
}

if("dagene" %in% ExtraAnalysis){
	seuset = AddModuleScore(seuset,list(dagenes),name="DissociationAssociated.Score")
	colnames(seuset@meta.data)[colnames(seuset@meta.data)=="DissociationAssociated.Score1"]="DissociationAssociated.Score"
}

if(!is.null(regressOut)){
    regressOut = unlist(strsplit(regressOut, "[,]"))
    if(length(unique(seuset@meta.data$orig.ident))==1){
        regressOut = regressOut[!regressOut%in%"orig.ident"]
    }
    if(!"CellCycle" %in% ExtraAnalysis){
        regressOut = regressOut[!regressOut%in%c("S.Score","G2M.Score","CC.Difference")]
    }
	if(!"dagene" %in% ExtraAnalysis){
        regressOut = regressOut[!regressOut%in%c("percent.DissociationAssociated","DissociationAssociated.Score")]
    }
    if(NormType=="sctransform"){
        regressOut = regressOut[!regressOut%in%c("nFeature_RNA","nCount_RNA")]
    }
    print(regressOut)
}

if(NormType=="scale"){
	top10 <- head(VariableFeatures(seuset), 10)

	plot1 <- VariableFeaturePlot(seuset)
	gg <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
	ggsave("1.Normalization/VariableGenesPlot.pdf",gg,width = 10,height = 7)
	ggsave("1.Normalization/VariableGenesPlot.png",gg,width = 10,height = 7)

	#write.table(data.frame(Gene=rownames(seuset@hvg.info),seuset@hvg.info),file="1.Normalization/hvgInfo.txt",sep="\t",row.names=F,quote=F)
	write.table(VariableFeatures(seuset),file="1.Normalization/varGene.txt",sep="\t",row.names=F,quote=F)

	if(RunFast){
        	scalegene = VariableFeatures(seuset)
        	print(paste0("Use Variable genes to scale data"))
	}else{
        	scalegene = rownames(seuset)

	}
	seuset <- ScaleData(object = seuset, features = scalegene, vars.to.regress = regressOut, model.use = ScaleModel)
}

if(NormType=="sctransform"){
	seuset <- SCTransform(seuset, vars.to.regress = regressOut, variable.features.n = varGeneNum, return.only.var.genes=RunFast)
	top10 <- head(VariableFeatures(seuset), 10)

    plot1 <- VariableFeaturePlot(seuset)
    gg <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    ggsave("1.Normalization/VariableGenesPlot.pdf",gg,width = 10,height = 7)
    ggsave("1.Normalization/VariableGenesPlot.png",gg,width = 10,height = 7)
	write.table(VariableFeatures(seuset),file="1.Normalization/varGene.txt",sep="\t",row.names=F,quote=F)
}

write.table(data.frame(CellName=rownames(seuset@meta.data),seuset@meta.data),file="1.Normalization/CellInfosFilter.txt",sep="\t",quote=F,row.names=F)
write.table(data.frame(CellName=rownames(seuset@meta.data),seuset@meta.data),file="CellInfosFilter.txt",sep="\t",quote=F,row.names=F)

pcadim = getPCDim(pcadim)

seuset <- RunPCA(object = seuset, features = VariableFeatures(seuset), npcs = computePCs, do.print = TRUE, ndims.print = 1:computePCs, nfeatures.print = 10)

gg = VizDimLoadings(seuset, dims = pcadim, reduction = "pca",balanced = TRUE,nfeatures = 20)
ggsave("2.PCAAnalysis/VizPCAPlot.pdf",gg,width = 10,height = 10)
ggsave("2.PCAAnalysis/VizPCAPlot.png",gg,width = 10,height = 10)

gg = DimHeatmap(seuset, dims = 1:6, cells = min(500,length(colnames(seuset))),nfeatures = 20, balanced = TRUE, fast = FALSE)
ggsave("2.PCAAnalysis/PCHeatmap.pdf",gg,width = 10,height = 10)
ggsave("2.PCAAnalysis/PCHeatmap.png",gg,width = 10,height = 10)

if(NormType=="scale"){
	seuset <- JackStraw(seuset, num.replicate = 100, dims = computePCs)
	seuset <- ScoreJackStraw(seuset, dims = 1:computePCs)

	print(seuset@reductions$pca@jackstraw@overall.p.values)
	write.table(seuset@reductions$pca@jackstraw@overall.p.values,file="2.PCAAnalysis/pcapvalue.txt",sep="\t",quote=F,row.names=F)

}

pcgene <- c()
for(i in 1:computePCs){
    topgenes = TopFeatures(object = seuset[["pca"]], dim = i,balanced=TRUE, nfeatures=60)
    topgenes = data.frame(PC=paste0("PC",i),topgenes$positive,topgenes$negative)
    pcgene = rbind(pcgene,topgenes)
}
write.table(pcgene,file="2.PCAAnalysis/pcGene.txt",sep="\t",quote=F,row.names=F,col.names=T)

if(NormType=="scale"){
	gg = JackStrawPlot(seuset, dims = pcadim)
	ggsave("2.PCAAnalysis/JackStrawPlot.pdf",gg,width = 8,height = 8)
	ggsave("2.PCAAnalysis/JackStrawPlot.png",gg,width = 8,height = 8)
}
gg = ElbowPlot(seuset,ndims = computePCs)
ggsave("2.PCAAnalysis/PCElbowPlot.pdf",gg,width = 8,height = 8)
ggsave("2.PCAAnalysis/PCElbowPlot.png",gg,width = 8,height = 8)

gg = DimPlot(seuset, reduction = "pca")
ggsave("2.PCAAnalysis/PCAPlot.pdf",gg,width = 10,height = 8)
ggsave("2.PCAAnalysis/PCAPlot.png",gg,width = 10,height = 8)

write.table(data.frame(CellName=rownames(seuset@reductions$pca@cell.embeddings),seuset@reductions$pca@cell.embeddings),file="2.PCAAnalysis/pca.txt",sep="\t",quote=F,row.names=F)
write.table(data.frame(CellName=rownames(seuset@reductions$pca@cell.embeddings),seuset@reductions$pca@cell.embeddings),file="pca.txt",sep="\t",quote=F,row.names=F)

print("RunUMAP")

seuset <- RunUMAP(seuset, dims = pcadim, n.neighbors = nneighbors, reduction.name = "umap3d", n.components = 3, reduction.key = "umap3d_")
seuset <- RunUMAP(seuset, dims = pcadim, n.neighbors = nneighbors)
#umap
gg = DimPlot(seuset, reduction = "umap")
ggsave("4.UMAPAnalysis/UMAPPlot.pdf",gg,width = 9,height = 8)
ggsave("4.UMAPAnalysis/UMAPPlot.png",gg,width = 9,height = 8)
#ncounts
gg = metadataplot(seuset, features="nCount_RNA", cols = rainbow(8)[6:1], reduction = "umap")
ggsave("4.UMAPAnalysis/UMAPPlot_nCount_RNA.pdf",gg,width = 9,height = 8)
ggsave("4.UMAPAnalysis/UMAPPlot_nCount_RNA.png",gg,width = 9,height = 8)
#nFeature
gg = metadataplot(seuset, features="nFeature_RNA", cols = rainbow(8)[6:1], reduction = "umap")
ggsave("4.UMAPAnalysis/UMAPPlot_nFeature_RNA.pdf",gg,width = 9,height = 8)
ggsave("4.UMAPAnalysis/UMAPPlot_nFeature_RNA.png",gg,width = 9,height = 8)
#mt
gg = metadataplot(seuset, features="percent.mito", cols = rainbow(8)[6:1], reduction = "umap")
ggsave("4.UMAPAnalysis/UMAPPlot_percent.mito.pdf",gg,width = 9,height = 8)
ggsave("4.UMAPAnalysis/UMAPPlot_percent.mito.png",gg,width = 9,height = 8)

if("PT" %in% ExtraAnalysis){
	gg = metadataplot(seuset, features="percent.pt", cols = rainbow(8)[6:1], reduction = "umap")
	ggsave("4.UMAPAnalysis/UMAPPlot_percent.pt.pdf",gg,width = 9,height = 8)
	ggsave("4.UMAPAnalysis/UMAPPlot_percent.pt.png",gg,width = 9,height = 8)
}
if("CellCycle" %in% ExtraAnalysis){
    gg = DimPlot(seuset, reduction = "umap", group.by = "Phase")
    ggsave("4.UMAPAnalysis/UMAPPlot_CellCyclePhase.pdf",gg,width = 9,height = 8)
    ggsave("4.UMAPAnalysis/UMAPPlot_CellCyclePhase.png",gg,width = 9,height = 8)
    gg = metadataplot(seuset, features="S.Score", cols = rainbow(8)[6:1], reduction = "umap")
    ggsave("4.UMAPAnalysis/UMAPPlot_S.Score.pdf",gg,width = 9,height = 8)
    ggsave("4.UMAPAnalysis/UMAPPlot_S.Score.png",gg,width = 9,height = 8)
    gg = metadataplot(seuset, features="G2M.Score", cols = rainbow(8)[6:1], reduction = "umap")
    ggsave("4.UMAPAnalysis/UMAPPlot_G2M.Score.pdf",gg,width = 9,height = 8)
    ggsave("4.UMAPAnalysis/UMAPPlot_G2M.Score.png",gg,width = 9,height = 8)
    gg = metadataplot(seuset, features="CC.Difference", cols = rainbow(8)[6:1], reduction = "umap")
    ggsave("4.UMAPAnalysis/UMAPPlot_CC.Difference.pdf",gg,width = 9,height = 8)
    ggsave("4.UMAPAnalysis/UMAPPlot_CC.Difference.png",gg,width = 9,height = 8)
}

if("HouseKeeping" %in% ExtraAnalysis){
	gg = metadataplot(seuset, features="percent.HouseKeeping", cols = rainbow(8)[6:1], reduction = "umap")
    ggsave("4.UMAPAnalysis/UMAPPlot_percent.HouseKeeping.pdf",gg,width = 9,height = 8)
    ggsave("4.UMAPAnalysis/UMAPPlot_percent.HouseKeeping.png",gg,width = 9,height = 8)
}

if("dagene" %in% ExtraAnalysis){
	gg = metadataplot(seuset, features="percent.DissociationAssociated", cols = rainbow(8)[6:1], reduction = "umap")
    ggsave("4.UMAPAnalysis/UMAPPlot_percent.DissociationAssociated.pdf",gg,width = 9,height = 8)
    ggsave("4.UMAPAnalysis/UMAPPlot_percent.DissociationAssociated.png",gg,width = 9,height = 8)
	gg = metadataplot(seuset, features="DissociationAssociated.Score", cols = rainbow(8)[6:1], reduction = "umap")
    ggsave("4.UMAPAnalysis/UMAPPlot_DissociationAssociated.Score.pdf",gg,width = 9,height = 8)
    ggsave("4.UMAPAnalysis/UMAPPlot_DissociationAssociated.Score.png",gg,width = 9,height = 8)
}

if(NormType=="sctransform"){
    gg = metadataplot(seuset, features="nCount_SCT", cols = rainbow(8)[6:1], reduction = "umap")
    ggsave("4.UMAPAnalysis/UMAPPlot_nCount_SCT.pdf",gg,width = 9,height = 8)
    ggsave("4.UMAPAnalysis/UMAPPlot_nCount_SCT.png",gg,width = 9,height = 8)
    gg = metadataplot(seuset, features="nFeature_SCT", cols = rainbow(8)[6:1], reduction = "umap")
    ggsave("4.UMAPAnalysis/UMAPPlot_nFeature_SCT.pdf",gg,width = 9,height = 8)
    ggsave("4.UMAPAnalysis/UMAPPlot_nFeature_SCT.png",gg,width = 9,height = 8)
}

dir.create("umapCoordsForBrowser")
write.table(data.frame(CellName=rownames(seuset@reductions$umap@cell.embeddings),seuset@reductions$umap@cell.embeddings),file="4.UMAPAnalysis/umap.txt",sep="\t",quote=F,row.names=F)
write.table(data.frame(CellName=rownames(seuset@reductions$umap@cell.embeddings),seuset@reductions$umap@cell.embeddings),file="umapCoordsForBrowser/coord.2d.txt",sep="\t",quote=F,row.names=F)
write.table(data.frame(CellName=rownames(seuset@reductions$umap3d@cell.embeddings),seuset@reductions$umap3d@cell.embeddings),file="4.UMAPAnalysis/umap3d.txt",sep="\t",quote=F,row.names=F)
write.table(data.frame(CellName=rownames(seuset@reductions$umap3d@cell.embeddings),seuset@reductions$umap3d@cell.embeddings),file="umapCoordsForBrowser/coord.3d.txt",sep="\t",quote=F,row.names=F)

print("RunTSNE")
seuset <- RunTSNE(object = seuset,dims = pcadim,reduction.name = "tsne3d",dim.embed=3, reduction.key = "tSNE3d_",check_duplicates = FALSE)
seuset <- RunTSNE(object = seuset,dims = pcadim,check_duplicates = FALSE)
gg = DimPlot(seuset, reduction = "tsne")
ggsave("3.TSNEAnalysis/TSNEPlot.pdf",gg,width = 9,height = 8)
ggsave("3.TSNEAnalysis/TSNEPlot.png",gg,width = 9,height = 8)
gg = metadataplot(seuset, features="nCount_RNA", cols = rainbow(8)[6:1], reduction = "tsne")
ggsave("3.TSNEAnalysis/TSNEPlot_nCount_RNA.pdf",gg,width = 9,height = 8)
ggsave("3.TSNEAnalysis/TSNEPlot_nCount_RNA.png",gg,width = 9,height = 8)
gg = metadataplot(seuset, features="nFeature_RNA", cols = rainbow(8)[6:1], reduction = "tsne")
ggsave("3.TSNEAnalysis/TSNEPlot_nFeature_RNA.pdf",gg,width = 9,height = 8)
ggsave("3.TSNEAnalysis/TSNEPlot_nFeature_RNA.png",gg,width = 9,height = 8)
gg = metadataplot(seuset, features="percent.mito", cols = rainbow(8)[6:1], reduction = "tsne")
ggsave("3.TSNEAnalysis/TSNEPlot_percent.mito.pdf",gg,width = 9,height = 8)
ggsave("3.TSNEAnalysis/TSNEPlot_percent.mito.png",gg,width = 9,height = 8)


if("PT" %in% ExtraAnalysis){
	gg = metadataplot(seuset, features="percent.pt", cols = rainbow(8)[6:1], reduction = "tsne")
	ggsave("3.TSNEAnalysis/TSNEPlot_percent.pt.pdf",gg,width = 9,height = 8)
	ggsave("3.TSNEAnalysis/TSNEPlot_percent.pt.png",gg,width = 9,height = 8)
}

if("CellCycle" %in% ExtraAnalysis){
    gg = DimPlot(seuset, reduction = "tsne", group.by = "Phase")
    ggsave("3.TSNEAnalysis/TSNEPlot_CellCyclePhase.pdf",gg,width = 9,height = 8)
    ggsave("3.TSNEAnalysis/TSNEPlot_CellCyclePhase.png",gg,width = 9,height = 8)
    gg = metadataplot(seuset, features="S.Score", cols = rainbow(8)[6:1], reduction = "tsne")
    ggsave("3.TSNEAnalysis/TSNEPlot_S.Score.pdf",gg,width = 9,height = 8)
    ggsave("3.TSNEAnalysis/TSNEPlot_S.Score.png",gg,width = 9,height = 8)
    gg = metadataplot(seuset, features="G2M.Score", cols = rainbow(8)[6:1], reduction = "tsne")
    ggsave("3.TSNEAnalysis/TSNEPlot_G2M.Score.pdf",gg,width = 9,height = 8)
    ggsave("3.TSNEAnalysis/TSNEPlot_G2M.Score.png",gg,width = 9,height = 8)
    gg = metadataplot(seuset, features="CC.Difference", cols = rainbow(8)[6:1], reduction = "tsne")
    ggsave("3.TSNEAnalysis/TSNEPlot_CC.Difference.pdf",gg,width = 9,height = 8)
    ggsave("3.TSNEAnalysis/TSNEPlot_CC.Difference.png",gg,width = 9,height = 8)
}
if("HouseKeeping" %in% ExtraAnalysis){
    gg = metadataplot(seuset, features="percent.HouseKeeping", cols = rainbow(8)[6:1], reduction = "tsne")
    ggsave("3.TSNEAnalysis/TSNEPlot_percent.HouseKeeping.pdf",gg,width = 9,height = 8)
    ggsave("3.TSNEAnalysis/TSNEPlot_percent.HouseKeeping.png",gg,width = 9,height = 8)
}
if("dagene" %in% ExtraAnalysis){
    gg = metadataplot(seuset, features="percent.DissociationAssociated", cols = rainbow(8)[6:1], reduction = "tsne")
    ggsave("3.TSNEAnalysis/TSNEPlot_percent.DissociationAssociated.pdf",gg,width = 9,height = 8)
    ggsave("3.TSNEAnalysis/TSNEPlot_percent.DissociationAssociated.png",gg,width = 9,height = 8)
	gg = metadataplot(seuset, features="DissociationAssociated.Score", cols = rainbow(8)[6:1], reduction = "tsne")
    ggsave("3.TSNEAnalysis/TSNEPlot_DissociationAssociated.Score.pdf",gg,width = 9,height = 8)
    ggsave("3.TSNEAnalysis/TSNEPlot_DissociationAssociated.Score.png",gg,width = 9,height = 8)
}
if(NormType=="sctransform"){
    gg = metadataplot(seuset, features="nCount_SCT", cols = rainbow(8)[6:1], reduction = "tsne")
    ggsave("3.TSNEAnalysis/TSNEPlot_nCount_SCT.pdf",gg,width = 9,height = 8)
    ggsave("3.TSNEAnalysis/TSNEPlot_nCount_SCT.png",gg,width = 9,height = 8)
    gg = metadataplot(seuset, features="nFeature_SCT", cols = rainbow(8)[6:1], reduction = "tsne")
    ggsave("3.TSNEAnalysis/TSNEPlot_nFeature_SCT.pdf",gg,width = 9,height = 8)
    ggsave("3.TSNEAnalysis/TSNEPlot_nFeature_SCT.png",gg,width = 9,height = 8)
}

dir.create("tsneCoordsForBrowser")
write.table(data.frame(CellName=rownames(seuset@reductions$tsne@cell.embeddings),seuset@reductions$tsne@cell.embeddings),file="3.TSNEAnalysis/tsne.txt",sep="\t",quote=F,row.names=F)
write.table(data.frame(CellName=rownames(seuset@reductions$tsne@cell.embeddings),seuset@reductions$tsne@cell.embeddings),file="tsneCoordsForBrowser/coord.2d.txt",sep="\t",quote=F,row.names=F)
write.table(data.frame(CellName=rownames(seuset@reductions$tsne3d@cell.embeddings),seuset@reductions$tsne3d@cell.embeddings),file="3.TSNEAnalysis/tsne3d.txt",sep="\t",quote=F,row.names=F)
write.table(data.frame(CellName=rownames(seuset@reductions$tsne3d@cell.embeddings),seuset@reductions$tsne3d@cell.embeddings),file="tsneCoordsForBrowser/coord.3d.txt",sep="\t",quote=F,row.names=F)
#保存RDS	
print("SaveRDS")
saveRDS(seuset, file = paste0("../",prefix,".seuset.rds"))

if(writeMatrix){
	outfilename = paste0(tmppath,"/",prefix,"/AllNormalizedCounts.txt")

	print("writeNormMatrix")
	if(NormType=="sctransform"){
        normdata = seuset@assays$SCT@data
	}else{
        normdata = seuset@assays$RNA@data
	}

	tryCatch({
 		normdata2 = data.frame(GeneID = row.names(normdata),as.matrix(normdata));
 		colnames(normdata2) = c("GeneID",colnames(normdata));
	    print("Write Matrix");
		write.table(normdata2,file=outfilename,sep="\t",quote=F,row.names=F)},
		error=function(e) {
		print(e)
		print("Write Lines")
		txtwrite = file(outfilename,"w");
		writeLines(paste("Gene",str_c(as.character(colnames(normdata)),collapse='\t'),sep="\t"), txtwrite)
		for (i in 1:nrow(normdata)) {
			writeLines( paste(rownames(normdata)[i],str_c(as.character(normdata[i,]),collapse='\t'),sep="\t"), txtwrite);
		}
		close(txtwrite);
	}
	)
}


version = seuset@version
activeassay = seuset@active.assay
info = rbind(paste0("SeuratVersion:",packageVersion("Seurat")),
	paste0("SeuratObjectVersion:",packageVersion("SeuratObject")),
	paste0("DefaultAssay:",activeassay),
	paste0("Norm&ScaleMethod:",NormType),
	paste0("GeneMinCellNum:",geneExpMinCell),
	paste0("CellMinGeneNum:",cellMinGene),
	paste0("CellMaxGeneNum:",cellMaxGene),
	paste0("MaxMTPercent:",MaxMTPercent),
	paste0("computePCs:",computePCs),
	paste0("pcadim:",paste(pcadim,collapse=",")),
	paste0("nneighbors:",nneighbors),
	paste0("VarGeneNum:",varGeneNum),
	paste0("RegressOut:",paste(regressOut,collapse=","))
	)
write.table(info,file=paste0("../",prefix,".parameterInfo.txt"),sep="\t",quote=F,row.names=F,col.names=F)
quit(status=returnValue)







