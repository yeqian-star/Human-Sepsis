library(Seurat)
library(getopt)
library(stringr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(scales)


command = matrix(c('inputFile','i',1,"character",
		'rangetab','r',2,"character",
		'rangeOrderAndColor','R',2,"character",
		'clustertab','c',2,"character",
		'clusterOrderAndColor','C',2,"character",
		'compare','v',1,"character",
		'anno','a',2,"character",
		'method','M',1,"character",
		'minpct','m',1,"numeric",
		'logfc','l',1,"numeric",
		'cutoffType','P',1,"character",
		'cutoff','V',1,"numeric",
		'outpath','o',1,"character",
		'base','e',1,"character")
		,byrow=TRUE, ncol=4)
args=getopt(command)
options(stringsAsFactors = FALSE)
rdsFile = args$inputFile;
clustertab = args$clustertab;
rangetab = args$rangetab;
rangeOrderAndColor = args$rangeOrderAndColor;
outpath = args$outpath;
cutoffType = args$cutoffType;
cutoff = args$cutoff;
minpct = args$minpct;
logfc = args$logfc;
method = args$method;
base = args$base;
compareFile = args$compare
annoFile = args$anno;

seuset=readRDS(rdsFile);
allGene = rownames(seuset)

print("Start")
dir.create(outpath)

rangeColor = NULL

setwd(outpath)

if("SCT" %in% names(seuset@assays)){
	DefaultAssay(seuset) = "SCT"
}else if(length(grep("Spatial",names(seuset@assays)))>0){
    assayloc <- min(grep("Spatial",names(seuset@assays)))
    spatialassay = names(seuset@assays)[assayloc]
    DefaultAssay(seuset) = spatialassay 
}else{
	DefaultAssay(seuset) = "RNA"
}

if(base=="e"){
	base=exp(1)
	colFC = "avg_logFC"
}else if(base=="2"){
	base=2
	colFC = "avg_log2FC"
}

isFileExist = function(file){
  if(is.null(file)){
    return(FALSE)
  }
  if(!file.exists(file)){
    return(FALSE)
  }
  return(TRUE)
}

getNum = function(vector = vector, value = value){
    num = table(vector)[value]
    if(is.na(num)){num=0}
    return(as.numeric(num))
}

theme_classic_new <- function(base_size = 11, base_family = "", 
                              base_line_size = base_size/22, 
                              base_rect_size = base_size/22)
{
  theme_classic(base_size = base_size, base_family = base_family, 
                base_line_size = base_line_size, base_rect_size = base_rect_size) %+replace%
    theme(axis.ticks = element_line(colour = "black", size = 0.8),
          axis.line = element_line(size = 0.8),
          axis.ticks.length = unit(1.5, "mm"),
          axis.text = element_text(colour = "black", size = 12),
          axis.title = element_text(size = 20),
          legend.text = element_text(colour = "black", size = 12),
          legend.title = element_text(colour = "black", size = 12))
}

getTopGene <- function(sigDETest,colFC){
  upGene = top_n(sigDETest[sigDETest[,colFC]>0,],5,get(colFC))[,"Gene"]
  downGene = top_n(sigDETest[sigDETest[,colFC]<0,],5,-get(colFC))[,"Gene"]
  topGene = c(upGene[1:min(5,length(upGene))],downGene[1:min(5,length(downGene))])
  return(topGene)
}

createPlotData <- function(DETest, colFC, logfc, cutoffType, cutoff, topGene){
  plotGene = DETest[,c("Gene",colFC,cutoffType)]
  colnames(plotGene) <- c("gene", "plotfc", "FDR_pvalue")
  dat_df <- mutate(plotGene, label = ifelse(gene %in% topGene, gene, NA))
  dat_df <- mutate(dat_df, FDR_type = ifelse(FDR_pvalue >= cutoff, paste0(">=", cutoff), paste0("<", cutoff))) %>% 
    mutate(Style = ifelse(plotfc >= logfc & FDR_pvalue < cutoff, 
                          "Upregulated", ifelse(plotfc < -logfc & FDR_pvalue < cutoff, 
                                                "Downregulated", "Normal")))
  dat_df$FDR_type <- factor(dat_df$FDR_type, levels = c(paste0("<", cutoff),paste0(">=", cutoff)))
 
  return(dat_df) 
}

draw_volcano <- function(dat_df, colFC, logfc, cutoffType, cutoff){
  point_color <- c("red","grey","blue")
  my_plot <- ggplot(dat_df, aes(x = plotfc, y = -log10(FDR_pvalue), color = Style, label = label))
  
  my_plot <- my_plot + 
    geom_point() + scale_color_manual(values = c(point_color[1],point_color[3],point_color[2]), limits = c("Upregulated", "Downregulated", "Normal")) +
    labs(x = colFC, y = paste0("-log10(", cutoffType, ")")) +
    theme_classic_new() + theme(axis.title = element_text(size = 13)) + 
    geom_vline(xintercept = c(-logfc, logfc), linetype = "longdash", size = 0.2) +
    geom_hline(yintercept = c(-log10(cutoff)), linetype = "longdash", size = 0.2)
  
  my_plot <- my_plot +
    geom_text_repel(color = "black", segment.color = "black", 
                    min.segment.length = 0, 
                    force = 10, size = 4,
                    max.overlaps = 10)
  
  
  return(my_plot)
}


if(!is.null(clustertab)){
  newInfo = read.table(clustertab,sep="\t",header=F,comment.char = "")
  realInfo = newInfo[is.finite(match(newInfo[,1],colnames(seuset))),]
  seuset = subset(seuset,cells=as.character(realInfo[,1]))
  extraInfo = realInfo[,2]
  names(extraInfo) = realInfo[,1]
  extraInfo = factor(extraInfo,levels=unique(extraInfo))
  seuset <- AddMetaData(object = seuset, metadata = extraInfo, col.name = "CompareCluster")
}	

if(isFileExist(rangetab)){
	newInfo = read.table(rangetab,sep="\t",header=F,comment.char = "")
	realInfo = newInfo[is.finite(match(newInfo[,1],colnames(seuset))),]
  seuset = subset(seuset,cells=as.character(realInfo[,1]))
  extraInfo = realInfo[,2]
	names(extraInfo) = realInfo[,1]
	extraInfo = factor(extraInfo,levels=unique(extraInfo))
	seuset <- AddMetaData(object = seuset, metadata = extraInfo, col.name = "CompareRange")
	table = table(as.character(seuset@meta.data$CompareCluster),as.character(seuset@meta.data$CompareRange))
}else{
  table = table(as.character(seuset@meta.data$CompareCluster))
}

print(table)

if(isFileExist(rangeOrderAndColor)){
  rangeColor = read.table(rangeOrderAndColor,sep="\t",header=F,comment.char = "")
}

if(isFileExist(annoFile)){
  geneAnno = read.delim(annoFile,he=T, sep="\t", check.names = FALSE,comment.char = "")
  x = match(allGene,geneAnno$AccID)
  z = which(is.na(x))
  geneuid = geneAnno$GeneUID[x]
  geneuid[z] = allGene[z]
}else{
  geneuid = allGene 
}
geneAnno = data.frame(Gene = allGene, GeneUID = geneuid)
rownames(geneAnno) = allGene

allCluster = unique(as.character(seuset@meta.data$CompareCluster))

compare = read.table(compareFile,sep="\t",header=T, comment.char = "",colClasses="character")

allStatistics = c()

for (i in 1:nrow(compare)) {
  case = compare[i,1]
  ctrl = compare[i,2]
  outprefix = compare[i,3]
  compareStatistics = c()
  
  if(!case %in% allCluster | !ctrl %in% allCluster){
    print(paste0("SkipCompare ",outprefix))
    next
  }
  
  caseCells = WhichCells(seuset,expression = CompareCluster==case,return.null =TRUE)
  ctrlCells = WhichCells(seuset,expression = CompareCluster==ctrl,return.null =TRUE)
  
  if(length(caseCells)<3 | length(ctrlCells)<3){
    print(paste0("SkipCompare ",outprefix))
    compareStatistics = rbind(compareStatistics, c(outprefix,"",length(caseCells),length(ctrlCells),minpct,logfc,cutoff,0,0,0))
    allStatistics = rbind(allStatistics, compareStatistics)
    next
  }
  
  dir.create(outprefix)
  print(paste0("Compare:",outprefix))
  DETest = FindMarkers(seuset,ident.1 = case, ident.2 = ctrl, group.by = 'CompareCluster',min.pct = minpct, logfc.threshold = 0, test.use =method,base=base)
  DETest = cbind(geneAnno[rownames(DETest),,drop=F],DETest);
  write.table(DETest,paste0(outprefix,"/",outprefix,".pctFilter.txt"),quote=F,sep="\t",row.name=F);
  sigDETest = DETest[DETest[,cutoffType]<cutoff&abs(DETest[,colFC])>=logfc&(DETest$pct.1>=minpct|DETest$pct.2>=minpct),];
  if(nrow(sigDETest)>0){
    sigDETest[sigDETest[,colFC]>=0,"up_down"] = "upregulated";
    sigDETest[sigDETest[,colFC]<0,"up_down"] = "downregulated";
    upGeneNum = getNum(sigDETest$up_down,"upregulated")
    downGeneNum = getNum(sigDETest$up_down,"downregulated")
    write.table(sigDETest,paste0(outprefix,"/",outprefix,".minpct",minpct,".logfc",logfc,".",cutoffType,cutoff,".txt"),quote=F,sep="\t",row.name=F)
  }else{
    upGeneNum = 0;
    downGeneNum = 0;
  }
  
  compareStatistics = rbind(compareStatistics, c(outprefix,"",length(caseCells),length(ctrlCells),minpct,logfc,cutoff,nrow(DETest),upGeneNum,downGeneNum))
  
  topGene = getTopGene(sigDETest,colFC)
  dat_df = createPlotData(DETest,colFC,logfc,cutoffType,cutoff,topGene)
  gg = draw_volcano(dat_df,colFC,logfc,cutoffType,cutoff)
  
  ggsave(gg, filename = paste0(outprefix,"/",outprefix,".Volcano.png"), 
         width = 8, height = 7)
  ggsave(gg, filename = paste0(outprefix,"/",outprefix,".Volcano.pdf"), 
         width = 8, height = 7)
  
  if(!"CompareRange" %in% colnames(seuset@meta.data)){
    colnames(compareStatistics) = c("Compare","Range","CaseCellNum","CtrlCellNum","MinPct",colFC,paste0(cutoffType,"_Cutoff"),"PctFilteredGene","UpGene","DownGene")
    write.table(compareStatistics,paste0(outprefix,"/",outprefix,".Statistics.txt"),quote=F,sep="\t",row.name=F);
    allStatistics = rbind(allStatistics, compareStatistics)
    next
  }
  
  range = unique(seuset@meta.data$CompareRange)
  
  multiPlot = c()
  
  for(j in range){
    caseCells = WhichCells(seuset,expression = CompareRange==j & CompareCluster==case,return.null =TRUE)
    ctrlCells = WhichCells(seuset,expression = CompareRange==j & CompareCluster==ctrl,return.null =TRUE)
    if(length(caseCells)<3 | length(ctrlCells)<3){
      print(paste0("SkipRangeCompare ",j,"@",outprefix))
      compareStatistics = rbind(compareStatistics, c(outprefix,j,length(caseCells),length(ctrlCells),minpct,logfc,cutoff,0,0,0))
      next
    }
    print(paste0("Compare:",j,"@",outprefix))
    
    subset = subset(seuset,subset = CompareRange==j)
    DETest = FindMarkers(subset,ident.1 = case, ident.2 = ctrl, group.by = 'CompareCluster',min.pct = minpct, logfc.threshold = 0, test.use =method,base=base)
    DETest = cbind(geneAnno[rownames(DETest),,drop=F],DETest);
    write.table(DETest,paste0(outprefix,"/",j,"_Range_",outprefix,".pctFilter.txt"),quote=F,sep="\t",row.name=F);
    sigDETest = DETest[DETest[,cutoffType]<cutoff&abs(DETest[,colFC])>=logfc&(DETest$pct.1>=minpct|DETest$pct.2>=minpct),];
    if(nrow(sigDETest)>0){
      sigDETest[sigDETest[,colFC]>=0,"up_down"] = "upregulated";
      sigDETest[sigDETest[,colFC]<0,"up_down"] = "downregulated";
      upGeneNum = getNum(sigDETest$up_down,"upregulated")
      downGeneNum = getNum(sigDETest$up_down,"downregulated")
      write.table(sigDETest,paste0(outprefix,"/",j,"_Range_",outprefix,".minpct",minpct,".logfc",logfc,".",cutoffType,cutoff,".txt"),quote=F,sep="\t",row.name=F)
    }else{
      upGeneNum = 0;
      downGeneNum = 0;
    }
    compareStatistics = rbind(compareStatistics, c(outprefix,j,length(caseCells),length(ctrlCells),minpct,logfc,cutoff,nrow(DETest),upGeneNum,downGeneNum))
    topGene = getTopGene(sigDETest,colFC)
    dat_df = createPlotData(DETest,colFC,logfc,cutoffType,cutoff,topGene)
    gg = draw_volcano(dat_df,colFC,logfc,cutoffType,cutoff)
    
    ggsave(gg, filename = paste0(outprefix,"/",j,"_Range_",outprefix,".Volcano.png"), 
           width = 8, height = 7)
    ggsave(gg, filename = paste0(outprefix,"/",j,"_Range_",outprefix,".Volcano.pdf"), 
           width = 8, height = 7)
    dat_df$Range = j
    multiPlot = rbind(multiPlot,dat_df)
  }
  plotRange = unique(multiPlot$Range)
  plotColor = hue_pal()(length(plotRange))
  if(!is.null(rangeColor)){
    plotRangeColor = rangeColor[rangeColor$V1 %in% plotRange,]
    if(nrow(plotRangeColor)==length(plotRange)){
      plotRange = plotRangeColor[,1]
      plotColor = plotRangeColor[,2]
    }
  }
  if("" %in% plotColor){
    plotColor = hue_pal()(length(plotRange))
  }
  multiPlot$Range = factor(multiPlot$Range,levels = plotRange)
  
  coldata = aggregate(multiPlot[,"plotfc"],by=multiPlot[,"Range",drop=F],max)
  coldata2 = aggregate(multiPlot[,"plotfc"],by=multiPlot[,"Range",drop=F],min)
  coldata$Range = factor(coldata$Range,levels = plotRange)
  coldata2$Range = factor(coldata2$Range,levels = plotRange)
  coldata$y = 0
  
  multiPlot = multiPlot[order(multiPlot$FDR_type,decreasing = T),]
  
  gg = ggplot() +geom_col(data=coldata,aes(x=Range,y=x),fill="#d3d3d3",alpha=0.5)+geom_col(data=coldata2,aes(x=Range,y=x),fill="grey",alpha=0.5)
  gg = gg + geom_jitter(data = multiPlot, aes(x=Range,y=plotfc,col=FDR_type),position = position_jitter(seed=1))+geom_text_repel(data=multiPlot,aes(x=Range,y=plotfc,label=label),color = "black", segment.color = "black", min.segment.length = 0, force = 10, size = 4,max.overlaps = 10,position = position_jitter(seed=1))
  gg = gg + scale_color_manual(values = c("red","#999999"))
  gg = gg + geom_tile(data = coldata,aes(x=Range,y=y,fill=Range),color="black",height=max(0.25,max(coldata$x)/10),show.legend = F)+scale_fill_manual(values = plotColor)
  gg = gg + geom_text(data = coldata,aes(x=Range,y=y,label=Range),size = 5)
  gg = gg + theme_classic_new() + theme(axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank())
  gg = gg + xlab("")+ylab(colFC)+labs(color=cutoffType)
  
  ggsave(gg, filename = paste0(outprefix,"/Range_",outprefix,".jitter.png"), 
         width = 4+1.2*length(plotRange), height = 8)
  ggsave(gg, filename = paste0(outprefix,"/Range_",outprefix,".jitter.pdf"), 
         width = 4+1.2*length(plotRange), height = 8)
  colnames(compareStatistics) = c("Compare","Range","CaseCellNum","CtrlCellNum","MinPct",colFC,paste0(cutoffType,"_Cutoff"),"PctFilteredGene","UpGene","DownGene")
  write.table(compareStatistics,paste0(outprefix,"/",outprefix,".Statistics.txt"),quote=F,sep="\t",row.name=F);
  allStatistics = rbind(allStatistics, compareStatistics)
}

write.table(allStatistics,"All.statistics.txt",quote=F,sep="\t",row.name=F)
