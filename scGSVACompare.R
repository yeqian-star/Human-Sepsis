library(limma)
library(getopt)
library(stringr)
library(ggplot2)
library(ggthemes)
library(scales)
library(dplyr)
command = matrix(c('inputFile','i',1,"character",
                   'prefix','P',1,"character",
                   'rangetab','r',2,"character",
                   'rangeOrderAndColor','R',2,"character",
                   'clustertab','c',2,"character",
                   'clusterOrderAndColor','C',2,"character",
                   'compare','v',1,"character",
		                'cuttype','t',1,"character",
		                'cutoff','V',1,"numeric",
                   'plottopN','N',1,"numeric",
		                'outpath','o',1,"character")
		                ,byrow=TRUE, ncol=4)
args=getopt(command)
options(stringsAsFactors = FALSE)
inputFile = args$inputFile;
prefix = args$prefix;
clustertab = args$clustertab;
clusterOrderAndColor = args$clusterOrderAndColor;
rangetab = args$rangetab;
rangeOrderAndColor = args$rangeOrderAndColor;
compareFile = args$compare
outpath = args$outpath;
cuttype = args$cuttype;
cutoff = args$cutoff;
topN = args$plottopN;

print("Start")
dir.create(outpath)
setwd(outpath)

datacutoff = function(result=result,cuttype=c("pvalue","fdr"),cutoff=0.05){
	if(cuttype=="pvalue"){
                TF = result$P.Value<cutoff
        }else if(cuttype=="fdr"){
                TF = result$adj.P.Val<cutoff
        }
	data = result[TF,]
	return(data)
}

getTopData = function(result = result,topN = 10){
  result = result[order(result$t,decreasing = T),,drop=F]
  resultUp = result[result$t>0,,drop=F]
  resultUp = top_n(resultUp,topN,t)
  if(nrow(resultUp)>0){
    resultUp = resultUp[1:min(topN,nrow(resultUp)),,drop=F]
  }
  result = result[order(result$t),,drop=F]
  resultdown = result[result$t<0,,drop=F]
  resultdown = top_n(resultdown,topN,-t)
  if(nrow(resultdown)>0){
    resultdown = resultdown[1:min(topN,nrow(resultdown)),,drop=F]
  }
  result = rbind(resultUp,resultdown)
  return(result)
}

ssgseabarplot = function(result=result,cuttype=c("pvalue","fdr"),cutoff=0.05){
	if(cuttype=="pvalue"){
		TF = result$P.Value<cutoff
	}else if(cuttype=="fdr"){	
		TF = result$adj.P.Val<cutoff
	}
	forplot = data.frame(GeneSet = rownames(result),Tvalue=result$t,TF)
	forplot = forplot[order(forplot$Tvalue),]
	forplot$GeneSet = factor(forplot$GeneSet,levels = as.character(forplot$GeneSet))
	limits = 5*ceiling(max(abs(forplot$Tvalue))/5)
	gg = ggplot(data=forplot,mapping = aes(y=Tvalue,x=GeneSet,fill=TF))+geom_col()+ylab("t value")+ylim(-limits,limits)+theme_few()+coord_flip()+scale_fill_discrete(name=paste0(cuttype,"<",cutoff),limits = c('FALSE','TRUE'))
	return(gg)
}

ssgseasigbarplot = function(result=result){
	TF = result$t<=0
        forplot = data.frame(GeneSet = rownames(result),Tvalue=result$t,TF)
        forplot = forplot[order(forplot$Tvalue),]
        forplot$GeneSet = factor(forplot$GeneSet,levels = as.character(forplot$GeneSet))
        limits = 5*ceiling(max(abs(forplot$Tvalue))/5)
        gg = ggplot(data=forplot,mapping = aes(y=Tvalue,x=GeneSet,fill=TF))+geom_col()+ylab("t value")+ylim(-limits,limits)+theme_few()+coord_flip()+theme(legend.position = "none")
        return(gg)
}


if(is.null(inputFile)){
	stop("please check your inputfiles")
}

compare = read.table(compareFile,sep="\t",header=T, comment.char = "",colClasses="character")
clusterInfo = read.table(clustertab,sep="\t",header=T,comment.char = "",row.names = 1)

rangeInfo = NULL
if(!is.null(rangetab)){
    rangeInfo = read.table(rangetab,sep="\t",header=T,comment.char = "",row.names = 1)
}

inputFiles = unlist(strsplit(inputFile,","))
prefixes = unlist(strsplit(prefix,","))

for (num in 1:length(prefixes)) {
    inFile = inputFiles[num]
    fileName = prefixes[num]
    dir.create(fileName)
    gsva_matrix = read.delim(inFile, he=T, sep="\t", row.names=1,comment.char = "",check.names = F)
    gsva_matrix = t(gsva_matrix)
    cellList = colnames(gsva_matrix)
  
    cellList = intersect(cellList,rownames(clusterInfo))
    compareClusterInfo = clusterInfo[cellList,,drop=F]
    allCluster = as.character(unique(compareClusterInfo[,1]))
    allRange = NULL
    
    if(!is.null(rangeInfo)){
      cellList = intersect(cellList,rownames(rangeInfo))
      compareRangeInfo = rangeInfo[cellList,,drop=F]
      allRange = as.character(unique(compareRangeInfo[,1]))
    }
    
    for (i in 1:nrow(compare)) {
      case = compare[i,1]
      ctrl = compare[i,2]
      outprefix = compare[i,3]
      
      if(!case %in% allCluster | !ctrl %in% allCluster){
        print(paste0("SkipCompare ",outprefix))
        next
      }
      
      caseCells = rownames(compareClusterInfo[compareClusterInfo[,1]==case,,drop=F])
      ctrlCells = rownames(compareClusterInfo[compareClusterInfo[,1]==ctrl,,drop=F])
      
      if(length(caseCells)<3 | length(ctrlCells)<3){
        print(paste0("SkipCompare ",outprefix))
        next
      }
      
      dir.create(paste0(fileName,"/",outprefix))
      print(paste0("Compare:",outprefix))
      
      data = gsva_matrix[,c(caseCells,ctrlCells),drop=F]
      design <- cbind(Case=1, CasevsCtrl=c(rep(0, length(caseCells)), rep(1, length(ctrlCells))))
      fit <- lmFit(data, design)
      fit <- eBayes(fit)
      result = topTable(fit, coef="CasevsCtrl",number =nrow(data))
      if(nrow(result)==1){
        result = data.frame(GeneSet=rownames(data),result)
      }else{
        result = data.frame(GeneSet=rownames(result),result)
      }
      write.table(result,paste0(fileName,"/",outprefix,"/",outprefix,".txt"),quote=F,sep="\t",row.name=F)
      gg = ssgseabarplot(getTopData(result,topN),cuttype=cuttype,cutoff=cutoff)
      ggsave(paste0(fileName,"/",outprefix,"/",outprefix,".bar.png"),gg,width=8,height=min(100,max(2,0.2*nrow(result))),limitsize = FALSE)
      ggsave(paste0(fileName,"/",outprefix,"/",outprefix,".bar.pdf"),gg,width=8,height=min(100,max(2,0.2*nrow(result))),limitsize = FALSE)
      result = datacutoff(result,cuttype=cuttype,cutoff=cutoff)
      write.table(result,paste0(fileName,"/",outprefix,"/",outprefix,".",cuttype,cutoff,".txt"),quote=F,sep="\t",row.name=F)
      gg = ssgseasigbarplot(getTopData(result,topN))
      ggsave(paste0(fileName,"/",outprefix,"/",outprefix,".sigbar.png"),gg,width=8,height=min(100,max(2,0.2*nrow(result))),limitsize = FALSE)
      ggsave(paste0(fileName,"/",outprefix,"/",outprefix,".sigbar.pdf"),gg,width=8,height=min(100,max(2,0.2*nrow(result))),limitsize = FALSE)
      
      if(is.null(allRange)){
        next
      }
      
      for(j in allRange){
        rangeCell = rownames(compareRangeInfo[compareRangeInfo[,1]==j,,drop=F])
        rangeCaseCells = intersect(rangeCell,caseCells)
        rangeCtrlCells = intersect(rangeCell,ctrlCells)
        
        if(length(rangeCaseCells)<3 | length(rangeCtrlCells)<3){
          print(paste0("SkipRangeCompare ",j,"@",outprefix))
          next
        }
        print(paste0("Compare:",j,"@",outprefix))
        data = gsva_matrix[,c(rangeCaseCells,rangeCtrlCells),drop=F]
        design <- cbind(Case=1, CasevsCtrl=c(rep(0, length(rangeCaseCells)), rep(1, length(rangeCtrlCells))))
        fit <- lmFit(data, design)
        fit <- eBayes(fit)
        result = topTable(fit, coef="CasevsCtrl",number =nrow(data))
        if(nrow(result)==1){
          result = data.frame(GeneSet=rownames(data),result)
        }else{
          result = data.frame(GeneSet=rownames(result),result)
        }
        write.table(result,paste0(fileName,"/",outprefix,"/",j,"_Range_",outprefix,".txt"),quote=F,sep="\t",row.name=F)
        gg = ssgseabarplot(getTopData(result,topN),cuttype=cuttype,cutoff=cutoff)
        ggsave(paste0(fileName,"/",outprefix,"/",j,"_Range_",outprefix,".bar.png"),gg,width=8,height=min(100,max(2,0.2*nrow(result))),limitsize = FALSE)
        ggsave(paste0(fileName,"/",outprefix,"/",j,"_Range_",outprefix,".bar.pdf"),gg,width=8,height=min(100,max(2,0.2*nrow(result))),limitsize = FALSE)
        result = datacutoff(result,cuttype=cuttype,cutoff=cutoff)
        write.table(result,paste0(fileName,"/",outprefix,"/",j,"_Range_",outprefix,".",cuttype,cutoff,".txt"),quote=F,sep="\t",row.name=F)
        gg = ssgseasigbarplot(getTopData(result,topN))
        ggsave(paste0(fileName,"/",outprefix,"/",j,"_Range_",outprefix,".sigbar.png"),gg,width=8,height=min(100,max(2,0.2*nrow(result))),limitsize = FALSE)
        ggsave(paste0(fileName,"/",outprefix,"/",j,"_Range_",outprefix,".sigbar.pdf"),gg,width=8,height=min(100,max(2,0.2*nrow(result))),limitsize = FALSE)
        
      }
    }
}

