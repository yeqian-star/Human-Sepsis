library(Seurat)
library(GSVA)
library(limma)
library(getopt)
library(stringr)
command = matrix(c('inputFile','r',1,"character",
		'prefix','p',1,"character",
		'gmtFile','g',1,"character",
		'lineNum',"N",1,"numeric",
		'isRDS',"R",0,"logical",
		'outpath','o',1,"character")
		,byrow=TRUE, ncol=4)
args=getopt(command)
options(stringsAsFactors = FALSE)
inputFile = args$inputFile;
outpath = args$outpath;
gmtFile = args$gmtFile;
prefix = args$prefix;
lineNum = args$lineNum;

if(is.null(args$isRDS)){args$isRDS = FALSE}

print("Start")
dir.create(outpath)
setwd(outpath)


gmtread = function (file)
{
     if (!grepl("\\.gmt$", file)[1]) {
        data = read.table(file,header = T,sep="\t")
	types = as.character(unique(data[,1]))
	geneSetDB = list()
	for(i in 1:length(types)){
		type = types[i]
		geneset = as.character(unique(data[data[,1]==type,2]))
		geneSetDB[[i]]=geneset
	}
	names(geneSetDB) = types 
     }else{
     	geneSetDB = readLines(file)
     	geneSetDB = strsplit(geneSetDB, "\t")
     	names(geneSetDB) = sapply(geneSetDB, "[", 1)
    	geneSetDB = lapply(geneSetDB, "[", -1:-2)
     	geneSetDB = lapply(geneSetDB, function(x) {
        	x[which(x != "")]
     	})
     }
     return(geneSetDB)
}

geneset = gmtread(gmtFile)

if(args$isRDS){
	seuset = readRDS(inputFile)
	data = as.matrix(seuset@data)
	gsva_matrix<- gsva(data, geneset,method="ssgsea")
	gsva_matrix = t(gsva_matrix)
	forwrite = data.frame(CellName=rownames(gsva_matrix),gsva_matrix)
	colnames(forwrite) = c("CellName",colnames(gsva_matrix))
	write.table(forwrite,file=paste0(prefix,".ssgsea.txt"),sep="\t",row.names=F,quot=F)

}else{
	gsva.matrix <- data.frame()
	txtread = file(inputFile,"r")
	title = readLines(txtread,n=1)
	genenames = as.character(unlist(strsplit(title, "[\t]")))
	genenames = genenames[2:length(genenames)]

	line = readLines(txtread,n=lineNum)

	while( length(line) != 0) {
		lineInfos = data.frame(strsplit(line, "[\t]"),stringsAsFactors = F)
		data.use = lineInfos[2:nrow(lineInfos),,drop=F]
		colnames(data.use)=lineInfos[1,]
		rownames(data.use)=genenames

		data = as.matrix(data.use)
		gsva_matrix<- gsva(data, geneset,method="ssgsea",ssgsea.norm=FALSE)

		newdata = t(gsva_matrix)
		gsva.matrix = rbind(gsva.matrix,newdata)
		print(dim(gsva.matrix))
		line=readLines(txtread,n=lineNum)
	}
	close(txtread)
	dif = max(gsva.matrix)-min(gsva.matrix)
	gsva.matrix = gsva.matrix/dif
	forwrite = data.frame(CellName=rownames(gsva.matrix),gsva.matrix)
	colnames(forwrite) = c("CellName",colnames(gsva.matrix))
	write.table(forwrite,file=paste0(prefix,".ssgsea.txt"),sep="\t",row.names=F,quot=F)
}
