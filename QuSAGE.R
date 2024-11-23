## Copyright (C) 2019-2020 Novelbio
## Running Platform: Ubuntu / Novelbrain
## Updator: LI, Yadang          Update time: 2020.01.08
## Author:  LI, Yadang          Create time: 2019.05.28
## Instruction: The function of this script is to run QuSAGE algorithm
## Input file1: file_eset, short for expression matrix file. In the matrix, each row specifies a gene and each column is a cell.
## Input file2: file_ccls, short for cell-cluster list file. In the list, each row specifies a cell with the cell name in the first column and the cluster name in the second.
## Input file3: file_gmt, short for GeneSets files with GMT format. Users can also import tsv files, which contain at least two columns, a geneset name, and a gene, as the first two columns respectively. This program will convert the tsv files to GMT format automatically.
## Output: The program will create file folders corresponding to each input gmt file. In each file folders, the result of running QuSAGE algorithm will be saved in ClusterID.qs.rds files.

## Function list:
##   SetGlobalVariable
##   LoadInputPara
##   CreateEsetAndClist
##   CreateGeneSets
##   SetResultFileName
##     multi_conv <- multi_conv
##     aggregateGeneSet.M <- aggregateGeneSet
##     qusage.single.M <- qusage.single
##     qusage.M <- qusage
##   genesetCheck
##   CreateQSresults
##   main

## Global Variable list:
##   crtwd     single character (current work directory)
##   species   single character (human/mouse/rat)
##   file_eset single character (expression matrix file)
##   file_ccls single character (cell-cluster list file)
##   file_gmt  single character (geneset file)
##   keepQSrds          logical (if keep QS rds file)
##   eset                matrix (expression matrix)
##   ccls                  list (cell-cluster list)
##   cset                       (cluster list)
##   clist                 list (clusterID list)
##   genesets              list (geneset list)
##   gslist                list (GeneSets file name list)
##   label_gmt single character (GeneSets file name)
##   ERRgset   character vector (error gene list)
##   errlist               list (clusterID-errGeneSetName list)
##   labelset  character vector (labels("CLUSTER"/"BACKGROUND") for each column in eset)
##   contrast  single character ("CLUSTER-BACKGROUND")

library(qusage)
library(getopt)
library(optparse)

SetGlobalVariable <- function () {
  clist <<- list()
  gslist <<- list()
  errlist <<- list("Cluster__GeneSet")
  Ctrl_stop <<- 0
}

## Func: Load input arguments
LoadInputPara <- function () {
  command = matrix(c('species'       , 's', 1, "character",
                     'file_eset'     , 'e', 1, "character",
                     'file_ccls'     , 'l', 1, "character",
                     'outpath'       , 'o', 1, "character",
                     'outinpath'     , 'i', 1, "character",
                     'binpath'       , 'b', 2, "character",
                     'input_table'   , 't', 2, "character",
                     'input_MSigDB'  , 'm', 2, "character",
                     'input_Novelbio', 'n', 2, "character",
                     'keepQuSAGErds' , 'k', 2, "logical"
                    ), byrow = TRUE, ncol = 4)
  args = getopt(command)
  species <<- args$species
  file_eset <<- args$file_eset
  file_ccls <<- args$file_ccls
  keepQSrds <<- T
  if (!is.null(args$keepQuSAGErds)) keepQSrds <<- as.logical(args$keepQuSAGErds); print(paste("keep QuSAGE rds file:", keepQSrds))
  setwd(args$outpath); crtwd <<- getwd(); print(paste("OutPath: ", crtwd, sep = " "))
  task.path <<- args$outinpath; print(paste("TaskPath:", task.path, sep = " "))
  if (!is.null(args$input_table)) INgmt <- unlist(strsplit(args$input_table, "[,]"))
  if (!is.null(args$input_MSigDB) & species == "human") {
    DBgmt <- paste(args$binpath, "/MSigDB_human/", unlist(strsplit(args$input_MSigDB, "[,]")), sep = "")
  }
  if (!is.null(args$input_MSigDB) & species == "mouse") {
    DBgmt <- paste(args$binpath, "/MSigDB_mouse/", unlist(strsplit(args$input_MSigDB, "[,]")), sep = "")
  }
  if (!is.null(args$input_MSigDB) & species == "rat") {
    DBgmt <- paste(args$binpath, "/MSigDB_rat/", unlist(strsplit(args$input_MSigDB, "[,]")), sep = "")
  }
  if (!is.null(args$input_Novelbio)) {
    NBgmt <- paste(args$binpath, "/NovelbioDB/", unlist(strsplit(args$input_Novelbio, "[,]")), sep = "")
  }
  list_gmt <<- c()
  if (exists("INgmt")) list_gmt <<- c(list_gmt, INgmt); #print(INgmt)
  if (exists("DBgmt")) list_gmt <<- c(list_gmt, DBgmt); #print(DBgmt)
  if (exists("NBgmt")) list_gmt <<- c(list_gmt, NBgmt); #print(NBgmt)
  print('GeneSet File list:'); print(list_gmt); #for (i in list_gmt) print(i)
}

## Function: Get expression matrix and cell-cluster list from input files. If the amount of the cells in the expression matrix and the cell-cluster list are not same, the non-overlapping cells will be removed both in eset (expression matrix) and ccls (list of each cell with its cluster).
CreateEsetAndClist <- function (
    file_eset, ## expression matrix file
    file_ccls  ## cell-cluster list file
  ) {
  print(paste("Loading ESet file start ", date()))
  eset <<- as.matrix(read.table(file_eset, he = T, sep = "\t", row.names = 1, quote = ""))
  coln <<- read.table(file_eset, nrow = 1, he = F, sep = "\t", row.names = 1, colClasses = 'character'); print("Eset size:");print(dim(eset));
  colnames(eset) <- coln
  print(paste("Loading CCLs file start ", date()))
  ccls <<- read.table(file_ccls, he = T, sep = "\t", quote = ""); print(paste("CCls size:", length(ccls[,1]))); #print("ccls:");print(ccls[,2])
  overlap <- intersect(coln, ccls[,1]); print(paste("overlap:  ", length(overlap)))
  e_index <- match(overlap, coln)
  c_index <- match(overlap, ccls[,1])
  eset <<- eset[,e_index]
  cset <<- as.character(ccls[c_index,2]); #print(cset) ## class(cset): character
  if (length(clist)==0) clist <<- unique(cset); #print(paste("clist:", clist)); ## class(clist): character
  clist <<- unlist(clist); print(paste("Cluster Num:", length(clist))); #print("Cluster:");print(clist); ## class(clist): character
  clist <<- rev(clist[order(gsub("^([0-9])$", "0\\1", clist, perl = T))]); print("Cluster:");print(clist);
}

## Function: Write gene set list into genesets based on input gmt / tsv file with GeneSets info
CreateGeneSets <- function (
    file_gmt ## geneset file
  ) {
  #print(paste('file_gmt:', file_gmt))
  if ( grepl("gmt$", file_gmt, perl = T, ignore.case = T)) {
    genesets <<- read.gmt(file_gmt); #print(genesets)
  } else {
    info <- read.table(file_gmt, header = T, sep = "\t")[,1:2]
    colnames(info) <- c("set", "gene")
    genesets <<- list()
    for (i in unique(info$set)) {
      set <- list(as.character(info[info$set==i,2]))
      genesets <<- c(genesets, set)
    }
    names(genesets) <<- unique(info$set)
    #print(genesets)
  }
  #genesets
}


## Function: Set output file names
SetResultFileName <- function (
    cluster ## cluster name
  ) {
  label_gmt <<- sub("\\.\\w+$", "", file_gmt, perl = T, ignore.case = T)
  label_gmt <<- sub(".symbols$", "", label_gmt, perl = T, ignore.case = T)
  label_gmt <<- sub(".v6.2$", "", label_gmt, perl = T, ignore.case = T)
  label_gmt <<- sub("^.+/", "", label_gmt, perl = T, ignore.case = T)  #print(paste('label_gmt:', label_gmt))
  label <- paste("Cluster", cluster, label_gmt, sep = "_")
  #outwd <- paste(crtwd, "QuSAGE_result", label_gmt, sep = "/")
  tempd <- paste(crtwd, "tmp", label_gmt, sep = "/")
  #if (!file.exists(outwd)) dir.create(outwd, recursive = T)
  if (!file.exists(tempd)) dir.create(tempd, recursive = T)
  gslist <<- c (gslist, label_gmt)
}

## Instruction: This function is the original code from qusage library
################################FFT method for convolution
multi_conv<-function(x){ #print(length(x))#print(paste("x:", x))
  x_fft<-apply(x,2,function(x)fft(x, inverse = FALSE)); ## x_fft class: matrix
  M<-max(Mod(x_fft)); #print(table(is.na(Mod(x_fft)))); #print(paste("M:", M))
  x_fft<-x_fft/M; #print(paste("x_fft:", x_fft))
  Prod_fft<-apply(x_fft,1,prod); #print(paste("Prod_fft:", Prod_fft))
  p1<-Re(fft(Prod_fft,inverse=TRUE)); #print(paste("p1:", p1)) #print(fft(Prod_fft,inverse=TRUE))
  N<-nrow(x); #print(paste("N:", N))
  #     if(ncol(x)%%2==0)p1<-c(p1[(N/2+1):N],p1[1:(N/2)])
  Mp1<-which.max(p1); #print(paste("Mp1", Mp1))
  Delta<-N/2-Mp1; #print(paste("Delta:", Delta))
  if(Delta>0){
    p1<-c(p1[(N-Delta):N],p1[1:(N-Delta-1)])
  }
  if(Delta<0){
    p1<-c(p1[(1-Delta):N],p1[1:(1-Delta-1)])
  }
  p2 = abs(p1)
  return(p2/sum(p2))
}

## Instruction: This function is modified code based on aggregateGeneSet from qusage library
aggregateGeneSet.M<-function(geneResults,  ##A QSarray object, as generated by makeComparison
                           geneSets,     ##a list of pathways to be compared, each item in the list is a vector of names that correspond to the gene names from Baseline/PostTreatment
                           n.points=2^12, ##the number of points to sample the convoluted t-distribution at.
                           silent=TRUE   ##If false, print a "." every fifth pathway, as a way to keep track of progress
                          ){

#   NumSDs<-c(20,20,20,10,5,2,2,2,2,1,1,rep(.5,220))*30
  NumSDs<-abs(qt(10^-10,1:250))
  NumSDs[NumSDs>750] = 750
#   ,rep(abs(qt(10^-8,50)),220))

  #saveRDS(geneResults,"generesults.rds")
  Means = geneResults$mean
  SD = geneResults$SD
  DOF=geneResults$dof
  COLS = names(Means)

  if(is.vector(geneSets) & !is.list(geneSets)){
    n = deparse(substitute(geneSets))
    geneSets = list(geneSets)
    names(geneSets) = n
  }
  if(is.null(names(geneSets))){names(geneSets) = 1:length(geneSets)}
  geneSets = lapply(geneSets,function(x){
    if(is.numeric(x)){
      if(any(!(x %in% 1:length(COLS)))){stop("Numeric gene set indices out of bounds")}
      return(x)
    }
    which(COLS%in%x)
  })

  #########First set MaxDiff to adjust to data:
  ##calculate standard deviation
  SumSigma<-sapply(names(geneSets),function(i){
      Indexes = geneSets[[i]]
      x<-sqrt(sum((SD^2*(DOF/(DOF-2)))[Indexes])); #print(paste("SumSigma:", x, "i", i))
      return(x)
  }); #print(paste("SumSigma:", SumSigma))

  MinDof<-sapply(names(geneSets),function(i){
      Indexes = geneSets[[i]]
      if(length(Indexes)==0){return(NA)}
      return(floor(min(DOF[Indexes])))
      }); #print(paste("MinDof:", MinDof))
#   MaxDiff<-pmax(NumSDs[floor(min(DOF))]*SumSigma,1)
  MaxDiff<-pmax(NumSDs[MinDof]*SumSigma,1,na.rm=TRUE); #print(paste("MaxDiff:", MaxDiff))
  #print(paste("NumSDs[MinDof]", NumSDs[MinDof])); print(paste("SumSigma:", SumSigma)); print(paste("qt:", abs(qt(10^-10,1000)))); print(paste("calc:",abs(qt(10^-10,1000))*SumSigma))
# MaxDiff[is.na(NumSDs[MinDof])] = abs(qt(10^-10,1000))*SumSigma; print(paste("MaxDiff:", MaxDiff)) ## Bug created by the original author of QuSAGE
  PDFs = pathMeans = Sizes = NULL
  for(i in 1:length(geneSets)){ #print(paste("geneSets[i]:", geneSets[i]))
    if(!silent & i%%5==0){cat(".")}
    Indexes = geneSets[[i]]; #print(paste("Indexes:", Indexes)); print(paste("length of Indexes:", length(Indexes)))
    if(length(Indexes)!=0){
      Norm<-(2*MaxDiff[i]/{n.points-1}); #print(paste("i:", i, " MaxDiff[i]:", MaxDiff[i], " Norm:", Norm)) #normalize
      PDF<-sapply(Indexes,function(j){
          x = SD[j]; ## x class: numeric
          MaxDiffInt<-MaxDiff[i]/x; #print(paste("MaxDiffInt:", MaxDiffInt))
          #y<-dt(seq(-MaxDiffInt,MaxDiffInt,length.out=n.points),DOF)
          x1<-seq(-MaxDiffInt,MaxDiffInt,length.out=n.points); ## x1 class: numeric
          y<-dt(x1[1:(n.points/2)],DOF[j]); #print(c(paste(names(MaxDiffInt), j), names(DOF[j]), sum(y), DOF[j]))## y class: numeric
					if(sum(y) == 0) { ## return error gene name
						ERRgene <- names(DOF[j])
						names(ERRgene) <- names(MaxDiffInt)
						ERRgset <<- c(ERRgset, ERRgene)
						#print(y)
					}
          y<-c(y,rev(y)); #print(table(is.na(y/sum(y)))) ## y class: numeric
          y/sum(y)/Norm;
      }); #print(table(is.na(PDF))); print(dim(PDF)); saveRDS(PDF,"PDF.rds")
			PDF[is.na(PDF)] = 0; #print(table(is.na(PDF))) ## replace error PDF ## PDF class: matrix
      Tmp<-multi_conv(PDF); #print(paste(length(PDF), dim(PDF)))
      Tmp = Tmp*(n.points-1)/MaxDiff[i]/2
    }else{
      warning(paste("Gene set: (index ",i, ") has 0 overlap with eset.",sep=""))
      #print(paste("Gene set: (index ",i, ") has 0 overlap with eset.",sep=""))
      Tmp = rep(NA, n.points)
    }
    PDFs = cbind(PDFs,Tmp)
    pathMeans = c(pathMeans, mean(Means[Indexes]))
    Sizes = c(Sizes, length(Indexes))
  }

  colnames(PDFs) = names(pathMeans) = names(Sizes) = names(geneSets)
  ##add the new data to the existing QSarray object
  geneResults$pathways = geneSets
  geneResults$path.mean = pathMeans
  geneResults$path.size = Sizes
  geneResults$ranges = MaxDiff/Sizes
  geneResults$n.points=n.points
  geneResults$path.PDF = PDFs
  return(geneResults)
}

## Instruction: This function is modified code based on qusage.single from qusage library
qusage.single.M = function(eset,              ##a matrix of log2(expression values), with rows of features and columns of samples. OR an object of class ExpressionSet
                           labels,            ##vector of labels representing each column of eset.
                           contrast,          ##a string describing which of the groups in 'labels' we want to compare. This is usually of the form 'trt-ctrl', where 'trt' and 'ctrl' are groups represented in 'labels'.
                           geneSets,          ##a list of pathways to be compared. Each item in the list is a vector of names that correspond to the row names of eset.
                           pairVector=NULL,   ##A vector of factors (usually just 1,2,3,etc.) describing the sample pairings. This is often just a vector of patient IDs or something similar. If not provided, all samples are assumed to be independent.
                           var.equal=FALSE,   ##a logical variable indicating whether to treat the two variances as being equal. If TRUE then the pooled variance is used to estimate the variance otherwise the Welch approximation is used.
                           filter.genes=FALSE,##a boolean indicating whether the genes in eset should be filtered to remove genes with low mean and sd.
                           n.points=2^12      ##The number of points to sample the convolution at. Passed to aggregateGeneSet
                 ){
  cat("Calculating gene-by-gene comparisons...")
  results = makeComparison(eset, labels, contrast, pairVector=pairVector,var.equal=var.equal)
  if(filter.genes){
    results = filterGenes(results)
  }

  cat("Done.\nAggregating gene data for gene sets.")
  nu = floor(min(results$dof,na.rm=T))
  if(nu<5){cat("\nLow sample size detected. Increasing n.points in aggregateGeneSet.")}
  results = aggregateGeneSet.M(results, geneSets, silent=F, n.points=n.points)
  cat("Done.\nCalculating variance inflation factors...")
  results = calcVIF(eset, results)
  #cat("Done.\nCalculating homogeneity scores...")
  #results = calcHomogeneity(results, silent=TRUE, addVIF=FALSE)
  #cat("Done.\nCalculating correlation matrix...")
  #results = calcPCor(eset, results)
  cat("Done.\n")
  results
}

## Instruction: This function is modified code based on qusage from qusage library
qusage.M = function(eset,              ##a matrix of log2(expression values), with rows of features and columns of samples. OR an object of class ExpressionSet
                    labels,            ##vector of labels representing each column of eset.
                    contrast,          ##a string describing which of the groups in 'labels' we want to compare. This is usually of the form 'trt-ctrl', where 'trt' and 'ctrl' are groups represented in 'labels'.
                    geneSets,          ##a list of pathways to be compared. Each item in the list is a vector of names that correspond to the row names of eset.
                    pairVector=NULL,   ##A vector of factors (usually just 1,2,3,etc.) describing the sample pairings. This is often just a vector of patient IDs or something similar. If not provided, all samples are assumed to be independent.
                    var.equal=FALSE,   ##a logical variable indicating whether to treat the two variances as being equal. If TRUE then the pooled variance is used to estimate the variance otherwise the Welch approximation is used.
                    filter.genes=FALSE,##a boolean indicating whether the genes in eset should be filtered to remove genes with low mean and sd.
                    n.points=2^12      ##The number of points to sample the convolution at. Passed to aggregateGeneSet
                 ){
  ##determine the type of comparison (i.e. single, double, or more complex?)

  ##if double
  if(grepl("\\s*\\(.+-.+\\)\\s*-\\s*\\(.+-.+\\)\\s*",contrast, perl=TRUE)){

    ##split into two single comparisons
    temp = unlist(strsplit(contrast, "\\)\\s*-\\s*\\("))
    contrast1 = gsub("\\s*\\(","", x=temp[1], perl=TRUE)
    contrast2 = gsub("\\)\\s*","", x=temp[2], perl=TRUE)
    ##flip contrast2
    temp = unlist(strsplit(contrast2, "\\s*\\-\\s*"))
    contrast2 = paste0(temp[2:1], collapse="-")

    # run qusage for contrast1 and contrast2,
    qs.results.1 = qusage.single.M(eset, labels, contrast1, geneSets, pairVector, var.equal, filter.genes, n.points)
    qs.results.2 = qusage.single.M(eset, labels, contrast2, geneSets, pairVector, var.equal, filter.genes, n.points)

    # set the number of samples to be equal (i.e. don't weight the combined pdf)
    n.s = c(qs.results.1$n.samples, qs.results.2$n.samples)
    qs.results.1$n.samples = qs.results.2$n.samples = 1

    qs.results.comb = combinePDFs(list(qs.results.1, qs.results.2), n.points=n.points*2)

    ##double the means and ranges (because we're doing addition, not average)
    qs.results.comb$path.mean = 2 * qs.results.comb$path.mean
    qs.results.comb$ranges = 2 * qs.results.comb$ranges

    ##remove the QSlist (so this object can be considered a normal QSarray)
    qs.results.comb$QSlist = NULL

    ##fix the n.samples
    qs.results.comb$n.samples = mean(n.s)


    ##add other missing bits to the qs.results
    qs.results.comb$contrast = contrast
    for(i in c("var.method","labels","pairVector","pathways","path.size")){
      if(!is.null(qs.results.1[[i]])){qs.results.comb[[i]] = qs.results.1[[i]]}
    }

    return(qs.results.comb)

  ##else, hopefully limma can handle the contrast definition
  }else{
    return(qusage.single.M(eset, labels, contrast, geneSets, pairVector, var.equal, filter.genes, n.points))
  }

}

## Function: remove the unexpressed or nonexistent gene in a certain cluster from the genesets
## Input:    cluster (cluster name)
genesetCheck <- function (cluster) {	#print("Start to check geneset")
  gset <<- genesets
#  for (i in 1:length(gset)) { print(paste(length(gset[[i]]), "genes in", names(gset[[i]])))
#    for (j in 1:length(gset[[i]])) {
#      t = sum(seuset@data[gset[[i]][j], names(seuset@ident[seuset@ident==cluster])])
#      check = c();print(t)
#      if (t == 0) {
#				check[j] = FALSE
#				print(gset[[i]][j])
#			}
#      if (t != 0) check[j] = TRUE
#    }
#    gset[[i]] = gset[[i]][check]; print(length(gset[[i]]))
#  }; print("Finish")

	for (i in 1:length(gset)) {
		#l0 = length(gset[[i]])
		#gene_none <- setdiff(c(gset[[i]]), row.names(seuset@data)); #print(paste("Nonexist:", gene_none))
		#gset[[i]] <<- intersect(c(gset[[i]]), row.names(seuset@data))
		#l1 = length(gset[[i]])
		#print(paste(l0, "->", l1, names(gset[i])))

		#sum = rowSums(seuset@data[,names(seuset@ident[seuset@ident==cluster])])
		#gene_zero <- intersect(c(gset[[i]]), names(sum[sum==0])); print(gene_zero)
    #gset[[i]] <<- setdiff(c(gset[[i]]), gene_zero);
	}; #print("Finish")
}

## Function: Create qs.results and save it into a RDS file
CreateQSresults <- function (
    cluster ## cluster name
  ) {
  outwd <- paste0(task.path, "tmp/", label_gmt)
#  if (!file.exists(outwd) & keepQSrds==T) dir.create(outwd, recursive = T)
  if (!file.exists(outwd)) dir.create(outwd, recursive = T)
  file_qs_result <- paste0(outwd, "/Cluster_", cluster, ".qs.rds") #print(file_qs_result)
  if (!file.exists(file_qs_result)) {
    qs.results <<- qusage.M(eset, labelset, contrast, gset); #print(class(qs.results)); #print(length(gset[[1]]))
    if (!is.null(ERRgset)) {
      print("ERROR gene:")
      print(ERRgset)
      task.color <<- 223
    }
#		if (keepQSrds==T) saveRDS(qs.results, file = file_qs_result)
    saveRDS(qs.results, file = file_qs_result); print(paste("Output:", file_qs_result))
  }
  else {
      print("QSarray file already exists") #print(paste(file_qs_result, "already exists"))
      #qs.results <<- readRDS(file_qs_result)
  }
  #head(qs.results)
  #print(qs.results$path)
}

write_file <- function (file, info) {
  if (file.exists(file)) {
    print(paste("file already exists:", file))
  }
  else {
    for (i in unique(info)) {
      write(i, file = file, ap = T)
    }
    print(paste("finish writing", file))
  }
}

main <- function () {
  print(date())
  print("Running QuSAGE_20200108.R")
  task.color <<- 0
  LoadInputPara()
  CreateEsetAndClist (file_eset, file_ccls)
  for (f in list_gmt) {
    file_gmt <<- f
    print(file_gmt)
    CreateGeneSets(file_gmt)
    for (i in clist) {
      ERRgset <<- c()
      print(paste("Cluster", i, date(), sep = " "))
      labelset <<- cset; #print(labelset)
      labelset[labelset!=i] <<- "BACKGROUND"; #print(labelset)
      labelset[labelset==i] <<- "CLUSTER"; #print(labelset)
      contrast <<- "CLUSTER-BACKGROUND"
      SetResultFileName(i)
      genesetCheck(i)
      CreateQSresults(i)
      if (length(ERRgset) > 0) errlist <<- c(errlist, paste0(i, "__", names(ERRgset)))
    }
    date()
  }
  write_file(paste0(task.path, "/tmp/clist.txt"), clist)
  write_file(paste0(task.path, "/tmp/gslist.txt"), gslist)
  write_file(paste0(task.path, "/tmp/errlist.txt"), errlist)
}

warnings()
SetGlobalVariable()
main();
quit(status = task.color)
