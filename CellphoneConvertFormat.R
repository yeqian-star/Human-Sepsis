#将cellphonev4结果改成cellphonev1的格式
#pvalue文件中interacting_pair含有多段下划线，影响后续genepair拆分
#由于interacting_pair是两个gene名称由_连接而成，reorder过程中按照下划线切分会出现bug
#现计划根据cellphone拼接interacting_pair的逻辑，有gene_a/b用gene_a/b，没有则用partner_a/b
#interacting_pair	id_cp_interaction	 partner_a	 partner_b	 gene_a	gene_b annotation_strategy secreted is_integrin
#
library(getopt)
command = matrix(c('pvaluefile','p',1,"character",
	'meanfile','m',1,"character",
	'outpath','o',1,"character",
	'dbgenefile','d',1,"character")
	,byrow=TRUE,ncol=4)
args=getopt(command)
options(stringsAsFactors = FALSE)
pvaluefile <- args$pvaluefile
meanfile <- args$meanfile
outpath <- args$outpath
dbgenefile <- args$dbgenefile
print("Start converting")
setwd(outpath)
pvalue <- read.table(pvaluefile,sep="\t",header=T)
title1 = read.table(pvaluefile,sep="\t", he=F, nrow = 1,colClasses='character')
colnames(pvalue) <- title1	
clustertitle1 = title1[12:ncol(pvalue)]
sigmean <- read.table(meanfile,sep="\t",header=T)
title2 = read.table(meanfile,sep="\t", he=F, nrow = 1,colClasses='character')
colnames(sigmean) <- title2	
clustertitle2 = title2[13:ncol(sigmean)]
dbgene <- read.table(dbgenefile,sep=",",header=T)

#将pvalue文件拆成info+矩阵，重新组合info信息
pvalue1= pvalue[,1:11]
pvalue2= pvalue[,-1:-11]
pvalue1 <- pvalue[,c(2,1,3,4,5,6,10,7,11)]

ptest <- pvalue1
namecut <- gsub(c("complex:"),"",ptest[,3])
namecut <- gsub(c("simple:"),"",namecut)
ptest[,3]=namecut
namecut2 <- gsub(c("complex:"),"",ptest[,4])
namecut2 <- gsub(c("simple:"),"",namecut2)
ptest[,4]=namecut2

for(i in 1:nrow(ptest)){
	if(ptest[i,5]!=""){
		ptest[i,5] <- dbgene[match(ptest[i,5],dbgene[,4]),3]
	}
	if(ptest[i,6]!=""){
		ptest[i,6] <- dbgene[match(ptest[i,6],dbgene[,4]),3]
	}
}
#将ensembl ID转回gene symbol

for(i in 1:nrow(ptest)){
	name1=""
	name2=""
	if(ptest[i,5]!=""){
		name1=ptest[i,5]
	}else{
		name1=ptest[i,3]
	}
	if(ptest[i,6]!=""){
		name2=ptest[i,6]
	}else{
		name2=ptest[i,4]
	}
	ptest[i,10] <-paste0(name1,"&",name2)
}
pvalue1[,1]<- ptest[,10]
# 20231124 clusterpair可能出现相同的可能： 1：Ins1  2：Gcg_Sst_Ppy   3: Ins1_Gcg 4: Sst_Ppy  
#12 和34 都是Ins1_Gcg_Sst_Ppy,在R里做dataframe操作会在第二个名字后加 .1 做区分
#所以需要在最后写出时重新赋值title去掉.1 的情况
colnames(pvalue2) = clustertitle1
newpvalue <- cbind(pvalue1,pvalue2)

write.table(newpvalue,"pvalues_converted.txt",sep="\t",col.names=T,row.names=F,quot=F)

sigmean1 <- sigmean[,1:12]
sigmean2<- sigmean[,-1:-12]
sigmean1 <-sigmean[,c(2,1,3,4,5,6,10,7,11,12)]
sigtest <- sigmean1
namecut <-gsub(c("complex:"),"",sigtest[,3])
namecut <- gsub(c("simple:"),"",namecut)
sigtest[,3]=namecut
namecut2 <- gsub(c("complex:"),"",sigtest[,4])
namecut2 <- gsub(c("simple:"),"",namecut2)
sigtest[,4]=namecut2
#将ensembl ID转回gene symbol
for(i in 1:nrow(sigtest)){
	if(sigtest[i,5]!=""){
		sigtest[i,5] <- dbgene[match(sigtest[i,5],dbgene[,4]),3]
	}
	if(sigtest[i,6]!=""){
		sigtest[i,6] <- dbgene[match(sigtest[i,6],dbgene[,4]),3]
	}
}

for(i in 1:nrow(sigtest)){
	name1=""
	name2=""
	if(sigtest[i,5]!=""){
		name1=sigtest[i,5]
	}else{
		name1=sigtest[i,3]
	}
	if(sigtest[i,6]!=""){
		name2=sigtest[i,6]
	}else{
		name2=sigtest[i,4]
	}
	sigtest[i,11] <-paste0(name1,"&",name2)
}
sigmean1[,1]<- sigtest[,11]
# 20231124 clusterpair可能出现相同的可能： 1：Ins1  2：Gcg_Sst_Ppy   3: Ins1_Gcg 4: Sst_Ppy  
#12 和34 都是Ins1_Gcg_Sst_Ppy,在R里做dataframe操作会在第二个名字后加 .1 做区分
#所以需要在最后写出时重新赋值title去掉.1 的情况
colnames(sigmean2) = clustertitle2

newsigmean <- cbind(sigmean1,sigmean2)
write.table(newsigmean,"significant_means_converted.txt",sep="\t",col.names=T,row.names=F,quot=F,na="")