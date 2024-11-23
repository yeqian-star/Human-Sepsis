args <- commandArgs(T)

options(stringsAsFactors = FALSE)
inFile = args[1];
prefixes = args[2];
outFile = args[3];

file = unlist(strsplit(inFile, "[,]"))
prefix = unlist(strsplit(prefixes, "[,]"))

result = data.frame()
for(i in 1:length(file)){
	data = read.csv(file[i],header = T)
	rownames(data) = prefix[i]
	title = read.csv(file[i],colClasses = 'character',nrows = 1,header = F)
	colnames(data)=title
	result = rbind(result,data)
}

finalresult = t(result)
finalresult = data.frame(Sample=rownames(finalresult),finalresult)
colnames(finalresult) = c("Sample",rownames(result))

write.table(finalresult,outFile,sep="\t",quote = F,row.names = F)
