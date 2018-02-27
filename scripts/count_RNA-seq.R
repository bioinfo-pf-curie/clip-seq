#!/bin/Rscript
rm(list=ls())

args <- commandArgs(TRUE)

gene_file <- args[2]
file_RNA <- args[3]
outfile <- args[4]

print(gene_file)
print(file_RNA)
print(outfile)
print(class(file_RNA))

require(Rsubread)

##Create annotation file
whole_gene <- read.table(gene_file, sep="\t", header=F)
whole_gene_2 <- whole_gene[,-5]
whole_gene_4_featureCounts <- cbind(whole_gene[,4],whole_gene_2[,-4])
colnames(whole_gene_4_featureCounts) <- c("geneid", "chr", "start", "end", "strand")

##Bam files structure
file_RNA<-strsplit(file_RNA, " ")
if(file_RNA[[1]][1]==""){
	file_RNA<-file_RNA[[1]][-1]
}


countTab <- featureCounts(files=file_RNA, annot.ext=whole_gene_4_featureCounts, strandSpecific=1, reportReads=FALSE)

CountTable <- cbind(countTab$annotation, countTab$counts)

write.table(CountTable, file=outfile, sep="\t", col.names=F, row.names=F, quote=F)
