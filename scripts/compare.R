rm(list = ls())

args <- commandArgs(TRUE)
la <- length(args)
if (la > 0){
  for (i in 1:la)
    eval(parse(text = args[[i]]))
}

require(gplots)
require(dplyr)
require(magrittr)
require(DESeq2)
source(polyA_lib)


Total_file <- read.csv(peakfile, header = TRUE, sep = "\t", check.names = FALSE)
rownames(Total_file) <- Total_file$peak

## Get condition
samples_id <- colnames(Total_file)[9:(ncol(Total_file) - 1)]

design <- read.csv(input_list, sep = "\t", check.names = FALSE, header = FALSE, colClasses = c(rep("character",3)))

if (ncol(design == 4)) {
  condition <- design[match(samples_id, design[,1]), 3]
  group <- design[match(samples_id, design[,1]), 4]
}else{
  print("the input_bed_list.txt file is false")
}

cond <- c(unique(condition[group==1]), unique(condition[group==0]))
condition <- as.factor(condition)
print(samples_id)
print(group)
print(condition)
print(cond)

## Filering
### coverage sum by condition > MIN_COUNT_PER_COND
Total_file <- Total_file[which(rowSums(Total_file[,samples_id[group == 0]]) >= as.numeric(min_count_cond) | rowSums(Total_file[,samples_id[group == 1]]) >= as.numeric(min_count_cond)),]
print(dim(Total_file))

## keep gene with at least 2 peaks
Total_file <- Total_file[Total_file$gene %in% Total_file$gene[duplicated(Total_file$gene)],]
Total_file$gene <- factor(Total_file$gene)
print(dim(Total_file))

## Estimate size factor
ddsFirstC <- DESeqDataSetFromMatrix(countData=as.matrix(rbind(Total_file[Total_file$status=="CLIP|",samples_id],Total_file[Total_file$status=="CLIP",samples_id])), colData=data.frame(condition),design = ~ condition)
ddsNormC <- estimateSizeFactors(ddsFirstC)
print(ddsNormC)

ddsFirstR <- DESeqDataSetFromMatrix(countData=as.matrix(rbind(Total_file[Total_file$status=="RNA|",samples_id],Total_file[Total_file$status=="RNA",samples_id])), colData=data.frame(condition),design = ~ condition)
ddsNormR <- estimateSizeFactors(ddsFirstR)
print(ddsNormR)


## DESEQ
## peakIntron vs LastExon (LE sum)
Data4DESEQ <- group_by(Total_file, gene) %>% do(res=PeakIntron_LastPeakSum(.)) %>% extract2(2) %>% rbind_all
rownames(Data4DESEQ) <- paste0(Data4DESEQ$chr,":",Data4DESEQ$start,"-",Data4DESEQ$end,":",Data4DESEQ$peak)

LE <- c(rep("NLE",length(condition)),rep("LE",length(condition)))
colData <- data.frame(t(rbind("condition"=as.vector(condition),LE)))
myData <- as.matrix(Data4DESEQ[,7:(ncol(Data4DESEQ))])
rownames(myData) <- rownames(Data4DESEQ)
ddsT <- DESeqDataSetFromMatrix(countData=myData, colData=colData, design = ~ LE + condition + LE:condition)
sizeFactors(ddsT) <- c(sizeFactors(ddsNormC),sizeFactors(ddsNormR))
print(sizeFactors(ddsT))
myddsT <- DESeq(ddsT)
resT <- results(myddsT,contrast=c(0,0,0,1))

rld <- rlogTransformation(myddsT, blind=TRUE)
png(file=paste(dirname(peakfile),"/MAplot_Total.png",sep=""))
plotMA(resT, alpha=0.05, main="MAplot_Total\nalpha=0.05")
dev.off()
png(file=paste(dirname(peakfile),"/PCAplot.png",sep=""))
print(plotPCA(rld,intgroup="condition"))
dev.off()

## Histogramm of pvalue & padj
pdf(file = paste(dirname(peakfile), "/hist_pval.pdf", sep=""))
hist(resT$pvalue, breaks=100, main="Histogramm of pvalue", xlab="pvalue")
dev.off()

## Hierarchical clustering
cdslog <- log(counts(ddsT)+1)
dist.cor <- 1-cor(cdslog, use="pairwise.complete.obs", method="spearman")
hist.cor <- hclust(as.dist(dist.cor), method="ward.D")
pdf(file = paste(dirname(peakfile), "/Ascending_hierarchical_classification.pdf", sep=""))
plot(hist.cor, main=paste("samples classification by \n", nrow(cdslog), "peaks", sep=""), sub="Distance= 1-correlation")
dev.off()

##################################
NLE_cond1 <- which(colData[,"condition"] == cond[1] & colData[,"LE"] == "NLE")
NLE_cond2 <- which(colData[,"condition"] == cond[2] & colData[,"LE"] == "NLE")
LE_cond1 <- which(colData[,"condition"] == cond[1] & colData[,"LE"] == "LE")
LE_cond2 <- which(colData[,"condition"] == cond[2] & colData[,"LE"] == "LE")
logFC_CLIP <- log2(rowMeans(counts(ddsT,normalized=TRUE)[,NLE_cond1])/rowMeans(counts(ddsT,normalized=TRUE)[,NLE_cond2]))
logFC_RNA <- log2(rowMeans(counts(ddsT,normalized=TRUE)[,LE_cond1])/rowMeans(counts(ddsT,normalized=TRUE)[,LE_cond2]))
Difference_logFC <- logFC_CLIP-logFC_RNA

resO <- cbind(as.data.frame(resT),counts(ddsT,normalized=T), "logFC_CLIP"=logFC_CLIP, "logFC_RNA"=logFC_RNA, "logFC_CLIP-logFC_RNA"=Difference_logFC)
colnames(resO)[6+c(1:ncol(myData))] <- colnames(myData)
outfile <- sub(".bed", "_peakClip_SumRNA.txt", peakfile)
write.table(resO,file=outfile,col.names=TRUE,row.names=TRUE,sep="\t",quote=FALSE)

resP <- cbind(data.frame(Data4DESEQ$chr),Data4DESEQ$start,Data4DESEQ$end,data.frame(Data4DESEQ$strand),resO)
colnames(resP)[1:4] <- c("chr", "start", "end", "strand")
outfile2 <- sub(".bed", "_peakClip_SumRNA_total.txt", peakfile)
write.table(resP,file=outfile2,col.names=TRUE,row.names=TRUE,sep="\t",quote=FALSE)


nameUP_padj <- resO[resO$padj<0.05 & resO$logFC_CLIP-logFC_RNA>0,]
write.table(nameUP_padj, file=paste(dirname(peakfile),"/up_expressed_peaks_padj.txt",sep=""), col.names=T, row.names=T, quote=F, sep="\t")

nameDOWN_padj <- resO[resO$padj<0.05 & resO$logFC_CLIP-logFC_RNA<0,]
write.table(nameDOWN_padj, file=paste(dirname(peakfile),"/down_expressed_peaks_padj.txt",sep=""), col.names=T, row.names=T, quote=F, sep="\t")



## Make plots

pdf(file = paste(dirname(peakfile), "/diffan_heatmap.pdf", sep=""))
par(font.lab=2, cex=.8)

plot(x=logFC_RNA, y=logFC_CLIP,pch=20,main=paste("CLIP-seq vs RNA-seq",paste(cond[1],cond[2],sep="-"), "(FDR 5%)"),xlab="log2 FC RNA-seq",ylab="log2 FC CLIP-seq", frame=FALSE)
points(x=logFC_RNA[which(resT$padj<0.05 & logFC_CLIP-logFC_RNA > 0)],y=logFC_CLIP[which(resT$padj<0.05 & logFC_CLIP-logFC_RNA > 0)],pch=20,col="red")
points(x=logFC_RNA[which(resT$padj<0.05 & logFC_CLIP-logFC_RNA < 0)],y=logFC_CLIP[which(resT$padj<0.05 & logFC_CLIP-logFC_RNA < 0)],pch=20,col="blue")
up <- length(logFC_RNA[which(resT$padj<0.05 & logFC_CLIP-logFC_RNA > 0)])
down <-  length(logFC_RNA[which(resT$padj<0.05 & logFC_CLIP-logFC_RNA < 0)])
noDiff <- length(logFC_RNA[which(resT$padj > 0.05)])
legend("topleft",legend=up,text.col="red",bty="n",cex=1.5)
legend("bottomright",legend=down,text.col="blue",bty="n",cex=1.5)
legend("topright", legend = noDiff, text.col = "black", bty = "n", cex = 1.5)
dev.off()

plot(x=logFC_RNA, y=logFC_CLIP,pch=20,main=paste("CLIP-seq vs RNA-seq",paste(cond[1],cond[2],sep="-"), "(FDR 10%)"),xlab="log2 FC RNA-seq",ylab="log2 FC CLIP-seq", frame=FALSE)
points(x=logFC_RNA[which(resT$padj<0.1 & logFC_CLIP-logFC_RNA > 0)],y=logFC_CLIP[which(resT$padj<0.1 & logFC_CLIP-logFC_RNA > 0)],pch=20,col="red")
points(x=logFC_RNA[which(resT$padj<0.1 & logFC_CLIP-logFC_RNA < 0)],y=logFC_CLIP[which(resT$padj<0.1 & logFC_CLIP-logFC_RNA < 0)],pch=20,col="blue")
up <- length(logFC_RNA[which(resT$padj<0.1 & logFC_CLIP-logFC_RNA > 0)])
down <-  length(logFC_RNA[which(resT$padj<0.1 & logFC_CLIP-logFC_RNA < 0)])
noDiff <- length(logFC_RNA[which(resT$padj > 0.1)])
legend("topleft",legend=up,text.col="red",bty="n",cex=1.5)
legend("bottomright",legend=down,text.col="blue",bty="n",cex=1.5)
legend("topright", legend = noDiff, text.col = "black", bty = "n", cex = 1.5)
dev.off()

plot(x=logFC_RNA, y=logFC_CLIP,pch=20,main=paste("CLIP-seq vs RNA-seq",paste(cond[1],cond[2],sep="-"), "(pvalue 5%)"),xlab="log2 FC RNA-seq",ylab="log2 FC CLIP-seq", frame=FALSE)
points(x=logFC_RNA[which(resT$pvalue<0.05 & Difference_logFC > 0)],y=logFC_CLIP[which(resT$pvalue<0.05 & Difference_logFC > 0)],pch=20,col="red")
points(x=logFC_RNA[which(resT$pvalue<0.05 & Difference_logFC < 0)],y=logFC_CLIP[which(resT$pvalue<0.05 & Difference_logFC < 0)],pch=20,col="blue")
up <- length(logFC_RNA[which(resT$pvalue<0.05 & Difference_logFC> 0)])
down <-  length(logFC_RNA[which(resT$pvalue<0.05 & Difference_logFC< 0)])
noDiff <- length(logFC_RNA[which(resT$pvalue > 0.05)])
legend("topleft",legend=up,text.col="red",bty="n",cex=1.5)
legend("bottomright",legend=down,text.col="blue",bty="n",cex=1.5)
legend("topright", legend = noDiff, text.col = "black", bty = "n", cex = 1.5)

dev.off()

## R pictures
png(file=paste(dirname(peakfile),"/Mean_normalized_counts.png",sep=""))
plot(resT$baseMean+1, -log10(resT$pvalue),log="x", xlab="mean of normalized counts",ylab=expression(-log[10](pvalue)),ylim=c(0,30), cex=.4, col=rgb(0,0,0,.3))
dev.off()

