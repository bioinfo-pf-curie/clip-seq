rm(list=ls())
args <- commandArgs(TRUE)
la <- length(args)
if (la > 0){
for (i in 1:la)
eval(parse(text=args[[i]]))
}

peaks <- read.table(peakFile, header=FALSE, sep="\t")
newCol <- cbind("name_temp_", row.names(peaks))
newCol <- data.frame(name_temp=paste(newCol[,1], newCol[,2],sep=""))
peaks[,4] <- newCol
peaks <- peaks[,-7]
outfile=sub(".bed$", "_intersect.bed", peakFile)
write.table(as.data.frame(peaks), file=outfile, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
