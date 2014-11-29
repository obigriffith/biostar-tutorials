library("GenomicRanges")
library("IRanges")

datadir="/Users/ogriffit/Desktop/BioStars/"
setwd(datadir)

#Import regions to exclude
Data2Mergefile="IrangesTestData.txt"
Data2Merge=read.table(file=Data2Mergefile, sep="\t", header=TRUE)

#Create GenomicRanges object for manipulation
ranges2merge = GRanges(seqnames = Rle(Data2Merge[,"Id"]),
ranges = IRanges(Data2Merge[,"Start"], end = Data2Merge[,"End"]), 
strand = Rle(strand("*"), length(rownames(Data2Merge) ) ) )

#Simplify ranges to merge overlaps
reduced=reduce(ranges2merge)

#Reformat and print to file
results=data.frame(Id=as.vector(seqnames(reduced)),start=start(reduced),end=end(reduced))
write.table(results, file="merged.txt", sep="\t", row.names=FALSE, quote=FALSE)







