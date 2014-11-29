setwd("/Users/ogriffit/git/biostar-tutorials/Heatmaps") #Put heatmap.3 code here
library("gplots")

#Create a fake dataset
prob_matrix=replicate(100, rnorm(20)) 
drug_names=paste("drug",letters[1:20],sep="_")
patient_ids=paste("patient",c(1:100),sep="_")
rownames(prob_matrix)=drug_names
colnames(prob_matrix)=patient_ids

#Create fake color side bars
drugclass_colors=sample(c("darkorchid","darkred"), length(drug_names), replace = TRUE, prob = NULL)
drugcategory_colors=sample(c("green","darkgreen"), length(drug_names), replace = TRUE, prob = NULL)
subtype_colors=sample(c("red","blue","cyan","pink","yellow","green"), length(patient_ids), replace = TRUE, prob = NULL)
Mcolors=sample(c("black","white","grey"), length(patient_ids), replace = TRUE, prob = NULL)
Ncolors=sample(c("black","white","grey"), length(patient_ids), replace = TRUE, prob = NULL)
Tcolors=sample(c("black","white","grey"), length(patient_ids), replace = TRUE, prob = NULL)
HER2colors=sample(c("black","white","grey"), length(patient_ids), replace = TRUE, prob = NULL)
PRcolors=sample(c("black","white","grey"), length(patient_ids), replace = TRUE, prob = NULL)
ERcolors=sample(c("black","white","grey"), length(patient_ids), replace = TRUE, prob = NULL)
rlab=t(cbind(drugclass_colors,drugcategory_colors))
clab=cbind(subtype_colors,Mcolors,Ncolors,Tcolors,HER2colors,PRcolors,ERcolors)
rownames(rlab)=c("Class","Category")
colnames(clab)=c("Subtype","M","N","T","HER2","PR","ER")

#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}

#Create heatmap using custom heatmap.3 source code
source("heatmap.3.R")
pdf(file="heatmap3_example.pdf")
main_title="Drug Response Predictions"
par(cex.main=1)
heatmap.3(prob_matrix, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", dendrogram="both", margins=c(6,12), 
Rowv=TRUE, Colv=TRUE, ColSideColors=clab, RowSideColors=rlab, symbreaks=FALSE, key=TRUE, symkey=FALSE, 
density.info="none", trace="none", main=main_title, labCol=FALSE, labRow=drug_names, cexRow=1, col=rev(heat.colors(75)), 
ColSideColorsSize=7, RowSideColorsSize=2, KeyValueName="Prob. Response")
legend("topright",legend=c("Basal","LumA","LumB","Her2","Claudin","Normal","","Positive","Negative","NA","","Targeted","Chemo","","Approved","Experimental"), 
fill=c("red","blue","cyan","pink","yellow","green","white","black","white","grey","white","darkorchid","darkred","white","green","darkgreen"), border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
dev.off()


# example with just single column, single row
mat <- matrix(1:100, byrow=T, nrow=10)
column_annotation <- sample(c("red", "blue", "green"), 10, replace=T)
column_annotation <- as.matrix(column_annotation)
colnames(column_annotation) <- c("Variable X")
 
row_annotation <- sample(c("red", "blue", "green"), 10, replace=T)
row_annotation <- as.matrix(t(row_annotation))
rownames(row_annotation) <- c("Variable Y")
 
heatmap.3(mat, RowSideColors=row_annotation, ColSideColors=column_annotation) 