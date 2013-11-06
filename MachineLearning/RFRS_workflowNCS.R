library(randomForest)
library(ROCR)
require(Hmisc)
library(genefilter)
library(mclust)
library(heatmap.plus)
library(fBasics)

#Set working directory and filenames for Input/output
setwd("/Users/nspies/Documents/Genome Center/R Scripts")

datafile="~/Documents/Genome Center/R Scripts/R Data/GSE2034/GSE2034_gcrma.txt" 
clindatafile="~/Documents/Genome Center/R Scripts/R Data/GSE2034/GSE2034_clindetails.txt"

outfile="Cepheid_RFoutput.txt"
case_pred_outfile="Cepheid_CasePredictions.txt"
varimp_pdffile="Cepheid_varImps.pdf"
MDS_pdffile="Cepheid_MDS.pdf"
case_margins_file="Cepheid_Margins.pdf"
ROC_pdffile="Cepheid_ROC.pdf"
vote_dist_pdffile="Cepheid_vote_dist.pdf"
vote_density_pdffile="density_RF.pdf"
margins_mixed_model_clustering_pdffile="Margins_MM.pdf"
top100heatmap_pdffile="top100heatmap.pdf"
new_case_pred_outfile="Cepheid_NewCasePredictions.txt"
combined_case_pred_outfile="Cepheid_CasePredictions_combined.txt"
combined_case_pred_1000downsamp_10yrRFS_outfile="Cepheid_CasePredictions_combined_100downsamp_10yrRFS.txt"

#Read in data (expecting a tab-delimited file with header line and rownames)
data_import=read.table(datafile, header = TRUE, na.strings = "NA", sep="\t")
clin_data_import=read.table(clindatafile, header = TRUE, na.strings = "NA", sep="\t")

#NEED TO GET CLINICAL DATA IN SAME ORDER AS GCRMA DATA
clin_data_order=order(clin_data_import[,"GSM"])
clindata=clin_data_import[clin_data_order,]
data_order=order(colnames(data_import)[4:length(colnames(data_import))])+3 #Order data without first three columns, then add 3 to get correct index in original file
rawdata=data_import[,c(1:3,data_order)] #grab first three columns, and then remaining columns in order determined above

header=colnames(rawdata)

#If there are predictor variables that are constant/invariant, consider removing them
#Preliminary gene filtering
X=rawdata[,4:length(header)]
#Take values and un-log2 them, then filter out any genes according to following criteria (recommended in multtest/MTP documentation): 
#At least 20% of samples should have raw intensity greater than 100 
#The coefficient of variation (sd/mean) is between 0.7 and 10
ffun=filterfun(pOverA(p = 0.2, A = 100), cv(a = 0.7, b = 10))
filt=genefilter(2^X,ffun)
filt_Data=rawdata[filt,] 

#Get potential predictor variables
predictor_data=t(filt_Data[,4:length(header)]) #Filtered
predictor_names=paste(filt_Data[,3]," (",filt_Data[,1],")", sep="") #Filtered, gene symbol + probe ids
colnames(predictor_data)=predictor_names

#Filter down to just cases which have 10yr FU (i.e., exclude NAs)
cases_10yr = !is.na(clindata[,"X10yr_relapse"])
clindata_10yr=clindata[cases_10yr,]
predictor_data_10yr=predictor_data[cases_10yr,]

#Also put aside cases without 10yr FU
cases_no10yr = is.na(clindata[,"X10yr_relapse"])
clindata_no10yr=clindata[cases_no10yr,]
predictor_data_no10yr=predictor_data[cases_no10yr,]

#Get target variable and specify as factor/categorical
target=clindata_10yr[,"X10yr_relapse"] #recurrences after 10yrs not considered events
target[target==0]="NoRelapse"
target[target==1]="Relapse"
target=as.factor(target)

#Run RandomForests
#NOTE: use an ODD number for ntree. When the forest/ensembl is used on test data, ties are broken randomly.
#Having an odd number of trees avoids this issue and makes the model fully deterministic
rf_output=randomForest(x=predictor_data_10yr, y=target, importance = TRUE, ntree = 100001, proximity=TRUE)

#######################################
#Save RF classifier with save()
save(rf_output, file="RF_model")

#Load saved model - this will save time if re-running and you want to skip RF run
load("RF_model")

#Get importance measures
rf_importances=importance(rf_output, scale=FALSE)

#Determine performance statistics
confusion=rf_output$confusion
sensitivity=(confusion[2,2]/(confusion[2,2]+confusion[2,1]))*100
specificity=(confusion[1,1]/(confusion[1,1]+confusion[1,2]))*100
overall_error=rf_output$err.rate[length(rf_output$err.rate[,1]),1]*100
overall_accuracy=1-overall_error
class1_error=paste(rownames(confusion)[1]," error rate= ",confusion[1,3], sep="")
class2_error=paste(rownames(confusion)[2]," error rate= ",confusion[2,3], sep="")
overall_accuracy=100-overall_error

#Prepare stats for output to file
sens_out=paste("sensitivity=",sensitivity, sep="")
spec_out=paste("specificity=",specificity, sep="")
err_out=paste("overall error rate=",overall_error,sep="")
acc_out=paste("overall accuracy=",overall_accuracy,sep="")
misclass_1=paste(confusion[1,2], rownames(confusion)[1],"misclassified as", colnames(confusion)[2], sep=" ")
misclass_2=paste(confusion[2,1], rownames(confusion)[2],"misclassified as", colnames(confusion)[1], sep=" ")

#Prepare confusion table for writing to file
confusion_out=confusion[1:2,1:2]
confusion_out=cbind(rownames(confusion_out), confusion_out)

#Print results to file
write.table(rf_importances[,4],file=outfile, sep="\t", quote=FALSE, col.names=FALSE)
write("confusion table", file=outfile, append=TRUE)
write.table(confusion_out,file=outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE, append=TRUE)
write(c(sens_out,spec_out,acc_out,err_out,class1_error,class2_error,misclass_1,misclass_2), file=outfile, append=TRUE)

#Produce graph of variable importances for top 30 markers
pdf(file=varimp_pdffile)
varImpPlot(rf_output, type=2, n.var=30, scale=FALSE, main="Variable Importance (Gini) for top 30 predictors")
dev.off()

#Produce MDS plot
pdf(file=MDS_pdffile)
target_labels=as.vector(target)
target_labels[target_labels=="NoRelapse"]="N"
target_labels[target_labels=="Relapse"]="R"
MDSplot(rf_output, target, k=2, xlab="", ylab="", pch=target_labels, palette=c("red", "blue"), main="MDS plot")
dev.off()

#Determine which cases were misclassified based on OOB testing and by how much
#Margins represent the fraction of votes for the correct class minus the fraction for the incorrect class
#If class is determined by majority rules, then cases with positive margins are correct and vice versa
margins=margin(rf_output, target)
#combine the patient ids with target/known class, predicted class, votes, and margins
case_predictions=cbind(clindata_10yr[,c(1,3,4,16,17,19)],target,rf_output$predicted,rf_output$votes,as.vector(margins))
misclass_cases=case_predictions[case_predictions[,"target"]!=case_predictions[,"rf_output$predicted"],]

#Write case predictions to file
write.table(case_predictions,file=case_pred_outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#Produce graph of margins for all cases
#Note, when plot is called on a 'margin' object, it uses brewer.pal for color selection.
pdf(file=case_margins_file)
plot(margins, ylab="margins (fraction votes true - fraction votes false)", xlab="patient cases", main="Margins plot")
legend("bottomright", legend=c("NoRelapse","Relapse"), pch=c(20,20), col=c(brewer.pal(3,"Set1")[1], brewer.pal(3,"Set1")[2]))
dev.off()

#Attempt to define cutoffs by mixed-model clustering of votes/probabilities
#Force mclust to break into 3 clusters 
x=as.numeric(case_predictions[,"Relapse"])
mclust_margins=Mclust(x, G=3)
summary(mclust_margins, x) #Gives you list of values returned
classification_margins=mclust_margins$classification
num_clusters=mclust_margins$G

pdf(file=margins_mixed_model_clustering_pdffile)
par(mfrow=c(3,1), oma=c(2,2,2,2))
hist(x[classification_margins==1], xlim=c(0,1), col="blue", xlab="Log2 GCRMA value", main="component 1")
hist(x[classification_margins==2], xlim=c(0,1), col="red", xlab="Log2 GCRMA value", main="component 2")
hist(x[classification_margins==3], xlim=c(0,1), col="green", xlab="Log2 GCRMA value", main="component 3")
title(main="Separation of relapse probabilities by model-based clustering", outer=TRUE)
dev.off()

#Choose cutoffs - use max and min of "middle" cluster (#2) to break into three components
mm_cutoff1=min(x[classification_margins==2])
mm_cutoff2=max(x[classification_margins==2])

#Produce back-to-back histogram of vote distributions for Relapse and NoRelapse
options(digits=2) 
pdf(file=vote_dist_pdffile)
out <- histbackback(split(rf_output$votes[,"Relapse"], target), probability=FALSE, xlim=c(-50,50), main = 'Vote distributions for patients classified by RF', axes=TRUE, ylab="Fraction votes (Relapse)")
#add color
barplot(-out$left, col="red" , horiz=TRUE, space=0, add=TRUE, axes=FALSE) 
barplot(out$right, col="blue", horiz=TRUE, space=0, add=TRUE, axes=FALSE) 
dev.off()

#Create ROC curve plot and calculate AUC
#Can use Malignant vote fractions or Malignant-Benign vote fractions as predictive variable
#The ROC curve will be generated by stepping up through different thresholds for calling malignant vs benign
predictions=as.vector(rf_output$votes[,2])
pred=prediction(predictions,target)
#First calculate the AUC value
perf_AUC=performance(pred,"auc")
AUC=perf_AUC@y.values[[1]]
#Then, plot the actual ROC curve
perf_ROC=performance(pred,"tpr","fpr")
pdf(file=ROC_pdffile)
plot(perf_ROC, main="ROC plot")
text(0.5,0.5,paste("AUC = ",format(AUC, digits=5, scientific=FALSE)))
dev.off()

#Plot density of RF scores for two classes - allows comparison to Oncotype results
dx=case_predictions[,"Relapse"]
classes=case_predictions[,"X10yr_relapse"]
pdf(file=vote_density_pdffile)
plot(density(dx[classes==0],na.rm=T),main='Density of RF scores by 10 year relapse');lines(density(dx[classes==1],na.rm=T),col='red');legend('topright',legend=c('nonrecurrent','recurrent'),col=c('black','red'),pch=20)
dev.off()

#Try running data through forest
RF_predictions_responses=predict(rf_output, predictor_data_no10yr, type="response")
RF_predictions_probs=predict(rf_output, predictor_data_no10yr, type="prob")
RF_predictions_votes=predict(rf_output, predictor_data_no10yr, type="vote")

#Write new case predictions to file
new_case_predictions=cbind(clindata_no10yr[,c(1,3,4,16,17,19)],RF_predictions_responses,RF_predictions_probs)
write.table(new_case_predictions,file=new_case_pred_outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#Combine all case predictions into a single file. Only combine shared columns
case_predictions_for_combine=case_predictions[,c("Study","GSE","GSM","t_rfs","e_rfs","X10yr_relapse","rf_output$predicted","NoRelapse","Relapse")]
colnames(case_predictions_for_combine)=c("Study","GSE","GSM","t_rfs","e_rfs","X10yr_relapse","RF_predictions_responses","NoRelapse","Relapse")
case_predictions_all_combined=rbind(case_predictions_for_combine,new_case_predictions)

#Define risk groups
#risk group 1: Based on cutoffs chosen by mixed-model clustering of margins
case_predictions_all_combined[case_predictions_all_combined[,"Relapse"] <= mm_cutoff1,"RF_group1"]="low"
case_predictions_all_combined[case_predictions_all_combined[,"Relapse"] > mm_cutoff1 & case_predictions_all_combined[,"Relapse"] < mm_cutoff2,"RF_group1"]="int"
case_predictions_all_combined[case_predictions_all_combined[,"Relapse"] >= mm_cutoff2,"RF_group1"]="high"

write.table(case_predictions_all_combined,file=combined_case_pred_outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#Down-sample relapses to create dataset more comparable to Oncotype publication (i.e., 15% overall relapse rate)
#Downsample Relapse cases so that they represent only 15% of total cases:  x / (429+x) = 0.15 [solving for x, ~76]
#Keep all NoRelapses
NoRelapseCases=which(is.na(case_predictions_all_combined[,"X10yr_relapse"]) | case_predictions_all_combined[,"X10yr_relapse"]==0)
RelapseCases=which(case_predictions_all_combined[,"X10yr_relapse"]==1)

#Do multiple downsampling and get N and 10yr relapse rates for each risk group (RF_group1), then average
I=1000
downsampledata=matrix(data=NA, nrow=I, ncol=6)
for (i in 1:I){
 random_RelapseCases=sample(x=RelapseCases, size=76, replace = FALSE, prob = NULL)
 case_predictions_all_combined_down=case_predictions_all_combined[c(NoRelapseCases,random_RelapseCases),]
 low_10yr_relapses=case_predictions_all_combined_down[case_predictions_all_combined_down[,"RF_group1"]=="low","X10yr_relapse"]
 int_10yr_relapses=case_predictions_all_combined_down[case_predictions_all_combined_down[,"RF_group1"]=="int","X10yr_relapse"]
 high_10yr_relapses=case_predictions_all_combined_down[case_predictions_all_combined_down[,"RF_group1"]=="high","X10yr_relapse"]
 perc_10yr_relapse_low=sum(low_10yr_relapses, na.rm=TRUE)/length(low_10yr_relapses)*100
 perc_10yr_relapse_int=sum(int_10yr_relapses, na.rm=TRUE)/length(int_10yr_relapses)*100
 perc_10yr_relapse_high=sum(high_10yr_relapses, na.rm=TRUE)/length(high_10yr_relapses)*100
 downsampledata[i,1:6]=c(perc_10yr_relapse_low,length(low_10yr_relapses),perc_10yr_relapse_int,length(int_10yr_relapses),perc_10yr_relapse_high,length(high_10yr_relapses))
}
colnames(downsampledata)=c("low_perc","low_N","int_perc","int_N","high_perc","high_N")
write.table(downsampledata,file=combined_case_pred_1000downsamp_10yrRFS_outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


#Cluster top 100 genes to look for redundancy
#Manually select best probe (where redundancy exists) from Cepheid_RFoutput.txt
top100file="C:/Users/Obi/Documents/My Dropbox/Projects/Cepheid/analyzing/analysis_final2/RandomForests/train_survival/unbalanced/final_100k_trees/top100_Vars_nonRedundant.txt" #combined (standardCDF + customCDF)
#Read in data (expecting a tab-delimited file with header line and rownames)
top100_data_import=read.table(top100file, header = FALSE, na.strings = "NA", sep="\t")
rownames(top100_data_import)=top100_data_import[,4] #Set gene/probe names as rownames
top100=as.vector(top100_data_import[,4]) #gene/probe labels
top100_predictor_data_10yr=predictor_data_10yr[,top100]


#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {
  dist(c,method="euclidian")
}
myclust=function(c) { 
  hclust(c,method="average")
}

#Perform k-means clustering to break into desired number of clusters
kmeans_clusters=kmeans(t(top100_predictor_data_10yr), 20)

#cluster kmeans clusters to get order to display the clusters
kmeans_cluster_order=heatmap.2(as.matrix(kmeans_clusters$centers), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", dendrogram="both", labCol=NA, RowSideColors=rainbow(20), cexRow=1, col=rev(heat.colors(75)))$rowInd

#Now get order of genes
genes_kmeans_order_new=rep(NA,length(kmeans_clusters$cluster))
topXgenes=rep(NA,length(kmeans_cluster_order))
for (k in 1:length(kmeans_cluster_order)){
 cluster=kmeans_cluster_order[k]
 cluster_genes=kmeans_clusters$cluster[kmeans_clusters$cluster==cluster]
 i=length(which(!is.na(genes_kmeans_order_new)))+1
 j=i+length(cluster_genes)-1
 genes_kmeans_order_new[i:j]=names(cluster_genes)
 topXgenes[k]=names(cluster_genes[1])
}

#Create color side bar for heatmap based on kmeans clusters
colors=rainbow(20)
kmeans_colors=colors[as.vector(kmeans_clusters$cluster)]
#kmeans_order=order(kmeans_clusters$cluster)
#kmeans_colors_ordered=colors[as.vector(kmeans_clusters$cluster)[kmeans_order]]
kmeans_colors_ordered=colors[as.vector(kmeans_clusters$cluster[genes_kmeans_order_new])]

#top100_ordered=top100[kmeans_order] #gene/probe labels - reordered according to kmeans
top100_ordered=genes_kmeans_order_new #gene/probe labels - reordered according to kmeans

#Get actual variable importances, and create colors
top100VarImp=as.vector(top100_data_import[,3]) #ordered already from large to small
top100VarImp_colors=seqPalette(100, name="Blues")[round(top100VarImp*100)]

top100VarImp_ordered=as.vector(top100_data_import[genes_kmeans_order_new,3]) #kmeans order
top100VarImp_colors_ordered=seqPalette(100, name="Blues")[round(top100VarImp_ordered*100)]

#Set target colors for patients
target_colors=as.vector(target)
target_colors[target_colors=="NoRelapse"]="yellow"
target_colors[target_colors=="Relapse"]="red"

risk_group_colors=case_predictions_all_combined[1:325,"RF_group1"]
risk_group_colors[risk_group_colors=="low"]="lightblue"
risk_group_colors[risk_group_colors=="int"]="blue"
risk_group_colors[risk_group_colors=="high"]="darkblue"

#Get actual relapse risk scores so that heatmap can be ordered with this
risk_group_scores=case_predictions_all_combined[1:325,"Relapse"]
risk_order=order(risk_group_scores)
risk_group_colors_ordered=risk_group_colors[risk_order]
target_colors_ordered=target_colors[risk_order]


#Reformat with other heatmap function to allow multiple color side bars
clab=cbind(target_colors,risk_group_colors)
clab_riskorder=cbind(target_colors_ordered,risk_group_colors_ordered)
colnames(clab)=c("Relapse status","Risk group")
colnames(clab_riskorder)=c("Relapse status","Risk group")
rlab=cbind(kmeans_colors, top100VarImp_colors)
rlab_kmeansorder=cbind(kmeans_colors_ordered, top100VarImp_colors_ordered)
colnames(rlab)=c("Kmeans cluster","Var. Imp.")
colnames(rlab_kmeansorder)=c("Kmeans cluster","Var. Imp.")


pdf(file=top100heatmap_pdffile)

#Create separate heatmaps. One with columns/rows ordered by hclust, the other by risk score and kmeans cluster
#Set par margins for heatmap.plus (heatmap.2 has to be handled differently)?
par(oma=c(3,4.5,1,2)) #bottom, left, top, right

heatmap.2(as.matrix(t(top100_predictor_data_10yr)), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", dendrogram="both", labRow=top100, labCol=NA, RowSideColors=kmeans_colors, ColSideColors=risk_group_colors, cexRow=0.40, col=rev(heat.colors(75)))
heatmap.plus(as.matrix(t(top100_predictor_data_10yr)), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", labRow=top100, labCol=NA, ColSideColors=clab, RowSideColors=rlab, cexRow=0.50, col=rev(heat.colors(75)))

#Will need to rearrange the actual patient/gene data and set Colv/Rowv=FALSE because reordering with a vector doesn't seem to work as expected.
#top100_predictor_data_10yr_riskorder=top100_predictor_data_10yr[risk_order,kmeans_order]
top100_predictor_data_10yr_riskorder=top100_predictor_data_10yr[risk_order,genes_kmeans_order_new]
heatmap.2(as.matrix(t(top100_predictor_data_10yr_riskorder)), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", Colv=FALSE, Rowv=FALSE, symbreaks=FALSE, key=TRUE, symkey=FALSE, density.info="none", trace="none", dendrogram="none", labRow=top100_ordered, labCol=NA, RowSideColors=kmeans_colors_ordered, ColSideColors=risk_group_colors_ordered, cexRow=0.40, col=rev(heat.colors(75)))
heatmap.plus(as.matrix(t(top100_predictor_data_10yr_riskorder)), hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", Colv=NA, Rowv=NA, labRow=top100_ordered, labCol=NA, ColSideColors=clab_riskorder, RowSideColors=rlab_kmeansorder, cexRow=0.50, col=rev(heat.colors(75)))
dev.off()

