#Survival Analysis - Kaplan Meier Curve

#Install required libraries if necessary
install.packages("survival")

#Load required libraries
library(survival)

#Start with case predictions
datadir="/Users/obigriffith/git/biostar-tutorials/MachineLearning"
setwd(datadir)
case_pred_outfile="testset_CasePredictions.txt"
KMplotfile="KaplanMeier_TestSet_RFRS.pdf"

#read RF results and clinical data from file
clindata_plusRF=read.table(case_pred_outfile, header = TRUE, na.strings = "NA", sep="\t")

#Create new risk grouping with additional groups
#Add column for new grouping
quantiles=quantile(clindata_plusRF[,"Relapse"], probs=c(0.33333,0.66667))
clindata_plusRF[,"RF_Group2"]=clindata_plusRF[,"Relapse"]
clindata_plusRF[which(clindata_plusRF[,"Relapse"]<=quantiles[1]),"RF_Group2"]="low"
clindata_plusRF[which(clindata_plusRF[,"Relapse"]>quantiles[1] &  clindata_plusRF[,"Relapse"]<=quantiles[2]),"RF_Group2"]="int"
clindata_plusRF[which(clindata_plusRF[,"Relapse"]>quantiles[2]),"RF_Group2"]="high"

#Rename time column for easy scripting
clindata_plusRF[,"t_rfs"]=clindata_plusRF[,"time.rfs"]

#Add column of 10yr censored data
clindata_plusRF[,"e_rfs_10yrcens"]=clindata_plusRF[,"event.rfs"]
clindata_plusRF[which(clindata_plusRF[,"t_rfs"]>10),"e_rfs_10yrcens"]=0

#First, perform survival analysis for entire patient cohort without down-sampling
#Create survival plot and statistics
#Calculate logrank survival statistic between groups
#Create new dataframe with just necessary data
surv_data=clindata_plusRF[,c("t_rfs","e_rfs_10yrcens","RF_Group2")]

#create a survival object using data
surv_data.surv = with(surv_data, Surv(t_rfs, e_rfs_10yrcens==1))
#Calculate p-value
survdifftest=survdiff(surv_data.surv ~ RF_Group2, data = surv_data)
survpvalue = 1 - pchisq(survdifftest$chisq, length(survdifftest$n) - 1)
survpvalue = format(as.numeric(survpvalue), digits=3)

#Linear test p-value 
#Using the "Score (logrank) test" pvalue from coxph with riskgroup coded as ordinal variable
#See http://r.789695.n4.nabble.com/Trend-test-for-survival-data-td857144.html
#recode  risk groups as 1,2,3
surv_data_lin=clindata_plusRF[,c("t_rfs","e_rfs_10yrcens","RF_Group2")]
surv_data_lin[,"RF_Group2"]=as.vector(surv_data_lin[,"RF_Group2"])
surv_data_lin[which(surv_data_lin[,"RF_Group2"]=="low"),"RF_Group2"]=1
surv_data_lin[which(surv_data_lin[,"RF_Group2"]=="int"),"RF_Group2"]=2
surv_data_lin[which(surv_data_lin[,"RF_Group2"]=="high"),"RF_Group2"]=3
surv_data_lin[,"RF_Group2"]=as.numeric(surv_data_lin[,"RF_Group2"])
survpvalue_linear=summary(coxph(Surv(t_rfs, e_rfs_10yrcens)~RF_Group2, data=surv_data_lin))$sctest[3]
survpvalue_linear = format(as.numeric(survpvalue_linear), digits=3)

##Plot KM curve
krfit.by_RFgroup = survfit(surv_data.surv ~ RF_Group2, data = surv_data)
pdf(file=KMplotfile)
colors = rainbow(5)
title="Survival by RFRS - Test Set"
plot(krfit.by_RFgroup, col = colors, xlab = "Time (Years)", ylab = "Relapse Free Survival", main=title, cex.axis=1.3, cex.lab=1.4)
abline(v = 10, col = "black", lty = 3)

#Set order of categories, categories are by default assigned colors alphabetically by survfit
groups=sort(unique(surv_data[,"RF_Group2"])) #returns unique factor levels sorted alphabetically
names(colors)=groups
groups_custom=c("low","int","high")
colors_custom=colors[groups_custom]
group_sizes_custom=table(surv_data[,"RF_Group2"])[groups_custom]
groups_custom=c("Low","Intermediate","High") #Reset names for consistency with manuscript
legend_text=c(paste(groups_custom, " ", "(", group_sizes_custom, ")", sep=""),paste("p =", survpvalue_linear, sep=" "))
legend(x = "bottomleft", legend = legend_text, col = c(colors_custom,"white"), lty = "solid", bty="n", cex=1.2)
dev.off()