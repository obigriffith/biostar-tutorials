#Survival Analysi - Kaplan Meier Curve

library(dplyr)
library(survival)

#Start with case predictions
#datadir="/Users/ogriffit/Dropbox/LBNL/Projects/Cepheid/analyzing/analysis_final2/RandomForests/train_survival/unbalanced/final_100k_trees/"
#datadir="/Users/nspies/biostar-tutorials/MachineLearning/"
#setwd(datadir)
case_pred_outfile="testset_CasePredictions.txt"
KMplotfile="KaplanMeier_TestSet_RFRS.pdf"

#read RF results and clinical data from file
print('reading data')
clindata_plusRF=read.table(case_pred_outfile, header = TRUE, na.strings = "NA", sep="\t")

#Create new risk grouping with additional groups
clindata_plusRF = clindata_plusRF %>% 
    rename(t_rfs = time.rfs) %>%                                     # Rename time column for easy scripting
    mutate(
        RF_Group2 = c('low', 'int', 'high')[ntile(Relapse, 3)],      # Add column for new grouping
        e_rfs_10yrcens = ifelse (t_rfs>10, 0, event.rfs)             # Add column of 10yr censored data
        ) 

#First, perform survival analysis for entire patient cohort without down-sampling
#Create survival plot and statistics
#Calculate logrank survival statistic between groups
#Create new dataframe with just necessary data
#surv_data=clindata_plusRF[,c("t_rfs","e_rfs_10yrcens","RF_Group2")]
surv_data = clindata_plusRF %>% 
    select(t_rfs, e_rfs_10yrcens, RF_Group2)

#create a survival object using data
print('calculating survival')
surv_data.surv = with(surv_data, Surv(t_rfs, e_rfs_10yrcens==1))
#Calculate p-value
survdifftest=survdiff(surv_data.surv ~ RF_Group2, data = surv_data)
survpvalue = 1 - pchisq(survdifftest$chisq, length(survdifftest$n) - 1)
survpvalue = format(as.numeric(survpvalue), digits=3)

#Linear test p-value 
#Using the "Score (logrank) test" pvalue from coxph with riskgroup coded as ordinal variable
#See http://r.789695.n4.nabble.com/Trend-test-for-survival-data-td857144.html
#recode  risk groups as 1,2,3
print('Linear test p-value')
surv_data_lin = clindata_plusRF %>%
    select (t_rfs, e_rfs_10yrcens, RF_Group2) %>%
    mutate(RF_Group2 = factor(RF_Group2, levels=c('low', 'int', 'high'))) %>%
    mutate(RF_Group2 = as.numeric(RF_Group2))

survpvalue_linear = summary(coxph(Surv(t_rfs, e_rfs_10yrcens)~RF_Group2, data=surv_data_lin))$sctest[3]
survpvalue_linear = format(as.numeric(survpvalue_linear), digits=3)

##Plot KM curve
print('plotting KM curve')
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
