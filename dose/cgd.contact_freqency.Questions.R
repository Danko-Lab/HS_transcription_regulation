# Contact frequency analysis
#
# Questions:
# 0) Can we predict which genes are regulated by HSF?
# 1) Do we get a better predictor based on HiC data (vs linear distance)?
# 2) Does the peak strength also play a role (in addition to the CF)?
# 3) Does it matter if HSF1 binding site is transcibed?
#
# Columns in data file:
# V1: Gene name
# V2: HSF-dep (0/1)
# V3: HSF-indep (0/1)
# V4: Down reg (0/1)
# V5: Unreg (0/1)
# V6: Chromo
# V7: TSS start
# V8: TSS end
# V9: Log-2 fold change for that gene
# V10: Strand
# V11: Gene class
# V12: Distance to nearest transcribed HSF peak
# V13: Distance to nearest non-transcribed HSF peak
# V14: CF to nearest transcribed
# V15: CF to nearest non-transcribed
# V16: CF to before bin transcribed
# V17: CF to before bin non-transcribed
# V18: CF to after bin transcribed
# V19: CF to after bin non-transcribed
# V20: CF to all transcribed (comma separated list)
# V21: CF to all non-transcribed (comma separated list)
# V22: Peak strength of nearest transcribed peak
# V23: Peak strength of nearest non-transcribed peak
# V24: Peak strengths of all transcribed peaks (comma separated list)
# V25: Peak strengths of all non-transcribed peaks (comma separated list)
# V26: Distance to transcribed site (comma separated list)
# V27: Distance to non-transcribed site (comma separated list)
# V28: PRO-seq reads in the first 100bp from the TSS. Only for up-regulated genes

#####################
## Dependencies
## 
library(DMwR)
library(ROCR)
library(caTools)
require(ggplot2)
library(reshape2)

auc_pr <- function(obs, pred) {
  xx.df <- prediction(pred, obs)
  perf  <- performance(xx.df, "prec", "rec")
  xy    <- data.frame(recall=perf@x.values[[1]], precision=perf@y.values[[1]])
  
  # take out division by 0 for lowest threshold
  xy <- subset(xy, !is.nan(xy$precision))
  
  res   <- trapz(xy$recall, xy$precision)
  res
}

#####################
## Read in data.
##
fileName = "K562_gene_class_contact_frequency_3K_w_PROseq_counts_table.txt"

cf_table = read.table(fileName, sep="\t")
head(cf_table)

#####################
## Sanity check...
largestContactFreq <- sapply(1:NROW(cf_table), function(i) {max(cf_table$V14[i], cf_table$V15[i])})
shortestDistance   <- sapply(1:NROW(cf_table), function(i) {min(cf_table$V12[i], cf_table$V13[i])})
summary(largestContactFreq[cf_table$V2 == 1])
summary(shortestDistance[cf_table$V2 == 1])

#####################
## Define closest
##
indx <- as.double(cf_table[,12] > cf_table[,13])
distClosest <- sapply(1:NROW(cf_table), function(i) {cf_table[i,12+indx[i]]})
cfClosest   <- sapply(1:NROW(cf_table), function(i) {cf_table[i,14+indx[i]]})
psClosest   <- sapply(1:NROW(cf_table), function(i) {cf_table[i,22+indx[i]]})

#####################
## Define dose
##
txDoseNaive <- sapply(1:NROW(cf_table), function(i) {
 txcf <- as.double(strsplit(as.character(cf_table[i,20]), ",")[[1]])
 txps <- as.double(strsplit(as.character(cf_table[i,24]), ",")[[1]])
 
  ntxcf <- as.double(strsplit(as.character(cf_table[i,21]), ",")[[1]])
  ntxps <- as.double(strsplit(as.character(cf_table[i,25]), ",")[[1]])

 sum(txcf*txps)+sum(ntxcf*ntxps)
})

txDoseOptim <- sapply(1:NROW(cf_table), function(i) {
 txcf <- as.double(strsplit(as.character(cf_table[i,20]), ",")[[1]])
 txps <- as.double(strsplit(as.character(cf_table[i,24]), ",")[[1]])
 
 ntxcf <- as.double(strsplit(as.character(cf_table[i,21]), ",")[[1]])
 ntxps <- as.double(strsplit(as.character(cf_table[i,25]), ",")[[1]])
 
 sum((1/ (1 + exp(-0.1*(txcf-84.07369)))* txps))+0.4637866*sum(1/ (1 + exp(-0.1*(ntxcf-84.07369)))* ntxps)
# sum((1/ (1 + exp(-0.06638191*(txcf-40.97855576)))* txps))+0*sum(1/ (1 + exp(-0.06638191*(ntxcf-40.97855576)))* ntxps)
})

plot(txDoseOptim[cf_table$V2 == 1], cf_table$V9[cf_table$V2 == 1], ylab= "HS Fold Change", xlab= "HSF1 Dose", pch=19, col=alpha("gray", 0.8))
cor.test(txDoseOptim[cf_table$V2 == 1], cf_table$V9[cf_table$V2 == 1], method= "pearson")
#cor.test(cf_table$V22[cf_table$V2 == 1], cf_table$V9[cf_table$V2 == 1], method= "pearson") ## About the same.

## Create optimized dataset.
txDose <- txDoseOptim
cf_table <- cbind(cf_table[,c(1:28)], NearestHSF1= cfClosest*psClosest, txDoseNaive, txDose= txDoseOptim)
upreg <- cf_table[cf_table$V2 == 1 | cf_table$V3 == 1,]


################################################################################
## Compare total dose quantity across distinct classes of HS regulated genes.
boxplot(txDose[cf_table$V2 == 1], txDose[cf_table$V3 == 1], txDose[cf_table$V4 == 1], txDose[cf_table$V5 == 1], ylab="HSF-1 Dose")

dat <- melt(list(HSF1dep= txDose[cf_table$V2 == 1], Upindep= txDose[cf_table$V3 == 1], Down= txDose[cf_table$V4 == 1], Unreg= txDose[cf_table$V5 == 1]))
names(dat) <- c("HSF1Dose", "GeneType")

p<-ggplot(dat, aes(x=GeneType, y=HSF1Dose, fill=GeneType)) +
        #geom_violin(trim=FALSE) +
        geom_boxplot(width=0.1, fill="white") +
        labs(title="Pausing index by stage",x="Stage", y = "Pausing index [log-2]")#+
        #ylim(0,7)
p

##############################################################################
## Do multiple HSF-1 binding sites contribute substantially to dose?

## Plot the number of TFBSs near HSF1 depdendent up-regulated genes.
N_HSF1_bs <- sapply(1:NROW(cf_table), function(i) {
 txcf <- as.double(strsplit(as.character(cf_table[i,20]), ",")[[1]])
 ntxcf <- as.double(strsplit(as.character(cf_table[i,21]), ",")[[1]])

 NROW(txcf) #+ NROW(ntxcf)
})

boxplot(N_HSF1_bs[cf_table$V2 == 1], N_HSF1_bs[cf_table$V3 == 1], N_HSF1_bs[cf_table$V4 == 1], N_HSF1_bs[cf_table$V5 == 1], ylab="# HSF-1 binding sites", ylim=c(0,10))

## YES! And this is especially true for transcribed HSF-1 binding sites.


txN_hiD <- sapply(1:NROW(cf_table), function(i) {
 txcf <- as.double(strsplit(as.character(cf_table[i,20]), ",")[[1]])
 txps <- as.double(strsplit(as.character(cf_table[i,24]), ",")[[1]])

 ntxcf <- as.double(strsplit(as.character(cf_table[i,21]), ",")[[1]])
 ntxps <- as.double(strsplit(as.character(cf_table[i,25]), ",")[[1]])

 sum(doses/ sum(doses) > 0.05)
})


## First, visualize all points - nearest vs. total.
plot(upreg$NearestHSF1+1, upreg$txDoseNaive+1, log="xy",pch=19, xlab= "Nearest HSF1 dose", ylab="Total HSF1 dose", col=alpha("gray", 0.2))
points((upreg$NearestHSF1+1)[upreg$V2 == 1], (upreg$txDoseNaive+1)[upreg$V2 == 1], pch=19, col=alpha("blue", 0.2))

txN_hiD <- sapply(1:NROW(cf_table), function(i) {
 txcf <- as.double(strsplit(as.character(cf_table[i,20]), ",")[[1]])
 txps <- as.double(strsplit(as.character(cf_table[i,24]), ",")[[1]])

 ntxcf <- as.double(strsplit(as.character(cf_table[i,21]), ",")[[1]])
 ntxps <- as.double(strsplit(as.character(cf_table[i,25]), ",")[[1]])

# sum(1/ (1 + exp(-0.1*(txcf-84.07369)))* txps > 1)+sum(0.4637866*(1/ (1 + exp(-0.1*(ntxcf-84.07369)))* ntxps) > 1)

 doses <- c(c(1/ (1 + exp(-0.1*(txcf-84.07369)))* txps), c(0.4637866*(1/ (1 + exp(-0.1*(ntxcf-84.07369)))* ntxps)))
 
 sum(doses/ sum(doses) > 0.05)
 #max(doses/ sum(doses)) < 0.8
})

#boxplot(txN[cf_table$V2 == 1], txN[cf_table$V2 == 0])
dose_n_summary <- summary(as.factor(txN_hiD[cf_table$V2 == 1]))
dose_n_summary/ sum(dose_n_summary)
cf_table[cf_table$V2 == 1 & txN_hiD > 1,] ## Which genes does this happen at?!


#############################################
## Do multiple sites contribute to changes after HSF1?!
cf_table <- cf_table[!((txDose_- NearestHSF1_)<0) & cf_table$V9 > 0,] ## for now...

NearestHSF1_ <- cf_table$NearestHSF1 #cf_table$V14*cf_table$V22# #1/ (1 + exp(-0.03479813*(cfClosest-35.46118269)))* psClosest #cfClosest*psClosest
txDose_      <- cf_table$txDoseNaive

indx <- cf_table$V2 == 1

cor.test(txDose_[indx], cf_table$V9[indx], method="pearson")
cor.test(NearestHSF1_[indx], cf_table$V9[indx], method="pearson")
 
cor.test(txDose_[indx], cf_table$V9[indx], method="spearman")
cor.test(NearestHSF1_[indx], cf_table$V9[indx], method="spearman")

## If multiple sites are contributing, residuals for genes w/ large contributions of other sites should have residuals that are systematically negative in  glm(FC ~ Nearest).
model <- glm(cf_table$V9[indx]~NearestHSF1_[indx])
summary(model$residuals)
summary(model$residuals[((txDose_- NearestHSF1_) == 0)[indx]])
summary(model$residuals[((txDose_- NearestHSF1_) > 0)[indx]])

summary(model$residuals[((txDose_- NearestHSF1_) > 1 & (txDose_- NearestHSF1_)<2)[indx]])
summary(model$residuals[((txDose_- NearestHSF1_) > 2 & (txDose_- NearestHSF1_)<50)[indx]])
summary(model$residuals[((txDose_- NearestHSF1_) > 50 & (txDose_- NearestHSF1_)<100)[indx]])
summary(model$residuals[((txDose_- NearestHSF1_) > 100 & (txDose_- NearestHSF1_)<500)[indx]])
summary(model$residuals[((txDose_- NearestHSF1_) > 500 & (txDose_- NearestHSF1_)<1000)[indx]])
summary(model$residuals[((txDose_- NearestHSF1_) > 1000)[indx]])

wilcox.test(model$residuals[((txDose_- NearestHSF1_) == 0)[indx]], model$residuals[((txDose_- NearestHSF1_) > 0)[indx]])

cor.test(model$residuals, ((txDose_- NearestHSF1_)/ (NearestHSF1_+1))[indx], method="spearman")
cor.test(model$residuals/(cf_table$V9[indx]), ((txDose_- NearestHSF1_)/ (NearestHSF1_+1))[indx], method="spearman")

## If we estimate using genes w/ only 1 HSF1 site contributing, we should underestimate residuals for genes where multimple HSF1 sites contribute.
model <- glm(V9~NearestHSF1, data= cf_table[indx & (txDose_- NearestHSF1_) == 0,])
cf_table$NearestHSF1 <- cf_table$txDoseNaive
residuals_ <- cf_table$V9[indx] - predict(model, cf_table[indx,])

summary(residuals_)
summary(residuals_[((txDose_- NearestHSF1_) == 0)[indx]])
summary(residuals_[((txDose_- NearestHSF1_) > 0)[indx]])

wilcox.test(residuals_[((txDose_- NearestHSF1_) == 0)[indx]], residuals_[((txDose_- NearestHSF1_) > 0)[indx]])
cor.test(residuals_, (txDose_- NearestHSF1_)[indx], method="spearman")

summary(txDose_- NearestHSF1_) ## Sanity check...

##### YAY! I'm pretty convinced by this analysis, so long as the fit does not really suck badly.

#####################
## Run a basic model
## 

# Put this code in a for loop and plot distribution
totalNum=NROW(upreg)
trainNum=ceiling(totalNum*0.90)

N <- 200
cf_list <- rep(0, N)
distance_list <- rep(0, N)
dose_list <- rep(0, N)
peak_strength_list <- rep(0, N)
cf_distance_list <- rep(0, N)
cf_distance_peak_strength_list <- rep(0, N)
cf_distance_peak_strength_dose_list <- rep(0, N)


for (i in 1:N){
  #print(paste("i is", i))

  trainSample.idx=sample(1:NROW(upreg), trainNum)  #80% of indices to train on

  model_cf=glm(V2~V14, family="binomial", data=upreg[trainSample.idx,])
  scores_cf=predict(model_cf, upreg[-trainSample.idx,])

  model_distance=glm(V2~V12, family="binomial", data=upreg[trainSample.idx,])
  scores_distance=predict(model_distance, upreg[-trainSample.idx,])
  
  model_dose=glm(V2~txDose, family="binomial", data=upreg[trainSample.idx,])
  scores_dose=predict(model_dose, upreg[-trainSample.idx,])

  model_peak_strength=glm(V2~V22, family="binomial", data=upreg[trainSample.idx,])
  scores_peak_strength=predict(model_peak_strength, upreg[-trainSample.idx,])

  # Combine cf * peak strength
  model_cf_distance=glm(V2~NearestHSF1, family="binomial", data=upreg[trainSample.idx,])
  scores_cf_distance=predict(model_cf_distance, upreg[-trainSample.idx,])

  # Combine cf and distance and peak strength
  model_cf_distance_peak_strength=glm(V2~NearestHSF1+txDoseNaive, family="binomial", data=upreg[trainSample.idx,])
  scores_cf_distance_peak_strength=predict(model_cf_distance_peak_strength, upreg[-trainSample.idx,])

  # Combine cf and distance and peak strength and dosage
  model_cf_distance_peak_strength_dose=glm(V2~txDose+V28, family="binomial", data=upreg[trainSample.idx,])
  scores_cf_distance_peak_strength_dose=predict(model_cf_distance_peak_strength_dose, upreg[-trainSample.idx,])

  # Store results
  cf_list[[i]] <- auc_pr(upreg[-trainSample.idx,2],scores_cf)
  distance_list[[i]] <- auc_pr(upreg[-trainSample.idx,2],scores_distance)
  dose_list[[i]] <- auc_pr(upreg[-trainSample.idx,2],scores_dose)
  peak_strength_list[[i]] <- auc_pr(upreg[-trainSample.idx,2],scores_peak_strength)
  cf_distance_list[[i]] <- auc_pr(upreg[-trainSample.idx,2],scores_cf_distance)
  cf_distance_peak_strength_list[[i]] <- auc_pr(upreg[-trainSample.idx,2],scores_cf_distance_peak_strength)
  cf_distance_peak_strength_dose_list[[i]] <- auc_pr(upreg[-trainSample.idx,2],scores_cf_distance_peak_strength_dose)
}

boxplot(cf_list, distance_list, dose_list, peak_strength_list, cf_distance_list, cf_distance_peak_strength_list, cf_distance_peak_strength_dose_list)



