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
# V29: 'body_count_HS_br1'
# V30: 'body_count_HS_br2'
# V31: 'body_count_NHS_br1'
# V32: 'body_count_NHS_br2'
# V33: 'pause_count_HS_br1'
# V34: 'pause_count_HS_br2'
# V35: 'pause_count_NHS_br1'
# V36: 'pause_count_NHS_br2'
# V37: 'postcps_count_HS_br1'
# V38: 'postcps_count_HS_br2'
# V39: 'postcps_count_NHS_br1'
# V40: 'postcps_count_NHS_br2'
# V41: 'body_rpkm_HS_br1'
# V42: 'body_rpkm_HS_br2'
# V43: 'body_rpkm_NHS_br1'
# V44: 'body_rpkm_NHS_br2'
# V45: 'pause_rpkm_HS_br1'
# V46: 'pause_rpkm_HS_br2'
# V47: 'pause_rpkm_NHS_br1'
# V48: 'pause_rpkm_NHS_br2'
# V49: 'postcps_rpkm_HS_br1'
# V50: 'postcps_rpkm_HS_br2'
# V51: 'postcps_rpkm_NHS_br1'
# V52: 'postcps_rpkm_NHS_br2'

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
fileName = "K562_gene_class_contact_frequency_3K_w_PROseq_counts_and_rpkm_minus_HS_table.txt"

cf_table = read.table(fileName, sep="\t")
head(cf_table)

## Get training/ test definitions that Paul is using.
tset <- read.table("training_and_test_data_w_geneNames.txt", header=TRUE, sep="\t")
#sum(as.character(cf_table$V1[cf_table$V2 == 1 | cf_table$V3 == 1]) == tset$X0)/NROW(tset) ## = 1; suggests we don't have to reorder.

#####################
## Sanity check...
largestContactFreq <- sapply(1:NROW(cf_table), function(i) {max(cf_table$V14[i], cf_table$V15[i])})
shortestDistance   <- sapply(1:NROW(cf_table), function(i) {min(cf_table$V12[i], cf_table$V13[i])})
summary(largestContactFreq[cf_table$V2 == 1])
summary(shortestDistance[cf_table$V2 == 1])

cf_table[(txDose < 0.1 & cf_table$V2 == 1 & shortestDistance < 1000),]

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

#############################################
## Optimize a and b by gradient descent.

#sample(which(cf_table$V2 == 1 | cf_table$V3 == 1), 0.2*sum(cf_table$V2 == 1 | cf_table$V3 == 1), replace=FALSE)
indx <- which(cf_table$V2 == 1 | cf_table$V3 == 1)[tset$is_train == "True"] 
indx_tst <- which(cf_table$V2 == 1 | cf_table$V3 == 1)[tset$is_train == "False"]

dose <- function(x) {
 b <- x[1] #r2= 81.28158 #r1= 71.6659 #d= 20 # 
 a <- x[2] #r2= 0.326084 #r1= 0.6096329 #d= 0.1  #
 z <- x[3] #r2= 0.6083511 #r1= 0.6177003 # # d= 0 #

 cf_small <- cf_table[indx,] #[cf_table$V2 == 1 | cf_table$V3 == 1,]

 txDose <- sapply(1:NROW(cf_small), function(i) {
  txcf <- as.double(strsplit(as.character(cf_small[i,20]), ",")[[1]])
  txps <- as.double(strsplit(as.character(cf_small[i,24]), ",")[[1]])
  
  ntxcf <- as.double(strsplit(as.character(cf_small[i,21]), ",")[[1]])
  ntxps <- as.double(strsplit(as.character(cf_small[i,25]), ",")[[1]])

  sum(1/ (1 + exp(-a*(txcf-b)))* txps) + z*sum(1/ (1 + exp(-a*(ntxcf-b)))* ntxps)
 })

 1- auc_pr(cf_small$V2 ==1, txDose)
}
optim(c(0), dose, method= "Brent", lower= c(0), upper= c(1))
optim(c(20), dose, method= "Brent", lower= c(0), upper= c(500))


optim(c(20), dose, method= "CG")
#optim(c(80, 0.1, 0.5), dose, method= "CG")
optim(c(150, 0.1, 1), dose, method= "L-BFGS-B", lower= c(0, 0, 0), upper= c(500, 1, 1))


####################
## Incorporate order... And down-weight strenth of far out sites.

## This is the functino that we want. The 2 and 5 are free parameters to tune.
sr <- c(1:10)
plot((1/(1+exp(-2*(5-sr))))~sr)

dose <- function(x) {
 k <- x[1]
 sl<- x[2]

 b <- 81.28158 #
 a <- 0.326084
 z <- 0.6083511 # 

 cf_small <- cf_table[cf_table$V2 == 1 | cf_table$V3 == 1,]

 txDose <- sapply(1:NROW(cf_small), function(i) {
  txcf <- as.double(strsplit(as.character(cf_small[i,20]), ",")[[1]])
  txps <- as.double(strsplit(as.character(cf_small[i,24]), ",")[[1]])
  txdists <- as.double(strsplit(as.character(cf_small[i,26]), ",")[[1]])
  
  ntxcf <- as.double(strsplit(as.character(cf_small[i,21]), ",")[[1]])
  ntxps <- as.double(strsplit(as.character(cf_small[i,25]), ",")[[1]])
  ntxdists <- as.double(strsplit(as.character(cf_small[i,27]), ",")[[1]])
  
  OD <- order(c(txdists, ntxdists))
  
  sum(1/ (1 + exp(-a*(txcf-b))) * 1/ (1 + exp(-sl*(k-OD))) * txps) + z*sum(1/ (1 + exp(-a*(ntxcf-b))) * 1/ (1 + exp(-sl*(k-OD))) * ntxps)
 })

 1- auc_pr(cf_small$V2 ==1, txDose)
}

optim(c(5, 1), dose, method= "L-BFGS-B", lower= c(0, 0), upper= c(20, 1))


optim(c(1, 0.5), dose, method= "Brent", lower= c(0, 0), upper= c(500, 1))
optim(c(1), dose, method= "L-BFGS-B", lower= c(0), upper= c(500))
#optim(c(85, 0.1), dose, method= "L-BFGS-B", lower= c(0, 0), upper= c(500, 1))

dose <- function(x) { ## Optimize recovery of fold-change.
 b <- x[1] # 40.97855576 #
 a <- x[2] # 0.06638191
# z <- x[3] # 0 # 
 cf_small <- cf_table[cf_table$V2 == 1,]

 txDose <- sapply(1:NROW(cf_small), function(i) {
  txcf <- as.double(strsplit(as.character(cf_small[i,20]), ",")[[1]])
  txps <- as.double(strsplit(as.character(cf_small[i,24]), ",")[[1]])
  
  ntxcf <- as.double(strsplit(as.character(cf_small[i,21]), ",")[[1]])
  ntxps <- as.double(strsplit(as.character(cf_small[i,25]), ",")[[1]])

  sum(1/ (1 + exp(-a*(txcf-b)))* txps) #+ z* sum(1/ (1 + exp(-a*(ntxcf-b)))* ntxps)
 })

 -cor(cf_small$V9, txDose, method="pearson")
}
optim(c(80, 0.1), dose, method= "L-BFGS-B", lower= c(0, 0), upper= c(500, 1))

optim(c(80, 0.1, 0.5), dose, method= "L-BFGS-B", lower= c(0, 0, 0), upper= c(500, 1, 1))

closetss <- function(x) {
 a <- x[2]
 b <- x[1]
 cf_small <- cf_table[cf_table$V2 == 1,]
 
 model2 <- 1/ (1 + exp(-a*(cfClosest[cf_table$V2 == 1]-b)))* psClosest[cf_table$V2 == 1]
 
 -cor(cf_small$V9, model2, method="pearson")
}
optim(c(85, 0.1), closetss, method= "L-BFGS-B", lower= c(0, 0), upper= c(500, 1))

closetss <- function(x) {
 a <- x[2]
 b <- x[1]
 cf_small <- cf_table[cf_table$V2 == 1 | cf_table$V3 == 1,]
 
 model2 <- 1/ (1 + exp(-a*(cfClosest[cf_table$V2 == 1 | cf_table$V3 == 1]-b)))* psClosest[cf_table$V2 == 1 | cf_table$V3 == 1]
 
 1- auc_pr(cf_small$V2 ==1, model2)
}
optim(c(85, 0.1), closetss, method= "L-BFGS-B", lower= c(0, 0), upper= c(500, 1))

#######################
## Compute the important values

txDoseOptim <- sapply(1:NROW(cf_table), function(i) {

 b <- 81.28158
 a <- 0.326084	
 z <- 0.6083511 

 txcf <- as.double(strsplit(as.character(cf_table[i,20]), ",")[[1]])
 txps <- as.double(strsplit(as.character(cf_table[i,24]), ",")[[1]])
 
 ntxcf <- as.double(strsplit(as.character(cf_table[i,21]), ",")[[1]])
 ntxps <- as.double(strsplit(as.character(cf_table[i,25]), ",")[[1]])
 
 sum((1/ (1 + exp(-a*(txcf-b)))* txps))+z*sum(1/ (1 + exp(-a*(ntxcf-b)))* ntxps)
})

plot(txDoseOptim[cf_table$V2 == 1], cf_table$V9[cf_table$V2 == 1], ylab= "HS Fold Change", xlab= "HSF1 Dose", pch=19, col=alpha("gray", 0.8))
cor.test(txDoseOptim[cf_table$V2 == 1], cf_table$V9[cf_table$V2 == 1], method= "pearson")
#cor.test(cf_table$V22[cf_table$V2 == 1], cf_table$V9[cf_table$V2 == 1], method= "pearson") ## About the same.

txDose <- txDoseOptim
cf_table <- cbind(cf_table[,c(1:28)], NearestHSF1= cfClosest*psClosest, txDoseNaive, txDose= txDoseOptim)
upreg <- cf_table[cf_table$V2 == 1 | cf_table$V3 == 1,]
plot(upreg$NearestHSF1+1, upreg$txDoseNaive+1, log="xy",pch=19, col=alpha("gray", 0.2))
points((upreg$NearestHSF1+1)[upreg$V2 == 1], (upreg$txDoseNaive+1)[upreg$V2 == 1], pch=19, col=alpha("blue", 0.2))

plot(upreg$NearestHSF1+1, upreg$txDose+1, log="xy", pch=19, col=alpha("gray", 0.2))
points((upreg$NearestHSF1+1)[upreg$V2 == 1], (upreg$txDose+1)[upreg$V2 == 1], pch=19, col=alpha("blue", 0.2))

###############################
## Evaluate.
boxplot(txDose[cf_table$V2 == 1], txDose[cf_table$V3 == 1], txDose[cf_table$V4 == 1], txDose[cf_table$V5 == 1], ylab="HSF-1 Dose")

dat <- melt(list(HSF1dep= txDose[cf_table$V2 == 1], Upindep= txDose[cf_table$V3 == 1], Down= txDose[cf_table$V4 == 1], Unreg= txDose[cf_table$V5 == 1]))
names(dat) <- c("HSF1Dose", "GeneType")

p<-ggplot(dat, aes(x=GeneType, y=HSF1Dose, fill=GeneType)) +
        #geom_violin(trim=FALSE) +
        geom_boxplot(width=0.1, fill="white") +
        labs(title="Pausing index by stage",x="Stage", y = "Pausing index [log-2]")#+
        #ylim(0,7)
p

## auPRC.
auc_pr(cf_table$V2 == 1,txDoseNaive)
auc_pr(cf_table$V2 == 1,txDoseOptim)
auc_pr(cf_table$V2 == 1,-cf_table$V12)
auc_pr(cf_table$V2 == 1,cfClosest*psClosest)
auc_pr(cf_table$V2 == 1,cf_table$V14*cf_table$V22)
auc_pr(cf_table$V2 == 1,cf_table$V15*cf_table$V23)
auc_pr(cf_table$V2 == 1,1/ (1 + exp(-0.326084*(cfClosest-81.28158)))* psClosest)

PRcurve(-cf_table$V12, cf_table[,2] == 1, ylim=c(0,1), xlim=c(0,1), col="gray")
par(new=TRUE)
PRcurve(txDoseNaive, cf_table[,2] == 1, ylim=c(0,1), xlim=c(0,1), col="blue")
par(new=TRUE)
PRcurve(txDoseOptim, cf_table[,2] == 1, ylim=c(0,1), xlim=c(0,1), col="red")

## auPRC on holdout.
auc_pr((cf_table$V2 == 1)[indx_tst],txDoseNaive[indx_tst])
auc_pr((cf_table$V2 == 1)[indx_tst],txDoseNaive[indx_tst])
auc_pr((cf_table$V2 == 1)[indx_tst],-cf_table$V12[indx_tst])

PRcurve(txDoseNaive[indx_tst], cf_table[indx_tst,2] == 1, ylim=c(0,1), xlim=c(0,1), col="black")

## Plot distance against dose: 
plot(txDose+1, shortestDistance+1, log="xy", ylab= "Distance (+1)", xlab= "Dose (+1)", pch=19, col=alpha("gray", 0.2))
points((txDose+1)[cf_table$V2 == 1], (shortestDistance+1)[cf_table$V2 == 1], pch=19, col=alpha("blue", 0.2))
points((txDose+1)[cf_table$V3 == 1], (shortestDistance+1)[cf_table$V3 == 1], pch=19, col=alpha("red", 0.2))
points((txDose+1)[cf_table$V4 == 1], (shortestDistance+1)[cf_table$V4 == 1], pch=19, col=alpha("green", 0.2))

## Violin plot for Dose among all of the groups.
boxplot(upreg$txDose[upreg$V2 == 1], upreg$txDose[upreg$V2 == 0])#, ylim=c(0,5000))
auc_pr(upreg$V2 == 1,upreg$txDose)
auc_pr(upreg$V2 == 1,-upreg$V12)
auc_pr(upreg$V2 == 1,1/ (1 + exp(-0.0452*(cfClosest[cf_table$V2 == 1 | cf_table$V3 == 1]-85)))* psClosest[cf_table$V2 == 1 | cf_table$V3 == 1])

PRcurve(-upreg$V12, upreg[,2] == 1, ylim=c(0,1), xlim=c(0,1), col="gray")
par(new=TRUE)
PRcurve(upreg$txDose, upreg[,2] == 1, ylim=c(0,1), xlim=c(0,1), col="red")

PRcurve(-cf_table$V12, cf_table[,2] == 1, ylim=c(0,1), xlim=c(0,1), col="gray")
par(new=TRUE)
PRcurve(txDoseNaive, cf_table[,2] == 1, ylim=c(0,1), xlim=c(0,1), col="blue")
par(new=TRUE)
PRcurve(txDoseOptim, cf_table[,2] == 1, ylim=c(0,1), xlim=c(0,1), col="red")

#####################
## Run a basic model
## 

# Put this code in a for loop and plot distribution
totalNum=NROW(upreg)
trainNum=ceiling(totalNum*0.90)

N <- 200
cf_list <- rep(0, N)
distance_list <- rep(0, N)
peak_strength_list <- rep(0, N)
dose_list <- rep(0, N)
distance_peak_strength_list <- rep(0, N)
cf_distance_peak_strength_list <- rep(0, N)
dose_cf_distance_peak_strength_list <- rep(0, N)

for (i in 1:N){
  #print(paste("i is", i))

  trainSample.idx=sample(1:NROW(upreg), trainNum)  #80% of indices to train on

  model_cf=glm(V2~V14, family="binomial", data=upreg[trainSample.idx,])
  scores_cf=predict(model_cf, upreg[-trainSample.idx,])

  model_distance=glm(V2~V12, family="binomial", data=upreg[trainSample.idx,])
  scores_distance=predict(model_distance, upreg[-trainSample.idx,])
  
  model_peak_strength=glm(V2~V22, family="binomial", data=upreg[trainSample.idx,])
  scores_peak_strength=predict(model_peak_strength, upreg[-trainSample.idx,])

  # Full dose model.
  model_dose=glm(V2~txDose, family="binomial", data=upreg[trainSample.idx,])
  scores_dose=predict(model_dose, upreg[-trainSample.idx,])

  # Combine distance * peak strength
  model_distance_peak_strength=glm(V2~V12+V22, family="binomial", data=upreg[trainSample.idx,]) 
  scores_distance_peak_strength=predict(model_distance_peak_strength, upreg[-trainSample.idx,])

  # Combine cf and distance and peak strength
  model_cf_distance_peak_strength=glm(V2~NearestHSF1+txDoseNaive, family="binomial", data=upreg[trainSample.idx,])
  scores_cf_distance_peak_strength=predict(model_cf_distance_peak_strength, upreg[-trainSample.idx,])

  # Combine cf and distance and peak strength and dose
  model_cf_distance_peak_strength_dose=glm(V2~txDose, family="binomial", data=upreg[trainSample.idx,])
  scores_cf_distance_peak_strength_dose=predict(model_cf_distance_peak_strength_dose, upreg[-trainSample.idx,])

  # Store results
  cf_list[[i]] <- auc_pr(upreg[-trainSample.idx,2],scores_cf)
  distance_list[[i]] <- auc_pr(upreg[-trainSample.idx,2],scores_distance)
  peak_strength_list[[i]] <- auc_pr(upreg[-trainSample.idx,2],scores_peak_strength)
  dose_list[[i]] <- auc_pr(upreg[-trainSample.idx,2],scores_dose)
  distance_peak_strength_list[[i]] <- auc_pr(upreg[-trainSample.idx,2],scores_distance_peak_strength)
  cf_distance_peak_strength_list[[i]] <- auc_pr(upreg[-trainSample.idx,2],scores_cf_distance_peak_strength)
  dose_cf_distance_peak_strength_list[[i]] <- auc_pr(upreg[-trainSample.idx,2],scores_cf_distance_peak_strength_dose)
}

boxplot(cf_list, distance_list, peak_strength_list, dose_list, distance_peak_strength_list, cf_distance_peak_strength_list, dose_cf_distance_peak_strength_list)

wilcox.test(dose_list, dose_cf_distance_peak_strength_list) ## NS
wilcox.test(dose_list, cf_distance_peak_strength_list) ## NS

