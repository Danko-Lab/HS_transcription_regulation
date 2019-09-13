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
require(beeswarm)

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
cf_table <- cbind(cf_table[,c(1:52)], NearestHSF1= cfClosest*psClosest, txDoseNaive, txDose= txDoseOptim)
cf_table <- cbind(cf_table[,c(1:52)], NearestHSF1= cfClosest*psClosest, txDoseNaive, txDose= txDoseOptim)


################################################################################
## Compare total dose quantity across distinct classes of HS regulated genes.
#boxplot(txDose[cf_table$V2 == 1], txDose[cf_table$V3 == 1], txDose[cf_table$V4 == 1], txDose[cf_table$V5 == 1], ylab="HSF-1 Dose")
boxplot(txDoseNaive[cf_table$V2 == 1], txDoseNaive[cf_table$V3 == 1], txDoseNaive[cf_table$V4 == 1], txDoseNaive[cf_table$V5 == 1], ylab="HSF-1 Dose")

dat <- melt(list(HSF1dep= txDose[cf_table$V2 == 1], Upindep= txDose[cf_table$V3 == 1], Down= txDose[cf_table$V4 == 1], Unreg= txDose[cf_table$V5 == 1]))
names(dat) <- c("HSF1Dose", "GeneType")

p<-ggplot(dat, aes(x=GeneType, y=HSF1Dose, fill=GeneType)) +
        geom_violin(trim=TRUE, scale="width") +
        geom_boxplot(width=0.1, fill="white") +
        labs(title="HSF1 Dose",x="Stage", y = "Total HSF1 dose.")#+
        #ylim(0,7)
p

##############################################################################
## Do multiple HSF-1 binding sites contribute substantially to dose?

## Plot the number of TFBSs near HSF1 depdendent up-regulated genes.
N_HSF1_bs <- sapply(1:NROW(cf_table), function(i) {
 txcf <- as.double(strsplit(as.character(cf_table[i,20]), ",")[[1]])
 ntxcf <- as.double(strsplit(as.character(cf_table[i,21]), ",")[[1]])

 NROW(ntxcf)
 #NROW(txcf) + NROW(ntxcf)
})

boxplot(N_HSF1_bs[cf_table$V2 == 1], N_HSF1_bs[cf_table$V3 == 1], N_HSF1_bs[cf_table$V4 == 1], N_HSF1_bs[cf_table$V5 == 1], ylab="# HSF-1 binding sites", ylim=c(0,10))

## YES! And this is especially true for transcribed HSF-1 binding sites.

dosen <- list(integer(0), integer(0),integer(0),integer(0),integer(0),integer(0),integer(0),integer(0),integer(0),integer(0))
for(i in which(cf_table$V2 == 1)) {
 txcf <- as.double(strsplit(as.character(cf_table[i,20]), ",")[[1]])
 txps <- as.double(strsplit(as.character(cf_table[i,24]), ",")[[1]])

 ntxcf <- as.double(strsplit(as.character(cf_table[i,21]), ",")[[1]])
 ntxps <- as.double(strsplit(as.character(cf_table[i,25]), ",")[[1]])
 
 doses <- rev(sort(c(txcf*txps, ntxcf*ntxps)))
 if(NROW(doses) > 0 & sum(doses) > 0) {
 for(j in 1:min(10, NROW(doses))) {
   dosen[[j]] <- c(dosen[[j]], doses[j]/ sum(doses))
 }
 }
}

## Now plot a beeswarm plot of this list...
pdf("Dose_n.pdf")
  beeswarm(dosen, corral="random")
dev.off()

## Is the closest typically the strongest?
dose_closest <- integer(0)
for(i in which(cf_table$V2 == 1)) {
 txcf <- as.double(strsplit(as.character(cf_table[i,20]), ",")[[1]])
 txps <- as.double(strsplit(as.character(cf_table[i,24]), ",")[[1]])

 ntxcf <- as.double(strsplit(as.character(cf_table[i,21]), ",")[[1]])
 ntxps <- as.double(strsplit(as.character(cf_table[i,25]), ",")[[1]])
 
 doses <- c(txcf*txps, ntxcf*ntxps)
 dists <- c(as.double(strsplit(as.character(cf_table[i,26]), ",")[[1]]), as.double(strsplit(as.character(cf_table[i,27]), ",")[[1]]))
 doses <- doses[order(dists)]

 if(NROW(doses) > 0 & sum(doses) > 0) {
  dose_closest <- c(dose_closest, which.max(doses))
 }
}

## Now plot a beeswarm plot.
summary(as.factor(dose_closest))
summary(as.factor(dose_closest))/NROW(dose_closest)

#############################################################
## Find genes where multiple sites contribute to the total.
plot(cf_table$NearestHSF1+1, cf_table$txDoseNaive+1, log="xy",pch=19, xlab= "Nearest HSF1 dose", ylab="Total HSF1 dose", col=alpha("gray", 0.2))
points((cf_table$NearestHSF1+1)[cf_table$V2 == 1], (cf_table$txDoseNaive+1)[cf_table$V2 == 1], pch=19, col=alpha("blue", 0.2))

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
##
## If so, we expect that a glm based on points where a primary site contributes will 
## under-estimate fold-changes.

#cf_table <- cf_table[!((txDose_- NearestHSF1_)<0) & cf_table$V9 > 0,] ## for now...

NearestHSF1_ <- cf_table$NearestHSF1 #cf_table$V14*cf_table$V22# #1/ (1 + exp(-0.03479813*(cfClosest-35.46118269)))* psClosest #cfClosest*psClosest
txDose_      <- cf_table$txDoseNaive

indx <- cf_table$V2 == 1

cor.test(cf_table$txDose[indx], cf_table$V9[indx], method="pearson") ## Optimum dose.
cor.test(txDose_[indx], cf_table$V9[indx], method="pearson") ## Dose naive. 
cor.test(NearestHSF1_[indx], cf_table$V9[indx], method="pearson") ## Dose for nearest HSF1 site. 
cor.test(cf_table$V22[indx], cf_table$V9[indx], method="pearson") ## HSF1 binding quantity for HSF1 site. 

cor.test(cf_table$txDose[indx], cf_table$V9[indx], method="spearman")
cor.test(txDose_[indx], cf_table$V9[indx], method="spearman")
cor.test(NearestHSF1_[indx], cf_table$V9[indx], method="spearman")
cor.test(cf_table$V22[indx], cf_table$V9[indx], method="spearman")

plot(cf_table$V22[indx], cf_table$V9[indx], method="spearman")

## If multiple sites are contributing, residuals for genes w/ large contributions of other sites should have residuals that are systematically negative in  glm(FC ~ Nearest).
model <- glm(cf_table$V9[indx]~NearestHSF1_[indx])
summary(model$residuals)
summary(model$residuals[((txDose_ - NearestHSF1_) == 0)[indx]])
summary(model$residuals[((txDose_ - NearestHSF1_) > 0)[indx]])

wilcox.test(model$residuals[((txDose_- NearestHSF1_) == 0)[indx]], model$residuals[((txDose_- NearestHSF1_) > 0)[indx]])

cor.test(model$residuals, ((txDose_- NearestHSF1_)/ (NearestHSF1_+1))[indx], method="spearman")
cor.test(model$residuals/(cf_table$V9[indx]), ((txDose_- NearestHSF1_)/ (NearestHSF1_+1))[indx], method="spearman")

## If we estimate using genes w/ only 1 HSF1 site contributing, we should underestimate residuals for genes where multimple HSF1 sites contribute.
model <- glm(V9~txDoseNaive, data= cf_table[indx & (txDose_- NearestHSF1_) == 0,])
residuals_ <- cf_table$V9[indx] - predict(model, cf_table[indx,])

summary(residuals_)
summary(residuals_[((txDose_- NearestHSF1_) == 0)[indx]])
summary(residuals_[((txDose_- NearestHSF1_) > 0)[indx]])

wilcox.test(residuals_[((txDose_- NearestHSF1_) == 0)[indx]], residuals_[((txDose_- NearestHSF1_) > 0)[indx]])
cor.test(residuals_, (txDose_- NearestHSF1_)[indx], method="pearson")
cor.test(residuals_, (txDose_- NearestHSF1_)[indx], method="spearman")

plot(V9~txDoseNaive, data=cf_table[indx,], xlab="HSF1 dose", ylab="Log-2 fold-change, HS")
abline(model)
points(V9~txDoseNaive, data=cf_table[indx & (txDose_- NearestHSF1_) > 0,], col="red", pch=19)
points(cf_table$V9[indx]~txDose_[indx], col="blue", pch=19)

summary(txDose_- NearestHSF1_) ## Sanity check...
plot(txDose_[indx]~NearestHSF1_[indx])

#############################################
## Repeat multiple site analysis using pausing indices.
##
## HSF1 releases Pol II into productive elongation.

# V41: 'body_rpkm_HS_br1'
# V42: 'body_rpkm_HS_br2'
# V43: 'body_rpkm_NHS_br1'
# V44: 'body_rpkm_NHS_br2'
# V45: 'pause_rpkm_HS_br1'
# V46: 'pause_rpkm_HS_br2'
# V47: 'pause_rpkm_NHS_br1'
# V48: 'pause_rpkm_NHS_br2'

indx <- cf_table$V2 == 1

PIhs <- rowMeans(cf_table[,c(45:46)])/rowMeans(cf_table[,c(41:42)])
PInhs<- rowMeans(cf_table[,c(47:48)])/rowMeans(cf_table[,c(43:44)])

PI <- PInhs/PIhs; cf_table$PI <- PI
FC <- rowMeans(cf_table[,c(45:46)])/ rowMeans(cf_table[,c(47:48)])
FCb<- rowMeans(cf_table[,c(41:42)])/ rowMeans(cf_table[,c(43:44)])

cor.test(txDose_[indx], PI[indx], method="spearman")
cor.test(NearestHSF1_[indx], PI[indx], method="spearman")

cor.test(txDose_[indx], FC[indx], method="spearman")
cor.test(NearestHSF1_[indx], FC[indx], method="spearman")

cor.test(txDose_[indx], FCb[indx], method="spearman")
cor.test(NearestHSF1_[indx], FCb[indx], method="spearman")

plot(PI~txDoseNaive, data=cf_table[indx,], xlab="HSF1 Dose", ylab="Fold change in pausing index", log="y")
plot(txDose_[indx], FCb[indx], log="y")

#############################################
## Do we see a higher contact frequency at distal sites than expected by chance?!
##
## If so, we expect that a glm based on points where a primary site contributes will 
## under-estimate fold-changes.

# V12: Distance to nearest transcribed HSF peak
# V13: Distance to nearest non-transcribed HSF peak
# V14: CF to nearest transcribed
# V15: CF to nearest non-transcribed
# V16: CF to before bin transcribed
# V17: CF to before bin non-transcribed
# V18: CF to after bin transcribed
# V19: CF to after bin non-transcribed

## Are points closer than expected by chance?
plot((cf_table$V14+1), (cf_table$V12+1), log="xy", ylab="Linear distance", xlab="contact frequency", col=alpha("gray", 0.5))
points((cf_table$V14+1)[indx], (cf_table$V12+1)[indx], pch=19, col=alpha("red", 0.5))

## Are points closer than the surrounding windows for HSF1 dependent genes?
indx <- cf_table$V12 > 5000 & cf_table$V2 == 1

cf_diff <- c(cf_table$V14-rowMeans(cf_table[,c(16,18)]))
boxplot(cf_diff[indx], cf_diff[cf_table$V2 == 0], names=c("HSF1 dep", "All"), ylab="CF center - avg(CF flank)")
wilcox.test(cf_diff[indx])
