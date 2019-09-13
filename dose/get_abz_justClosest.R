# Given a list of genes, estimate a, b, and z.
#

## Example: 
## R --no-save --args training_and_test_data_w_geneNames.txt
##
## Note that this expects a column name: "is_train" where values of "True" will be used to train.

#####################
## Dependencies
## 
library(ROCR)
library(caTools)

auc_pr <- function(obs, pred) {
  xx.df <- prediction(pred, obs)
  perf  <- performance(xx.df, "prec", "rec")
  xy    <- data.frame(recall=perf@x.values[[1]], precision=perf@y.values[[1]])
  
  # take out division by 0 for lowest threshold
  xy <- subset(xy, !is.nan(xy$precision))
  
  res   <- trapz(xy$recall, xy$precision)
  res
}

args = commandArgs(trailingOnly=TRUE)

#####################
## Read in data.
##
fileName = "K562_gene_class_contact_frequency_3K_w_PROseq_counts_for_HSF1_peaks_table.txt"
cf_table = read.table(fileName, sep="\t")

## Get training/ test definitions that Paul is using.
tset <- read.table(args[1], header=TRUE, sep="\t") ## "training_and_test_data_w_geneNames.txt"
stopifnot(sum(as.character(cf_table$V1[cf_table$V2 == 1 | cf_table$V3 == 1]) == tset$X0)/NROW(tset) == 1)

indx <- which(cf_table$V2 == 1 | cf_table$V3 == 1)[tset$is_train == "True"] 

#############################################
## Optimize a and b by gradient descent.

## Starting parameters...
b <- 20
a <- 0.1
z <- 0

doseb <- function(x) {
 b <- x[1]

 cf_small <- cf_table[indx,]

 txDose <- sapply(1:NROW(cf_small), function(i) {
  txcf <- as.double(cf_small[i,14])
  txps <- as.double(cf_small[i,22])
  
  ntxcf <- as.double(cf_small[i,15])
  ntxps <- as.double(cf_small[i,23])

  sum(1/ (1 + exp(-a*(txcf-b)))* txps) + z*sum(1/ (1 + exp(-a*(ntxcf-b)))* ntxps)
 })

 1- auc_pr(cf_small$V2 ==1, txDose)
}


dosea <- function(x) {
 a <- x[1]

 cf_small <- cf_table[indx,]

 txDose <- sapply(1:NROW(cf_small), function(i) {
  txcf <- as.double(cf_small[i,14])
  txps <- as.double(cf_small[i,22])
  
  ntxcf <- as.double(cf_small[i,15])
  ntxps <- as.double(cf_small[i,23])

  sum(1/ (1 + exp(-a*(txcf-b)))* txps) + z*sum(1/ (1 + exp(-a*(ntxcf-b)))* ntxps)
 })

 1- auc_pr(cf_small$V2 ==1, txDose)
}

dosez <- function(x) {
 z <- x[1]

 cf_small <- cf_table[indx,]

 txDose <- sapply(1:NROW(cf_small), function(i) {
  txcf <- as.double(cf_small[i,14])
  txps <- as.double(cf_small[i,22])
  
  ntxcf <- as.double(cf_small[i,15])
  ntxps <- as.double(cf_small[i,23])

  sum(1/ (1 + exp(-a*(txcf-b)))* txps) + z*sum(1/ (1 + exp(-a*(ntxcf-b)))* ntxps)
 })

 1- auc_pr(cf_small$V2 ==1, txDose)
}

for(i in 1:3) {
 b <- optim(c(b), doseb, method= "Brent", lower= c(0), upper= c(500))$par 
 a <- optim(c(a), dosea, method= "Brent", lower= c(0), upper= c(1))$par 
 z <- optim(c(z), dosez, method= "Brent", lower= c(0), upper= c(1))$par 
}

#############################################
## Write out a, b, and z to a file.
write.table(list(b= b, a= a, z= z), "tmp.abz.out", row.names=FALSE, quote=FALSE, sep="\t")
