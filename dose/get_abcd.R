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
d <- 20
c <- 0.1

doseb <- function(x) {
 b <- x[1]

 cf_small <- cf_table[indx,]

 txDose <- sapply(1:NROW(cf_small), function(i) {
  cf <- c(as.double(strsplit(as.character(cf_small[i,20]), ",")[[1]]), as.double(strsplit(as.character(cf_small[i,21]), ",")[[1]]))
  ps <- c(as.double(strsplit(as.character(cf_small[i,24]), ",")[[1]]), as.double(strsplit(as.character(cf_small[i,25]), ",")[[1]]))
  tx <- c(as.double(strsplit(as.character(cf_small[i,57]), ",")[[1]]), as.double(strsplit(as.character(cf_small[i,59]), ",")[[1]]))

  sum(1/ (1 + exp(-a*(cf-b))) * (1 + exp(-c*(tx-d))) * ps)
 })

 1- auc_pr(cf_small$V2 ==1, txDose)
}


dosea <- function(x) {
 a <- x[1]

 cf_small <- cf_table[indx,]

 txDose <- sapply(1:NROW(cf_small), function(i) {
  cf <- c(as.double(strsplit(as.character(cf_small[i,20]), ",")[[1]]), as.double(strsplit(as.character(cf_small[i,21]), ",")[[1]]))
  ps <- c(as.double(strsplit(as.character(cf_small[i,24]), ",")[[1]]), as.double(strsplit(as.character(cf_small[i,25]), ",")[[1]]))
  tx <- c(as.double(strsplit(as.character(cf_small[i,57]), ",")[[1]]), as.double(strsplit(as.character(cf_small[i,59]), ",")[[1]]))

  sum(1/ (1 + exp(-a*(cf-b))) * (1 + exp(-c*(tx-d))) * ps)
 })

 1- auc_pr(cf_small$V2 ==1, txDose)
}


dosed <- function(x) {
 d <- x[1]

 cf_small <- cf_table[indx,]

 txDose <- sapply(1:NROW(cf_small), function(i) {
  cf <- c(as.double(strsplit(as.character(cf_small[i,20]), ",")[[1]]), as.double(strsplit(as.character(cf_small[i,21]), ",")[[1]]))
  ps <- c(as.double(strsplit(as.character(cf_small[i,24]), ",")[[1]]), as.double(strsplit(as.character(cf_small[i,25]), ",")[[1]]))
  tx <- c(as.double(strsplit(as.character(cf_small[i,57]), ",")[[1]]), as.double(strsplit(as.character(cf_small[i,59]), ",")[[1]]))

  sum(1/ (1 + exp(-a*(cf-b))) * (1 + exp(-c*(tx-d))) * ps)
 })

 1- auc_pr(cf_small$V2 ==1, txDose)
}


dosec <- function(x) {
 c <- x[1]

 cf_small <- cf_table[indx,]

 txDose <- sapply(1:NROW(cf_small), function(i) {
  cf <- c(as.double(strsplit(as.character(cf_small[i,20]), ",")[[1]]), as.double(strsplit(as.character(cf_small[i,21]), ",")[[1]]))
  ps <- c(as.double(strsplit(as.character(cf_small[i,24]), ",")[[1]]), as.double(strsplit(as.character(cf_small[i,25]), ",")[[1]]))
  tx <- c(as.double(strsplit(as.character(cf_small[i,57]), ",")[[1]]), as.double(strsplit(as.character(cf_small[i,59]), ",")[[1]]))

  sum(1/ (1 + exp(-a*(cf-b))) * (1 + exp(-c*(tx-d))) * ps)
 })

 1- auc_pr(cf_small$V2 ==1, txDose)
}


for(i in 1:3) {
 b <- optim(c(b), doseb, method= "Brent", lower= c(0), upper= c(500))$par 
 a <- optim(c(a), dosea, method= "Brent", lower= c(0), upper= c(1))$par 
 d <- optim(c(d), dosed, method= "Brent", lower= c(0), upper= c(1500))$par 
 c <- optim(c(c), dosec, method= "Brent", lower= c(0), upper= c(1))$par 
}

#############################################
## Write out a, b, and z to a file.
write.table(list(b= b, a= a, d= d, c= c), "tmp.abcd.out", row.names=FALSE, quote=FALSE, sep="\t")
