# Given a, b, c, and d, compute dose for each gene.
#

## Example: 
## R --no-save --args 81.2706904282568 0.326164077657388 81.2706904282568 0.326164077657388
##

#####################
## Dependencies
## 
args = commandArgs(trailingOnly=TRUE)

## Starting parameters...
b <- as.double(args[1])
a <- as.double(args[2])
d <- as.double(args[3])
c <- as.double(args[4])

#####################
## Read in data.
##
fileName = "K562_gene_class_contact_frequency_3K_w_PROseq_counts_for_HSF1_peaks_table.txt"
cf_table = read.table(fileName, sep="\t")

## Use only the up-reg genes.
indx <- which(cf_table$V2 == 1 | cf_table$V3 == 1)
cf_table <- cf_table[indx,]

#############################################
## Compute dose.

txDoseOptim <- sapply(1:NROW(cf_small), function(i) {
  cf <- c(as.double(strsplit(as.character(cf_small[i,20]), ",")[[1]]), as.double(strsplit(as.character(cf_small[i,21]), ",")[[1]]))
  ps <- c(as.double(strsplit(as.character(cf_small[i,24]), ",")[[1]]), as.double(strsplit(as.character(cf_small[i,25]), ",")[[1]]))
  tx <- c(as.double(strsplit(as.character(cf_small[i,57]), ",")[[1]]), as.double(strsplit(as.character(cf_small[i,59]), ",")[[1]]))

  sum(1/ (1 + exp(-a*(cf-b))) * (1 + exp(-c*(tx-d))) * ps)
 })
 
#############################################
## Write out a, b, and z to a file.
write.table(data.frame(geneID= cf_table$V1, Dose= txDoseOptim), "tmp.dose.out", row.names=FALSE, quote=FALSE, sep="\t")
