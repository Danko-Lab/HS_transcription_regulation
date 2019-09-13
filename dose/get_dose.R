# Given a, b, and z, compute dose for each gene.
#

## Example: 
## R --no-save --args 81.2706904282568 0.326164077657388 0.612254252326722
##

#####################
## Dependencies
## 
args = commandArgs(trailingOnly=TRUE)

## Starting parameters...
b <- as.double(args[1])
a <- as.double(args[2])
z <- as.double(args[3])

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

txDoseOptim <- sapply(1:NROW(cf_table), function(i) {

 txcf <- as.double(strsplit(as.character(cf_table[i,20]), ",")[[1]])
 txps <- as.double(strsplit(as.character(cf_table[i,24]), ",")[[1]])
 
 ntxcf <- as.double(strsplit(as.character(cf_table[i,21]), ",")[[1]])
 ntxps <- as.double(strsplit(as.character(cf_table[i,25]), ",")[[1]])
 
 sum((1/ (1 + exp(-a*(txcf-b)))* txps))+z*sum(1/ (1 + exp(-a*(ntxcf-b)))* ntxps)
})

#############################################
## Write out a, b, and z to a file.
write.table(data.frame(geneID= cf_table$V1, Dose= txDoseOptim), "tmp.dose.out", row.names=FALSE, quote=FALSE, sep="\t")
