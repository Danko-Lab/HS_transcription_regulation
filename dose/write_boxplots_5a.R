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
# V53: Count (plus strand + minus strand) for closest transcribed peak (HS)
# V54: Count for closest transcribed peak (NHS)
# V55: Count for closest non-transcribed peak (HS)
# V56: Count for closest non-transcribed peak (NHS)
# V57: Comma separated list of counts for all transcribed peaks (HS)
# V58: Comma separated list of counts for all transcribed peaks (NHS)
# V59: Comma separated list of counts for all non-transcribed peaks (HS)
# V60: Comma separated list of counts for all non-transcribed peaks (NHS)


#####################
## Dependencies
## 
library(DMwR)
library(ROCR)
library(caTools)
require(ggplot2)
library(reshape2)
library(gridExtra)

#####################
## Read in data.
##
fileName = "K562_gene_class_contact_frequency_3K_w_PROseq_counts_for_HSF1_peaks_table_HC.txt" # "K562_gene_class_contact_frequency_3K_w_PROseq_counts_and_rpkm_minus_HS_table.txt"

cf_table = read.table(fileName, sep="\t")
head(cf_table)

#####################
## Define closest
##
indx <- as.double(cf_table[,12] > cf_table[,13])

distClosest <- sapply(1:NROW(cf_table), function(i) {cf_table[i,12+indx[i]]})
cfClosest   <- sapply(1:NROW(cf_table), function(i) {cf_table[i,14+indx[i]]})
psClosest   <- sapply(1:NROW(cf_table), function(i) {cf_table[i,22+indx[i]]})
txHSClosest   <- sapply(1:NROW(cf_table), function(i) {cf_table[i,53+(2*indx[i])]})
txNHSClosest   <- sapply(1:NROW(cf_table), function(i) {cf_table[i,54+(2*indx[i])]})


#############################################
## How many sites found within the distance range...
nbs <- sapply(1:NROW(cf_table), function(i) {
 txcf <- as.double(strsplit(as.character(cf_table[i,20]), ",")[[1]])
 txps <- as.double(strsplit(as.character(cf_table[i,24]), ",")[[1]])

 NROW(txcf)+NROW(txps)
})

#############################################
## Now plot all of these...

datcf <- melt(list("HSF1-dep"= cfClosest[cf_table$V2 == 1], 
				"Up Reg"= cfClosest[cf_table$V3 == 1], 
				"Down Reg"= cfClosest[cf_table$V4 == 1], 
				"Unreg"= cfClosest[cf_table$V5 == 1]))		
names(datcf) <- c("ContactFrequency", "GeneType")
datcf$GeneType <- factor(datcf$GeneType, levels= c("HSF1-dep", "Up Reg", "Down Reg", "Unreg"), ordered=TRUE)


datdist <- melt(list("HSF1-dep"= distClosest[cf_table$V2 == 1], 
				"Up Reg"= distClosest[cf_table$V3 == 1], 
				"Down Reg"= distClosest[cf_table$V4 == 1], 
				"Unreg"= distClosest[cf_table$V5 == 1]))		
names(datdist) <- c("Distance", "GeneType")
datdist$GeneType <- factor(datdist$GeneType, levels= c("HSF1-dep", "Up Reg", "Down Reg", "Unreg"), ordered=TRUE)

datps <- melt(list("HSF1-dep"= psClosest[cf_table$V2 == 1], 
				"Up Reg"= psClosest[cf_table$V3 == 1], 
				"Down Reg"= psClosest[cf_table$V4 == 1], 
				"Unreg"= psClosest[cf_table$V5 == 1]))		
names(datps) <- c("PeakStrnegth", "GeneType")
datps$GeneType <- factor(datps$GeneType, levels= c("HSF1-dep", "Up Reg", "Down Reg", "Unreg"), ordered=TRUE)

datnbs <- melt(list("HSF1-dep"= nbs[cf_table$V2 == 1], 
				"Up Reg"= nbs[cf_table$V3 == 1], 
				"Down Reg"= nbs[cf_table$V4 == 1], 
				"Unreg"= nbs[cf_table$V5 == 1]))		
names(datnbs) <- c("NumHSF1BS", "GeneType")
datnbs$GeneType <- factor(datnbs$GeneType, levels= c("HSF1-dep", "Up Reg", "Down Reg", "Unreg"), ordered=TRUE)

dattxhs <- melt(list("HSF1-dep"= txHSClosest[cf_table$V2 == 1], 
				"Up Reg"= txHSClosest[cf_table$V3 == 1], 
				"Down Reg"= txHSClosest[cf_table$V4 == 1], 
				"Unreg"= txHSClosest[cf_table$V5 == 1]))		
names(dattxhs) <- c("TxHS", "GeneType")
dattxhs$GeneType <- factor(dattxhs$GeneType, levels= c("HSF1-dep", "Up Reg", "Down Reg", "Unreg"), ordered=TRUE)

dattxnhs <- melt(list("HSF1-dep"= txNHSClosest[cf_table$V2 == 1], 
				"Up Reg"= txNHSClosest[cf_table$V3 == 1], 
				"Down Reg"= txNHSClosest[cf_table$V4 == 1], 
				"Unreg"= txNHSClosest[cf_table$V5 == 1]))		
names(dattxnhs) <- c("TxNHS", "GeneType")
dattxnhs$GeneType <- factor(dattxnhs$GeneType, levels= c("HSF1-dep", "Up Reg", "Down Reg", "Unreg"), ordered=TRUE)


theme_set(theme_classic(base_size = 24)) 
cols = c("#8b0000", "#ff0000", "#0000ff", "#ffa500")

pcf<-ggplot(datcf, aes(x=GeneType, y=ContactFrequency, fill=GeneType)) +
        geom_violin(trim=FALSE, scale="width") +
        geom_boxplot(width=0.1, fill="white", outlier.size=0.5) +
        labs(x="Gene class", y = "Contact frequency to HSF1")+
        ylim(0,1500)+ guides(fill=FALSE)+ theme_classic()+ scale_fill_manual(values=cols)+ scale_colour_manual(values=cols)

pdist<-ggplot(datdist, aes(x=GeneType, y=Distance, fill=GeneType)) +
        geom_violin(trim=FALSE, scale="width") +
        geom_boxplot(width=0.1, fill="white", outlier.size=0.5) +
        labs(x="Gene class", y = "Distance to HSF1")+
        ylim(0,250000)+ guides(fill=FALSE)+ theme_classic()+ scale_fill_manual(values=cols)+ scale_colour_manual(values=cols)

pps <-ggplot(datps, aes(x=GeneType, y=PeakStrnegth, fill=GeneType)) +
        geom_violin(trim=FALSE, scale="width") +
        geom_boxplot(width=0.1, fill="white", outlier.size=0.5) +
        labs(x="Gene class", y = "HSF1 peak strength")+
        ylim(0,35)+ guides(fill=FALSE)+ theme_classic()+ scale_fill_manual(values=cols)+ scale_colour_manual(values=cols)
	
pnbs <-ggplot(datnbs, aes(x=GeneType, y=NumHSF1BS, fill=GeneType)) +
        geom_violin(trim=FALSE, scale="width") +
        geom_boxplot(width=0.1, fill="white", outlier.size=0.5) +
        labs(x="Gene class", y = "# HSF1 peaks")+
        ylim(0,35)+ guides(fill=FALSE)+ theme_classic()+ scale_fill_manual(values=cols)+ scale_colour_manual(values=cols)

pntxhs <-ggplot(dattxhs, aes(x=GeneType, y=TxHS, fill=GeneType)) +
        geom_violin(trim=FALSE, scale="width") +
        geom_boxplot(width=0.1, fill="white", outlier.size=0.5) +
        labs(x="Gene class", y = "Txn HS at HSF1 peak")+
        ylim(0,1500)+ guides(fill=FALSE)+ theme_classic()+ scale_fill_manual(values=cols)+ scale_colour_manual(values=cols)

pntxnhs <-ggplot(dattxnhs, aes(x=GeneType, y=TxNHS, fill=GeneType)) +
        geom_violin(trim=FALSE, scale="width") +
        geom_boxplot(width=0.1, fill="white", outlier.size=0.5) +
        labs(x="Gene class", y = "Txn HS at HSF1 peak")+
        ylim(0,1500)+ guides(fill=FALSE)+ theme_classic()+ scale_fill_manual(values=cols)+ scale_colour_manual(values=cols)

		
pdf("Fig5a.pdf")
 grid.arrange(pcf, pdist, pnbs, pps, pntxhs, pntxnhs, nrow = 2)
dev.off()
		
pdf("FullPDF.pdf")
 pcf
 pdist
 pps
 pnbs
dev.off()


