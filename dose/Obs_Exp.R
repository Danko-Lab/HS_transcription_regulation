## 
## 1. Chromosome, 
## 2. HSF1_binding_site, 
## 3. Chromosome, 
## 4. TSS_being_tested, 
## 5. Distance_from_HSF1_binding_site_to_TSS, 
## 6. HS_Test_p-value, 
## 7. HS_Observed_count, 
## 8. HS_Expected_count, 
## 9. HS_adjusted_p-value, 
## 10.Chromosome, 
## 11.HSF1_binding_site, 
## 12.Chromosome, 
## 13.TSS_being_tested, 
## 14.Distance_from_HSF1_binding_site_to_TSS, 
## 15.NHS_Test_p-value, 
## 16.NHS_Observed_count, 
## 17.NHS_Expected_count, 
## 18.NHS_adjusted_p-value
## 

cc <- read.table("Pasted.K562_HA_NHS_combined.chr1-22.cc")

HS_oe <- cc[,7]/ cc[,8]
NHS_oe <- cc[,16]/ cc[,17]

## Compare HS to NHS observed / expected (all pairs).
boxplot(HS_oe, NHS_oe, ylim=c(0,10))
abline(h=1, lty="dotted")

wilcox.test(HS_oe, NHS_oe)

boxplot(HS_oe-NHS_oe, ylim=c(-10,10)) ## NOTE: Expecting >0 if observed:expected contact ratio is gained.
summary(HS_oe-NHS_oe)
abline(h=0, lty="dotted")

## Compare to HS and NHS as a function of distance.
EP_dist <- abs(cc[,2] - cc[,4])
plot(EP_dist, HS_oe-NHS_oe)
abline(h=0, lty="dotted", col="dark red")

plot(EP_dist, HS_oe-NHS_oe, xlim=c(0,5000))
abline(h=0, lty="dotted", col="dark red")

plot(EP_dist, HS_oe-NHS_oe, ylim=c(-5,5))
abline(h=0, lty="dotted", col="dark red")

## Condition on significance in James' model.
which(cond <- cc[,9] < 0.15)
EP_dist <- abs(cc[,2] - cc[,4])
plot(EP_dist[cond], (HS_oe-NHS_oe)[cond], ylim=c(-5,5))
abline(h=0, lty="dotted", col="dark red")
plot(EP_dist[cond], (HS_oe-NHS_oe)[cond], xlim=c(0,5000), ylim=c(-5,5))
abline(h=0, lty="dotted", col="dark red")


which(cond <- cc[,18] < 0.15)
EP_dist <- abs(cc[,2] - cc[,4])
plot(EP_dist[cond], (HS_oe-NHS_oe)[cond], ylim=c(-5,5))
abline(h=0, lty="dotted", col="dark red")
plot(EP_dist[cond], (HS_oe-NHS_oe)[cond], xlim=c(0,15000), ylim=c(-5,5))
abline(h=0, lty="dotted", col="dark red")








