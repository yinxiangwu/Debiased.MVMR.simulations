# Summarize simulation results for adIVW, 10, 20, 30 exposures
# created by Yinxiang Wu on Jan 02, 2023
# last update: Dec 30, 2023

res.tmp1 <- NULL
K=10 # number of exposures
k=1 # index of beta to extract
for (i in 1:100) {
  tryCatch({
    res.tmp1 <- rbind(res.tmp1, read.csv(paste0('Res/for_pub/10_20_30_exp_example/storeallresults/K',K,'/res_dmvmr_beta',k,'_K',K,'_job',i,'.csv')))
  }, error=function(e){})
}

# IVW
res.ivw1 <- mean(res.tmp1[,1])
res.ivw.sd1 <- sd(res.tmp1[,1])
res.ivw.se1 <- mean(res.tmp1[,2])
res.ivw.cp1 <- mean(res.tmp1[,1] + qnorm(0.025) * res.tmp1[,2] < 0.4 & 0.4 < res.tmp1[,1] + qnorm(0.975) * res.tmp1[,2])

# MR-Egger
res.egger1 <- mean(res.tmp1[,3])
res.egger.sd1 <- sd(res.tmp1[,3])
res.egger.se1 <- mean(res.tmp1[,4])
res.egger.cp1 <- mean(res.tmp1[,3] + qnorm(0.025) * res.tmp1[,4] < 0.4 & 0.4 < res.tmp1[,3] + qnorm(0.975) * res.tmp1[,4])

# MR-Median
res.median1 <- mean(res.tmp1[,5])
res.median.sd1 <- sd(res.tmp1[,5])
res.median.se1 <- mean(res.tmp1[,6])
res.median.cp1 <- mean(res.tmp1[,5] + qnorm(0.025) * res.tmp1[,6] < 0.4 & 0.4 < res.tmp1[,5] + qnorm(0.975) * res.tmp1[,6])

# Grapple
res.grapple1 <- mean(res.tmp1[,7])
res.grapple.sd1 <- sd(res.tmp1[,7])
res.grapple.se1 <- mean(res.tmp1[,8])
res.grapple.cp1 <- mean(res.tmp1[,7] + qnorm(0.025) * res.tmp1[,8] < 0.4 & 0.4 < res.tmp1[,7] + qnorm(0.975) * res.tmp1[,8])

# MRBEE
res.mrbee1 <- mean(res.tmp1[,9])
res.mrbee.sd1 <- sd(res.tmp1[,9])
res.mrbee.se1 <- mean(res.tmp1[,10])
res.mrbee.cp1 <- mean(res.tmp1[,9] + qnorm(0.025) * res.tmp1[,10] < 0.4 & 0.4 < res.tmp1[,9] + qnorm(0.975) * res.tmp1[,10])

# dIVW
res.divw1 <- mean(res.tmp1[,11])
res.divw.sd1 <- sd(res.tmp1[,11])
res.divw.se1 <- mean(res.tmp1[,12])
res.divw.cp1 <- mean(res.tmp1[,11] + qnorm(0.025) * res.tmp1[,12] < 0.4 & 0.4 < res.tmp1[,11] + qnorm(0.975) * res.tmp1[,12])

# adIVW
res.adivw1 <- mean(res.tmp1[,13])
res.adivw.sd1 <- sd(res.tmp1[,13])
res.adivw.se1 <- mean(res.tmp1[,14])
res.adivw.cp1 <- mean(res.tmp1[,13] + qnorm(0.025) * res.tmp1[,14] < 0.4 & 0.4 < res.tmp1[,13] + qnorm(0.975) * res.tmp1[,14])

res.tmp2 <- NULL
K=20
for (i in 1:100) {
  tryCatch({
    res.tmp2 <- rbind(res.tmp2, read.csv(paste0('Res/for_pub/10_20_30_exp_example/storeallresults/K',K,'/res_dmvmr_beta',k,'_K',K,'_job',i,'.csv')))
  }, error=function(e){})
}

# IVW
res.ivw2 <- mean(res.tmp2[,1])
res.ivw.sd2 <- sd(res.tmp2[,1])
res.ivw.se2 <- mean(res.tmp2[,2])
res.ivw.cp2 <- mean(res.tmp2[,1] + qnorm(0.025) * res.tmp2[,2] < 0.4 & 0.4 < res.tmp2[,1] + qnorm(0.975) * res.tmp2[,2])

# MR-Egger
res.egger2 <- mean(res.tmp2[,3])
res.egger.sd2 <- sd(res.tmp2[,3])
res.egger.se2 <- mean(res.tmp2[,4])
res.egger.cp2 <- mean(res.tmp2[,3] + qnorm(0.025) * res.tmp2[,4] < 0.4 & 0.4 < res.tmp2[,3] + qnorm(0.975) * res.tmp2[,4])

# MR-Median
res.median2 <- mean(res.tmp2[,5])
res.median.sd2 <- sd(res.tmp2[,5])
res.median.se2 <- mean(res.tmp2[,6])
res.median.cp2 <- mean(res.tmp2[,5] + qnorm(0.025) * res.tmp2[,6] < 0.4 & 0.4 < res.tmp2[,5] + qnorm(0.975) * res.tmp2[,6])

# Grapple
res.grapple2 <- mean(res.tmp2[,7])
res.grapple.sd2 <- sd(res.tmp2[,7])
res.grapple.se2 <- mean(res.tmp2[,8])
res.grapple.cp2 <- mean(res.tmp2[,7] + qnorm(0.025) * res.tmp2[,8] < 0.4 & 0.4 < res.tmp2[,7] + qnorm(0.975) * res.tmp2[,8])

# MRBEE
res.mrbee2 <- mean(res.tmp2[,9])
res.mrbee.sd2 <- sd(res.tmp2[,9])
res.mrbee.se2 <- mean(res.tmp2[,10])
res.mrbee.cp2 <- mean(res.tmp2[,9] + qnorm(0.025) * res.tmp2[,10] < 0.4 & 0.4 < res.tmp2[,9] + qnorm(0.975) * res.tmp2[,10])

# dIVW
res.divw2 <- mean(res.tmp2[,11])
res.divw.sd2 <- sd(res.tmp2[,11])
res.divw.se2 <- mean(res.tmp2[,12])
res.divw.cp2 <- mean(res.tmp2[,11] + qnorm(0.025) * res.tmp2[,12] < 0.4 & 0.4 < res.tmp2[,11] + qnorm(0.975) * res.tmp2[,12])

# adIVW
res.adivw2 <- mean(res.tmp2[,13])
res.adivw.sd2 <- sd(res.tmp2[,13])
res.adivw.se2 <- mean(res.tmp2[,14])
res.adivw.cp2 <- mean(res.tmp2[,13] + qnorm(0.025) * res.tmp2[,14] < 0.4 & 0.4 < res.tmp2[,13] + qnorm(0.975) * res.tmp2[,14])


res.tmp3 <- NULL
K=30
for (i in 1:100) {
  tryCatch({
    res.tmp3 <- rbind(res.tmp3, read.csv(paste0('Res/for_pub/10_20_30_exp_example/storeallresults/K',K,'/res_dmvmr_beta',k,'_K',K,'_job',i,'.csv')))
  }, error=function(e){})
}

# IVW
res.ivw3 <- mean(res.tmp3[,1])
res.ivw.sd3 <- sd(res.tmp3[,1])
res.ivw.se3 <- mean(res.tmp3[,2])
res.ivw.cp3 <- mean(res.tmp3[,1] + qnorm(0.025) * res.tmp3[,2] < 0.4 & 0.4 < res.tmp3[,1] + qnorm(0.975) * res.tmp3[,2])

# MR-Egger
res.egger3 <- mean(res.tmp3[,3])
res.egger.sd3 <- sd(res.tmp3[,3])
res.egger.se3 <- mean(res.tmp3[,4])
res.egger.cp3 <- mean(res.tmp3[,3] + qnorm(0.025) * res.tmp3[,4] < 0.4 & 0.4 < res.tmp3[,3] + qnorm(0.975) * res.tmp3[,4])

# MR-Median
res.median3 <- mean(res.tmp3[,5])
res.median.sd3 <- sd(res.tmp3[,5])
res.median.se3 <- mean(res.tmp3[,6])
res.median.cp3 <- mean(res.tmp3[,5] + qnorm(0.025) * res.tmp3[,6] < 0.4 & 0.4 < res.tmp3[,5] + qnorm(0.975) * res.tmp3[,6])


# Grapple
res.grapple3 <- mean(res.tmp3[,7])
res.grapple.sd3 <- sd(res.tmp3[,7])
res.grapple.se3 <- mean(res.tmp3[,8])
res.grapple.cp3 <- mean(res.tmp3[,7] + qnorm(0.025) * res.tmp3[,8] < 0.4 & 0.4 < res.tmp3[,7] + qnorm(0.975) * res.tmp3[,8])

# MRBEE
res.mrbee3 <- mean(res.tmp3[,9])
res.mrbee.sd3 <- sd(res.tmp3[,9])
res.mrbee.se3 <- mean(res.tmp3[,10])
res.mrbee.cp3 <- mean(res.tmp3[,9] + qnorm(0.025) * res.tmp3[,10] < 0.4 & 0.4 < res.tmp3[,9] + qnorm(0.975) * res.tmp3[,10])

# dIVW
res.divw3 <- mean(res.tmp3[,11])
res.divw.sd3 <- sd(res.tmp3[,11])
res.divw.se3 <- mean(res.tmp3[,12])
res.divw.cp3 <- mean(res.tmp3[,11] + qnorm(0.025) * res.tmp3[,12] < 0.4 & 0.4 < res.tmp3[,11] + qnorm(0.975) * res.tmp3[,12])

# adIVW
res.adivw3 <- mean(res.tmp3[,13])
res.adivw.sd3 <- sd(res.tmp3[,13])
res.adivw.se3 <- mean(res.tmp3[,14])
res.adivw.cp3 <- mean(res.tmp3[,13] + qnorm(0.025) * res.tmp3[,14] < 0.4 & 0.4 < res.tmp3[,13] + qnorm(0.975) * res.tmp3[,14])

res.summary <- cbind(rbind(c(res.ivw1,res.ivw.sd1,res.ivw.se1,res.ivw.cp1),
                           c(res.egger1,res.egger.sd1,res.egger.se1,res.egger.cp1),
                           c(res.median1,res.median.sd1,res.median.se1,res.median.cp1),
                           c(res.grapple1,res.grapple.sd1,res.grapple.se1,res.grapple.cp1),
                           c(res.mrbee1,res.mrbee.sd1,res.mrbee.se1,res.mrbee.cp1),
                           c(res.divw1,res.divw.sd1,res.divw.se1,res.divw.cp1),
                           c(res.adivw1,res.adivw.sd1,res.adivw.se1,res.adivw.cp1)),
                     rbind(c(res.ivw2,res.ivw.sd2,res.ivw.se2,res.ivw.cp2),
                           c(res.egger2,res.egger.sd2,res.egger.se2,res.egger.cp2),
                           c(res.median2,res.median.sd2,res.median.se2,res.median.cp2),
                           c(res.grapple2,res.grapple.sd2,res.grapple.se2,res.grapple.cp2),
                           c(res.mrbee2,res.mrbee.sd2,res.mrbee.se2,res.mrbee.cp2),
                           c(res.divw2,res.divw.sd2,res.divw.se2,res.divw.cp2),
                           c(res.adivw2,res.adivw.sd2,res.adivw.se2,res.adivw.cp2)),
                     rbind(c(res.ivw3,res.ivw.sd3,res.ivw.se3,res.ivw.cp3),
                           c(res.egger3,res.egger.sd3,res.egger.se3,res.egger.cp3),
                           c(res.median3,res.median.sd3,res.median.se3,res.median.cp3),
                           c(res.grapple3,res.grapple.sd3,res.grapple.se3,res.grapple.cp3),
                           c(res.mrbee3,res.mrbee.sd3,res.mrbee.se3,res.mrbee.cp3),
                           c(res.divw3,res.divw.sd3,res.divw.se3,res.divw.cp3),
                           c(res.adivw3,res.adivw.sd3,res.adivw.se3,res.adivw.cp3))
)

write.csv(res.summary, "tables/res.summary.main.table2.csv",row.names = FALSE)
