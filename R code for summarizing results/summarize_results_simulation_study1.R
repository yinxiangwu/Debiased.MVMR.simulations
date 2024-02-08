# Summarize simulation results for adIVW, three exposures example
# created by Yinxiang Wu on Jan 02, 2023
# last update: Dec 30, 2023
res.summary <- NULL
K = 3
all_divided_by_D = TRUE # change this to FALSE, if only first exposure's SNP-exposure associations are divided by D
if (all_divided_by_D) {
  D_vec <- c(3, 5.5, 6.5)
} else {
  D_vec <- c(4.5, 7, 9.25)
}
beta_l <- list(b1 = c(-0.5,-0.7,0.3),b2 = c(0,0.5,-0.5))
rho_vec <- c(0.3, 0.7)
sim_settings <- expand.grid(1:length(beta_l),D_vec,rho_vec)
if (all_divided_by_D) {
  sim_settings[7:12,2] <-  c(2.5, 2.5, 4.5, 4.5, 5.5, 5.5) 
} else {
  sim_settings[7:12,2] <-  c(4, 4, 6.5, 6.5, 9, 9) 
}
colnames(sim_settings) <- c('beta','D','rho')

tbl_to_present <- list("main_table1" = c(7,9,11), # the row numbers corresponding to main table
                       "supp_tables1" = c(1,3,5),
                       "supp_tables2" = c(2,4,6),
                       "supp_tables3" = c(8,10,12))

tbl_num <- "main_table1" # to summarize the results for the other tables, change this 
for (m in tbl_to_present[[tbl_num]]) {
  D = sim_settings[m,]$D
  rho = sim_settings[m,]$rho
  beta = beta_l[[sim_settings[m,]$beta]]
  # beta 1
  res.beta1 <- NULL
  for (i in c(1:100)) {
    tryCatch({
      res.beta1 <- rbind(res.beta1, read.csv(file = paste0("Res/for_pub/3exp_example/all_divided_by_D/beta",sim_settings[m,]$beta,'/rho',sim_settings[m,]$rho,'/res_beta1_D',sim_settings[m,]$D,'_job',i,'.csv')))
    }, error=function(e){})
  }
  
  ivw.est1 <- mean(res.beta1[,1])
  ivw.sd1 <- sd(res.beta1[,1])
  ivw.se1 <- mean(res.beta1[,2])
  ivw.cp1 <- mean(res.beta1[,1] + qnorm(0.025) * res.beta1[,2] < beta[1] & beta[1] < 
                     res.beta1[,1] + qnorm(0.975) * res.beta1[,2])
  
  egger.est1 <- mean(res.beta1[,3])
  egger.sd1 <- sd(res.beta1[,3])
  egger.se1 <- mean(res.beta1[,4])
  egger.cp1 <- mean(res.beta1[,3] + qnorm(0.025) * res.beta1[,4] < beta[1] & beta[1] < 
                       res.beta1[,3] + qnorm(0.975) * res.beta1[,4])
  
  median.est1 <- mean(res.beta1[,5])
  median.sd1 <- sd(res.beta1[,5])
  median.se1 <- mean(res.beta1[,6])
  median.cp1 <- mean(res.beta1[,5] + qnorm(0.025) * res.beta1[,6] < beta[1] & beta[1] < 
                       res.beta1[,5] + qnorm(0.975) * res.beta1[,6])
  
  grapple.est1 <- mean(res.beta1[,7])
  grapple.sd1 <- sd(res.beta1[,7])
  grapple.se1 <- mean(res.beta1[,8])
  grapple.cp1 <- mean(res.beta1[,7] + qnorm(0.025) * res.beta1[,8] < beta[1] & beta[1] < 
                       res.beta1[,7] + qnorm(0.975) * res.beta1[,8])
  
  mrbee.est1 <- mean(res.beta1[,9])
  mrbee.sd1 <- sd(res.beta1[,9])
  mrbee.se1 <- mean(res.beta1[,10])
  mrbee.cp1 <- mean(res.beta1[,9] + qnorm(0.025) * res.beta1[,10] < beta[1] & beta[1] < 
                        res.beta1[,9] + qnorm(0.975) * res.beta1[,10])
  
  divw.est1 <- mean(res.beta1[,11])
  divw.sd1 <- sd(res.beta1[,11])
  divw.se1 <- mean(res.beta1[,12])
  divw.cp1 <- mean(res.beta1[,11] + qnorm(0.025) * res.beta1[,12] < beta[1] & beta[1] < 
                      res.beta1[,11] + qnorm(0.975) * res.beta1[,12])
  
  adivw.est1 <- mean(res.beta1[,13])
  adivw.sd1 <- sd(res.beta1[,13])
  adivw.se1 <- mean(res.beta1[,14])
  adivw.cp1 <- mean(res.beta1[,13] + qnorm(0.025) * res.beta1[,14] < beta[1] & beta[1] < 
                     res.beta1[,13] + qnorm(0.975) * res.beta1[,14])
  
  # beta 2
  res.beta2 <- NULL
  for (i in c(1:100)) {
    tryCatch({
      res.beta2 <- rbind(res.beta2, read.csv(file = paste0("Res/for_pub/3exp_example/all_divided_by_D/beta",sim_settings[m,]$beta,'/rho',sim_settings[m,]$rho,'/res_beta2_D',sim_settings[m,]$D,'_job',i,'.csv')))
    }, error=function(e){})
  }
  
  ivw.est2 <- mean(res.beta2[,1])
  ivw.sd2 <- sd(res.beta2[,1])
  ivw.se2 <- mean(res.beta2[,2])
  ivw.cp2 <- mean(res.beta2[,1] + qnorm(0.025) * res.beta2[,2] < beta[2] & beta[2] < 
                    res.beta2[,1] + qnorm(0.975) * res.beta2[,2])
  
  egger.est2 <- mean(res.beta2[,3])
  egger.sd2 <- sd(res.beta2[,3])
  egger.se2 <- mean(res.beta2[,4])
  egger.cp2 <- mean(res.beta2[,3] + qnorm(0.025) * res.beta2[,4] < beta[2] & beta[2] < 
                      res.beta2[,3] + qnorm(0.975) * res.beta2[,4])
  
  median.est2 <- mean(res.beta2[,5])
  median.sd2 <- sd(res.beta2[,5])
  median.se2 <- mean(res.beta2[,6])
  median.cp2 <- mean(res.beta2[,5] + qnorm(0.025) * res.beta2[,6] < beta[2] & beta[2] < 
                       res.beta2[,5] + qnorm(0.975) * res.beta2[,6])
  
  grapple.est2 <- mean(res.beta2[,7])
  grapple.sd2 <- sd(res.beta2[,7])
  grapple.se2 <- mean(res.beta2[,8])
  grapple.cp2 <- mean(res.beta2[,7] + qnorm(0.025) * res.beta2[,8] < beta[2] & beta[2] < 
                        res.beta2[,7] + qnorm(0.975) * res.beta2[,8])
  
  mrbee.est2 <- mean(res.beta2[,9])
  mrbee.sd2 <- sd(res.beta2[,9])
  mrbee.se2 <- mean(res.beta2[,10])
  mrbee.cp2 <- mean(res.beta2[,9] + qnorm(0.025) * res.beta2[,10] < beta[2] & beta[2] < 
                      res.beta2[,9] + qnorm(0.975) * res.beta2[,10])
  
  divw.est2 <- mean(res.beta2[,11])
  divw.sd2 <- sd(res.beta2[,11])
  divw.se2 <- mean(res.beta2[,12])
  divw.cp2 <- mean(res.beta2[,11] + qnorm(0.025) * res.beta2[,12] < beta[2] & beta[2] < 
                     res.beta2[,11] + qnorm(0.975) * res.beta2[,12])
  
  adivw.est2 <- mean(res.beta2[,13])
  adivw.sd2 <- sd(res.beta2[,13])
  adivw.se2 <- mean(res.beta2[,14])
  adivw.cp2 <- mean(res.beta2[,13] + qnorm(0.025) * res.beta2[,14] < beta[2] & beta[2] < 
                      res.beta2[,13] + qnorm(0.975) * res.beta2[,14])
  
  
  # beta 3
  res.beta3 <- NULL
  for (i in c(1:100)) {
    tryCatch({
      res.beta3 <- rbind(res.beta3, read.csv(file = paste0("Res/for_pub/3exp_example/all_divided_by_D/beta",sim_settings[m,]$beta,'/rho',sim_settings[m,]$rho,'/res_beta3_D',sim_settings[m,]$D,'_job',i,'.csv')))
      }, error=function(e){})
  }
  
  ivw.est3 <- mean(res.beta3[,1])
  ivw.sd3 <- sd(res.beta3[,1])
  ivw.se3 <- mean(res.beta3[,2])
  ivw.cp3 <- mean(res.beta3[,1] + qnorm(0.025) * res.beta3[,2] < beta[3] & beta[3] < 
                    res.beta3[,1] + qnorm(0.975) * res.beta3[,2])
  
  egger.est3 <- mean(res.beta3[,3])
  egger.sd3 <- sd(res.beta3[,3])
  egger.se3 <- mean(res.beta3[,4])
  egger.cp3 <- mean(res.beta3[,3] + qnorm(0.025) * res.beta3[,4] < beta[3] & beta[3] < 
                      res.beta3[,3] + qnorm(0.975) * res.beta3[,4])
  
  median.est3 <- mean(res.beta3[,5])
  median.sd3 <- sd(res.beta3[,5])
  median.se3 <- mean(res.beta3[,6])
  median.cp3 <- mean(res.beta3[,5] + qnorm(0.025) * res.beta3[,6] < beta[3] & beta[3] < 
                       res.beta3[,5] + qnorm(0.975) * res.beta3[,6])
  
  grapple.est3 <- mean(res.beta3[,7])
  grapple.sd3 <- sd(res.beta3[,7])
  grapple.se3 <- mean(res.beta3[,8])
  grapple.cp3 <- mean(res.beta3[,7] + qnorm(0.025) * res.beta3[,8] < beta[3] & beta[3] < 
                        res.beta3[,7] + qnorm(0.975) * res.beta3[,8])
  
  mrbee.est3 <- mean(res.beta3[,9])
  mrbee.sd3 <- sd(res.beta3[,9])
  mrbee.se3 <- mean(res.beta3[,10])
  mrbee.cp3 <- mean(res.beta3[,9] + qnorm(0.025) * res.beta3[,10] < beta[3] & beta[3] < 
                      res.beta3[,9] + qnorm(0.975) * res.beta3[,10])
  
  divw.est3 <- mean(res.beta3[,11])
  divw.sd3 <- sd(res.beta3[,11])
  divw.se3 <- mean(res.beta3[,12])
  divw.cp3 <- mean(res.beta3[,11] + qnorm(0.025) * res.beta3[,12] < beta[3] & beta[3] < 
                     res.beta3[,11] + qnorm(0.975) * res.beta3[,12])
  
  adivw.est3 <- mean(res.beta3[,13])
  adivw.sd3 <- sd(res.beta3[,13])
  adivw.se3 <- mean(res.beta3[,14])
  adivw.cp3 <- mean(res.beta3[,13] + qnorm(0.025) * res.beta3[,14] < beta[3] & beta[3] < 
                      res.beta3[,13] + qnorm(0.975) * res.beta3[,14])
  
  
  kappa <- mean(res.beta3[,16])
  res.summary <- rbind(res.summary, 
                       rbind(c(D,rho,kappa,ivw.est1,ivw.sd1,ivw.se1,ivw.cp1,
                               ivw.est2,ivw.sd2,ivw.se2,ivw.cp2,
                               ivw.est3,ivw.sd3,ivw.se3,ivw.cp3),
                             c(D,rho,kappa,egger.est1,egger.sd1,egger.se1,egger.cp1,
                               egger.est2,egger.sd2,egger.se2,egger.cp2,
                               egger.est3,egger.sd3,egger.se3,egger.cp3),
                             c(D,rho,kappa,median.est1,median.sd1,median.se1,median.cp1,
                               median.est2,median.sd2,median.se2,median.cp2,
                               median.est3,median.sd3,median.se3,median.cp3),
                             c(D,rho,kappa,grapple.est1,grapple.sd1,grapple.se1,grapple.cp1,
                               grapple.est2,grapple.sd2,grapple.se2,grapple.cp2,
                               grapple.est3,grapple.sd3,grapple.se3,grapple.cp3),
                             c(D,rho,kappa,mrbee.est1,mrbee.sd1,mrbee.se1,mrbee.cp1,
                               mrbee.est2,mrbee.sd2,mrbee.se2,mrbee.cp2,
                               mrbee.est3,mrbee.sd3,mrbee.se3,mrbee.cp3),
                             c(D,rho,kappa,divw.est1,divw.sd1,divw.se1,divw.cp1,
                               divw.est2,divw.sd2,divw.se2,divw.cp2,
                               divw.est3,divw.sd3,divw.se3,divw.cp3),
                             c(D,rho,kappa,adivw.est1,adivw.sd1,adivw.se1,adivw.cp1,
                               adivw.est2,adivw.sd2,adivw.se2,adivw.cp2,
                               adivw.est3,adivw.sd3,adivw.se3,adivw.cp3)
                       ))
  
}

write.csv(res.summary, paste0("tables/res.summary.",tbl_num,".csv"),row.names = FALSE)