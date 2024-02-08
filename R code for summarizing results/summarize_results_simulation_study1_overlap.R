# Summarize simulation results for adIVW, three exposures example with overlapping exposure and outcome datasets
# created by Yinxiang Wu on Jan 02, 2023
# last update: Dec 30, 2023
res.summary <- NULL
K = 3
D_vec <- c(1,1,1,2.5,4.5,5.5) 
beta_list <- list(b1 = c(0.6,0.1,-0.2), b2 = c(-0.2,0.1,-0.2), b3 = c(0,0,0), b4 = c(-0.5,-0.7,0.3), b5 = c(-0.5,-0.7,0.3), b6 = c(-0.5,-0.7,0.3))
rho <- c(0.7)
n_rep = 100

tbl_to_present <- list("supp_tables13" = 1, # the row numbers corresponding to main table
                       "supp_tables14" = 2,
                       "supp_tables15" = 3,
                       "supp_tables16" = 4:6)

tbl_num <- "supp_tables16" # to summarize the results for the other tables, change this 

for (m in tbl_to_present[[tbl_num]]) {
  D = D_vec[m]
  # beta 1
  res.beta1 <- NULL
  for (i in c(1:100)) {
    tryCatch({
      res.beta1 <- rbind(res.beta1, read.csv(file = paste0('Res/for_pub/3exp_example/with_overlap/res_beta1_b',m,'_job',i,'.csv')))
    }, error=function(e){})
  }
  
  adivw.est1 <- mean(res.beta1[,1])
  adivw.sd1 <- sd(res.beta1[,1])
  adivw.se1 <- mean(res.beta1[,2])
  adivw.cp1 <- mean(res.beta1[,1] + qnorm(0.025) * res.beta1[,2] < beta[1] & beta[1] < 
                      res.beta1[,1] + qnorm(0.975) * res.beta1[,2])
  
  adivw2.est1 <- mean(res.beta1[,3])
  adivw2.sd1 <- sd(res.beta1[,3])
  adivw2.se1 <- mean(res.beta1[,4])
  adivw2.cp1 <- mean(res.beta1[,3] + qnorm(0.025) * res.beta1[,4] < beta[1] & beta[1] < 
                      res.beta1[,3] + qnorm(0.975) * res.beta1[,4])
  
  grapple.est1 <- mean(res.beta1[,5])
  grapple.sd1 <- sd(res.beta1[,5])
  grapple.se1 <- mean(res.beta1[,6])
  grapple.cp1 <- mean(res.beta1[,5] + qnorm(0.025) * res.beta1[,6] < beta[1] & beta[1] < 
                       res.beta1[,5] + qnorm(0.975) * res.beta1[,6])
  

  # beta 2
  res.beta2 <- NULL
  for (i in c(1:100)) {
    tryCatch({
      res.beta2 <- rbind(res.beta2, read.csv(file = paste0('Res/for_pub/3exp_example/with_overlap/res_beta2_b',m,'_job',i,'.csv')))
    }, error=function(e){})
  }
  
  adivw.est2 <- mean(res.beta2[,1])
  adivw.sd2 <- sd(res.beta2[,1])
  adivw.se2 <- mean(res.beta2[,2])
  adivw.cp2 <- mean(res.beta2[,1] + qnorm(0.025) * res.beta2[,2] < beta[2] & beta[2] < 
                      res.beta2[,1] + qnorm(0.975) * res.beta2[,2])
  
  adivw2.est2 <- mean(res.beta2[,3])
  adivw2.sd2 <- sd(res.beta2[,3])
  adivw2.se2 <- mean(res.beta2[,4])
  adivw2.cp2 <- mean(res.beta2[,3] + qnorm(0.025) * res.beta2[,4] < beta[2] & beta[2] < 
                      res.beta2[,3] + qnorm(0.975) * res.beta2[,4])
  
  grapple.est2 <- mean(res.beta2[,5])
  grapple.sd2 <- sd(res.beta2[,5])
  grapple.se2 <- mean(res.beta2[,6])
  grapple.cp2 <- mean(res.beta2[,5] + qnorm(0.025) * res.beta2[,6] < beta[2] & beta[2] < 
                        res.beta2[,5] + qnorm(0.975) * res.beta2[,6])
  
  
  
  # beta 3
  res.beta3 <- NULL
  for (i in c(1:100)) {
    tryCatch({
      res.beta3 <- rbind(res.beta3, read.csv(file = paste0('Res/for_pub/3exp_example/with_overlap/res_beta3_b',m,'_job',i,'.csv')))
      }, error=function(e){})
  }
  
  adivw.est3 <- mean(res.beta3[,1])
  adivw.sd3 <- sd(res.beta3[,1])
  adivw.se3 <- mean(res.beta3[,2])
  adivw.cp3 <- mean(res.beta3[,1] + qnorm(0.025) * res.beta3[,2] < beta[3] & beta[3] < 
                      res.beta3[,1] + qnorm(0.975) * res.beta3[,2])
  
  adivw2.est3 <- mean(res.beta3[,3])
  adivw2.sd3 <- sd(res.beta3[,3])
  adivw2.se3 <- mean(res.beta3[,4])
  adivw2.cp3 <- mean(res.beta3[,3] + qnorm(0.025) * res.beta3[,4] < beta[3] & beta[3] < 
                      res.beta3[,3] + qnorm(0.975) * res.beta3[,4])
  
  grapple.est3 <- mean(res.beta3[,5])
  grapple.sd3 <- sd(res.beta3[,5])
  grapple.se3 <- mean(res.beta3[,6])
  grapple.cp3 <- mean(res.beta3[,5] + qnorm(0.025) * res.beta3[,6] < beta[3] & beta[3] < 
                        res.beta3[,5] + qnorm(0.975) * res.beta3[,6])
  
  
  kappa <- mean(res.beta3[,8])
  res.summary <- rbind(res.summary, 
                       rbind(c(D,rho,kappa,adivw.est1,adivw.sd1,adivw.se1,adivw.cp1,
                               adivw.est2,adivw.sd2,adivw.se2,adivw.cp2,
                               adivw.est3,adivw.sd3,adivw.se3,adivw.cp3),
                             c(D,rho,kappa,adivw2.est1,adivw2.sd1,adivw2.se1,adivw2.cp1,
                               adivw2.est2,adivw2.sd2,adivw2.se2,adivw2.cp2,
                               adivw2.est3,adivw2.sd3,adivw2.se3,adivw2.cp3),
                             c(D,rho,kappa,grapple.est1,grapple.sd1,grapple.se1,grapple.cp1,
                               grapple.est2,grapple.sd2,grapple.se2,grapple.cp2,
                               grapple.est3,grapple.sd3,grapple.se3,grapple.cp3)
                       ))
  
}

write.csv(res.summary, paste0("tables/res.summary.",tbl_num,".csv"),row.names = FALSE)