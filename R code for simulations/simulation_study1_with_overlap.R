# Simulations for adIVW, three exposures example with overlapping exposure and outcome datasets
# created by Yinxiang Wu on Jan 02, 2023
# updated on: 
# last update: Dec 17, 2023
task_id <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

source('F_stats_calculator.R')

library(mvtnorm)
library(dplyr)
library(GRAPPLE)
library(MendelianRandomization)
library(MVMR)
library(MRBEE)
library(mr.divw)

data("rawdat_mvmr")
p<-dim(rawdat_mvmr)[1]

K = 3

par.beta.exposure<-matrix(nrow=p,ncol=K)
se.exposure<-par.beta.exposure

D_vec <- c(1,1,1,2.5,4.5,5.5) 
beta_list <- list(b1 = c(0.6,0.1,-0.2), b2 = c(-0.2,0.1,-0.2), b3 = c(0,0,0), b4 = c(-0.5,-0.7,0.3), b5 = c(-0.5,-0.7,0.3), b6 = c(-0.5,-0.7,0.3))
rho <- c(0.7)
n_rep = 100

for (m in 1:6) {
  res.beta1 <- matrix(nrow = n_rep, ncol = 9)
  res.beta2 <- matrix(nrow = n_rep, ncol = 9)
  res.beta3 <- matrix(nrow = n_rep, ncol = 9)
  # set.seed(12512341 + task_id + m)
  # setting
  par.beta.exposure[,1]<-rawdat_mvmr$LDL_beta/D_vec[m]
  par.beta.exposure[,2]<-rawdat_mvmr$HDL_beta/D_vec[m]
  par.beta.exposure[,3]<-rawdat_mvmr$Trg_beta/D_vec[m]
  se.exposure[,1]<-rawdat_mvmr$LDL_se
  se.exposure[,2]<-rawdat_mvmr$HDL_se
  se.exposure[,3]<-rawdat_mvmr$Trg_se
  beta <- beta_list[[m]]
  par.beta.outcome <- par.beta.exposure %*% beta
  se.outcome<-rawdat_mvmr$SBP_se
  P_full <- matrix(rho,nrow = K+1, ncol = K+1)
  diag(P_full) <- 1 
  P <- matrix(rho,nrow = K, ncol = K)
  diag(P) <- 1
        
  full.par <- cbind(par.beta.exposure, par.beta.outcome)
  full.se <- cbind(se.exposure,se.outcome)
  V_full <- lapply(1:p, function(j) diag(full.se[j,]) %*% P_full %*% diag(full.se[j,]))
        
  Vj <- lapply(1:p, function(j) diag(se.exposure[j,]) %*% P %*% diag(se.exposure[j,]))
  P_eigen <- eigen(P)
  P_root_inv <- P_eigen$vectors %*% diag(1/sqrt(P_eigen$values)) %*% t(P_eigen$vectors)
  Vj_root_inv <- lapply(1:p, function(j) {
    P_root_inv %*% diag(1/se.exposure[j,])
  })
  # calcualte IV strength parameter
  IV_strength_matrix <- Reduce("+",lapply(1:p, function(j) {
    beta.exposure.V <- Vj_root_inv[[j]] %*% par.beta.exposure[j,]
    beta.exposure.V %*% t(beta.exposure.V)}))
  min_kappa_value <- min(eigen(IV_strength_matrix/sqrt(p))$values)
        
  for (i in 1:n_rep) {
      dat <- lapply(1:p, function(j) {rmvnorm(1,mean=full.par[j,],sigma=V_full[[j]])}) %>% Reduce(rbind,.)
      beta.exposure <- dat[,1:3]
      beta.outcome <- dat[,4]
      # adIVW assuming independent exposure and outcome datasets
      res.mvmr.adivw <- mvmr.divw(beta.exposure, se.exposure, beta.outcome, se.outcome, gen_cor = P,phi_cand = NULL)
      res.dmvmr.adivw<- res.mvmr.adivw$beta.hat
      res.dmvmr.adivw.se <- res.mvmr.adivw$beta.se
      # adIVW allowing for overlapping exposure and outcome datasets
      res.mvmr.adivw.overlap <- mvmr.divw(beta.exposure, se.exposure, beta.outcome, se.outcome, gen_cor = P_full,phi_cand = NULL,overlap = TRUE)
      res.dmvmr.adivw.overlap <- res.mvmr.adivw.overlap$beta.hat
      res.dmvmr.adivw.overlap.se <- res.mvmr.adivw.overlap$beta.se
      lambda_selected <- res.mvmr.adivw.overlap$phi_selected
      # GRAPPLE
      grapple_dat <- data.frame(SNP = 1:p, beta.exposure,se.exposure,beta.outcome,se.outcome)
      colnames(grapple_dat) <- c('SNP',paste0('gamma_exp',1:3),paste0('se_exp',1:3),'gamma_out1','se_out1')
      P_grapple <- P_full
      res.grapple <- grappleRobustEst(grapple_dat,
                                      cor.mat = P_grapple,
                                      plot.it = FALSE)
      res.mvmr.grapple <- res.grapple$beta.hat
      res.mvmr.grapple.se <- sqrt(diag(res.grapple$beta.var))
      
      # Conditional F-stats
      F.data <- format_mvmr(BXGs = beta.exposure,
                            BYG = beta.outcome,
                            seBXGs = se.exposure,
                            seBYG = se.outcome,
                            RSID = 1:p)
      fres <- strength_mvmr2(r_input = F.data, gencov = lapply(1:p, function(j) {Vj[[j]] }))
      F_stat <- as.numeric(fres)
      
      tmp.res <- cbind(
                       res.dmvmr.adivw,
                       res.dmvmr.adivw.se,
                       res.dmvmr.adivw.overlap,
                       res.dmvmr.adivw.overlap.se,
                       res.mvmr.grapple,
                       res.mvmr.grapple.se,
                       lambda_selected,
                       min_kappa_value,
                       F_stat)
      
      res.beta1[i,] <- tmp.res[1,]
      res.beta2[i,] <- tmp.res[2,]
      res.beta3[i,] <- tmp.res[3,]
  }
  write.csv(res.beta1, file = paste0('Res/for_pub/3exp_example/with_overlap/res_beta1_b',m,'_job',task_id,'.csv'), row.names = FALSE)
  write.csv(res.beta2, file = paste0('Res/for_pub/3exp_example/with_overlap/res_beta2_b',m,'_job',task_id,'.csv'), row.names = FALSE)
  write.csv(res.beta3, file = paste0('Res/for_pub/3exp_example/with_overlap/res_beta3_b',m,'_job',task_id,'.csv'), row.names = FALSE)
}
