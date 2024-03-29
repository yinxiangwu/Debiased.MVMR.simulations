# Simulations for adIVW, three exposures example
# created by Yinxiang Wu on Jan 02, 2023
# updated on: 
# last update: Aug 15, 2023
task_id <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

source('F_stats_calculator.R') # adapted from MVMR package

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

all_divided_by_D = TRUE # change this to FALSE, if only first exposure's SNP-exposure associations are divided by D
folder_name <- ifelse(all_divided_by_D,"all_divided_by_D","first_divided_by_D")
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
n_rep = 100

for (m in 1:nrow(sim_settings)) {
      res.beta1 <- matrix(nrow = n_rep, ncol = 17)
      res.beta2 <- matrix(nrow = n_rep, ncol = 17)
      res.beta3 <- matrix(nrow = n_rep, ncol = 17)
      set.seed(123 + task_id + m*10000 + ifelse(all_divided_by_D,0,5234))
      # setting
      if (all_divided_by_D) {
        D = rep(sim_settings[m,]$D,3)
      } else {
        D = c(sim_settings[m,]$D,1,1)
      }
      rho = sim_settings[m,]$rho
      beta = beta_l[[sim_settings[m,]$beta]]
      par.beta.exposure[,1]<-rawdat_mvmr$LDL_beta/D[1]
      par.beta.exposure[,2]<-rawdat_mvmr$HDL_beta/D[2]
      par.beta.exposure[,3]<-rawdat_mvmr$Trg_beta/D[3]
      se.exposure[,1]<-rawdat_mvmr$LDL_se
      se.exposure[,2]<-rawdat_mvmr$HDL_se
      se.exposure[,3]<-rawdat_mvmr$Trg_se
      se.outcome<-rawdat_mvmr$SBP_se
      P <- matrix(rho,nrow = K, ncol = K)
      diag(P) <- 1
      
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
      beta.exposure <- lapply(1:p, function(j) {rmvnorm(1,mean=par.beta.exposure[j,],sigma=
                                                          Vj[[j]])}) %>% Reduce(rbind,.)
      # generating true and sample outcome
      par.beta.outcome <- par.beta.exposure %*% beta
      beta.outcome<-rnorm(p,mean=par.beta.outcome,sd=se.outcome)
      # kappa values
      IV_strength_matrix <- Reduce("+",lapply(1:p, function(j) {
        beta.exposure.V <- Vj_root_inv[[j]] %*% beta.exposure[j,]
        beta.exposure.V %*% t(beta.exposure.V)})) - p*diag(K)
      kappa_values <- min(eigen(IV_strength_matrix/sqrt(p))$values)
      # dIVW
      res.mvmr.divw <- mvmr.divw(beta.exposure, se.exposure, beta.outcome, se.outcome, gen_cor = P)
      res.dmvmr<- res.mvmr.divw$beta.hat
      res.dmvmr.se <- res.mvmr.divw$beta.se
      # adIVW
      res.mvmr.adivw <- mvmr.divw(beta.exposure, se.exposure, beta.outcome, se.outcome, gen_cor = P,phi_cand = NULL)
      res.dmvmr.adivw<- res.mvmr.adivw$beta.hat
      res.dmvmr.adivw.se <- res.mvmr.adivw$beta.se
      lambda_selected <- res.mvmr.adivw$phi_selected
      # IVW
      res.mvmr.ivw <- mvmr.ivw(beta.exposure, se.exposure, beta.outcome, se.outcome, gen_cor = P)
      res.dmvmr.IVW <- res.mvmr.ivw$beta.hat
      res.dmvmr.IVW.se <- res.mvmr.ivw$beta.se
      # MVMR median
      mr_median_obj <- mr_mvmedian(mr_mvinput(bx = beta.exposure,bxse = se.exposure,by = beta.outcome,byse = se.outcome),iterations = 10000)
      res.mvmr.median <- mr_median_obj$Estimate
      res.mvmr.median.se <- mr_median_obj$StdError
      # MVMR Egger
      mr_eggr_obj <- mr_mvegger(mr_mvinput(bx = beta.exposure,bxse = se.exposure,by = beta.outcome,byse = se.outcome),orientate = 1)
      res.mvmr.egger <- mr_eggr_obj$Estimate
      res.mvmr.egger.se <- mr_eggr_obj$StdError.Est
      # GRAPPLE
      grapple_dat <- data.frame(SNP = 1:p, beta.exposure,se.exposure,beta.outcome,se.outcome)
      colnames(grapple_dat) <- c('SNP',paste0('gamma_exp',1:3),paste0('se_exp',1:3),'gamma_out1','se_out1')
      P_grapple <- diag(K+1)
      P_grapple[1:K,1:K] <- P
      res.grapple <- grappleRobustEst(grapple_dat,
                                      cor.mat = P_grapple,
                                      plot.it = FALSE)
      res.mvmr.grapple <- res.grapple$beta.hat
      res.mvmr.grapple.se <- sqrt(diag(res.grapple$beta.var))
      # MRBEE
      R <- diag(K+1)
      R[1:3,1:3] <- P
      fit=MRBEE.IMRP(by=beta.outcome,bX=beta.exposure,byse=se.outcome,bXse=se.exposure,Rxy=R)
      res.mrbee.est <- fit$theta
      res.mrbee.se <- sqrt(diag(fit$covtheta))
      
      # Conditional F-stats
      F.data <- format_mvmr(BXGs = beta.exposure,
                            BYG = beta.outcome,
                            seBXGs = se.exposure,
                            seBYG = se.outcome,
                            RSID = 1:p)
      fres <- strength_mvmr2(r_input = F.data, gencov = lapply(1:p, function(j) {Vj[[j]] }))
      F_stat <- as.numeric(fres)
      
      tmp.res <- cbind(res.dmvmr.IVW,
                       res.dmvmr.IVW.se,
                       res.mvmr.egger,
                       res.mvmr.egger.se,
                       res.mvmr.median,
                       res.mvmr.median.se,
                       res.mvmr.grapple,
                       res.mvmr.grapple.se,
                       res.mrbee.est,
                       res.mrbee.se,
                       res.dmvmr,
                       res.dmvmr.se,
                       res.dmvmr.adivw,
                       res.dmvmr.adivw.se,
                       lambda_selected,
                       min_kappa_value,
                       F_stat)
      
      res.beta1[i,] <- tmp.res[1,]
      res.beta2[i,] <- tmp.res[2,]
      res.beta3[i,] <- tmp.res[3,]
      }
      write.csv(res.beta1, file = paste0("Res/final/3exp_example/",folder_name,"/beta",sim_settings[m,]$beta,'/rho',sim_settings[m,]$rho,'/res_beta1_D',sim_settings[m,]$D,'_job',task_id,'.csv'), row.names = FALSE)
      write.csv(res.beta2, file = paste0("Res/final/3exp_example/",folder_name,"/beta",sim_settings[m,]$beta,'/rho',sim_settings[m,]$rho,'/res_beta2_D',sim_settings[m,]$D,'_job',task_id,'.csv'), row.names = FALSE)
      write.csv(res.beta3, file = paste0("Res/final/3exp_example/",folder_name,"/beta",sim_settings[m,]$beta,'/rho',sim_settings[m,]$rho,'/res_beta3_D',sim_settings[m,]$D,'_job',task_id,'.csv'), row.names = FALSE)
      cat("Finish simulation with m =", m)
}
