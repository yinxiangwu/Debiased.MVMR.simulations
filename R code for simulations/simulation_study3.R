# Individual-level simulation for adIVW, three exposures
# created by Yinxiang Wu on Jan 02, 2023
# updated on: 
# last update: Jan 17, 2023
task_id <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

source('F_stats_calculator.R')
source('data_gen_individual.R')

library(mvtnorm)
library(dplyr)
library(GRAPPLE)
library(MendelianRandomization)
library(MVMR)
library(mr.divw)
library(MRBEE)

K = 3
n_rep = 100
p.val.threshold <- 10^{-2}

res.beta1 <- matrix(nrow = n_rep, ncol = 18)
res.beta2 <- matrix(nrow = n_rep, ncol = 18)
res.beta3 <- matrix(nrow = n_rep, ncol = 18)
lambda_choice <- rep(NA, n_rep)

for (i in 1:n_rep) {
      set.seed(4327+(task_id-1)*n_rep+i)
      # generate data
      dat <- data_gen_individual(hsquare = 0.1)
      # p-value thresholding selection of SNPs
      sel.SNPs <- apply(dat[,paste0('pval.selection',1:3)],1,function(x) any(x < p.val.threshold/K))
      sel.SNPs.cor <- apply(dat[,paste0('pval.selection',1:3)],1,function(x) all(x > 0.5))
      
      # exposure and outcome datasets used to compute causal effects
      beta.exposure <- dat[sel.SNPs,paste0('beta.exposure',1:3)] %>% as.matrix()
      se.exposure <- dat[sel.SNPs,paste0('se.exposure',1:3)] %>% as.matrix()
      beta.outcome <- dat[sel.SNPs,'beta.outcome']
      se.outcome <- dat[sel.SNPs,'se.outcome']
      p <- sum(sel.SNPs)
      
      # exposure and outcome datasets used to compute correlation matrix
      beta.exposure.cor <- dat[sel.SNPs.cor,paste0('beta.exposure',1:3)]
      se.exposure.cor <- dat[sel.SNPs.cor,paste0('se.exposure',1:3)]
      beta.outcome.cor <- dat[sel.SNPs.cor,'beta.outcome']
      se.outcome.cor <- dat[sel.SNPs.cor,'se.outcome']
      
      
      zz <- cbind(beta.exposure.cor, se.exposure.cor, beta.outcome.cor, se.outcome.cor)
      z.values <- cbind(beta.exposure.cor/se.exposure.cor,
                        beta.outcome.cor/se.outcome.cor)
      colnames(z.values) <- c(paste0("exposure", 1:K),
                              'outcome')
      
      z.values <- as.matrix(z.values)
      z.values <- z.values[rowSums(is.na(z.values)) == 0,  , drop = F]
      covv <- t(z.values) %*% z.values / nrow(z.values)
      varr <- colMeans(z.values^2, na.rm = T)
      corr <- t(covv / sqrt(varr))/sqrt(varr)
      P <- corr[1:K, 1:K] # estimated correlation matrix
      
      Vj <- lapply(1:p, function(j) diag(se.exposure[j,]) %*% P %*% diag(se.exposure[j,]))
      P_eigen <- eigen(P)
      P_root_inv <- P_eigen$vectors %*% diag(1/sqrt(P_eigen$values)) %*% t(P_eigen$vectors)
      Vj_root_inv <- lapply(1:p, function(j) {
        P_root_inv %*% diag(1/se.exposure[j,])
      })
      # calcualte IV strength parameter
      IV_strength_matrix <- Reduce("+",lapply(1:p, function(j) {
        beta.exposure.V <- Vj_root_inv[[j]] %*% beta.exposure[j,]
        beta.exposure.V %*% t(beta.exposure.V)})) - p*diag(K)
      min_kappa_value <- min(eigen(IV_strength_matrix/sqrt(p))$values)
      
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
      P_grapple <- corr
      res.grapple <- grappleRobustEst(grapple_dat,
                                      cor.mat = P_grapple,
                                      plot.it = FALSE)
      res.mvmr.grapple <- res.grapple$beta.hat
      res.mvmr.grapple.se <- sqrt(diag(res.grapple$beta.var))
      # MRBEE
      R <- diag(K+1)
      R[1:K,1:K] <- P
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
                       F_stat,
                       p)
      
      res.beta1[i,] <- tmp.res[1,]
      res.beta2[i,] <- tmp.res[2,]
      res.beta3[i,] <- tmp.res[3,]
      }
      write.csv(res.beta1, file = paste0('Res/final/individual_lvl_sim/res_beta1_h01_job',task_id,'.csv'), row.names = FALSE)
      write.csv(res.beta2, file = paste0('Res/final/individual_lvl_sim/res_beta2_h01_job',task_id,'.csv'), row.names = FALSE)
      write.csv(res.beta3, file = paste0('Res/final/individual_lvl_sim/res_beta3_h01_job',task_id,'.csv'), row.names = FALSE)
      cat("Finish simulation with m =", m)

