# Simulations for adIVW, 10, 20, 30 exposures example
# created by Yinxiang Wu on Jan 02, 2023
# updated on: 
# last update: June 27, 2023
task_id <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

source('F_stats_calculator.R')

library(mvtnorm)
library(dplyr)
library(MVMR)
library(mr.divw)
library(parallel)
library(mvnfast)
library(GRAPPLE)
library(MendelianRandomization)
library(MRBEE)

data("bmi.cad")
p<-dim(bmi.cad)[1]

D_vec <- c(5.9,5.9,5.9)
K_vec = c(10,20,30)
set.seed(1423)
beta0 <- c(0.4, rnorm(29,mean = 0,sd = 0.15))
beta0 <- round(beta0,1)
pii <- c(0.3)
n_rep = 100

for (m in 1:3) {
      K = K_vec[m]
      beta <- beta0[1:K]
      
      dat.gen <- function() {
        set.seed(388)
        par.beta.exposure <- matrix(nrow=p,ncol=K)
        se.exposure<-par.beta.exposure
        par.beta.exposure[,1]<-bmi.cad$beta.exposure/D_vec[m]
        se.exposure[,1]<-bmi.cad$se.exposure
        tau <- sum(par.beta.exposure[,1]^2)/10
        for (k in 2:K) {
          ita_k <- rbinom(p,1,prob = pii)
          delta_k <- rnorm(p,mean = 0, sd = sqrt(tau))
          par.beta.exposure[,k] <- ita_k * par.beta.exposure[,1] + (1 - ita_k)*delta_k
          se.exposure[,k] <- bmi.cad$se.exposure
        }
        return(list(dat1 = par.beta.exposure, dat2 = se.exposure))
      }
      
      generated.dat <- dat.gen()
      
      par.beta.exposure <-  generated.dat$dat1
      se.exposure <- generated.dat$dat2
      
      par.beta.outcome <- par.beta.exposure %*% beta
      se.outcome <- bmi.cad$se.outcome
    
      P <- matrix(0.3,nrow = K, ncol = K)
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
      kappa_value_true <- min(eigen(IV_strength_matrix/sqrt(p))$values)
      
      res.list <- vector(mode = "list", length = K)
      for (k in 1:K) res.list[[k]] <- matrix(nrow = n_rep, ncol = 18)
      
      set.seed(436 + task_id + m*10000)
      for (i in 1:n_rep) {
      beta.exposure <- lapply(1:p, function(j) {rmvnorm(1,mean=par.beta.exposure[j,],sigma=
                                                          Vj[[j]])}) %>% Reduce(rbind,.)
      # generating true and sample outcome
      beta.outcome<-rnorm(p,mean=par.beta.outcome,sd=se.outcome)
      # IV strength parameter
      IV_strength_matrix <- Reduce("+",lapply(1:p, function(j) {
        beta.exposure.V <- Vj_root_inv[[j]] %*% beta.exposure[j,]
        beta.exposure.V %*% t(beta.exposure.V)})) - p*diag(K)
      kappa_values <- eigen(IV_strength_matrix/sqrt(p))$values
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
      mr_median_obj <- mr_mvmedian(mr_mvinput(bx = beta.exposure,bxse = se.exposure,by = beta.outcome,byse = se.outcome))
      res.mvmr.median <- mr_median_obj$Estimate
      res.mvmr.median.se <- mr_median_obj$StdError
      # MVMR Egger
      mr_eggr_obj <- mr_mvegger(mr_mvinput(bx = beta.exposure,bxse = se.exposure,by = beta.outcome,byse = se.outcome),orientate = 1)
      res.mvmr.egger <- mr_eggr_obj$Estimate
      res.mvmr.egger.se <- mr_eggr_obj$StdError.Est
      # GRAPPLE
      grapple_dat <- data.frame(SNP = 1:p, beta.exposure,se.exposure,beta.outcome,se.outcome)
      colnames(grapple_dat) <- c('SNP',paste0('gamma_exp',1:K),paste0('se_exp',1:K),'gamma_out1','se_out1')
      P_grapple <- diag(K+1)
      P_grapple[1:K,1:K] <- P
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
      
      F_hist <- as.numeric(fres)
      
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
                       kappa_value_true,
                       kappa_values,
                       F_hist)
      
      for (k in 1:K) {res.list[[k]][i,] <- tmp.res[k,]}
      }
      for (k in 1:K) {write.csv(res.list[[k]], file = paste0('Res/final/10_20_30_exp_example/K',K,'/res_dmvmr_beta',k,'_K',K,'_job',task_id,'.csv'), row.names = FALSE)
      }
      cat("Finish simulation with m =", m)
}
