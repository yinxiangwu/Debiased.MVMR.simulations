data_gen_individual <- function(hsquare = 0.2) {
    s<-1000
    n<-1e4
    p<-2000
    h<-sqrt(hsquare) # 20% total heritability
    beta0<-c(1,-0.5, 0.5)
    standard_normal_random_effect<-mr.divw:::standard_normal_random_effect[1:p]
    gamma_coef1<-standard_normal_random_effect*(h/sqrt(s))*sqrt(2)
    if(s<p){
      gamma_coef1[(s+1):p]<-0
    }
    standard_normal_random_effect<-mr.divw:::standard_normal_random_effect[(p+1):(2*p)]
    gamma_coef2<-standard_normal_random_effect*(h/sqrt(s))*sqrt(2)
    if(s<p){
      gamma_coef2[(s+1):p]<-0
    }
    standard_normal_random_effect<-mr.divw:::standard_normal_random_effect[(2*p+1):(3*p)]
    gamma_coef3<-standard_normal_random_effect*(h/sqrt(s))*sqrt(2)
    if(s<p){
      gamma_coef3[(s+1):p]<-0
    }
    gamma_coef <- cbind(gamma_coef1, gamma_coef2, gamma_coef3)
    data_gen_onesample<-function(){
        z<-matrix(nrow=n,ncol=p)
        for(j in 1:p){
          tmp<-rmultinom(n,1,c(0.25,0.5,0.25)) # z's are independent
          z[,j]<-0*(tmp[1,]==1)+1*(tmp[2,]==1)+2*(tmp[3,]==1)
        }
        u<-rnorm(n,0,sqrt((1-h^2)*0.6))
        x<-z[,1:s] %*% gamma_coef[1:s, ] + u + rmvnorm(n,mean=c(0,0,0),sigma = diag(rep((1-h^2)*0.4,3)))
        y<-10+x%*%beta0+u+rnorm(n)
        return(list(z=z,x=x,y=y))
    }
    tmp1<-data_gen_onesample() # exposure dataset
    tmp2<-data_gen_onesample() # outcome dataset
    tmp3<-data_gen_onesample() # selection dataset
    full_df<-data.frame(beta.exposure1=numeric(p))
    full_df[,c(paste0('beta.exposure',1:3),
               paste0('se.exposure',1:3))] <- sapply(1:p, function(i) {
        m <- lm(tmp1$x ~ tmp1$z[,i])
        m_summary <- summary(m)
        return(c(m$coef[2,],
                 m_summary$`Response gamma_coef1`$coefficients[2,2],
                 m_summary$`Response gamma_coef2`$coefficients[2,2],
                 m_summary$`Response gamma_coef3`$coefficients[2,2]))
      }) %>% t(.)
    full_df[,c('beta.outcome','se.outcome')]<-sapply(1:p, function(i) {
        m <- lm(tmp2$y ~ tmp2$z[,i])
        m_summary <- summary(m)
        return(c(m$coefficients[-1], m_summary$coefficients[2,2]))
      }) %>% t(.)
    full_df[,c(paste0('beta.exposure.selection',1:3),
               paste0('se.exposure.selection',1:3))]<-sapply(1:p, function(i) {
          m <- lm(tmp3$x ~ tmp3$z[,i])
          m_summary <- summary(m)
          return(c(m$coef[2,],
                   m_summary$`Response gamma_coef1`$coefficients[2,2],
                   m_summary$`Response gamma_coef2`$coefficients[2,2],
                   m_summary$`Response gamma_coef3`$coefficients[2,2]))
    }) %>% t(.)
    full_df$relevant.ind<-0
    full_df$relevant.ind[1:s]<-1
    #####################
    for (i in 1:3) {
    full_df[,paste0('pval.exposure',i)]<-2*pnorm(abs(full_df[,paste0('beta.exposure',i)])/full_df[,paste0('se.exposure',i)],lower.tail = FALSE)
    full_df[,paste0('pval.selection',i)]<-2*pnorm(abs(full_df[,paste0('beta.exposure.selection',i)])/full_df[,paste0('se.exposure.selection',i)],lower.tail = FALSE)
    full_df[,paste0('z.exposure',i)]<-full_df[,paste0('beta.exposure',i)]/full_df[,paste0('se.exposure',i)]
    }
    return(full_df)
}
