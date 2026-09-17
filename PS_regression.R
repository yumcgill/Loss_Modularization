library(gtools)
library(parallel)
library(rstanarm)
options(mc.cores = 23)

stick.breaking<-function(av,Nv=1000){
  u<-rbeta(Nv,1,av)	
  v<-u
  w<-c(1,cumprod(1-u[-Nv]))
  return(v*w)
}
al<-5
Nv <- 1000

causal_mod<-function(seed,N,alp_true){
  set.seed(seed)
  x1<-rnorm(n=N,mean=1,sd=1)
  x2<-rnorm(n=N,mean=1,sd=1) 
  x3<-rnorm(n=N,mean=-1,sd=1)
  x4<-rnorm(n=N,mean=-1,sd=1)
  X<-cbind(1,x1,x2,x3,x4)
  
  px<-1/(1+exp(-c(X%*%alp_true)))
  D<-rbinom(n=N,size=1,prob=px)
  # Now generate the continuous outcome
  beta0<-0.25
  beta1<-0.25
  beta2<-0.25
  beta3<-0.25
  beta4<-1.5
  linp<-beta0*x1+beta1*x2+beta2*x3+beta3*x4+beta4*x3*x4
  sigma_y<-1
  Y<-rnorm(N,mean=linp,sd=sigma_y)
  sim.data<-data.frame(x1=x1,x2=x2,x3=x3,x4=x4,D=as.integer(D),Y=Y)
  ps_dp <- glm(D~x1+x2+x3+x4,data=sim.data,family = binomial(link = "logit"))$fitted.values
  sim.data$pred<-predict(lm(Y ~ D+ps_dp))
  gen_ps <- stan_glm(D~x1+x2+x3+x4,data=sim.data,family = binomial(link = "logit"),chain = 1, warmup = 500, iter = 1500)
  gen_ps_post <- posterior_epred(gen_ps)
  post_theta_BB<-NULL
  post_theta_DP<-NULL
  post_theta_GEN<-NULL
  for (i in 1:1000){
   
    ps_i      <- gen_ps_post[i, ]
    gen_design <- cbind(1, D, ps_i)
    gen_fit    <- lm.fit(x = gen_design, y = Y)
    post_theta_GEN[i] <- gen_fit$coefficients[2]
    
    w<-as.vector(rdirichlet(1,rep(1,N)))
    ps_bb <- glm(D~x1+x2+x3+x4,weights=w,data= sim.data,family = binomial(link = "logit"))$fitted.values
    mod<-lm(Y ~ D+ps_bb, weights=w,data= sim.data)
    post_theta_BB[i] <- mod$coef[2]
    u<-runif(Nv)
    datasetnew<-sim.data
    for (nv in 1:Nv) {
      if(u[nv] > al/(al+N+nv-1)){
        ind<-sample(1:(N+nv-1),size=1)
        datasetnew<-rbind(datasetnew,datasetnew[ind,])
      }else{
        ind<-sample(1:(N+nv-1),size=1)
        newdata<-datasetnew[ind,]
        newdata$Y<-rnorm(1,newdata$pred,1)
        datasetnew<-rbind(datasetnew,newdata)
      }
    }
    
    
    resamp.data <- datasetnew[-c(1:N),]
    resamp.data$ps <-glm(D ~ x1+x2+x3+x4,data=resamp.data,family = binomial(link = "logit"))$fitted.values
    mod_dp <- lm(Y ~ D+ps,data=resamp.data)
    post_theta_DP[i] <- mod_dp$coef[2]
    
    
  } 
  
  

  postmean_theta_GEN<-mean(post_theta_GEN)
  ci_GEN<-0
  if(quantile(post_theta_GEN,prob=c(0.025))<0 && quantile(post_theta_GEN,prob=c(0.975))>0) {ci_GEN<-1}
  postmean_theta_BB<-median(post_theta_BB)
  ci_BB<-0
  if(quantile(post_theta_BB,prob=c(0.025))<0 && quantile(post_theta_BB,prob=c(0.975))>0) {ci_BB<-1}
  postmean_theta_DP<-median(post_theta_DP)
  ci_DP<-0
  if(quantile(post_theta_DP,prob=c(0.025))<0 && quantile(post_theta_DP,prob=c(0.975))>0) {ci_DP<-1}
  return(list(ate_gen=postmean_theta_GEN,ci_GEN=ci_GEN, ate_bb=postmean_theta_BB,ci_BB=ci_BB,ate_dp=postmean_theta_DP,ci_DP=ci_DP))
}

alp<-c(0.00,0.30,0.80,0.30,0.80)
set.seed(12)
seeds <- sample(c(1:100000000),1000)
modular_causal_50<-mclapply(seeds,function(x) causal_mod(x,50,alp))
ate_gen_50<-unlist(mclapply(modular_causal_50, '[[', "ate_gen"))
coverge_gen_50<-sum(unlist(mclapply(modular_causal_50, '[[', "ci_GEN")))/1000
mean(ate_gen_50);var(ate_gen_50);coverge_gen_50


ate_BB_50<-unlist(mclapply(modular_causal_50, '[[', "ate_bb"))
coverge_BB_50<-sum(unlist(mclapply(modular_causal_50, '[[', "ci_BB")))/1000
mean(ate_BB_50);var(ate_BB_50);coverge_BB_50


ate_DP_50<-unlist(mclapply(modular_causal_50, '[[', "ate_dp"))
coverge_DP_50<-sum(unlist(mclapply(modular_causal_50, '[[', "ci_DP")))/1000
mean(ate_DP_50);var(ate_DP_50);coverge_DP_50



modular_causal_100<-mclapply(seeds,function(x) causal_mod(x,100,alp))
ate_gen_100<-unlist(mclapply(modular_causal_100, '[[', "ate_gen"))
coverge_gen_100<-sum(unlist(mclapply(modular_causal_100, '[[', "ci_GEN")))/1000
mean(ate_gen_100);var(ate_gen_100);coverge_gen_100


ate_BB_100<-unlist(mclapply(modular_causal_100, '[[', "ate_bb"))
coverge_BB_100<-sum(unlist(mclapply(modular_causal_100, '[[', "ci_BB")))/1000
mean(ate_BB_100);var(ate_BB_100);coverge_BB_100


ate_DP_100<-unlist(mclapply(modular_causal_100, '[[', "ate_dp"))
coverge_DP_100<-sum(unlist(mclapply(modular_causal_100, '[[', "ci_DP")))/1000
mean(ate_DP_100);var(ate_DP_100);coverge_DP_100



modular_causal_30<-mclapply(seeds,function(x) causal_mod(x,30,alp))
ate_gen_30<-unlist(mclapply(modular_causal_30, '[[', "ate_gen"))
coverge_gen_30<-sum(unlist(mclapply(modular_causal_30, '[[', "ci_GEN")))/1000
mean(ate_gen_30);var(ate_gen_30);coverge_gen_30


ate_BB_30<-unlist(mclapply(modular_causal_30, '[[', "ate_bb"))
coverge_BB_30<-sum(unlist(mclapply(modular_causal_30, '[[', "ci_BB")))/1000
mean(ate_BB_30);var(ate_BB_30);coverge_BB_30


ate_DP_30<-unlist(mclapply(modular_causal_30, '[[', "ate_dp"))
coverge_DP_30<-sum(unlist(mclapply(modular_causal_30, '[[', "ci_DP")))/1000
mean(ate_DP_30);var(ate_DP_30);coverge_DP_30
