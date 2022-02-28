source("prediction_functions.R")

S = 25000

#########################################
## vary theta-mu; tau2 = 1/2,2
#########################################

alpha = .25
## MODEL PARAMS
s2 = 1
MU.VAL = seq(from = -2.5, to = 2.5, by = .05); n_M = length(MU.VAL)
NS = c(3,7,11,15,19); n_N = length(NS)

# grand sim to collect MC estimates of mean for each n, t2 pair for 2 methods
mc_mean_out = array(NA,dim=c(n_N,n_M,3))
for ( j in 1:n_N ){
  
  n = NS[j]
  
  interval_out = array(NA,dim=c(2,S,n_M,3)) # fab 1/2, fab 2, avg
  for ( s in 1:S ){
    
    for ( T in 1:n_M){
      
      theta = MU.VAL[T]
      Y = rnorm(n,theta,sqrt(s2))
      
      # get interval
      interval_out[,s,T,1] = fab_conf_pred(Y,0,1/2,alpha)$bounds
      interval_out[,s,T,2] = fab_conf_pred(Y,0,2,alpha)$bounds
      interval_out[,s,T,3] = avg_conf_pred(Y,alpha)$bounds
    }
    
    
  }
  
  temp = apply(interval_out,c(2,3,4),diff)
  mc_mean_out[j,,] = apply(temp,c(2,3),mean)# rows are t2; columns are method
  
  if ( j == 1 ){
    bounds = apply(interval_out,c(1,3,4),mean)
  }
  
}


cols = gray.colors(n_N)

### TAU2 = 1/2
pdf("./single_pop_varymu_tau_HALF.pdf",family="Times",height = 7,width = 6.6)
for ( j in 1:n_N){
  
  mean.mat = mc_mean_out[j,,1]/mc_mean_out[j,,3]
  
  if(j==1){
    plot(MU.VAL,mean.mat,
         xlab = expression(theta), ylab = "Expected Interval Width Ratio",
         ylim = c(range(mean.mat,1.01)),
         type="l",cex.axis=1.5,cex.lab=1.5)
    title("(c)", adj = .05, line = 1,cex.main=2)
    abline(h=1,lty="dashed")
    
    min.val = min(mean.mat)
  }
  
  lines(MU.VAL,(mean.mat),lwd=2,col=cols[j])
  
  
}
dev.off()


### TAU2 = 2
pdf("./single_pop_varymu_tau_2.pdf",family="Times",height = 7,width = 6.6)
for ( j in 1:n_N){
  
  mean.mat = mc_mean_out[j,,2]/mc_mean_out[j,,3]
  
  
  if(j==1){
    plot(MU.VAL,mean.mat,
         xlab = expression(theta), ylab = "Expected Interval Width Ratio",
         ylim = c(range(mean.mat,1.01)),
         type="l",cex.axis=1.5,cex.lab=1.5)
    title("(d)", adj = .05, line = 1,cex.main=2)
    abline(h=1,lty="dashed")
  }
  
  lines(MU.VAL,(mean.mat),lwd=2,col=cols[j])
  
  
}
dev.off()

## plot pred interval bounds
pdf("./single_pop_varymu_bounds.pdf",family="Times",height = 7,width = 6.6)
make.mat = cbind(bounds[1,,-2],bounds[2,,-2])
cols = c("grey","black")
matplot(MU.VAL,(make.mat),
        xlab = expression(theta), ylab = "Prediction Interval Endpoints",
        col = rep(cols,2), type="l",lwd=c(3,1.5,3,1.5),
        lty = 1,cex.axis=1.5,cex.lab=1.5)
title("(a)", adj = .05, line = 1,cex.main=2)
abline(a=0,b=1,lty="dashed")
abline(v=0,lty="dashed")
dev.off()


##############################
##############################
##############################
##############################

#########################################
## vary tau2; bayes risk; NOT a function of MU if mu=MU
#########################################


## model params
MU = 0; mu = MU
T2 = seq(from = .05, to = 6, by = .05); n_T = length(T2)

# grand sim to collect MC estimates of mean for each n, t2 pair for 2 methods
mc_mean_out = array(NA,dim=c(n_N,n_T,5))
for ( j in 1:n_N ){
  
  n = NS[j]
  
  interval_out = array(NA,dim=c(2,S,n_T,5)) # fab, avg
  for ( s in 1:S ){
    
    # hyperparams
    
    for ( T in 1:n_T){
      
      tau2 = T2[T]
      
      theta = rnorm(1,MU,sqrt(tau2)*sqrt(s2))
      Y = rnorm(n,theta,sqrt(s2))
      
      # get interval 
      interval_out[,s,T,1] = fab_conf_pred(Y,mu,tau2,alpha)$bounds
      interval_out[,s,T,2] = avg_conf_pred(Y,alpha)$bounds
    }
    
  }
  
  temp = apply(interval_out,c(2,3,4),diff)
  mc_mean_out[j,,] = apply(temp,c(2,3),mean)# rows are t2; columns are method
  
  
}


cols = gray.colors(n_N)
pdf("./single_pop_bayesrisk_varytau.pdf",family="Times",height = 7,width = 6.6)
for ( j in 1:n_N){ #n_N
  
  mean.mat = mc_mean_out[j,,1]/mc_mean_out[j,,2]
  
  if(j==1){
    plot(T2,mean.mat,
         xlab = expression(tau^2), ylab = "Bayes Expected Interval Width Ratio",
         ylim = c(range(mean.mat,1.01)),
         pch="",cex.axis=1.5,cex.lab=1.5)
    title("(b)", adj = .05, line = 1,cex.main=2)
    abline(h=1,lty="dashed")
  }
  lo = smooth.spline(T2,mean.mat)
  #lines(T2,(mean.mat),lwd=3,col=cols[j])
  lines(lo,lwd=3,col=cols[j])
  
}
dev.off()

