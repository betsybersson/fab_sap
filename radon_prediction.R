source("prediction_functions.R")
##############################################
## Prediction on Radon data set
##############################################

## load data
source("org_radon_data.R")

## error rate
alpha = 1/3

##############################################
#### obtain intervals
##############################################

## get hyperparamter estimates
spatial_ests = sapply(1:p,function(j) EB_values(y[j!=group],rep(1:(p-1),n[-j]),
                                                             W[-j,-j],u[-j]))


## get prediction intervals
opt_conf_spatial = array(NA,dim = c(2,p))
for (j in 1:p){
  
  ## SPATIAL MODEL WITH COVARIATES
  cov.out = get_spatial_moments(spatial_ests,j,W.stand,cbind(1,u))
  
  muj = cov.out$mu
  tau2j = cov.out$tau2

  ## GET OPT CONF INTERVALS
  opt_conf_spatial[,j] = fab_conf_pred(y[group == j],
                                       muj,tau2j,
                                       alpha)$bounds
}


## conformal prediction, no info sharing
avg_conf = sapply(1:p,function(j)
  avg_conf_pred(y[group == j],alpha)$bounds )

## plot 
means = tapply(y,group,mean)

plot(rep(means,each = 3),rbind(avg_conf,NA),
     col = "gray",lwd = 5,type = "l",
     ylab = expression(y[n+1]),
     xlab = expression(bar(y)[j]),
     ylim = range(opt_conf_spatial,avg_conf))
lines(rep(means,each = 3),rbind(opt_conf_spatial,NA),col = "black")
abline(v=mean(means),col="black",lty="dashed")
abline(h=mean(means),col="black",lty="dashed")
abline(b=1,a=0,col="black",lty="dashed")

## overall comparison
opt.spat.width = diff(opt_conf_spatial)
avg.width = diff(avg_conf)

mean(avg.width>opt.spat.width)




