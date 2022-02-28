################################
### Frequentist and Bayesian Coverage under accurate and inaccurate parametric assumptions
################################
## pred and EB functions

oracle.pred = function(theta,sd,n,alpha){
  
  zval = qnorm(1-alpha/2)
  
  ci_help = zval * sd * sqrt(1/n+1)
  
  out = c(theta - ci_help,theta + ci_help)
  
  return(out)
}

oracle.bayes = function(theta,sd,n,mu,tau2,alpha){
  
  sig2 = sd^2
  
  varj = 1 / (1/tau2 + n/sig2)
  thetaj = (mu/tau2 + theta * n/sig2) * varj
  
  zval = qnorm(1-alpha/2)
  
  ci_help = zval * sqrt(varj + sig2)
  
  out = c(thetaj - ci_help,thetaj + ci_help)
  
  return(out)
}
pivot_pred_knownvar = function(z, alpha, sd = 1){
  # standard prediction for single group
  n = length(z)
  zbar = mean(z)
  
  tval = qnorm(1-alpha/2)
  
  ci_help = tval * sd * sqrt(1/n+1)
  
  out = c(zbar - ci_help,zbar + ci_help)
  
  return(out)
}
bayes_pred_knownvar = function(Y,mu,tau2,ALPHA,sd = 1){
  
  n = length(Y)
  ybar = mean(Y)
  
  sig2 = sd^2
  
  varj = 1 / (1/tau2 + n/sig2)
  thetaj = (mu/tau2 + ybar * n/sig2) * varj
  
  tval = qnorm(1-ALPHA/2)
  
  ci_help = tval * sqrt(varj + sig2)
  
  out = c(thetaj - ci_help,thetaj + ci_help)
  
  return(out)
}
################################
## obtain coverage percentages for 
## classical pivot and bayesian prediction methods
################################

## model parameters
theta.range = seq(from=-8,to=8,by=.1);n_T = length(theta.range)
sig2 = 1; sd = sqrt(sig2)
mu = 0
tau2 = 1/2

n = 3
alpha = .25

## get exact coverage when Normal model assumptions are accurate

# frequentist coverage- constant over theta.range
loc = 1
pred.b = oracle.pred(theta.range[loc],1,n,alpha)
freq.cons.cov = diff(pnorm(pred.b,theta.range[loc],sd*sqrt(1/n+1)))

# bayes coverage
bayes.var = 1 / (1/tau2 + n/sig2)
helper.var = (n/sig2)*bayes.var^2 + sig2

bayes.cov = array(NA,dim=c(n_T))
for ( j in 1:n_T){
  theta = theta.range[j]
  pred.b = oracle.bayes(theta,1,n,mu,tau2,alpha)
  bayes.cov[j] = diff(pnorm(pred.b,theta,sqrt(helper.var)))
}

## obtain coverage under incorrect param assumption

# model assumptions- two normal distns with mean theta, variance 1
fx = function(n,theta = 0,c = 1){
  wt = rbinom(n,1,.5)
  tight.sd = .000001
  wt * rnorm(n,(theta - c),tight.sd) + (1-wt) * rnorm(n,(theta + c),tight.sd)
}

# MC approx coverage

S = 1000000
bounds.b = bounds.f = array(NA,dim=c(S,2))
bayes.wrong.cov = freq.wrong.cov = array(NA,dim=c(n_T))
for ( j in 1:n_T){
  theta = theta.range[j]
  for (s in 1:S){
    y = fx(n,theta)
    bounds.f[s,] = pivot_pred_knownvar(y,alpha)
    bounds.b[s,] = bayes_pred_knownvar(y,mu,tau2,alpha,sd = 1)
  }
  ystar = fx(S,theta)
  freq.wrong.cov[j] = mean((ystar>bounds.f[,1])&(ystar<bounds.f[,2]))
  bayes.wrong.cov[j] = mean((ystar>bounds.b[,1])&(ystar<bounds.b[,2]))
  
}

#pdf("./plots/coverage_muline.pdf",family="Times",height = 7,width = 6.6)
plot(theta.range,bayes.cov*100,
     type= "l", col = "red",
     ylab = "Frequentist Coverage Percent",
     xlab = expression(theta),
     ylim=range(bayes.cov*100,freq.wrong.cov*100),
     cex.axis=1.5,cex.lab=1.5,
     yaxt = "n")
axis(side = 2, at = seq(from = 15, to = 75, by = 15),
     cex.axis=1.5,cex.lab=1.5)
abline(h=freq.cons.cov*100,col="blue")
abline(v = 0,lwd = 3)
lines(theta.range,freq.wrong.cov*100,col="blue",lty="dashed")
lines(theta.range,bayes.wrong.cov*100,col="red",lty="dashed")
#dev.off()

