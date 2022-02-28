source("prediction_functions.R")
library(dplyr)

## problem set-up 

# params
ALPHA = .2
n = 4
mu = 0
tau2 = 4

# get dummy data
set.seed(1234)
y = rnorm(n,0,1)
set.seed(Sys.time())

## get conformity scores over yn+1

# obtain conformity score of c n+1
params = t_post_pred(y,mu,tau2,1,1)

y_star = seq(from = params$mu-5, 
             to = params$mu+5, by = .005)

c_star = sapply(y_star,function(j)
  dt_ncns(j,params$nu,params$mu,params$sig2) )


ci_func = function(potential_y_star,I){
  temp_Y = c(y[-I],potential_y_star)
  params2 = t_post_pred(temp_Y,mu,tau2,1,1)
  ci_out = dt_ncns(y[I],params2$nu,params2$mu,params2$sig2)
  return(ci_out)
}

cis = sapply(1:n, function(i)
  sapply(y_star,function(j)ci_func(j,i)))


#pdf("./plots/computation_ex.pdf",family="Times",height = 7,width = 6.6)
plot(y_star,c_star,type="l",
     lwd = 3,
     xlab = expression('Y'[n+1]),ylab="",
     xlim = c(min(y_star)-.8,max(y_star)),
     ylim = c(0,.75), yaxt="n",cex.axis=1.5,cex.lab=1.5)
title("(a)", adj = .05, line = 1,cex.main=2)
for( i in 1:n){
  lines(y_star,cis[,i],lty=i+1,col="red")
}
text(min(y_star)-c(.6,.4,.4,.4,.4),c(c_star[1],cis[1,]+c(.02,0,0,-.02)),
     c(expression(c[n+1]),expression(c[1]),
       expression(c[2]),expression(c[3]),
       expression(c[4])),cex=1.8)
abline(v=params$mu,lwd = 3)
#dev.off()


## get num of ci < cn+1 graph

# obtain critical values
tau2_theta = 1/(1/tau2 + n + 1)
sol1s = y
sol2s = (2*(mu/tau2 + sum(y))*tau2_theta-y)/
  (1-2*tau2_theta)
S.uns = c(sol1s,sol2s); S.uns.mat = cbind(sol1s,sol2s)[n:1,]


#pdf("./plots/computation_cs_ex.pdf",family="Times",height = 7,width = 6.6)
plot(y_star,c_star,
     pch="",ylab="",lwd = 3,
     xlab = expression('Y'[n+1]),
     xlim = c(min(y_star)-.1,max(y_star)),
     ylim = c(-.2,4.2),cex=1.5,
     yaxt = "n",cex.axis=1.5,cex.lab = 1.5)
title("(b)", adj = .05, line = 1.,cex.main=2)
lines(c(min(y_star)-.1,max(y_star)),c(0,0))
for ( j in 1:4){
  lines(S.uns.mat[j,],rep(j,2))
  
}
points(c(S[1,1],S[9,2]),c(0,0),pch=c("<",">"),cex=1.5)
points(c(S.uns.mat),rep(1:4,2),pch=16,cex=1.5)
abline(v=params$mu,lwd = 3)
axis(2,at = 0:4,labels = c(expression(S[n+1]),expression(S[4]),
                           expression(S[3]),expression(S[2]),
                           expression(S[1])),
     las = 2,cex.axis=1.5)
#dev.off()

## get py graph

# order critical values
S = rep(sort(S.uns),each=2)
S = matrix(c(min(y_star)-.1,S,max(y_star)),
           ncol = 2, byrow = T )
NUM = c(1:(n+1),(n:1))/(n+1)


## plot
#pdf("./plots/computation_num_ex.pdf",family="Times",height = 7,width = 6.6)
plot(y_star,c_star,
     pch="",lwd = 3,
     xlab = expression('Y'[n+1]),
     ylab = expression(p[y]),
     xlim = c(min(y_star)-.1,max(y_star)),
     ylim = c(0.9/5,1),cex.axis=1.5,cex.lab=1.5)
title("(c)", adj = .05, line = 1,cex.main=2)
for ( j in 1:9){
  lines(S[j,],rep(NUM[j],2))
  if (j<5){
    if (j > 1){
      points(S[j,1],NUM[j],pch=16,cex=1.5) }
    if (j < 9){
      points(S[j,2],NUM[j],pch=1,cex=1.5)
      lines(rep(S[j,2],2),c(NUM[j],NUM[j+1]),lty=3) }
  } else if (j == 5){
    points(S[j,1],NUM[j],pch=16,cex=1.5)
    points(S[j,2],NUM[j],pch=16,cex=1.5)
    lines(rep(S[j,2],2),c(NUM[j],NUM[j+1]),lty=3)
  } else {
    points(S[j,1],NUM[j],pch=1,cex=1.5)
    if (j < 9){
      points(S[j,2],NUM[j],pch=16,cex=1.5)
      lines(rep(S[j,2],2),c(NUM[j],NUM[j+1]),lty=3) }
  }
  
}
# end points
points(S[1,1],NUM[1],pch="<")
points(S[9,2],NUM[9],pch=">")
abline(v=params$mu,lwd = 3)
#dev.off()


