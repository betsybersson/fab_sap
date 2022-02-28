#library(sae)
##################################
## fab conformal algorithm
##################################
fab_conf_pred = function(Y,mu,tau2,ALPHA){
  # input data vector Y
  # input prior parameters mu, tau2
  # input desired confidence level alpha
  # return alpha level prediction interval
  
  # parameters
  N = length(Y)
  tau2_theta = 1/(1/tau2 + N + 1)
  
  # critical values
  sol1s = Y
  sol2s = (2*(mu/tau2 + sum(Y))*tau2_theta-Y)/
    (1-2*tau2_theta)
  
  # obtain final solution
  S = sort(c(sol1s,sol2s))
  k = floor(ALPHA*(N+1))
  int = c(S[k],S[2*N-k+1])
  
  if ( k == 0 ) {
    stop("ALPHA is set too low.")
  }
  
  return(list("bounds" = int, "prob" = (1-k/(N+1))*100))
  
}
##################################
##################################  
# distance to avg conformal algorithm
##################################
avg_conf_pred = function(z,ALPHA){
  # z is training data, alpha is error rate
  #
  # |a + bz*| <= |ci + di z^*| where i indexes training data
  #
  ## get a, b, ci, di
  N = length(z)
  sum_z = sum(z)
  avg_z = (N+1) * z
  
  a = sum_z
  b = -N
  ci = sum_z-avg_z # constant
  di = rep(1,N) # multiplied by zstar
  
  b1 = (a-ci)/(di-b)
  b2 = (-a-ci)/(di+b)
  
  # obtain final solution
  S = sort(c(b1,b2))
  k = floor(ALPHA*(N+1))
  int = c(S[k],S[2*N-k+1])
  
  if ( k == 0 ) {
    stop("ALPHA is set too low.")
  }
  
  return(list("bounds" = int, "prob" = (1-k/(N+1))*100))
}
##################################
##################################  
# obtain emp bayes estimates
##################################
EB_values = function(y,group,W = NA,X = NA){
  
  ybar = tapply(y,group,mean)
  n = table(group)
  ss = (n-1)*c(tapply(y,group,var))
  
  ## obtain mles of s2
  mll = function(XX){
    a = XX[1]
    b = XX[2]
    
    ap = (a+n-1)/2
    bp = (b+ss)/2
    
    out = -sum( ( (a/2)*log(b/2)-lgamma(a/2) ) - ( ap*log(bp)-lgamma(ap) ) )
    
    return(out)
  }
  init = c(1,1)
  ab_out = optim(init, mll)$par
  
  ## posterior modes for each s2
  a = ab_out[1]
  b = ab_out[2]
  
  ap = (a+n-1)/2
  bp = (b+ss)/2
  
  # est of var
  s2 = bp/(ap+1)
  
  ## spatial ests of mu, t2, rho, theta using plugin s2
  W_t = diag(1/rowSums(W)) %*% W #row standardize
  
  df = data.frame("radon" = c(ybar),"uran" = c(X),"var" = c(s2/n))
  resultML = eblupSFH(radon ~ uran, var, as.data.frame(W_t),
                           method="ML",data = df)
  
  theta.eb = resultML$eblup
  mu.eb = resultML$fit$estcoef[1]
  rho.eb = resultML$fit$spatialcorr
  tau2.eb = resultML$fit$refvar
  
  out = list(mu = mu.eb,
             tau2 = tau2.eb,
             rho = rho.eb,
             theta = theta.eb,
             s20 = ab_out[2]/(ab_out[1]+1),
             s2 = s2)

}
##################################  
##################################  
# get conditional moments from eb output
##################################
get_spatial_moments = function(ebayes.output,j,W.stand,X){
  
  # dimension
  p = ncol(W.stand)
  
  # extract params
  mu = unlist(ebayes.output[,j]$mu)
  rho = c(ebayes.output[,j]$rho)
  omega2 = c(ebayes.output[,j]$tau2)
  theta = c(ebayes.output[,j]$theta)
  
  # helper
  G = omega2 * solve((eye(p)-rho*t(W.stand))%*%(eye(p)-rho*(W.stand))) 
  R = G[j,-j] %*% solve(G[-j,-j])
  
  # output
  MU = c(X[j,,drop=F] %*% mu + R %*% (theta - X[-j,,drop=F] %*% mu))
  TAU2 = c(G[j,j] - R %*% G[-j,j])/ebayes.output[,j]$s20 
  
  return(list("mu" = MU,"tau2" = TAU2))
}
##################################  
##################################  
# EBLUP SAE CODE
# AUTHORS: Isabel Molina, Yolanda Marhuenda
## SOURCE: https://rdrr.io/cran/sae/src/R/eblupSFH.R
##################################
eblupSFH <-
  function(formula,vardir,proxmat,method="REML",MAXITER=150,PRECISION=0.0001,data)
  {
    result <- list(eblup=NA, 
                   fit=list(method=method, convergence=TRUE, iterations=0, estcoef=NA, 
                            refvar=NA, spatialcorr=NA, goodness=NA)
    ) 
    
    if (method!="REML" & method!="ML")
      stop(" method=\"",method, "\" must be \"REML\" or \"ML\".")
    
    namevar     <- deparse(substitute(vardir))
    if (!missing(data))
    {
      formuladata <- model.frame(formula,na.action = na.omit,data)
      X           <- model.matrix(formula,data)        
      vardir      <- data[,namevar]
    } else
    {
      formuladata <- model.frame(formula,na.action = na.omit)
      X           <- model.matrix(formula)        
    }
    y <- formuladata[,1]            
    
    if (attr(attributes(formuladata)$terms,"response")==1){
      textformula <- paste(formula[2],formula[1],formula[3])
      } else{
      textformula <- paste(formula[1],formula[2])}
    
    if (length(na.action(formuladata))>0){
      stop("Argument formula=",textformula," contains NA values.")}
    if (any(is.na(vardir))){
      stop("Argument vardir=",namevar," contains NA values.")}
    
    proxmatname <- deparse(substitute(proxmat))
    if (any(is.na(proxmat))){
      stop("Argument proxmat=",proxmatname," contains NA values.")}
    
    if (!is.matrix(proxmat)){
      proxmat <- as.matrix(proxmat)}
    
    nformula  <- nrow(X)   
    nvardir   <- length(vardir) 
    nproxmat  <- nrow(proxmat)
    if (nformula!=nvardir | nformula!=nproxmat){
      stop("   formula=",textformula," [rows=",nformula,"],\n", 
           "     vardir=",namevar," [rows=",nvardir,"] and \n",
           "     proxmat=",proxmatname," [rows=",nproxmat,"]\n",
           "  must be the same length.")}
    if (nproxmat!=ncol(proxmat)){
      stop(" Argument proxmat=",proxmatname," is not a square matrix [rows=",nproxmat,",columns=",ncol(proxmat),"].")}
    
    m <-length(y)  # Sample size or number of areas
    p <-dim(X)[2]  # Num. of X columns of num. of auxiliary variables (including intercept)
    Xt <-t(X)
    yt <-t(y)
    proxmatt<-t(proxmat)
    I <-diag(1,m)
    
    # Initialize vectors containing estimators of variance and spatial correlation
    par.stim <-matrix(0,2,1)
    stime.fin<-matrix(0,2,1)
    
    # Initialize scores vector and Fisher information matrix
    s<-matrix(0,2,1)
    Idev<-matrix(0,2,2)
    
    # Initial value of variance set to the mean of sampling variances vardir
    # Initial value of spatial correlation set to 0.5
    sigma2.u.stim.S<-0
    rho.stim.S<-0
    
    sigma2.u.stim.S[1]<-median(vardir)
    rho.stim.S[1]<-0.5
    
    if (method=="REML")   # Fisher-scoring algorithm for REML estimators start
    { 
      k<-0
      diff.S<-PRECISION+1
      while ((diff.S>PRECISION)&(k<MAXITER))
      {
        k<-k+1
        
        # Derivative of covariance matrix V with respect to variance
        derSigma<-solve((I-rho.stim.S[k]*proxmatt)%*%(I-rho.stim.S[k]*proxmat))
        
        # Derivative of covariance matrix V with respect to spatial correlation parameter
        derRho<-2*rho.stim.S[k]*proxmatt%*%proxmat-proxmat-proxmatt
        derVRho<-(-1)*sigma2.u.stim.S[k]*(derSigma%*%derRho%*%derSigma)
        
        # Covariance matrix and inverse covariance matrix
        V<-sigma2.u.stim.S[k]*derSigma+I*vardir
        Vi<-solve(V)
        
        # Matrix P and coefficients'estimator beta
        XtVi<-Xt%*%Vi
        Q<-solve(XtVi%*%X)
        P<-Vi-t(XtVi)%*%Q%*%XtVi
        b.s<-Q%*%XtVi%*%y
        
        # Terms involved in scores vector and Fisher information matrix
        PD<-P%*%derSigma
        PR<-P%*%derVRho
        Pdir<-P%*%y
        
        # Scores vector
        s[1,1]<-(-0.5)*sum(diag(PD))+(0.5)*(yt%*%PD%*%Pdir)
        s[2,1]<-(-0.5)*sum(diag(PR))+(0.5)*(yt%*%PR%*%Pdir)
        
        # Fisher information matrix
        Idev[1,1]<-(0.5)*sum(diag(PD%*%PD))
        Idev[1,2]<-(0.5)*sum(diag(PD%*%PR))
        Idev[2,1]<-Idev[1,2]
        Idev[2,2]<-(0.5)*sum(diag(PR%*%PR))
        
        # Updating equation
        par.stim[1,1]<-sigma2.u.stim.S[k]
        par.stim[2,1]<-rho.stim.S[k]
        
        stime.fin<-par.stim+solve(Idev)%*%s
        
        # Restricting the spatial correlation to (-0.999,0.999)
        if (stime.fin[2,1]<=-1)
          stime.fin[2,1] <- -0.999
        if (stime.fin[2,1]>=1)
          stime.fin[2,1] <- 0.999
        
        # Restricting the spatial correlation to (-0.999,0.999) and the variance to (0.0001,infty)
        #if ((stime.fin[2,1]<=0.999)&(stime.fin[2,1]>=-0.999)&(stime.fin[1,1]>0.0001)){
        sigma2.u.stim.S[k+1]<-stime.fin[1,1]
        rho.stim.S[k+1]<-stime.fin[2,1]
        diff.S<-max(abs(stime.fin-par.stim)/par.stim)
        #}else
        #{
        #   sigma2.u.stim.S[k+1]<-stime.fin[1,1]
        #    rho.stim.S[k+1]<-stime.fin[2,1]
        #    diff.S<-PRECISION/10
        #}
      } # End of while
      
    } else       # Fisher-scoring algorithm for ML estimators start
    {
      k<-0
      diff.S<-PRECISION+1
      while ((diff.S>PRECISION)&(k<MAXITER))
      {
        k<-k+1
        
        # Derivative of covariance matrix V with respect to variance
        derSigma<-solve((I-rho.stim.S[k]*proxmatt)%*%(I-rho.stim.S[k]*proxmat))
        
        # Derivative of covariance matrix V with respect to spatial correlation
        derRho<-2*rho.stim.S[k]*proxmatt%*%proxmat-proxmat-proxmatt
        derVRho<-(-1)*sigma2.u.stim.S[k]*(derSigma%*%derRho%*%derSigma)
        
        # Covariance matrix and inverse covariance matrix
        V<-sigma2.u.stim.S[k]*derSigma+I*vardir
        Vi<-solve(V)
        
        # Coefficients'estimator beta and matrix P
        XtVi<-Xt%*%Vi
        Q<-solve(XtVi%*%X)
        P<-Vi-t(XtVi)%*%Q%*%XtVi
        b.s<-Q%*%XtVi%*%y
        
        # Terms involved in scores vector and Fisher information matrix
        PD<-P%*%derSigma
        PR<-P%*%derVRho
        Pdir<-P%*%y
        ViD<-Vi%*%derSigma
        ViR<-Vi%*%derVRho
        
        # Scores vector
        s[1,1]<-(-0.5)*sum(diag(ViD))+(0.5)*(yt%*%PD%*%Pdir)
        s[2,1]<-(-0.5)*sum(diag(ViR))+(0.5)*(yt%*%PR%*%Pdir)
        
        # Fisher information matrix
        Idev[1,1]<-(0.5)*sum(diag(ViD%*%ViD))
        Idev[1,2]<-(0.5)*sum(diag(ViD%*%ViR))
        Idev[2,1]<-Idev[1,2]
        Idev[2,2]<-(0.5)*sum(diag(ViR%*%ViR))
        
        # Updating equation      
        par.stim[1,1]<-sigma2.u.stim.S[k]
        par.stim[2,1]<-rho.stim.S[k]
        
        stime.fin<-par.stim+solve(Idev)%*%s
        
        # Restricting the spatial correlation to (-0.999,0.999)
        if (stime.fin[2,1]<=-1)
          stime.fin[2,1] <- -0.999
        if (stime.fin[2,1]>=1)
          stime.fin[2,1] <- 0.999
        
        # Restricting the spatial correlation to (-0.999,0.999) and the variance to (0.0001,infty)
        #if ((stime.fin[2,1]<=0.999)&(stime.fin[2,1]>=-0.999)&(stime.fin[1,1]>0.0001)){
        sigma2.u.stim.S[k+1]<-stime.fin[1,1]
        rho.stim.S[k+1]<-stime.fin[2,1]
        diff.S<-max(abs(stime.fin-par.stim)/par.stim)
        #}else
        #{
        #   sigma2.u.stim.S[k+1]<-stime.fin[1,1]
        #   rho.stim.S[k+1]<-stime.fin[2,1]
        #   diff.S<-PRECISION/10
        #}
        
      } # End of while
      
    }   
    
    # Final values of estimators
    if (rho.stim.S[k+1]==-0.999){
      rho.stim.S[k+1] <- -1
      }else if (rho.stim.S[k+1]==0.999){
      rho.stim.S[k+1] <- 1}
    rho <-rho.stim.S[k+1]
    
    sigma2.u.stim.S[k+1]<-max(sigma2.u.stim.S[k+1],0)
    sigma2u <- sigma2.u.stim.S[k+1]
    
    #print(rho.stim.S)
    #print(sigma2.u.stim.S)
    
    # Indicator of convergence
    result$fit$iterations  <- k  
    if(k>=MAXITER && diff.S>=PRECISION) 
    {
      result$fit$convergence <- FALSE
      #return(result)
    }
    
    result$fit$refvar       <- sigma2u
    result$fit$spatialcorr <- rho
    if (sigma2u<0 || rho<(-1) || rho>1 )  # COMPROBAR
    {
      print("eblupSFH: este mensaje no debe salir")
      return(result)
    }
    # Computation of the coefficients'estimator (Bstim)
    A   <-solve((I-rho*proxmatt)%*%(I-rho*proxmat))    
    G   <-sigma2u*A
    V   <-G+I*vardir
    Vi  <-solve(V)
    XtVi<-Xt%*%Vi
    Q   <-solve(XtVi%*%X)
    Bstim<-Q%*%XtVi%*%y
    
    # Significance of the regression coefficients
    std.errorbeta<-sqrt(diag(Q))
    tvalue       <-Bstim/std.errorbeta
    pvalue       <-2*pnorm(abs(tvalue),lower.tail=FALSE)
    coef         <-data.frame(beta=Bstim,std.error=std.errorbeta,tvalue,pvalue)
    
    # Goodness of fit measures: loglikelihood, AIC, BIC
    Xbeta <-X%*%Bstim
    resid <-y-Xbeta
    loglike<-(-0.5)*(m*log(2*pi)+determinant(V,logarithm=TRUE)$modulus+t(resid)%*%Vi%*%resid)
    AIC    <-(-2)*loglike+2*(p+2)
    BIC    <-(-2)*loglike+(p+2)*log(m)
    goodness<-c(loglike=loglike,AIC=AIC,BIC=BIC)
    
    # Computation of the Spatial EBLUP
    res<-y-X%*%Bstim
    thetaSpat<-X%*%Bstim+G%*%Vi%*%res
    
    result$fit$estcoef  <- coef
    result$fit$goodness <- goodness
    result$eblup       <- thetaSpat
    
    return(result)
  }
##########################################################
## Identity
# N: length of diagonal
##########################################################
eye = function(N){
  diag(rep(1,N))
}
##########################################################
## noncentral, non standard t density evaluated at a point X
# X: point to evaluate density at
# nu: number of degrees of freedom
# mu: mean of t distribution
# sig2: variance of t distribution
## density is written down on STA 721 midterms in my classes
##########################################################
dt_ncns = function(X,nu,mu,sig2,propto = TRUE){
  # density at X for t distribution with 
  # nu degrees of freedom
  # mean mu
  # variance sig2
  # if propto = true, then return density up to normalizing constant
  
  out = gamma((nu+1)/2) * (nu * sig2 * pi)^(-1/2) / gamma(nu/2) *
    (1 + (X - mu)^2/(sig2 * nu))^(-(nu+1)/2) 
  
  if (propto == TRUE){
    out = out/(gamma((nu+1)/2)/ gamma(nu/2) /sqrt(nu*pi))
  }
  
  return(out)
}
##########################################################
## posterior parameters of NIG conjugate model 
##########################################################
NIG_post = function(Y,mu,tau2,alpha,beta,Y2){
  # for data Y and hyperparatemets mu, tau2, beta, alpha,
  # return posterior parameters of NIG distribution
  n = length(Y)
  tau2_theta = 1/(1/tau2 + n)
  mu_theta = (mu/tau2 + n*mean(Y)) * tau2_theta
  alpha_sigma = (alpha + n)/2
  beta_sigma = (beta + sum(Y2) + mu^2/tau2 - mu_theta^2/tau2_theta)/2
  
  return(list("mu" = mu_theta,"tau2" = tau2_theta,
              "alpha" = alpha_sigma, "beta" = beta_sigma))
  
}
##########################################################
## posterior predictive parameters of NIG conjugate model 
##########################################################
t_post_pred = function(Y,mu,tau2,alpha,beta,Y2 = Y^2){
  # for data Y and hyperparatemets mu, tau2, beta, alpha,
  # return parameters ofposterior predictive t distn
  
  params = NIG_post(Y,mu,tau2,alpha,beta,Y2)
  
  nu_t = params$alpha * 2
  mu_t = params$mu
  var_t = params$beta / params$alpha *(1 + params$tau2)
  
  return(list("nu" = nu_t,"mu" = mu_t, "sig2" = var_t))
  
}
