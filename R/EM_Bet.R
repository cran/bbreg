#############################################################################################
#' @title infmat_bet
#' @description Function to compute standard errors based on the Fisher information matrix for the beta regression.
#' @param theta vector of parameters (all coefficients: kappa and lambda).
#' @param z response vector with 0 < z[i] < 1.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @return Vector of standard errors.
infmat_bet = function(theta,z,x,v)
{
  n = length(z)
  nkap = ncol(x)
  kap = theta[1:nkap]
  lam = theta[-c(1:nkap)]; nlam = length(lam)
  mu = exp(x%*%kap)/(1+exp(x%*%kap))
  phi = exp(v%*%lam)
  # E[-d2/d.kap.kap]
  d2.k.k = array(0,c(nkap,nkap))
  aux1 = trigamma(mu*phi)+trigamma((1-mu)*phi)
  aux1 = (phi^2)*aux1*(mu^2)*((1-mu)^2)
  aux2 = log(z)-log(1-z)-digamma(mu*phi)+digamma((1-mu)*phi)
  aux2 = aux2*phi*mu*(1-mu)*(1-2*mu)
  for(j in 1:nkap){
    for(l in 1:nkap){
        d2.k.k[j,l] = sum(aux1*x[,j]*x[,l]-aux2*x[,j]*x[,l])
    }
  }
  # E[-d2/d.lam.lam]
  d2.l.l = array(0,c(nlam,nlam))
  aux1 = (mu^2)*trigamma(mu*phi)+((1-mu)^2)*trigamma((1-mu)*phi)
  aux1 = aux1*(phi^2)
  aux2 = mu*(log(z)-log(1-z))+digamma(phi)+log(1-z)-mu*digamma(mu*phi)
  aux2 = aux2 -(1-mu)*digamma((1-mu)*phi)
  aux2 = aux2*phi
  for(j in 1:nlam){
    for(l in 1:nlam){
        d2.l.l[j,l] = sum(aux1*v[,j]*v[,l]-aux2*v[,j]*v[,l])
    }
  }
  # E[-d2/d.kap.lam]
  d2.k.l = array(0,c(nkap,nlam));
  aux1 = -log(z)+log(1-z)+digamma(mu*phi)-digamma((1-mu)*phi)
  aux1 = aux1 +mu*phi*trigamma(mu*phi) -(1-mu)*phi*trigamma((1-mu)*phi)
  aux1 = aux1*mu*(1-mu)*phi
  for(j in 1:nkap){
    for(l in 1:nlam){
        d2.k.l[j,l] = sum(aux1*x[,j]*v[,l])
    }
  }
  d2.Q1 = cbind(d2.k.k, d2.k.l)
  d2.Q2 = cbind(t(d2.k.l), d2.l.l)
  d2.Q = rbind(d2.Q1, d2.Q2)
  #
  # E[d/d.kap.j d/d.kap.l] and E[d/d.kap.j d/d.lam.l]
  dk.dk = array(0,c(nkap,nkap))
  dk.dl = array(0,c(nkap,nlam))
  aux1 = log(z)-log(1-z)-digamma(mu*phi)+digamma((1-mu)*phi)
  aux1 = aux1*phi*mu*(1-mu)
  aux2 = mu*(log(z)-log(1-z))+digamma(phi)+log(1-z)
  aux2 = aux2-mu*digamma(mu*phi)-(1-mu)*digamma((1-mu)*phi)
  aux2 = aux2*phi
  for(j in 1:nkap){
    for(l in 1:nkap){ dk.dk[j,l] = sum(aux1*x[,j])*sum(aux1*x[,l]) }
    for(l in 1:nlam){ dk.dl[j,l] = sum(aux1*x[,j])*sum(aux2*v[,l]) }
  }
  # E[d/d.lam.j d/d.lam.l]
  dl.dl = array(0,c(nlam,nlam))
  aux1 = mu*(log(z)-log(1-z))+digamma(phi)+log(1-z)-mu*digamma(mu*phi)-(1-mu)*digamma((1-mu)*phi)
  aux1 = aux1*phi
  aux2 = digamma(phi)+log(1-z)+mu*(log(z)-log(1-z))
  aux2 = aux2-mu*digamma(mu*phi)-(1-mu)*digamma((1-mu)*phi)
  aux2 = (trigamma(phi)+aux2^2)*(phi^2)
  for(j in 1:nlam){
    for(l in 1:nlam){
      dl.dl[j,l] = sum(aux2*v[,j]*v[,l])
      for(i in 1:n){
        for(k in 1:n){
          if(k != i){ dl.dl[j,l] = dl.dl[j,l] + aux1[i]*aux1[k]*v[i,j]*v[k,l] }
        }
      }
    }
  }
  #
  dd.Q1 = cbind(dk.dk, dk.dl)
  dd.Q2 = cbind(t(dk.dl), dl.dl)
  dd.Q = rbind(dd.Q1, dd.Q2)
  #
  aux = d2.Q-dd.Q # Fisher Information Matrix.
  inv.aux = tryCatch(solve(aux), error = function(e) rep(NA,nrow(aux)))
  out = inv.aux
  if(is.matrix(inv.aux)){ out = sqrt(diag(inv.aux)) } # Standard error.
  return(out)
}

#############################################################################################
#' @title Qf_bet
#' @description Auxiliary function required in the Expectation-Maximization algorithm (Q-function related to the beta model).
#' @param theta vector of parameters (all coefficients).
#' @param phiold previous value of the precision parameter (phi).
#' @param z response vector with 0 < z[i] < 1.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @return Scalar representing the output of this auxiliary function for the beta case.
Qf_bet = function(theta,phiold,z,x,v)
{
  n = length(z)
  nkap = ncol(x)
  kap = theta[1:nkap];
  lam = theta[-c(1:nkap)];
  mu = exp(x%*%kap)/(1+exp(x%*%kap)); # mean parameter.
  phi = exp(v%*%lam); # precision parameter.
  #
  mu0phi = mu*phi
  mu1phi = (1-mu)*phi
  mu0phi[which(mu0phi <=0)] = 10^(-10)
  mu1phi[which(mu1phi <=0)] = 10^(-10)
  #
  out = phi*(mu*log(z/(1-z))+digamma(phiold)+log(1-z));
  out = out -lgamma(mu0phi)-lgamma(mu1phi);
  out = out -log(z*(1-z)) -digamma(phiold) -log(1-z) -phiold;
  return(sum(out))
}

#############################################################################################
#' @title gradtheta_bet
#' @description Function to calculate the gradient required for optimization via \code{optim}.
#' This option is related to the beta regression.
#' @param theta vector of parameters (all coefficients).
#' @param phiold previous value of the precision parameter (phi).
#' @param z response vector with 0 < z[i] < 1.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @return Scalar representing the output of this auxiliary gradient function for the beta case.
gradtheta_bet = function(theta,phiold,z,x,v)
{
  n = length(z)
  nkap = ncol(x)
  kap = theta[1:nkap]
  lam = theta[-c(1:nkap)]
  nlam = length(lam)
  mu = exp(x%*%kap)/(1+exp(x%*%kap))
  phi = exp(v%*%lam)
  #
  Ukap = array(0,c(1,nkap))
  aux = log(z)-log(1-z)-digamma(mu*phi)+digamma((1-mu)*phi)
  aux = aux*mu*(1-mu)*phi
  for(j in 1:nkap){ Ukap[j] = sum(aux*x[,j]) }
  #
  Ulam = array(0,c(1,nlam))
  aux = mu*(log(z)-log(1-z)) +digamma(phiold)
  aux = aux +log(1-z)-mu*digamma(mu*phi)
  aux = aux -(1-mu)*digamma((1-mu)*phi)
  aux = aux*phi
  for(j in 1:nlam){ Ulam[j] = sum(aux*v[,j]) }
  return(c(Ukap,Ulam))
}

#############################################################################################
#' @title EMrun_bet
#' @description Function to run the Expectation-Maximization algorithm for the beta regression.
#' @param kap initial values for the coefficients in kappa related to the mean parameter.
#' @param lam initial values for the coefficients in lambda related to the precision parameter.
#' @param z response vector with 0 < z[i] < 1.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param epsilon tolerance to control convergence criterion.
#' @return Vector containing the estimates for kappa and lambda.
EMrun_bet = function(kap,lam,z,x,v,epsilon)
{
  n = nrow(x); nkap = ncol(x); nlam = ncol(v)
  mu = exp(x%*%kap)/(1+exp(x%*%kap))
  phi = exp(v%*%lam)
  theta = c(kap,lam)
  count = 0
  repeat{
    kap_r = kap; lam_r = lam; theta_r = theta
    mu_r = mu; phi_r = phi
    ### E step ------------------------------
    ### M step ------------------------------
    M = stats::optim(par = theta, fn = Qf_bet, gr = gradtheta_bet, phiold=phi_r, z=z ,x=x, v=v, control=list(fnscale=-1), method = 'BFGS');
    theta = M$par; kap = theta[1:nkap]; lam = theta[-c(1:nkap)]
    mu = exp(x%*%kap)/(1+exp(x%*%kap)); phi = exp(v%*%lam)
    # Compute Q -----------------------------
    Q_r = Qf_bet(theta_r,phi_r,z,x,v)
    Q = Qf_bet(theta,phi_r,z,x,v)
    ### Convergence criterion ---------------
    term1 = sqrt(sum((theta-theta_r)^2))
    term2 = abs(Q - Q_r)
    ### -------------------------------------
    count = count+1
    if(max(term1,term2) < epsilon){ break }
    if(count >= 10000){epsilon = 10^(-3)}
    ### -------------------------------------
  }
  # g(phi)
  gphi = 1/(1+phi)
  out = list()
  out[[1]] = c(kap,lam)
  out[[2]] = cbind(mu,gphi)
  out[[3]] = count
  names(out) = c("coeff","mu_gphi","n_iter")
  return(out)
}

#############################################################################################
#' @title simdata_bet
#' @description Function to generate synthetic data from the beta regression.
#' @param kap coefficients kappa related to the mean parameter.
#' @param lam coefficients lambda related to the precision parameter.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @return Response vector z (with 0 < z[i] < 1).
#' @seealso
#' \code{\link{simdata_bes}}, \code{\link{dbessel}}, \code{\link{dbbtest}}
#' @examples
#' n = 100; x = cbind(rbinom(n, 1, 0.5), runif(n, -1, 1)); v = runif(n, -1, 1);
#' z = simdata_bet(kap = c(1, -1, 0.5), lam = c(0.5,- 0.5), x, v)
#' hist(z, xlim = c(0, 1), prob = TRUE)
#' @export
simdata_bet <- function(kap,lam,x,v)
{
  ncolx = ncol(x); if(is.null(ncolx)==TRUE){ncolx=1}
  ncolv = ncol(v); if(is.null(ncolv)==TRUE){ncolv=1}
  nkap = length(kap); nlam = length(lam)
  if((nkap-ncolx)==1){x = cbind(1,x)}
  if((nlam-ncolv)==1){v = cbind(1,v)}
  if(abs(nkap-ncolx)>1){stop("check dimension of kappa and x")}
  if(abs(nlam-ncolv)>1){stop("check dimension of lambda and v")}
  #
  xk = x%*%kap
  mu = exp(xk)/(1+exp(xk))
  phi = exp(v%*%lam)
  s1 = mu*phi
  s2 = phi*(1-mu)
  n = length(s1)
  Z = rep(0,n)
  for(i in 1:n){
   Z[i] = stats::rbeta(1, shape1=s1[i], shape2=s2[i])
   if(Z[i] < 0.00001){Z[i] = 0.00001}
   if(Z[i] > 0.99999){Z[i] = 0.99999}
  }
  return(Z)
}

#############################################################################################
#' @title envelope_bet
#' @description Function to calculate envelopes based on residuals for the beta regression.
#' @param residual character indicating the type of residual ("pearson" or "score").
#' @param kap coefficients in kappa related to the mean parameter.
#' @param lam coefficients in lambda related to the precision parameter.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param nsim_env number of synthetic data sets to be generated.
#' @param prob confidence level of the envelope (a number between 0 and 1).
#' @param n sample size.
#' @param epsilon tolerance parameter used in the Expectation-Maximization algorithm applied to the synthetic data.
#' @return Matrix with dimension 2 x n (1st row = upper bound, second row = lower bound).
envelope_bet <- function(residual,kap,lam,x,v,nsim_env,prob,n,epsilon)
{
  Res = array(0,c(nsim_env,n))
  bar = utils::txtProgressBar(min=0,max=nsim_env,style=3)
  for(i in 1:nsim_env){
    zsim = simdata_bet(kap,lam,x,v)
    EM = EMrun_bet(kap,lam,z=zsim,x,v,epsilon)
    if(residual=="pearson"){
      musim = as.numeric(EM[[2]][,1])
      gpsim = as.numeric(EM[[2]][,2])
      Res[i,] = (zsim-musim)/sqrt(gpsim*musim*(1-musim)) }
    if(residual=="score"){
      nkap = length(kap)
      kapsim = EM[[1]][1:nkap]
      lamsim = EM[[1]][-(1:nkap)]
      Res[i,] = score_residual_bet(kapsim,lamsim,zsim,x,v) }
    if(residual=="quantile"){
      nkap = length(kap)
      kapsim = EM[[1]][1:nkap]
      lamsim = EM[[1]][-(1:nkap)]
      Res[i,] = quantile_residual_bet(kapsim,lamsim,zsim,x,v) }
    utils::setTxtProgressBar(bar,i)
  }
  Res = t(apply(Res,1,sort))
  Res = apply(Res,2,sort)
  id1 = max(1,round(nsim_env*(1-prob)/2))
  id2 = round(nsim_env*(1+prob)/2)
  Env = rbind(Res[id2,],apply(Res,2,mean),Res[id1,])
  rownames(Env) = c("upper","mean","lower")
  return(Env)
}

#############################################################################################
#' @title pred_accuracy_bet
#' @description Function to calculate the Residual Sum of Suqares for partitions (training and test sets) of
#' the data set. Residuals are calculated based on the beta regression.
#' @param residual character indicating the type of residual ("pearson", "score" or "quantile").
#' @param kap coefficients in kappa related to the mean parameter.
#' @param lam coefficients in lambda related to the precision parameter.
#' @param z response vector with 0 < z[i] < 1.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param ntest number observations in the test set for prediction.
#' @param predict number of partitions (training and test sets) to be evaluated.
#' @param epsilon tolerance parameter used in the Expectation-Maximization algorithm for the training data set.
#' @return Vector containing the RSS for each partition of the full data set.
pred_accuracy_bet = function(residual,kap,lam,z,x,v,ntest,predict,epsilon){
  n = length(z)
  nkap = ncol(x)
  nlam = ncol(v)
  RSS_pred = rep(0,predict)
  if(predict>=50){ bar = utils::txtProgressBar(min=0,max=predict,style=3) }
  for(i in 1:predict){
    id_pred = sample(1:n,ntest,replace=FALSE)
    ztr = z[-id_pred]
    xtr = as.matrix(x[-id_pred,])
    vtr = as.matrix(v[-id_pred,])
    zte = z[id_pred]
    xte = as.matrix(x[id_pred,])
    vte = as.matrix(v[id_pred,])
    EMtr = EMrun_bet(kap,lam,ztr,xtr,vtr,epsilon)
    kaptr = EMtr[[1]][1:nkap]
    lamtr = EMtr[[1]][-(1:nkap)]
    mupred = exp(xte%*%kaptr)/(1+exp(xte%*%kaptr))
    phipred = exp(vte%*%lamtr)
    gphi_pred = 1/(1+phipred)
    if(residual=="pearson"){ respred = (zte-mupred)/sqrt(gphi_pred*mupred*(1-mupred)) }
    if(residual=="score"){ respred = score_residual_bet(kaptr,lamtr,zte,xte,vte) }
    if(residual=="quantile"){ respred = quantile_residual_bet(kaptr,lamtr,zte,xte,vte) }
    RSS_pred[i] = sum(respred^2)
    if(predict>=50){ utils::setTxtProgressBar(bar,i) }
  }
  return(RSS_pred)
}

#############################################################################################
#' @title score_residual_bet
#' @description Function to calculate empirical score residuals based on the bessel regression.
#' @param kap coefficients in kappa related to the mean parameter.
#' @param lam coefficients in lambda related to the precision parameter.
#' @param z response vector with 0 < z[i] < 1.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param nsim_score number synthetic data sets (default = 100) to be generated as a support to estime the mean and s.d. of log(z)-log(1-z).
#' @return Vector containing the score residuals.
score_residual_bet = function(kap,lam,z,x,v,nsim_score=100){
  n = length(z)
  zs = array(0,c(nsim_score,n))
  mu.sim = exp(x%*%kap)/(1+exp(x%*%kap))
  phi.sim = exp(v%*%lam)
  for(j in 1:nsim_score){
    aux = numeric(n)
    for(i in 1:n){
      s1 = mu.sim[i]*phi.sim[i]
      s2 = (1-mu.sim[i])*phi.sim[i]
      aux[i] = stats::rbeta(1, shape1=s1,shape2=s2)
      if(aux[i] < 0.00001){aux[i] = 0.00001}
      if(aux[i] > 0.99999){aux[i] = 0.99999} }
    zs[j,] = log(aux)-log(1-aux)
  }
  me_zs = apply(zs,2,mean)
  sd_zs = apply(zs,2,stats::sd)
  out = log(z)-log(1-z)
  out = (out-me_zs)/sd_zs
  return(out)
}

#############################################################################################
#' @title quantile_residual_bet
#' @description Function to calculate quantile residuals based on the beta regression. Details about this type of residual can be found in \emph{Pereira (2019)}.
#' @param kap coefficients in kappa related to the mean parameter.
#' @param lam coefficients in lambda related to the precision parameter.
#' @param z response vector with 0 < z[i] < 1.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @references DOI:10.1080/03610918.2017.1381740 (\href{https://www.tandfonline.com/doi/abs/10.1080/03610918.2017.1381740}{Pereira; 2019})
#' @return Vector containing the quantile residuals.
quantile_residual_bet = function(kap,lam,z,x,v){
  n = length(z)
  xk = x%*%kap
  mu = exp(xk)/(1+exp(xk))
  phi = exp(v%*%lam)
  out = rep(0,n)
  for(i in 1:n){
    s1 = mu[i]*phi[i]
    s2 = (1-mu[i])*phi[i]
    aux = stats::pbeta(z[i], shape1=s1, shape2=s2)
    out[i] = stats::qnorm(aux,mean=0,sd=1)
  }
  return(out)
}
