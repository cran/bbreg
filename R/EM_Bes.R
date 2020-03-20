#############################################################################################
#' @title infmat_bes
#' @description Function to compute standard errors based on the Fisher information matrix for the bessel regression.
#' @param theta vector of parameters (all coefficients: kappa and lambda).
#' @param z response vector with 0 < z_i < 1.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @return Vector of standard errors.
infmat_bes = function(theta,z,x,v)
{
  n = length(z)
  nkap = ncol(x)
  kap = theta[1:nkap]
  lam = theta[-c(1:nkap)]
  nlam = length(lam)
  mu = exp(x%*%kap)/(1+exp(x%*%kap))
  phi = exp(v%*%lam)
  wz = Ew1z(z,mu,phi)
  chi = Ew2z(z,mu,phi)
  xi2 = ((mu^2)/z) + (((1-mu)^2)/(1-z))
  xit = (mu-z)/(z*(1-z))
  # E[-d2/d.kap.kap]
  d2.k.k = array(0,c(nkap,nkap))
  aux1 = (1-2*mu)*xit + mu*(1-mu)/(z*(1-z))
  aux1 = aux1*wz*(phi^2) + 2
  aux1 = aux1*mu*(1-mu)
  for(j in 1:nkap){
    for(l in 1:nkap){
      d2.k.k[j,l] = sum(aux1*x[,j]*x[,l])
    }
  }
  # E[-d2/d.lam.lam]
  d2.l.l = array(0,c(nlam,nlam))
  aux1 = (2*wz*phi*xi2-1)*phi
  for(j in 1:nlam){
    for(l in 1:nlam){
      d2.l.l[j,l] = sum(aux1*v[,j]*v[,l])
    }
  }
  # E[-d2/d.kap.lam]
  d2.k.l = array(0,c(nkap,nlam))
  aux1 = 2*(phi^2)*wz*mu*(1-mu)*xit
  for(j in 1:nkap){
    for(l in 1:nlam){
      d2.k.l[j,l] = sum(aux1*x[,j]*v[,l])
    }
  }
  #
  d2.Q1 = cbind(d2.k.k, d2.k.l)
  d2.Q2 = cbind(t(d2.k.l), d2.l.l)
  d2.Q = rbind(d2.Q1, d2.Q2)
  #
  # E[d/d.kap.j d/d.kap.l]
  dk.dk = array(0,c(nkap,nkap))
  aux1 = ((1-2*mu)^2)-2*(1-2*mu)*wz*(phi^2)*mu*(1-mu)*xit
  aux2 = (mu^2)*((1-mu)^2)*(phi^4)*chi*(xit^2)
  aux = 1-2*mu-mu*(1-mu)*(phi^2)*wz*xit
  for(j in 1:nkap){
    for(l in 1:nkap){
      dk.dk[j,l] = sum((aux1+aux2)*x[,j]*x[,l])
      for(i in 1:n){
        for(k in 1:n){
          if(k != i){ dk.dk[j,l] = dk.dk[j,l] + aux[i]*aux[k]*x[i,j]*x[k,l] }
        }
      }
    }
  }
  # E[d/d.lam.j d/d.lam.l]
  dl.dl = array(0,c(nlam,nlam))
  aux1 = ((2+phi)^2)-2*(2+phi)*(phi^2)*wz*xi2
  aux2 = (phi^4)*chi*(xi2^2)
  aux = 2+phi-(phi^2)*wz*xi2
  for(j in 1:nlam){
    for(l in 1:nlam){
      dl.dl[j,l] = sum((aux1+aux2)*v[,j]*v[,l])
      for(i in 1:n){
        for(k in 1:n){
          if(k != i){ dl.dl[j,l] = dl.dl[j,l] + aux[i]*aux[k]*v[i,j]*v[k,l] }
        }
      }
    }
  }
  # E[d/d.kap.j d/d.lam.l]
  dk.dl = array(0,c(nkap,nlam))
  aux1 = (1-2*mu)*(2+phi)
  aux2 = -(1-2*mu)*wz*(phi^2)*xi2
  aux3 = -(2+phi)*wz*(phi^2)*mu*(1-mu)*xit
  aux4 = chi*(phi^4)*mu*(1-mu)*xit*xi2
  aux5 = 1-2*mu-wz*(phi^2)*mu*(1-mu)*xit
  aux6 = 2+phi-wz*(phi^2)*xi2
  for(j in 1:nkap){
    for(l in 1:nlam){
      dk.dl[j,l] = sum((aux1+aux2+aux3+aux4)*x[,j]*v[,l])
      for(i in 1:n){
        for(k in 1:n){
          if(k != i){ dk.dl[j,l] = dk.dl[j,l] + aux5[i]*aux6[k]*x[i,j]*v[k,l] }
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
#' @title Ew1z
#' @description Auxiliary function required in the Expectation-Maximization algorithm (E-step) and in the calculation of the
#' Fisher information matrix. It represents the conditional expected value E(W[i]^s|Z[i]), with s = -1; i.e.,
#' latent W[i]^(-1) given the observed Z[i].
#' @param z response vector with 0 < z[i] < 1.
#' @param mu mean parameter (vector having the same size of z).
#' @param phi precision parameter (vector having the same size of z).
#' @return Vector of expected values.
Ew1z = function(z,mu,phi)
{
  xi = sqrt( ((mu^2)/z) + (((1-mu)^2)/(1-z)) )
  prod = phi*xi
  bK.num = besselK(prod,-2,expon.scaled=TRUE)
  bK.den = besselK(prod,-1,expon.scaled=TRUE)
  out = (1/prod)*(bK.num/bK.den)
  return(out)
}

#############################################################################################
#' @title Ew2z
#' @description Auxiliary function required in the calculation of the
#' Fisher information matrix. It represents the conditional expected value E(W[i]^s|Z[i]), with s = -2; i.e.,
#' latent W[i]^(-2) given the observed Z[i].
#' @param z response vector with 0 < z_i < 1.
#' @param mu mean parameter (vector having the same size of z).
#' @param phi precision parameter (vector having the same size of z).
#' @return vector of expected values.
Ew2z = function(z,mu,phi)
{
  xi = sqrt( ((mu^2)/z) + (((1-mu)^2)/(1-z)) )
  prod = phi*xi
  prod2 = prod^2
  bK.num = besselK(prod,-3,expon.scaled=TRUE)
  bK.den = besselK(prod,-1,expon.scaled=TRUE)
  out = (1/prod2)*(bK.num/bK.den)
  return(out)
}

#############################################################################################
#' @title Qf_bes
#' @description Auxiliary function required in the Expectation-Maximization algorithm (Q-function related to the bessel model).
#' @param theta vector of parameters (all coefficients: kappa and lambda).
#' @param wz parameter representing E(1/W[i]|Z[i] = z[i], theta).
#' @param z response vector with 0 < z[i] < 1.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @return Scalar representing the output of this auxiliary function for the bessel case.
Qf_bes = function(theta,wz,z,x,v)
{
  n = length(z)
  nkap = ncol(x)
  kap = theta[1:nkap]
  lam = theta[-c(1:nkap)]
  mu = exp(x%*%kap)/(1+exp(x%*%kap)) # mean parameter.
  phi = exp(v%*%lam) # phi precision parameter.
  #
  out1 = log(mu) + log(1-mu) + 2*log(phi) + phi
  out2 = 0.5*wz*(phi^2)*( ((mu^2)/z) + (((1-mu)^2)/(1-z)) )
  return( sum(out1-out2) )
}

#############################################################################################
#' @title gradtheta_bes
#' @description Function to calculate the gradient required for optimization via \code{optim}.
#' This option is related to the bessel regression.
#' @param theta vector of parameters (all coefficients: kappa and lambda).
#' @param wz parameter representing E(1/W[i]|Z[i] = z[i], theta).
#' @param z response vector with 0 < z[i] < 1.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @return Scalar representing the output of this auxiliary gradient function for the bessel case.
gradtheta_bes = function(theta,wz,z,x,v)
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
  aux = (1/mu)-(1/(1-mu))
  aux = aux -wz*(phi^2)*( (mu/z)-((1-mu)/(1-z)) )
  aux = aux*mu*(1-mu)
  for(j in 1:nkap){ Ukap[j] = sum(aux*x[,j]) }
  #
  Ulam = array(0,c(1,nlam))
  aux = (2/phi)+1-wz*phi*( ((mu^2)/z) + (((1-mu)^2)/(1-z)) )
  aux = aux*phi
  for(j in 1:nlam){ Ulam[j] = sum(aux*v[,j]) }
  return(c(Ukap,Ulam))
}

#############################################################################################
#' @title EMrun_bes
#' @description Function to run the Expectation-Maximization algorithm for the bessel regression.
#' @param kap initial values for the coefficients in kappa related to the mean parameter.
#' @param lam initial values for the coefficients in lambda related to the precision parameter.
#' @param z response vector with 0 < z_i < 1.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param epsilon tolerance to control the convergence criterion.
#' @return Vector containing the estimates for kappa and lambda in the bessel regression.
EMrun_bes = function(kap,lam,z,x,v,epsilon){
  n = length(z)
  nkap = ncol(x)
  nlam = ncol(v)
  mu = exp(x%*%kap)/(1+exp(x%*%kap))
  phi = exp(v%*%lam)
  theta = c(kap,lam)
  count = 0
  repeat{
    kap_r = kap; lam_r = lam; theta_r = theta
    mu_r = mu; phi_r = phi
    ### E step ------------------------------
    wz_r = Ew1z(z,mu_r,phi_r)
    ### M step ------------------------------
    M = stats::optim(par = theta, fn = Qf_bes, gr = gradtheta_bes, wz=wz_r, z=z ,x=x, v=v, control=list(fnscale=-1), method = 'BFGS');
    theta = M$par; kap = theta[1:nkap]; lam = theta[-c(1:nkap)]
    mu = exp(x%*%kap)/(1+exp(x%*%kap)); phi = exp(v%*%lam)
    # Compute Q -----------------------------
    Q_r = Qf_bes(theta_r,wz_r,z,x,v)
    Q = Qf_bes(theta,wz_r,z,x,v)
    ### Convergence criterion ---------------
    term1 = sqrt(sum((theta-theta_r)^2))
    term2 = abs(Q - Q_r)
    ### -------------------------------------
    count = count+1
    if(max(term1,term2) < epsilon){ break }
    if(count >= 10000){epsilon = 10^(-3)}
    ### -------------------------------------
  }
  gphi = (1-phi+(phi^2)*exp(phi)*expint_En(phi,order=1))/2
  out = list()
  out[[1]] = c(kap,lam)
  out[[2]] = cbind(mu,gphi)
  out[[3]] = count
  names(out) = c("coeff","mu_gphi","n_iter")
  return(out)
}

#############################################################################################
#' @title simdata_bes
#' @description Function to generate synthetic data from the bessel regression.
#' Requires the R package "statmod" to deal with the Inverse Gaussian distribution (\emph{Giner and Smyth, 2016}).
#' @param kap coefficients in kappa related to the mean parameter.
#' @param lam coefficients in lambda related to the precision parameter.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @return Response vector z (with 0 < z[i] < 1).
#' @references DOI:10.32614/RJ-2016-024 (\href{https://journal.r-project.org/archive/2016/RJ-2016-024/index.html}{Giner and Smyth; 2016})
#' @seealso
#' \code{\link{dbessel}}, \code{\link{dbbtest}}, \code{\link{simdata_bet}}
#' @examples
#' n = 100; x = cbind(rbinom(n, 1, 0.5), runif(n, -1, 1)); v = runif(n, -1, 1);
#' z = simdata_bes(kap = c(1, -1, 0.5), lam = c(0.5, -0.5), x, v)
#' hist(z, xlim = c(0, 1), prob = TRUE)
#' @export
simdata_bes <- function(kap,lam,x,v)
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
  a = mu*phi
  b = phi*(1-mu)
  n = length(a)
  Y1 = rep(0,n)
  Y2 = rep(0,n)
  for(i in 1:n){
    # Y1 ~ Inv.Gaussian s.t. E(Y1) = a[i], Var(Y1) = a[i]
    Y1[i] = rinvgauss(1, mean = a[i], shape = a[i]^2)
    # Y2 ~ Inv.Gaussian s.t. E(Y2) = b[i], Var(Y2) = b[i]
    Y2[i] = rinvgauss(1, mean = b[i], shape = b[i]^2)
  }
  Z = Y1/(Y1+Y2) # Z[i] ~ Bessel(mu[i],phi[i])
  return(Z)
}

#############################################################################################
#' @title envelope_bes
#' @description Function to calculate envelopes based on residuals for the bessel regression.
#' @param residual character indicating the type of residual ("pearson", "score" or "quantile").
#' @param kap coefficients in kappa related to the mean parameter.
#' @param lam coefficients in lambda related to the precision parameter.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param nsim_env number of synthetic data sets to be generated.
#' @param n sample size.
#' @param prob confidence level of the envelope (number between 0 and 1).
#' @param epsilon tolerance parameter used in the Expectation-Maximization algorithm applied to the synthetic data.
#' @return Matrix with dimension 2 x n (1st row = upper bound, second row = lower bound).
envelope_bes <- function(residual,kap,lam,x,v,nsim_env,prob,n,epsilon)
{
 Res = array(0,c(nsim_env,n))
 bar = utils::txtProgressBar(min=0,max=nsim_env,style=3)
 for(i in 1:nsim_env){
   zsim = simdata_bes(kap,lam,x,v)
   EM = EMrun_bes(kap,lam,z=zsim,x,v,epsilon)
   if(residual=="pearson"){
     musim = as.numeric(EM[[2]][,1])
     gpsim = as.numeric(EM[[2]][,2])
     Res[i,] = (zsim-musim)/sqrt(gpsim*musim*(1-musim)) }
   if(residual=="score"){
     nkap = length(kap)
     kapsim = EM[[1]][1:nkap]
     lamsim = EM[[1]][-(1:nkap)]
     Res[i,] = score_residual_bes(kapsim,lamsim,zsim,x,v) }
   if(residual=="quantile"){
     nkap = length(kap)
     kapsim = EM[[1]][1:nkap]
     lamsim = EM[[1]][-(1:nkap)]
     Res[i,] = quantile_residual_bes(kapsim,lamsim,zsim,x,v) }
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
#' @title pred_accuracy_bes
#' @description Function to calculate the Residual Sum of Suqares for partitions (training and test sets) of
#' the data set. Residuals are calculated here based on the bessel regression.
#' @param residual Character indicating the type of residual ("pearson", "score" or "quantile").
#' @param kap coefficients in kappa related to the mean parameter.
#' @param lam coefficients in lambda related to the precision parameter.
#' @param z response vector with 0 < z[i] < 1.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param ntest number of observations in the test set for prediction.
#' @param predict number of partitions (training and test sets) to be evaluated.
#' @param epsilon tolerance parameter used in the Expectation-Maximization algorithm for the training data set.
#' @return Vector containing the RSS for each partition of the full data set.
pred_accuracy_bes = function(residual,kap,lam,z,x,v,ntest,predict,epsilon){
  n = length(z)
  nkap = ncol(x)
  nlam = ncol(v)
  RSS_pred = rep(0,predict)
  if(predict >= 50){ bar = utils::txtProgressBar(min=0,max=predict,style=3) }
  for(i in 1:predict){
    id_pred = sample(1:n,ntest,replace=FALSE)
    ztr = z[-id_pred]
    xtr = as.matrix(x[-id_pred,])
    vtr = as.matrix(v[-id_pred,])
    zte = z[id_pred]
    xte = as.matrix(x[id_pred,])
    vte = as.matrix(v[id_pred,])
    EMtr = EMrun_bes(kap,lam,ztr,xtr,vtr,epsilon)
    kaptr = EMtr[[1]][1:nkap]
    lamtr = EMtr[[1]][-(1:nkap)]
    mupred = exp(xte%*%kaptr)/(1+exp(xte%*%kaptr))
    phipred = exp(vte%*%lamtr)
    gphi_pred = (1-phipred+(phipred^2)*exp(phipred)*expint_En(phipred,order=1))/2
    if(residual=="pearson"){ respred = (zte-mupred)/sqrt(gphi_pred*mupred*(1-mupred)) }
    if(residual=="score"){ respred = score_residual_bes(kaptr,lamtr,zte,xte,vte) }
    if(residual=="quantile"){ respred = quantile_residual_bes(kaptr,lamtr,zte,xte,vte) }
    RSS_pred[i] = sum(respred^2)
    if(predict >= 50){ utils::setTxtProgressBar(bar,i) }
  }
  return(RSS_pred)
}

#############################################################################################
#' @title score_residual_bes
#' @description Function to calculate the empirical score residuals based on the bessel regression.
#' @param kap coefficients in kappa related to the mean parameter.
#' @param lam coefficients in lambda related to the precision parameter.
#' @param z response vector with 0 < z[i] < 1.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @param nsim_score number synthetic data sets (default = 100) to be generated as a support to estime mean and s.d. of log(z)-log(1-z).
#' @return Vector containing the score residuals.
score_residual_bes = function(kap,lam,z,x,v,nsim_score=100){
  n = length(z)
  zs = array(0,c(nsim_score,n))
  for(j in 1:nsim_score){
    aux = simdata_bes(kap,lam,x,v)
    zs[j,] = log(aux)-log(1-aux) }
  me_zs = apply(zs,2,mean)
  sd_zs = apply(zs,2,stats::sd)
  out = log(z)-log(1-z)
  out = (out-me_zs)/sd_zs
  return(out)
}

#############################################################################################
#' @title dbessel
#' @description Function to calculate the probability density of the bessel distribution.
#' @param z scalar (0 < z < 1) for which the p.d.f. is to be evaluated.
#' @param mu scalar representing the mean parameter.
#' @param phi scalar representing the precision parameter.
#' @return scalar expressing the value of the density at z.
#' @seealso
#' \code{\link{simdata_bes}}, \code{\link{dbbtest}}, \code{\link{simdata_bet}}
#' @examples
#' z = seq(0.01, 0.99, 0.01); np = length(z);
#' density = rep(0, np)
#' for(i in 1:np){ density[i] = dbessel(z[i], 0.5, 1) }
#' plot(z, density, type = "l", lwd = 2, cex.lab = 2, cex.axis = 2)
#' @export
dbessel = function(z,mu,phi){
  zeta = sqrt( (z*(1-2*mu)+mu^2)/(z-z^2) )
  out = mu*(1-mu)*phi*exp(phi)*besselK((phi*zeta),1)
  out = out/(zeta*pi*(z*(1-z))^(3/2))
  return(out)
}

#############################################################################################
#' @title quantile_residual_bes
#' @description Function to calculate quantile residuals based on the bessel regression. Details about this type of residual can be found in \emph{Pereira (2019)}.
#' @param kap coefficients in kappa related to the mean parameter.
#' @param lam coefficients in lambda related to the precision parameter.
#' @param z response vector with 0 < z[i] < 1.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
#' @references DOI:10.1080/03610918.2017.1381740 (\href{https://www.tandfonline.com/doi/abs/10.1080/03610918.2017.1381740}{Pereira; 2019})
#' @return Vector containing the quantile residuals.
quantile_residual_bes = function(kap,lam,z,x,v){
  n = length(z)
  xk = x%*%kap
  mu = exp(xk)/(1+exp(xk))
  phi = exp(v%*%lam)
  out = rep(0,n)
  for(i in 1:n){
   aux = stats::integrate(f = dbessel, lower = 0, upper = z[i], mu[i], phi[i])
   aux = aux$value
   if(aux > 0.99999){ aux = 0.99999 }
   if(aux < 0.00001){ aux = 0.00001 }
   out[i] = stats::qnorm(aux,mean=0,sd=1)
  }
  return(out)
}
