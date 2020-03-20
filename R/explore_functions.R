## Import key packages
#' @import Formula
#' @import betareg
#' @import expint
#' @import statmod

#############################################################################################
#' @title eplot
#' @description Function to build the envelope graph of residuals against the theoretical quantiles from
#' the standard normal distribution. Envelopes will not be displayed when the data is fitted with argument
#' envelope = 0.
#' @param M object of class "bbreg" containing results from the fitted model.
#' If the model is fitted with envelope = 0, the graph cannot be created.
#' @param pch symbol to use. This can either be a single character or an integer code for one of a set of graphics symbols (default: solid circle).
#' @param xlim x limits (x1, x2) of the plot. Note that x1 > x2 is allowed and leads to a 'reversed axis'.
#' @param ylim y limits of the plot.
#' @param xlab label for the x axis.
#' @param ylab label for the y axis.
#' @param cex.lab magnification to be used for x and y labels relative to the current setting of cex.
#' @param cex.axis magnification to be used for axis annotation relative to the current setting of cex.
#' @seealso
#' \code{\link{bbsummary}}, \code{\link{dbbtest}}
#' @examples
#' \donttest{
#' n = 100; x = cbind(rbinom(n, 1, 0.5), runif(n, -1, 1)); v = runif(n, -1, 1);
#' z = simdata_bes(kap = c(1, 1, -0.5), lam = c(0.5, -0.5), x, v)
#' M = bbreg(z ~ x | v, envelope = 10)
#' eplot(M)}
#' @export
eplot = function(M,pch=19,xlim=NULL,ylim=NULL,xlab=NULL,ylab=NULL,cex.lab=NULL,cex.axis=NULL)
{
  res = M$residuals
  env = M$envelope
  name = M$modelname
  n = length(res)
  residualname = paste0(toupper(substring(M$residualname,1,1)),substring(M$residualname,2))
  if(is.null(ylim)==TRUE){ yli = range(res,env) }else{ yli = ylim }
  if(is.null(xlim)==TRUE){ xli = c(stats::qnorm(0.5/n),stats::qnorm(1-0.5/n)) }else{ xli = xlim }
  if(is.null(ylab)==TRUE){ yla = paste0(residualname," residual") }else{ yla = ylab }
  if(is.null(xlab)==TRUE){ xla = "Theoretical quantiles" }else{ xla = xlab }
  if(is.null(cex.lab)==TRUE){ cl = 1.5 }else{ cl = cex.lab }
  if(is.null(cex.axis)==TRUE){ ca = 1.5 }else{ ca = cex.axis }
  RR = stats::qqnorm(res,xlab=xla,ylab=yla,xlim=xli,ylim=yli,pch="",main=name,cex.lab=cl,cex.axis=ca)
  if(is.null(env)==FALSE){
    aux = sort(RR$x)
    graphics::lines(aux,env[1,],col=grDevices::rgb(0.7,0.7,0.7))
    graphics::lines(aux,env[3,],col=grDevices::rgb(0.7,0.7,0.7))
    graphics::polygon(c(aux,rev(aux)),c(env[3,],rev(env[1,])),col=grDevices::rgb(0.7,0.7,0.7),border=NA)
    graphics::lines(aux,env[2,],lty=2,lwd=2) }
  graphics::points(RR$x,RR$y,pch=pch)
}

#############################################################################################
#' @title bbsummary
#' @description Function providing a summary of results related to the regression model (bessel or beta).
#' @param M object of class "bbreg" containing results from the fitted model.
#' @seealso
#' \code{\link{eplot}}, \code{\link{dbbtest}}
#' @examples
#' \donttest{
#' fit1 = bbreg(agreement ~ priming + eliciting, data = WT)
#' bbsummary(fit1)}
#' @export
bbsummary = function(M)
{
  nkap = length(M$kappa)
  nlam = length(M$lambda)
  #
  itc = M$intercept
  varnames = NULL
  if(itc[1]==TRUE){ varnames = "intercept" }
  varnames = c(varnames,labels(stats::terms(stats::formula(M$call,rhs=1))))
  if(itc[2]==TRUE){ varnames = c(varnames,"intercept") }
  varnames = c(varnames,labels(stats::terms(stats::formula(M$call,rhs=2))))
  #
  if(length(varnames) < (nkap+nlam)){ varnames = names(M$start) }
  #
  Est = c(M$kappa,M$lambda)
  SEr = M$std_errors
  tab = cbind(Est, SEr, Est/SEr, 2*stats::pnorm(-abs(Est/SEr)))
  colnames(tab) = c("Estimate", "Std.error", "z-value", "Pr(>|z|)")
  rownames(tab) = varnames
  tab = list(mean = tab[seq.int(length.out = nkap), , drop = FALSE], precision = tab[seq.int(length.out = nlam) + nkap, , drop = FALSE])
  #
  digits = max(3,getOption("digits")-3)
  #
  message("\nCall:", deparse(M$call, width.cutoff = floor(getOption("width") * 0.85)),"",sep="\n")
  message(sprintf("%s:\n",M$message))
  message(sprintf("%s\n",paste0("Number of iterations of the EM algorithm = ",M$niter)))
  #
  if(is.null(M$DBB)==FALSE){
    cat(sprintf("\n %s:\n","Results of the discrimination test DBB"))
    print(structure(round(M$DBB,digits=digits))) }
  #
  residualname = paste0(toupper(substring(M$residualname,1,1)),substring(M$residualname,2))
  cat(sprintf("\n %s:\n",paste0(residualname," residuals")))
  print(structure(round(c(M$RSS,as.vector(stats::quantile(M$residuals))), digits=digits), .Names = c("RSS","Min", "1Q", "Median", "3Q", "Max")))
  #
  if(NROW(tab$mean)){
    cat(paste0("\n Coefficients modeling the mean:\n"))
    stats::printCoefmat(tab$mean, digits=digits, signif.legend = FALSE)
  }else{ message("\n No coefficients modeling the mean. \n") }
  #
  if(NROW(tab$precision)) {
    cat(paste0("\n Coefficients modeling the precision:\n"))
    stats::printCoefmat(tab$precision, digits=digits, signif.legend = FALSE)
  }else{ message("\n No coefficients modeling the precision. \n") }
  #
  gp = unique(M$gphi)
  if(length(gp)==1){ cat(sprintf("%s\n",paste0("g(phi) = ",round(gp,digits=digits)))) }
  #
  if(getOption("show.signif.stars")){
    cat("---\n Signif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n\n") }
  #
  if(is.null(M$RSS_pred)==FALSE){
    message(sprintf("%s\n",paste0("Average RSS from the predictive accuracy simulation = ",round(mean(M$RSS_pred),digits=digits))))}
  #
  if(is.null(M$envelope)==FALSE){
    message(sprintf("%s\n",paste0("Percentage of residual within the envelope = ",round(M$envelope_prop,digits=digits))))}
}

#############################################################################################
#' @title startvalues
#' @description Function providing initial values for the Expectation-Maximization algorithm.
#' @param z response vector with 0 < z[i] < 1.
#' @param x matrix containing the covariates for the mean submodel. Each column is a different covariate.
#' @param v matrix containing the covariates for the precision submodel. Each column is a different covariate.
startvalues = function(z,x,v)
{
  n = length(z)
  nkap = ncol(x)
  nlam = ncol(v)
  idx = which(apply(x==1,2,sum)!=n) # find non-intercept columns in x.
  idv = which(apply(v==1,2,sum)!=n) # find non-intercept columns in v.
  #
  if(length(idx) > 0 & length(idv) > 0){
    if(length(idx) == nkap & length(idv) == nlam){
      aux = betareg(z ~ 0 + x[,idx] | 0 + v[,idv], link ="logit",link.phi = "log", type = "ML")
      intercept = c(FALSE,FALSE) }
    if(length(idx) < nkap & length(idv) == nlam){
      aux = betareg(z ~ x[,idx] | 0 + v[,idv], link ="logit",link.phi = "log", type = "ML")
      intercept = c(TRUE,FALSE) }
    if(length(idx) == nkap & length(idv) < nlam){
      aux = betareg(z ~ 0 + x[,idx] | v[,idv], link ="logit",link.phi = "log", type = "ML")
      intercept = c(FALSE,TRUE) }
    if(length(idx) < nkap & length(idv) < nlam){
      aux = betareg(z ~ x[,idx] | v[,idv], link ="logit",link.phi = "log", type = "ML")
      intercept = c(TRUE,TRUE) } }
  #
  if(length(idx) > 0 & length(idv) == 0){
    if(length(idx) == nkap){
      aux = betareg(z ~ 0 + x[,idx], link ="logit",link.phi = "log", type = "ML")
      intercept = c(FALSE,TRUE) }
    if(length(idx) < nkap){
      aux = betareg(z ~ x[,idx], link ="logit",link.phi = "log", type = "ML")
      intercept = c(TRUE,TRUE) } }
  #
  if(length(idx) == 0 & length(idv) > 0){
    if(length(idv)==nlam){
      aux = betareg(z ~ 1 | 0 + v[,idv], link ="logit",link.phi = "log", type = "ML")
      intercept = c(TRUE,FALSE) }
    if(length(idv) <nlam){
      aux = betareg(z ~ 1 | v[,idv], link ="logit", link.phi = "log", type = "ML")
      intercept = c(TRUE,TRUE) } }
  #
  if(length(idx) == 0 & length(idv) == 0){
    aux = betareg(z ~ 1, link ="logit", link.phi = "log", type = "ML")
    intercept = c(TRUE,TRUE) }
  #
  kap_start = as.numeric(aux$start[1:nkap])
  lam_start = as.numeric(aux$start[-c(1:nkap)])
  lam_start[1] = log(2) + log(1+exp(lam_start[1]))
  out = list()
  out[[1]] = kap_start
  out[[2]] = lam_start
  out[[3]] = intercept
  return(out)
}
