#' ANOVA in kernel distribution estimation with binned data.
#'
#' @param n Vector of positive integers. Sizes of the complete samples corresponding to each treatment.
#' @param y Vector. Observed values. They define the extremes of the sequence of intervals in which data is binned.
#' @param trt.w Matrix. Proportion of observations within each interval. Each column corresponds to a different treatment.
#' @param abs.values Logical. Indicates if the values of trt.w are given in absolute (TRUE) or relative (FALSE) format.
#' @param alpha Real number between 0 and 1. Significance level of the test.
#' @param B Positive integer. Number of bootstrap replicates used to compute the confidence bands.
#' 
#' @details 
#' Constructs bootstrap confidence bands for each treatment and checks whether they overlap or not.
#' 
#' @return TRUE if the null hypothesis is accepted and FALSE otherwise.
#' 
#' @references 
#' \insertRef{TesisMiguel2015}{binnednp}
#'
#' @useDynLib binnednp
#' @importFrom Rcpp sourceCpp
#'
#' @export
anv.binned <- function(n,y,trt.w,abs.values=FALSE,alpha=0.05,B=500)
{
  
  ntrt <- ncol(trt.w)
  alpha2 <- alpha/ntrt
  if(abs.values)
  {
    for(i in 1:ntrt)
    {
      trt.w[,i] <- trt.w[,i]/n[i]
    }
  }
  
  ncores <- parallel::detectCores()
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  parout <- foreach::foreach(i=1:ntrt,.combine=rbind,.packages=c("Rcpp","binnednp")) %dopar%{
    band.boot <- bw.dist.binned.boot(n[i],y,trt.w[,i],confband=TRUE,alpha=alpha2,plot=FALSE,print=FALSE,parallel=FALSE,B=B)$confband
    return(band.boot)
  }
  
  parallel::stopCluster(cl)
  H0.boot <- TRUE
  maxmin.i <- numeric(ntrt)
  minmax.i <- numeric(ntrt)
  kk <- length(y)
  for(i in 1:kk)
  {
    for(j in 1:ntrt)
    {
      maxmin.i[j] <- max(parout[i+kk*(j-1),2])
      minmax.i[j] <- min(parout[i+kk*(j-1),1])
    }
    H0.boot <- ifelse(any(maxmin.i>minmax.i),FALSE,TRUE)
    if(H0.boot == FALSE) break
  }
  
  return(H0.boot)
  
  
}
