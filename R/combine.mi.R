# This function applies Rubin's rule to combine results from different imputations
combine.mi <- function(mi.param.est,mi.param.cov){
  # Input variables:
  # mi.param: a matrix containing estimated parameters in each imputation
  # as its rows
  # mi.param.cov: a list contaning the covariance matrix estimated within 
  # each imputation.
  M=length(mi.param.cov)
  if (M>1){
    param.est=apply(mi.param.est,2,mean)
    est.diff=scale(mi.param.est , center = TRUE, scale = FALSE)
    Sample.cov.param1=NULL
    for (i in 1:M){
      Sample.cov.param1[[i]]=(est.diff[i,])%*%t(est.diff[i,])
    }
    sample.cov.param1=Reduce('+',Sample.cov.param1)/ (M-1)
    Factor=(M+1)/M
    sample.cov.param = Factor* sample.cov.param1
    mean.cov.param=Reduce(`+`, mi.param.cov)/M
    param.cov=mean.cov.param+ sample.cov.param
  }
  if (M==1){
    param.est=mi.param.est
    param.cov=mi.param.cov
    sample.cov.param1=NULL
    mean.cov.param=NULL
  }
  
  return(list(param.est=param.est,param.cov=param.cov,between.cov=sample.cov.param1,within.cov=mean.cov.param))
}