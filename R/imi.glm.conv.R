imi.glm.conv <- function(data.imp,family=binomial(link='logit'),epsilon,resp,regressors,
                         conv.plot=TRUE,
                         dis.method='mahalanobis',mah.scale='combined',successive.valid=3,max.iter.glm=1000){
  # possible methods for distance
  # dis.method='euclidean'
  # dis.method='inf.norm'
  # dis.method='mahalanobis'
  # # possible methods for scale in Mahalanobis
  # mah.scale='within'
  # mah.scale='between'
  # mah.scale='combined'
  # successive.valid: length of distances successively smaller than epsilon to stop, minimum is 1

  M=length(data.imp)
  dis=rep(0,M-1)
  param.cov=NULL
  # for the first one
  data.imp1=data.imp[[1]]
  model.expression<-as.formula(paste(paste(resp,"~"),paste(regressors,collapse="+")))
  model.fit=glm(model.expression,family=family,data=data.imp1,control = list(maxit = max.iter.glm))
  param.est=model.fit$coefficients
  param.cov[[1]]=vcov(model.fit)
  comb.mi0=combine.mi(param.est,param.cov)
  for (i in 2:M){
    data.imp1=data.imp[[i]]
    model.fit=glm(model.expression,family=family,data=data.imp1,control = list(maxit = max.iter.glm))
    param.est=rbind(param.est,model.fit$coefficients)
    param.cov[[i]]=vcov(model.fit)
    comb.mi1=combine.mi(param.est,param.cov)
    # now compute the distance between those two based on the specified method
    diff.est=comb.mi1$param.est-comb.mi0$param.est
    if (dis.method=='euclidean'){
      dis[i-1]=sqrt(t(diff.est)%*%diff.est)
    }
    if (dis.method=='inf.norm'){
      dis[i-1]=max(abs(diff.est))
    }
    if (dis.method=='mahalanobis'){
      # to compute generalized inverse
      require('MASS')
      if (mah.scale=='within'){
        #S.mah=comb.mi0$within.cov + comb.mi1$within.cov
        S.mah= comb.mi1$within.cov

      }
      if (mah.scale=='between'){
        #S.mah=comb.mi0$between.cov + comb.mi1$between.cov
        S.mah= comb.mi1$between.cov

      }
      if (mah.scale=='combined'){
        #S.mah=comb.mi0$param.cov + comb.mi1$param.cov
        S.mah=comb.mi1$param.cov
      }
      dis[i-1]=sqrt((t(diff.est)%*%ginv(S.mah))%*%diff.est)
    }
    comb.mi0=comb.mi1
    #dis.all=c(dis.all,dis)
    #dis.extra=tail(dis.all,n=successive.valid)
  }
  crit.vec=dis<epsilon
  diff.cumsum=diff(cumsum(crit.vec),lag=successive.valid)
  sufficient.M='Not sufficient!'
  if (successive.valid %in% diff.cumsum){
    sufficient.M=which(diff.cumsum==successive.valid)[1]+(successive.valid-1)
  }
  if (conv.plot==TRUE){
    plot(2:M,dis,xlab='Number of imputations',ylab='Distance')
  }
  return(list(dis.steps=dis,sufficient.M=sufficient.M))}
