imi.t.test.conv <- function(data.imp,epsilon,x, y = NULL, alternative='two.sided',mu,paired = FALSE, var.equal = FALSE,
                       conf.level = 0.95,conv.plot=TRUE,successive.valid=3){

  y.var=NULL
  M=length(data.imp)
  dis=rep(0,M-1)
  param.cov=NULL
  # for the first one
  data.imp1=data.imp[[1]]
  x.num=which(names(data.imp1)==x)
  if (is.null(y)==0){
    y.num=which(names(data.imp1)==y)
    y.var=data.imp1[,y.num]
  }
  
  comb.mi0=t.test(x=data.imp1[,x.num], y = y.var, alternative = alternative,
                     mu = mu, paired = paired, var.equal = var.equal, conf.level = conf.level)
   for (i in 2:M){
    data.imp1=data.imp[1:i]
    comb.mi1=mi.t.test(data.imp1, x=x, y = y, alternative = alternative,
                       mu = mu, paired = paired, var.equal = var.equal, conf.level = conf.level)
    # now compute the distance between those two based on the specified method
    dis[i-1]=sqrt((comb.mi1$p.value-comb.mi0$p.value)^2)
    comb.mi0=comb.mi1
    #dis.all=c(dis.all,dis)
    #dis.extra=tail(dis.all,n=successive.valid)
  }
  crit.vec=dis<epsilon
  diff.cumsum=diff(cumsum(crit.vec),lag=successive.valid)
  sufficient.M='Not sufficient!'
  if (successive.valid %in% diff.cumsum){
    sufficient.M=which(diff.cumsum==successive.valid)[1]+(successive.valid+1)
  }
  if (conv.plot==TRUE){
    plot(2:M,dis,xlab='Number of imputations',ylab='Distance')
  }
  return(list(dis.steps=dis,sufficient.M=sufficient.M))

 
}
