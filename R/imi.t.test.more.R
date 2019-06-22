imi.t.test.more <- function(data.miss,data.imp0,max.M=500,epsilon,method='mvn',
                       x, y = NULL, alternative='two.sided',mu,paired = FALSE, var.equal = FALSE,
                       conf.level = 0.95,conv.plot=TRUE,successive.valid=3,
                       print.progress=TRUE){
  # dis.measure= "p.value", "half.ci", "df".
  # possible values for method:
  # all possible things in MICE,
  # in case mvn is chosen, then it switchs to amelia.
  # possibilities for dis.method and mah.scale
  # possible methods for distance
  # dis.method='euclidean'
  # dis.method='inf.norm'
  # dis.method='mahalanobis'
  # # possible methods for scale in Mahalanobis
  # mah.scale='within'
  # mah.scale='between'
  # mah.scale='combined'
  # successive.valid: length of distances successively smaller than epsilon to stop, minimum is zero
  if (length(mu)>1 ){
    stop('Please give a scalar mu')
  }


  # successive.valid: length of distances successively smaller than epsilon to stop, minimum is 1
  if (is.character(successive.valid)==FALSE){
    if (successive.valid<1 | floor(successive.valid)!=successive.valid){
      stop('Please enter a postive integer successive.valid')
    }
  }
  if (is.character(successive.valid)==TRUE){
    stop('Please enter a postive integer successive.valid')
  }
  # selecting the x part
  #miss.x=data.miss[-which(names(data.miss)==resp),]

  M=length(data.imp0)
  if (print.progress==TRUE){
    cat('We are working on performing the test on the initial imputed datasets.',fill = TRUE)
  }

  comb.mi0=mi.t.test(data.imp0, x=x, y = y, alternative = alternative,
                     mu = mu, paired = paired, var.equal = var.equal, conf.level = conf.level)
  data.imp=data.imp0
  #max.M=500
  dis.extra=epsilon+1
  dis.all=0
  while(sum(dis.extra>epsilon)>0 & M<max.M){
    M=M+1
    if (print.progress=='TRUE'){
      cat(paste('We are working on imputed datasets number ',M),fill = TRUE)
    }

    if (method=='mvn'){
      data.imp1=amelia(data.miss,m=1,p2s =0)$imputations$imp1
      data.imp[[M]]=data.imp1
    }
    if (method!='mvn' & method!='auto'){
      data.imp1_ini=mice(data.miss,m=1,method=method,print=FALSE)
      data.imp1=complete(data.imp1_ini, action = 1)
      data.imp[[M]]=data.imp1
    }
    if (method=='auto'){
      data.imp1_ini=mice(data.miss,m=1,print=FALSE)
      data.imp1=complete(data.imp1_ini, action = 1)
      data.imp[[M]]=data.imp1
    }
    # now fit the model to the new imputed data
    comb.mi1=mi.t.test(data.imp, x=x, y = y, alternative = alternative,
                       mu = mu, paired = paired, var.equal = var.equal, conf.level = conf.level)

    # now compute the distance between those two based on the specified method

    dis=sqrt((comb.mi1$p.value-comb.mi0$p.value)^2)
    comb.mi0=comb.mi1
    dis.all=c(dis.all,dis)
    dis.extra=tail(dis.all,n=successive.valid)
  }
  dis.all=dis.all[-1]
  if (conv.plot==TRUE){
    plot(dis.all,xlab='Number of imputations',ylab='Distance')
  }
  conv.status=0
  if (sum(dis.extra<epsilon)>0){
    conv.status=1
  }
  if (print.progress==TRUE){
    if (conv.status==1){
      cat(paste('We are done! The convergence is acheived and the sufficient number of imputations is ',M),
          fill = TRUE)
    }
    if (conv.status==0){
      cat('The convergence could not be achieved, please increase max.M or deacrese espsilon or
          number of successive validations steps.',
          fill = TRUE)
    }
  }
  return(list(test.result=comb.mi1,data.imp=data.imp,dis.steps=dis.all,conv.status=conv.status,M=M))

}
