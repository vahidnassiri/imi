imi.t.test <- function(data.miss,M0=2,max.M=500,epsilon,method='mvn',
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


  if (is.character(M0)==FALSE ){
    if (M0<2){
      stop('Please select an initial number of imputations larger than 1')
    }
  }
  if (is.character(M0)==FALSE ){
    if (floor(M0)!=M0){
      stop('Please select an integer initial number of imputations')
    }
  }
  if (is.character(successive.valid)==FALSE){
    if (successive.valid<1 | floor(successive.valid)!=successive.valid){
      stop('Please enter a postive integer successive.valid')
    }
  }
  if (M0!='manual' & successive.valid!='manual'){
    if (print.progress==TRUE){
      cat('We are working on the initial imputation.',fill = TRUE)
    }

    if (method=='mvn'){
      require(Amelia)
      data.imp0=amelia(data.miss,m=M0,p2s =0)$imputations
    }
    if (method!='mvn' & method!='auto'){
      require(mice)
      data.imp0_ini=mice(data.miss,m=M0,method=method,print=FALSE)
      data.imp0=lapply(seq(M0), function(i) complete(data.imp0_ini, action = i))
    }
    if (method=='auto'){
      require(mice)
      data.imp0_ini=mice(data.miss,m=M0,print=FALSE)
      data.imp0=lapply(seq(M0), function(i) complete(data.imp0_ini, action = i))
    }
    data.imp=data.imp0

    # first testing multivariate normality
    require(MKmisc)
    comb.mi0=mi.t.test(data.imp, x=x, y = y, alternative = alternative,
                       mu = mu, paired = paired, var.equal = var.equal, conf.level = conf.level)

    # from here

    # Now we should add imputations till we reach that point
    # variable to define
    #epsilon=0.05

    M=length(data.imp)
  }
  ## interactive part
  if (M0=='manual' & successive.valid=='manual'){
    t1=proc.time()
    if (method=='mvn'){
      require(Amelia)
      data.imp0=amelia(data.miss,m=1,p2s =0)$imputations
    }
    if (method!='mvn' & method!='auto'){
      require(mice)
      data.imp0_ini=mice(data.miss,m=1,method=method,print=FALSE)
      data.imp0=lapply(seq(1), function(i) complete(data.imp0_ini, action = i))
    }
    if (method=='auto'){
      require(mice)
      data.imp0_ini=mice(data.miss,m=1,print=FALSE)
      data.imp0=lapply(seq(1), function(i) complete(data.imp0_ini, action = i))
    }
    # first testing multivariate normality
    require(MKmisc)
    comb.mi0=mi.t.test(data.imp0, x=x, y = y, alternative = alternative,
                       mu = mu, paired = paired, var.equal = var.equal, conf.level = conf.level)

    t2=proc.time()



    print(paste('The time it takes (in seconds) to imput the data once and fit the model to it is:',round((t2-t1)[3],3)))

    M0.ini=readline('What is your choice of initial number of imputations?')
    successive.valid.ini=readline('What is your choice for successive steps validation?')
    if (print.progress==TRUE){
      cat('We are working on the initial imputation.',fill = TRUE)
    }

    M0 = as.numeric(unlist(strsplit(M0.ini, ",")))-1
    successive.valid=as.numeric(unlist(strsplit(successive.valid.ini, ",")))
    if (successive.valid<1 | floor(successive.valid)!=successive.valid){
      stop('Please enter a postive integer successive.valid')
    }
    data.imp=data.imp0
    if (floor(M0+1)!=M0+1){
      stop('Please select an integer initial number of imputations')
    }
    if ((M0+1)<2){
      stop('Please select an initial number of imputations larger than 1')
    }
    if (method=='mvn'){
      require(Amelia)
      data.imp0=amelia(data.miss,m=M0,p2s =0)$imputations
    }
    if (method!='mvn' & method!='auto'){
      require(mice)
      data.imp0_ini=mice(data.miss,m=M0,method=method,print=FALSE)
      data.imp0=lapply(seq(M0), function(i) complete(data.imp0_ini, action = i))
    }
    if (method=='auto'){
      require(mice)
      data.imp0_ini=mice(data.miss,m=M0,print=FALSE)
      data.imp0=lapply(seq(M0), function(i) complete(data.imp0_ini, action = i))
    }
    data.imp[2:(M0+1)]=data.imp0

    # first testing multivariate normality
    require(MKmisc)
    comb.mi0=mi.t.test(data.imp, x=x, y = y, alternative = alternative,
                       mu = mu, paired = paired, var.equal = var.equal, conf.level = conf.level)

    # from here

    # Now we should add imputations till we reach that point
    # variable to define
    #epsilon=0.05

    M=length(data.imp)

  }
  # Now for one of them
  if (M0=='manual' & successive.valid!='manual'){
    t1=proc.time()
    if (method=='mvn'){
      require(Amelia)
      data.imp0=amelia(data.miss,m=1,p2s =0)$imputations
    }
    if (method!='mvn' & method!='auto'){
      require(mice)
      data.imp0_ini=mice(data.miss,m=1,method=method,print=FALSE)
      data.imp0=lapply(seq(1), function(i) complete(data.imp0_ini, action = i))
    }
    if (method=='auto'){
      require(mice)
      data.imp0_ini=mice(data.miss,m=1,print=FALSE)
      data.imp0=lapply(seq(1), function(i) complete(data.imp0_ini, action = i))
    }
    # first testing multivariate normality
    require(MKmisc)
    comb.mi0=mi.t.test(data.imp0, x=x, y = y, alternative = alternative,
                       mu = mu, paired = paired, var.equal = var.equal, conf.level = conf.level)

    t2=proc.time()

    print(paste('The time it takes (in seconds) to imput the data once and fit the model to it is:',round((t2-t1)[3],3)))

    M0.ini=readline('What is your choice of initial number of imputations?')
    if (print.progress==TRUE){
      cat('We are working on the initial imputation.',fill = TRUE)
    }

    M0 = as.numeric(unlist(strsplit(M0.ini, ",")))-1
    data.imp=data.imp0
    if (floor(M0+1)!=M0+1){
      stop('Please select an integer initial number of imputations')
    }
    if ((M0+1)<2){
      stop('Please select an initial number of imputations larger than 1')
    }
    if (method=='mvn'){
      require(Amelia)
      data.imp0=amelia(data.miss,m=M0,p2s =0)$imputations
    }
    if (method!='mvn' & method!='auto'){
      require(mice)
      data.imp0_ini=mice(data.miss,m=M0,method=method,print=FALSE)
      data.imp0=lapply(seq(M0), function(i) complete(data.imp0_ini, action = i))
    }
    if (method=='auto'){
      require(mice)
      data.imp0_ini=mice(data.miss,m=M0,print=FALSE)
      data.imp0=lapply(seq(M0), function(i) complete(data.imp0_ini, action = i))
    }
    data.imp[2:(M0+1)]=data.imp0

    # first testing multivariate normality
    require(MKmisc)
    comb.mi0=mi.t.test(data.imp, x=x, y = y, alternative = alternative,
                       mu = mu, paired = paired, var.equal = var.equal, conf.level = conf.level)

    # from here

    # Now we should add imputations till we reach that point
    # variable to define
    #epsilon=0.05

    M=length(data.imp)
  }
  # successive
  if (successive.valid=='manual' & M0!='manual'){
    t1=proc.time()
    if (method=='mvn'){
      require(Amelia)
      data.imp0=amelia(data.miss,m=1,p2s =0)$imputations
    }
    if (method!='mvn' & method!='auto'){
      require(mice)
      data.imp0_ini=mice(data.miss,m=1,method=method,print=FALSE)
      data.imp0=lapply(seq(1), function(i) complete(data.imp0_ini, action = i))
    }
    if (method=='auto'){
      require(mice)
      data.imp0_ini=mice(data.miss,m=1,print=FALSE)
      data.imp0=lapply(seq(1), function(i) complete(data.imp0_ini, action = i))
    }
    # first testing multivariate normality
    require(MKmisc)
    comb.mi0=mi.t.test(data.imp0, x=x, y = y, alternative = alternative,
                       mu = mu, paired = paired, var.equal = var.equal, conf.level = conf.level)

    t2=proc.time()

    print(paste('The time it takes (in seconds) to imput the data once and fit the model to it is:',round((t2-t1)[3],3)))

    successive.valid.ini=readline('What is your choice for successive steps validation?')
    if (print.progress==TRUE){
      cat('We are working on the initial imputation.',fill = TRUE)
    }

    successive.valid=as.numeric(unlist(strsplit(successive.valid.ini, ",")))
    if (successive.valid<1 | floor(successive.valid)!=successive.valid){
      stop('Please enter a postive integer successive.valid')
    }

    data.imp=data.imp0


    if (method=='mvn'){
      require(Amelia)
      data.imp0=amelia(data.miss,m=M0,p2s =0)$imputations
    }
    if (method!='mvn' & method!='auto'){
      require(mice)
      data.imp0_ini=mice(data.miss,m=M0,method=method,print=FALSE)
      data.imp0=lapply(seq(M0), function(i) complete(data.imp0_ini, action = i))
    }
    if (method=='auto'){
      require(mice)
      data.imp0_ini=mice(data.miss,m=M0,print=FALSE)
      data.imp0=lapply(seq(M0), function(i) complete(data.imp0_ini, action = i))
    }
    data.imp[2:(M0+1)]=data.imp0

    # first testing multivariate normality
    require(MKmisc)
    comb.mi0=mi.t.test(data.imp, x=x, y = y, alternative = alternative,
                       mu = mu, paired = paired, var.equal = var.equal, conf.level = conf.level)

    # from here

    # Now we should add imputations till we reach that point
    # variable to define
    #epsilon=0.05

    M=length(data.imp)


  }



  #method='mice'
  #method='mvn'

  # selecting the x part
  #miss.x=data.miss[-which(names(data.miss)==resp),]


  #max.M=500
  dis.extra=epsilon+1
  dis.all=0
  while(sum(dis.extra>epsilon)>0 & M<max.M){
    M=M+1
    if (print.progress=='TRUE'){
      cat(paste('We are working on imputed dataset number ',M),fill = TRUE)
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
