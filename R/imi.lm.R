imi.lm <-
  function(data.miss,
           M0 = 2,
           max.M = 500,
           epsilon,
           method = 'pmm',
           resp,
           regressors,
           conv.plot = TRUE,
           dis.method = 'mahalanobis',
           mah.scale = 'combined',
           successive.valid = 3,
           print.progress=TRUE) {
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

    if (is.character(successive.valid) == FALSE) {
      if (successive.valid < 1 |
          floor(successive.valid) != successive.valid) {
        stop('Please enter a postive integer successive.valid')
      }
    }
    if (is.character(M0) == FALSE) {
      if (M0 < 2) {
        stop('Please select an initial number of imputations larger than 1')
      }
    }
    if (is.character(M0) == FALSE) {
      if (floor(M0) != M0) {
        stop('Please select an integer initial number of imputations')
      }
    }
    ## nothing is manual
    if (M0 != 'manual' & successive.valid != 'manual') {
      if (print.progress==TRUE){
        cat('We are working on the initial imputation.',fill = TRUE)
      }

      if (method == 'mvn') {
        require(Amelia)
        data.imp0 = amelia(data.miss, m = M0,p2s =0)$imputations
      }
      if (method != 'mvn' & method!='auto') {
        require(mice)
        data.imp0_ini = mice(data.miss, m = M0, method = method,print=FALSE)
        data.imp0 = lapply(seq(M0), function(i)
          complete(data.imp0_ini, action = i))
      }
      if (method=='auto') {
        require(mice)
        data.imp0_ini = mice(data.miss, m = M0,print=FALSE)
        data.imp0 = lapply(seq(M0), function(i)
          complete(data.imp0_ini, action = i))
      }

      # Now fitting the modle on these
      # first creating the modle formula
      model.expression <-
        as.formula(paste(paste(resp, "~"), paste(regressors, collapse = "+")))
      # fit the model to the first one to find the length of parameters
      model.fit = lm(model.expression, data = data.imp0[[1]])
      param.est1 = model.fit$coefficients
      param.cov1 = vcov(model.fit)
      param.est = matrix(0, M0, length(param.est1))
      param.cov = NULL
      param.est[1, ] = param.est1
      param.cov[[1]] = param.cov1
      # now find the estimate for the rest
      for (i in 2:M0) {
        model.fit = lm(model.expression, data = data.imp0[[i]])
        param.est[i, ] = model.fit$coefficients
        param.cov[[i]] = vcov(model.fit)
      }
      comb.mi0 = combine.mi(param.est, param.cov)
      data.imp = data.imp0
      # Now we should add imputations till we reach that point
      # variable to define
      #epsilon=0.05

      M = M0
    }
    ## interactive part: both manual
    if (M0 == 'manual' & successive.valid == 'manual') {
      t1 = proc.time()
      if (method == 'mvn') {
        require(Amelia)
        data.imp0 = amelia(data.miss, m = 2,p2s =0)$imputations
      }
      if (method != 'mvn' & method!='auto') {
        require(mice)
        data.imp0_ini = mice(data.miss, m = 2, method = method,print=FALSE)
        data.imp0 = lapply(seq(2), function(i)
          complete(data.imp0_ini, action = i))
      }
      if (method=='auto') {
        require(mice)
        data.imp0_ini = mice(data.miss, m = 2,print=FALSE)
        data.imp0 = lapply(seq(2), function(i)
          complete(data.imp0_ini, action = i))
      }

      # Now fitting the modle on these
      # first creating the modle formula
      model.expression <-
        as.formula(paste(paste(resp, "~"), paste(regressors, collapse = "+")))
      # fit the model to the first one to find the length of parameters
      model.fit = lm(model.expression, data = data.imp0[[1]])
      param.est1 = model.fit$coefficients
      param.cov1 = vcov(model.fit)
      param.est = matrix(0, 2, length(param.est1))
      param.cov = NULL
      param.est[1, ] = param.est1
      param.cov[[1]] = param.cov1
      model.fit = lm(model.expression, data = data.imp0[[2]])
      param.est[2, ] = model.fit$coefficients
      param.cov[[2]] = vcov(model.fit)
      comb.mi0 = combine.mi(param.est, param.cov)
      data.imp = data.imp0
      M = length(data.imp)
      t2 = proc.time()

      print(
        paste(
          'The time it takes (in seconds) to imput the data two times and fit the model to them is:',
          round((t2 - t1)[3], 3)
        )
      )

      M0.ini = readline('What is your choice of initial number of imputations?')
      successive.valid.ini = readline('What is your choice for successive steps validation?')

      if (print.progress==TRUE){
        cat('We are working on the initial imputation.',fill = TRUE)
      }

      M0 = as.numeric(unlist(strsplit(M0.ini, ","))) - 2
      successive.valid = as.numeric(unlist(strsplit(successive.valid.ini, ",")))
      if (successive.valid < 1 |
          floor(successive.valid) != successive.valid) {
        stop('Please enter a postive integer successive.valid')
      }
      data.imp = data.imp0
      if (floor(M0 + 1) != M0 + 1) {
        stop('Please select an integer initial number of imputations')
      }
      if ((M0 + 2) < 2) {
        stop('Please select an initial number of imputations larger than 1')
      }
      if (M0 > 0) {
        if (method == 'mvn') {
          require(Amelia)
          data.imp0 = amelia(data.miss, m = M0,p2s =0)$imputations
        }
        if (method != 'mvn' & method!='auto') {
          require(mice)
          data.imp0_ini = mice(data.miss, m = M0, method = method,print=FALSE)
          data.imp0 = lapply(seq(M0), function(i)
            complete(data.imp0_ini, action = i))
        }
        if (method=='auto') {
          require(mice)
          data.imp0_ini = mice(data.miss, m = M0,print=FALSE)
          data.imp0 = lapply(seq(M0), function(i)
            complete(data.imp0_ini, action = i))
        }

        # Now fitting the modle on these
        # first creating the modle formula
        #model.expression<-as.formula(paste(paste(resp,"~"),paste(regressors,collapse="+")))
        # fit the model to the first one to find the length of parameters
        model.fit = lm(model.expression, data = data.imp0[[1]])
        param.est1 = model.fit$coefficients
        param.cov1 = vcov(model.fit)
        #param.est=matrix(0,M0,length(param.est1))
        #param.cov=NULL
        param.est = rbind(param.est, param.est1)
        param.cov[[1 + 2]] = param.cov1
        # now find the estimate for the rest
        if (M0 > 1) {
          for (i in 2:M0) {
            model.fit = lm(model.expression, data = data.imp0[[i]])
            param.est = rbind(param.est, model.fit$coefficients)
            param.cov[[i + 2]] = vcov(model.fit)
          }
        }

        comb.mi0 = combine.mi(param.est, param.cov)
        data.imp[3:(M0 + 2)] = data.imp0
        # Now we should add imputations till we reach that point
        # variable to define
        #epsilon=0.05

        M = length(data.imp)
      }

    }
    ## interactive part: M0 manual
    if (M0 == 'manual' & successive.valid != 'manual') {
      t1 = proc.time()
      if (method == 'mvn') {
        require(Amelia)
        data.imp0 = amelia(data.miss, m = 2,p2s =0)$imputations
      }
      if (method != 'mvn' & method!='auto') {
        require(mice)
        data.imp0_ini = mice(data.miss, m = 2, method = method,print=FALSE)
        data.imp0 = lapply(seq(2), function(i)
          complete(data.imp0_ini, action = i))
      }
      if (method=='auto') {
        require(mice)
        data.imp0_ini = mice(data.miss, m = 2,print=FALSE)
        data.imp0 = lapply(seq(2), function(i)
          complete(data.imp0_ini, action = i))
      }

      # Now fitting the modle on these
      # first creating the modle formula
      model.expression <-
        as.formula(paste(paste(resp, "~"), paste(regressors, collapse = "+")))
      # fit the model to the first one to find the length of parameters
      model.fit = lm(model.expression, data = data.imp0[[1]])
      param.est1 = model.fit$coefficients
      param.cov1 = vcov(model.fit)
      param.est = matrix(0, 2, length(param.est1))
      param.cov = NULL
      param.est[1, ] = param.est1
      param.cov[[1]] = param.cov1
      model.fit = lm(model.expression, data = data.imp0[[2]])
      param.est[2, ] = model.fit$coefficients
      param.cov[[2]] = vcov(model.fit)
      comb.mi0 = combine.mi(param.est, param.cov)
      data.imp = data.imp0
      M = length(data.imp)
      t2 = proc.time()

      print(
        paste(
          'The time it takes (in seconds) to imput the data two times and fit the model to them is:',
          round((t2 - t1)[3], 3)
        )
      )

      M0.ini = readline('What is your choice of initial number of imputations?')
      if (print.progress==TRUE){
        cat('We are working on the initial imputation.',fill = TRUE)
      }

      M0 = as.numeric(unlist(strsplit(M0.ini, ","))) - 2
      data.imp = data.imp0
      if (floor(M0 + 1) != M0 + 1) {
        stop('Please select an integer initial number of imputations')
      }
      if ((M0 + 2) < 2) {
        stop('Please select an initial number of imputations larger than 1')
      }
      if (M0 > 0) {
        if (method == 'mvn') {
          require(Amelia)
          data.imp0 = amelia(data.miss, m = M0,p2s =0)$imputations
        }
        if (method != 'mvn' & method!='auto') {
          require(mice)
          data.imp0_ini = mice(data.miss, m = M0, method = method,print=FALSE)
          data.imp0 = lapply(seq(M0), function(i)
            complete(data.imp0_ini, action = i))
        }
        if (method=='auto') {
          require(mice)
          data.imp0_ini = mice(data.miss, m = M0,print=FALSE)
          data.imp0 = lapply(seq(M0), function(i)
            complete(data.imp0_ini, action = i))
        }

        # Now fitting the modle on these
        # first creating the modle formula
        #model.expression<-as.formula(paste(paste(resp,"~"),paste(regressors,collapse="+")))
        # fit the model to the first one to find the length of parameters
        model.fit = lm(model.expression, data = data.imp0[[1]])
        param.est1 = model.fit$coefficients
        param.cov1 = vcov(model.fit)
        #param.est=matrix(0,M0,length(param.est1))
        #param.cov=NULL
        param.est = rbind(param.est, param.est1)
        param.cov[[1 + 2]] = param.cov1
        # now find the estimate for the rest
        if (M0 > 1) {
          for (i in 2:M0) {
            model.fit = lm(model.expression, data = data.imp0[[i]])
            param.est = rbind(param.est, model.fit$coefficients)
            param.cov[[i + 2]] = vcov(model.fit)
          }
        }

        comb.mi0 = combine.mi(param.est, param.cov)
        data.imp[3:(M0 + 2)] = data.imp0
        # Now we should add imputations till we reach that point
        # variable to define
        #epsilon=0.05

        M = length(data.imp)
      }

    }
    ## interactive part: successive.valid is manual
    if (M0 != 'manual' & successive.valid == 'manual') {
      t1 = proc.time()
      if (method == 'mvn') {
        require(Amelia)
        data.imp0 = amelia(data.miss, m = 2,p2s =0)$imputations
      }
      if (method != 'mvn' & method!='auto') {
        require(mice)
        data.imp0_ini = mice(data.miss, m = 2, method = method,print=FALSE)
        data.imp0 = lapply(seq(2), function(i)
          complete(data.imp0_ini, action = i))
      }
      if (method=='auto') {
        require(mice)
        data.imp0_ini = mice(data.miss, m = 2,print=FALSE)
        data.imp0 = lapply(seq(2), function(i)
          complete(data.imp0_ini, action = i))
      }

      # Now fitting the modle on these
      # first creating the modle formula
      model.expression <-
        as.formula(paste(paste(resp, "~"), paste(regressors, collapse = "+")))
      # fit the model to the first one to find the length of parameters
      model.fit = lm(model.expression, data = data.imp0[[1]])
      param.est1 = model.fit$coefficients
      param.cov1 = vcov(model.fit)
      param.est = matrix(0, 2, length(param.est1))
      param.cov = NULL
      param.est[1, ] = param.est1
      param.cov[[1]] = param.cov1
      model.fit = lm(model.expression, data = data.imp0[[2]])
      param.est[2, ] = model.fit$coefficients
      param.cov[[2]] = vcov(model.fit)
      comb.mi0 = combine.mi(param.est, param.cov)
      data.imp = data.imp0
      M = length(data.imp)
      t2 = proc.time()

      print(
        paste(
          'The time it takes (in seconds) to imput the data two times and fit the model to them is:',
          round((t2 - t1)[3], 3)
        )
      )

      successive.valid.ini = readline('What is your choice for successive steps validation?')

      if (print.progress==TRUE){
        cat('We are working on the initial imputation.',fill = TRUE)
      }

      successive.valid = as.numeric(unlist(strsplit(successive.valid.ini, ",")))
      if (successive.valid < 1 |
          floor(successive.valid) != successive.valid) {
        stop('Please enter a postive integer successive.valid')
      }
      data.imp = data.imp0
      if (floor(M0 + 1) != M0 + 1) {
        stop('Please select an integer initial number of imputations')
      }
      if ((M0 + 2) < 2) {
        stop('Please select an initial number of imputations larger than 1')
      }
      if (M0 > 0) {
        if (method == 'mvn') {
          require(Amelia)
          data.imp0 = amelia(data.miss, m = M0,p2s =0)$imputations
        }
        if (method != 'mvn' & method!='auto') {
          require(mice)
          data.imp0_ini = mice(data.miss, m = M0, method = method,print=FALSE)
          data.imp0 = lapply(seq(M0), function(i)
            complete(data.imp0_ini, action = i))
        }
        if (method=='auto') {
          require(mice)
          data.imp0_ini = mice(data.miss, m = M0,print=FALSE)
          data.imp0 = lapply(seq(M0), function(i)
            complete(data.imp0_ini, action = i))
        }

        # Now fitting the modle on these
        # first creating the modle formula
        #model.expression<-as.formula(paste(paste(resp,"~"),paste(regressors,collapse="+")))
        # fit the model to the first one to find the length of parameters
        model.fit = lm(model.expression, data = data.imp0[[1]])
        param.est1 = model.fit$coefficients
        param.cov1 = vcov(model.fit)
        #param.est=matrix(0,M0,length(param.est1))
        #param.cov=NULL
        param.est = rbind(param.est, param.est1)
        param.cov[[1 + 2]] = param.cov1
        # now find the estimate for the rest
        if (M0 > 1) {
          for (i in 2:M0) {
            model.fit = lm(model.expression, data = data.imp0[[i]])
            param.est = rbind(param.est, model.fit$coefficients)
            param.cov[[i + 2]] = vcov(model.fit)
          }
        }

        comb.mi0 = combine.mi(param.est, param.cov)
        data.imp[3:(M0 + 2)] = data.imp0
        # Now we should add imputations till we reach that point
        # variable to define
        #epsilon=0.05

        M = length(data.imp)
      }

    }






    # Doing the rest of initial imputations

    #max.M=500
    dis.extra = epsilon + 1
    dis.all = 0
    while (sum(dis.extra > epsilon) > 0 & M < max.M) {
      M = M + 1
      if (print.progress=='TRUE'){
        cat(paste('We are working on imputed dataset number ',M),fill = TRUE)
      }

      if (method == 'mvn') {
        data.imp1 = amelia(data.miss, m = 1,p2s =0)$imputations$imp1
        data.imp[[M]] = data.imp1
      }
      if (method != 'mvn' & method!='auto') {
        data.imp1_ini = mice(data.miss, m = 1, method = method,print=FALSE)
        data.imp1 = complete(data.imp1_ini, action = 1)
        data.imp[[M]] = data.imp1
      }
      if (method=='auto') {
        data.imp1_ini = mice(data.miss, m = 1,print=FALSE)
        data.imp1 = complete(data.imp1_ini, action = 1)
        data.imp[[M]] = data.imp1
      }
      # now fit the model to the new imputed data
      model.fit = lm(model.expression, data = data.imp1)
      param.est = rbind(param.est, model.fit$coefficients)
      param.cov[[M]] = vcov(model.fit)
      comb.mi1 = combine.mi(param.est, param.cov)
      # now compute the distance between those two based on the specified method
      diff.est = comb.mi1$param.est - comb.mi0$param.est
      if (dis.method == 'euclidean') {
        dis = sqrt(t(diff.est) %*% diff.est)
      }
      if (dis.method == 'inf.norm') {
        dis = max(abs(diff.est))
      }
      if (dis.method == 'mahalanobis') {
        require('MASS')
        if (mah.scale == 'within') {
          #S.mah=comb.mi0$within.cov + comb.mi1$within.cov
          S.mah = comb.mi1$within.cov
        }
        if (mah.scale == 'between') {
          #S.mah=comb.mi0$between.cov + comb.mi1$between.cov
          S.mah = comb.mi1$between.cov
        }
        if (mah.scale == 'combined') {
          #S.mah=comb.mi0$param.cov + comb.mi1$param.cov
          S.mah = comb.mi1$param.cov
        }
        dis = sqrt((t(diff.est) %*% ginv(S.mah)) %*% diff.est)
      }
      comb.mi0 = comb.mi1
      dis.all = c(dis.all, dis)
      dis.extra = tail(dis.all, n = successive.valid)
    }
    dis.all = dis.all[-1]
    if (conv.plot == TRUE) {
      plot(dis.all, xlab = 'Number of imputations', ylab = 'Distance')
    }
    conv.status = 0
    if (sum(dis.extra < epsilon) > 0) {
      conv.status = 1
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
    return(
      list(
        mi.param = comb.mi1,
        data.imp = data.imp,
        dis.steps = dis.all,
        conv.status = conv.status,
        M = M
      )
    )
  }

