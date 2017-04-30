### this function is modified from gwr.predict in GWmodel package ###

gwr.predict.mod = function(formula, data, predictdata, bw, kernel, adaptive = FALSE, generate.prediction.var=F) {
  
  fd.locat <- coordinates(data) ### extract training co-ordinate
  fd.n <- nrow(fd.locat) ### number of training co-ordinate
  
  pd.locat <- coordinates(predictdata) ### extract test co-ordinate
  pd.n <- nrow(pd.locat) ### number of test co-ordinate
  
  #####################################################################################################################################
  ### THIS WHOLE CHUNK IS TO EXTRACT OUT THE PROPER X (independent variable) AND Y (dependent variable) VALUES from the training data #
  #####################################################################################################################################
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.extract(mf, "response")
  x <- model.matrix(mt, mf)
  
  var.n <- ncol(x) ### number of predictor (including intercept)
  idx1 <- match("(Intercept)", colnames(x)) ### identify intercept column
  if (!is.na(idx1)) {colnames(x)[idx1] <- "Intercept"} ### some renaming
  inde_vars <- colnames(x)[-1] ### names of independent variables (excluding intercept column)

  
  ########################################
  ### THIS IS TO EXTRACT THE TEST DATA ###
  ########################################
  predictdata <- as(predictdata, "data.frame") ### convert the test data to dataframe
  x.p <- predictdata[, inde_vars] ### extract the test indepedent variables
  x.p <- cbind(rep(1, pd.n), x.p) ### add intercept ("1")
  x.p <- as.matrix(x.p)
  
  
  ###########################################
  #### PART 1: FOCUS ON TEST DATA ###########
  #### THIS IS THE ACTUAL PREDICTION STEP ###
  ###########################################
  
  ### matrix to hold the weights
  wt <- matrix(nrow = fd.n, ncol = pd.n)
  
  ### matrix to hold xtxinv (i.e. inverse of t(x) %*% t)
  xtxinv <- as.double(rep(0, pd.n * var.n * var.n))
  dim(xtxinv) <- c(pd.n, var.n, var.n)
  
  ### matrix to hold beta
  betas1 <- matrix(nrow = pd.n, ncol = var.n)
  
  ### for each set of test co-ordinate
  for (i in 1:pd.n) {
    
    ### get GW distance between train co-ordinate and test co-ordinate
    dist.vi <- gw.dist(fd.locat, pd.locat, focus = i)
    
    ### assign corresponding GW weights
    W.i <- gw.weight(dist.vi, bw, kernel, adaptive)
    
    ### record the weight (W.i) in wt
    wt[, i] <- W.i
    
    ### GW regression on x and y using weight (W.i)
    gw.resi <- gw.reg1(x, y, W.i, i)
    
    ### record the beta values
    betas1[i, ] <- gw.resi[[1]]
    
    ### record the xtxinv
    xtxinv[i, , ] <- gw.resi[[2]]
  }
  
  ### THIS IS THE PREDICTION STEP!!!
  ### predict using x.p (predictor independent variables) and betas1 (previously recorded)
  gw.predict <- gwr.fitted(x.p, betas1)
  
  
  if (generate.prediction.var) {
    
    ##########################################################################
    ### Everything below are for generating the regression statistics only ###
    ##########################################################################
    
    #######################################
    #### PART 2: FOCUS ON TRAINING DATA ###
    #######################################
    
    ### matrix to hold beta
    betas2 <- matrix(nrow = fd.n, ncol = var.n)
    
    ### S is analogous to xtxinv above, but the form is much simpler
    S <- matrix(nrow = fd.n, ncol = fd.n)
    
    for (i in 1:fd.n) {
      dist.vi <- gw.dist(fd.locat, fd.locat, focus = i)
      W.i <- gw.weight(dist.vi, bw, kernel, adaptive)
      gw.resi <- gw.reg(x, y, W.i, T, i)
      betas2[i, ] <- gw.resi[[1]]
      S[i, ] <- gw.resi[[2]]
    }
    
    tr.S <- sum(diag(S))
    tr.StS <- sum(S^2)
    Q <- t(diag(fd.n) - S) %*% (diag(fd.n) - S)
    RSS.gw <- t(y) %*% Q %*% y
    sigma.hat <- RSS.gw/(fd.n - 2 * tr.S + tr.StS)
    
    
    predict.var <- numeric(pd.n)
    for (i in 1:pd.n) {
      w2 <- wt[, i] * wt[, i]
      w2x <- x * w2
      xtw2x <- t(x) %*% w2x
      xtxinvp <- c()
      for (j in 1:var.n) xtxinvp <- rbind(xtxinvp, xtxinv[i, , j])
      s0 <- xtxinvp %*% xtw2x %*% xtxinvp
      x.pi <- x.p[i, ]
      dim(x.pi) <- c(1, length(x.pi))
      s1 <- x.pi %*% s0 %*% t(x.pi)
      pse <- (sqrt(sigma.hat)) * (sqrt(1 + s1))
      pvar <- (pse)^2
      predict.var[i] <- pvar
    }
    
    gwr.pred.df <- data.frame(betas1, gw.predict, predict.var)
  } else {
    gwr.pred.df <- data.frame(betas1, gw.predict, NA)
  }
  
  

  
  ### clean up and output data
  colnames(gwr.pred.df) <- c(paste(colnames(x), "coef", sep = "_"), "prediction", "prediction_var")
  
  form =  as.character(formula)
  form = form[2]
  gwr.pred.df = cbind(gwr.pred.df, actual=predictdata[,form])
  rownames(pd.locat) <- rownames(gwr.pred.df)
  SDF = SpatialPointsDataFrame(coords = pd.locat, data = gwr.pred.df)
  return(SDF)

  
}
