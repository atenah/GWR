gwr.bootstrap = function(formula.reg, data, output.folder, bootstrap, percentage.test, kernel.type, bandwidth=NULL) {
  
  data.full = data
  
  out.fold = file.path(output.folder)

  ### calculate GW matrix
  print("Calculate GW distance matrix...")
  DM = gw.dist(dp.locat=coordinates(data.full))
  
  ### select best bw
  if (is.null(bandwidth)) {
    print("Select best bandwidth...")
    bw.best = bw.gwr(formula(formula.reg), data=data.full, kernel = kernel.type, dMat = DM)
  } else {
    bw.best = bandwidth
  }
  write(paste("Bandwith:", bw.best), paste(out.fold, "/bootstrap_log.txt", sep=""), append = T)
  
  ### bootstrap
  print("Start Bootstrap procedure...")
  num.data = nrow(data.full)
  write(paste("Number of total data:", num.data), paste(out.fold, "/bootstrap_log.txt", sep=""), append = T)
  write(paste("Percentage of test data:", percentage.test, "%"), paste(out.fold, "/bootstrap_log.txt", sep=""), append = T)
  num.test = round(num.data*percentage.test/100, 0)
  write(paste("Number of test data:", num.test), paste(out.fold, "/bootstrap_log.txt", sep=""), append = T)
  
  row.numbering = 1:num.data
  
  bootstrap.raw = list()
  
  progbar = utils::txtProgressBar(style = 3)
  for (i in 1:bootstrap) {
    
    test.numbering = sort(sample(row.numbering, num.test))
    train.numbering = row.numbering[!row.numbering %in% test.numbering]
    
    ### training
    data.train = data.full[train.numbering,]
    
    ### test
    data.test = data.full[test.numbering,]
    
    ### use GW regression (simplified function) for prediction
    temp =  gwr.predict.mod.quick(formula(formula.reg), data = data.train, predictdata = data.test, bw = bw.best, kernel = kernel.type, DM, train.numbering, test.numbering)
    
    bootstrap.raw[[i]] = as.matrix(data.frame(temp)[,c("prediction", "actual")])
    
    utils::setTxtProgressBar(progbar, i/bootstrap)
  }
  close(progbar)
  
  print("Saving bootstrap results...")
  save(bootstrap.raw, file = paste(out.fold, "/bootstrap_raw.RObj", sep=""))
  
  print("Calculate summary statistics...")
  bootstrap.result = do.call(rbind, lapply(bootstrap.raw, function(x) {
    t.ssq = sum((x[,2]-mean(x[,2]))^2)
    r.ssq = sum(abs(x[,1]-x[,2])^2)
    R2.model = 1 - r.ssq/t.ssq
    temp = cor.test(x[,1], x[,2])
    return(c(p=temp$p.value, cor=as.numeric(temp$est), R2.model=R2.model))
  }))
  bootstrap.result = data.frame(bootstrap.result)
  write.table(bootstrap.result, file = paste(out.fold, "/bootstrap_result.txt", sep=""), sep="\t", quote=F, row.names=F, col.names=T)
  save(bootstrap.result, file = paste(out.fold, "/bootstrap_result.RObj", sep=""))
  
  bootstrap.summary = summary(bootstrap.result)
  write.table(bootstrap.summary, file = paste(out.fold, "/bootstrap_summary.txt", sep=""), quote=F)

  print("DONE!!!")
  
}

