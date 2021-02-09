

ptm <- proc.time()
RE=matrix(0,1,ntrees)
for (k in 1:ntrees){
  # print(k)
  set.seed(k)
  sam_rep=sample(c(1:n),n, replace = TRUE)
  sam_feature=sample(c(1:p),m,replace = FALSE)
  x.new=xx[sam_rep,sam_feature];y.new=y.train[sam_rep]
  ###gb###gb(x.new,y.new,x.test[,sam_feature],v=0.1,M=100)$predicted.prop
  p0=mean(y.new);odd0=p0/(1-p0);beta0=log(odd0)
  beta=beta0;pp=rep(p0,length(y.new));n=length(y.new)
  prediction1=rep(beta0,nrow(xx.test[,sam_feature]))
  ### step 2: iterations
  for (s in 1:M){
    tp=(k-1)*M+s
    #####for prediction
    test_pred = predict(prune_tree1[[tp]],xx.test[,sam_feature])
    for (i in 1: len1[tp]){
      stump_test=as.numeric(which(test_pred==leaf_node1[[tp]][i]))
      prediction1[stump_test]=prediction1[stump_test]+v*gama1[[tp]][i]
    }
    
  }
  #pred.prop1=exp(prediction1)/(1+exp((prediction1)))
  ##
  RE[,k]=exp(prediction1)/(1+exp((prediction1)))
}

tes=apply(RE, 1,mean)
tes

proc.time() - ptm
