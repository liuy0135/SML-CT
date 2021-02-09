
compa_step2=function(x,gap=0.0025){
  se=seq(0,0.02,gap)
  q1=quantile(x,se)
  t1=NULL
  for (k in 1:length(se)){
    tmp1=abs(x-q1[k])
    t1[k]=which(tmp1==min(tmp1))[1]
  }
  t1
}

a1_s2=a2_s2=list()
for(k in 1:82){
  dat_ran=as.matrix(Pancreas_CT[[k]])
  S=t(dat_ran)%*%dat_ran
  q1=svd(S)$u[,1:2]
  sign01=ifelse(cal_area(sort(q1[,1]))>0.5,1,-1)
  sign02=ifelse(cal_area(sort(q1[,2]))>0.5,1,-1)
  a1_s2[[k]]=compa_step2(sign01*q1[,1])
  a2_s2[[k]]=compa_step2(sign02*q1[,2])
  print(k)
}

b1_s2=b2_s2=list()
for(k in 1:281){
  dat_ran=as.matrix(Cancer[[k]])
  S=t(dat_ran)%*%dat_ran
  q1=svd(S)$u[,1:2]
  sign01=ifelse(cal_area(sort(q1[,1]))>0.5,1,-1)
  sign02=ifelse(cal_area(sort(q1[,2]))>0.5,1,-1)
  b1_s2[[k]]=compa_step2(sign01*q1[,1])
  b2_s2[[k]]=compa_step2(sign02*q1[,2])
  print(k)
}



HtoC=CtoH=CR=NULL

for (tt in 1:50){
  print(tt)
p=128*128;n1=82;n2=281;
n=n1+n2
sam01=sample(1:n1,50);sam02=sample(1:n2,200)


###alpha we need
id_B=which(rank(ll2)<6)#c(41)#31#c(1:51)#51#

p_or=128*128;n2_or=281
data_tr1=matrix(0,p_or,n1)
for (k in 1:n1){
  temp1=Pancreas_CT[[k]]
  id=a2_s2[[k]][id_B]#sample(1:ncol(temp1),3)#
  data_tr1[,k]=apply(temp1[,id],1,mean)
  #data_tr1[,k]=apply(temp1,1,mean)
  #data_tr1[,k]=temp1[,id]
}
data_tr2=matrix(0,p_or,n2_or)
for (k in 1:n2){
  temp1=Cancer[[k]]
  id=b2_s2[[k]][id_B]#sample(1:ncol(temp1),3)#
  data_tr2[,k]=apply(temp1[,id],1,mean)
  #data_tr2[,k]=apply(temp1,1,mean)
  #data_tr2[,k]=temp1[,id]
}


idx_var1=which(va1==0)

data_tr=cbind(data_tr1[,sam01], data_tr2[,sam02])
va1=apply( data_tr,1,va)
idx_var1=which(va1==0)
uuf=NULL
for(ss in 1:p){
  uuf[ss]=two_test2(data_tr1[ss,sam01],data_tr2[ss,sam02])
  #uuf[ss]=length(which(c(data_tr1[ss,sam01],data_tr2[ss,sam02])==-1024))
}
reA2_ff=setdiff(which(uuf>1/log(250)),idx_var1)

length(idx_var1)/p_or;1-length(idx_var1)/p_or-length(reA2_ff)/p_or;length(reA2_ff)/p_or;

data_te=cbind(data_tr1[,-sam01],data_tr2[,-sam02])
y.train=rep(0,250)
y.train[1:50]=1
y.test=rep(0,113)
y.test[1:32]=1


xx=data.frame(t(data_tr[reA2_ff,]))
xx.test=data.frame(t(data_te[reA2_ff,]))


########RF-GB
m=20;ntrees=1000;v=0.1;M=100
##tes=random_forest(xx,y.train,xx.test,m=15,ntrees=2000)
####
n=length(y.train);p=ncol(xx);n2=nrow(xx.test)
#sam_rep=matrix(0,n,ntrees)
#sam_feature=matrix(0,m,ntrees)
RE=matrix(0,n2,ntrees)

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
    
    pred.prop=pp
    res1=y.new-exp(beta)/(1+exp(beta))
    new_observation=data.frame(x.new,res1)
    base_tree<-tree(res1~., data=new_observation)
    #if(length(base_tree$frame$var))
    ###choose the best stumps of trees
    # print(s)
    prune_tree=base_tree
    try(
      if( length(which(base_tree$where==2))>2){
        #try(cv_base = cv.tree(base_tree, FUN = prune.tree))
        cv_base = cv.tree(base_tree, FUN = prune.tree)
        #if (length(cv_base$dev)>3){
        ratio=cv_base$dev[1:8]/cv_base$dev[2:9]
        prune_tree = prune.tree(base_tree, best = max(which(ratio==max(ratio)),2))#}
      } )
    
    tree.pred = predict(prune_tree, new_observation)
    
    ###extract the leaf nodes of a tree
    ta1=prune_tree$frame
    leaf_node=ta1$yval[which(ta1$var=="<leaf>")]
    len=length(leaf_node)
    ###calculate gamma
    gama=NULL
    gam1=rep(0,length(y.new))
    for (j in 1: len){
      stump=as.numeric(which(tree.pred==leaf_node[j]))
      #obs=new_observation[stump,]
      gama[j]=sum(res1[stump])/sum(pp[stump]*(1-pp[stump]))
      gam1[stump]=sum(res1[stump])/sum(pp[stump]*(1-pp[stump]))
    }
    beta=beta+v*gam1
    pp=exp(beta)/(1+exp(beta))
    #####for prediction
    test_pred = predict(prune_tree,xx.test[,sam_feature])
    for (i in 1: len){
      stump_test=as.numeric(which(test_pred==leaf_node[i]))
      prediction1[stump_test]=prediction1[stump_test]+v*gama[i]
    }
    
  }
  #pred.prop1=exp(prediction1)/(1+exp((prediction1)))
  ##
  RE[,k]=exp(prediction1)/(1+exp((prediction1)))
}



tes=apply(RE, 1,mean)
y.test=c(rep(1,32),rep(0,81))
ttt=ifelse(tes>0.5,1,0)-y.test
CR[tt]=1-length(which((ifelse(tes>0.5,1,0)-y.test)==0))/length(y.test)
HtoC[tt]=length(which(ttt[1:31]!=0))
CtoH[tt]=length(which(ttt[32:113]!=0))
}
mean(CR)
mean( HtoC)
mean(CtoH)
