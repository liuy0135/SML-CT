mce=function(x,n1=50,n2=200){
  x1=x[1:n1];x2=x[(n1+1):(n1+n2)]
  s1=getmode(x1)
  po1=1-length(which(x1==s1))/n1
  
  s2=c(1:2)[-s1]
  po2=1-length(which(x2==s2))/n2
  (po1+po2)/2
  #a1=mean(abs(x-c(rep(1,n1),rep(2,n2))))
  #a2=mean(abs(x-c(rep(2,n1),rep(1,n2))))
  #min(a1,a2)
}


va=function(x){
  OutVals = boxplot(x,plot=FALSE)$out
  if (length(OutVals)==0){
    ifelse(var(x)==0,var(x), var(x)/sd((x-mean(x))^2))
    #var(x)/sd((x-mean(x))^2)
  }else
  { id=which(x %in% OutVals)
  y=x[-id]
  ifelse(var(y)==0,var(y), var(y)/sd((x-mean(y))^2))}
}

two_test2=function(x,y){
  OutVals1 = boxplot(x,plot=FALSE)$out
  OutVals2 = boxplot(y,plot=FALSE)$out
  x1=x;y1=y
  if (length(OutVals1)>0){
    id1=which(x %in% OutVals1)
    x1=x[-id1]        
  }
  if (length(OutVals2)>0){
    id2=which(y %in% OutVals2)
    y1=y[-id2]        
  }
  n1=length(x1);n2=length(y1)
  S=((n1-1)*var(x1)+(n2-1)*var(y1))/(n1+n2-2)
  T=mean(x1)-mean(y1)
  T*S^(-1)*T
}

p=128*128;n1=82;n2=281;
n=n1+n2
ll1=ll2=NULL
sam01=sample(1:n1,50);sam02=sample(1:n2,200)

####remove the noises
for(l in 1:51){
  dat_ran1=  dat_ran12=matrix(0,p,n1)
  for (k in 1:n1){
    temp=Pancreas_CT[[k]]
    sam=a1[[k]][l]
    sam12=a2[[k]][l]
    dat_ran1[,k]=temp[,sam]
    dat_ran12[,k]=temp[,sam12]
  }
  dat_ran2=dat_ran22=matrix(0,p,n2)
  for (k in 1:n2){
    temp2=Cancer[[k]]
    sam2=b1[[k]][l]
    sam22=b2[[k]][l]
    dat_ran2[,k]=temp2[,sam2]
    dat_ran22[,k]=temp2[,sam22]
  }
  dat_ran=cbind(dat_ran1,dat_ran2)
  dat_ran_02=cbind(dat_ran12,dat_ran22)
  
  x.train1=cbind(dat_ran1[,sam01],dat_ran2[,sam02])
  x.train2=cbind(dat_ran12[,sam01],dat_ran22[,sam02])
  
  va1=apply( x.train1,1,va)
  va2=apply( x.train2,1,va)
  
  idx_var1=which(va1==0)
  idx_var2=which(va2==0)
  
  
  uu1=uu2=NULL
  for(ss in 1:p){
    uu1[ss]=two_test2(dat_ran1[ss,sam01],dat_ran2[ss,sam02])
    uu2[ss]=two_test2(dat_ran12[ss,sam01],dat_ran22[ss,sam02])
  }
  reA2_01=setdiff(which(uu1>1/log(250)),idx_var1)
  reA2_02=setdiff(which(uu2>1/log(250)),idx_var1)
  
  S1=t(x.train1[reA2_01,])%*%x.train1[reA2_01,]
  S2=t(x.train2[reA2_02,])%*%x.train2[reA2_02,]
  
  q1=svd(S1)$u[,1:2]
  q2=svd(S2)$u[,1:2]
  k1=kmeans(q1[,1],2)$cluster
  k2=kmeans(q1[,2],2)$cluster
  k12=kmeans(q1,2)$cluster
  ll1[l]=min(mce(k1),mce(k2),mce(k12))
  
  k2_1=kmeans(q2[,1],2)$cluster
  k2_2=kmeans(q2[,2],2)$cluster
  k2_12=kmeans(q2,2)$cluster
  ll2[l]=min(mce(k2_1),mce(k2_2),mce(k2_12))
  print(l)
  #plot.ts(q2)
}
###alpha we need
id_B=which(rank(ll2)<3)