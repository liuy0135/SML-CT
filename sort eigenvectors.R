setwd("/Users/lym/Google Drive/cancer data")
#load("Pancreas_CT_health")
#load("Pancreas_CT_Cancer")
cal_area=function(x){
  library(MESS)
  n=length(x)
  gap=max(x)-min(x)
  x_ax=c(1:n)/n;y_ax=(x-min(x))/gap
  auc(x_ax,y_ax)
}

compa=function(x,gap=0.02){
  se=seq(0,1,gap)
  q1=quantile(x,se)
  t1=NULL
  for (k in 1:length(se)){
    tmp1=abs(x-q1[k])
    t1[k]=which(tmp1==min(tmp1))[1]
  }
  t1
}

a1=a2=list()
for(k in 1:82){
  dat_ran=as.matrix(Pancreas_CT[[k]])
  S=t(dat_ran)%*%dat_ran
  q1=svd(S)$u[,1:2]
  sign01=ifelse(cal_area(sort(q1[,1]))>0.5,1,-1)
  sign02=ifelse(cal_area(sort(q1[,2]))>0.5,1,-1)
  a1[[k]]=compa(sign01*q1[,1])
  a2[[k]]=compa(sign02*q1[,2])
  print(k)
}

b1=b2=list()
for(k in 1:281){
  dat_ran=as.matrix(Cancer[[k]])
  S=t(dat_ran)%*%dat_ran
  q1=svd(S)$u[,1:2]
  sign01=ifelse(cal_area(sort(q1[,1]))>0.5,1,-1)
  sign02=ifelse(cal_area(sort(q1[,2]))>0.5,1,-1)
  b1[[k]]=compa(sign01*q1[,1])
  b2[[k]]=compa(sign02*q1[,2])
  print(k)
}
