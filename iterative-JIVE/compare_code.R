setwd("D:/git-project/multi-source/iterative-JIVE")
para = matrix(0,nrow=0,ncol=8)
para = rbind(para,c(80,100,10,1,1,100,1,5)) 
para = rbind(para,c(80,100,10,1,1,1000,1,5)) 
para = rbind(para,c(80,100,10,10,1,10,1,10))
para = rbind(para,c(80,100,10,10,1,100,1,10))
para = rbind(para,c(80,100,10,100,1,100,1,100))
para = rbind(para,c(80,100,10,100,1,1000,1,100))

para = rbind(para,c(80,100,10,1,1,100,1,1))
para = rbind(para,c(80,100,10,1,1,1000,1,1))
para = rbind(para,c(80,100,10,100,1,100,1,1))
para = rbind(para,c(80,100,10,100,1,1000,1,1))
para = rbind(para,c(80,100,10,1000,1,100,1,1))
para = rbind(para,c(80,100,10,1000,1,1000,1,1))
para = rbind(para,c(80,100,10,1000,1,100,100,1))
para = rbind(para,c(80,100,10,1000,1,1000,1,100))
para = rbind(para,c(80,100,10,1,1,100,100,1))
para = rbind(para,c(80,100,10,1,1,1000,1,100))
para = rbind(para,c(80,100,10,1,1,100,1000,1))
para = rbind(para,c(80,100,10,1,1,1000,1,1000))
total.result = NULL
total.result.sd = NULL
for(case in seq(1,nrow(para))){
  G = 2
  n = para[case,1]
  p1 = para[case,2]
  p2 = para[case,3]
  p=p1+p2
  r = 1
  r1 = 1
  r2 = 1
  r.list = list(r1, r2)
  L = 20
  alpha = rep(para[case,4], r)
  alpha1 = rep(para[case,5], r1) #OLS: 0
  alpha2 = rep(para[case,6], r2) #OLS: 0
  file.name = paste0("multi_source_result_pls_n=",n,"_p1=",p1,"_p2=",p2,"_a=",alpha[1],"_a1=",alpha1[1],"_a2=",alpha2[1],"_scale1=",para[case,7],"_scale2=",para[case,8])
  load(file = paste0(file.name,".rData"))
  total.result = rbind(total.result,colMeans(result))
  total.result.sd = rbind(total.result.sd,apply(result,2,sd))
}

total.result.order = NULL
for(case in seq(1,nrow(para))){
  total.result.order = rbind(total.result.order,colnames(total.result)[order(total.result[case,])])
}
compare.result = total.result[,c("continuum.pls","pls","pls.scaled")]
compare.result.sd = total.result.sd[,c("continuum.pls","pls","pls.scaled")]
compare.result.order = NULL
for(case in seq(1,nrow(para))){
  compare.result.order = rbind(compare.result.order,colnames(compare.result)[order(compare.result[case,])])
}
View(total.result.order)
View(compare.result.sd)
View(compare.result)
View(compare.result.order)





para = matrix(0,nrow=0,ncol=8)
para = rbind(para,c(80,100,10,1,1,100,1,1))
para = rbind(para,c(80,100,10,1,1,1000,1,1))
para = rbind(para,c(80,100,10,100,1,100,1,1))
para = rbind(para,c(80,100,10,100,1,1000,1,1))
para = rbind(para,c(80,100,10,1000,1,100,1,1))
para = rbind(para,c(80,100,10,1000,1,1000,1,1))
para = rbind(para,c(80,100,10,1000,1,100,100,1))
para = rbind(para,c(80,100,10,1000,1,1000,1,100))
para = rbind(para,c(80,100,10,1,1,100,100,1))
para = rbind(para,c(80,100,10,1,1,1000,1,100))
para = rbind(para,c(80,100,10,1,1,100,1000,1))
para = rbind(para,c(80,100,10,1,1,1000,1,1000))
total.est.result = NULL
for(case in seq(1,nrow(para))){
  G = 2
  n = para[case,1]
  p1 = para[case,2]
  p2 = para[case,3]
  p=p1+p2
  r = 1
  r1 = 1
  r2 = 1
  r.list = list(r1, r2)
  L = 20
  alpha = rep(para[case,4], r)
  alpha1 = rep(para[case,5], r1) #OLS: 0
  alpha2 = rep(para[case,6], r2) #OLS: 0
  file.name = paste0("multi_source_est_result_pls_n=",n,"_p1=",p1,"_p2=",p2,"_a=",alpha[1],"_a1=",alpha1[1],"_a2=",alpha2[1],"_scale1=",para[case,7],"_scale2=",para[case,8])
  load(file = paste0(file.name,".rData"))
  total.est.result = rbind(total.est.result,colMeans(est.result))
}

compare.result.i1 = total.est.result[,c("continuum.pls.i1","AJIVE.i1")]
compare.result.i2 = total.est.result[,c("continuum.pls.i2","AJIVE.i2")]
total.result.i1 = total.est.result[,c("continuum.pls.i1","AJIVE.i1","continuum.pls.stand.i1")]
total.result.i2 = total.est.result[,c("continuum.pls.i2","AJIVE.i2","continuum.pls.stand.i2")]
total.est.result.i1.order = NULL
total.est.result.i2.order = NULL
compare.est.result.i1.order = NULL
compare.est.result.i2.order = NULL
for(case in seq(1,nrow(para))){
  total.est.result.i1.order = rbind(total.est.result.i1.order,colnames(total.result.i1)[order(total.result.i1[case,])])
  total.est.result.i2.order = rbind(total.est.result.i2.order,colnames(total.result.i2)[order(total.result.i2[case,])])
  compare.est.result.i1.order = rbind(compare.est.result.i1.order,colnames(compare.result.i1)[order(compare.result.i1[case,])])
  compare.est.result.i2.order = rbind(compare.est.result.i2.order,colnames(compare.result.i2)[order(compare.result.i2[case,])])
}
View(compare.est.result.i1.order)
View(compare.est.result.i2.order)












para = matrix(0,nrow=0,ncol=8)
para = rbind(para,c(80,50,50,1,1,1,1,1))
para = rbind(para,c(80,50,50,10,1,1,1,1))
para = rbind(para,c(80,50,50,1,10,1,1,1))
para = rbind(para,c(80,50,50,1,1,10,1,1))
para = rbind(para,c(80,50,50,1,1,1,10,1))
para = rbind(para,c(80,50,50,1,1,1,1,10))
para = rbind(para,c(80,100,10,1,1,1,1,1))
para = rbind(para,c(80,100,10,10,1,1,1,1))
para = rbind(para,c(80,100,10,100,1,1,1,1))
para = rbind(para,c(80,100,10,1,1,1,10,1))
para = rbind(para,c(80,100,10,10,1,1,10,1))
para = rbind(para,c(80,100,10,100,1,1,10,1))
para = rbind(para,c(80,100,10,1,1,1,1,10))
para = rbind(para,c(80,100,10,10,1,1,1,10))
para = rbind(para,c(80,100,10,100,1,1,1,10))
para = rbind(para,c(80,100,10,1,10,1,1,1))
para = rbind(para,c(80,100,10,1,100,1,1,1))
para = rbind(para,c(80,100,10,1,1,10,1,1))
para = rbind(para,c(80,100,10,1,1,100,1,1))
para = rbind(para,c(80,100,10,1,1,10,1,10))
para = rbind(para,c(80,100,10,1,1,100,1,10))
para = rbind(para,c(80,100,10,1,1,10,10,1))
para = rbind(para,c(80,100,10,1,1,100,10,1))
total.result = NULL
for(case in seq(1:nrow(para))){
  G = 2
  n = para[case,1]
  p1 = para[case,2]
  p2 = para[case,3]
  p=p1+p2
  r = 1
  r1 = 1
  r2 = 1
  r.list = list(r1, r2)
  L = 20
  alpha = rep(para[case,4], r)
  alpha1 = rep(para[case,5], r1) #OLS: 0
  alpha2 = rep(para[case,6], r2) #OLS: 0
  file.name = paste0("multi_source_result_pls_n=",n,"_p1=",p1,"_p2=",p2,"_a=",alpha[1],"_a1=",alpha1[1],"_a2=",alpha2[1],"_scale1=",para[case,7],"_scale2=",para[case,8],"_sigma")
  load(file = paste0(file.name,".rData"))
  total.result = rbind(total.result,colMeans(result))
}
compare.result = total.result[,c("continuum.pls","continuum.pls.stand","pls","pls.scaled")]
total.result.order = NULL
compare.result.order = NULL
for(case in seq(1,nrow(para))){
  total.result.order = rbind(total.result.order,colnames(total.result)[order(total.result[case,])])
  compare.result.order = rbind(compare.result.order,colnames(compare.result)[order(compare.result[case,])])
}
total.est.result = NULL
for(case in seq(1,nrow(para))){
  file.name = paste0("multi_source_est_result_pls_n=",para[case,1],"_p1=",para[case,2],"_p2=",para[case,3],"_a=",para[case,4],"_a1=",para[case,5],"_a2=",para[case,6],"_scale1=",para[case,7],"_scale2=",para[case,8],"_sigma")
  load(file = paste0(file.name,".rData"))
  total.est.result = rbind(total.est.result,colMeans(est.result))
}
compare.result.i1 = total.est.result[,c("continuum.pls.i1","AJIVE.i1")]
compare.result.i2 = total.est.result[,c("continuum.pls.i2","AJIVE.i2")]
total.result.i1 = total.est.result[,c("continuum.pls.i1","AJIVE.i1","continuum.pls.stand.i1")]
total.result.i2 = total.est.result[,c("continuum.pls.i2","AJIVE.i2","continuum.pls.stand.i2")]
total.est.result.i1.order = NULL
total.est.result.i2.order = NULL
compare.est.result.i1.order = NULL
compare.est.result.i2.order = NULL
for(case in seq(1,nrow(para))){
  total.est.result.i1.order = rbind(total.est.result.i1.order,colnames(total.result.i1)[order(total.result.i1[case,])])
  total.est.result.i2.order = rbind(total.est.result.i2.order,colnames(total.result.i2)[order(total.result.i2[case,])])
  compare.est.result.i1.order = rbind(compare.est.result.i1.order,colnames(compare.result.i1)[order(compare.result.i1[case,])])
  compare.est.result.i2.order = rbind(compare.est.result.i2.order,colnames(compare.result.i2)[order(compare.result.i2[case,])])
}
View(compare.est.result.i1.order)
View(compare.est.result.i2.order)

View(compare.result.order)



n.rep=30
para = matrix(0,nrow=0,ncol=9)
para = rbind(para,c(40,80,80,1,1,1,1,1,1))
para = rbind(para,c(40,80,80,100,1,1,1,1,1))
para = rbind(para,c(40,80,80,300,1,1,1,1,1))
para = rbind(para,c(40,80,80,1,100,1,1,1,1))
para = rbind(para,c(40,80,80,1,300,1,1,1,1))
para = rbind(para,c(40,80,80,1,1,100,1,1,1))
para = rbind(para,c(40,80,80,1,1,300,1,1,1))
para = rbind(para,c(40,80,80,1,1,1,100,1,1))
para = rbind(para,c(40,80,80,1,1,1,1,100,1))
para = rbind(para,c(40,80,80,1,1,1,1,1,100))
para = rbind(para,c(40,200,40,1,1,1,1,1,1))
para = rbind(para,c(40,200,40,1,1,1,1,100,1))
para = rbind(para,c(40,200,40,1,1,100,1,1,100))
para = rbind(para,c(40,200,40,1,1,300,1,1,300))
para = rbind(para,c(40,200,40,1,1,100,1,1,1))
para = rbind(para,c(40,200,40,1,1,300,1,1,1))
para = rbind(para,c(40,200,40,1,1,100,1,50,1))
para = rbind(para,c(40,200,40,1,1,300,1,50,1))
para = rbind(para,c(40,200,40,1,1,1,1,1,100))
para = rbind(para,c(40,200,40,1,1,1,1,1,300))
para = rbind(para,c(40,200,40,100,1,1,1,1,1))
total.result = NULL
for(case in seq(1,nrow(para))){
  G = 2
  n = para[case,1]
  p1 = para[case,2]
  p2 = para[case,3]
  p=p1+p2
  r = 1
  r1 = 1
  r2 = 1
  r.list = list(r1, r2)
  L = 20
  alpha = rep(para[case,4], r)
  alpha1 = rep(para[case,5], r1) #OLS: 0
  alpha2 = rep(para[case,6], r2) #OLS: 0
  file.name = paste0("multi_source_result_pls_n=",n,"_p1=",p1,"_p2=",p2,"_a=",alpha[1],"_a1=",alpha1[1],"_a2=",alpha2[1],"_scale1=",para[case,7],"_scale2=",para[case,8],"_scalej=",para[case,9],"_sigma3")
  load(file = paste0(file.name,".rData"))
  total.result = rbind(total.result,colMeans(result))
}
compare.result = total.result[,c("continuum.pls","continuum.pls.stand","pls","pls.scaled")]
total.result.order = NULL
compare.result.order = NULL
for(case in seq(1,nrow(para))){
  total.result.order = rbind(total.result.order,colnames(total.result)[order(total.result[case,])])
  compare.result.order = rbind(compare.result.order,colnames(compare.result)[order(compare.result[case,])])
}
total.est.result = NULL
for(case in seq(1,nrow(para))){
  file.name = paste0("multi_source_est_result_pls_n=",para[case,1],"_p1=",para[case,2],"_p2=",para[case,3],"_a=",para[case,4],"_a1=",para[case,5],"_a2=",para[case,6],"_scale1=",para[case,7],"_scale2=",para[case,8],"_scalej=",para[case,9],"_sigma3")
  load(file = paste0(file.name,".rData"))
  total.est.result = rbind(total.est.result,colMeans(est.result))
}
compare.result.i1 = total.est.result[,c("continuum.pls.i1","AJIVE.i1","continuum.pcr.i1")]
compare.result.i2 = total.est.result[,c("continuum.pls.i2","AJIVE.i2","continuum.pcr.i2")]
compare.result.j = total.est.result[,c("continuum.pls.j","AJIVE.J","continuum.pcr.j")]
total.result.i1 = total.est.result[,c("continuum.pls.i1","AJIVE.i1","continuum.pls.stand.i1","continuum.pcr.i1","continuum.pcr.stand.i1")]
total.result.i2 = total.est.result[,c("continuum.pls.i2","AJIVE.i2","continuum.pls.stand.i2","continuum.pcr.i2","continuum.pcr.stand.i2")]
total.result.j = total.est.result[,c("continuum.pls.j","AJIVE.J","continuum.pls.stand.j","continuum.pcr.j","continuum.pcr.stand.j")]

total.est.result.i1.order = NULL
total.est.result.i2.order = NULL
total.est.result.j.order = NULL
compare.est.result.i1.order = NULL
compare.est.result.i2.order = NULL
compare.est.result.j.order = NULL
total.est.result.j.order = NULL
compare.est.result.i1.order = NULL
for(case in seq(1,nrow(para))){
  total.est.result.i1.order = rbind(total.est.result.i1.order,colnames(total.result.i1)[order(total.result.i1[case,])])
  total.est.result.i2.order = rbind(total.est.result.i2.order,colnames(total.result.i2)[order(total.result.i2[case,])])
  compare.est.result.i1.order = rbind(compare.est.result.i1.order,colnames(compare.result.i1)[order(compare.result.i1[case,])])
  total.est.result.j.order = rbind(total.est.result.j.order,colnames(total.result.j)[order(total.result.j[case,])])
  compare.est.result.j.order = rbind(compare.est.result.j.order,colnames(compare.result.j)[order(compare.result.j[case,])])
  compare.est.result.i2.order = rbind(compare.est.result.i2.order,colnames(compare.result.i2)[order(compare.result.i2[case,])])
}
View(total.est.result.i1.order)
View(total.est.result.i2.order)
View(total.est.result.j.order)
View(compare.result.order)
View(total.result.order)