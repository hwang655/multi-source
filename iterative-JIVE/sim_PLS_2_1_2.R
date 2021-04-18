# ---------------------- reading shell command --------------------- 
# args = (commandArgs(TRUE))
# cat(args, "\n")
# for (i in 1:length(args)) {
#   eval(parse(text = args[[i]]))
# }
# ------------------------------------------------------------------ 
case.seq=6*3+1:3
noise=1
library(MASS)
library(rlist)
library(glmnet)
library(methods)
library(pls)
library(ajive)
current = getwd()
# setwd("/nas/longleaf/home/peiyao/continuum/")
setwd("D:/git-project/multi-source/iterative-JIVE")
source("./function/jive_continuum.R")
# setwd(current)
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
seed.seq = 1:n.rep
for(case in case.seq){
  
  result = NULL
  est.result = NULL
  for(myseed in seed.seq){
    # print(myseed*10)
    set.seed(myseed*10)
    
    MSE = NULL
    est.mse = NULL
    
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
    
    # Q = randortho(p)
    # V = matrix(Q[,1:r], ncol = r)
    # V1 = matrix(Q[,r+(1:r1)], ncol = r1)
    # V2 = matrix(Q[,r+r1+(1:r2)], ncol = r2)
    
    Siga = diag(p)
    for(i in 1:p)
      for(j in 1:p)
        Siga[i,j] = (0.6)^{abs(i-j)}
    X = mvrnorm(n,rep(0,p),Siga)
    X1 = X[,1:p1]
    X2 = X[,(p1+1):p]
    
    # X1 = mvrnorm(n, rep(0, p1), diag(p1))
    # X2 = mvrnorm(n, rep(0, p2), diag(p2))
    
    q = min(n, p)/2
    q1 = min(n, p1)/2
    q2 = min(n, p2)/2
    V = svd(X)$v[,1:q]%*%rep(1/sqrt(q), q)
    V1 = svd(X1- X%*%V%*%t(V[1:p1,]))$v[,1:q1]%*%rep(1/sqrt(q1), q1)
    V2 = svd(X2- X%*%V%*%t(V[(p1+1):p,]))$v[,1:q2]%*%rep(1/sqrt(q2), q2)
    X1 = (X1-X1%*%V1%*%t(V1)-X%*%V%*%t(matrix(V[1:p1,],ncol=1)))*noise+X1%*%V1%*%t(V1)*para[case,7]+X%*%V%*%t(matrix(V[1:p1,],ncol=1))*para[case,9]
    V1 = V1*para[case,7]
    X2 = (X2-X2%*%V2%*%t(V2)-X%*%V%*%t(matrix(V[(p1+1):p,],ncol=1)))*noise+X2%*%V2%*%t(V2)*para[case,8]+X%*%V%*%t(matrix(V[(p1+1):p,],ncol=1))*para[case,9]
    V2 = V2*para[case,8]
    X = cbind(X1, X2)
    e = rnorm(n)*.2
    Y = X%*%V%*%alpha + X1%*%V1%*%alpha1 + X2%*%V2%*%alpha2 + e
    X.list = list(t(X1), t(X2))
    I1 = X1%*%V1%*%t(V1)
    I2 = X2%*%V2%*%t(V2)
    J = X%*%V%*%t(V)
    ajive.result = ajive::ajive(list(X1,X2),initial_signal_ranks = c(r1,r2),joint_rank = r)
    est.mse[["AJIVE.i1"]]= norm(ajive.result$block_decomps[[1]]$individual$full-X1%*%V1%*%t(V1),"2")
    est.mse[["AJIVE.i2"]]= norm(ajive.result$block_decomps[[2]]$individual$full-X2%*%V2%*%t(V2),"2")
    est.mse[["AJIVE.J"]]= norm(cbind(ajive.result$block_decomps[[1]]$joint$full,ajive.result$block_decomps[[2]]$joint$full)-J,"2")
    
    temp = cbind(ajive.result$block_decomps[[1]]$joint$full,ajive.result$block_decomps[[2]]$joint$full)
    ajive.x = svd(temp)$u[,1:r]%*%diag(svd(temp)$d[1:r],r)
    temp = ajive.result$block_decomps[[1]]$individual$full
    ajive.x = cbind(ajive.x,svd(temp)$u[,1:r1]%*%diag(svd(temp)$d[1:r1],r1))
    temp = ajive.result$block_decomps[[2]]$individual$full
    ajive.x = cbind(ajive.x,svd(temp)$u[,1:r2]%*%diag(svd(temp)$d[1:r2],r2))
    ajive.beta = ginv(t(ajive.x)%*%ajive.x)%*%t(ajive.x)%*%Y
    
    
    X = mvrnorm(n,rep(0,p),Siga)
    X1 = X[,1:p1]
    X2 = X[,(p1+1):p]
    X1 = (X1-X1%*%V1%*%t(V1)/para[case,7]-X%*%V%*%t(matrix(V[1:p1,],ncol=1)/para[case,9]))*noise+X1%*%V1%*%t(V1)+X%*%V%*%t(matrix(V[1:p1,],ncol=1))
    X2 = (X2-X2%*%V2%*%t(V2)/para[case,8]-X%*%V%*%t(matrix(V[(p1+1):p,],ncol=1)/para[case,9]))*noise+X2%*%V2%*%t(V2)+X%*%V%*%t(matrix(V[(p1+1):p,],ncol=1))
    
    e = rnorm(n)*.2
    X.test = cbind(X1, X2)
    Y.test = X.test%*%V%*%alpha + X1%*%V1%*%alpha1 + X2%*%V2%*%alpha2 + e
    
    X.test = cbind(X1, X2)
    
    X.test.list = list(t(X1), t(X2))
    ajive.result = ajive::ajive(list(X1,X2),initial_signal_ranks = c(r1,r2),joint_rank = r)
    temp = cbind(ajive.result$block_decomps[[1]]$joint$full,ajive.result$block_decomps[[2]]$joint$full)
    ajive.x = svd(temp)$u[,1:r]%*%diag(svd(temp)$d[1:r],r)
    temp = ajive.result$block_decomps[[1]]$individual$full
    ajive.x = cbind(ajive.x,svd(temp)$u[,1:r1]%*%diag(svd(temp)$d[1:r1],r1))
    temp = ajive.result$block_decomps[[2]]$individual$full
    ajive.x = cbind(ajive.x,svd(temp)$u[,1:r2]%*%diag(svd(temp)$d[1:r2],r2))
    MSE["AJIVE"] = mean((Y.test - ajive.x%*%ajive.beta)^2)
    
    X.list.scaled = list(X.list[[1]]/norm(X.list[[1]],"F"),X.list[[2]]/norm(X.list[[2]],"F"))
    X.scaled = cbind(t(X.list.scaled[[1]]),t(X.list.scaled[[2]]))
    X.test.list.scaled = list(X.test.list[[1]]/norm(X.test.list[[1]],"F"),X.test.list[[2]]/norm(X.test.list[[2]],"F"))
    X.test.scaled = cbind(t(X.test.list.scaled[[1]]),t(X.test.list.scaled[[2]]))
    X = cbind(t(X.list[[1]]),t(X.list[[2]]))
    
    
    ml.continuum.pcr = continuum.multisource.iter.v1(X.list, Y, lambda = 0, gam = 1e10, rankJ = r, rankA = c(r1, r2),center.X = F, scale.X = F, center.Y = F, scale.Y = F,  orthIndiv = F)
    ml.continuum = ml.continuum.pcr
    # ml = decomposeX(X.test.list, ml.continuum$U[[ml.continuum$nrun]], ml.continuum$W[[ml.continuum$nrun]], ml.continuum$centerValues.X, ml.continuum$scaleValues.X)
    beta.Cind = ml.continuum$beta.Cind
    Yhat.heter = (X1)%*%beta.Cind[[1]]+(X2)%*%beta.Cind[[2]]
    # 
    MSE[["continuum.pcr"]]= mean((Y.test - as.numeric(ml.continuum$intercept) - X.test%*%ml.continuum$beta.C - Yhat.heter)^2)
    est.mse[["continuum.pcr.i1"]]= norm(t(ml.continuum$I[[1]])-I1,"2")
    est.mse[["continuum.pcr.i2"]]= norm(t(ml.continuum$I[[2]])-I2,"2")
    est.mse[["continuum.pcr.j"]]= norm(t(ml.continuum$J)-(J),"2")
    
    ml.continuum.pls = continuum.multisource.iter.v1(X.list, Y, maxiter=1000,lambda = 0, gam = 1, rankJ = r, rankA = c(r1, r2),center.X = F, scale.X = F, center.Y = F, scale.Y = F,  orthIndiv = F)
    ml.continuum = ml.continuum.pls
    # ml = decomposeX(X.test.list, ml.continuum$U[[ml.continuum$nrun]], ml.continuum$W[[ml.continuum$nrun]], ml.continuum$centerValues.X, ml.continuum$scaleValues.X)
    beta.Cind = ml.continuum$beta.Cind
    Yhat.heter = (X1)%*%beta.Cind[[1]]+(X2)%*%beta.Cind[[2]]
    # 
    MSE[["continuum.pls"]]= mean((Y.test - as.numeric(ml.continuum$intercept) - X.test%*%ml.continuum$beta.C - Yhat.heter)^2)
    est.mse[["continuum.pls.i1"]]= norm(t(ml.continuum$I[[1]])-I1,"2")
    est.mse[["continuum.pls.i2"]]= norm(t(ml.continuum$I[[2]])-I2,"2")
    est.mse[["continuum.pls.j"]]= norm(t(ml.continuum$J)-(J),"2")
    
    ml.continuum.pcr = continuum.multisource.iter.v1(X.list, Y, lambda = 0, gam = 1e10, rankJ = r, rankA = c(r1, r2),center.X = T, scale.X = T, center.Y = T, scale.Y = T,  orthIndiv = F)
    ml.continuum = ml.continuum.pcr
    # ml = decomposeX(X.test.list, ml.continuum$U[[ml.continuum$nrun]], ml.continuum$W[[ml.continuum$nrun]], ml.continuum$centerValues.X, ml.continuum$scaleValues.X)
    beta.Cind = ml.continuum$beta.Cind
    Yhat.heter = (X1)%*%beta.Cind[[1]]+(X2)%*%beta.Cind[[2]]
    # 
    MSE[["continuum.pcr.stand"]]= mean((Y.test - as.numeric(ml.continuum$intercept) - X.test%*%ml.continuum$beta.C - Yhat.heter)^2)
    est.mse[["continuum.pcr.stand.i1"]]= norm(t(ml.continuum$I[[1]])-I1,"2")
    est.mse[["continuum.pcr.stand.i2"]]= norm(t(ml.continuum$I[[2]])-I2,"2")
    est.mse[["continuum.pcr.stand.j"]]= norm(t(ml.continuum$J)-(J),"2")
    
    ml.continuum.pls = continuum.multisource.iter.v1(X.list, Y, lambda = 0, gam = 1, rankJ = r, rankA = c(r1, r2),center.X = T, scale.X = T, center.Y = T, scale.Y = T,  orthIndiv = F)
    ml.continuum = ml.continuum.pls
    # ml = decomposeX(X.test.list, ml.continuum$U[[ml.continuum$nrun]], ml.continuum$W[[ml.continuum$nrun]], ml.continuum$centerValues.X, ml.continuum$scaleValues.X)
    beta.Cind = ml.continuum$beta.Cind
    Yhat.heter = (X1)%*%beta.Cind[[1]]+(X2)%*%beta.Cind[[2]]
    # 
    MSE[["continuum.pls.stand"]]= mean((Y.test - as.numeric(ml.continuum$intercept) - X.test%*%ml.continuum$beta.C - Yhat.heter)^2)
    est.mse[["continuum.pls.stand.i1"]]= norm(t(ml.continuum$I[[1]])-I1,"2")
    est.mse[["continuum.pls.stand.i2"]]= norm(t(ml.continuum$I[[2]])-I2,"2")
    est.mse[["continuum.pls.stand.j"]]= norm(t(ml.continuum$J)-(J),"2")
    
    
    
    
    
    ml.ridge = cv.glmnet(x = X, y = Y, lambda = 10^seq(3,-10,length=80),alpha = 0, standardize = T, intercept = T)
    MSE[["ridge"]]= mean((Y.test -predict(ml.ridge, newx = X.test, s = ml.ridge$lambda.min) )^2)
    
    ml.pls = plsr(Y~X, validation = "CV", center = T, scale = T)
    ncomp.pls = selectNcomp(ml.pls, method = "randomization", plot = F)
    MSE[["pls"]]= mean((predict(ml.pls, newdata = X.test, ncomp = r+r1+r2)[,,1] - Y.test)^2)
    
    ml.pcr = pcr(Y~X, validation = "CV", center = T, scale = T)
    ncomp.pcr = selectNcomp(ml.pcr, method = "randomization", plot = F)
    MSE[["pcr"]]= mean((predict(ml.pcr, newdata = X.test, ncomp = r+r1+r2)[,,1] - Y.test)^2)
    

    ml.ridge = cv.glmnet(x = X.scaled, y = Y, alpha = 0, standardize = T, intercept = T)
    MSE[["ridge.scaled"]]= mean((Y.test -predict(ml.ridge, newx = X.test.scaled, s = ml.ridge$lambda.min) )^2)

    ml.pls = plsr(Y~X.scaled, validation = "CV", center = T, scale = T)
    ncomp.pls = selectNcomp(ml.pls, method = "randomization", plot = F)
    MSE[["pls.scaled"]]= mean((predict(ml.pls, newdata = X.test.scaled, ncomp = r+r1+r2)[,,1] - Y.test)^2)

    ml.pcr = pcr(Y~X.scaled, validation = "CV", center = T, scale = T)
    ncomp.pcr = selectNcomp(ml.pcr, method = "randomization", plot = F)
    MSE[["pcr.scaled"]]= mean((predict(ml.pcr, newdata = X.test.scaled, ncomp = r+r1+r2)[,,1] - Y.test)^2)

    # MSE[["noise"]] = mean(e^2)
    
    if (myseed==1){
      result = MSE
      est.result = est.mse
    }
    else{
      result = rbind(result,MSE)
      est.result = rbind(est.result,est.mse)
      #   result = mapply("+", result, MSE, SIMPLIFY = FALSE)
    }
    # print(MSE)
  }
  # colnames(result) = c("cr.pcr","cr.pls","cr.pls.n.orth","ridge","pls","pcr")
  # result = mapply("/", result,length(seed.seq) , SIMPLIFY = FALSE)
  file.name = paste0("multi_source_result_pls_n=",n,"_p1=",p1,"_p2=",p2,"_a=",alpha[1],"_a1=",alpha1[1],"_a2=",alpha2[1],"_scale1=",para[case,7],"_scale2=",para[case,8],"_scalej=",para[case,9],"_sigma3")
  # write.table(result, file = file.name, sep = ',', append = T, col.names = F, row.names = F)
  save(result,file = paste0("2_1_noise1/",file.name,".rData"))
  file.name = paste0("multi_source_est_result_pls_n=",n,"_p1=",p1,"_p2=",p2,"_a=",alpha[1],"_a1=",alpha1[1],"_a2=",alpha2[1],"_scale1=",para[case,7],"_scale2=",para[case,8],"_scalej=",para[case,9],"_sigma3")
  # print(c(case,colMeans(result)))
  save(est.result,file = paste0("2_1_noise1/",file.name,".rData"))
}




# a = seq(0, 1, length.out = L+1)
# gam.list = a/(1-a)
# gam.list[L+1] = 1e10
# 
# ml.100.list = list()
# for (gam in gam.list){
#   # print(gam)
#   ml.100.list = list.append(ml.100.list, continuum.multisource.iter.v1(X.list, Y, lambda = 0, maxiter = 300, gam = gam, rankJ = 1, rankA = c(0, 0)))
# }
# 
# ml.200.list = list()
# for (gam in gam.list){
#   # print(gam)
#   ml.200.list = list.append(ml.200.list, continuum.multisource.iter.v1(X.list, Y, lambda = 0, maxiter = 300, gam = gam, rankJ = 2, rankA = c(0, 0)))
# }
# 
# ml.111.list = list()
# for (gam in gam.list){
#   ml.111.list = list.append(ml.111.list, continuum.multisource.iter.v1(X.list, Y, lambda = 0, maxiter = 300, gam = gam, rankJ = 1, rankA = c(1, 1)))
# }
# 
# ml.011.list = list()
# for (gam in gam.list){
#   # print(gam)
#   ml.011.list = list.append(ml.011.list, continuum.multisource.iter.v1(X.list, Y, lambda = 0, maxiter = 300, gam = gam, rankJ = 0, rankA = c(1, 1)))
# }
# 
# ml.022.list = list()
# for (gam in gam.list){
#   ml.022.list = list.append(ml.022.list, continuum.multisource.iter.v1(X.list, Y, lambda = 0, maxiter = 300, gam = gam, rankJ = 0, rankA = c(2, 2)))
# }
# 
# ml.222.list = list()
# for (gam in gam.list){
#   ml.222.list = list.append(ml.222.list, continuum.multisource.iter.v1(X.list, Y, lambda = 0, maxiter = 300, gam = gam, rankJ = 2, rankA = c(2, 2)))
# }
# ml.ridge = cv.glmnet(x = X, y = Y, alpha = 0, standardize = F, intercept = F)
# 
# ml.pls = plsr(Y~X, validation = "CV", center = F, scale = F)
# ncomp.pls = selectNcomp(ml.pls, method = "randomization", plot = F)
# 
# ml.pcr = pcr(Y~X, validation = "CV", center = F, scale = F)
# ncomp.pcr = selectNcomp(ml.pcr, method = "randomization", plot = F)
# 
# ml.ridge.list = lapply(1:G, function(g) cv.glmnet(x = X.list[[g]], y = Y.list[[g]], alpha = 0,
#                                                   standardize = F, intercept = F))
# 
# ml.pls.list = lapply(1:G, function(g) plsr(Y.list[[g]] ~ X.list[[g]], validation = "CV", center = F, scale = F))
# ncomp.pls.list = lapply(1:G, function(g) selectNcomp(ml.pls.list[[g]], method = "randomization", plot = F))
# 
# ml.pcr.list = lapply(1:G, function(g) pcr(Y.list[[g]]~X.list[[g]], validation = "CV", center = F, scale = F))
# ncomp.pcr.list = lapply(1:G, function(g) selectNcomp(ml.pcr.list[[g]], method = "randomization", plot = F))


# MSE = list()
# for (ml in ml.100.list){
#   beta.Cind = ml$beta.Cind
#   Yhat.heter = (X1)%*%beta.Cind[[1]]+(X2)%*%beta.Cind[[2]]
#   MSE = list.append(MSE, mean((Y.test - as.numeric(ml$intercept) - X.test%*%ml$beta.C - Yhat.heter)^2))
# }
# 
# for (ml in ml.200.list){
#   beta.Cind = ml$beta.Cind
#   Yhat.heter = (X1)%*%beta.Cind[[1]]+(X2)%*%beta.Cind[[2]]
#   MSE = list.append(MSE, mean((Y.test - as.numeric(ml$intercept) - X.test%*%ml$beta.C - Yhat.heter)^2))
# }
# 
# for (ml in ml.111.list){
#   beta.Cind = ml$beta.Cind
#   Yhat.heter = (X1)%*%beta.Cind[[1]]+(X2)%*%beta.Cind[[2]]
#   MSE = list.append(MSE, mean((Y.test - as.numeric(ml$intercept) - X.test%*%ml$beta.C - Yhat.heter)^2))
# }
# 
# for (ml in ml.011.list){
#   beta.Cind = ml$beta.Cind
#   Yhat.heter = (X1)%*%beta.Cind[[1]]+(X2)%*%beta.Cind[[2]]
#   MSE = list.append(MSE, mean((Y.test - as.numeric(ml$intercept) - X.test%*%ml$beta.C - Yhat.heter)^2))
# }
# 
# for (ml in ml.022.list){
#   beta.Cind = ml$beta.Cind
#   Yhat.heter = (X1)%*%beta.Cind[[1]]+(X2)%*%beta.Cind[[2]]
#   MSE = list.append(MSE, mean((Y.test - as.numeric(ml$intercept) - X.test%*%ml$beta.C - Yhat.heter)^2))
# }
# 
# for (ml in ml.222.list){
#   beta.Cind = ml$beta.Cind
#   Yhat.heter = (X1)%*%beta.Cind[[1]]+(X2)%*%beta.Cind[[2]]
#   MSE = list.append(MSE, mean((Y.test - as.numeric(ml$intercept) - X.test%*%ml$beta.C - Yhat.heter)^2))
# }
# 
# MSE = do.call(c, MSE)
# 
# ml.ridge = cv.glmnet(x = X, y = Y, alpha = 0, standardize = F, intercept = F)
# MSE= c(MSE, mean((Y.test -predict(ml.ridge, newx = X.test, s = ml.ridge$lambda.min) )^2))
# 
# ml.pls = plsr(Y~X, validation = "CV", center = F, scale = F)
# ncomp.pls = selectNcomp(ml.pls, method = "randomization", plot = F)
# MSE= c(MSE,mean((predict(ml.pls, newdata = X.test, ncomp = r+r1+r2)[,,1] - Y.test)^2))
# 
# ml.pcr = pcr(Y~X, validation = "CV", center = F, scale = F)
# ncomp.pcr = selectNcomp(ml.pcr, method = "randomization", plot = F)
# MSE= c(MSE,mean((predict(ml.pcr, newdata = X.test, ncomp = r+r1+r2)[,,1] - Y.test)^2))
# print(MSE)
# result = rbind(result,MSE)
# }
# # group models
# ml = ml.ridge.list
# MSE = list.append(MSE, sapply(1:G, function(g)
#   mean((Y.test.list[[g]] - predict(ml[[g]], newx = X.test.list[[g]], s = ml[[g]]$lambda.min))^2)))
# 
# ml = ml.pls.list
# MSE = list.append(MSE, sapply(1:G, function(g) mean((Y.test.list[[g]] - predict(ml[[g]], newdata = X.test.list[[g]], ncomp = r + r.list[[g]])[,,1]
# )^2)))
# 
# ml = ml.pcr.list
# MSE = list.append(MSE, sapply(1:G, function(g) mean((Y.test.list[[g]] - predict(ml[[g]], newdata = X.test.list[[g]], ncomp = r + r.list[[g]])[,,1]
# )^2)))
# 
# row.names(MSE) = c("iter.OLS", "iter.PLS", "iter.PCR", "global.ridge", "global.PLS",
#                    "global.PCR", "group.ridge", "group.PLS",
#                    "group.PCR")
# file.name = "result.csv"
# write.table(t(c(myseed, as.vector(MSE))), file = file.name, sep = ',', append = T, col.names = F, row.names = F)
# 
