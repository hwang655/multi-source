# ---------------------- reading shell command --------------------- 
# args = (commandArgs(TRUE))
# cat(args, "\n")
# for (i in 1:length(args)) {
#   eval(parse(text = args[[i]]))
# }
# ------------------------------------------------------------------ 

library(MASS)
library(rlist)
library(glmnet)
library(methods)
library(pls)

current = getwd()
# setwd("/nas/longleaf/home/peiyao/continuum/")
source("./function/jive_continuum.R")
# setwd(current)
para = matrix(0,nrow=10,ncol=6)
para[1,] = c(100,50,50,1,1,4)
para[2,] = c(100,90,50,1,1,4)
para[3,] = c(100,90,50,1,4,1)
para[4,] = c(100,50,50,1,1,1)
para[5,] = c(100,90,50,1,4,1)
para[6,] = c(100,90,50,1,1,4)
seed.seq = 1:20
for(case in 1:6){
  result = NULL
  for(myseed in seed.seq){
    print(myseed)
    set.seed(myseed)
    
    MSE = NULL
    
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
    
    X1 = mvrnorm(n, rep(0, p1), diag(p1))
    X2 = mvrnorm(n, rep(0, p2), diag(p2))
    X = cbind(X1, X2)
    
    q = min(n, p)/2
    q1 = min(n, p1)/2
    q2 = min(n, p2)/2
    V = svd(X)$v[,1:q]%*%rep(1/sqrt(q), q)
    V1 = svd(X1- X%*%V%*%t(V[1:p1,]))$v[,1:q1]%*%rep(1/sqrt(q1), q1)
    V2 = svd(X2- X%*%V%*%t(V[(p1+1):p,]))$v[,1:q2]%*%rep(1/sqrt(q2), q2)
    # V = svd(X)$v[,1:q]
    # V1 = svd(X1- X%*%V%*%t(V[1:p1,]))$v[,1:q1]
    # V2 = svd(X2- X%*%V%*%t(V[(p1+1):p,]))$v[,1:q2]
    
    e = rnorm(n)*.2
    Y = X%*%V%*%alpha + X1%*%V1%*%alpha1 + X2%*%V2%*%alpha2 + e
    
    
    sum(e^2)/sum(Y^2)
    
    X.list = list(t(X1), t(X2))
    
    X1 = mvrnorm(n, rep(0, p1), diag(p1))
    X2 = mvrnorm(n, rep(0, p2), diag(p2))
    
    e = rnorm(n)*.2
    Y.test = X%*%V%*%alpha + X1%*%V1%*%alpha1 + X2%*%V2%*%alpha2 + e
    
    X.test = cbind(X1, X2)
    
    X.test.list = list(t(X1), t(X2))
    {
      a = seq(0, 1, length.out = L+1)
      gam.list = a/(1-a)
      gam.list[L+1] = 1e10
      
      rr = matrix(0,nrow=6,ncol=3)
      rr[1,] = c(1,0,0)
      rr[2,] = c(2,0,0)
      rr[3,] = c(1,1,1)
      rr[4,] = c(2,2,2)
      rr[5,] = c(0,1,1)
      rr[6,] = c(0,2,2)
      
      MSE = list()
      for (ind in 1:nrow(rr)){
        MSE0 = NULL
        for (gam in gam.list){
          r0=rr[ind,1]
          r11=rr[ind,2]
          r22=rr[ind,3]
          ml = continuum.multisource.iter.v1(X.list, Y, lambda = 0, gam = gam, rankJ = r0, rankA = c(r11, r22),
                                                                               center.X = F, scale.X = F, center.Y = F, scale.Y = F, orthIndiv = T)
          beta.Cind = ml$beta.Cind
          Yhat.heter = (X1)%*%beta.Cind[[1]]+(X2)%*%beta.Cind[[2]]
          MSE0= c(MSE0,mean((Y.test - as.numeric(ml$intercept) - X.test%*%ml$beta.C - Yhat.heter)^2))
        }
        MSE[[paste0("MSE",r0,r11,r22)]] = MSE0
      }

    }
    
    # ml.continuum.pcr = continuum.multisource.iter.v1(X.list, Y, lambda = 0, gam = 1e10, rankJ = r, rankA = c(r1, r2))
    # ml.continuum = ml.continuum.pcr
    # # ml = decomposeX(X.test.list, ml.continuum$U[[ml.continuum$nrun]], ml.continuum$W[[ml.continuum$nrun]], ml.continuum$centerValues.X, ml.continuum$scaleValues.X)
    # beta.Cind = ml.continuum$beta.Cind
    # Yhat.heter = (X1)%*%beta.Cind[[1]]+(X2)%*%beta.Cind[[2]]
    # 
    # MSE= c(MSE,mean((Y.test - as.numeric(ml.continuum$intercept) - X.test%*%ml.continuum$beta.C - Yhat.heter)^2))
    # 
    # ml.continuum.pls = continuum.multisource.iter.v1(X.list, Y, lambda = 0, gam = 1, rankJ = r, rankA = c(r1, r2))
    # ml.continuum = ml.continuum.pls
    # beta.Cind = ml.continuum$beta.Cind
    # Yhat.heter = (X1)%*%beta.Cind[[1]]+(X2)%*%beta.Cind[[2]]
    # 
    # MSE= c(MSE,mean((Y.test - as.numeric(ml.continuum$intercept) - X.test%*%ml.continuum$beta.C - Yhat.heter)^2))


    # ml.continuum.pcr = continuum.multisource.iter.v1(X.list, Y, lambda = 0, gam = 1e10, rankJ = r, rankA = c(r1, r2),center.X = F, scale.X = F, center.Y = F, scale.Y = F,  orthIndiv = T)
    # ml.continuum = ml.continuum.pcr
    # # ml = decomposeX(X.test.list, ml.continuum$U[[ml.continuum$nrun]], ml.continuum$W[[ml.continuum$nrun]], ml.continuum$centerValues.X, ml.continuum$scaleValues.X)
    # beta.Cind = ml.continuum$beta.Cind
    # Yhat.heter = (X1)%*%beta.Cind[[1]]+(X2)%*%beta.Cind[[2]]
    # 
    # MSE= c(MSE,mean((Y.test - as.numeric(ml.continuum$intercept) - X.test%*%ml.continuum$beta.C - Yhat.heter)^2))
    # 
    # ml.continuum.pls = continuum.multisource.iter.v1(X.list, Y, lambda = 0, gam = 1, rankJ = r, rankA = c(r1, r2),center.X = F, scale.X = F, center.Y = F, scale.Y = F,  orthIndiv = T)
    # ml.continuum = ml.continuum.pls
    # beta.Cind = ml.continuum$beta.Cind
    # Yhat.heter = (X1)%*%beta.Cind[[1]]+(X2)%*%beta.Cind[[2]]
    # 
    # MSE= c(MSE,mean((Y.test - as.numeric(ml.continuum$intercept) - X.test%*%ml.continuum$beta.C - Yhat.heter)^2))
    # ml.continuum.pls = continuum.multisource.iter.v1(X.list, Y, lambda = 0, gam = 1, rankJ = r, rankA = c(r1, r2),center.X = F, scale.X = F, center.Y = F, scale.Y = F,  orthIndiv = F)
    # ml.continuum = ml.continuum.pls
    # beta.Cind = ml.continuum$beta.Cind
    # Yhat.heter = (X1)%*%beta.Cind[[1]]+(X2)%*%beta.Cind[[2]]
    # 
    # MSE= c(MSE,mean((Y.test - as.numeric(ml.continuum$intercept) - X.test%*%ml.continuum$beta.C - Yhat.heter)^2))
    
    
    
    ml.ridge = cv.glmnet(x = X, y = Y, alpha = 0, standardize = F, intercept = F)
    MSE[["ridge"]]= mean((Y.test -predict(ml.ridge, newx = X.test, s = ml.ridge$lambda.min) )^2)
    
    ml.pls = plsr(Y~X, validation = "CV", center = F, scale = F)
    ncomp.pls = selectNcomp(ml.pls, method = "randomization", plot = F)
    MSE[["pls"]]= mean((predict(ml.pls, newdata = X.test, ncomp = r+r1+r2)[,,1] - Y.test)^2)
    
    ml.pcr = pcr(Y~X, validation = "CV", center = F, scale = F)
    ncomp.pcr = selectNcomp(ml.pcr, method = "randomization", plot = F)
    MSE[["pcr"]]= mean((predict(ml.pcr, newdata = X.test, ncomp = r+r1+r2)[,,1] - Y.test)^2)
    if (myseed==1){
      result = MSE
    }
    else{
      result = mapply("+", result, MSE, SIMPLIFY = FALSE)
    }
    # print(MSE)
    # result = rbind(result,MSE)
  }
  # colnames(result) = c("cr.pcr","cr.pls","cr.pls.n.orth","ridge","pls","pcr")
  result = mapply("/", result,length(seed.seq) , SIMPLIFY = FALSE)
  file.name = paste0("multi_source_result_pls_n",n,"_p1",p1,"_p2",p2,"_a",alpha[1],"_a1",alpha1[1],"_a2",alpha2[1])
  # write.table(result, file = file.name, sep = ',', append = T, col.names = F, row.names = F)
  save(result,file = paste0(file.name,".rData"))
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
