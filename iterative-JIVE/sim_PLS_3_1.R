# ---------------------- reading shell command --------------------- 
# args = (commandArgs(TRUE))
# cat(args, "\n")
# for (i in 1:length(args)) {
#   eval(parse(text = args[[i]]))
# }
# ------------------------------------------------------------------ 
noise=1
para.index=6
library(MASS)
library(rlist)
library(glmnet)
library(methods)
library(pls)
library(ajive)
library(pracma)
current = getwd()
# setwd("/nas/longleaf/home/peiyao/continuum/")
setwd("D:/git-project/multi-source/iterative-JIVE")
source("./function/jive_continuum.R")
# setwd(current)
n.rep=10
para = matrix(0,nrow=0,ncol=9)
para = rbind(para,c(30,50,50,1,0.5,1,0.01,0.1,0.1))
para = rbind(para,c(30,50,50,1,1,1,100,100,100))
para = rbind(para,c(30,50,50,1,5,1,0.5,1,1))
para = rbind(para,c(30,50,50,1,25,1,0.5,1,1))
para = rbind(para,c(30,50,50,1,0.5,1,1,1,1))
para = rbind(para,c(30,50,50,1,1,1,1,1,1))
para = rbind(para,c(30,50,50,1,5,1,1,1,1))
para = rbind(para,c(30,50,50,1,25,1,1,1,1))
para = rbind(para,c(30,50,50,1,0.5,1,5,1,1))
para = rbind(para,c(30,50,50,1,1,1,5,1,1))
para = rbind(para,c(30,50,50,1,5,1,5,1,1))
para = rbind(para,c(30,50,50,1,25,1,5,1,1))
para = rbind(para,c(30,50,50,1,0.5,1,25,1,1))
para = rbind(para,c(30,50,50,1,1,1,25,1,1))
para = rbind(para,c(30,50,50,1,5,1,25,1,1))
para = rbind(para,c(30,50,50,1,25,1,25,1,1))


para = rbind(para,c(30,50,50,5,0.5,1,0.5,1,1))
para = rbind(para,c(30,50,50,5,1,1,0.5,1,1))
para = rbind(para,c(30,50,50,5,5,1,0.5,1,1))
para = rbind(para,c(30,50,50,5,25,1,0.5,1,1))
para = rbind(para,c(30,50,50,5,0.5,1,1,1,1))
para = rbind(para,c(30,50,50,5,1,1,1,1,1))
para = rbind(para,c(30,50,50,5,5,1,1,1,1))
para = rbind(para,c(30,50,50,5,25,1,1,1,1))
para = rbind(para,c(30,50,50,5,0.5,1,5,1,1))
para = rbind(para,c(30,50,50,5,1,1,5,1,1))
para = rbind(para,c(30,50,50,5,5,1,5,1,1))
para = rbind(para,c(30,50,50,5,25,1,5,1,1))
para = rbind(para,c(30,50,50,5,0.5,1,25,1,1))
para = rbind(para,c(30,50,50,5,1,1,25,1,1))
para = rbind(para,c(30,50,50,5,5,1,25,1,1))
para = rbind(para,c(30,50,50,5,25,1,25,1,1))


para = rbind(para,c(30,20,80,1,0.5,1,0.5,1,1))
para = rbind(para,c(30,20,80,1,1,1,0.5,1,1))
para = rbind(para,c(30,20,80,1,5,1,0.5,1,1))
para = rbind(para,c(30,20,80,1,25,1,0.5,1,1))
para = rbind(para,c(30,20,80,1,0.5,1,1,1,1))
para = rbind(para,c(30,20,80,1,1,1,1,1,1))
para = rbind(para,c(30,20,80,1,5,1,1,1,1))
para = rbind(para,c(30,20,80,1,25,1,1,1,1))
para = rbind(para,c(30,20,80,1,0.5,1,5,1,1))
para = rbind(para,c(30,20,80,1,1,1,5,1,1))
para = rbind(para,c(30,20,80,1,5,1,5,1,1))
para = rbind(para,c(30,20,80,1,25,1,5,1,1))
para = rbind(para,c(30,20,80,1,0.5,1,25,1,1))
para = rbind(para,c(30,20,80,1,1,1,25,1,1))
para = rbind(para,c(30,20,80,1,5,1,25,1,1))
para = rbind(para,c(30,20,80,1,25,1,25,1,1))


para = rbind(para,c(30,20,80,5,0.5,1,0.5,1,1))
para = rbind(para,c(30,20,80,5,1,1,0.5,1,1))
para = rbind(para,c(30,20,80,5,5,1,0.5,1,1))
para = rbind(para,c(30,20,80,5,25,1,0.5,1,1))
para = rbind(para,c(30,20,80,5,0.5,1,1,1,1))
para = rbind(para,c(30,20,80,5,1,1,1,1,1))
para = rbind(para,c(30,20,80,5,5,1,1,1,1))
para = rbind(para,c(30,20,80,5,25,1,1,1,1))
para = rbind(para,c(30,20,80,5,0.5,1,5,1,1))
para = rbind(para,c(30,20,80,5,1,1,5,1,1))
para = rbind(para,c(30,20,80,5,5,1,5,1,1))
para = rbind(para,c(30,20,80,5,25,1,5,1,1))
para = rbind(para,c(30,20,80,5,0.5,1,25,1,1))
para = rbind(para,c(30,20,80,5,1,1,25,1,1))
para = rbind(para,c(30,20,80,5,5,1,25,1,1))
para = rbind(para,c(30,20,80,5,25,1,25,1,1))



para = rbind(para,c(30,80,20,1,0.5,1,0.5,1,1))
para = rbind(para,c(30,80,20,1,1,1,0.5,1,1))
para = rbind(para,c(30,80,20,1,5,1,0.5,1,1))
para = rbind(para,c(30,80,20,1,25,1,0.5,1,1))
para = rbind(para,c(30,80,20,1,0.5,1,1,1,1))
para = rbind(para,c(30,80,20,1,1,1,1,1,1))
para = rbind(para,c(30,80,20,1,5,1,1,1,1))
para = rbind(para,c(30,80,20,1,25,1,1,1,1))
para = rbind(para,c(30,80,20,1,0.5,1,5,1,1))
para = rbind(para,c(30,80,20,1,1,1,5,1,1))
para = rbind(para,c(30,80,20,1,5,1,5,1,1))
para = rbind(para,c(30,80,20,1,25,1,5,1,1))
para = rbind(para,c(30,80,20,1,0.5,1,25,1,1))
para = rbind(para,c(30,80,20,1,1,1,25,1,1))
para = rbind(para,c(30,80,20,1,5,1,25,1,1))
para = rbind(para,c(30,80,20,1,25,1,25,1,1))

para = rbind(para,c(30,80,20,5,0.5,1,0.5,1,1))
para = rbind(para,c(30,80,20,5,1,1,0.5,1,1))
para = rbind(para,c(30,80,20,5,5,1,0.5,1,1))
para = rbind(para,c(30,80,20,5,25,1,0.5,1,1))
para = rbind(para,c(30,80,20,5,0.5,1,1,1,1))
para = rbind(para,c(30,80,20,5,1,1,1,1,1))
para = rbind(para,c(30,80,20,5,5,1,1,1,1))
para = rbind(para,c(30,80,20,5,25,1,1,1,1))
para = rbind(para,c(30,80,20,5,0.5,1,5,1,1))
para = rbind(para,c(30,80,20,5,1,1,5,1,1))
para = rbind(para,c(30,80,20,5,5,1,5,1,1))
para = rbind(para,c(30,80,20,5,25,1,5,1,1))
para = rbind(para,c(30,80,20,5,0.5,1,25,1,1))
para = rbind(para,c(30,80,20,5,1,1,25,1,1))
para = rbind(para,c(30,80,20,5,5,1,25,1,1))
para = rbind(para,c(30,80,20,5,25,1,25,1,1))



case.seq=(para.index*16-15):(para.index*16)
seed.seq = 1:n.rep
n.test=50
for(case in case.seq){
  result = NULL
  est.result = NULL
  for(myseed in seed.seq)
    {
    set.seed(myseed*10)
    
    MSE = NULL
    est.mse = NULL
    
    para[case,1] = para[case,1]
    G = 2
    n = para[case,1]
    p1 = para[case,2]
    p2 = para[case,3]
    p=p1+p2
    # r = min(p,n)/5
    # r1 = min(p,n)/5
    # r2 = min(p,n)/5
    r = r1=r2=1
    r.list = list(r1, r2)
    L = 20
    
    alpha = rep(para[case,4], r)
    alpha1 = rep(para[case,5], r1) #OLS: 0
    alpha2 = rep(para[case,6], r2) #OLS: 0 
    
    
    Siga = diag(p)
    for(i in 1:p)
      for(j in 1:p)
        Siga[i,j] = (0.6)^{abs(i-j)}
    
    Q = randortho(n)
    U = matrix(Q[,1:r], ncol = r)
    U1 = matrix(Q[,r+(1:r1)], ncol = r1)
    U2 = matrix(Q[,r+r1+(1:r2)], ncol = r2)
    J = mvrnorm(n,rep(0,p),Siga)
    J = svd(J,nu=1,nv=1)
    j.scale = para[case,7]
    i1.scale = para[case,8]
    i2.scale = para[case,9]
    J = U%*%(J$d[1:r])%*%t(J$v[,1:r])*j.scale
    # J = U%*%J$d[1:r]%*%t(J$v[,1:r])*j.scale
    J1 = J[,1:p1]
    J2 = J[,(p1+1):p]
    
    I1 = mvrnorm(n,rep(0,p1),Siga[(1:p1),(1:p1)])*i1.scale
    I2 = mvrnorm(n,rep(0,p2),Siga[(1:p2),(1:p2)])*i2.scale
    I1 = svd(I1)
    # I1 = U1%*%diag(I1$d[1:r1])%*%t(I1$v[,1:r1])
    I1 = U1%*%I1$d[1:r1]%*%t(I1$v[,1:r1])
    I2 = svd(I2)
    # I2 = U2%*%diag(I2$d[1:r2])%*%t(I2$v[,1:r2])
    I2 = U2%*%(I2$d[1:r2])%*%t(I2$v[,1:r2])
    # 
    E1 = mvrnorm(n,rep(0,p1),Siga[(1:p1),(1:p1)])
    E2 = mvrnorm(n,rep(0,p2),Siga[(1:p2),(1:p2)])
    X1=J1+I1+E1
    X2=J2+I2+E2
    X = cbind(X1,X2)
    
    # X1 = mvrnorm(n, rep(0, p1), diag(p1))
    # X2 = mvrnorm(n, rep(0, p2), diag(p2))
    
    q = min(n, p)/2
    q1 = min(n, p1)/2
    q2 = min(n, p2)/2
    V = matrix(1,nrow=p,ncol=1)
    V1 = matrix(1,nrow=p1,ncol=1)
    V2 = matrix(1,nrow=p2,ncol=1)
    # V = svd(J)$v[,1:q]%*%rep(1/sqrt(q), q)
    # V1 = svd(I1)$v[,1:q1]%*%rep(1/sqrt(q1), q1)
    # V2 = svd(I2)$v[,1:q2]%*%rep(1/sqrt(q2), q2)
    # X1 = (X1-X1%*%V1%*%t(V1)-X%*%V%*%t(matrix(V[1:p1,],ncol=1)))*noise+X1%*%V1%*%t(V1)*para[case,7]+X%*%V%*%t(matrix(V[1:p1,],ncol=1))*para[case,9]
    # V1 = V1*para[case,7]
    # X2 = (X2-X2%*%V2%*%t(V2)-X%*%V%*%t(matrix(V[(p1+1):p,],ncol=1)))*noise+X2%*%V2%*%t(V2)*para[case,8]+X%*%V%*%t(matrix(V[(p1+1):p,],ncol=1))*para[case,9]
    # V2 = V2*para[case,8]
    # X = cbind(X1, X2)
    e = rnorm(n)*.2
    Y = J%*%V%*%alpha + I1%*%V1%*%alpha1 + I2%*%V2%*%alpha2 + e
    X.list = list(t(X1), t(X2))
    # I1 = X1%*%V1%*%t(V1)
    # I2 = X2%*%V2%*%t(V2)
    # J = X%*%V%*%t(V)
    ajive.result = ajive::ajive(list(X1,X2),initial_signal_ranks = c(r1,r2),joint_rank = r)
    est.mse[["AJIVE.i1"]]= norm(ajive.result$block_decomps[[1]]$individual$full-I1,"2")
    est.mse[["AJIVE.i2"]]= norm(ajive.result$block_decomps[[2]]$individual$full-I2,"2")
    est.mse[["AJIVE.J"]]= norm(cbind(ajive.result$block_decomps[[1]]$joint$full,ajive.result$block_decomps[[2]]$joint$full)-J,"2")
    
    ajive.x = cbind(ajive.result$block_decomps[[1]]$joint$full,ajive.result$block_decomps[[2]]$joint$full,ajive.result$block_decomps[[1]]$individual$full,ajive.result$block_decomps[[2]]$individual$full)
    ajive.beta = ginv(t(ajive.x)%*%ajive.x)%*%t(ajive.x)%*%Y
    
    
    Q = randortho(n.test)
    U = matrix(Q[,1:r], ncol = r)
    U1 = matrix(Q[,r+(1:r1)], ncol = r1)
    U2 = matrix(Q[,r+r1+(1:r2)], ncol = r2)
    J.test = mvrnorm(n.test,rep(0,p),Siga)
    J.test = svd(J.test)
    # J = U%*%diag(J$d[1:r])%*%t(J$v[,1:r])
    
    J.test = U%*%J.test$d[1:r]%*%t(J.test$v[,1:r])*j.scale
    J1.test = J.test[,1:p1]
    J2.test = J.test[,(p1+1):p]
    
    I1.test = mvrnorm(n.test,rep(0,p1),Siga[(1:p1),(1:p1)])
    I2.test = mvrnorm(n.test,rep(0,p2),Siga[(1:p2),(1:p2)])
    I1.test = svd(I1.test)
    # I1 = U1%*%diag(I1$d[1:r1])%*%t(I1$v[,1:r1])
    I1.test = U1%*%I1.test$d[1:r1]%*%t(I1.test$v[,1:r1])*i1.scale
    I2.test = svd(I2.test)
    # I2 = U2%*%diag(I2$d[1:r2])%*%t(I2$v[,1:r2])
    I2.test = U2%*%(I2.test$d[1:r2])%*%t(I2.test$v[,1:r2])*i2.scale
    
    E1 = mvrnorm(n.test,rep(0,p1),Siga[(1:p1),(1:p1)])
    E2 = mvrnorm(n.test,rep(0,p2),Siga[(1:p2),(1:p2)])
    X1=J.test[,1:p1]+I1.test+E1
    X2=J.test[,(p1+1):p]+ I2.test+E2
    e = rnorm(n.test)*.2
    X.test = cbind(X1, X2)
    Y.test = J.test%*%V%*%alpha + I1.test%*%V1%*%alpha1 + I2.test%*%V2%*%alpha2 + e
    
    X.test.list = list(t(X1), t(X2))
    ajive.result = ajive::ajive(list(X1,X2),initial_signal_ranks = c(r1,r2),joint_rank = r)
    ajive.x = cbind(ajive.result$block_decomps[[1]]$joint$full,ajive.result$block_decomps[[2]]$joint$full,ajive.result$block_decomps[[1]]$individual$full,ajive.result$block_decomps[[2]]$individual$full)
    
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
    ml.continuum.pcr = continuum.multisource.iter.v1(X.list, Y, maxiter=1000,lambda = 0, gam = 1e10, rankJ = r, rankA = c(r1, r2),center.X = F, scale.X = F, center.Y = F, scale.Y = F,  orthIndiv = F)
    ml.continuum = ml.continuum.pcr
    # ml = decomposeX(X.test.list, ml.continuum$U[[ml.continuum$nrun]], ml.continuum$W[[ml.continuum$nrun]], ml.continuum$centerValues.X, ml.continuum$scaleValues.X)
    beta.Cind = ml.continuum$beta.Cind
    Yhat.heter = (X1)%*%beta.Cind[[1]]+(X2)%*%beta.Cind[[2]]
    # Yhat.heter = t(SimData$x.test[[1]])%*%beta.Cind[[1]]+t(SimData$x.test[[2]])%*%beta.Cind[[2]]
    # x.test = cbind(t(SimData$x.test[[1]]),t(SimData$x.test[[2]]))
    # mean((SimData$y.test - as.numeric(ml.continuum$intercept) - x.test%*%ml.continuum$beta.C - Yhat.heter)^2)
    
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
    
    
    ml.continuum.pls = continuum.multisource.iter.v2(X.list, Y, maxiter=1000,lambda = 0, gam = 1, rankJ = r, rankA = c(r1, r2),center.X = F, scale.X = F, center.Y = F, scale.Y = F,  orthIndiv = F)
    ml.continuum = ml.continuum.pls
    # ml = decomposeX(X.test.list, ml.continuum$U[[ml.continuum$nrun]], ml.continuum$W[[ml.continuum$nrun]], ml.continuum$centerValues.X, ml.continuum$scaleValues.X)
    beta.Cind = ml.continuum$beta.Cind
    Yhat.heter = (X1)%*%beta.Cind[[1]]+(X2)%*%beta.Cind[[2]]
    #
    MSE[["continuum.pls2"]]= mean((Y.test - as.numeric(ml.continuum$intercept) - X.test%*%ml.continuum$beta.C - Yhat.heter)^2)
    est.mse[["continuum.pls2.i1"]]= norm(t(ml.continuum$I[[1]])-I1,"2")
    est.mse[["continuum.pls2.i2"]]= norm(t(ml.continuum$I[[2]])-I2,"2")
    est.mse[["continuum.pls2.j"]]= norm(t(ml.continuum$J)-(J),"2")

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

    
    
    
    
    ml.ridge = cv.glmnet(x = X, y = Y, lambda = 10^seq(1,-4,length=100),alpha = 0, standardize = F, intercept = F)
    MSE[["ridge"]]= mean((Y.test -predict(ml.ridge, newx = X.test, s = ml.ridge$lambda.min) )^2)
    
    
    ml.ridge = cv.glmnet(x = cbind(J,I1,I2), y = Y, lambda = 10^seq(1,-14,length=100),alpha = 0, standardize = F, intercept = F)
    MSE[["ridge.seperate"]]= mean((Y.test -predict(ml.ridge, newx = cbind(J.test,I1.test,I2.test), s = ml.ridge$lambda.min) )^2)
    
    
    ml.pls = plsr(Y~X, validation = "CV", center = T, scale = T)
    ncomp.pls = selectNcomp(ml.pls, method = "randomization", plot = F)
    MSE[["pls"]]= mean((predict(ml.pls, newdata = X.test, ncomp = r+r1+r2)[,,1] - Y.test)^2)

    ml.pcr = pcr(Y~X, validation = "CV", center = T, scale = T)
    ncomp.pcr = selectNcomp(ml.pcr, method = "randomization", plot = F)
    MSE[["pcr"]]= mean((predict(ml.pcr, newdata = X.test, ncomp = r+r1+r2)[,,1] - Y.test)^2)


    # ml.ridge = cv.glmnet(x = X.scaled, y = Y, alpha = 0, standardize = T, intercept = T)
    # MSE[["ridge.scaled"]]= mean((Y.test -predict(ml.ridge, newx = X.test.scaled, s = ml.ridge$lambda.min) )^2)
    # 
    # ml.pls = plsr(Y~X.scaled, validation = "CV", center = T, scale = T)
    # ncomp.pls = selectNcomp(ml.pls, method = "randomization", plot = F)
    # MSE[["pls.scaled"]]= mean((predict(ml.pls, newdata = X.test.scaled, ncomp = r+r1+r2)[,,1] - Y.test)^2)
    # 
    # ml.pcr = pcr(Y~X.scaled, validation = "CV", center = T, scale = T)
    # ncomp.pcr = selectNcomp(ml.pcr, method = "randomization", plot = F)
    # MSE[["pcr.scaled"]]= mean((predict(ml.pcr, newdata = X.test.scaled, ncomp = r+r1+r2)[,,1] - Y.test)^2)
    # 
    # 
    # ml.continuum.pcr = continuum.multisource.iter.v1(X.list.scaled, Y, lambda = 0, gam = 1e10, rankJ = r, rankA = c(r1, r2),center.X = F, scale.X = F, center.Y = F, scale.Y = F,  orthIndiv = F)
    # ml.continuum = ml.continuum.pcr
    # # ml = decomposeX(X.test.list, ml.continuum$U[[ml.continuum$nrun]], ml.continuum$W[[ml.continuum$nrun]], ml.continuum$centerValues.X, ml.continuum$scaleValues.X)
    # beta.Cind = ml.continuum$beta.Cind
    # Yhat.heter = t(X.test.list.scaled[[1]])%*%beta.Cind[[1]]+t(X.test.list.scaled[[2]])%*%beta.Cind[[2]]
    # #
    # MSE[["continuum.pcr.scaled"]]= mean((Y.test - as.numeric(ml.continuum$intercept) - X.test.scaled%*%ml.continuum$beta.C - Yhat.heter)^2)
    # 
    # ml.continuum = continuum.multisource.iter.v1(X.list.scaled, Y, maxiter=1000, lambda = 0, gam = 1, rankJ = r, rankA = c(r1, r2),center.X = F, scale.X = F, center.Y = F, scale.Y = F,  orthIndiv = F)
    # beta.Cind = ml.continuum$beta.Cind
    # Yhat.heter = t(X.test.list.scaled[[1]])%*%beta.Cind[[1]]+t(X.test.list.scaled[[2]])%*%beta.Cind[[2]]
    # #
    # MSE[["continuum.pls.scaled"]]= mean((Y.test - as.numeric(ml.continuum$intercept) - X.test.scaled%*%ml.continuum$beta.C - Yhat.heter)^2)
    # # print(sort(MSE))
    # 
    MSE = MSE/norm(Y.test)
    
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
  print(case)
  print(para[case,])
  print(cbind(sort(apply(result,2,mean)),apply(result,2,sd)[order(apply(result,2,mean))]))
  # colnames(result) = c("cr.pcr","cr.pls","cr.pls.n.orth","ridge","pls","pcr")
  # result = mapply("/", result,length(seed.seq) , SIMPLIFY = FALSE)
  file.name = paste0("multi_source_result_pls_n=",n,"_p1=",p1,"_p2=",p2,"_a=",alpha[1],"_a1=",alpha1[1],"_a2=",alpha2[1],"_scale1=",para[case,7],"_scale2=",para[case,8],"_scalej=",para[case,9],"_sigma3")
  # write.table(result, file = file.name, sep = ',', append = T, col.names = F, row.names = F)
  save(result,file = paste0("3_1_noise2/",file.name,".rData"))
  file.name = paste0("multi_source_est_result_pls_n=",n,"_p1=",p1,"_p2=",p2,"_a=",alpha[1],"_a1=",alpha1[1],"_a2=",alpha2[1],"_scale1=",para[case,7],"_scale2=",para[case,8],"_scalej=",para[case,9],"_sigma3")
  # print(c(case,colMeans(result)))
  save(est.result,file = paste0("3_1_noise2/",file.name,".rData"))
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
