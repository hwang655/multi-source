library(MASS)
library(pls)
library(purrr)
library(pracma)

r = 2
r1 = 2
r2 = 2
n = 50
p1 = 50
p2 = 50
L = c(50, 10)
L1 = c(1, 1)
L2 = c(1, 1)
p = p1 + p2
G = 2
alpha = rep(1, r)
alpha1 = rep(1, r1)
alpha2 = rep(1, r2)

Q = randortho(p, type = "orthonormal")
U = Q[,1:r]%*%diag(L)
Q1 = randortho(p1, type = "orthonormal")
Q2 = randortho(p2, type = "orthonormal")
W1 = Q1[,1:r1]%*%DIAG(L1)
W2 = Q2[,1:r2]%*%DIAG(L2)

P = randortho(n, type = "orthonormal")
S = matrix(P[1:r,], ncol = n)
S1 = matrix(P[(r+1):(r+r1),], ncol = n)
S2 = matrix(P[(r+r1+1):(r+r1+r2),], ncol = n)

E1 = t(mvrnorm(n, rep(0, p1), diag(p1))*10)
E2 = t(mvrnorm(n, rep(0, p2), diag(p2))*10)

J = U%*%S
J1 = J[1:p1,]
J2 = J[(p1+1):(p1+p2),]
I1 = W1%*%S1
I2 = W2%*%S2

X1 = J1 + I1 + E1
X2 = J2 + I2 + E2

e = rnorm(n)*0.01
Y = t(S)%*%alpha + t(S1)%*%alpha1 + t(S2)%*%alpha2 + e

X = cbind(t(X1), t(X2))
X.list = list(X1, X2)

ajive.result = ajive::ajive(list(t(X1),t(X2)),initial_signal_ranks = c(r1,r2),joint_rank = r)
est.mse[["AJIVE.i1"]]= norm(ajive.result$block_decomps[[1]]$individual$full-t(I1),"2")
est.mse[["AJIVE.i2"]]= norm(ajive.result$block_decomps[[2]]$individual$full-t(I2),"2")
est.mse[["AJIVE.J"]]= norm(cbind(ajive.result$block_decomps[[1]]$joint$full,ajive.result$block_decomps[[2]]$joint$full)-t(J),"2")
temp = cbind(ajive.result$block_decomps[[1]]$joint$full,ajive.result$block_decomps[[2]]$joint$full)
ajive.x = svd(temp)$u[,1:r]%*%diag(svd(temp)$d[1:r],r)
temp = ajive.result$block_decomps[[1]]$individual$full
ajive.x = cbind(ajive.x,svd(temp)$u[,1:r1]%*%diag(svd(temp)$d[1:r1],r1))
temp = ajive.result$block_decomps[[2]]$individual$full
ajive.x = cbind(ajive.x,svd(temp)$u[,1:r2]%*%diag(svd(temp)$d[1:r2],r2))
ajive.beta = ginv(t(ajive.x)%*%ajive.x)%*%t(ajive.x)%*%Y

# ----------------------------------------------- testing ----------------------------------------------------
P = randortho(n, type = "orthonormal")
S = matrix(P[1:r,], ncol = n)
S1 = matrix(P[(r+1):(r+r1),], ncol = n)
S2 = matrix(P[(r+r1+1):(r+r1+r2),], ncol = n)

E1 = t(mvrnorm(n, rep(0, p1), diag(p1))*0.01)
E2 = t(mvrnorm(n, rep(0, p2), diag(p2))*0.01)

J = U%*%S
J1 = J[1:p1,]
J2 = J[(p1+1):(p1+p2),]
I1 = W1%*%S1
I2 = W2%*%S2

X1 = J1 + I1 + E1
X2 = J2 + I2 + E2

e = rnorm(n)*0.01
Y.test = t(S)%*%alpha + t(S1)%*%alpha1 + t(S2)%*%alpha2 + e

X.test = cbind(t(X1), t(X2))
X.test.list = list(X1, X2)
ajive.result = ajive::ajive(list(t(X1),t(X2)),initial_signal_ranks = c(r1,r2),joint_rank = r)
temp = cbind(ajive.result$block_decomps[[1]]$joint$full,ajive.result$block_decomps[[2]]$joint$full)
ajive.x = svd(temp)$u[,1:r]%*%diag(svd(temp)$d[1:r],r)
temp = ajive.result$block_decomps[[1]]$individual$full
ajive.x = cbind(ajive.x,svd(temp)$u[,1:r1]%*%diag(svd(temp)$d[1:r1],r1))
temp = ajive.result$block_decomps[[2]]$individual$full
ajive.x = cbind(ajive.x,svd(temp)$u[,1:r2]%*%diag(svd(temp)$d[1:r2],r2))
MSE["AJIVE"] = mean((Y.test - ajive.x%*%ajive.beta)^2)

# ----------------------------------------------- run model ----------------------------------------------------
result = list()
# ml.jive = jive.multisource(X.list, rankJ = r, rankA = c(r1, r2), method = "given", orthIndiv = T)

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
est.mse[["continuum.pcr.i1"]]= norm(t(ml.continuum$I[[1]])-t(I1),"2")
est.mse[["continuum.pcr.i2"]]= norm(t(ml.continuum$I[[2]])-t(I2),"2")
est.mse[["continuum.pcr.j"]]= norm(t(ml.continuum$J)-t(J),"2")

ml.continuum.pls = continuum.multisource.iter.v1(X.list, Y, maxiter=1000,lambda = 0, gam = 1, rankJ = r, rankA = c(r1, r2),center.X = F, scale.X = F, center.Y = F, scale.Y = F,  orthIndiv = F)
ml.continuum = ml.continuum.pls
# ml = decomposeX(X.test.list, ml.continuum$U[[ml.continuum$nrun]], ml.continuum$W[[ml.continuum$nrun]], ml.continuum$centerValues.X, ml.continuum$scaleValues.X)
beta.Cind = ml.continuum$beta.Cind
Yhat.heter = (X1)%*%beta.Cind[[1]]+(X2)%*%beta.Cind[[2]]
# 
MSE[["continuum.pls"]]= mean((Y.test - as.numeric(ml.continuum$intercept) - X.test%*%ml.continuum$beta.C - Yhat.heter)^2)
est.mse[["continuum.pls.i1"]]= norm(t(ml.continuum$I[[1]])-t(I1),"2")
est.mse[["continuum.pls.i2"]]= norm(t(ml.continuum$I[[2]])-t(I2),"2")
est.mse[["continuum.pls.j"]]= norm(t(ml.continuum$J)-t(J),"2")

ml.continuum.pcr = continuum.multisource.iter.v1(X.list, Y, lambda = 0, gam = 1e10, rankJ = r, rankA = c(r1, r2),center.X = T, scale.X = T, center.Y = T, scale.Y = T,  orthIndiv = F)
ml.continuum = ml.continuum.pcr
# ml = decomposeX(X.test.list, ml.continuum$U[[ml.continuum$nrun]], ml.continuum$W[[ml.continuum$nrun]], ml.continuum$centerValues.X, ml.continuum$scaleValues.X)
beta.Cind = ml.continuum$beta.Cind
Yhat.heter = (X1)%*%beta.Cind[[1]]+(X2)%*%beta.Cind[[2]]
# 
MSE[["continuum.pcr.stand"]]= mean((Y.test - as.numeric(ml.continuum$intercept) - X.test%*%ml.continuum$beta.C - Yhat.heter)^2)
est.mse[["continuum.pcr.stand.i1"]]= norm(t(ml.continuum$I[[1]])-t(I1),"2")
est.mse[["continuum.pcr.stand.i2"]]= norm(t(ml.continuum$I[[2]])-t(I2),"2")
est.mse[["continuum.pcr.stand.j"]]= norm(t(ml.continuum$J)-t(J),"2")

ml.continuum.pls = continuum.multisource.iter.v1(X.list, Y, lambda = 0, gam = 1, rankJ = r, rankA = c(r1, r2),center.X = T, scale.X = T, center.Y = T, scale.Y = T,  orthIndiv = F)
ml.continuum = ml.continuum.pls
# ml = decomposeX(X.test.list, ml.continuum$U[[ml.continuum$nrun]], ml.continuum$W[[ml.continuum$nrun]], ml.continuum$centerValues.X, ml.continuum$scaleValues.X)
beta.Cind = ml.continuum$beta.Cind
Yhat.heter = (X1)%*%beta.Cind[[1]]+(X2)%*%beta.Cind[[2]]
# 
MSE[["continuum.pls.stand"]]= mean((Y.test - as.numeric(ml.continuum$intercept) - X.test%*%ml.continuum$beta.C - Yhat.heter)^2)
est.mse[["continuum.pls.stand.i1"]]= norm(t(ml.continuum$I[[1]])-t(I1),"2")
est.mse[["continuum.pls.stand.i2"]]= norm(t(ml.continuum$I[[2]])-t(I2),"2")
est.mse[["continuum.pls.stand.j"]]= norm(t(ml.continuum$J)-t(J),"2")





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
