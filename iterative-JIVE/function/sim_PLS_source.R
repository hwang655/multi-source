library(MASS)
library(pls)
library(purrr)
library(pracma)
library(glmnet)
library(R.matlab)
library(r.jive)
library(ajive)
source("./function/jive_continuum.R")

Y.pcr.list = NULL
Y.pls.list = NULL

result = NULL

myseed=1
seed.seq = 1:20
for(myseed in seed.seq){
  MSE = NULL
set.seed(134+myseed)
r = 2
r1 = 2
r2 = 2
n = 50
p1 = 10
p2 = 100
L = c(1, 1)
L1 = c(1, 1)
L2 = c(1, 1)
p = p1 + p2
G = 2
alpha = c(5, 3)*10
alpha1 = c(1, 0.5)*100
alpha2 = c(2, 0.5)*10
# beta = c(rep(1, p))*10
# beta1 = c(rep(0, 5), rep(1, p1-5))*5
# beta2 = c(rep(0, 5), rep(-1, p2-5))*5

Q = randortho(n, type = "orthonormal")
U = Q[,1:r]%*%DIAG(L)
W1 = Q[,(r+1):(r+r1)]%*%DIAG(L1)
W2 = Q[,(r+r1+1):(r+r1+r2)]%*%DIAG(L2)

P = randortho(p1, type = "orthonormal")
S = P[1:r,]
S1 = P[(r+1):(r+r1),]
P = randortho(p2, type = "orthonormal")
S = cbind(S,P[1:r,])/2
S2 = P[(r+1):(r+r2),]

# V = randortho(n, type = "orthonormal")
# Z = matrix(V[,1:r], nrow = n)
# V = randortho(n, type = "orthonormal")
# Z1 = matrix(V[,1:r1], nrow = n)
# V = randortho(n, type = "orthonormal")
# Z2 = matrix(V[,1:r2], nrow = n)
# a = t(P)%*%DIAG(c(1/L, rep(0, n - r)))%*%Z
# a1 = t(P)%*%DIAG(c(rep(0, r1), 1/L1, rep(0, n - 2*r1)))%*%Z1
# a2 = t(P)%*%DIAG(c(rep(0, r1+r2), 1/L1, rep(0, n - 2*r1-r2)))%*%Z2

E1 = (mvrnorm(n, rep(0, p1), diag(p1))*0.01)*10
E2 = (mvrnorm(n, rep(0, p2), diag(p2))*0.01)*10

J = U%*%S
J1 = J[,1:p1]*5
J2 = J[,(p1+1):(p1+p2)]
I1 = W1%*%S1
I2 = W2%*%S2

X1 = (J1 + I1 + E1)
X2 = J2 + I2 + E2
# X1 = J1  + E1
# X2 = J2  + E2

e = rnorm(n)
#Y = t(J)%*%beta + t(I1)%*%(beta1) + t(I2)%*%(beta2) + e
#Y = t(S)%*%SOLVE(DIAG(L))^(1/2)%*%alpha + t(S1)%*%SOLVE(DIAG(L1))^(1/2)%*%alpha1 + t(S2)%*%SOLVE(DIAG(L1))^(1/2)%*%alpha2 + e

#Y = t(S)%*%alpha + t(S1)%*%alpha1 + t(S2)%*%alpha2 + e
#Y = t(S)%*%DIAG(L)%*%Z[1:r,]%*%alpha + t(S1)%*%DIAG(L1)%*%Z1[1:r1,]%*%alpha1 +
#  t(S2)%*%DIAG(L2)%*%Z2[1:r2,]%*%alpha2 + e
# Y = t(J)%*%J%*%a%*%alpha + t(I1)%*%I1%*%a1%*%alpha1 + t(I2)%*%I2%*%a2%*%alpha2 + e
# Y = J%*%t(S)%*%alpha + I1%*%t(S1)%*%alpha1+I2%*%t(S2)%*%alpha2+e
X = cbind(X1, X2)
Y = X%*%t(S)%*%alpha + X1%*%t(S1)%*%alpha1+X2%*%t(S2)%*%alpha2+e
# t(a)%*%t(J)%*%(J)%*%Y
# t(a1)%*%t(I1)%*%I1%*%Y
# t(a2)%*%t(I2)%*%I2%*%Y

X.list = list(t(X1), t(X2))
# ----------------------------------------------- testing ----------------------------------------------------
# P = randortho(n, type = "orthonormal")
# S = matrix(P[1:r,], ncol = n)
# S1 = matrix(P[(r+1):(r+r1),], ncol = n)
# S2 = matrix(P[(r+r1+1):(r+r1+r2),], ncol = n)

Q = randortho(n, type = "orthonormal")
U = Q[,1:r]%*%DIAG(L)
W1 = Q[,(r+1):(r+r1)]%*%DIAG(L1)
W2 = Q[,(r+r1+1):(r+r1+r2)]%*%DIAG(L2)

#Z = randortho(r, type = "orthonormal")
#Z1 = randortho(r1, type = "orthonormal")
#Z2 = randortho(r2, type = "orthonormal")
# a = t(S)%*%DIAG(1/L)%*%Z
# a1 = t(S1)%*%DIAG(1/L1)%*%Z1
# a2 = t(S2)%*%DIAG(1/L2)%*%Z2

# V = randortho(n, type = "orthonormal")
# P = randortho(n, type = "orthonormal")
# Z = matrix(P[1:r,], ncol = n)
# S = Z%*%V
# Z1 = matrix(P[(r+1):(r+r1),], ncol = n)
# S1 = Z1%*%V
# Z2 = matrix(P[(r+r1+1):(r+r1+r2),], ncol = n)
# S2 = Z2%*%V

E1 = (mvrnorm(n, rep(0, p1), diag(p1))*0.01)*10
E2 = (mvrnorm(n, rep(0, p2), diag(p2))*0.01)*10

J = U%*%S
J1 = J[,1:p1]*5
J2 = J[,(p1+1):(p1+p2)]
I1 = W1%*%S1
I2 = W2%*%S2
# 
X1 = (J1 + I1 + E1)
X2 = J2 + I2 + E2
# X1 = J1  + E1
# X2 = J2  + E2


e = rnorm(n)
#Y.test = t(S)%*%SOLVE(DIAG(L))^(1/2)%*%alpha + t(S1)%*%SOLVE(DIAG(L1))^(1/2)%*%alpha1 + t(S2)%*%SOLVE(DIAG(L1))^(1/2)%*%alpha2 + e
#Y.test = t(S)%*%alpha + t(S1)%*%alpha1 + t(S2)%*%alpha2 + e
#Y.test = t(J)%*%beta + t(I1)%*%(beta1) + t(I2)%*%(beta2) + e
# Y.test = t(J)%*%J%*%a%*%aY = J%*%t(S)%*%alpha + I1%*%t(S1)%*%alpha1+I2%*%t(S2)%*%alpha2+elpha + t(I1)%*%I1%*%a1%*%alpha1 + t(I2)%*%I2%*%a2%*%alpha2 + e
# Y.test = J%*%t(S)%*%alpha + I1%*%t(S1)%*%alpha1+I2%*%t(S2)%*%alpha2+e
X.test = cbind(X1, X2)
Y.test = X.test%*%t(S)%*%alpha + X1%*%t(S1)%*%alpha1+X2%*%t(S2)%*%alpha2+e
X.test.list = list(t(X1), t(X2))


# writeMat(paste0("./data/joint",myseed,".mat"),trainx = X,trainy = Y,testx = X.test,testy = Y.test)

# ----------------------------------------------- run model ----------------------------------------------------
# ml.jive = jive(X.list, rankJ = r, rankA = c(r1, r2), method = "given", orthIndiv = F,scale=F,center=F)
# ml.continuum.pcr = continuum.multisource.iter.v1(X.list, Y, lambda = 0, gam = 1e10, rankJ = r, rankA = c(0,0))
# r1=r2=0
ml.continuum.pcr = continuum.multisource.iter.v1(X.list, Y, lambda = 0, gam = 1e10, rankJ = r, rankA = c(r1, r2),orthIndiv = F,center.X = F,center.Y = F,scale.X=F,scale.Y=F)
ml.continuum = ml.continuum.pcr
# ml = decomposeX(X.test.list, ml.continuum$U[[ml.continuum$nrun]], ml.continuum$W[[ml.continuum$nrun]], ml.continuum$centerValues.X, ml.continuum$scaleValues.X)
beta.Cind = ml.continuum$beta.Cind
Yhat.heter = (X1)%*%beta.Cind[[1]]+(X2)%*%beta.Cind[[2]]
Y.pcr.list = list.append(Y.pcr.list,Yhat.heter)
MSE= c(MSE,mean((Y.test - as.numeric(ml.continuum$intercept) - X.test%*%ml.continuum$beta.C - Yhat.heter)^2))

# ml.continuum.pls = continuum.multisource.iter.v1(X.list, Y, lambda = 0, gam = 1, rankJ = r, rankA = c(r1, r2),orthIndiv = F,center.X=F, center.Y = F,scale.X=F,scale.Y=F)
# # ml.continuum.pls = continuum.multisource.iter.v1(X.list, Y, lambda = 0, gam = 1, rankJ = r, rankA = c(r1, r2))
# ml.continuum = ml.continuum.pls
# beta.Cind = ml.continuum$beta.Cind
# Yhat.heter = (X1)%*%beta.Cind[[1]]+(X2)%*%beta.Cind[[2]]

# Y.pls.list = list(Y.pls.list,Yhat.heter)
# MSE= c(MSE,mean((Y.test - as.numeric(ml.continuum$intercept) - X.test%*%ml.continuum$beta.C - Yhat.heter)^2))


ml.ridge = cv.glmnet(x = X, y = Y, alpha = 0, standardize = F, intercept = F)
MSE= c(MSE, mean((Y.test -predict(ml.ridge, newx = X.test, s = ml.ridge$lambda.min) )^2))

ml.pls = plsr(Y~X, validation = "CV", center = F, scale = F)
ncomp.pls = selectNcomp(ml.pls, method = "randomization", plot = F)
MSE= c(MSE,mean((predict(ml.pls, newdata = X.test, ncomp = r+r1+r2)[,,1] - Y.test)^2))

ml.pcr = pcr(Y~X, validation = "CV", center = F, scale = F)
ncomp.pcr = selectNcomp(ml.pcr, method = "randomization", plot = F)
MSE= c(MSE,mean((predict(ml.pcr, newdata = X.test, ncomp = r+r1+r2)[,,1] - Y.test)^2))

# ajive_result = ajive(list(t(X.list[[1]]),t(X.list[[2]])),c(r1,r2),joint_rank = r)
# ajive_result$block_decomps

print(MSE)
print(ml.continuum.pcr$nrun)
result = rbind(result,MSE)
}
# colnames(result) = c("cr.pcr","cr.pls")
# rownames(result) = seed.seq
row.names(result) = c("cr.pcr","cr.pls","ridge","pls","pcr")
# file.name = "result.csv"
# write.table(t(c(myseed, as.vector(MSE))), file = file.name, sep = ',', append = T, col.names = F, row.names = F)
