set.seed(0)
n <- 60
d <- 2
n2 <- 20
f1 <- function(x) {sin(2*pi*x[1]) + sin(2*pi*x[2])}
#f1 <- function(x) {abs(sin(2*pi*x[1])) + 10*x[2]}
#f1 <- TestFunctions::RFF_get(D=d)
X1 <- matrix(runif(n*d),n,d)
Z1 <- apply(X1,1,f1) + rnorm(n, 0, 1e-3)
X2 <- matrix(runif(n2*d),n2,d)
Z2 <- apply(X2,1,f1)
Xall <- rbind(X1, X2)
Zall <- c(Z1, Z2)
XX1 <- matrix(runif(5*d),nrow=5)
ZZ1 <- apply(XX1, 1, f1)
system.time(u <- IGP(package='GauPro',X=X1,Z=Z1))
u$theta()
cbind(u$predict(XX1), ZZ1)
u$predict.se(XX1)
cf::cf_func(u$predict,batchmax = 100, pts=X1)
u$update(Xnew=X2,Znew=Z2)
u$predict(XX1)
cf::cf_func(u$predict,batchmax = 100, pts=Xall)
u$delete()
cf::cf_func(function(xx)sin(8*pi*xx[1]) + sin(8*pi*xx[2]))
cf::cf_data(X1,Z1)
cf::cf_data(X2,Z2)
cf::cf_data(rbind(X1,X2),c(Z1,Z2))
cf::cf_data(u$X,u$Z)
cf::cf_func(u$grad_norm)

u <- UGP$new(package='laGP')
u$update(Xnew=X1, Znew=Z1)
cf::cf_func(u$predict,batchmax = 100);points(X1)
u$update(Xnew=X2,Znew=Z2)
cf::cf_func(u$predict,batchmax = 100);points(rbind(X1,X2))
u$delete()



# Test 1 D
set.seed(0)
n <- 10
d <- 1
n2 <- 5
f1 <- function(x) {sin(2*pi*x[1]) + rnorm(1,0,.005)}
X1 <- matrix(sort(runif(n*d)),n,d)
Z1 <- apply(X1,1,f1) #+ rnorm(n, 0, 1e-3)
X2 <- matrix(runif(n2*d),n2,d)
Z2 <- apply(X2,1,f1)
XX1 <- matrix(runif(10),ncol=1)
ZZ1 <- apply(XX1, 1, f1)
system.time(u <- IGP(package='mlegp',X=X1,Z=Z1, estimate.nugget=T, nugget=c(10,10,10,10,10,1,1,1,1,1)))
plot(X1, Z1)
cbind(u$predict(XX1), ZZ1)
u$predict.se(XX1)
curve(u$predict(matrix(x, ncol=1)));points(X1,Z1, col=2, cex=2, pch=19)
curve(u$predict(matrix(x, ncol=1)) + 2*u$predict.se(matrix(x, ncol=1)), add=T, col=3)
curve(u$predict(matrix(x, ncol=1)) - 2*u$predict.se(matrix(x, ncol=1)), add=T, col=3)
u$update(Xnew=X2,Znew=Z2)
curve(u$predict(matrix(x, ncol=1)));points(rbind(X1,X2),c(Z1,Z2), col=2, cex=2, pch=19)
curve(u$predict(matrix(x, ncol=1)) + 2*u$predict.se(matrix(x, ncol=1)), add=T, col=3)
curve(u$predict(matrix(x, ncol=1)) - 2*u$predict.se(matrix(x, ncol=1)), add=T, col=3)
u$delete()



# Test 1 D reversion to mean
set.seed(0)
n <- 10
d <- 1
n2 <- 5
f1 <- function(x) {0 *sin(2*pi*x[1]) + rnorm(1,0,.05)}
X1 <- matrix(((runif(n*d)-.5) * .1 + 1) %% 1 ,n,d)
Z1 <- apply(X1,1,f1) #+ rnorm(n, 0, 1e-3)
X2 <- matrix(runif(n2*d),n2,d)
Z2 <- apply(X2,1,f1)
XX1 <- matrix(runif(10),ncol=1)
ZZ1 <- apply(XX1, 1, f1)
system.time(u <- UGP$new(package='mlegp',X=X1,Z=Z1, estimate.nugget=T))
cbind(u$predict(XX1), ZZ1)
u$predict.se(XX1)
curve(u$predict(matrix(x, ncol=1)), ylim=c(-.1,.1));points(X1,Z1, col=2, cex=2, pch=19)
curve(u$predict(matrix(x, ncol=1)) + 2*u$predict.se(matrix(x, ncol=1)), add=T, col=3)
curve(u$predict(matrix(x, ncol=1)) - 2*u$predict.se(matrix(x, ncol=1)), add=T, col=3)
u$update(Xnew=X2,Znew=Z2)
curve(u$predict(matrix(x, ncol=1)));points(rbind(X1,X2),c(Z1,Z2), col=2, cex=2, pch=19)
curve(u$predict(matrix(x, ncol=1)) + 2*u$predict.se(matrix(x, ncol=1)), add=T, col=3)
curve(u$predict(matrix(x, ncol=1)) - 2*u$predict.se(matrix(x, ncol=1)), add=T, col=3)
u$delete()




# try tgp
library(tgp)
mod <- tgp::btgpllm(X2,Z2)
predict(mod,X1)

# laGP not working
## prior and inits for theta (d) and nugget (g)
library(laGP)
X3 <- X2
Z3 <- Z2
da <- darg(list(mle=TRUE), X=X3)
ga <- garg(list(mle=TRUE), y=Z3)

## initialize new separable GP
gpi <- newGPsep(X=X3, Z=Z3, d=da$start, g=ga$start, dK=TRUE)

## joint inference for lengthscale and nugget under prior calculated above
mle <- jmleGPsep(gpi, drange=c(da$min, da$max), grange=c(ga$min, ga$max), dab=da$ab, gab=ga$ab, verb=1)

## predict at a new set of XX locations, use lite=TRUE if you don't need a full predictive covariance matrics
p <- predGPsep(gpi, XX=XX1, lite=T)
cf::cf_func(function(xx){predGPsep(gpi,XX=xx,lite=T)$mean},batchmax = 100)

## clean up
deleteGPsep(gpi)







# test laGP updateGP here

# laGP not working
## prior and inits for theta (d) and nugget (g)
da <- darg(list(mle=TRUE), X=X1)
ga <- garg(list(mle=TRUE), y=Z1)

## initialize new separable GP
gpi <- newGPsep(X=X1, Z=Z1, d=da$start, g=ga$start, dK=TRUE)

## joint inference for lengthscale and nugget under prior calculated above
mle <- jmleGPsep(gpi, drange=c(da$min, da$max), grange=c(ga$min, ga$max), dab=da$ab, gab=ga$ab, verb=1)
cf::cf_func(function(xx){predGPsep(gpi,XX=xx,lite=T)$mean},batchmax = 100)
points(X1,pch=19)

updateGPsep(gpsepi=gpi, X=X2, Z=Z2)
cf::cf_func(function(xx){predGPsep(gpi,XX=xx,lite=T)$mean},batchmax = 100)
points(Xall,pch=19)

da <- darg(list(mle=TRUE), X=Xall)
ga <- garg(list(mle=TRUE), y=Zall)
jmleGPsep(gpi, drange=c(da$min, da$max), grange=c(ga$min, ga$max), dab=da$ab, gab=ga$ab, verb=1)



predGPsep(gpi, XX=XX1, lite=T)
cf::cf_func(function(xx){predGPsep(gpi,XX=xx,lite=T)$mean},batchmax = 100)
points(Xall,pch=19)

## clean up
deleteGPsep(gpi)



# Try GPy
library(rPython)
python.exec('import sys')
python.exec("sys.path.insert(0, '/Users/collin/anaconda/lib/python2.7/site-packages/')")
python.exec('print sys.path')
python.exec('import numpy')
python.exec('print numpy')
python.exec('import numpy as np')
python.exec('import scipy')
python.exec('print scipy.__version__')
python.exec('import GPy')
python.exec('print GPy')

python.assign("inputdim", ncol(X1))
python.assign("X1", (X1))
python.assign("y1", Z1)
python.exec('X =  np.matrix(X1)')
python.exec('y = np.matrix(y1).reshape((-1,1))')
python.exec("kernel = GPy.kern.RBF(input_dim=inputdim, variance=1., lengthscale=[1. for iii in range(inputdim)],ARD=True)")
python.exec("gp = GPy.models.GPRegression(X,y,kernel,normalizer=True)")
python.exec("gp.likelihood.variance = 1e-8")
python.exec("gp.optimize(messages=False)")
python.exec("gp.optimize_restarts(num_restarts = 5,  verbose=False)")
python.assign("xp1", XX1)
python.exec("xp = np.asmatrix(xp1)")
python.exec("y_pred, sigma2_pred = gp.predict(np.asarray(xp))")
unlist(python.get("y_pred.tolist()"))



# Try GPy with PythonInR since rPython not working on Windows
library(PythonInR)
pyConnect()
pyExecp('print sys.path')
pyExecp('import numpy')
pyExecp('print numpy')
pyExecp('import numpy as np')
pyExecp('import scipy')
pyExecp('print scipy.__version__')
pyExecp('import GPy')
pyExecp('print GPy')

pySet("inputdim", ncol(X1))
pySet("X1", X1, useNumpy=T)
pySet("y1", Z1, useSetPoly = F)
pyExecp('X =  np.asmatrix(X1)')
pyExecp('y = np.asarray(y1).reshape((-1,1))')
pyExecp("kernel = GPy.kern.RBF(input_dim=inputdim, variance=1., lengthscale=[1. for iii in range(inputdim)],ARD=True)")
pyExecp("gp = GPy.models.GPRegression(X,y,kernel,normalizer=True)")
pyExecp("gp.likelihood.variance = 1e-8")
pyExecp("gp.optimize(messages=False)")
pyExecp("gp.optimize_restarts(num_restarts = 5,  verbose=False)")
python.assign("xp1", XX1)
python.exec("xp = np.asmatrix(xp1)")
python.exec("y_pred, sigma2_pred = gp.predict(np.asarray(xp))")
unlist(python.get("y_pred.tolist()"))

# update with new data

python.assign("X1", (Xall))
python.assign("y1", Zall)
python.exec('X =  np.matrix(X1)')
python.exec('y = np.matrix(y1).reshape((-1,1))')
#try(python.exec("gp.set_Y(Y = y)"))
#python.exec("gp.set_X(X = X)")
#python.exec("gp.set_Y(Y = y)")
python.exec("gp.set_XY(X = X, Y = y)")
python.exec("gp.optimize(messages=False)")
python.exec("gp.optimize_restarts(num_restarts = 5,  verbose=False)")

pyExit()


# Dice kriging test
mod <- DiceKriging::km(design = X1, response = Z1, covtype = "gauss", nugget.estim = T)
predict.km(object = mod, newdata = XX1, type = "SK")



# Test flat surface
set.seed(0)
n <- 8
d <- 1
n2 <- 2
f1 <- function(x) {0}#sin(2*pi*x[1]) + sin(2*pi*x[2])}
#f1 <- branin
X1 <- matrix(runif(n*d),n,d)
Z1 <- apply(X1,1,f1) #+ rnorm(n, 0, 1e-6)
X2 <- matrix(runif(n2*d),n2,d)
Z2 <- apply(X2,1,f1)
Xall <- rbind(X1, X2)
Zall <- c(Z1, Z2)
XX1 <- matrix(runif(10),5,2)
ZZ1 <- apply(XX1, 1, f1)
system.time(u <- UGP$new(package='laGP',X=X1,Z=Z1, corr.power=2))
cbind(u$predict(XX1), ZZ1)
u$predict.se(XX1)


# Test update function
u <- UGP$new(package="laGP")
n <- 3
d <- 2
f <- TestFunctions::gaussian1
for (i in 1:10) {
  x <- matrix(runif(n*d), ncol=d)
  y <- f(x)
  u$update(Xnew=x, Znew=y)
  cf::cf(u$predict, pts=u$X)
}


# Test function with all near zero values
set.seed(0)
n <- 3
d <- 2
n2 <- 5
f1 <- function(x) {rnorm(1,0,1e-16)}
f1 <- TestFunctions::banana
X1 <- matrix(runif(n*d),n,d)
Z1 <- apply(X1,1,f1) #+ rnorm(n, 0, 1e-3)
X2 <- matrix(runif(n2*d),n2,d)
Z2 <- apply(X2,1,f1)
XX1 <- matrix(runif(10),ncol=1)
ZZ1 <- apply(XX1, 1, f1)
system.time(u <- UGP$new(package='laGP',X=X1,Z=Z1, estimate.nugget=T))
cf::cf(u$predict)
u$update(Xnew=X2,Znew=Z2)
cf::cf(u$predict)
u$delete()





# Try Matlab???
library(R.matlab)

Matlab$startServer()
matlab <- Matlab()
isOpen <- open(matlab)
if (!isOpen) throw("MATLAB server is not running: waited 30 seconds.")

# set a variable in R and send to MATLB
x <- 10
setVariable(matlab, x = x)
evaluate(matlab, "x")
evaluate(matlab, "y=20; z=x+y")
z <- getVariable(matlab, "z")
z
close(matlab)

setVariable(matlab, X1=X1)
evaluate(matlab, "X1 .* X1", capture=TRUE)->temp
temp
