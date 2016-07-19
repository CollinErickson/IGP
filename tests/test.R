n <- 10
d <- 2
n2 <- 100
X1 <- matrix(runif(n*d),n,d)
Z1 <- sin(8*pi*X1[,1]) + sin(8*pi*X1[,2])
X2 <- matrix(runif(n2*d),n2,d)
Z2 <- sin(8*pi*X2[,1]) + sin(8*pi*X2[,2])
XX1 <- matrix(runif(10),5,2)
u <- UGP$new(package='laGP',X=X1,Z=Z1)
u$predict(XX1)
#apply(XX1,1,u$predict.se)#(XX1)
u$predict.se(XX1)
#u$predict.se(XX1[1,])
contourfilled::contourfilled.func(u$predict,batchmax = 100);points(X1)
u$update(Xnew=X2,Znew=Z2)
u$predict(XX1)
#u$update()
#laGP::updateGPsep(gpsepi=u$mod[[1]], X=u$X, Z=u$Z)
contourfilled::contourfilled.func(u$predict,batchmax = 100);points(rbind(X1,X2))
u$delete()
contourfilled::contourfilled.func(function(xx)sin(8*pi*xx[1]) + sin(8*pi*xx[2]))
contourfilled::contourfilled.data(X1,Z1)
contourfilled::contourfilled.data(X2,Z2)
contourfilled::contourfilled.data(rbind(X1,X2),c(Z1,Z2))
contourfilled::contourfilled.data(u$X,u$Z)


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
contourfilled::contourfilled.func(function(xx){predGPsep(gpi,XX=xx,lite=T)$mean},batchmax = 100)

## clean up
deleteGPsep(gpi)
