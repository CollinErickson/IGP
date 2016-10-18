packs <- c('mlegp', 'GPfit', 'laGP', 'GauPro', 'DiceKriging')[-c(4)]
#i <- i+1;set.seed(i)
set.seed(13)
x <- matrix(lhs::maximinLHS(18,1), ncol=1)
xp <- matrix(seq(0,1,len=500), ncol=1)
f <- function(xx)abs(sin(xx^.9*2*pi))^.9
f <- function(xx)ifelse(xx<.5, abs(sin(xx^.9*2*pi))^.9, 1*abs(sin(xx^.9*2*pi))^.7)
f <- function(xx)abs(ifelse(xx<.5, abs(sin(xx^.9*2*pi))^.9, 1*abs(sin(xx^.9*2*pi))^.7) -.5)# n=18, set.seed= 3, 8, 13, 20, 31
#f <- function(xx)xx*(sin(xx*6*pi))
#f <- function(xx) pmin(abs(xx-1/3) , 2*abs(xx-2/3))
#f <- function(xx) (2 * xx) %% 1
#f <- function(xx) log((2 * xx) %% 1 +.1)
#f <- function(xx)xx*(sin(2*pi/(xx+.3))) + xx*(sin(2*pi/(1-xx^.77+.3)))
y <- f(x)
yp <- f(xp)

single.run <- function(pac) {
  u <- UGP::UGP$new(X=x,Z=y,package=pac)
  to <- u$predict(xp)
  u$delete()
  to
}

out <- lapply(packs, single.run)
packs <- packs[-2]
out <- out[-2]

plot(xp,yp, type='l', xlab='x', ylab='y')
legend(x='bottomright', legend=packs, fill=1+1:length(packs))
lapply(1:length(out),function(outi){outj <- out[[outi]];points(xp,outj, type='l', col=outi+1, lwd=6-outi)})
points(x,y,pch=19, cex=2)
