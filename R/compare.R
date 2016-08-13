compare.UGP <- function(packages, func, D, N, Npred=1000, reps=1, debug=F) {
  if (debug) browser()
  out <- data.frame()
  for (rep in 1:reps) {
    X <- lhs::maximinLHS(n=N, k=D)
    Y <- apply(X, 1, func)
    Xpred <- lhs::maximinLHS(n=Npred, k=D)
    Ypred <- apply(Xpred, 1, func)
    for (ipackage in seq_along(packages)) {
      package <- packages[ipackage]
      if(ipackage==1) fit.time <- system.time(u <- GauPro::GauPro$new(X=X, Z=Y, useOptim2=F))[3]
      else if(ipackage==2) fit.time <- system.time(u <- GauPro::GauPro$new(X=X, Z=Y, useOptim2=T))[3]
      else fit.time <- system.time(u <- UGP::UGP$new(X=X, Z=Y, package=package))[3]
      #fit.time <- system.time({
      #  u <- UGP::UGP$new(X=matrix(NA,ncol=0,nrow=0), Z=numeric(0), package=package)
      #  u$init(X=X, Z=Y)
      #})[3]
      predict.time <- system.time(up <- u$predict(Xpred, se.fit=T))[3]
      if (ipackage < 3) mse <- mean((up$me - Ypred)^2)
      else mse <- mean((up$fit - Ypred)^2)
      pmse <- mean((up$se)^2)
      rmse <- sqrt(mse)
      prmse <- sqrt(pmse)

      out.new <- data.frame(package=package, rep=rep, fit.time=fit.time, predict.time=predict.time,
                            mse=mse, pmse=pmse, rmse=rmse, prmse=prmse)
      out <- rbind(out, out.new)

      if(ipackage>3)u$delete()
    }
  }
  stripchart(rmse ~ package, data=out)
  out
}
