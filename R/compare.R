compare.UGP <- function(packages, func, D, N, Npred=1000, reps=1, debug=F, init_list=list(), pred_list=list()) {
  if (debug) browser()
  out <- data.frame()
  for (rep in 1:reps) {
    X <- lhs::maximinLHS(n=N, k=D)
    Y <- apply(X, 1, func)
    Xpred <- lhs::maximinLHS(n=Npred, k=D)
    Ypred <- apply(Xpred, 1, func)
    for (ipackage in seq_along(packages)) {
      package <- packages[ipackage]
      fit.time <- system.time({
        u <- do.call(UGP::UGP$new,
                     c(list(X=X, Z=Y, package=package),
                       if (!is.null(init_list[[as.character(ipackage)]])) init_list[[as.character(ipackage)]]))
      })[3]

      predict.time <- system.time(up <- do.call(u$predict, list(Xpred, se.fit=T)))[3]

      mse <- mean((up$fit - Ypred)^2)
      pmse <- mean((up$se)^2)
      rmse <- sqrt(mse)
      prmse <- sqrt(pmse)

      out.new <- data.frame(package=package, rep=rep, fit.time=fit.time, predict.time=predict.time,
                            mse=mse, pmse=pmse, rmse=rmse, prmse=prmse)
      out <- rbind(out, out.new)

      u$delete()
    }
  }
  stripchart(rmse ~ package, data=out)
  out
}
