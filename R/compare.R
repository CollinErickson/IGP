#' Compare multiple UGPs on both prediction accuracy and run time
#'
#' @param packages vector of packages to compare
#' @param func The function to test on
#' @param D Number of dimensions for the data
#' @param N Number of data points
#' @param Npred Number of prediction points
#' @param reps Number of replicates
#' @param debug If true will start browser
#' @param init_list Additional params to pass to init
#' @param pred_list Additional params to pass to pred
#'
#' @return Dataframe of results
#' @export
#' @importFrom graphics stripchart
#' @examples
#' compare.UGP(packages=c("laGP", "mlegp"), func=function(x){2*pi*x}, D=1, N=10, Npred = 1e3, reps = 3)
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
      package.use <- strsplit(package, '-')[[1]][1] # lets you add identifying name after a hyphen
      fit.time <- system.time({
        u <- do.call(UGP::UGP$new,
                     c(list(X=X, Z=Y, package=package.use),
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
  out$total.time <- out$fit.time + out$predict.time

  stripchart(rmse ~ package, data=out)
  #browser()
  com2 <- reshape::melt.data.frame(out, measure.var=c('rmse', 'prmse'), id.vars=c('package','rep'), variable_name="rmseprmse")

  #rmseprmse_plot <- (ggplot2::ggplot(com2, ggplot2::aes(x=value, y=rmseprmse, color=as.factor(rep)))
  #                   + ggplot2::geom_point(ggplot2::aes(shape=rmseprmse),size=3) + ggplot2::facet_grid(package ~ .)
  #browser()
  rmseprmse_plot <- (ggplot2::ggplot(com2, ggplot2::aes_string(x="value", y="rmseprmse", color="package", shape="as.factor(rep)"))
                     + ggplot2::geom_point(ggplot2::aes_string(shape="as.factor(rep)"),size=3) + ggplot2::facet_grid(package ~ .)
                     + ggplot2::guides(shape=F,color=F)
                     + ggplot2::ylab(NULL) + ggplot2::xlab(NULL)
                     +  ggplot2::geom_line(ggplot2::aes_string(x="value", y="rmseprmse", group="rep")))

  #time_plot <- (ggplot2::ggplot(out, ggplot2::aes(x=total.time,y=package, color=as.factor(rep)))
  time_plot <- (ggplot2::ggplot(out, ggplot2::aes_(x=~total.time,y=~package, color=~as.factor(rep)))
                + ggplot2::geom_point(size=3)
                + ggplot2::ylab(NULL) + ggplot2::xlab("Run time (s)")
                + ggplot2::guides(color=FALSE))
  print(time_plot)
  print(rmseprmse_plot)

  out
}
