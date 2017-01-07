library(R6)
GPcompare <- R6::R6Class(classname="GPcompare",
  public=list(
    D=NULL,
    reps=NULL,
    input.ss=NULL,
    test.ss=NULL,
    func=NULL,
    Xs=NULL,
    Zs=NULL,
    initialize = function(D, reps, input.ss, test.ss, func, packages) {

      self$D <- D
      self$reps <- reps
      self$input.ss <- input.ss
      self$test.ss <- test.ss
      self$func <- func
      self$xs <- list()
      self$ys <- list()
      if (is.character(packages)) {
        self$packages <- packages
        self$package_names <- packages
      } else { # list input, first of each is package, second is name, after that is options

        self$packages <- sapply(packages, function(xx) xx[1])
        self$package_names <- sapply(packages, function(xx) if (length(xx)>1) xx[2] else xx[1])
      }
    },
    create_data = function() {

      for (rep in 1:self$reps) {
        Xnew <- lhs::maximinLHS(self$input.ss, self$D)
        Znew <- apply(Xnew, 1, self$func)
        self$Xs <- c(self$Xs, Xnew)
        self$Zs <- c(self$Zs, Znew)
        self$outputlist <- list()
      }
    },
    run_fits = function() {
      for (rep in 1:reps) {
        for (ipackage in seq_along(self$packages)) {
          package <- self$packages[ipackage]
          package_name <- self$package_names[ipackage]
          #package.use <- strsplit(package, '-')[[1]][1] # lets you add identifying name after a hyphen
          fit.time <- system.time({
            u <- do.call(UGP::UGP$new,
                         c(list(X=self$Xs[[ipackage]], Z=self$Zs[[ipackage]], package=package),
                           if (!is.null(init_list[[as.character(ipackage)]])) init_list[[as.character(ipackage)]]))
          })[3]

          predict.time <- system.time(up <- do.call(u$predict, list(Xpred, se.fit=T)))[3]

          mse <- mean((up$fit - Ypred)^2)
          pmse <- mean((up$se)^2)
          rmse <- sqrt(mse)
          prmse <- sqrt(pmse)

          out.new <- data.frame(package=package, rep=rep, fit.time=fit.time, predict.time=predict.time,
                                mse=mse, pmse=pmse, rmse=rmse, prmse=prmse)
          #out <- rbind(out, out.new)
          self$outputlist[[package_name]][[rep]] <- out.new
          u$delete()
        }
      }
    },
    plot_output = function() {

      stripchart(rmse ~ package, data=out)
    }
  )
) # End GPcompare R6 class
