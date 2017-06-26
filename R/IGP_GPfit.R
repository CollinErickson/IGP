#' UGP
#' Class providing object with methods for fitting a GP model
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords data, kriging, Gaussian process, regression
#' @return Object of \code{\link{R6Class}} with methods for fitting GP model.
#' @format \code{\link{R6Class}} object.
#' @examples
#' n <- 20
#' d <- 2
#' n2 <- 20
#' f1 <- function(x) {sin(2*pi*x[1]) + sin(2*pi*x[2])}
#' X1 <- matrix(runif(n*d),n,d)
#' Z1 <- apply(X1,1,f1) + rnorm(n, 0, 1e-3)
#' X2 <- matrix(runif(n2*d),n2,d)
#' Z2 <- apply(X2,1,f1)
#' XX1 <- matrix(runif(10),5,2)
#' ZZ1 <- apply(XX1, 1, f1)
#' u <- IGP(package='laGP',X=X1,Z=Z1)
#' cbind(u$predict(XX1), ZZ1)
#' u$predict.se(XX1)
#' u$update(Xnew=X2,Znew=Z2)
#' u$predict(XX1)
#' u$delete()
#' @field X Design matrix
#' @field Z Responses
#' @field N Number of data points
#' @field D Dimension of data
#' @section Methods:
#' \describe{
#'   \item{Documentation}{For full documentation of each method go to https://github.com/CollinErickson/UGP/}
#'   \item{\code{new(X=NULL, Z=NULL, package=NULL,
#'   estimate.nugget=T, set.nugget=F, ...)}}{This method
#'   is used to create object of this class with \code{X} and \code{Z} as the data.
#'   The package tells it which package to fit the GP model.}
#'   \item{\code{Xall=NULL, Zall=NULL, Xnew=NULL, Znew=NULL, ...}}{This method
#'   updates the model, adding new data if given, then running optimization again.}}
IGP_GPfit <- R6::R6Class(classname = "IGP_GPfit", inherit = IGP_base,
                          public = list(
                            newthing = 33,
                            .init = function(...) {
                              if (!is.null(self$estimate.nugget) || self$set.nugget) {
                                warning("GPfit cannot estimate or set the nugget, it picks a stable value")
                              }
                              if (self$corr[[1]] == "gauss") {
                                self$mod <- GPfit::GP_fit(self$X, self$Z, corr = list(type="exponential",power=2))
                              } else if (self$corr[[1]] == "powerexp") {
                                self$mod <- GPfit::GP_fit(self$X, self$Z, corr = list(type="exponential",power=self$corr.power))
                              } else if (self$corr[[1]] == "matern") {
                                self$mod <- GPfit::GP_fit(self$X, self$Z, corr = list(type="matern",nu=self$corr[[2]]))
                              } else if (self$corr[[1]] == "matern5_2") {
                                self$mod <- GPfit::GP_fit(self$X, self$Z, corr = list(type="matern",nu=5/2))
                              } else if (self$corr[[1]] == "matern3_2") {
                                self$mod <- GPfit::GP_fit(self$X, self$Z, corr = list(type="matern",nu=3/2))
                              } else {
                                stop("corr not recognized for GPfit")
                              }
                            },
                            .update = function(...){
                              self$.init()
                            }, #"function",
                            .predict = function(XX, se.fit, ...){#browser()
                              if (se.fit) {
                                preds <- GPfit::predict.GP(self$mod, XX, se.fit=se.fit)
                                list(fit=preds$Y_hat, se.fit=sqrt(preds$MSE))
                              } else {
                                GPfit::predict.GP(self$mod, XX)$Y_hat
                              }
                            }, #"function",
                            .predict.se = function(XX, ...) {sqrt(GPfit::predict.GP(object=self$mod, xnew=XX, se.fit=T)$MSE)}, #"function",
                            .predict.var = function(XX, ...) {GPfit::predict.GP(object=self$mod, xnew=XX, se.fit=T)$MSE}, #"function",
                            #.grad = NULL,
                            .delete = function(...){self$mod <- NULL}, #"function",
                            .theta = function() {10^(self$mod$beta)}, #"function",
                            .nugget = function() {self$mod$delta} #"function",
                            #.mean = NULL # function that gives mean

                          )
)

#' UGP
#' Class providing object with methods for fitting a GP model
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords data, kriging, Gaussian process, regression
#' @return Object of \code{\link{R6Class}} with methods for fitting GP model.
#' @format \code{\link{R6Class}} object.
#' @examples
#' n <- 40
#' d <- 2
#' n2 <- 20
#' f1 <- function(x) {sin(2*pi*x[1]) + sin(2*pi*x[2])}
#' X1 <- matrix(runif(n*d),n,d)
#' Z1 <- apply(X1,1,f1) + rnorm(n, 0, 1e-3)
#' X2 <- matrix(runif(n2*d),n2,d)
#' Z2 <- apply(X2,1,f1)
#' XX1 <- matrix(runif(10),5,2)
#' ZZ1 <- apply(XX1, 1, f1)
#' u <- IGP_laGP$new(X=X1,Z=Z1)
#' cbind(u$predict(XX1), ZZ1)
#' u$predict.se(XX1)
#' u$update(Xnew=X2,Znew=Z2)
#' u$predict(XX1)
#' u$delete()
#' @field X Design matrix
#' @field Z Responses
#' @field N Number of data points
#' @field D Dimension of data
#' @section Methods:
#' \describe{
#'   \item{Documentation}{For full documentation of each method go to https://github.com/CollinErickson/UGP/}
#'   \item{\code{new(X=NULL, Z=NULL, package=NULL,
#'   estimate.nugget=T, set.nugget=F, ...)}}{This method
#'   is used to create object of this class with \code{X} and \code{Z} as the data.
#'   The package tells it which package to fit the GP model.}
#'   \item{\code{Xall=NULL, Zall=NULL, Xnew=NULL, Znew=NULL, ...}}{This method
#'   updates the model, adding new data if given, then running optimization again.}}
IGP_laGP <- R6::R6Class(classname = "IGP_laGP", inherit = IGP_base,
                            public = list(
                              .init = function(..., d=NULL, g=NULL, theta=NULL, nugget=NULL, no_update=FALSE) {#browser()
                                if (self$corr[[1]] != "gauss") {
                                  stop("laGP only uses Gaussian correlation")
                                }

                                if (is.null(d) & !is.null(theta)) {d <- 1/theta}
                                if (is.null(g) && is.null(nugget) && !is.null(self$set.nugget)) {g <- self$set.nugget}
                                if (is.null(g) & !is.null(nugget)) {g <- nugget}

                                da <- laGP::darg(list(mle=TRUE), X=self$X)
                                ga.try <- try(ga <- laGP::garg(list(mle=TRUE), y=self$Z), silent = T)
                                if (inherits(ga.try, "try-error")) {
                                  warning("Adding noise to ga in laGP");
                                  ga <- laGP::garg(list(mle=TRUE), y=self$Z+rnorm(length(self$Z),0,1e-2))
                                }

                                # Follow recommendations for small samples, otherwise use bigger range
                                drange <- if (nrow(self$X)<20) c(da$min, da$max) else c(1e-3,1e4) #c(da$min, da$max), # Don't like these small ranges
                                grange <- c(ga$min, ga$max)
                                da_start <- if (!is.null(d)) d else da$start
                                ga_start <- if (!is.null(g)) g else ga$start
                                mod1 <- laGP::newGPsep(X=self$X, Z=self$Z, d=da_start, g=ga$start, dK = TRUE)
                                #mod1 <- laGP::newGPsep(X=X, Z=Z, d=da$start, g=1e-6, dK = TRUE)
                                if (no_update) { #using d and g given
                                  self$mod.extra$theta = as.numeric(1 / d) # store theta params
                                  self$mod.extra$nugget = as.numeric(g) # store nugget
                                } else if (!no_update && !self$estimate.nugget) { # Only estimate d/theta
                                  mle.out <- laGP::mleGPsep(gpsepi = mod1,
                                                            param="d",
                                                             tmin=drange[1], tmax=drange[2],
                                                             verb=0, maxit=1000)
                                  self$mod.extra$theta <- as.numeric(1 / mle.out$d) # store theta params
                                  self$mod.extra$nugget <- g
                                  # Leave nugget as is
                                } else if (!no_update) { # Update both
                                  mle.out <- laGP::jmleGPsep(gpsepi = mod1,
                                                             drange=drange,
                                                             grange=grange,
                                                             #dab=da$ab, gab=ga$ab, # Will use MLE without these
                                                             verb=0, maxit=1000)
                                  self$mod.extra$theta = as.numeric(1 / mle.out[1,1:ncol(self$X)]) # store theta params
                                  self$mod.extra$nugget = as.numeric(mle.out[1,ncol(self$X) + 1]) # store nugget
                                } else {stop("Shouldn't be here IGP_laGP #32097555")}
                                self$mod <- mod1
                              }, #"function to initialize model with data
                              .update = function(..., no_update=FALSE) {#browser()
                                if (no_update) { # just add data and return
                                  laGP::updateGPsep(gpsepi=self$mod,
                                                    X=self$X[-(1:self$n.at.last.update), , drop=FALSE],
                                                    Z=self$Z[-(1:self$n.at.last.update)])
                                  return()
                                }

                                # Start over if not many points, had problems getting stuck in bad spots early
                                if (self$n.at.last.update < 20) {
                                  self$.delete()
                                  self$.init(...)
                                  return()
                                }

                                da <- laGP::darg(list(mle=TRUE), X=self$X)
                                # Had problem with garg that didn't make sense, trying to fix by adding noise to Z
                                ga.try <- try(ga <- laGP::garg(list(mle=TRUE), y=self$Z), silent = TRUE)
                                if (inherits(ga.try, "try-error")) {
                                  print("Adding noise to Z so garg doesn't give error")
                                  noise.sd <- 1e-8
                                  while (TRUE) {
                                    ga.try <- try(ga <- laGP::garg(list(mle=TRUE), y=self$Z + rnorm(length(self$Z), 0, noise.sd)), silent = TRUE)
                                    if (!inherits(ga.try, "try-error")) {break}
                                    noise.sd <- 10 * noise.sd
                                  }
                                }
                                n.since.last.update <- nrow(self$X) - self$n.at.last.update
                                if (n.since.last.update < 1) {
                                  warning("Can't update, no new X rows, but can optimize again")
                                } else {
                                  if (self$n.at.last.update < 10 || n.since.last.update > .25 * self$n.at.last.update) {
                                    # start over if too many
                                    self$.delete(...=...)
                                    self$.init(...=...)
                                  } else {
                                    lagpupdate.try <- try(
                                      laGP::updateGPsep(gpsepi=self$mod,
                                                        X=self$X[-(1:self$n.at.last.update), , drop=FALSE],
                                                        Z=self$Z[-(1:self$n.at.last.update)])
                                    )
                                    if (inherits(lagpupdate.try, "try-error")) {browser()}
                                  }
                                }
                                drange <- c(1e-3,1e4)
                                grange <- c(min(sqrt(.Machine$double.eps),self$mod.extra$nugget), max(1,self$mod.extra$nugget))
                                if (no_update) {
                                  stop("This is covered above at beginning of .update")
                                } else if (!no_update && !self$estimate.nugget) { # update d/theta but not nugget/g
                                  mle.try <- try(mle.out <- laGP::mleGPsep(gpsepi = self$mod,
                                                                           param = "d",
                                                                            tmin=drange[1], tmax=drange[2],
                                                                            verb=0, maxit=1000))
                                  if (inherits(mle.try, "try-error")) {
                                    # Sometimes gives error: L-BFGS-B needs finite values of 'fn'
                                    browser()
                                    warning('Restarting laGP model')
                                    self$delete()
                                    self$init(..., no_update=no_update)
                                    return()
                                  }
                                  # Update stored parameters for when user calls $theta() or $nugget()
                                  self$mod.extra$theta = as.numeric(1 / mle.out$d) # store theta params
                                  # leave nugget as it was
                                } else if (!no_update && self$estimate.nugget) {
                                  mle.try <- try(mle.out <- laGP::jmleGPsep(gpsepi = self$mod,
                                                                            #drange=c(da$min, da$max), # Getting rid of these here too
                                                                            #grange=c(ga$min, ga$max),
                                                                            drange=drange,
                                                                            grange=grange, # Had error of nugget starting outside bound
                                                                            #dab=da$ab, gab=ga$ab,
                                                                            verb=0, maxit=1000))
                                  if (inherits(mle.try, "try-error")) {
                                    # Sometimes gives error: L-BFGS-B needs finite values of 'fn'
                                    browser()
                                    warning('Restarting laGP model')
                                    self$delete()
                                    self$init(...)
                                    return()
                                  }
                                  # Update stored parameters for when user calls $theta() or $nugget()
                                  self$mod.extra$theta = as.numeric(1 / mle.out[1,1:ncol(self$X)]) # store theta params
                                  self$mod.extra$nugget = as.numeric(mle.out[1,ncol(self$X) + 1]) # store nugget
                                }
                              }, #"function to add data to model or reestimate params
                              .predict = function(XX, se.fit, ...){
                                if (se.fit) {
                                  preds <- laGP::predGPsep(self$mod, XX, lite=TRUE)
                                  list(fit=preds$mean, se.fit=sqrt(preds$s2))
                                } else {
                                  laGP::predGPsep(self$mod, XX, lite=TRUE)$mean
                                }
                              }, #"function to predict at new values
                              .predict.se = function(XX, ...) {sqrt(laGP::predGPsep(self$mod, XX, lite=TRUE)$s2)}, #"function predict the standard error/dev
                              .predict.var = function(XX, ...) {laGP::predGPsep(self$mod, XX, lite=TRUE)$s2}, #"function to predict the variance
                              .grad = NULL, # function to calculate the gradient
                              .delete = function(...) {
                                if (!is.null(self$mod)) {
                                  laGP::deleteGPsep(self$mod)
                                  self$mod <- NULL
                                }
                              }, #"function to delete model beyond simple deletion
                              .theta = function() {self$mod.extra$theta}, #"function to get theta, exp(-theta*(x-x)^2)
                              .nugget = function() {self$mod.extra$nugget}, #"function to get nugget
                              .mean = NULL # function that gives mean (constant, other functions not implemented)

                            )
)

#' UGP
#' Class providing object with methods for fitting a GP model
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords data, kriging, Gaussian process, regression
#' @return Object of \code{\link{R6Class}} with methods for fitting GP model.
#' @format \code{\link{R6Class}} object.
#' @examples
#' \dontrun{
#' n <- 40
#' d <- 2
#' n2 <- 20
#' f1 <- function(x) {sin(2*pi*x[1]) + sin(2*pi*x[2])}
#' X1 <- matrix(runif(n*d),n,d)
#' Z1 <- apply(X1,1,f1) + rnorm(n, 0, 1e-3)
#' X2 <- matrix(runif(n2*d),n2,d)
#' Z2 <- apply(X2,1,f1)
#' XX1 <- matrix(runif(10),5,2)
#' ZZ1 <- apply(XX1, 1, f1)
#' u <- IGP_tgp$new(X=X1,Z=Z1)
#' cbind(u$predict(XX1), ZZ1)
#' u$predict.se(XX1)
#' u$update(Xnew=X2,Znew=Z2)
#' u$predict(XX1)
#' u$delete()
#' }
#' @field X Design matrix
#' @field Z Responses
#' @field N Number of data points
#' @field D Dimension of data
#' @section Methods:
#' \describe{
#'   \item{Documentation}{For full documentation of each method go to https://github.com/CollinErickson/UGP/}
#'   \item{\code{new(X=NULL, Z=NULL, package=NULL,
#'   estimate.nugget=T, set.nugget=F, ...)}}{This method
#'   is used to create object of this class with \code{X} and \code{Z} as the data.
#'   The package tells it which package to fit the GP model.}
#'   \item{\code{Xall=NULL, Zall=NULL, Xnew=NULL, Znew=NULL, ...}}{This method
#'   updates the model, adding new data if given, then running optimization again.}}
IGP_tgp <- R6::R6Class(classname = "IGP_tgp", inherit = IGP_base,
                            public = list(
                              .init = function(package=self$package, ...) {#browser()
                                if (self$corr[[1]] != "gauss") {
                                  stop("tgp only uses Gaussian correlation")
                                }
                                modfunc <-  if (self$package == "blm") tgp::blm
                                else if (package == "btlm") tgp::btlm
                                else if (package == "bcart") tgp::bcart
                                else if (package == "bgp") tgp::bgp
                                else if (package == "bgpllm") tgp::bgpllm
                                else if (package == "btgp") tgp::btgp
                                else if (package == "btgpllm") tgp::btgpllm
                                capture.output(mod1 <- modfunc(self$X, self$Z))
                                self$mod <- mod1
                              }, #"function to initialize model with data
                              .update = function(...) {#browser()
                                self$.init(...=...)
                              }, #"function to add data to model or reestimate params
                              .predict = function(XX, se.fit, ...){#browser()
                                Xsplit0 <- self$mod$Xsplit # Need to add this Xsplit stuff so it can predict outside the box of original data
                                self$mod$Xsplit <- rbind(as.matrix(Xsplit0), XX) # This is the workaround mentioned on p37 https://cran.r-project.org/web/packages/tgp/tgp.pdf, have to as.matrix to avoid col name error
                                capture.output(preds <- with(globalenv(), predict)(self$mod, XX))
                                self$mod$Xsplit <- Xsplit0 # Reset value
                                if (se.fit) {
                                  list(fit=preds$ZZ.km, se.fit=sqrt(preds$ZZ.ks2))
                                } else {
                                  preds$ZZ.km
                                }
                              }, #"function to predict at new values
                              .predict.se = function(XX, ...) {
                                Xsplit0 <- self$mod$Xsplit # Need to add this Xsplit stuff so it can predict outside the box of original data
                                self$mod$Xsplit <- rbind(as.matrix(Xsplit0), XX) # This is the workaround mentioned on p37 https://cran.r-project.org/web/packages/tgp/tgp.pdf, have to as.matrix to avoid col name error
                                toreturn <- sqrt(with(globalenv(), predict)(self$mod, XX)$ZZ.ks2)
                                self$mod$Xsplit <- Xsplit0 # Reset value
                                toreturn
                              }, #"function predict the standard error/dev
                              .predict.var = function(XX, ...) {
                                Xsplit0 <- self$mod$Xsplit # Need to add this Xsplit stuff so it can predict outside the box of original data
                                self$mod$Xsplit <- rbind(as.matrix(Xsplit0), XX) # This is the workaround mentioned on p37 https://cran.r-project.org/web/packages/tgp/tgp.pdf, have to as.matrix to avoid col name error
                                toreturn <- with(globalenv(), predict)(self$mod, XX)$ZZ.ks2
                                self$mod$Xsplit <- Xsplit0 # Reset value
                                toreturn
                              }, #"function to predict the variance
                              .grad = NULL, # function to calculate the gradient
                              .delete = function(...) {self$mod <- NULL}, #"function to delete model beyond simple deletion
                              .theta = function() {rep(NA, ncol(self$X))}, #"function to get theta, exp(-theta*(x-x)^2)
                              .nugget = function() {NA}, #"function to get nugget
                              .mean = NULL # function that gives mean (constant, other functions not implemented)

                            )
)

#' UGP
#' Class providing object with methods for fitting a GP model
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords data, kriging, Gaussian process, regression
#' @return Object of \code{\link{R6Class}} with methods for fitting GP model.
#' @format \code{\link{R6Class}} object.
#' @examples
#' n <- 40
#' d <- 2
#' n2 <- 20
#' f1 <- function(x) {sin(2*pi*x[1]) + sin(2*pi*x[2])}
#' X1 <- matrix(runif(n*d),n,d)
#' Z1 <- apply(X1,1,f1) + rnorm(n, 0, 1e-3)
#' X2 <- matrix(runif(n2*d),n2,d)
#' Z2 <- apply(X2,1,f1)
#' XX1 <- matrix(runif(10),5,2)
#' ZZ1 <- apply(XX1, 1, f1)
#' u <- IGP_mlegp$new(X=X1,Z=Z1)
#' cbind(u$predict(XX1), ZZ1)
#' u$predict.se(XX1)
#' u$update(Xnew=X2,Znew=Z2)
#' u$predict(XX1)
#' u$delete()
#' @field X Design matrix
#' @field Z Responses
#' @field N Number of data points
#' @field D Dimension of data
#' @section Methods:
#' \describe{
#'   \item{Documentation}{For full documentation of each method go to https://github.com/CollinErickson/UGP/}
#'   \item{\code{new(X=NULL, Z=NULL, package=NULL,
#'   estimate.nugget=T, set.nugget=F, ...)}}{This method
#'   is used to create object of this class with \code{X} and \code{Z} as the data.
#'   The package tells it which package to fit the GP model.}
#'   \item{\code{Xall=NULL, Zall=NULL, Xnew=NULL, Znew=NULL, ...}}{This method
#'   updates the model, adding new data if given, then running optimization again.}}
IGP_mlegp <- R6::R6Class(classname = "IGP_mlegp", inherit = IGP_base,
                            public = list(
                              .init = function(...) {
                                if (self$corr[[1]] != "gauss") {
                                  stop("mlegp only uses Gaussian correlation")
                                }
                                # Have to do this for nugget or else it won't work
                                temp_nug <- if (is.null(self$estimate.nugget) || self$estimate.nugget == FALSE) NULL
                                            else if (self$estimate.nugget == TRUE) 1e-6
                                temp_nug_known <- if (is.null(self$set.nugget)) 0 else self$set.nugget

                                # Having trouble when nugget passed as ... since two match, so trying other method
                                if (FALSE) {
                                  # Use these in case they are passed as ...
                                  if (!exists("nugget")) temp_nug <- NULL
                                  if (!exists("nugget.known")) temp_nug_known <- NULL
                                  co <- capture.output(m <- mlegp::mlegp(X=self$X, Z=self$Z, verbose=0,
                                                                         nugget = nugget, #temp_nug,
                                                                         nugget.known = nugget.known, #temp_nug_known,
                                                                         #nugget = self$set.nugget,
                                                                         #nugget.known=as.integer(!self$estimate.nugget),
                                                                         ...)
                                                       )
                                } else { # use args and do.call to avoid multiple matching
                                  #browser()
                                  argss <- list(...)
                                  argss$X <- self$X
                                  argss$Z <- self$Z
                                  argss$verbose <- 0
                                  if (is.null(argss$nugget)) argss$nugget <- temp_nug
                                  if (is.null(argss$nugget.known)) argss$nugget.known <- temp_nug_known
                                  co <- capture.output(m <- do.call(mlegp::mlegp, argss))
                                }

                                self$mod <- m
                              }, #"function to initialize model with data
                              .update = function(...) {
                                self$.init(...)
                              }, #"function to add data to model or reestimate params
                              .predict = function(XX, se.fit, ...) {
                                mlegp::predict.gp(object=self$mod, newData=XX, se.fit = se.fit)
                              }, #"function to predict at new values
                              .predict.se = function(XX, ...) {mlegp::predict.gp(object=self$mod, newData=XX, se.fit=T)$se.fit}, #"function predict the standard error/dev
                              .predict.var = function(XX, ...) {mlegp::predict.gp(object=self$mod, newData=XX, se.fit=T)$se.fit^2}, #"function to predict the variance
                              .grad = NULL, # function to calculate the gradient
                              .delete = function(...){self$mod <- NULL}, #"function to delete model beyond simple deletion
                              .theta = function() {self$mod$beta}, #"function to get theta, exp(-theta*(x-x)^2)
                              .nugget = function() {self$mod$nugget}, #"function to get nugget
                              .mean = NULL # function that gives mean (constant, other functions not implemented)

                            )
)

#' UGP
#' Class providing object with methods for fitting a GP model
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords data, kriging, Gaussian process, regression
#' @return Object of \code{\link{R6Class}} with methods for fitting GP model.
#' @format \code{\link{R6Class}} object.
#' @examples
#' n <- 40
#' d <- 2
#' n2 <- 20
#' f1 <- function(x) {sin(2*pi*x[1]) + sin(2*pi*x[2])}
#' X1 <- matrix(runif(n*d),n,d)
#' Z1 <- apply(X1,1,f1) + rnorm(n, 0, 1e-3)
#' X2 <- matrix(runif(n2*d),n2,d)
#' Z2 <- apply(X2,1,f1)
#' XX1 <- matrix(runif(10),5,2)
#' ZZ1 <- apply(XX1, 1, f1)
#' u <- IGP_GauPro$new(X=X1,Z=Z1, parallel=FALSE)
#' cbind(u$predict(XX1), ZZ1)
#' u$predict.se(XX1)
#' u$update(Xnew=X2,Znew=Z2)
#' u$predict(XX1)
#' u$delete()
#' @field X Design matrix
#' @field Z Responses
#' @field N Number of data points
#' @field D Dimension of data
#' @section Methods:
#' \describe{
#'   \item{Documentation}{For full documentation of each method go to https://github.com/CollinErickson/UGP/}
#'   \item{\code{new(X=NULL, Z=NULL, package=NULL,
#'   estimate.nugget=T, set.nugget=F, ...)}}{This method
#'   is used to create object of this class with \code{X} and \code{Z} as the data.
#'   The package tells it which package to fit the GP model.}
#'   \item{\code{Xall=NULL, Zall=NULL, Xnew=NULL, Znew=NULL, ...}}{This method
#'   updates the model, adding new data if given, then running optimization again.}}
IGP_GauPro <- R6::R6Class(classname = "IGP_GauPro", inherit = IGP_base,
                            public = list(
                              .init = function(...) {#browser()
                                if (self$corr[[1]] != "gauss") {
                                  stop("GauPro only uses Gaussian correlation")
                                }
                                #m <- GauPro::GauPro$new(X=self$X, Z=self$Z, ...)
                                #m <- GauPro::GauPr_Gauss_par$new(X=self$X, Z=self$Z, ...)
                                m <- GauPro::GauPro(X=self$X, Z=self$Z, nug.est=self$estimate.nugget, nug=self$set.nugget, ...)
                                self$mod <- m
                              }, #"function to initialize model with data
                              .update = function(...) {
                                self$mod$update(Xall=self$X, Zall=self$Z, ...)
                              }, #"function to add data to model or reestimate params
                              .predict = function(XX, se.fit, ...) {
                                if (se.fit) {
                                  preds <- self$mod$pred(XX=XX, se.fit=T)
                                  list(fit=preds$mean, se.fit=preds$se)
                                } else {
                                  self$mod$pred(XX=XX)
                                }
                              }, #"function to predict at new values
                              .predict.se = function(XX, ...) {self$mod$pred(XX=XX, se.fit=T)$se}, #"function predict the standard error/dev
                              .predict.var = function(XX, ...) {self$mod$pred(XX=XX, se.fit=T)$s2}, #"function to predict the variance
                              .grad = function(XX) {self$mod$grad(XX=XX)}, # function to calculate the gradient
                              .delete = function(...){self$mod <- NULL}, #"function to delete model beyond simple deletion
                              .theta = function() {self$mod$theta}, #"function to get theta, exp(-theta*(x-x)^2)
                              .nugget = function() {self$mod$nug}, #"function to get nugget
                              .mean = NULL # function that gives mean (constant, other functions not implemented)

                            )
)

#' UGP
#' Class providing object with methods for fitting a GP model
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords data, kriging, Gaussian process, regression
#' @return Object of \code{\link{R6Class}} with methods for fitting GP model.
#' @format \code{\link{R6Class}} object.
#' @examples
#' n <- 40
#' d <- 2
#' n2 <- 20
#' f1 <- function(x) {sin(2*pi*x[1]) + sin(2*pi*x[2])}
#' X1 <- matrix(runif(n*d),n,d)
#' Z1 <- apply(X1,1,f1) + rnorm(n, 0, 1e-3)
#' X2 <- matrix(runif(n2*d),n2,d)
#' Z2 <- apply(X2,1,f1)
#' XX1 <- matrix(runif(10),5,2)
#' ZZ1 <- apply(XX1, 1, f1)
#' u <- IGP_DiceKriging$new(X=X1,Z=Z1)
#' cbind(u$predict(XX1), ZZ1)
#' u$predict.se(XX1)
#' u$update(Xnew=X2,Znew=Z2)
#' u$predict(XX1)
#' u$delete()
#' @field X Design matrix
#' @field Z Responses
#' @field N Number of data points
#' @field D Dimension of data
#' @section Methods:
#' \describe{
#'   \item{Documentation}{For full documentation of each method go to https://github.com/CollinErickson/UGP/}
#'   \item{\code{new(X=NULL, Z=NULL, package=NULL,
#'   estimate.nugget=T, set.nugget=F, ...)}}{This method
#'   is used to create object of this class with \code{X} and \code{Z} as the data.
#'   The package tells it which package to fit the GP model.}
#'   \item{\code{Xall=NULL, Zall=NULL, Xnew=NULL, Znew=NULL, ...}}{This method
#'   updates the model, adding new data if given, then running optimization again.}}
IGP_DiceKriging <- R6::R6Class(classname = "IGP_DiceKriging", inherit = IGP_base,
                            public = list(
                              .init = function(...) {
                                if (!exists("covtype")) {
                                  if (self$corr[[1]] %in% c("gauss", "powexp", "exp", "matern5_2", "matern3_2")) {
                                    covtype <- self$corr[[1]]
                                  } else if (self$corr[[1]] == "powerexp") {
                                    covtype <- "powexp"
                                  } else if (self$corr[[1]] == "matern") {
                                    if (self$corr[[2]] == 3/2) {covtype <- "matern3_2"}
                                    else if (self$corr[[2]] == 5/2) {covtype <- "matern5_2"}
                                    else {stop("DiceKriging can only do Matern 3/2 and 5/2")}
                                  } else {
                                    print(self$corr)
                                    stop("corr not recognized for DiceKriging, printed above this line")
                                  }
                                }
                                #capture.output(mod1 <- DiceKriging::km(design=X, response=Z, covtype="gauss", nugget.estim=T))
                                capture.output(mod1 <- DiceKriging::km(design=self$X, response=self$Z, covtype=covtype, nugget.estim=self$estimate.nugget, nugget=self$set.nugget, ...))
                                self$mod <- mod1
                              }, #"function to initialize model with data
                              .update = function(...) {#browser()
                                n.since.last.update = nrow(self$X) - self$n.at.last.update
                                if (n.since.last.update < 1) {
                                  message("Can't update, no new X rows")
                                } else {
                                  if (self$n.at.last.update < 10 || n.since.last.update > .25 * self$n.at.last.update) {
                                    # start over if too many
                                    self$.delete(...=...)
                                    self$.init(...=...)
                                  } else {
                                    capture.output(DiceKriging::update(object=mod, newX=self$X[-(1:self$n.at.last.update),], newy=self$Z[-(1:self$n.at.last.update)], nugget.reestim=T))
                                  } #TRYING TO LET UPDATES BE BIG, ELSE UNCOMMENT THIS PART
                                }
                              }, #"function to add data to model or reestimate params
                              .predict = function(XX, se.fit, ...){
                                if (se.fit) {
                                  preds <- DiceKriging::predict.km(self$mod, XX, type = "SK", checkNames=F)
                                  list(fit=preds$mean, se.fit=sqrt(preds$sd))
                                } else {
                                  DiceKriging::predict.km(self$mod, XX, type = "SK", checkNames=F)$mean
                                }
                              }, #"function to predict at new values
                              .predict.se = function(XX, ...) {DiceKriging::predict.km(self$mod, XX, type = "SK", checkNames=F)$sd}, #"function predict the standard error/dev
                              .predict.var = function(XX, ...) {(DiceKriging::predict.km(self$mod, XX, type = "SK", checkNames=F)$sd) ^ 2}, #"function to predict the variance
                              .grad = NULL, # function to calculate the gradient
                              .delete = function(...) {self$mod <- NULL}, #"function to delete model beyond simple deletion
                              .theta = function() {self$mod@covariance@range.val}, #"function to get theta, exp(-theta*(x-x)^2)
                              .nugget = function() {self$mod@covariance@nugget}, #"function to get nugget
                              .mean = NULL # function that gives mean (constant, other functions not implemented)

                            )
)

#' UGP
#' Class providing object with methods for fitting a GP model
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords data, kriging, Gaussian process, regression
#' @return Object of \code{\link{R6Class}} with methods for fitting GP model.
#' @format \code{\link{R6Class}} object.
#' @examples
#' \dontrun{
#' n <- 40
#' d <- 2
#' n2 <- 20
#' f1 <- function(x) {sin(2*pi*x[1]) + sin(2*pi*x[2])}
#' X1 <- matrix(runif(n*d),n,d)
#' Z1 <- apply(X1,1,f1) + rnorm(n, 0, 1e-3)
#' X2 <- matrix(runif(n2*d),n2,d)
#' Z2 <- apply(X2,1,f1)
#' XX1 <- matrix(runif(10),5,2)
#' ZZ1 <- apply(XX1, 1, f1)
#' u <- IGP_sklearn$new(X=X1,Z=Z1)
#' cbind(u$predict(XX1), ZZ1)
#' u$predict.se(XX1)
#' u$update(Xnew=X2,Znew=Z2)
#' u$predict(XX1)
#' u$delete()
#' }
#' @field X Design matrix
#' @field Z Responses
#' @field N Number of data points
#' @field D Dimension of data
#' @section Methods:
#' \describe{
#'   \item{Documentation}{For full documentation of each method go to https://github.com/CollinErickson/UGP/}
#'   \item{\code{new(X=NULL, Z=NULL, package=NULL,
#'   estimate.nugget=T, set.nugget=F, ...)}}{This method
#'   is used to create object of this class with \code{X} and \code{Z} as the data.
#'   The package tells it which package to fit the GP model.}
#'   \item{\code{Xall=NULL, Zall=NULL, Xnew=NULL, Znew=NULL, ...}}{This method
#'   updates the model, adding new data if given, then running optimization again.}}
IGP_sklearn <- R6::R6Class(classname = "IGP_sklearn", inherit = IGP_base,
      public = list(
        pypack = "PythonInR",
        py = list(
          PythonInR = list(
            conn = function() {PythonInR::pyOptions("numpyAlias", "np");if (!PythonInR::pyIsConnected()) {PythonInR::pyConnect()}},
            exec = PythonInR::pyExec,
            assign = function(key, value) {PythonInR::pySet(key=key, value=value, useSetPoly = TRUE, useNumpy = TRUE)},
            get = PythonInR::pyGet,
            close = PythonInR::pyExit
          )#,
          # rPython = list(
          #   conn = function() {},
          #   exec = function() {},#rPython::python.exec,
          #   assign = function() {},#rPython::assign,
          #   get = function() {},#rPython::get
          #   close = function() {}
          # )
        ),
        .init = function(...) {#browser()
          #rPython::python.exec('import sys') # These first two lines need to go
          #rPython::python.exec("sys.path.insert(0, '/Users/collin/anaconda/lib/python2.7/site-packages/')")
          self$py[[self$pypack]]$conn()
          self$py[[self$pypack]]$exec('import numpy as np')
          #self$py[[self$pypack]]$exec('import numpy')
          #rPython::python.exec('import numpy as np')
          #rPython::python.exec('from sklearn import gaussian_process')
          self$py[[self$pypack]]$exec('from sklearn.gaussian_process import GaussianProcessRegressor')
          # rPython::python.exec('from sklearn.gaussian_process import GaussianProcessRegressor')
          self$py[[self$pypack]]$exec("import warnings")
          self$py[[self$pypack]]$exec("warnings.filterwarnings('ignore')")
          # rPython::python.exec("import warnings")
          # rPython::python.exec("warnings.filterwarnings('ignore')")

          self$py[[self$pypack]]$assign("inputdim", ncol(self$X))
          self$py[[self$pypack]]$assign("X1", (self$X))
          self$py[[self$pypack]]$assign("y1", matrix(self$Z, ncol=1))
          self$py[[self$pypack]]$exec('X =  np.matrix(X1)')
          self$py[[self$pypack]]$exec('y = np.matrix(y1).reshape((-1,1))')
          #rPython::python.exec("gp = gaussian_process.GaussianProcess(                      \
          #                     theta0=np.asarray([1e-1 for ijk in range(inputdim)]),       \
          #                     thetaL=np.asarray([1e-4 for ijk in range(inputdim)]),       \
          #                     thetaU=np.asarray([200 for ijk in range(inputdim)]),        \
          #                     optimizer='Welch') ")

          if (!is.null(self$estimate.nugget) || self$set.nugget) {
            warning("GPfit cannot estimate or set the nugget, it picks a stable value")
          }
          if (self$corr[[1]] == "gauss") {
            self$py[[self$pypack]]$exec('from sklearn.gaussian_process.kernels import RBF')
            kernline <- 'kernel = RBF(length_scale=np.asarray([1. for ijk in range(inputdim)]))'
          } else if (self$corr[[1]] == "matern") {
            self$py[[self$pypack]]$exec('from sklearn.gaussian_process.kernels import Matern')
            kernline <- paste0('kernel = Matern(length_scale=np.asarray([1 for ijk in range(inputdim)]), nu=', self$corr[[2]],')')
          } else if (self$corr[[1]] == "matern5_2") {
            self$py[[self$pypack]]$exec('from sklearn.gaussian_process.kernels import Matern')
            kernline <- 'kernel = Matern(length_scale=np.asarray([1 for ijk in range(inputdim)]), nu=2.5)'
          } else if (self$corr[[1]] == "matern3_2") {
            self$py[[self$pypack]]$exec('from sklearn.gaussian_process.kernels import Matern')
            kernline <- 'kernel = Matern(length_scale=np.asarray([1 for ijk in range(inputdim)]), nu=1.5)'
          } else {
            stop("corr not recognized for sklearn")
          }
          #rPython::python.exec('kernel = RBF(length_scale=np.asarray([1. for ijk in range(inputdim)]))') # This and line below added 1/10/17
          self$py[[self$pypack]]$exec(kernline)
          if (!is.null(self$set.nugget) & !self$estimate.nugget) { # set nug and don't estimate
            self$py[[self$pypack]]$exec('gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=10, alpha=',self$set.nugget,')')
          } else if (!is.null(self$set.nugget) & self$estimate.nugget) { # estimate nugget but not given
            self$py[[self$pypack]]$exec('gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=10)')
          } else if (self$estimate.nugget) { # nug not given but estimate it
            self$py[[self$pypack]]$exec('from sklearn.gaussian_process.kernels import WhiteKernel')
            self$py[[self$pypack]]$exec('kernel += WhiteKernel(noise_level=1e-6, noise_level_bounds=(1e-10, 1e5))')
            self$py[[self$pypack]]$exec('gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=10)')
          } else {
            stop("no sklearn option error #928248")
          }
          #rPython::python.exec('gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=10)')
          # Need to give it restarts, just predicted zero when this argument was left out

          self$py[[self$pypack]]$exec("gp.fit(X, y)")

          self$mod <- "sklearn model is in Python"
        }, #"function to initialize model with data
        .update = function(...) {
          self$py[[self$pypack]]$assign("X1", (self$X))
          self$py[[self$pypack]]$assign("y1", matrix(self$Z, ncol=1))
          self$py[[self$pypack]]$exec('X =  np.matrix(X1)')
          self$py[[self$pypack]]$exec('y = np.matrix(y1).reshape((-1,1))')
          self$py[[self$pypack]]$exec("gp.fit(X, y)")
        }, #"function to add data to model or reestimate params
        .predict = function(XX, se.fit, ...) {
          self$py[[self$pypack]]$assign("xp1", XX)
          self$py[[self$pypack]]$exec("xp = np.asmatrix(xp1)")
          self$py[[self$pypack]]$exec("y_pred, sigma2_pred = gp.predict(xp, return_std=True)")
          if (se.fit) {
            list(fit=unlist(self$py[[self$pypack]]$get("y_pred.tolist()")),
                 se.fit=unlist(self$py[[self$pypack]]$get("np.sqrt(sigma2_pred).tolist()")))
          } else {
            unlist(self$py[[self$pypack]]$get("y_pred.tolist()"))
          }
        }, #"function to predict at new values
        .predict.se = function(XX, ...) {
          self$py[[self$pypack]]$assign("xp1", XX)
          self$py[[self$pypack]]$exec("xp = np.asmatrix(xp1)")
          self$py[[self$pypack]]$exec("y_pred, sigma2_pred = gp.predict(xp, return_std=True)")
          unlist(self$py[[self$pypack]]$get("np.sqrt(sigma2_pred).tolist()"))
        }, #"function predict the standard error/dev
        .predict.var = function(XX, ...) {
          self$py[[self$pypack]]$assign("xp1", XX)
          self$py[[self$pypack]]$exec("xp = np.asmatrix(xp1)")
          self$py[[self$pypack]]$exec("y_pred, sigma2_pred = gp.predict(xp, return_std=True)")
          unlist(self$py[[self$pypack]]$get("sigma2_pred.tolist()"))
        }, #"function to predict the variance
        .grad = NULL, # function to calculate the gradient
        .delete = function(...){
          self$py[[self$pypack]]$exec('X =  None')
          self$py[[self$pypack]]$exec('y =  None')
          self$py[[self$pypack]]$exec('xp =  None')
          self$py[[self$pypack]]$exec('X1 =  None')
          self$py[[self$pypack]]$exec('y1 =  None')
          self$py[[self$pypack]]$exec('xp1 =  None')
          self$py[[self$pypack]]$exec('y_pred =  None')
          self$py[[self$pypack]]$exec('sigma2_pred =  None')
          self$py[[self$pypack]]$exec('gp =  None')
          self$py[[self$pypack]]$exec('inputdim =  None')
          self$py[[self$pypack]]$close()
          self$mod <- NULL
        }, #"function to delete model beyond simple deletion
        .theta = function() {rep(NA, ncol(self$X))}, #"function to get theta, exp(-theta*(x-x)^2)
        .nugget = function() {self$py[[self$pypack]]$get('gp.kernel.get_params()')}, #"function to get nugget
        .mean = NULL # function that gives mean (constant, other functions not implemented)

      )
)


#' UGP
#' Class providing object with methods for fitting a GP model
#'
#' @docType class
#' @importFrom R6 R6Class
#' @importFrom PythonInR pyOptions pySet
#' @export
#' @keywords data, kriging, Gaussian process, regression
#' @return Object of \code{\link{R6Class}} with methods for fitting GP model.
#' @format \code{\link{R6Class}} object.
#' @examples
#' \dontrun{
#' n <- 40
#' d <- 2
#' n2 <- 20
#' f1 <- function(x) {sin(2*pi*x[1]) + sin(2*pi*x[2])}
#' X1 <- matrix(runif(n*d),n,d)
#' Z1 <- apply(X1,1,f1) + rnorm(n, 0, 1e-3)
#' X2 <- matrix(runif(n2*d),n2,d)
#' Z2 <- apply(X2,1,f1)
#' XX1 <- matrix(runif(10),5,2)
#' ZZ1 <- apply(XX1, 1, f1)
#' u <- IGP_GPy$new(X=X1,Z=Z1)
#' cbind(u$predict(XX1), ZZ1)
#' u$predict.se(XX1)
#' u$update(Xnew=X2,Znew=Z2)
#' u$predict(XX1)
#' u$delete()
#' }
#' @field X Design matrix
#' @field Z Responses
#' @field N Number of data points
#' @field D Dimension of data
#' @section Methods:
#' \describe{
#'   \item{Documentation}{For full documentation of each method go to https://github.com/CollinErickson/UGP/}
#'   \item{\code{new(X=NULL, Z=NULL, package=NULL,
#'   estimate.nugget=T, set.nugget=F, ...)}}{This method
#'   is used to create object of this class with \code{X} and \code{Z} as the data.
#'   The package tells it which package to fit the GP model.}
#'   \item{\code{Xall=NULL, Zall=NULL, Xnew=NULL, Znew=NULL, ...}}{This method
#'   updates the model, adding new data if given, then running optimization again.}}
IGP_GPy <- R6::R6Class(classname = "IGP_GPy", inherit = IGP_base,
                            public = list(
                              pypack = "PythonInR",
                              py = list(
                                PythonInR = list(
                                  conn = function() {PythonInR::pyOptions("numpyAlias", "np");if (!PythonInR::pyIsConnected()) {PythonInR::pyConnect()}},
                                  exec = PythonInR::pyExec,
                                  assign = function(key, value) {PythonInR::pySet(key=key, value=value, useSetPoly = TRUE, useNumpy = TRUE)},
                                  get = PythonInR::pyGet,
                                  close = PythonInR::pyExit
                                )#,
                                # rPython = list(
                                #   conn = function() {},
                                #   exec = function() {},#rPython::python.exec,
                                #   assign = function() {},#rPython::assign,
                                #   get = function() {},#rPython::get
                                #   close = function() {}
                                # )
                              ),
                              .init = function(...) {browser()
                                if (self$corr[[1]] == "gauss") {
                                  kernline <- 'kernel = GPy.kern.RBF(input_dim=inputdim, ARD=True)'
                                } else if (self$corr[[1]] == "matern") {
                                  stop("GPy only does matern5_2 and matern3_2")
                                  #kernline <- 'kernel = Matern(length_scale=np.asarray([1 for ijk in range(inputdim)]), nu=2.5)'
                                } else if (self$corr[[1]] == "matern5_2") {
                                  kernline <- 'kernel = GPy.kern.Matern52(input_dim=inputdim, ARD=True)'
                                } else if (self$corr[[1]] == "matern3_2") {
                                  kernline <- 'kernel = GPy.kern.Matern52(input_dim=inputdim, ARD=True)'
                                } else {
                                  stop("corr not recognized for GPy")
                                }
                                #rPython::python.exec('import sys') # These first two lines need to go
                                #rPython::python.exec("sys.path.insert(0, '/Users/collin/anaconda/lib/python2.7/site-packages/')")
                                self$py[[self$pypack]]$conn()
                                self$py[[self$pypack]]$exec('import numpy as np')
                                self$py[[self$pypack]]$exec('import GPy')

                                self$py[[self$pypack]]$assign("inputdim", ncol(self$X))
                                self$py[[self$pypack]]$assign("X1", (self$X))
                                self$py[[self$pypack]]$assign("y1", matrix(self$Z, ncol=1))
                                self$py[[self$pypack]]$exec('X =  np.matrix(X1)')
                                self$py[[self$pypack]]$exec('y = np.matrix(y1).reshape((-1,1))')
                                #rPython::python.exec("kernel = GPy.kern.RBF(input_dim=inputdim, variance=1., lengthscale=[1. for iii in range(inputdim)],ARD=True)")
                                self$py[[self$pypack]]$exec(kernline)
                                self$py[[self$pypack]]$exec("gp = GPy.models.GPRegression(X,y,kernel)")
                                if (is.null(self$set.nugget)) {
                                  self$py[[self$pypack]]$exec("gp.likelihood.variance = 1e-8")
                                } else {
                                  self$py[[self$pypack]]$exec(paste0("gp.likelihood.variance = ",self$set.nugget,""))
                                }
                                if (!self$estimate.nugget) {
                                  self$py[[self$pypack]]$exec("gp.likelihood.variance.fix()")
                                }
                                self$py[[self$pypack]]$exec("gp.optimize(messages=False)")
                                self$py[[self$pypack]]$exec("gp.optimize_restarts(num_restarts = 5,  verbose=False)")

                                self$mod <- "GPy model is in Python"
                              }, #"function to initialize model with data
                              .update = function(...) {
                                self$py[[self$pypack]]$assign("X1", (self$X))
                                self$py[[self$pypack]]$assign("y1", matrix(self$Z, ncol=1))
                                self$py[[self$pypack]]$exec('X =  np.matrix(X1)')
                                self$py[[self$pypack]]$exec('y = np.matrix(y1).reshape((-1,1))')
                                self$py[[self$pypack]]$exec("gp.set_XY(X = X, Y = y)")
                                self$py[[self$pypack]]$exec("gp.optimize(messages=False)")
                                self$py[[self$pypack]]$exec("gp.optimize_restarts(num_restarts = 5,  verbose=False)")
                              }, #"function to add data to model or reestimate params
                              .predict = function(XX, se.fit, ...) {
                                self$py[[self$pypack]]$assign("xp1", XX)
                                self$py[[self$pypack]]$exec("xp = np.asmatrix(xp1)")
                                self$py[[self$pypack]]$exec("y_pred, sigma2_pred = gp.predict(np.asarray(xp))")
                                if (se.fit) {
                                  list(fit=unlist(self$py[[self$pypack]]$get("y_pred.tolist()")),
                                       se.fit=unlist(self$py[[self$pypack]]$get("np.sqrt(sigma2_pred).tolist()")))
                                } else {
                                  unlist(self$py[[self$pypack]]$get("y_pred.tolist()"))
                                }
                              }, #"function to predict at new values
                              .predict.se = function(XX, ...) {
                                self$py[[self$pypack]]$assign("xp1", XX)
                                self$py[[self$pypack]]$exec("xp = np.asmatrix(xp1)")
                                self$py[[self$pypack]]$exec("y_pred, sigma2_pred = gp.predict(np.asarray(xp))")
                                unlist(self$py[[self$pypack]]$get("np.sqrt(sigma2_pred).tolist()"))
                              }, #"function predict the standard error/dev
                              .predict.var = function(XX, ...) {
                                self$py[[self$pypack]]$assign("xp1", XX)
                                self$py[[self$pypack]]$exec("xp = np.asmatrix(xp1)")
                                self$py[[self$pypack]]$exec("y_pred, sigma2_pred = gp.predict(np.asarray(xp))")
                                unlist(self$py[[self$pypack]]$get("sigma2_pred.tolist()"))
                              }, #"function to predict the variance
                              .grad = NULL, # function to calculate the gradient
                              .delete = function(...){
                                self$py[[self$pypack]]$exec('X =  None')
                                self$py[[self$pypack]]$exec('y =  None')
                                self$py[[self$pypack]]$exec('xp =  None')
                                self$py[[self$pypack]]$exec('X1 =  None')
                                self$py[[self$pypack]]$exec('y1 =  None')
                                self$py[[self$pypack]]$exec('xp1 =  None')
                                self$py[[self$pypack]]$exec('y_pred =  None')
                                self$py[[self$pypack]]$exec('sigma2_pred =  None')
                                self$py[[self$pypack]]$exec('gp =  None')
                                self$py[[self$pypack]]$exec('kernel =  None')
                                self$py[[self$pypack]]$exec('inputdim =  None')
                                self$py[[self$pypack]]$close()
                                self$mod <- NULL
                              }, #"function to delete model beyond simple deletion
                              .theta = function() {rep(NA, ncol(self$X))}, #"function to get theta, exp(-theta*(x-x)^2)
                              .nugget = function() {self$py[[self$pypack]]$get('gp.likelihood.variance')}, #"function to get nugget
                              .mean = NULL # function that gives mean (constant, other functions not implemented)

                            )
)





#' UGP
#' Class providing object with methods for fitting a GP model
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords data, kriging, Gaussian process, regression
#' @return Object of \code{\link{R6Class}} with methods for fitting GP model.
#' @format \code{\link{R6Class}} object.
#' @examples
#' \dontrun{
#' n <- 40
#' d <- 2
#' n2 <- 20
#' f1 <- function(x) {sin(2*pi*x[1]) + sin(2*pi*x[2])}
#' X1 <- matrix(runif(n*d),n,d)
#' Z1 <- apply(X1,1,f1) + rnorm(n, 0, 1e-3)
#' X2 <- matrix(runif(n2*d),n2,d)
#' Z2 <- apply(X2,1,f1)
#' XX1 <- matrix(runif(10),5,2)
#' ZZ1 <- apply(XX1, 1, f1)
#' u <- IGP_DACE$new(X=X1,Z=Z1)
#' cbind(u$predict(XX1), ZZ1)
#' u$predict.se(XX1)
#' u$update(Xnew=X2,Znew=Z2)
#' u$predict(XX1)
#' u$delete()
#' }
#' @field X Design matrix
#' @field Z Responses
#' @field N Number of data points
#' @field D Dimension of data
#' @section Methods:
#' \describe{
#'   \item{Documentation}{For full documentation of each method go to https://github.com/CollinErickson/UGP/}
#'   \item{\code{new(X=NULL, Z=NULL, package=NULL,
#'   estimate.nugget=T, set.nugget=F, ...)}}{This method
#'   is used to create object of this class with \code{X} and \code{Z} as the data.
#'   The package tells it which package to fit the GP model.}
#'   \item{\code{Xall=NULL, Zall=NULL, Xnew=NULL, Znew=NULL, ...}}{This method
#'   updates the model, adding new data if given, then running optimization again.}}
IGP_DACE <- R6::R6Class(classname = "IGP_DACE", inherit = IGP_base,
                            public = list(
                              matlab_path = "C:\\Users\\cbe117\\School\\DOE\\GP_codes\\DACE\\dace",
                              .init = function(...) {#browser()

                                R.matlab::Matlab$startServer()
                                matlab <- R.matlab::Matlab()
                                self$mod <- matlab
                                isOpen <- open(matlab)
                                if (!isOpen) throw("MATLAB server is not running: waited 30 seconds.")

                                # set a variable in R and send to MATLAB
                                R.matlab::setVariable(matlab, X = self$X)
                                R.matlab::setVariable(matlab, Z = self$Z)
                                #R.matlab::setVariable(matlab, meanfunc = '')
                                R.matlab::setVariable(matlab, theta = 1)
                                R.matlab::setVariable(matlab, lob = 1e-4)
                                R.matlab::setVariable(matlab, upb = 1e4)
                                R.matlab::evaluate(matlab, 'meanfunc = @regpoly0')
                                R.matlab::evaluate(matlab, 'corrfunc = @corrgauss')
                                R.matlab::evaluate(matlab, paste0("addpath('",self$matlab_path,"');"))
                                R.matlab::evaluate(matlab, "[dmodel, perf] = dacefit(X, Z, meanfunc, corrfunc, theta, lob, upb);")
                                #R.matlab::evaluate(matlab, "y=20; z=x+y")
                                #z <- R.matlab::getVariable(matlab, "z")
                                #z
                                #close(matlab)

                                #R.matlab::setVariable(matlab, X1=X1)
                                #temp <- R.matlab::evaluate(matlab, "X1 .* X1", capture=TRUE)
                                #temp
                              }, #"function to initialize model with data
                              .update = function() { # function to add data to model or reestimate params
                                matlab <- self$mod
                                R.matlab::setVariable(matlab, X = self$X)
                                R.matlab::setVariable(matlab, Z = self$Z)
                                R.matlab::evaluate(matlab, "[dmodel, perf] = dacefit(X, Z, meanfunc, corrfunc, theta, lob, upb);")
                              },
                              .predict = function(XX, se.fit, ...) {#browser()
                                R.matlab::setVariable(self$mod, XX = XX)
                                R.matlab::evaluate(self$mod, '[YP, MSEP] = predictor(XX, dmodel);')
                                YP <- R.matlab::getVariable(self$mod, 'YP')
                                MSEP <- R.matlab::getVariable(self$mod, 'MSEP')
                                if (se.fit) {
                                  cbind(YP$YP, sqrt(MSEP$MSEP))
                                } else {
                                  YP$YP
                                }
                              }, #"function to predict at new values
                              .predict.se = function(XX, ...) {#browser()
                                R.matlab::setVariable(self$mod, XX = XX)
                                R.matlab::evaluate(self$mod, '[YP, MSEP] = predictor(XX, dmodel);')
                                YP <- R.matlab::getVariable(self$mod, 'YP')
                                MSEP <- R.matlab::getVariable(self$mod, 'MSEP')
                                #cbind(YP$YP, sqrt(MSEP$MSEP))
                                sqrt(MSEP$MSEP)
                              }, #"function predict the standard error/dev
                              .predict.var = function(XX, ...) {#browser()
                                R.matlab::setVariable(self$mod, XX = XX)
                                R.matlab::evaluate(self$mod, '[YP, MSEP] = predictor(XX, dmodel);')
                                # YP <- R.matlab::getVariable(self$mod, 'YP')
                                MSEP <- R.matlab::getVariable(self$mod, 'MSEP')
                                # cbind(YP$YP, MSEP$MSEP)
                                MSEP$MSEP
                              }, #"function to predict the variance
                              .grad = NULL, # function to calculate the gradient
                              .delete = function() {
                                close(self$mod)
                                #close(matlab)
                              }, #"function to delete model beyond simple deletion
                              .theta = function() {
                                R.matlab::getVariable(matlab, 'dmodel.theta')
                              }, #"function to get theta, exp(-theta*(x-x)^2)
                              .nugget = NULL, #"function to get nugget
                              .mean = function() {
                                R.matlab::getVariable(matlab, 'dmodel.beta')
                              } # function that gives mean (constant, other functions not implemented)

                            )
)






#' UGP
#' Class providing object with methods for fitting a GP model
#'
#' @docType class
#' @importFrom R6 R6Class
# @export
#' @keywords data, kriging, Gaussian process, regression
#' @return Object of \code{\link{R6Class}} with methods for fitting GP model.
#' @format \code{\link{R6Class}} object.
#' @examples
#' \dontrun{
#' n <- 40
#' d <- 2
#' n2 <- 20
#' f1 <- function(x) {sin(2*pi*x[1]) + sin(2*pi*x[2])}
#' X1 <- matrix(runif(n*d),n,d)
#' Z1 <- apply(X1,1,f1) + rnorm(n, 0, 1e-3)
#' X2 <- matrix(runif(n2*d),n2,d)
#' Z2 <- apply(X2,1,f1)
#' XX1 <- matrix(runif(10),5,2)
#' ZZ1 <- apply(XX1, 1, f1)
#' u <- IGP_ooDACE$new(X=X1,Z=Z1)
#' cbind(u$predict(XX1), ZZ1)
#' u$predict.se(XX1)
#' u$update(Xnew=X2,Znew=Z2)
#' u$predict(XX1)
#' u$delete()
#' }
#' @field X Design matrix
#' @field Z Responses
#' @field N Number of data points
#' @field D Dimension of data
#' @section Methods:
#' \describe{
#'   \item{Documentation}{For full documentation of each method go to https://github.com/CollinErickson/UGP/}
#'   \item{\code{new(X=NULL, Z=NULL, package=NULL,
#'   estimate.nugget=T, set.nugget=F, ...)}}{This method
#'   is used to create object of this class with \code{X} and \code{Z} as the data.
#'   The package tells it which package to fit the GP model.}
#'   \item{\code{Xall=NULL, Zall=NULL, Xnew=NULL, Znew=NULL, ...}}{This method
#'   updates the model, adding new data if given, then running optimization again.}}
IGP_ooDACE <- R6::R6Class(classname = "IGP_ooDACE", inherit = IGP_base,
                         public = list(
                           .init = function(...) {browser()

                             R.matlab::Matlab$startServer()
                             matlab <- R.matlab::Matlab()
                             self$mod <- matlab
                             isOpen <- open(matlab)
                             if (!isOpen) throw("MATLAB server is not running: waited 30 seconds.")

                             # set a variable in R and send to MATLAB
                             R.matlab::setVariable(matlab, X = self$X)
                             R.matlab::setVariable(matlab, inputdim = ncol(self$X))
                             R.matlab::setVariable(matlab, Z = self$Z)
                             #R.matlab::setVariable(matlab, meanfunc = '')
                             R.matlab::setVariable(matlab, theta = 1)
                             R.matlab::setVariable(matlab, lob = 1e-4)
                             R.matlab::setVariable(matlab, upb = 1e4)
                             R.matlab::evaluate(matlab, 'meanfunc = @regpoly0')
                             R.matlab::evaluate(matlab, 'corrfunc = @corrgauss')
                             R.matlab::evaluate(matlab, "addpath(genpath('C:\\Users\\cbe117\\School\\DOE\\GP_codes\\ooDACE'));")
                             # R.matlab::evaluate(matlab, "addpath('C:\\Users\\cbe117\\School\\DOE\\GP_codes\\ooDACE');")
                             # R.matlab::evaluate(matlab, "addpath('C:\\Users\\cbe117\\School\\DOE\\GP_codes\\ooDACE\\optimizers\');")
                             #R.matlab::evaluate(matlab, "[dmodel, perf] = DACEfit(X, Z, meanfunc, corrfunc, theta, lob, upb);")
                             browser()
                             R.matlab::evaluate(matlab, "opts = Kriging.getDefaultOptions();")
                             R.matlab::evaluate(matlab, "opts.hpBounds = [lob ; upb]; % hyperparameter optimization bounds")
                             R.matlab::evaluate(matlab, "opts.hpOptimizer = SQPLabOptimizer( inputdim, 1 );")
                             if (self$estimate.nugget) {
                               R.matlab::evaluate(matlab, "opts.lambda0 = 0;")
                               R.matlab::evaluate(matlab, "opts.lambdaBounds = [-5; 5];")
                             }
                             R.matlab::evaluate(matlab, "")
                             R.matlab::evaluate(matlab, "")
                             R.matlab::evaluate(matlab, "")
                             R.matlab::evaluate(matlab, "k = Kriging( opts, theta, meanfunc, corrfunc);")
                             R.matlab::evaluate(matlab, "k = k.fit( X, Z);")
                             #R.matlab::evaluate(matlab, "y=20; z=x+y")
                             #z <- R.matlab::getVariable(matlab, "z")
                             #z
                             #close(matlab)

                             #R.matlab::setVariable(matlab, X1=X1)
                             #temp <- R.matlab::evaluate(matlab, "X1 .* X1", capture=TRUE)
                             #temp
                           }, #"function to initialize model with data
                           .update = function (...) { #function to add data to model or reestimate params
                             matlab <- self$mod
                             R.matlab::setVariable(matlab, X = self$X)
                             R.matlab::setVariable(matlab, Z = self$Z)
                             R.matlab::evaluate(matlab, "k = k.fit( X, Z);")
                           },
                           .predict = function(XX, se.fit, ...) {browser()
                             R.matlab::setVariable(self$mod, XX = XX)
                             R.matlab::evaluate(self$mod, '[YP, MSEP] = k.predict(XX);')
                             YP <- R.matlab::getVariable(self$mod, 'YP')
                             MSEP <- R.matlab::getVariable(self$mod, 'MSEP')
                             if (se.fit) {
                               cbind(YP$YP, sqrt(MSEP$MSEP))
                             } else {
                               YP$YP
                             }
                           }, #"function to predict at new values
                           .predict.se = function(XX, ...) {browser()
                             R.matlab::setVariable(self$mod, XX = XX)
                             R.matlab::evaluate(self$mod, '[YP, MSEP] = k.predict(XX);')
                             YP <- R.matlab::getVariable(self$mod, 'YP')
                             MSEP <- R.matlab::getVariable(self$mod, 'MSEP')
                             YP$YP
                             cbind(YP$YP, sqrt(MSEP$MSEP))
                           }, #"function predict the standard error/dev
                           .predict.var = function(XX, ...) {browser()
                             R.matlab::setVariable(self$mod, XX = XX)
                             R.matlab::evaluate(self$mod, '[YP, MSEP] = k.predict(XX);')
                             YP <- R.matlab::getVariable(self$mod, 'YP')
                             MSEP <- R.matlab::getVariable(self$mod, 'MSEP')
                             cbind(YP$YP, MSEP$MSEP)
                           }, #"function to predict the variance
                           .grad = NULL, # function to calculate the gradient
                           .delete = function() {
                             close(self$mod)
                             #close(matlab)
                           }, #"function to delete model beyond simple deletion
                           .theta = function() {
                             R.matlab::getVariable(matlab, 'k.getHyperparameters')
                           }, #"function to get theta, exp(-theta*(x-x)^2)
                           .nugget = NULL, #"function to get nugget
                           .mean = function() {
                             R.matlab::getVariable(matlab, 'methods(k)') # can't do this yet
                           } # function that gives mean (constant, other functions not implemented)

                         )
)




#' UGP
#' Class providing object with methods for fitting a GP model
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords data, kriging, Gaussian process, regression
#' @return Object of \code{\link{R6Class}} with methods for fitting GP model.
#' @format \code{\link{R6Class}} object.
#' @examples
#' \dontrun{
#' n <- 40
#' d <- 2
#' n2 <- 20
#' f1 <- function(x) {sin(2*pi*x[1]) + sin(2*pi*x[2])}
#' X1 <- matrix(runif(n*d),n,d)
#' Z1 <- apply(X1,1,f1) + rnorm(n, 0, 1e-3)
#' X2 <- matrix(runif(n2*d),n2,d)
#' Z2 <- apply(X2,1,f1)
#' XX1 <- matrix(runif(10),5,2)
#' ZZ1 <- apply(XX1, 1, f1)
#' u <- IGP_sklearn$new(X=X1,Z=Z1)
#' cbind(u$predict(XX1), ZZ1)
#' u$predict.se(XX1)
#' u$update(Xnew=X2,Znew=Z2)
#' u$predict(XX1)
#' u$delete()
#' }
#' @field X Design matrix
#' @field Z Responses
#' @field N Number of data points
#' @field D Dimension of data
#' @section Methods:
#' \describe{
#'   \item{Documentation}{For full documentation of each method go to https://github.com/CollinErickson/UGP/}
#'   \item{\code{new(X=NULL, Z=NULL, package=NULL,
#'   estimate.nugget=T, set.nugget=F, ...)}}{This method
#'   is used to create object of this class with \code{X} and \code{Z} as the data.
#'   The package tells it which package to fit the GP model.}
#'   \item{\code{Xall=NULL, Zall=NULL, Xnew=NULL, Znew=NULL, ...}}{This method
#'   updates the model, adding new data if given, then running optimization again.}}
IGP_laGP_GauPro <- R6::R6Class(classname = "IGP_lagP_GauPro", inherit = IGP_base,
                           public = list(
                             .init = function(...) {#browser()
                               # Fit model to data with laGP
                               self$mod.extra$laGP <- IGP(X=self$X, Z=self$Z, package="laGP",
                                                          corr=self$corr,
                                                          estimate.nugget=self$estimate.nugget,
                                                          set.nugget=self$set.nugget)
                               #self$mod.extra$laGP$init(X=self$X, Z=self$Z, ...)

                               # Copy params to GauPro, don't fit, use this for predicting
                               self$mod.extra$GauPro <- IGP(X=self$X, Z=self$Z, package="GauPro",
                                                          corr=self$corr,
                                                          estimate.nugget=FALSE,
                                                          set.nugget=self$mod.extra$laGP$nugget(),
                                                          theta=self$mod.extra$laGP$theta(),
                                                          #nug=self$mod.extra$laGP$nug,
                                                          param.est=FALSE)
                               #self$mod.extra$GauPro$init(X=self$X, Z=self$Z,
                               #                          theta=self$mod.extra$laGP$theta,
                               #                           nug=self$mod.extra$laGP$nug)
                               self$mod <- "laGP model is mod.extra$laGP, GauPro model is mod.extra$GauPro. This fits with laGP but predicts with GauPro"
                             }, #"function to initialize model with data
                             .update = function(..., no_update=FALSE) {
                               # Update model in laGP
                               if (!no_update) { # won't have it update data even not updating params since I don't like when it gives issues, and we pass Xall anyways
                                 self$mod.extra$laGP$update(Xall=self$X, Zall=self$Z,
                                                            no_update=no_update,
                                                            ...)
                               }
                               # Pass GauPro new theta and nugget if it was updated
                               if (!no_update) {
                                 self$mod.extra$GauPro$mod$theta <- self$mod.extra$laGP$theta()
                                 self$mod.extra$GauPro$mod$nug <- self$mod.extra$laGP$nugget()
                               }
                               self$mod.extra$GauPro$update(Xall=self$X, Zall=self$Z,
                                                            no_update=TRUE)

                             }, #"function to add data to model or reestimate params
                             .predict = function(XX, se.fit, ...) {
                               self$mod.extra$GauPro$.predict(XX=XX, se.fit=se.fit, ...)
                             }, #"function to predict at new values
                             .predict.se = function(XX, ...) {
                               self$mod.extra$GauPro$.predict.se(XX=XX, ...)
                             }, #"function predict the standard error/dev
                             .predict.var = function(XX, ...) {
                               self$mod.extra$GauPro$.predict.var(XX=XX, ...)
                             }, #"function to predict the variance
                             .grad = function(XX) {self$mod.extra$GauPro$grad(XX=XX)}, # function to calculate the gradient
                             .delete = function(...){
                               if (!is.null(self$mod.extra)) {
                                 self$mod.extra$laGP$delete()
                                 self$mod.extra$GauPro$delete()
                                 self$mod.extra <- NULL
                               }
                              self$mod <- NULL
                             }, #"function to delete model beyond simple deletion
                             .theta = function() {self$mod.extra$GauPro$theta()}, #"function to get theta, exp(-theta*(x-x)^2)
                             .nugget = function() {self$mod.extra$GauPro$nugget()}, #"function to get nugget
                             .mean = NULL # function that gives mean (constant, other functions not implemented)

                           )
)
