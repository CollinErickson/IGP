addpath(genpath('C:/Users/cbe117/Documents/R/win-library/3.4/UGP/gpml-matlab-v4.0-2016-10-19'))
tout <- R.matlab::evaluate.Matlab(
  self$mod,
  "
  x = gpml_randn(0.8, 20, 1);                 % 20 training inputs
  y = sin(3*x) + 0.1*gpml_randn(0.9, 20, 1);  % 20 noisy training targets
  xs = linspace(-3, 3, 61)';                  % 61 test inputs
  meanfunc = [];                    % empty: don't use a mean function
  covfunc = @covSEiso;              % Squared Exponental covariance function
  likfunc = @likGauss;              % Gaussian likelihood
  hyp = struct('mean', [], 'cov', [0 0], 'lik', -1);
  hyp2 = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, x, y);
  ",
  capture=TRUE
)
print(tout)


# Start all here
set.seed(0)
n <- 40
d <- 2
n2 <- 20
f1 <- function(x) {sin(2*pi*x[1]) + sin(2*pi*x[2])}
X1 <- matrix(runif(n*d),n,d)
Z1 <- apply(X1,1,f1) + rnorm(n, 0, 1e-3)
# X2 <- matrix(runif(n2*d),n2,d)
# Z2 <- apply(X2,1,f1)
# XX1 <- matrix(runif(10),5,2)
# ZZ1 <- apply(XX1, 1, f1)

R.matlab::Matlab$startServer()
matlab <- R.matlab::Matlab()
isOpen <- open(matlab)
print(R.matlab::evaluate.Matlab(matlab, '1+2', capture=TRUE))
GPML_file_path <- system.file("gpml-matlab-v4.0-2016-10-19", package="UGP")
# addpath(genpath('C:/Users/cbe117/Documents/R/win-library/3.4/UGP/gpml-matlab-v4.0-2016-10-19'))
R.matlab::evaluate(matlab, paste0("addpath(genpath('", GPML_file_path, "'));"))

R.matlab::setVariable(matlab, X = X1)
R.matlab::setVariable(matlab, Z = Z1)

R.matlab::evaluate(matlab, 'meanfunc = @meanConst; hyp.mean = [0;];')
R.matlab::evaluate(matlab, 'corrfunc = @corrSEard; hyp.cov = [zeros(size(X, 2), 1); 0;]; hyp.lik = log(0.1);')
R.matlab::evaluate(matlab, 'likfunc = @likGauss')
R.matlab::evaluate(matlab, 'hyp = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, corrfunc, likfunc, X, Z);')

close(matlab)



filepath <- tempfile(fileext='.m')
filetext <- paste0("
1+2
addpath(genpath('", GPML_file_path, "'));
X = [0; .3; .5; .7; .9; 1;];
Z = [0; .4; .5; .6; .7; 1;];
meanfunc = @meanConst; hyp.mean = [0;];
corrfunc = @covSEard; hyp.cov = [zeros(size(X, 2), 1); 0;]; hyp.lik = log(0.1);
likfunc = @likGauss
hyp = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, corrfunc, likfunc, X, Z);
")
write(x=filetext, file=filepath)
# system(paste0("matlab ", filepath))
system(paste0("matlab -nodisplay -nosplash -minimize -r \"run('",filepath,"'); exit\""))
