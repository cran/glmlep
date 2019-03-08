\name{cv.glmlep}
\alias{cv.glmlep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cross-validation for \code{glmlep}.
}
\description{
 Does k-fold cross-validation for \code{glmlep}, produces a plot, and returns a value for \code{lambda}.
}
\usage{
cv.glmlep(x, y, family = c("gaussian", "binomial"), lambda = NULL, 
lambda.min = ifelse(n < p, 0.05, 0.001), nlambda = 100, lambda2 = 0, 
kappa = ifelse(n < p, 0.1, 0.05), pen.fac = rep(1, p), tol = 1e-06, 
max.ite = 1000, foldid, nfolds = 5, cv.seed = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{
	The design matrix, without an intercept.
}
\item{y}{
	The response vector. Quantitative for family="gaussian". For family="binomial" should be a vector with two levels.
}
\item{family}{
 Response type (see above)
}
  \item{lambda}{
A user supplied lambda sequence. Typical usage is to have the program compute its own lambda sequence based on nlambda and lambda.min.ratio. Supplying a value of lambda overrides this. WARNING: use with care. Do not supply a single value for lambda. Supply instead a decreasing sequence of lambda values. glmnet relies on its warms starts for speed, and its often faster to fit a whole path than compute a single fit.
}
  \item{lambda.min}{
Smallest value for lambda, as a fraction of lambda.max, the (data derived) entry value (i.e. the smallest value for which all coefficients are zero). The default depends on the sample size nobs relative to the number of variables nvars. If nobs > nvars, the default is 0.001, close to zero. If nobs < nvars, the default is 0.05. 
}
  \item{nlambda}{
 The number of \code{lambda} values; default is 100.
}
  \item{lambda2}{
 The tuning parameter for additional L_2 penalty. Use for better grouping effect. The default is 0.
}
  \item{kappa}{
  The scale tuning parameter of the LEP penalty. One can specify it to get the desired estimates because of the homotopy of LEP function to the L_0 function. If nobs > nvars, the default is 0.05, close to zero. If nobs < nvars, the default is 0.1.
}
  \item{pen.fac}{
 Separate penalty factors can be applied to each coefficient. This is a number that multiplies lambda to allow differential shrinkage. Can be 0 for some variables, which implies no shrinkage, and that variable is always included in the model. Default is 1 for all variables (and implicitly infinity for variables listed in exclude). Note: the penalty factors are internally rescaled to sum to nobs, and the lambda sequence will reflect this change.
}
  \item{tol}{
  Convergence tolerance for MCD. Each inner MCD loop continues until the change in the estimates is less than \code{tol}. default is 1E-6.
}
  \item{max.ite}{
  Maximum number of passes over the data for all lambda values; default is 10^3.
}
  \item{foldid}{
	An optional vector of values between 1 and \code{nfold} identifying what fold each observation is in. If supplied, nfold can be missing.
}
  \item{nfolds}{
 Number of folds - default is 5.
}
  \item{cv.seed}{
 The seed for cross-validation. This could be used for simulation replicability.
}
}
\details{
 The function runs glmlep nfolds+1 times; the first to get the lambda sequence and the final estimate, and then the remainder to compute the fit with each of the folds omitted. The loss is accumulated, and the average loss over the folds is computed. Note that cv.glmlep does NOT search for values for \code{kappa}. A specific value should be supplied, else \code{kappa}=0.05 is assumed by default. If users would like to cross-validate \code{kappa} as well, they should call cv.glmlep with a pre-computed vector foldid, and then use this same fold vector in separate calls to cv.glmlep with different values of \code{kappa}. Note that \code{n} is the sample size and \code{p} is the dimension of variables.
}
\value{
	An object of class "cv.glmlep" is returned, which is a list with the ingredients of the cross-validation fit.
    \item{beta }{A nrow(x) x length(\code{lambda}) matrix of estimated coefficient.}
    \item{lambda }{The sequence of regularization parameter values used}
	\item{df}{The degree of freedom for each value of \code{lambda}.}
	\item{loss}{The -2*log-likelihood value for each value of \code{lambda}.}
    \item{lambda.min }{The value of lambda with the minimum EBIC.}
	\item{beta.min }{The coefficient with the minimum EBIC.}
	\item{call}{The call that produces this object}
}
\references{
Wen, C., Wang, X., & Wang, S. (2013). Laplace Error Penalty based variable selection in ultra high-dimension. In press.
}
\author{
 Canhong Wen, Hao Lin, Shaoli Wang and Xueqin Wang.

Maintainer: Canhong Wen <wencanhong@gmail.com>
}
\examples{

## generate data from multivariate normal distribution
require(mvtnorm)
n = 100;
beta <- c(3,1.5,0,0,2,0,0,0)

set.seed(100)
p=length(beta);
corr_data=diag(rep(1,p));
x=as.matrix(rmvnorm(n,rep(0,p),corr_data))
noise=rnorm(n);

## Gaussian
y <- tcrossprod(x,t(beta)) + noise;
fit <- cv.glmlep(x,y,family="gaussian")


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ LEP }
\keyword{ GLM }
\keyword{ Cross-validation }
