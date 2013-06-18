\name{glmlep-package}
\alias{glmlep-package}
\docType{package}
\title{
 Regularization paths for LEP-penalized regression models
}
\description{
Efficient algorithms for fitting regularization paths for linear or logistic regression models penalized by LEP.
}
\details{
\tabular{ll}{
Package: \tab glmlep\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2013-06-05\cr
License: \tab GPL-2\cr
}
Accepts a design matrix X and vector of responses y, produces the regularization path over a grid of values for the tuning parameter lambda. Also provides methods for plotting and for determining locally convex regions of the coefficients paths.
}
\author{
 Canhong Wen, Hao Lin, Shaoli Wang and Xueqin Wang.

Maintainer: Canhong Wen <wencanhong@gmail.com>

}
\references{
Wen, C., Wang, X., & Wang, S. (2013). Laplace Error Penalty based variable selection in ultra high-dimension. In press.
}

\examples{

## generate data
require(mvtnorm)
n <- 100;
beta <- c(3,1.5,0,0,2,0,0,0)

set.seed(100)
p <- length(beta);
corr_data <- diag(rep(1,p));

x <- as.matrix(rmvnorm(n,rep(0,p),corr_data))
noise <- rnorm(n)

y <- tcrossprod(x,t(beta)) + noise;
fit <- glmlep(x,y,family="gaussian")
}
