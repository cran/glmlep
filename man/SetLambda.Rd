\name{SetLambda}
\alias{SetLambda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Internal glmlep functions
}
\description{
 Internal glmlep functions
}
\usage{
SetLambda(x, y, lambda.min, nlambda, penalty.factor)
}
\arguments{
  \item{x}{
The design matrix, without an intercept.
}
  \item{y}{
The response vector. Quantitative for family="gaussian". For family="binomial" should be a vector with two levels.
}
  \item{lambda.min}{
Smallest value for lambda, as a fraction of lambda.max, the (data derived) entry value (i.e. the smallest value for which all coefficients are zero). The default depends on the sample size nobs relative to the number of variables nvars. If nobs > nvars, the default is 0.001, close to zero. If nobs < nvars, the default is 0.05. 
}
  \item{nlambda}{
 The number of \code{lambda} values; default is 100.
}
  \item{penalty.factor}{
 Separate penalty factors can be applied to each coefficient. This is a number that multiplies lambda to allow differential shrinkage. 
}
}
\details{
These are not intended for use by users.
}

\author{
 Canhong Wen, Hao Lin, Shaoli Wang and Xueqin Wang.

Maintainer: Canhong Wen <wencanhong@gmail.com>
}