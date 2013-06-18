\name{loglike}
\alias{loglike}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Internal glmlep functions
}
\description{
 Internal glmlep functions
}
\usage{
loglike(x, y, beta, family = c("gaussian", "binomial"))
}
\arguments{
  \item{x}{
	The design matrix, without an intercept.
}
  \item{y}{
	The response vector. Quantitative for family="gaussian". For family="binomial" should be a vector with two levels.
}
  \item{beta}{
 The estimated coefficients.
}
  \item{family}{
 Response type (see above)
}
}
\details{
These are not intended for use by users.}

\author{
 Canhong Wen, Hao Lin, Shaoli Wang and Xueqin Wang.

Maintainer: Canhong Wen <wencanhong@gmail.com>
}

