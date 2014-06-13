\name{PCDHess.wf}
\alias{PCDHess.wf}

\title{
Calculating the Hessian matrix using Louis's method (1982)
}
\description{
Calculating the Hessian matrix using Louis's Method under the Gamma frailty non-homogeneous Poisson process model. This is a support function for \code{\link{PCDReg.wf}}.
}
\usage{
PCDHess.wf(DATA, beta, gamma, nu, order, knots)
}

\arguments{
  \item{DATA}{
use specified data structure.
}
  \item{beta}{
estimates of regression coefficients.
}
  \item{gamma}{
estimates of spline coefficients.
}
  \item{nu}{
estimates of gamma frailty variance parameter.
}
  \item{order}{
the order of basis functions.
}
  \item{knots}{
the equally spaced knots.
}
}
\details{
To obtain the Hessian matrix of the observed likelihood evaluated at the last step output of the EM algorithm.
}
\value{
\item{HESS}{Hessian matrix.}
}

\references{
Louis, T. (1982). Finding the observed information matrix when using the EM algorithm. \emph{Journal
of the Royal Statistical Society}, Series B 44, 226-233.

Yao,B., Wang, L., and He, X. (2014+).Semiparametric regression analysis of panel count data allowing for within-subject correlation.
}


\keyword{Hessian matrix}
