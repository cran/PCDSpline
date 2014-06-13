\name{Ispline}
\alias{Ispline}

\title{
Ispline
}
\description{
Generates the I-spline basis matrix associated with integrated spline basis functions. Created by Cai and Wang in October, 2009. Details can be found in Ramsay (1988).
}
\usage{
Ispline(x, order, knots)
}
\arguments{
  \item{x}{a vector of observations.}
  \item{order}{the order of the basis functions.}
  \item{knots}{a sequence of increasing points specifying the placement of the knots.}
}
\value{
An I-spline basis matrix of dimension c(length(knots)+order-2, length(x)).
}

\references{
Ramsay, J. (1988). Monotone regression splines in action. \emph{Statistical Science} 3, 425-441.       
}
