\name{PCDReg.nf}
\alias{PCDReg.nf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Regression analysis of panel count data under the non-homogeneous Poisson process model
}
\description{
Fits the nonhomogeneous Poisson process model to panel count data using EM algorithm.
}
\usage{
PCDReg.nf(DATA, order, placement, nknot, myknots, binit, ginit, t.seq, tol)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DATA}{
  use specified data structure.
}
  \item{order}{
  the order of basis functions.
}
  \item{placement}{
  logical, if TRUE knots are placed evenly across the observed intervals based on the input data set; if FALSE knots should be specified     by the user. see myknots.
  }
  \item{nknot}{
  the number of knots to be used.
}
  \item{myknots}{
  a sequence of increasing points whose length is nknot.
}
  \item{binit}{
  initial estimate of regression coefficients.
}
 \item{ginit}{
 initial estimate of spline coefficients whose length should be (order+nknot-2) or (order+length(myknots)-2).
}
 \item{t.seq}{
 an increasing sequence of points at which the baseline mean function is evaluated.
 }
 \item{tol}{
 the convergence criterion of the EM algorithm.
 }
 }
 
\details{
The above function fits the non-homogeneous Poisson process model to panel count data via EM algorithm. 
}
\value{
\item{beta}{estimates of regression coefficients.}
\item{gamma}{estimates of spline coefficients.}
\item{var.b}{the variance covariance matrix of regression coefficients.}
\item{Hess}{Hessian matrix.}
\item{knots}{the knots used in the analysis.}
\item{bmf}{estimated baseline mean function evaluated at the points t.seq; use pmf to plot the baseline mean fuction.}
\item{AIC}{the Akaike information criterion.}
\item{BIC}{the Bayesian information/Schwarz criterion.}
\item{flag}{the indicator whether the Hessian matrix is non-singular. When flag="TRUE",the variance estimate may not be accurate. }
}
\references{
Yao, B., Wang, L., and He, X. (2014+). Semiparametric regression analysis of panel count data allowing for within-subject correlation.
}

\note{
The non-homogeneous Poisson process model involves no Gamma frailty.
}

\seealso{
\code{\link{PCDReg.wf}}
}
\examples{
##Simulated Data

n=13; #number of subjects

##generate the number of observations for each subject
k=rpois(n,6)+1; K=max(k);
  
##generate random time gaps for each subject
y=matrix(,n,K);
for (i in 1:n){y[i,1:k[i]]=rexp(k[i],1)} 

##get observation time points for each subject
t=matrix(,n,K);
for (i in 1:n){
  for (j in 2:K){
    t[i,1] = y[i,1]
    t[i,j] = y[i,j]+t[i,j-1]
  }
}

##covariate x1 and x2 generated from Normal(0,0.5^2) and Bernoulli(0.5) respectively
x1=rnorm(n,0,0.5); x2=rbinom(n,1,0.5); x=cbind(x1,x2)

##true regression parameters and frailty variance parameter
beta1=1; beta2=-1; nu=0.5; 
parms=c(beta1,beta2)
phi=rgamma(n,nu,nu) 

##true baseline mean function
mu=function(t){2*t^(0.5)} 

##get the number of events between time intervals
z=matrix(,n,K);
xparms<-c();for (s in 1:nrow(x)){xparms[s]<-sum(x[s,]*parms)}
for (i in 1:n){
 z[i,1]<-rpois(1,mu(t[i,1])*exp(xparms[i])*phi[i]) 
 if (k[i]>1){
 z[i,2:k[i]]<-rpois(k[i]-1,(mu(t[i,2:k[i]])-mu(t[i,1:(k[i]-1)]))*exp(xparms[i])*phi[i])
 }
}

TestD<-list(t=t, x=x, z=z, k=k, K=K)

fit<-PCDReg.nf(DATA = TestD, order = 1, placement = TRUE, nknot=3,
               myknots, binit = c(-0.5,0.5) , ginit = seq(0.1,2),
               t.seq = seq(0,15,0.2), tol=10^(-3))
}

\keyword{Non-homogeneous Poisson process model}

