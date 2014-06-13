\name{PCD.Corr}
\alias{PCD.Corr}

\title{
Pearson's coefficient of correlation between panel counts within two non-overlapping intervals
}
\description{
Calculating Pearson's coefficient of correlation between panel counts within two non-overlapping time intervals (t1, t2] and (t3, t4]. The arguments of this function are the output of function \code{\link{PCDReg.wf}}. For the corresponding formula see Yao, Wang and He (2014+).
}
\usage{
PCD.Corr(x, beta, nu, gamma, t1, t2, t3, t4, order, knots)
}

\arguments{
  \item{x}{
 the covariate vector.
}
\item{beta}{estimates of regression coefficients.}
\item{nu}{estimate of gamma frailty variance parameter.}
\item{gamma}{estimates of spline coefficients.}
\item{t1, t2}{interval endpoints. t1 must be less than t2.}
\item{t3, t4}{interval endpoints.  t3 must be less than t4.}
\item{order}{the order of basis functions.}
\item{knots}{knots used in the analysis.}
}

\value{
\item{corr.rho}{Pearson's coefficient of correlation.}
}

\note{
The two intervals (t1, t2] and (t3, t4] must not be overlapped.
}

\seealso{
\code{\link{PCDReg.wf}}
}

\references{
Yao, B., Wang, L., and He, X. (2014+). Semiparametric regression analysis of panel count data allowing for within-subject correlation.     
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
xparms=c();for (s in 1:nrow(x)){xparms[s]<-sum(x[s,]*parms)}
for (i in 1:n){
 z[i,1]<-rpois(1,mu(t[i,1])*exp(xparms[i])*phi[i]) 
 if (k[i]>1){
 z[i,2:k[i]]<-rpois(k[i]-1,(mu(t[i,2:k[i]])-mu(t[i,1:(k[i]-1)]))*exp(xparms[i])*phi[i])
 }
}

TestD<-list(t=t, x=x, z=z, k=k, K=K)

fit<-PCDReg.wf(DATA = TestD, order = 1, placement = TRUE, nknot = 3, myknots, 
               binit = c(0.5,-0.5), ninit = 0.1, ginit = seq(0.1,2),
               t.seq = seq(0,15,0.2), tol = 10^(-3))

x1=c(1,1);
b1=fit$beta; n1=fit$nu; g1=fit$gamma;
t1=0; t2=6; t3=6; t4=12;
order=1; knots=fit$knots;

PCD.Corr(x1, b1, n1, g1, t1, t2, t3, t4, order, knots)
}
\keyword{Pearson correlation}

