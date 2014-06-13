\name{BladTumor}
\alias{BladTumor}
\docType{data}
\title{
Bladder Tumor Cancer Data
}
\description{
Bladder tumor data were from the bladder cancer study conducted by the Veterans Administration Cooperative Urological Research Group, and used by many people to demonstrate methodology for recurrent event modelling. In this study 118 patients who had superficial bladder tumors were randomized into one of three treatment groups: placebo (48), thiotepa (38), and pyridoxine (32). During the study at each follow-up visit, new tumors since the last visit were counted, measured, and removed transurethrally. For each patient the initial number of tumors and the size of largest initial tumors were also recorded. For more details about this study see Byar et al. (1977).
}

\usage{data(BladTumor)}
\format{
  A data frame with 116 observations on the following 8 variables.
 \tabular{ll}{
  subject:\tab patient ID\cr
  time:\tab observation time\cr
  count:\tab cumulative number of tumors\cr
  number:\tab initial number of tumors (8=8 or more)\cr
  size:\tab size(cm) of largest initial tumors\cr
  pyridoxine:\tab dummy variable for pyridoxine treatment\cr
  thiotepa:\tab dummy variable for thiotepa treatment\cr
  count1:\tab number of new tumors since last observation time\cr
  }
}
\details{
This data include 116 subjects who have at least one follow-up observation after the study enrollment.
}

\note{
To further use all the functions of this package one must convert the original data structure into the specified data structre which is in a list form. For more details please see the following example using bladder tumor data. 
}


\source{
Wang, X. and Yan, J. (2011). Fitting semiparametric regressions for panel count survival data with an R package spef. \emph{Computer Methods and Programs in Biomedicine} 104,2 278-285
}
\seealso{
\code{\link{BladTumor1}}
}
\references{
Byar, D.P., Blackard,C., and the VACURG. (1977). Comparisons of placebo, pyridoxine, and topical thiotepa in preventing recurrence of stage I bladder cancer. \emph{Urology} 10, 556-561.
}
\examples{
data(BladTumor)

n<-max(BladTumor$subject)
#record the number of observations for all patients
k<-as.numeric(table(BladTumor$subject))   
K<-max(k)
t<-matrix(,n,K)
z<-matrix(,n,K)

x1<-c();x2<-c();x3<-c();x4<-c();

for (r in 1:n){
rownum<-which(BladTumor$subject==r)
#record all observation times
t[r,][1:k[r]]<-BladTumor[rownum,]$time    
#record all panel counts from non-overlapping intervals
z[r,][1:k[r]]<-BladTumor[rownum,]$count1  
x1[r]<-BladTumor[which(BladTumor$subject==r),]$number[1]
x2[r]<-BladTumor[which(BladTumor$subject==r),]$size[1]
x3[r]<-BladTumor[which(BladTumor$subject==r),]$pyridoxine[1];
x4[r]<-BladTumor[which(BladTumor$subject==r),]$thiotepa[1]
}

x<-cbind(x1,x2,x3,x4)
BladTumor1<-list(t=t,x=x,z=z,k=k,K=K)
}

\keyword{datasets}
