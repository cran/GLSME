\name{GLSME-package}
\alias{GLSME-package}
\docType{package}
\title{Generalized Least Squares with Measurement Error}
\description{
The code fits the general linear model with correlated data and
observation error in both dependent and independent variables. The code fits the
model 
\deqn{y = D\beta + r, r \sim N(0,V), V = \sigma^{2} T + V_{e} + Var[U\beta|D],}
where \code{y} is a vector of observed response variables, \code{D} is an observed design matrix, 
\eqn{\beta}
is a vector of regression parameters to be estimated, 
\eqn{\sigma^{2}}\code{T} is a matrix representing the true residual
variance, \eqn{V_{e}} is a matrix
of known measurement variance in the response variable, and \code{Var[U}\eqn{\beta}\code{|D]} is a matrix
representing effects of measurement error in the predictor variables (see Hansen and Bartoszek 2012).
}
\details{
\tabular{ll}{
Package: \tab GLSME\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-05-19\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
The code fits the general linear model with correlated data and
observation error in both dependent and independent variables. The code fits the
model 
\deqn{y = D\beta + r, r \sim N(0,V), V = \sigma^{2}T + V_{e} + Var[U\beta|D],}
where \code{y} is a vector of observed response variables, \code{D} is an observed design matrix, 
\eqn{\beta}
is a vector of regression parameters to be estimated, 
\eqn{\sigma^{2}}\code{T} is a matrix representing the true residual
variance, \eqn{V_{e}} is a matrix
of known measurement variance in the response variable, and \code{Var[U}\eqn{\beta}\code{|D]} is a matrix
representing effects of measurement error in the predictor variables (see Hansen and Bartoszek 2012).

The estimation function is \code{GLSME}. It is an iterated (if the variance parameters are unknown)
generalized least squares estimation procedure. 

The motivation for the approach is that the observations and errors are correlated due
to an underlying phylogeny but the program allows for any dependence structure. 

In the \code{\link[mvSLOUCH]{mvSLOUCH}} package an alternative method of correcting 
for observation error is used. The error variance-covariance matrix enters
the likelihood function by being added to the biological variance-covariance matrix.
}
\author{
Krzysztof Bartoszek
Maintainer: <bartoszekkj@gmail.com>
}
\references{
Bartoszek, K. and Pienaar, J. and Mostad. P. and Andersson, S. and Hansen, T. F. (2012)
  A phylogenetic comparative method for studying  multivariate adaptation.
  Journal of Theoretical Biology 314:204-215.
  
Hansen, T.F. (1997) 
  Stabilizing selection and the comparative analysis of adaptation. 
  Evolution 51:1341-1351.

Hansen, T.F. and Bartoszek, K. (2012)
  Interpreting the evolutionary regression: the interplay between observational and biological errors in phylogenetic comparative studies.
  Systematic Biology 61(3):413-425.

Hansen, T.F. and Pienaar, J. and Orzack, S.H. (2008) 
  A comparative method for studying adaptation to randomly evolving environment. 
  Evolution 62:1965-1977.

}

\keyword{generalized least squares}
\keyword{regression}
\keyword{measurement error}

\seealso{
\code{\link[mvSLOUCH]{mvSLOUCH}} 
}
\examples{
library(mvSLOUCH)
library(ape)
library(ouch)
n<-5 ## number of species
apetree<-rtree(n)
phyltree<-ape2ouch(apetree) ##mvslouch requires ouch format
### Correct the names of the internal node labels.
phyltree@nodelabels[1:(phyltree@nnodes-phyltree@nterm)]<-
as.character(1:(phyltree@nnodes-phyltree@nterm))
### Define Brownian motion parameters to be able to simulate data under the Brownian motion model.
BMparameters<-list(vX0=matrix(0,nrow=2,ncol=1),Sxx=rbind(c(1,0),c(0.2,1)))
### Now simulate the data and remove the values corresponding to the internal nodes.
xydata<-simulBMProcPhylTree(phyltree,X0=BMparameters$vX0,Sigma=BMparameters$Sxx)
xydata<-xydata[(nrow(xydata)-n+1):nrow(xydata),]

x<-xydata[,1]
y<-xydata[,2]

yerror<-diag((rnorm(n,mean=0,sd=0.1))^2) #create error matrix
y<-rmvnorm(1,mean=y,sigma=yerror)[1,]
xerror<-diag((rnorm(n,mean=0,sd=0.1))^2) #create error matrix
x<-rmvnorm(1,mean=x,sigma=xerror)[1,]
GLSME(y=y, CenterPredictor=TRUE, D=cbind(rep(1, n), x), Vt=vcv(apetree), 
Ve=yerror, Vd=list("F",vcv(apetree)), Vu=list("F", xerror))
}
