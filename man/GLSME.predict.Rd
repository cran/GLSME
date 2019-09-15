\name{GLSME.predict}
\alias{GLSME.predict}
\title{Prediction for a new observation using parameters estimated by the \code{GLSME} function}
\description{
The function takes parameters estimated by the \code{GLSME} function and predicts the response
for a new observation of predictors. It also returns confidence intervals on the prediction.
The function is still under development.
}
\usage{
GLSME.predict(xo, glsme.estimate, vy, vx, alpha = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xo }{ The new observed predictors. In a intercept is in the model then a 1 has to be included for it.}
  \item{glsme.estimate }{ The output of the \code{GLSME} function. Has to have format \code{"long"}.}
  \item{vy }{ Residual variance, both biological and measurement error.}
  \item{vx }{ Biological variance in predictor, NOT observation variance of predictor. If there is a predictor in the model then a 0 row and column have to included for it.}
  \item{alpha }{ Level for confidence interval.}
}
%%\details{
%%}
\value{
  \item{BiasCorr }{ Prediction using the bias corrected estimate.
  \itemize{
	\item{prediction }{ Predicted value}
	\item{MSE }{ Estimate of mean square error. They are calculated by the formula
	\deqn{v_{y}+x_{o}^{T}(MSE[\beta | XO])x_{o}+ \beta'^{T} v_{x} \beta',}
	where \eqn{\beta'} is the bias corrected estimate of \eqn{\beta}.
	}
	\item{CI }{ \eqn{1-\alpha} level confidence intervals. They are calculated by the formula
	\deqn{\sqrt{1+1/n}*t_{\alpha}*(v_{y}+ \beta'^{T} v_{x} \beta'),}
	where \eqn{t_{\alpha}} is the \eqn{1-\alpha/2} level quantile of the t-distribution with n-k degrees of freedom,
	k is the number of regression parameters to estimate, 
	\eqn{\beta'} is the bias corrected estimate of \eqn{\beta} and n is the sample size used in 
	the estimation.}
    }
  }
  \item{BiasUncorr }{ Prediction using the bias uncorrected estimate.
  \itemize{
	\item{prediction }{ Predicted value}
	\item{MSE }{ Estimate of mean square error. They are calculated by the formula
	\deqn{v_{y}+x_{o}^{T}(MSE[\beta | XO])x_{o}+ \beta'^{T} v_{x} \beta',}
	where \eqn{\beta'} is the bias uncorrected estimate of \eqn{\beta}.
	}
	\item{CI }{ \eqn{1-\alpha} level confidence intervals. They are calculated by the formula
	\deqn{\sqrt{1+1/n}*t_{\alpha}*(v_{y}+ \beta'^{T} v_{x} \beta'),}
	where \eqn{t_{\alpha}} is the \eqn{1-\alpha/2} level quantile of the t-distribution with n-k degrees of freedom,
	k is the number of regression parameters to estimate, \eqn{\beta'} is the bias uncorrected estimate of \eqn{\beta} 
	and n is the sample size used in the estimation.}
    }
  }
 
}
\references{
Hansen, T.F. and Bartoszek, K. (2012)
  Interpreting the evolutionary regression: the interplay between observational and biological errors in phylogenetic comparative studies.
  Systematic Biology 61(3):413-425.
}

\author{Krzysztof Bartoszek}

\examples{
set.seed(12345)
n<-3 ## number of species
apetree<-ape::rtree(n)
### Define Brownian motion parameters to be able to simulate data under the Brownian motion model.
BMparameters<-list(vX0=matrix(0,nrow=2,ncol=1),Sxx=rbind(c(1,0),c(0.2,1)))
### Now simulate the data and remove the values corresponding to the internal nodes.
xydata<-mvSLOUCH::simulBMProcPhylTree(apetree,X0=BMparameters$vX0,Sigma=BMparameters$Sxx)
xydata<-xydata[(nrow(xydata)-n+1):nrow(xydata),]

x<-xydata[,1]
y<-xydata[,2]

yerror<-diag((stats::rnorm(n,mean=0,sd=0.1))^2) #create error matrix
y<-mvtnorm::rmvnorm(1,mean=y,sigma=yerror)[1,]
xerror<-diag((stats::rnorm(n,mean=0,sd=0.1))^2) #create error matrix
x<-mvtnorm::rmvnorm(1,mean=x,sigma=xerror)[1,]
glsme.res<-GLSME(y=y, CenterPredictor=TRUE, D=cbind(rep(1, n), x), Vt=ape::vcv(apetree), 
Ve=yerror, Vd=list("F",ape::vcv(apetree)), Vu=list("F", xerror),OutputType="long")
GLSME.predict(c(1,1), glsme.res, vy=1, vx=rbind(c(0,0),c(0,1)))
}
\keyword{prediction}
\keyword{generalized least squares}
\keyword{regression}
\keyword{measurement error}
