\name{GLSME}
\alias{GLSME}
\title{
Estimate regression parameters with correlated observations measurement errors.
}
\description{
The \code{GLSME} function estimates parameters of a linear model via generalized least squares. It allows 
for correlated predictors and responses. Furthermore it allows for correlated measurement errors
both in predictors and responses. The program specifically corrects for biase caused by these errors.
}
\usage{
GLSME(y, D, Vt, Ve, Vd, Vu, EstimateVariance = c(TRUE, TRUE), CenterPredictor = TRUE, 
InitialGuess = NULL, eps = 0.001, MaxIter = 50, MaxIterVar = 50, epsVar = 0.001, 
OutputType = "short", Vttype = NULL, Vetype = NULL, Vdtype = NULL, Vutype = NULL, 
ED = NULL, EDtype = "SingleValue")
}

\arguments{
  \item{y }{ A vector of observed response variables.}
  \item{D }{ 
a design matrix in which each column corresponds to a parameter to be estimated
in the \code{B}-vector/matrix. 
Each entry in these columns corresponds to a data point (e.g. a species
in comparative studies). The first column will typically be a column of ones, which
will estimate an intercept. Columns with indicators for categorical fixed effects can
also be added. Each regression variable is added as a column vector. The program
will automatically estimate one coefficient for each column in the design matrix and
these will be output in the order of the columns in the design matrix. Note that
columns corresponding to "random effects", indicated by nonzero entry in the \code{Vd}
matrix below, will be centered on their means unless the option 
\code{CenterPredictor = FALSE} is used to instruct the program to not do this.
If there is to be an INTERCEPT the user needs to put into \code{D} a constant column of 1s.}
  \item{Vt }{ The response biological residual covariance matrix (see Details).}
  \item{Ve }{ The response observation error covariance matrix (see Details).
observation errors in the response variable. In a comparative study in which the
response consists of species means, this will typically be a diagonal matrix with
squared standard errors of the means along the diagonal.}
  \item{Vd }{ Represents the true variance structure for the predictor variables. (see Details).}
  \item{Vu }{ The predictor observation variances (see Details)}
  \item{EstimateVariance }{ 
Option to turn off estimation of the variance parameters. This is a
vector of \code{TRUE} or \code{FALSE} values, 
where the first value corresponds to the true residual
variance, and the others correspond to the rest of the true predictor
variances. All the predictor variances can also be turned on or
off jointly by providing a single \code{TRUE} or \code{FALSE} value. 
The default is to estimate all
variance components. If a \code{FALSE} value is given the program assumes that the input
variance matrices are exact.}
  \item{CenterPredictor }{ \code{TRUE} or \code{FALSE} option to turn off automatic centering of predictors.}
  \item{InitialGuess }{ Starting value for the regression in the iterated GLS. The default is
\code{NULL}, in which case the program will use an OLS estimate. A specific starting
value can be given as a vector of numbers corresponding to each entry in the \code{B}-vector.
An additional number can also be given to specify the starting value of the residual
variance parameter.}
  \item{eps }{ tolerance for iterated GLS}
  \item{MaxIter }{ maximum number of iterations for iterated GLS}
  \item{MaxIterVar }{ maximum number of iterations for iterated GLS}
  \item{epsVar }{ tolerance for estimating variance parameters in predictors}
  \item{OutputType }{ should just the estimates be presented and their standard errors (\code{"short"}) or more detailed information (\code{"long"})}
  \item{Vttype }{ \code{Vt} matrix type (see Details)}
  \item{Vetype }{ \code{Ve} matrix type (see Details)}
  \item{Vdtype }{ \code{Vd} matrix type (see Details)}
  \item{Vutype }{ \code{Vu} matrix type (see Details)}
  \item{ED }{ the expected value of the design matrix, can be \code{NULL} then is estimated from the data}
  \item{EDtype }{ if \code{ED} is provided then specifies what is provided, allowed values are :
\itemize{
\item{\code{"constant"} }{ \code{ED} is a number and each value of \code{D} has mean equal to this number}
\item{\code{"variablemean"} }{ \code{ED} is a vector of length of number of variables, each value is a mean for the given predictor variable}
\item{\code{NULL} }{ \code{ED} is assumed to be calculated}}}
}

\details{
The code fits the general linear model with correlated data and
observation error in both dependent and independent variables. The code fits the
model 
\deqn{y = D\beta + r, r \sim N(0,V), V = \sigma^{2}T + V_{e} + Var[U\beta|D],}
where \code{y} is a vector of observed response variables, \code{D} is an observed design matrix, 
\eqn{\beta}
is a vector of parameters to be estimated, \code{T} is a matrix representing the true residual
variance up to a scale parameter, \eqn{\sigma^{2}}, that is estimated by the program, 
\eqn{V_{e}} is a matrix
of known measurement variance in the response variable, and \code{Var[U}\eqn{\beta}\code{|D]} 
is a matrix
representing effects of measurement error in the predictor variables (see appendix of
Hansen and Bartoszek 2012). To build the \code{Var[U}\eqn{\beta}\code{|D]} matrix, 
the program needs a known measurement
variance matrix \eqn{V_{u}} and a true variance matrix \eqn{V_{xt}} 
for each of the predictor variables
(these will be zero for fixed effects). The true variance matrices are assumed to be on
the form \eqn{V_{xt} = \sigma_{x}^{2}Sx} , where 
\code{Sx} is a matrix supplied by the user, and \eqn{\sigma_{x}^{2}} is a scale
parameter that the program estimates by maximum likelihood.

Note that this program cannot be used to fit parameters that enter nonlinearly
into the variance or the design matrix, as the \eqn{\alpha} 
in the adaptation-inertia model, but it
can be used to fit the other parameters in such models conditionally on given values
of the parameterized values of the matrices (and could hence be used as a subroutine
in a program for fitting such models).

Three important notes for the user :
\itemize{
\item{}{The program does NOT assume there will be an intercept -> hence the user needs to provide a column on \code{1}s in the design matrix if an intercept is desired.}
\item{}{The program by default centres predictors (controlled by \code{CenterPredictor}). This means that estimates of
fixed effects will be changed due to them absorbing the mean of the predictors. Using the centering has been
found to improve estimation especially of variance constants (\code{PredictorVarianceConstantEstimate} and \code{ResponseVarianceConstantEstimate} see Value).
The user should try out the option with \code{CenterPredictor} \code{TRUE} and \code{FALSE} (here fixed effects will not be effected) and compare results.}
\item{}{The program uses a Monte Carlo procedure as part of the estimation algorithm therefore the user should run the 
code a couple of times to see stability, and combine the results by e.g. a (weighted) average or
choose the best estimate according to e.g. the likelihood or \eqn{R^{2}}.}
}

The program tries to recognize the structure of the \code{Vt}, \code{Ve}, \code{Vd} and \code{Vu}
matrices passed (see the supplementary information to Hansen and Bartoszek 2012)
otherwise the user can specify how the matrix looks like in the appropriate matrix type variable, these can be 
in the respective \code{Vttype}, \code{Vetype}, \code{Vdtype} or \code{Vutype} parameter:
\itemize{
\item{\code{"SingleValue"} }{ the matrix variable is a single number that will be on the diagonal of the covariance matrix, used when the deviations are assumed to be uncorrelated and homoscedastic}
\item{\code{"Vector"} }{ the matrix variable is a vector each value corresponding to one of the variables and the covariance matrix will have that vector 
            appropriately on its diagonal, if an element of the list has the value \code{"F"} then this means that the variable is a fixed effect 
            and will get a 0 covariance matrix}
\item{\code{"CorrelatedPredictors"} }{ the matrix is a covariance matrix, it assumes that the observations are independent so the resulting covariance structure 
                is block diagonal, if some of the variables are fixed effects then in the matrix the values of the corresponding rows 
                and columns have to be 0 (this is a special case of BM with the second element equal to the identity matrix)}
\item{\code{"MatrixList"} }{ a list of length equal to the number of variables, each list element is the covariance structure 
                          for the given variable, if an element of the list has the value \code{"F"} then this means that the variable is a 
                         fixed effect and will get a 0 covariance matrix}
\item{\code{"BM"} }{ 
the matrix variable \code{Vx} is to be a list of two values, 
"Vx = Vx[[1]] %x%  Vx[[2]]" if e.g. this independent contrast 
        then the first value corresponds to the variable vector covariance while the second will be the matrix of distances between 
        species, if the first value is a number or vector then it is changed to a diagonal matrix, 
        if some of the variables are fixed effects then in the matrix of the first element of the list the values of the corresponding 
        rows and columns have to be 0
        }
\item{\code{NULL} or \code{"Matrix"} }{ the matrix is assumed calculated as given}
}
}

\value{
\itemize{
\item{GLSestimate }{ the GLS estimates without any correction (centering the predictors CHANGES fixed effects)}
\item{errorGLSestim }{ the estimates of their standard errors}
\item{BiasCorrectedGLSestimate }{ the bias corrected estimates (centering the predictors CHANGES fixed effects)}
\item{K }{ the bias attenuation factor matrix}
\item{R2 }{ \eqn{R^{2}} of the model with the GLS estimates not bias corrected}
\item{BiasCorrectedR2 }{ \eqn{R^{2}} of the model with the GLS estimates bias corrected}
\item{PredictorVarianceConstantEstimate }{ if \code{EstimateVariance[2]} is \code{TRUE} then the estimates of the unknown variance constants for the predictors otherwise not present}
\item{ResponseVarianceConstantEstimate }{ if \code{EstimateVariance[1]} is \code{TRUE} then the estimate of the unknown variance constant for the response otherwise not present}
if the \code{outputType} variable is set to \code{"long"} then the following additional fields will be in the output :
\item{CovarianceGLSestimate }{ estimate of the covariance matrix of the bias uncorrected GLS estimates}
\item{CovarianceBiasCorrectedGLSestimate }{ estimate of the covariance matrix of the bias corrected GLS estimates}
\item{response }{ the provided \code{y} vector}
\item{design }{ the provided design matrix \code{D}}
\item{Vt }{ the final used \code{Vt} matrix with the unknown variance constant incorporated (if estimated)}
\item{Ve }{ the final used \code{Ve} matrix} 
\item{Vd }{ the final used \code{Vd} matrix with the unknown variance constant(s) incorporated (if estimated)}
\item{Vu }{ the final used \code{Vu} matrix}
}
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

\author{Krzysztof Bartoszek}

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

\keyword{generalized least squares}
\keyword{regression}
\keyword{measurement error}
