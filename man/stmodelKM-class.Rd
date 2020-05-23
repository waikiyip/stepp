\name{stmodelKM-class}
\Rdversion{1.1}
\docType{class}
\alias{stmodelKM-class}
\alias{estimate,stmodelKM-method}
\alias{print,stmodelKM-method}
\alias{test,stmodelKM-method}

\title{Class \code{"stmodelKM"}}
\description{
	This is the S4 class for the stepp model of survival data using Kaplan-Meier method.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("stmodelKM", ...)} or by \cr
  the constructor function stmodel.KM.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{coltrt}:}{Object of class \code{"numeric"} \cr
	the treatment variable}
    \item{\code{survTime}:}{Object of class \code{"numeric"} \cr
	the time to event variable}
    \item{\code{censor}:}{Object of class \code{"numeric"} \cr
	the censor variable}
    \item{\code{trts}:}{Object of class \code{"numeric"} \cr
	a vector containing the codes for the 2 treatment groups, first and second treatment groups, respectively}
    \item{\code{timePoint}:}{Object of class \code{"numeric"} \cr
	timepoint to estimate survival}
  }
}
\section{Extends}{
Class \code{"\linkS4class{stmodel}"}, directly.
}
\section{Methods}{
  \describe{
    \item{estimate}{\code{signature(.Object = "stmodelKM")}:  \cr
	estimate the effect in absolute and relative scale of the overall population and each subpopulation.  }
    \item{print}{\code{signature(.Object = "stmodelKM")}: \cr
	print the estimate, covariance matrices and statistics. }
    \item{test}{\code{signature(.Object = "stmodelKM")}:  \cr
	perform the permutation tests or GEE and obtain various statistics. }
	 }
}

\value{
  The new method returns the stmodelKM object.

  The estimate method returns a list with the following fields:

  \item{model}{the stepp model - "KMe"}
  \item{sObs1}{a vector of effect estimates of all subpopulations based on the first treatment}
  \item{sSE1}{a vector of standard errors of effect estimates of all subpopulations based on the first treatment}
  \item{oObs1}{effect estimate of the entire population based on the first treatment}
  \item{oSE1}{the standard error of the effect estimate of the entire population based on the first treatment}
  \item{sObs2}{a vector of effect estimates of all subpopulations based on the group treatment}
  \item{sSE2}{a vector of standard errors of effect estimates of all subpopulations based on the first treatment}
  \item{oObs2}{effect estimate of the entire population based on the first treatment}
  \item{oSE2}{the standard error of the effect estimate of the entire population based on the first treatment}
  \item{skmw}{Wald's statistics for the effect estimate differences between the two treatments}
  \item{logHR}{a vector of log hazard ratio estimate of the subpopulations comparing first and second treatments}
  \item{logHRSE}{a vector of standard error of the log hazard ratio estimate of the subpopulations comparing first and second treatment}
  \item{ologHR}{the log hazard ratio estimate of the entire population comparing first and second treatment}
  \item{ologHRSE}{the standard error of the log hazard ratio estimate of the entire population comparing first and second treatment}
  \item{logHRw}{Wald's statistics for the log hazard ratio between the two treatment}

  The test method returns a list with the following fields:

  \item{model}{the stepp model - "KMt"}
  \item{sigma}{the covariance matrix for subpopulations based on effect differences}
  \item{hasigma}{the homogeneous association covariance matrix for subpopulations based on effect differences}
  \item{HRsigma}{the covariance matrix for the subpopulations based on hazard ratios}
  \item{haHRsigma}{the homogeneous association covariance matrix for subpopulations based on hazard ratios}
  \item{pvalue}{the supremum pvalue based on effect difference}
  \item{chi2pvalue}{the chisquare pvalue based on effect difference}
  \item{hapvalue}{the homogeneous association pvalue based on effect difference}

}

\author{
Wai-Ki YIp
}

\seealso{

  \code{\linkS4class{stwin}},     \code{\linkS4class{stsubpop}},  
  \code{\linkS4class{stmodelCI}}, \code{\linkS4class{stmodelGLM}}, 
  \code{\linkS4class{steppes}},   \code{\linkS4class{stmodel}}, 
  \code{\link{stepp.win}},        \code{\link{stepp.subpop}},     \code{\link{stepp.KM}},
  \code{\link{stepp.CI}},         \code{\link{stepp.GLM}},
  \code{\link{stepp.test}},       \code{\link{estimate}},         \code{\link{generate}}

}
\examples{
showClass("stmodelKM")

#GENERATE TREATMENT VARIABLE:
N <- 1000
Txassign <- sample(c(1,2), N, replace=TRUE, prob=c(1/2, 1/2))
n1 <- length(Txassign[Txassign==1])
n2 <- N - n1
#GENERATE A COVARIATE:
covariate <- rnorm(N, 55, 7)
#GENERATE SURVIVAL AND CENSORING VARIABLES ASSUMING A TREATMENT COVARIATE INTERACTION:
Entry <- sort( runif(N, 0, 5) )
SurvT1 <- .5
beta0 <-  -65 / 75
beta1 <- 2 / 75
Surv <- rep(0, N)
lambda1 <- -log(SurvT1) / 4
Surv[Txassign==1] <- rexp(n1, lambda1)
Surv[Txassign==2] <- rexp(n2, (lambda1*(beta0+beta1*covariate[Txassign==2])))
EventTimes <- rep(0, N)
EventTimes <- Entry + Surv
censor <- rep(0, N)
time <- rep(0,N)
for ( i in 1:N )
   {
     censor[i] <- ifelse( EventTimes[i] <= 7, 1, 0 )
     time[i] <- ifelse( EventTimes[i] < 7, Surv[i], 7 - Entry[i] ) 
   }

modKM <- new("stmodelKM", coltrt=Txassign, survTime=time, censor=censor, trts=c(1,2), timePoint=4)

}
\keyword{classes}
