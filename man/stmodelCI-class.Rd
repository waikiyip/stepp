\name{stmodelCI-class}
\Rdversion{1.1}
\docType{class}
\alias{stmodelCI-class}
\alias{estimate,stmodelCI-method}
\alias{print,stmodelCI-method}
\alias{test,stmodelCI-method}

\title{Class \code{"stmodelCI"}}
\description{
	This is the stepp model of survival data with competing risks.
}
\section{Objects from the Class}{
	Objects can be created by calls of the form \code{new("stmodelCI", ...)} or by \cr
	the constructor function stepp.CI.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{coltrt}:}{Object of class \code{"numeric"} \cr
	the treatment variable}
    \item{\code{coltime}:}{Object of class \code{"numeric"} \cr
	the time to event variable}
    \item{\code{coltype}:}{Object of class \code{"numeric"} \cr
	variable with distinct codes for different causes of failure where coltype=0 for censored observations; coltype=1 for event of interest; coltype=2 for other causes of failure}
    \item{\code{trts}:}{Object of class \code{"numeric"} \cr
	a vector containing the codes for the 2 treatment groups, 1st and 2nd treatment groups, respectively}
    \item{\code{timePoint}:}{Object of class \code{"numeric"} \cr
	timepoint to estimate survival}
  }
}

\section{Extends}{
Class \code{"\linkS4class{stmodel}"}, directly.
}

\section{Methods}{
  \describe{
    \item{estimate}{\code{signature(.Object = "stmodelCI")}: \cr
	estimate the effect in absolute and relative scale of the overall and each subpopulation }
    \item{print}{\code{signature(.Object = "stmodelCI")}: \cr
	print the estimate, covariance matrices and statistics }
    \item{test}{\code{signature(.Object = "stmodelCI")}: \cr
	perform the permutation tests or GEE and obtain various statistics }
	 }
}

\value{
  The new method returns the stmodelCI object.

  The estimate method returns a list with the following fields:

  \item{model}{the stepp model - "CIe"}
  \item{sObs1}{a vector of effect estimates of all subpopulations based on the 1st treatment}
  \item{sSE1}{a vector of standard errors of effect estimates of all subpopulations based on the 1st treatment}
  \item{oObs1}{effect estimate of the entire population based on the 1st treatment}
  \item{oSE1}{the standard error of the effect estimate of the entire population based on the 1st treatment}
  \item{sObs2}{a vector of effect estimates of all subpopulations based on the 1st treatment}
  \item{sSE2}{a vector of standard errors of effect estimates of all subpopulations based on the 1st treatment}
  \item{oObs2}{effect estimate of the entire population based on the 1st treatment}
  \item{oSE2}{the standard error of the effect estimate of the entire population based on the 1st treatment}
  \item{skmw}{Wald's statistics for the effect estimate differences between the two treatments}
  \item{logHR}{a vector of log hazard ratio estimate of the subpopulations comparing 1st and 2nd treatments}
  \item{logHRSE}{a vector of standard error of the log hazard ratio estimate of the subpopulations comparing 1st and 2nd treatments}
  \item{ologHR}{the log hazard ratio estimate of the entire population comparing 1st and 2nd treatments}
  \item{ologHRSE}{the standard error of the log hazard ratio estimate of the entire population comparing 1st and 2nd treatments}
  \item{logHRw}{Wald's statistics for the log hazard ratio between the two treatments}

  The test method returns a list with the following fields:

  \item{model}{the stepp model - "CIt"}
  \item{sigma}{the covariance matrix for subpopulations based on effect differences}
  \item{hasigma}{the homogeneous association covariance matrix for subpopulations based on effect differences}
  \item{HRsigma}{the covariance matrix for the subpopulations based on hazard ratio}
  \item{haHRsigma}{the homogeneous association covariance matrix for subpopulations based on hazard ratio}
  \item{pvalue}{the supremum pvalue based on effect difference}
  \item{chi2pvalue}{the chisquare pvalue based on effect difference}
  \item{hapvalue}{the homogeneous association pvalue based on effect difference}
  \item{HRpvalue}{the supremum pvalue based on hazard ratio}
  \item{haHRpvalue}{the homogeneous association pvalue based on hazard ratio}  

}

\author{
  Wai-Ki Yip
}

\seealso{
  \code{\linkS4class{stwin}},     \code{\linkS4class{stsubpop}},  \code{\linkS4class{stmodelKM}}, 
  \code{\linkS4class{stmodelCI}}, \code{\linkS4class{stmodelGLM}}, 
  \code{\linkS4class{steppes}},   \code{\linkS4class{stmodel}}, 
  \code{\link{stepp.win}},        \code{\link{stepp.subpop}},     \code{\link{stepp.KM}},
                                  \code{\link{stepp.GLM}},
  \code{\link{stepp.test}},       \code{\link{estimate}},         \code{\link{generate}}

}
\examples{
showClass("stmodelCI")

##
n <- 1000		# set the sample size
mu <- 0			# set the mean and sd of the covariate
sigma <- 1

beta0 <- log(-log(0.5)) # set the intercept for the log hazard
beta1 <- -0.2		# set the slope on the covariate
beta2 <- 0.5		# set the slope on the treatment indicator
beta3 <- 0.7		# set the slope on the interaction

prob2 <- 0.2		# set the proportion type 2 events
cprob <- 0.3		# set the proportion censored

set.seed(7775432) 	# set the random number seed
covariate <- rnorm(n,mean=mu,sd=sigma)	# generate the covariate values
Txassign <- rbinom(n,1,0.5)		# generate the treatment indicator
x3 <- covariate*Txassign		# compute interaction term
					# compute the hazard for type 1 event
lambda1 <- exp(beta0+beta1*covariate+beta2*Txassign+beta3*x3)
lambda2 <- prob2*lambda1/(1-prob2)	# compute the hazard for the type 2 event
					# compute the hazard for censoring time
lambda0 <- cprob*(lambda1+lambda2)/(1-cprob)
t1 <- rexp(n,rate=lambda1)		# generate the survival time for type 1 event
t2 <- rexp(n,rate=lambda2)		# generate the survival time for type 2 event
t0 <- rexp(n,rate=lambda0)		# generate the censoring time
time <- pmin(t0,t1,t2)			# compute the observed survival time
type <- rep(0,n)
type[(t1 < t0)&(t1 < t2)] <- 1
type[(t2 < t0)&(t2 < t1)] <- 2

# create the stepp model object to analyze the data using Cumulative Incidence approach
x <- new ("stmodelCI", coltrt=Txassign, trts=c(0,1), coltime=time, coltype=type, timePoint=1.0)

}
\keyword{classes}
