\name{stepp.CI}
\alias{stepp.CI}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
	The constructor to create the stmodelCI object
}

\description{
  This is the constructor function for the stmodelCI object. This object sets up the data with
  a stepp model using competing risks method for analysis. (CI stands for Cumulative Incidence.)\cr

  The model explores the treatment-effect interactions in competing risks data arising 
  from two or more treatment arms of a clinical trial. A permutation distribution approach to inference 
  is implemented that permutes covariate values within a treatment group. 
  The statistical significance of observed heterogeneity of treatment effects is calculated using 
  permutation tests:\cr\cr
  1) for the maximum difference between each subpopulation effect and the overall population 
     treatment effect or supremum based test statistic; \cr
  2) for the difference between each subpopulation effect and the overall population treatment 
     effect, which resembles the chi-square statistic. 
}

\usage{
	stepp.CI(coltrt, coltime, coltype, trts, timePoint)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coltrt}{the treatment variable}
  \item{coltime}{the time to event variable}
  \item{coltype}{variable with distinct codes for different causes of failure where coltype=0 for censored observations; coltype=1 for event of interest; coltype=2 for other causes of failure}
  \item{trts}{a vector containing the codes for the 2 treatment groups, 1st and 2nd treatment groups, respectively}
  \item{timePoint}{timepoint to estimate survival}
}

\value{
	It returns the stmodelCI object.
}

\author{
Wai-Ki Yip
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\references{
  Bonetti M, Gelber RD. Patterns of treatment effects in subsets of patients in clinical trials. Biostatistics 2004; 5(3):465-481.

  Bonetti M, Zahrieh D, Cole BF, Gelber RD. A small sample study of the STEPP approach to assessing treatment-covariate interactions in survival data. Statistics in Medicine 2009; 28(8):1255-68. 

  Lazar AA, Cole BF, Bonetti M, Gelber RD.  Evaluation of treatment-effect heterogeneity usiing biomarkers measured
  on a continuous scale: subpopulation treatment effect pattern plot.  Journal of Clinical Oncology, 2010; 28(29): 4539-4544.
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
x <- stepp.CI(coltrt=Txassign, trts=c(0,1), coltime=time, coltype=type, timePoint=1.0)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
