\name{stepp.KM}
\alias{stepp.KM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	The constructor to create the stmodelKM object
}
\description{
  This is the constructor function for the stmodelKM object.  This object sets up the data with 
  a stepp model using the Kaplan-Meier method for analysis.\cr

  The model explores the treatment-covariate interactions in survival data arising 
  from two treatment arms of a clinical trial. The treatment effects are measured using survival 
  functions at a specified time point estimated from the Kaplan-Meier method and the hazard ratio 
  based on observed-minus-expected estimation. A permutation distribution approach to inference 
  is implemented, based on permuting the covariate values within each treatment group. 
  The statistical significance of observed heterogeneity of treatment effects is calculated using 
  permutation tests:\cr\cr
  1) for the maximum difference between each subpopulation effect and the overall population 
     treatment effect or supremum based test statistic; \cr
  2) for the difference between each subpopulation effect and the overall population treatment 
     effect, which resembles the chi-square statistic. 
}
\usage{
	stepp.KM(coltrt, survTime, censor, trts, timePoint)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coltrt}{the treatment variable}
  \item{survTime}{the time to event variable}
  \item{censor}{the censor variable}
  \item{trts}{a vector containing the codes for the 2 treatment arms, 1st and 2nd treatment groups, respectively}
  \item{timePoint}{timepoint to estimate survival}
}

\value{
	It returns the stmodelKM object.
}

\author{
Wai-Ki Yip
}

\references{
  Bonetti M, Gelber RD. Patterns of treatment effects in subsets of patients in clinical trials. Biostatistics 2004; 5(3):465-481.

  Bonetti M, Zahrieh D, Cole BF, Gelber RD. A small sample study of the STEPP approach to assessing treatment-covariate interactions in survival data. Statistics in Medicine 2009; 28(8):1255-68. 

}

\seealso{
  \code{\linkS4class{stwin}},     \code{\linkS4class{stsubpop}},  \code{\linkS4class{stmodelKM}}, 
  \code{\linkS4class{stmodelCI}}, \code{\linkS4class{stmodelGLM}}, 
  \code{\linkS4class{steppes}},   \code{\linkS4class{stmodel}}, 
  \code{\link{stepp.win}},        \code{\link{stepp.subpop}},     
  \code{\link{stepp.CI}},         \code{\link{stepp.GLM}},
  \code{\link{stepp.test}},       \code{\link{estimate}},         \code{\link{generate}}

}
\examples{
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

modKM <- stepp.KM( coltrt=Txassign, survTime=time, censor=censor, trts=c(1,2), timePoint=4)

}