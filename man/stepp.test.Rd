\name{stepp.test}
\alias{stepp.test}
%- ALSO NEED an '\alias' for EACH other topic documented here.
\title{
        The constructor to generate a complete steppes object with effect estimates and test statistics
}
\description{
  This is a constructor function for the steppes object.   In addition, it estimates all the effects of
  each subpopulations, performs permutation tests and generates all covariance matrices and statistics.
}
\usage{
  stepp.test(subpop, model, nperm, showstatus = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{subpop}{the filled stepp subpopulation object}
  \item{model}{the model of the data to be used for stepp analysis}
  \item{nperm}{number of permutation used in the permutation test}
  \item{showstatus}{display the progress bar for the permutation test; default is TRUE}
}

\value{
  It returns a steppes object with all effects estimates, covariance matrices and statistics.
}

\details{
  Permutation tests for all the statistics are done (see ref below).  For best results, considering
  using 2500 permutations to obtain a rich distribution from which to draw inference.\cr
}

\author{
  Wai-ki Yip
}
\references{
  Bonetti M, Zahrieh D, Cole BF, Gelber RD. A small sample study of the STEPP approach to assessing treatment-covariate interactions in survival data. Statistics in Medicine 2009; 28(8):1255-68.

}
\seealso{
  \code{\linkS4class{stwin}},     \code{\linkS4class{stsubpop}},  \code{\linkS4class{stmodelKM}},
  \code{\linkS4class{stmodelCI}}, \code{\linkS4class{stmodelGLM}},
  \code{\linkS4class{steppes}},    \code{\linkS4class{stmodel}},
  \code{\link{stepp.win}},        \code{\link{stepp.subpop}},     \code{\link{stepp.KM}},
  \code{\link{stepp.CI}},         \code{\link{stepp.GLM}},
                                  \code{\link{estimate}},         \code{\link{generate}}
}
\examples{
data(bigCI)

rxgroup <- bigCI$trt
time    <- bigCI$time
evt     <- bigCI$event
cov     <- bigCI$ki67

#
# using constructor functions
swin    <- stepp.win(type="sliding", r1=50, r2=150)
subp    <- stepp.subpop(swin=swin, cov=cov)
summary(subp)

smodel  <- stepp.CI(coltrt=rxgroup, trts=c(1,2), coltime=time, coltype=evt, timePoint=4)
# Warning: In this example, the permutations have been set to 0 to allow the function
# to finish in a short amount of time.  IT IS RECOMMEND TO USE AT LEAST 2500 PERMUTATIONS TO 
# PROVIDE STABLE RESULTS.
statCI  <- stepp.test(subpop=subp, model=smodel, nperm=0)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ Keyword per line

