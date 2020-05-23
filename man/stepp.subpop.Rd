\name{stepp.subpop}
\alias{stepp.subpop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	The constructor to create the stsubpop object and generate the subpopulations based on the
	specified stepp window and covariate of interest
}
\description{
  This is the constructor function to create a stepp subpopulation object.  In addition, 
  it will also generate the stepp subpopulations based on the stepp window and covariate specified.
}
\usage{
stepp.subpop(swin, cov)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{swin}{the stepp window set up for the analysis}
  \item{cov}{the covariate of interest}
}

\value{
  It returns the stsubpop object with subpopulations generated. 
}

\author{
  Wai-Ki Yip
}

\seealso{
  \code{\linkS4class{stwin}},     \code{\linkS4class{stsubpop}},  \code{\linkS4class{stmodelKM}}, 
  \code{\linkS4class{stmodelCI}}, \code{\linkS4class{stmodelGLM}}, 
  \code{\linkS4class{steppes}},   \code{\linkS4class{stmodel}}, 
  \code{\link{stepp.win}},                                        \code{\link{stepp.KM}},
  \code{\link{stepp.CI}},         \code{\link{stepp.GLM}},
  \code{\link{stepp.test}},       \code{\link{estimate}},         \code{\link{generate}}

}
\examples{
# create a steppp window
win1 <- stepp.win(type="sliding", r1=5,r2=10)

# generate the covariate of interest
Y <- rnorm(100)

# create and generate the stepp subpopulation
sp   <- stepp.subpop(swin=win1, cov=Y)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
