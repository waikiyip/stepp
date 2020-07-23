\name{stepp.subpop}
\alias{stepp.subpop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	The constructor to create the stsubpop object and generate the subpopulations based on the
	specified stepp window and covariate of interest
}
\description{
  This is the constructor function to create a stepp subpopulation object. In addition, 
  it will also generate the stepp subpopulations based on the stepp window and covariate specified.
}
\usage{
stepp.subpop(swin, cov, coltype, coltrt, trts, minsubpops)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{swin}{the stepp window set up for the analysis}
  \item{cov}{the covariate of interest}
  \item{coltype}{variable providing the event type; required only for event-based windows}
  \item{coltrt}{variable providing the treatment indicators; required only for event-based windows}
  \item{trts}{vector containing the treatments list; required only for event-based windows}
  \item{minsubpops}{length-one numeric vector providing the minimum number of subpopulations to generate; required only for event-based windows}
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

# event-based windows using the BIG data set
data(bigKM)
rxgroup <- bigKM$trt
time    <- bigKM$time
evt     <- bigKM$event
cov     <- bigKM$ki67

swin_e <- new("stwin", type = "sliding_events", e1 = 10, e2 = 20)
subp_e <- new("stsubpop")
subp_e <- generate(subp_e, win = swin_e, covariate = cov, coltype = evt,
                   coltrt = rxgroup, trts = c(1, 2), minsubpops = 5)
summary(subp_e)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~subpopulation }
\keyword{ ~stepp }% __ONLY ONE__ keyword per line
