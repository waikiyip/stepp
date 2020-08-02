\name{stsubpop-class}
\Rdversion{1.3}
\docType{class}
\alias{stsubpop-class}
\alias{generate,stsubpop-method}
\alias{summary,stsubpop-method}

\title{Class \code{"stsubpop"}}
\description{
	This is the S4 class for stepp subpopulation object.  The subpopulations are generated
  based on the stepp windows and the covariate of interest.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("stsubpop")} or the constructor method stepp.subpop.
}
\section{Slots}{
  \describe{
    \item{\code{win}:}{Object of class \code{"stwin"} \cr
	the stepp window set up for the analysis }
    \item{\code{colvar}:}{Object of class \code{"numeric"} \cr
	the covariate of interest }
    \item{\code{nsubpop}:}{Object of class \code{"numeric"} \cr
	the number of subpopulations generated }
    \item{\code{subpop}:}{Object of class \code{"ANY"} \cr
	a matrix of subpopulations generated based on the stepp window and the specified covariate of interest }
    \item{\code{npatsub}:}{Object of class \code{"numeric"} 
	a vector of size of each subpopulation }
    \item{\code{medianz}:}{Object of class \code{"numeric"} \cr
	a vector of median value of the covariate of interest for each subpopulation }
    \item{\code{minz}:}{Object of class \code{"numeric"} \cr
	a vector of rounded minimum value of the covariate of interest for each subpopulation }
    \item{\code{minc}:}{Object of class \code{"numeric"} \cr
	a vector of the minimum value of the covariate of interest for each subpopulation }
     \item{\code{maxz}:}{Object of class \code{"numeric"} \cr
	a vector of rounded maximum value of the covariate of interest for each subpopulation }
    \item{\code{maxc}:}{Object of class \code{"numeric"} \cr
	a vector of the maximum value of the covariate of interest for each subpopulation }
    \item{\code{neventsubTrt0}:}{Object of class \code{"numeric"} or \code{NULL} \cr
  a vector containing the number of events in each subpopulation for the baseline treatment group}
    \item{\code{neventsubTrt1}:}{Object of class \code{"numeric"} or \code{NULL} \cr
  a vector containing the number of events in each subpopulation for the active treatment group}
    \item{\code{init}:}{Object of class \code{"logical"} \cr
	a logical value indicating if the subpopulations have already been generated or not }
  }
}
\section{Methods}{
  \describe{
    \item{generate}{\code{signature(.Object = "stsubpop", win, covariate, coltype, coltrt, trts, minsubpops)}: \cr
	a method to generate the subpopulations based on the stepp window object and the specified covariate of interest. For event-based windows, also the event type (\code{coltype}), treatment indicator (\code{coltrt}), treatments list (\code{trts}) and minimum number of subpopulations (\code{minsubpops}) must be provided }
    \item{summary}{\code{signature(.Object = "stsubpop")}: \cr
	a method to display the summary of the subpopulations generated }
	 }
}

\author{
  Wai-Ki Yip
}

\seealso{
  \code{\linkS4class{stwin}},                                     \code{\linkS4class{stmodelKM}}, 
  \code{\linkS4class{stmodelCI}}, \code{\linkS4class{stmodelGLM}}, 
  \code{\linkS4class{steppes}},   \code{\linkS4class{stmodel}}, 
  \code{\link{stepp.win}},        \code{\link{stepp.subpop}},     \code{\link{stepp.KM}},
  \code{\link{stepp.CI}},         \code{\link{stepp.GLM}},
  \code{\link{stepp.test}},       \code{\link{estimate}},         \code{\link{generate}}

}
\examples{
showClass("stsubpop")

# create a steppp window
win1 <- stepp.win(type="sliding", r1=5,r2=10)

# generate the covariate of interest
Y <- rnorm(100)

# create and generate the stepp subpopulation
sp <- new("stsubpop")
sp <- generate(sp, win=win1, cov=Y)
summary(sp)

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
\keyword{classes}
