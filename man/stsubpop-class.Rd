\name{stsubpop-class}
\Rdversion{1.1}
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
Objects can be created by calls of the form \code{new("stsubpop") or the constructor method stepp.subpop}.
%%  ~~ describe objects here ~~ 
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
    \item{\code{init}:}{Object of class \code{"logical"} \cr
	a logical value indicating if the subpopulations are generated or not }
  }
}
\section{Methods}{
  \describe{
    \item{generate}{\code{signature(.Object = "stsubpop", win, covariate)}: \cr
	a method to generate the subpopulations based on the stepp window object and the specified covariate of interest }
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

}
\keyword{classes}
