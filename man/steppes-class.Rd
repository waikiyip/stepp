\name{steppes-class}
\Rdversion{1.1}
\docType{class}
\alias{steppes-class}
\alias{estimate,steppes-method}
\alias{plot,steppes-method}
\alias{print,steppes-method}
\alias{summary,steppes-method}
\alias{test,steppes-method}

\title{Class \code{"steppes"}}
\description{
	This is the stepp object for stepp results.  It keeps track of all the estimates, covariance matrices and
	the test statistics
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("steppes", ...)} or by the constructor function
  stepp.test.
}
\section{Slots}{
  \describe{
    \item{\code{subpop}:}{Object of class \code{"stsubpop"} \cr
	the filled stepp subpopulation object }
    \item{\code{model}:}{Object of class \code{"stmodel"} \cr
	the model of the data to be used for stepp analysis }
    \item{\code{effect}:}{Object of class \code{"ANY"} \cr
	effect estimates of each subpopulation in absolute and relative scales \cr
	list return from the estimate methods of various stepp models \cr
	see documentation in various stepp models for details on what is in the list}
    \item{\code{result}:}{Object of class \code{"ANY"} \cr
	permutation and gee test results and various covariance matrices \cr
	list return from the test methods of various stepp models \cr
	see documentation in various stepp models for details on what is in the list}
    \item{\code{nperm}:}{Object of class \code{"numeric"} \cr
	number of permutation for the permutation test; default is 2500 }
  }
}
\section{Methods}{
  \describe{
    \item{estimate}{\code{signature(.Object = "steppes", sp, model)}: \cr
	estimate the effect in absolute and relative scale of the overall and each subpopulation for the specified subpopulation using the specified model. }
    \item{plot}{\code{signature(.Object = "steppes")}: \cr
	generate the three stepp plots}
    \item{print}{\code{signature(.Object = "steppes", estimate=TRUE, cov=TRUE, test=TRUE)}: \cr
	print the estimates, covariance matrices and statistics from the stepp analysis }
    \item{summary}{\code{signature(.Object = "steppes")}: \cr
	produce a summary of the steppes object }
    \item{test}{\code{signature(.Object = "steppes", nperm=100, showstatus=TRUE)}: \cr
	perform the permutation tests and obtain various statistics. Number of permutation is default to 100 and status bar to show.}
  }
}

\value{
  The new method returns the steppes object.
}

\note{
  The plot function for steppes class generates 3 plots at once.  If you are using R studio for development, 
  you can only get the last plot as R studio only allows one graphical device.  There is a simple workaround.
  By writing the plots to pdf files in your working directory, you can display them outside of R studio. 
}

\author{
  Wai-Ki Yip
}

\seealso{
  \code{\linkS4class{stwin}},     \code{\linkS4class{stsubpop}},  \code{\linkS4class{stmodelKM}}, 
  \code{\linkS4class{stmodelCI}}, \code{\linkS4class{stmodelGLM}},\code{\linkS4class{stmodel}}, 
  \code{\link{stepp.win}},        \code{\link{stepp.subpop}},     \code{\link{stepp.KM}},
  \code{\link{stepp.CI}},         \code{\link{stepp.GLM}},
  \code{\link{stepp.test}},       \code{\link{estimate}},         \code{\link{generate}}

}
\examples{
showClass("steppes")
}
\keyword{classes}
