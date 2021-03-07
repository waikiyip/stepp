\name{stwin-class}
\Rdversion{1.3}
\docType{class}
\alias{stwin-class}
\alias{summary,stwin-method}

\title{Class \code{"stwin"}}
\description{
  This the S4 class for the stepp window object.  The stepp window class describes the way
  to set up the subpopulation for a stepp analysis.  It specifies a pattern of 
  subpopulation you would like   to explore. There are three kinds of patterns: "sliding window", 
  "event-based sliding window" and "tail-oriented window".
 
  These patterns are specified through the following set of values (r1 and r2 for sliding and
  tail-oriented windows, e1 and e2 for event-based sliding windows):\cr
  \cr
  1. for sliding windows based on the number of patients ("sliding"), r1 is the minimum number
     of patients allowed in overlapping windows and r2 is the approximate size of subpopulation
     in each window;\cr
  \cr
  2. for sliding windows based on the number of events ("sliding_events"), e1 is the minimum
     number of events allowed in overlapping windows and e2 is the approximate size of
     subpopulation in each window in terms of number of events;\cr
  \cr
  3. for tail-oriented window, r1 is a vector of maximum covariate value for each subpopulation
     from the minimum value of the entire subpopulation, and r2 is a vector of minimum covariate
     values for each subpopulation from the maximum value of the window. The utility function 
     \code{\link{gen.tailwin}()} can be used to generate these vectors based on the number of desired 
     subpopulations.\cr
  \cr
  When "sliding_events" is chosen r1 and r2 are set to \code{NULL}, while when "sliding" or
  "tail-oriented" are chosen e1 and e2 are set to \code{NULL}.\cr
  \cr
  These pattern are based on all patients.
}

\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("stwin", type="sliding", r1=5, r2=20)} or
  the constructor function \code{stepp.win}.
}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"character"} \cr
		stepp window type; \code{"sliding"}, \code{"sliding_events"} or \code{"tail-oriented"}}
    \item{\code{r1}:}{Object of class \code{"numeric"} or \code{NULL} \cr 
		sliding window: minimum number of patients allowed in overlapping windows;\cr
	        tail-oriented window: a vector of maximum covariate value for each subpopulation}
    \item{\code{r2}:}{Object of class \code{"numeric"} or \code{NULL} \cr
		sliding window: size of subpopulation in each window;\cr
	    	tail-oriented window: a vector of minimum covariate value for each subpopulation}
    \item{\code{e1}:}{Object of class \code{"numeric"} or \code{NULL} \cr 
    event-based sliding window: minimum number of events allowed in overlapping windows}
    \item{\code{e2}:}{Object of class \code{"numeric"} or \code{NULL} \cr
    event-based sliding window: number of events in each subpopulation}
  }
}
\section{Methods}{
  \describe{
    \item{summary}{\code{signature(object = "stwin")}: \cr
		print a summary of the stepp windows object }
	 }
}

\author{
  Wai-Ki Yip
}

\seealso{
                                  \code{\linkS4class{stsubpop}},  \code{\linkS4class{stmodelKM}}, 
  \code{\linkS4class{stmodelCI}}, \code{\linkS4class{stmodelGLM}}, 
  \code{\linkS4class{steppes}},   \code{\linkS4class{stmodel}}, 
  \code{\link{stepp.win}},        \code{\link{stepp.subpop}},     \code{\link{stepp.KM}},
  \code{\link{stepp.CI}},         \code{\link{stepp.GLM}},
  \code{\link{stepp.test}},       \code{\link{estimate}},         \code{\link{generate}}

}
\examples{
  showClass("stwin")

  # create a stepp window of type "sliding" with (r2) size of subpopulation
  # in each window to be 200 and (r1) allows only 50 patients in the
  # overlapping windows
  mywin <- new("stwin", type="sliding", r1=50, r2=200)

  # print a summary of the stepp window object created
  summary(mywin)

  # create a stepp window object of type "sliding_events", 
  # (event-based) subpopulation size is 200 and allows
  # only 50 events between overlapping windows
  mywin <- stepp.win(type="sliding_events", e1=50, e2=200)

  # print a summary of the stepp window object
  summary(mywin) 
}
\keyword{classes}
