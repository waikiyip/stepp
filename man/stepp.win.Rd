\name{stepp.win}
\alias{stepp.win}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  The constructor to create the stepp window object
}

\description{
  This is the constructor function to create a stepp window object.
}

\usage{
  stepp.win(type = "sliding", r1 = 5, r2 = 20, e1 = NULL, e2 = NULL, basedon = "all")
}

\arguments{
  \item{type}{type of stepp window; "sliding", "sliding_events" or "tail-oriented"}
  \item{r1}{sliding window: minimum number of patients allowed in overlapping windows;\cr
	    tail-oriented window: a vector of maximum covariate values for each subpopulation}
  \item{r2}{sliding window: size of subpopulation in each window;\cr
	    tail-oriented window: a vector of minimum covariate values for each subpopulation}
  \item{e1}{sliding window: minimum number of events allowed in overlapping windows}
  \item{e2}{sliding window: number of events in each subpopulation}
  \item{basedon}{what the window is based on - "all" (default)"}
}

\details{
  This is the functional interface to construct a stepp window object besides
  using new("stwin", ...). A STEPP window specifies a pattern of subpopulation you would like
  to explore. There are three kinds of patterns: "sliding window", "event-based sliding window"
  and "tail-oriented window".\cr
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
     values for each subpopulation from the maximum value of the window. The utility function, 
     gen.tailwin, can be used to generate these two vectors based on the number of desired 
     subpopulations.\cr
  \cr
  When "sliding_events" is chosen r1 and r2 are set to \code{NULL}, while when "sliding" or
  "tail-oriented" are chosen e1 and e2 are set to \code{NULL}.\cr
  \cr
  These pattern are based on all patients.
}


\author{
 Wai-ki Yip
}


\seealso{
  \code{\linkS4class{stwin}},     \code{\linkS4class{stsubpop}},  \code{\linkS4class{stmodelKM}}, 
  \code{\linkS4class{stmodelCI}}, \code{\linkS4class{stmodelGLM}}, 
  \code{\linkS4class{steppes}},   \code{\linkS4class{stmodel}}, 
                                  \code{\link{stepp.subpop}},     \code{\link{stepp.KM}},
  \code{\link{stepp.CI}},         \code{\link{stepp.GLM}},
  \code{\link{stepp.test}},       \code{\link{estimate}},         \code{\link{generate}}

}

\examples{
  # create a stepp window object of type "sliding", 
  # subpopulation size is 200 and allows only 50 patients
  # between overlapping windows
  mywin <- stepp.win(type="sliding", r1=50, r2=200)

  # print a summary of the stepp window object
  summary(mywin) 

  # create a stepp window object of type "sliding_events", 
  # (event-based) subpopulation size is 200 and allows
  # only 50 events between overlapping windows
  mywin <- stepp.win(type="sliding_events", e1=50, e2=200)

  # print a summary of the stepp window object
  summary(mywin) 
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~event-based }
\keyword{ ~sliding }% __ONLY ONE__ keyword per line
