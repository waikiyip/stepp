\name{gen.tailwin}
\alias{gen.tailwin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Utility function to generate tail-oriented window 
}

\description{
   Utility function to generate tail-oriented windows inputs given the approximate number of subpopulations
   desired. 
}
\usage{
  gen.tailwin(covariate, nsub, dir="LE")
}
\arguments{
  \item{covariate}{covariate values}
  \item{nsub}{number of tail-oriented subpopulations to be generated}
  \item{dir}{"LE" (default) or "GE" - subpopulations with covariate values less than or equal/greater than or equal 
		to the generated values}
}
\value{
  It returns a list with fields: $v - vector of covariate values to be used in the constructor stepp.win. 
  $np - vector of subpopulation size associate with each tail-oriented window defined by $v.
}
\details{
  Use this together with the constructor, stepp.win, to generate tail-oriented windows.
}

\author{
 Wai-ki Yip
}


\seealso{
  \code{\linkS4class{stwin}},  \code{\link{stepp.win}} 
}
