\name{stepp_summary}
\alias{stepp_summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  The function to produce a summary of the size and various attributes of each subpopulation
}

\description{
  This function will be deprecated in the future.  Please use S4 classes and 
  corresponding generic functions e.g. summary, print and plot for future development.

  A method to print the summary of the size and various attributes for
  each subpopulation used in the STEPP analysis.
}

\usage{
  stepp_summary(x)
}

\arguments{
  \item{x}{a steppes object returned from stepp, analyze.KM.stepp or analyze.CumInc.stepp function}
}


\author{
 Wai-ki Yip, David Zahrieh, Marco Bonetti, Bernard Cole, Ann Lazar, Richard Gelber
}


\section{Warning }{
  This function together with other old functions will be depreciated in the future.  A new set of S4 classes are
  implemented to replace old interfaces.  Please use them for future development.
}


\seealso{
  The S4 classes: stwin, stsubpop, stmodelKM, stmodelCI, stmodelCOX, stmodelGLM, and steppes.

  Old functions to be deprecated: stepp, stepp_print, stepp_plot, analyze.KM.stepp, and analyze.CumInc.stepp 

}

\examples{
# set example in the documentation for the function stepp.

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
