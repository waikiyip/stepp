\name{stepp_print}
\alias{stepp_print}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  The function to print the estimate, covariance matrices and test statistics.
}

\description{
  This function will be deprecated in the future.  Please use S4 classes and generic
  functions summary, print and plot for future development.

  A method to print three essential information resulting from the STEPP analysis.
}

\usage{
  stepp_print(x, estimate=TRUE, cov=TRUE, test=TRUE)
}

\arguments{
  \item{x}{a steppes object returned from stepp, analyze.KM.stepp or analyze.CumInc.stepp function}
  \item{estimate}{whether to print the effect estimates; default to yes}
  \item{cov}{whether to print the covariance matrices; default to yes}
  \item{test}{whether to print the test statistics; default to yes}
}

\details{
  The STEPP analysis produces three important pieces of information.  User can decide to print them
  all or individually.  The three pieces of information are:\cr\cr
  1. effect estimates of each subpopulation for the two treatments, the differences in absolute 
     scale of the effects and the ratio in relative scale of the effects.\cr
  2. covariance matrics of the differences, logratios, differences in homogeneous association and 
     logratios in homogeneous association.\cr
  3. various permutation p values based on test statistics: supremum pvalues, homogeneous association pvalues, 
     chisquare pvalue. Not all statistics are available for all models.    

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

  Old functions to be deprecated: stepp, stepp_summary, stepp_plot, analyze.CumInc.stepp and analyze.KM.stepp. 

}

\examples{
# see example in the documentation for the function stepp.
}
