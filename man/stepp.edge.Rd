\name{stepp.edge}
\alias{stepp.edge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  The method performs an edge analysis on the STEPP GLM model estimate objects.
}

\description{
   Perform an edge analysis on the STEPP GLM model estimate objects.
}
\usage{
  stepp.edge(est, criteria, j=2, boot=0, seed=17, showstatus=TRUE, debug=0)
}
\arguments{
  \item{est}{a STEPP estimate object.}
  \item{criteria}{criteria to be used to identify the cut point; abs or rel scale and by how much.}
  \item{j}{number of treatment, default to 2}
  \item{boot}{perform a bootstrap analysis, default is none (0).}
  \item{seed}{seed used for bootstrap, default is 17.}
  \item{showstatus}{show the status of bootstrap, default is TRUE.}
  \item{debug}{internal debug flag, default is 0}
}


\value{
  It returns the result of the identified edge with the STEPP subpopulations.
}
\details{

  The criteria argument is a list with three elements:

  trtid - the treatment id;

  scale - "A" for absolute scale or "R" for relative scale;

  threshold - amount in either absolute or relative scale that would consider a jump.

  e.g. crit <- list(trtid=1, scale="A", threshold=-0.03)

  Only support for STEPP GLM models for now. Bootstrap is computational time intensive  but it will
  provide you with a quantification of the variability of the edge identified.
}

\author{
 Wai-ki Yip
}


\seealso{
  \code{\linkS4class{stwin}} 

}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
