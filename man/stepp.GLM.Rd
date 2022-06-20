\name{stepp.GLM}
\alias{stepp.GLM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	The constructor to create the stmodelGLM object
}
\description{
  	This is the constructor function for stmodelGLM object.  Three kinds of GLMs are supported at this time:
  	gaussian with identity link, binomial with logit link and Poisson with log link.
}
\usage{
	stepp.GLM(coltrt, trts, colY, MM, glm, link)
}

\arguments{
  \item{coltrt}{the treatment variable}
  \item{trts}{a vector containing the codes for the 2 treatment arms, 1st and 2nd treatment groups, respectively}
  \item{colY}{a vector containing the outcome}
  \item{MM}{a model matrix for additional covariates; default is NULL\cr
		no logical values or factors allowed\cr
		these values need to be converted to numeric with some encoding scheme}
  \item{glm}{the glm to be used for analysis: "gaussian", "binomial", "poisson"}
  \item{link}{the link function; reserved for future use}
}

\value{
  It returns a stmodelGLM object.
}

\author{
  Wai-Ki Yip
}

\references{
  Yip WK, Bonetti M, Cole BF, Barcella W, Wang XV, Lazar A and Gelber R (2016), 
  "STEPP- Subpopulation Analysis for Continuous, Binary and Count Outcomes",
  Clinical Trials 2016 August ; 13(4): 382-290, doi:10.1177/1740774516643297.
}

\seealso{
  \code{\linkS4class{stwin}},     \code{\linkS4class{stsubpop}},  \code{\linkS4class{stmodelKM}}, 
  \code{\linkS4class{stmodelCI}}, \code{\linkS4class{stmodelGLM}}, 
  \code{\linkS4class{steppes}},   \code{\linkS4class{stmodel}}, 
  \code{\link{stepp.win}},        \code{\link{stepp.subpop}},     \code{\link{stepp.KM}},
  \code{\link{stepp.CI}},         
  \code{\link{stepp.test}},       \code{\link{estimate}},         \code{\link{generate}}

}
\examples{
#
# see the example for the aspirin data for an example for a binomial glm model
#
}