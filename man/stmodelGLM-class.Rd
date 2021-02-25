\name{stmodelGLM-class}
\Rdversion{1.1}
\docType{class}
\alias{stmodelGLM-class}
\alias{estimate,stmodelGLM-method}
\alias{print,stmodelGLM-method}
\alias{test,stmodelGLM-method}

\title{Class \code{"stmodelGLM"}}
\description{
	This is the stepp model for data arising from the following Generalized Linear Models:
	1. gaussian with identity link,
	2. binomial with logit link, and
	3. Poisson with log link.

	One can specify additional covariates in the model using \code{model.matrix} in R.

}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("stmodelGLM", ...)} or by \cr 
  the construction function stmodel.GLM.
}
\section{Slots}{
  \describe{
    \item{\code{coltrt}:}{Object of class \code{"numeric"} \cr
	the treatment variable}
    \item{\code{colY}:}{Object of class \code{"numeric"} \cr
	a vector containing the outcome }
    \item{\code{trts}:}{Object of class \code{"numeric"} \cr
	a vector containing the codes for the 2 treatment groups, first and second treatment groups, respectively }
    \item{\code{MM}:}{Object of class \code{"ANY"} \cr
	a model matrix for extra adjustment covariates; default is NULL \cr
	currently, stepp can only support covariates with numeric values, no logical values or factors allowed \cr
	these values need to be converted to numeric with some encoding schemes}
    \item{\code{glm}:}{Object of class \code{"character"} \cr
	 the glm to be used for analysis: "gaussian", "binomial", "poisson"}
    \item{\code{link}:}{Object of class \code{"character"} \cr
	the link function; reserved for future use }
  }
}

\section{Extends}{
Class \code{"\linkS4class{stmodel}"}, directly.
}
\section{Methods}{
  \describe{
    \item{estimate}{\code{signature(.Object = "stmodelGLM")}: \cr
	estimate the effect in absolute and relative scale of the overall and each subpopulation }
    \item{print}{\code{signature(.Object = "stmodelGLM")}: \cr
	print the estimate, covariance matrices and statistics }
    \item{test}{\code{signature(.Object = "stmodelGLM")}: \cr
	perform the permutation tests or GEE and obtain various statistics }
  }
}

\value{
  The new method returns the stmodelGLM object.

  The estimate method returns a list with the following fields:

  \item{model}{the stepp model - "GLMGe" - Gaussian model, "GLMBe" - binomial model, and "GLMPe" - Poisson model}
  \item{sObs1}{a vector of effect estimates of all subpopulations based on the first treatment}
  \item{sSE1}{a vector of standard errors of effect estimates of all subpopulations based on the first treatment}
  \item{oObs1}{effect estimate of the entire population based on the first treatment}
  \item{oSE1}{the standard error of the effect estimate of the entire population based on the first treatment}
  \item{sObs2}{a vector of effect estimates of all subpopulations based on the first treatment}
  \item{sSE2}{a vector of standard errors of effect estimates of all subpopulations based on the first treatment}
  \item{oObs2}{effect estimate of the entire population based on the first treatment}
  \item{oSE2}{the standard error of the effect estimate of the entire population based on the first treatment}
  \item{sglmw}{Wald's statistics for the effect estimate differences between the two treatments}
  \item{RD}{a vector of effect estimates difference of all subpopulations between the two treatments}
  \item{RDSE}{a vector of the standard error effect estimates difference of all subpopulations between the two treatments}
  \item{ORD}{overall difference of effect estimates of the entire population between the two treatments}
  \item{ORDSE}{the standard error of the overall difference of the effect estimates of the entire population between the two treatments}
  \item{logR}{a vector of log ratio of effect estimates of all subpopulations between the two treatments}
  \item{logRSE}{a vector of standard error of log ratio of effect estimates of all subpopulations between the two treatments}
  \item{ologR}{log ratio of effect estimates of the entire population between the two treatments}
  \item{ologRSE}{the standard error of the log ratio of effect estimates of the entire population between the two treatments}
  \item{sglmlogrw}{Wald's statistics for the log ratio of effect estimates between the two treatments}

  The test method returns a list with the following fields:
  \item{model}{the stepp model - "GLMt"}
  \item{sigma}{the covariance matrix for subpopulations based on effect differences}
  \item{hasigma}{the homogeneous association covariance matrix for subpopulations based on effect differences}
  \item{HRsigma}{the covariance matrix for the subpopulations based on hazard ratio}
  \item{haHRsigma}{the homogeneous association covariance matrix for subpopulations based on hazard ratio}
  \item{pvalue}{the supremum pvalue based on effect difference}
  \item{chi2pvalue}{the chisquare pvalue based on effect difference}
  \item{hapvalue}{the homogeneous association pvalue based on effect difference}
  \item{HRpvalue}{the supremum pvalue based on hazard ratio}
  \item{haHRpvalue}{the homogeneous association pvalue based on hazard ratio}  

}


\author{
  Wai-Ki Yip
}

\seealso{
  \code{\linkS4class{stwin}},     \code{\linkS4class{stsubpop}},  \code{\linkS4class{stmodelKM}}, 
  \code{\linkS4class{stmodelCI}}, 
  \code{\linkS4class{steppes}},   \code{\linkS4class{stmodel}}, 
  \code{\link{stepp.win}},        \code{\link{stepp.subpop}},     \code{\link{stepp.KM}},
  \code{\link{stepp.CI}},         \code{\link{stepp.GLM}},
  \code{\link{stepp.test}},       \code{\link{estimate}},         \code{\link{generate}}

}
\examples{
showClass("stmodelGLM")
}
\keyword{classes}
