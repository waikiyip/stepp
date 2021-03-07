\name{balance_example}
\alias{balance_example}
\docType{data}
\title{
  Sample data to use with the \code{balance_patients()} function.
}

\description{
  A sample data set containing a single column which is used to show the usage of the \code{balance_patients()}
  function.
}

\usage{data(balance_example)}
\format{
  There is only one numeric column:
 	covar (continuous covariate).
}

\references{
  Lazar A, Cole B, Bonetti M, Gelber R (2010), "Evaluation of treatment-effect heterogeneity using biomarkers measured on a continuous scale: subpopulation treatment effect pattern plot."
  J Clin Oncol, 28(29), 4539-44.
}

\seealso{
  \code{\link{balance_patients}()}
}

\examples{
\dontrun{
data(balance_example, package = "stepp")
ranger2 <- c(950, 1050)
ranger1 <- c(300, 500)
maxnsubpops <- 50


\dontrun{res_bal <- balance_patients(ranger1, ranger2, maxnsubpops, balance_example$covar,
  plot = TRUE, verbose = TRUE, contour = TRUE, nlevels = 6, border = TRUE)}
  

}

\keyword{datasets}
