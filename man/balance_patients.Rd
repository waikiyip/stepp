\name{balance_patients}
\alias{balance_patients}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Utility function for determining the optimal values for generating the subpopulations.
}

\description{
  Utility function for determining the optimal values of the number of subpopulations and the corresponding r1 and r2
  values for creating subpopulations with the sliding window approach. The optimal values are those that make the
  subpopulations more balanced by minimizing the variance of the subpopulation sizes.
}

\usage{
  balance_patients(range_r1, range_r2, maxnsubpops, covar, verbose = FALSE,
    plot = FALSE, contour = FALSE, nlevels = 5)
}

\arguments{
  \item{range_r1}{numeric vector with two elements providing the range of values for the r1 parameter}
  \item{range_r2}{numeric vector with two elements providing the range of values for the r2 parameter}
  \item{maxnsubpops}{length-one numeric vector providing the maximum number of subpopulations to consider}
  \item{covar}{numeric vector containing the covariate values to use for generating the subpopulations}
  \item{verbose}{length-one logical vector; if \code{TRUE} prints a summary of the results in the console}
  \item{plot}{length-one logical vector; if \code{TRUE} produces a diagram showing the results of the calculations}
  \item{contour}{length-one logical vector; if \code{TRUE} adds to the plot the variance contour lines for each
  subpopulation number}
  \item{nlevels}{length-one numeric vector providing the number of contour lines to plot}
  \item{showstatus}{length-one logical vector; if \code{TRUE} displays a bar showing the progress of the calculations; default is TRUE}
}

\value{
  The \code{balance_patients()} function returns a list with the following items:
  \item{r1_best}{length-one numeric vector with overall best value of the r1 parameter}
  \item{r2_best}{length-one numeric vector with overall best value of the r2 parameter}
  \item{var_best}{length-one numeric vector with overall minimum value of the sizes variance}
  \item{nsubpops_best}{length-one numeric vector with overall best value for the number of subpopulations}
  \item{all_res}{numeric matrix with the details of all the calculations}
}

\references{
  Bonetti M, Gelber RD. Patterns of treatment effects in subsets of patients in clinical trials. Biostatistics 2004; 5(3):465-481.

  Bonetti M, Zahrieh D, Cole BF, Gelber RD. A small sample study of the STEPP approach to assessing treatment-covariate interactions in survival data. Statistics in Medicine 2009; 28(8):1255-68. 

  Lazar AA, Cole BF, Bonetti M, Gelber RD.  Evaluation of treatment-effect heterogeneity using biomarkers measured on a
  continuous scale: subpopulation treatment effect pattern plot.  Journal of Clinical Oncology, 2010; 28(29): 4539-4544.
}

\author{
 Marco Bonetti, Sergio Venturini
}

\seealso{
  \code{\linkS4class{stwin}},     \code{\linkS4class{stsubpop}},
  \code{\link{stepp.win}},        \code{\link{stepp.subpop}},     \code{\link{stepp.KM}}
}

\examples{
\dontrun{
data(balance_example, package = "stepp")
ranger2 <- c(950, 1050)
ranger1 <- c(300, 500)
maxnsubpops <- 50

res_bal <- balance_patients(ranger1, ranger2, maxnsubpops, balance_example$covar,
  plot = TRUE, verbose = TRUE, contour = TRUE, nlevels = 6)
}
}
