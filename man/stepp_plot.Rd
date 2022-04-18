\name{stepp_plot}
\alias{stepp_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  A function to generate the stepp plots
}

\description{
  This function will be deprecated in the future.  Please use S4 classes and its generic
  function e.g. summary, print and plot for future development.

  A method to generate the 3 stepp plots for users to examine the treatment-covariate interactions
  in various models data arising from two treatment arms of a clinical trial. Three plots are
  generated:\cr\cr
  1. effect estimates of both treatments against the median of each subpopulation\cr
  2. effect differences (in absolute scale) between the two treatments against the
     median of each subpopulation.\cr
  3. effect ratios (in relative scale) between the two treatments against the median
     of each subpopulation.

}

\usage{
  stepp_plot(x, legendy = 30, pline = -2.5, color = c("red", "black"), 
    ylabel= "Specify Timepoint & Endpoint", xlabel="Subpopulations by Median Covariate", 
    ncex = 0.7, tlegend=c("Specify 1st Treatment", "Specify 2nd Treatment"), 
    nlas = 0, alpha = 0.05, pointwise = FALSE, ci = TRUE, pv = TRUE, 
    showss = TRUE, ylimit=c(0,100,-100,100,0,3), which="", 
    noyscale=FALSE, at=NA, subplot=FALSE) 
}

\arguments{
  \item{x}{the steppes object returned from stepp, analyze.KM.stepp or analyze.CumInc.stepp}
  \item{legendy}{optional - the vertical location of the legend according to the units on the y-axis; default is 30}
  \item{pline}{optional - specify the vertical location of the p-value, that is, on which MARgin line, starting at 0 counting outwards; default is -2.5}
  \item{color}{optional - a vector containing the line colors for the 1st and 2nd treatment, respectively; default is "red" and "black"}
  \item{ylabel}{optional - specify the label for the y-axis; default is "Specify Timepoint & Endpoint"}
  \item{xlabel}{optional - specify the label for the x-axis; default is "Subpopulations by Median Covariate"}
  \item{ncex}{optional - specify the size of the text for the sample size annotation, that is, the character expansion factor; default is 0.7}
  \item{tlegend}{optional - a vector containing the treatment labels, 1st and 2nd treatment, respectively; default is c("Specify 1st Treatment", "Specify 2nd Treatment"}
  \item{nlas}{optional - specify the las paramter (0,1,2,3) to determine the orientation of the sample size annotation; default is 0}
  \item{alpha}{optional - specify the significance level; default is 0.05}
  \item{pointwise}{optional - specify pointwise confidence intervals (pointwise=TRUE), or confidence bands; default is FALSE}
  \item{ci}{optional - specify if you want to display the conf. interval or band; default is TRUE}
  \item{pv}{optional - pvalue will be displayed in the plots; FALSE: pvalue will not be displayed in the plots; default is TRUE}
  \item{showss}{optional - specify the label for the x-axisoptional - show the sample size; FALSE: sample size will not be shown in the plots; default is TRUE}
  \item{ylimit}{optional - pvalue will be displayed in the plots; FALSE: pvalue will not be displayed in the plots; default is TRUE}
  \item{which}{optional - if a subset of the plots is required, specify a subset of the numbers 1:3, which is the default}
  \item{noyscale}{optional - do not scale to the y-axis; default is FALSE}
  \item{at}{optional - specify the horizontal position of the p-value; use together with pline; default is at the minimal value of x}
  \item{subplot}{optional - generate a subpopulation distribution subplot on top of the first STEPP plot; default is FALSE (no subplot)}
}

\details{
  This is the function to generate three stepp plots for user to explore treatment-covariate 
  interactions in various models data arising from two treatment arms of a clinical trial. You can choose
  to produce a single figure with all the three plots or only a subset of these through the use of the `which' argument.

  Through these plots, the user can examine the differences in effects (either in absolute or relative scale)
  with respect to each subpopulation.

  This is mainly a graphics function.  The STEPP analysis is done.  Most of the parameters are for
  graphics control so that one can customize the plots for presentation purposes.
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

  Old functions to be deprecated: stepp, stepp_summary, stepp_print, analyze.KM.stepp and analyze.CumInc.stepp.
}

\examples{
  # see example in the documentation for the function stepp.
}
