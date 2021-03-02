\name{bigCI}
\alias{bigCI}
\docType{data}
\title{
  The BIG 1-98 trial dataset for cumulative incidence STEPP.
}

\description{
  This data set contains 2,685 patients in the Breast International Group (BIG) 1-98 randomized clinical trial.  The BIG 1-98 is a Phase III
  clinical trial of 8,010 post menopausal women with hormone-receptor-positive early invasive breast cancer who were randomly assigned adjuvant
  therapy of letrozole or tamoxifen.  Patterns of treatment effects for varying levels of the biomarker Ki-67 labeling index, a measure of 
  cell proliferation, were analyzed using STEPP.  The STEPP analysis showed that letrozole was more effective than tamoxifen for patients with 
  tumors expressing the highest levels of the Ki-67 labeling index.  The two treatment arms are letrozole and tamoxifen. 

}

\usage{data(bigCI)}
\format{
  There are four columns of numeric values: 
 	trt (treatment group),
	time (time to event),
	event (competing event types), and
	ki67 (continuous measurement of biomarker Ki-67).
}

\source{
  The Breast International Group (BIG) 1-98 Steering Committee and the International Breast Cancer Study Group (IBCSG) are acknowledged for 
  permission to use the data from the BIG 1-98 trial.
}

\references{
  Lazar A, Cole B, Bonetti M, Gelber R (2010), "Evaluation of treatment-effect heterogeneity using biomarkers measured on a continuous scale: subpopulation treatment effect pattern plot."
  J Clin Oncol, 28(29), 4539-44.

  Viale G et al (2008), "Prognostic and predictive value of cnetrally reviewed Ki-67 labeling index in postmenopausal women with endocrine-responsive breast cancer:
  results from Breast International Group Trial 1-98 comparing adjuvant tamoxifen and letrozole." J Clin Oncol, 28(34), 5569-75.
}


\examples{
data(bigCI)

rxgroup <- bigCI$trt
time    <- bigCI$time
evt     <- bigCI$event
cov     <- bigCI$ki67

# analyze using Cumulative Incidence method with
# sliding window size of 150 patients and a maximum of 50 patients in common
#
swin    <- new("stwin", type="sliding", r1=50, r2=150) # create a sliding window
subp    <- new("stsubpop")                             # create subpopulation object
subp    <- generate(subp, win=swin, covariate=cov) # generate the subpopulations
summary(subp)					   # summary of the subpopulations

# create a stepp model using Cumulative Incidences to analyze the data
#
smodel  <- new("stmodelCI", coltrt=rxgroup, trts=c(1,2), coltime=time, coltype=evt, timePoint=4)

statCI  <- new("steppes")		  # create a test object based on subpopulation and window
statCI  <- estimate(statCI, subp, smodel) # estimate the subpopulation results
# Warning: In this example, the permutations have been set to 0 to allow the function
# to finish in a short amount of time.  IT IS RECOMMEND TO USE AT LEAST 2500 PERMUTATIONS TO 
# PROVIDE STABLE RESULTS.
statCI  <- test(statCI, nperm=0)       # permutation test with 0 iterations

print(statCI)				  # print the estimates and test statistics
plot(statCI, ncex=0.65, legendy=30, pline=-15.5, color=c("blue","gold"),
        pointwise=FALSE, 
        xlabel="Median Ki-67 LI in Subpopulation (\% immunoreactivity)",
        ylabel="4-year Cumulative Incidence", 
        tlegend=c("Tamoxifen", "Letrozole"), nlas=3)

}

\keyword{datasets}

