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
\examples{
data(bigKM)

rxgroup <- bigKM$trt
time    <- bigKM$time
evt     <- bigKM$event
cov     <- bigKM$ki67

# analyze using Kaplan-Meier method with tail-oriented window
#
nsubpop_tmp <- 10
win_tmp <- gen.tailwin(cov, nsub = nsubpop_tmp, dir = "LE")
nsubpop <- length(win_tmp$v)
swin <- new("stwin", type = "tail-oriented", r1 = win_tmp$v, r2 = rep(min(cov), nsubpop)) # create a tail-oriented window
subp <- new("stsubpop")                             # create subpopulation object
subp <- generate(subp, win = swin, covariate = cov) # generate the subpopulations
summary(subp)                             # summary of the subpopulations

# create a stepp model using Kaplan Meier Method to analyze the data
#
smodel  <- new("stmodelKM", coltrt=rxgroup, trts=c(1,2), survTime=time, censor=evt, timePoint=4)

statKM  <- new("steppes")     # create a test object based on subpopulation and window
statKM  <- estimate(statKM, subp, smodel) # estimate the subpopulation results
# Warning: IT IS RECOMMEND TO USE AT LEAST 2500 PERMUTATIONS TO  PROVIDE STABLE RESULTS.
statKM <- test(statKM, nperm = 10)       # permutation test with 10 iterations

print(statKM)         # print the estimates and test statistics
plot(statKM, ncex=0.65, legendy=30, pline=-15.5, color=c("blue","gold"),
     pointwise=FALSE, 
     xlabel="Median Ki-67 LI in Subpopulation (\% immunoreactivity)",
     ylabel="4-year Disease Free Survival", 
     tlegend=c("Letrozole", "Tamoxifen"), nlas=3)
}
