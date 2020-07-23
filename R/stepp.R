#################################################################
#
# stepp.R
#
#############################
# old stepp interface       #
#############################
# These calls are maintained by backward compatibility reasons:
#
stepp <- function(trttype, coltrt, coltime, colcens=0, coltype=0, colvar, trts, 
			 patspop, minpatspop, timest, noperm)
{
  win.stepp    <- new("stwin", type="sliding", r1=minpatspop, r2=patspop)
  subp.stepp   <- new("stsubpop")

  if (trttype == "KM") {
    indata  <- cbind(coltrt, coltime, colcens, colvar)
    indata  <- indata[apply(indata, 1, function(x) !any(is.na(x))), 
        , drop = FALSE]
    coltrt  <- indata[, 1]
    coltime <- indata[, 2]
    colcens <- indata[, 3]
    colvar  <- indata[, 4]

    mod.stepp  <- new("stmodelKM", coltrt=coltrt, survTime=coltime, censor=colcens, 
				trts=trts, timePoint=timest)
  } else if (trttype == "CI") {
    indata  <- cbind(coltrt, coltime, coltype, colvar)
    indata  <- indata[apply(indata, 1, function(x) !any(is.na(x))), 
        , drop = FALSE]
    coltrt  <- indata[, 1]
    coltime <- indata[, 2]
    coltype <- indata[, 3]
    colvar  <- indata[, 4]

    mod.stepp  <- new("stmodelCI", coltrt=coltrt, coltime=coltime, coltype=coltype, 
				trts=trts, timePoint=timest)
  } else {
    stop("Treatment type can only be KM or CI !")
  }

  subp.stepp   <- generate(subp.stepp, win=win.stepp, covariate=colvar)

  result.stepp <- new("steppes")
  result.stepp <- estimate(result.stepp, subp.stepp, mod.stepp)
  result.stepp <- test(result.stepp, noperm)

  return(result.stepp)
}

stepp_summary <- function(x) {
  summary(x)
}

stepp_print <- function(x, estimate=TRUE, cov=TRUE, test=TRUE) {
  print(x@model, x, estimate, cov, test)
}

stepp_plot <- function(x, legendy = 30, pline = -2.5, color = c("red", "black"),
	ylabel= "Specify Timepoint & Endpoint", xlabel="Subpopulations by Median Covariate",
	ncex = 0.7, tlegend=c("Specify 1st Treatment", "Specify 2nd Treatment"), 
	nlas = 0, alpha = 0.05, pointwise = FALSE, diff = TRUE, ci = TRUE, pv = TRUE, 
	showss = TRUE, ylimit=c(0,100,-100,100,0,3), dev="", together=FALSE, noyscale=FALSE, at=NA, subplot=FALSE) {
  plot(x, legendy=legendy, pline=pline, color=color,
		ylabel=ylabel, xlabel=xlabel, ncex=ncex, tlegend=tlegend, nlas=nlas, alpha=alpha, 
		pointwise=pointwise, diff=diff, ci=ci, pv=pv, showss=showss, ylimit=ylimit, 
		dev=dev, together=together, noyscale=noyscale, at=at, subplot=subplot)
}

analyze.KM.stepp <- function ( coltrt, coltime, colcens, colvar, trts, patspop, minpatspop, 
				timest, noperm=2500,
			 	ncex = 0.70, legendy = 30, pline = -2.5, color = c("red", "black"),
			 	xlabel="Subpopulations by Median Covariate",
			 	ylabel = "?-year Disease-Free Survival", 
			 	tlegend = c("1st Treatment", "2nd Treatment"),
			 	nlas = 3, pointwise=FALSE) {
  stepp.KM <- stepp("KM", coltrt=coltrt, coltime=coltime, colcens=colcens, colvar=colvar,
			  trts=trts, patspop=patspop, minpatspop=minpatspop, timest=timest, noperm=noperm)
  stepp_summary(stepp.KM)
  stepp_print(stepp.KM)
  stepp_plot(stepp.KM, ncex=ncex,legendy=legendy, pline=pline, color=color,
		xlabel=xlabel, ylabel=ylabel, tlegend=tlegend, nlas=nlas,
		pointwise=pointwise)
  return(stepp.KM)
}

analyze.CumInc.stepp <- function(coltrt, coltime, coltype, colvar, trts, patspop, minpatspop, 
    				timest, noperm=2500,
				ncex = 0.7, legendy = 30, pline = -2.5, color = c("red", "black"),
    				xlabel = "Subpopulations by Median Covariate",
				ylabel = "?-year Disease-Free Survival", 
    				tlegend = c("1st Treatment", "2nd Treatment"), 
   				nlas = 3, pointwise = FALSE) {
  stepp.CI <- stepp("CI", coltrt=coltrt, coltime=coltime, coltype=coltype, colvar=colvar,
			  trts=trts, patspop=patspop, minpatspop=minpatspop, timest=timest, noperm=noperm)
  stepp_summary(stepp.CI)
  stepp_print(stepp.CI)
  stepp_plot(stepp.CI, ncex=ncex, legendy=legendy, pline=pline, color=color,
     		xlabel=xlabel, ylabel=ylabel, tlegend=tlegend, nlas=nlas,
		pointwise=pointwise)
  return(stepp.CI)
}

####
# Release notes for STEPP
#
stepp.rnote <- function() {
  cat("Release Note for STEPP version 3.2-0 (May 1, 2018)")
  cat("\nThis version of the STEPP package contains major additions.  Some of the features have some limitations.")

  cat("\n1. Multiple comparisons up to 8 groups (arms) - The user can specify up to 8 treatment groups in one STEPP")
  cat("\n   analysis. Treatment 0 is the baseline and all other treatments will be compared against the ")
  cat("\n   baseline treatment. The default is two treatment groups as before.")

  cat("\n2. Event-based sliding window - In addition to the sliding window pattern, the user can now specify")
  cat("\n   event-based sliding windows for the analysis. The event-based windows work for cumulative incidence")
  cat("\n   and Kaplan-Meier models. The first STEPP plot is modified to allow the user to specify an option to")
  cat("\n   show the subpopulations above the STEPP plot.")

  cat("\n3. Tail-oriented window - In addition to the sliding window pattern, the user can now specify")
  cat("\n   tail-oriented windows for the analysis. The tail-oriented windows work for all models. The first STEPP ")
  cat("\n   plot is modified to allow the user to specify an option to show the subpopulations above the STEPP plot.")

  cat("\n4. Updated BIG data and aspirin data - The original de-identified data sets are provided with this")
  cat("\n   version. The  bigKM data set is for the STEPP analysis based on the Kaplan-Meier method, and the bigCI ")
  cat("\n   data set is for the STEPP analysis based on the Cumulative Incidence method. The aspirin data set")
  cat("\n   is for the STEPP analysis based on the Bernoulli response using GLM.") 

  cat("\n5. Cut-Point (or edge identification) - This version contains an experimental feature to identify")
  cat("\n   the cut-point (or edge) for critical STEPP subgroups. A bootstrap confidence interval is ")
  cat("\n   provided for the cut-point identified. This feature is implemented only for the GLM models ")
  cat("\n   with 2 treatment-arms.")

  cat("\n")
  cat("\nAdditional STEPP related R software packages are available:")
  cat("\nSoftware for Meta-STEPP is available on  http://bcb.dfci.harvard.edu/~vwang/MetaSTEPP.html.")

  cat("\n")
}
