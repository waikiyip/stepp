setClassUnion("null_or_numeric", c("numeric", "NULL"))

#################################################################
#
# stwin.R
#
################
# stepp window #
################
setClass("stwin",
  representation(
    type    = "character",	        # stepp window type
  	r1      = "null_or_numeric",    # sliding window: largest number of patients in common
  					                        #                 among consecutive subpopulations
  					                        # tail-oriented window: < vector
  	r2      = "null_or_numeric",    # sliding window: minimum number of patients of each
                                    #                 subpopulation
    e1      = "null_or_numeric",    # event-based window: largest number of events in common
                                    #                     among consecutive subpopulations
    e2      = "null_or_numeric"),    # event-based window: minimum number of events of each
                                    #                     subpopulation
  prototype(type = "sliding", r1 = 5, r2 = 20, e1 = NULL, e2 = NULL)
)

setMethod("initialize", "stwin",
  function(.Object, type, r1, r2, e1, e2, ...) {
    if (missing(type)) {
      type <- "sliding"
    }
    if (missing(r1) & missing(e1)) {
      r1 <- 5
    } else if (missing(r1) & !missing(e1)) {
      r1 <- NULL
    }
    if (missing(r2) & missing(e2)) {
      r2 <- 20
    } else if (missing(r2) & !missing(e2)) {
      r2 <- NULL
    }
    if (missing(e1)) {
      e1 <- NULL
    }
    if (missing(e2)) {
      e2 <- NULL
    }
    if (type == "tail-oriented"){
      if (length(r1) != length(r2)) {
        errmsg2 <- "minpatspop (r1) and patspop (r2) must have the same length for 'tail-oriented' windows."
        stop(errmsg2)
      }
      r1 <- unique(sort(c(r1)))
      r2 <- unique(sort(c(r2)))
    }
    .Object@type <- type
    .Object@r1 <- r1
    .Object@r2 <- r2
    .Object@e1 <- e1
    .Object@e2 <- e2
    if (!validObject(.Object)) stop("")
    callNextMethod(.Object, ...)
  }
)

setValidity("stwin", 
	function(object) {
    status <- TRUE
	  if (is.na(match(object@type, c("sliding", "sliding_events", "tail-oriented")))) {
  		status <- FALSE
      errmsg1 <- "invalid stepp window type:"
  		errmsg1 <- paste(errmsg1, object@type)
  		print(errmsg1)
	  } else {
  	  if (object@type == "tail-oriented") {
        if (length(unique(object@r1)) > 1 & length(unique(object@r2)) > 1) {
          status <- FALSE
          errmsg2 <- "either minpatspop (r1) or patspop (r2) must contain a single unique value for 'tail-oriented' windows (type '?stepp.win' for more details)."
          print(errmsg2)
        }
        if (any(object@r1 <= object@r2)) {
          status <- FALSE
          errmsg2 <- "minpatspop (r1) must be larger (element-wise) than patspop (r2) for 'tail-oriented' windows."
          print(errmsg2)
        }
      } else if (object@type == "sliding") {
        if (is.na(object@r1) | is.na(object@r2)) {
          status <- FALSE
          errmsg2 <- "minpatspop (r1) and patspop (r2) must be both integers greater than or equal to 1."
          print(errmsg2)
        }
        if (object@r1 >= object@r2) {
          status <- FALSE
    		  errmsg2 <- "minpatspop (r1) MUST be less than patspop (r2)"
          print(errmsg2)
        }
      } else if (object@type == "sliding_events") {
        if (is.na(object@e1) | is.na(object@e2)) {
          status <- FALSE
          errmsg2 <- "mineventspop (e1) and eventspop (e2) must be both integers greater than or equal to 1."
          print(errmsg2)
        }
        if (object@e1 >= object@e2) {
          status <- FALSE
          errmsg2 <- "mineventspop (e1) MUST be less than eventspop (e2)"
          print(errmsg2)
        }
      }
    }
	  return(status)
	}
)

setMethod("summary", signature = "stwin",
	definition=function(object) {
	  write(paste("Window type:", object@type), file = "")
	  if (object@type == "tail-oriented") {
	    temp1 <- paste(paste(object@r1), collapse=" ")
	    temp2 <- paste(paste(object@r2), collapse=" ")
	    write(paste("Subpopulation for patients less than or equal to: ",    temp1), file = "")
	    write(paste("Subpopulation for patients greater than or equal to: ", temp2), file = "")
	  } else if (object@type == "sliding") {
	    write(paste("Number of patients per subpopulation (patspop r2):", object@r2), file = "")
    	    write(paste("Largest number of patients in common among consecutive subpopulations (minpatspop r1):", object@r1), file = "")
    } else if (object@type == "sliding_events") {
      write(paste("Number of events per subpopulation (eventspop e2):", object@e2), file = "")
          write(paste("Largest number of events in common among consecutive subpopulations (mineventspop e1):", object@e1), file = "")
	  }
	}
)

# constructor function for stepp window
stepp.win <- function(type = "sliding", r1 = 5, r2 = 20, e1 = NULL, e2 = NULL) {
	if (type == "tail-oriented"){
	  r1 <- unique(sort(c(r1)))
	  r2 <- unique(sort(c(r2)))
	} else if (type == "sliding") {
    e1 <- NULL
    e2 <- NULL
  } else if (type == "sliding_events") {
    r1 <- NULL
    r2 <- NULL
  }
	sw <- new("stwin", type = type, r1 = r1, r2 = r2, e1 = e1, e2 = e2)
  validObject(sw)
	return(sw)
}

#
# Utility function to generate tail subgroups
# - generate a vector of covariate values
#	for tail-oriented window
#	with each subgroup roughly about the specified percentage 
#	of the total cohort
#
# 1. covariate - continuous covariate of interest
# 2. nsub - number of subpopulation you would like (approx) in addition to the entire cohort
# 3. dir - "GE" or "LE"
#
gen.tailwin <- function(covariate, nsub, dir="LE") {
  if (nsub > length(unique(covariate))) stop("No of subpopulations can't be larger than covariate unique values")
  if (nsub <= 1) stop("No of subpopulations must be larger than 1")
  perc <- 1/(nsub + 1)

  cov.t <- sort(covariate)
  remain <- length(cov.t)
  lcov <- remain
  grpsize <- round(remain*perc,0)
  v <- rep(0, remain)
  np <- rep(0, remain)

  i <- 0
  if (dir == "LE"){
    while (remain >= 2*grpsize) {
    	cov.val <- cov.t[remain - grpsize]
    	# handle the case where the covariate value is the max.
    	if (cov.val == max(cov.t)) {
    	  pro.val <- which(cov.t < max(cov.t)) 
    	  cov.val <- cov.t[pro.val[length(pro.val)]]
    	}
    	mat <- which(cov.val == cov.t)
    	cov.r <- mat[length(mat)]
    	v[lcov - i] <- cov.val
    	cov.t <- cov.t[1:cov.r]
    	remain <- length(cov.t)
    	np[lcov - i] <- remain
    	i <- i + 1
    } 
    v <- v[(lcov - i + 1):lcov]
    np <- np[(lcov - i + 1):lcov]
  }
  else {
    while (remain >= 2*grpsize) {
    	cov.val <- cov.t[grpsize + 1]
    	# handle the case where the covariate value is the min.
    	if (cov.val == min(cov.t)) {
    	  pro.val <- which(cov.t > min(cov.t))
    	  cov.val <- cov.t[pro.val[1]]
    	}

    	temp <- which(cov.val == cov.t)
    	cov.l <- temp[1]
    	v[i + 1] <- cov.val
    	cov.t <- cov.t[(cov.l):length(cov.t)]
    	remain <- length(cov.t)
    	np[i + 1] <- length(cov.t)
    	i <- i + 1
    }
    v <- v[1:i]
    np <- np[1:i]
  }

  return(list(v = v, np = np))
}
