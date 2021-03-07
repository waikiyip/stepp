#################################################################
#
# stsubpop.R
#
#######################
# stepp subpopulation #
#######################

# supporting routine to generate subpopulation based just on the coviarate and window
generate.all <- function(sp, win, covariate, coltype = NULL, coltrt = NULL, trts = NULL, minsubpops = NULL) {
  r1 <- win@r1
  r2 <- win@r2
  e1 <- win@e1
  e2 <- win@e2

  if (win@type == "tail-oriented") {
    zvals   <- unique(sort(covariate))
    if (min(r1) < min(zvals) | max(r2) > max(zvals)) stop("Cannot create tail-oriented window")
    # nsubpop <- length(r1) + length(r2) + 1 ### [SV 27.02.2021: shouldn't it be length(r1)+length(r2)?]
    nsubpop <- length(r1) +length(r2)
    npatsub <- rep(0, nsubpop)
    I0      <- rep(1, nsubpop)
    I1      <- rep(length(zvals), nsubpop)
    # if (length(r1) != 0) { ### [SV 27.02.2021: shouldn't it be > 1?]
    if (length(r1) > 1) {
      for (i in 1:length(r1)) {
        npatsub[i] <- sum(covariate <= r1[i])
        sel <- which(zvals <= r1[i])
        I1[i] <- sel[length(sel)]
      }
      npatsub[length(r1) + 1] <- length(covariate)
    }
    # if (length(r2) != 0) { ### [SV 27.02.2021: shouldn't it be > 1?]
    if (length(r2) > 1) {
      npatsub[1] <- length(covariate)
      for (i in 1:length(r2)) {
        k <- i + length(r1)
        npatsub[k] <- sum(covariate >= r2[i])
        I0[k] <- which(zvals >= r2[i])[1]
      }
    }
  } else if (win@type == "sliding") {
    if (r2 >= length(covariate)) {
      errmsg <- "patspop (r2) must be strictly smaller than the number of unique covariate values."
      stop(errmsg)
    }
    zvals <- unique(sort(covariate))
    absfreq <- as.numeric(table(covariate))
    cumfreq <- cumsum(absfreq)

    I0 <- rep(0, 1000)  # max size of I0 and I1 is 1000; these limit the 
    I1 <- rep(0, 1000)  # max no. of stepp subpopulation to 1000
    I0[1] <- 1
    I1[1] <- sum(cumfreq < r2) + 1
    stopflag <- 0
    nsubpop <- 2
    while (stopflag == 0) {
      indinf <- I0[nsubpop - 1] + 1
      while ((cumfreq[I1[nsubpop - 1]] - cumfreq[indinf - 1]) > r1) {
        indinf <- indinf + 1
      }
      I0[nsubpop] <- indinf
      indsup <- I1[nsubpop - 1]
      while (((cumfreq[indsup] - cumfreq[I0[nsubpop]] + absfreq[I0[nsubpop]]) < r2) && (stopflag == 0)) {
        indsup <- indsup + 1
        stopflag <- 1 * (indsup == length(zvals))
      }
      I1[nsubpop] <- indsup
      nsubpop <- nsubpop + 1
    }

    nsubpop    <- nsubpop - 1
    npatsub    <- rep(0,  nsubpop)
    npatsub[1] <- cumfreq[I1[1]]
    for (i in 2:nsubpop) npatsub[i] <- cumfreq[I1[i]] - cumfreq[I0[i] - 1]
    I0 <- I0[1:nsubpop]
    I1 <- I1[1:nsubpop]
  } else if (win@type == "sliding_events") {
    #
    #   create 0/1 treatment assignment 
    #
    txassign <- coltrt
    for (i in 1:length(txassign)) {
      if (txassign[i] != trts[1] & txassign[i] != trts[2]) {
        txassign[i] <- NA
      }
      if (txassign[i] == trts[1]) {
        txassign[i] <- 1
      }
      else if (txassign[i] == trts[2]) {
        txassign[i] <- 0
      }
    }
    #
    # get the overlapping subgroups
    #
    # zvals is the unique values of the covariate that has been sorted at which there has been an event
    zvals <- unique(sort(covariate))
    # absfreq is the total number of events of that covariate value corresponding to zvals
    absfreq <- as.numeric(table(covariate))
    cumfreq <- cumsum(absfreq)
    #
    # Create vector of number of Type 1's recorded for each covariate
    # value for each treatment type
    # 
    type1valsTrt0 <- type1valsTrt1 <- rep(0, length(zvals))

    for (i in 1:length(covariate)) {
      if (coltype[i] == 1) {
        covariateIndices <- which(zvals == covariate[i])
        if (length(covariateIndices) != 1) {
          stop("Non-unique covariate detected in zvals vector or value missing.")
        }
        covariateIndex <- covariateIndices[1]
        #####NEW CODE
        if (txassign[i] == 0) {
          type1valsTrt0[covariateIndex] <- type1valsTrt0[covariateIndex] + 1
        }
        if (txassign[i] == 1) {
          type1valsTrt1[covariateIndex] <- type1valsTrt1[covariateIndex] + 1
        }
      }
    }

    # Create cummulative graph of typevals 
    cumtypevalsTrt0 <- cumsum(type1valsTrt0)
    cumtypevalsTrt1 <- cumsum(type1valsTrt1)
         
    I0 <- I1 <- rep(0, length(zvals))

    #
    # sets the boundaries of the first subpopulation
    #
    I0[1] <- 1

    #
    # Determines the minimum upper boundary that satisfies the value of e2 for both treatment types
    #
    upperTrt0 <- sum(cumtypevalsTrt0 < e2) + 1
    upperTrt1 <- sum(cumtypevalsTrt1 < e2) + 1
    I1[1] <- max(upperTrt0, upperTrt1)

    #
    # Check to see if the upper boundary of the first subpop is greater than
    # the size of the entire population.  This indicates that there are
    # insufficient events to satisfy E2.
    #
    if (I1[1] >= length(covariate)) {
      # obsInitialSubpopCreationFailures <- obsInitialSubpopCreationFailures + 1
      print(paste("Warning: Insufficient Type 1 Events found creating initial subpopulation (",
                  "Trt0=", cumtypevalsTrt0[length(cumtypevalsTrt0)],
                   " Trt1=", cumtypevalsTrt1[length(cumtypevalsTrt1)], ")", sep=""))
      stop("Throwing out dataset.")
    }

    #
    # Find the boundaries of each subpopulation
    #
    stopflag <- 0
    nsubpop <- 2
    while (stopflag == 0) {
      #
      # find the lower boundary of the next subpopulation
      #
      # indinf is the presumptive lower bound of the next subpopulation 
      indinf <- I0[nsubpop - 1] + 1

      while (((cumtypevalsTrt0[I1[nsubpop - 1]] - cumtypevalsTrt0[indinf - 1]) > e1) ||
        ((cumtypevalsTrt1[I1[nsubpop - 1]] - cumtypevalsTrt1[indinf - 1]) > e1)) {
        indinf <- indinf + 1
      }
      I0[nsubpop] <- indinf
      indsup <- I1[nsubpop - 1]
      #
      # find the upper boundary of the next subpopulation
      #
      while ((((cumtypevalsTrt0[indsup] - cumtypevalsTrt0[I0[nsubpop]] + type1valsTrt0[I0[nsubpop]]) < e2) ||
        ((cumtypevalsTrt1[indsup] - cumtypevalsTrt1[I0[nsubpop]] + type1valsTrt1[I0[nsubpop]]) < e2)) &&
        (stopflag == 0)) {
        indsup <- indsup + 1
        stopflag <- (indsup == length(zvals))
      }
      I1[nsubpop] <- indsup
      # increment total number of subpopulations for next loop iteration
      nsubpop <- nsubpop + 1
    }
    # decrement nsubpop to number of found subpopulations
    nsubpop <- nsubpop - 1

    #
    # If the last subpopulation is too small for either treatment arm (< e2)
    # include it in the previous subpopulation so that e2 is not violated
    #
    if (((cumtypevalsTrt0[I1[nsubpop]] - cumtypevalsTrt0[I0[nsubpop] - 1]) < e2) ||
      ((cumtypevalsTrt1[I1[nsubpop]] - cumtypevalsTrt1[I0[nsubpop] - 1]) < e2)) {
      I1[nsubpop - 1] <- I1[nsubpop]
      nsubpop <- nsubpop - 1
    }

    #
    # Make sure there are enough subpopulations
    #
    if (nsubpop < minsubpops) {
      # obsInsufficientSubpopsErrors <- obsInsufficientSubpopsErrors + 1
      print(paste("Warning: too few subpopulations (", nsubpop, ")", sep = ""))
      stop("Throwing out dataset.")
    }

    npatsub <- rep(0, nsubpop)
    npatsub[1] <- cumfreq[I1[1]]
    for (i in 2:nsubpop) npatsub[i] <- cumfreq[I1[i]] - cumfreq[I0[i] - 1]  
    neventsubTrt0 <- neventsubTrt1 <- rep(0, nsubpop)
    neventsubTrt0[1] <- cumtypevalsTrt0[I1[1]]
    neventsubTrt1[1] <- cumtypevalsTrt1[I1[1]]
    for (i in 2:nsubpop) {
      neventsubTrt0[i] <- cumtypevalsTrt0[I1[i]] - cumtypevalsTrt0[I0[i] - 1]
      neventsubTrt1[i] <- cumtypevalsTrt1[I1[i]] - cumtypevalsTrt1[I0[i] - 1]
    }
    I0 <- I0[1:nsubpop]
    I1 <- I1[1:nsubpop]
  }

  medians <- rep(0,  nsubpop)
  minz    <- rep(NA, nsubpop)
  minc    <- rep(NA, nsubpop)
  maxz    <- rep(NA, nsubpop)
  maxc    <- rep(NA, nsubpop)
  npats   <- length(covariate)
  subpop  <- matrix(rep(0, (npats*nsubpop)), ncol = nsubpop)
  for (i in 1:nsubpop) {
    subpop[, i] <- (covariate >= zvals[I0[i]]) * (covariate <= zvals[I1[i]])
    medians[i]  <- round((median(covariate[subpop[, i] == 1])), digits = 2)
    minz[i]     <- round(zvals[I0[i]], digits = 4)
    minc[i]     <- zvals[I0[i]]
    maxz[i]     <- round(zvals[I1[i]], digits = 4)
    maxc[i]     <- zvals[I1[i]]
  }

  # update the object
  sp@win <- win
  sp@colvar <- covariate
  sp@nsubpop <- nsubpop
  sp@subpop <- subpop
  sp@npatsub <- npatsub
  sp@medianz <- medians
  sp@minz <- minz
  sp@minc <- minc
  sp@maxz <- maxz
  sp@maxc <- maxc
  if (sp@win@type == "sliding_events") {
    sp@neventsubTrt0 <- neventsubTrt0
    sp@neventsubTrt1 <- neventsubTrt1
  } else {
    sp@neventsubTrt0 <- NULL
    sp@neventsubTrt1 <- NULL
  }
  sp@init <- TRUE

  return(sp)
}

setClass("stsubpop",
  representation(
    win      = "stwin",   # stepp window object
    colvar   = "numeric", # vector of covariate of interest (V) 
    nsubpop  = "numeric", # number of subpopulation generated
    subpop   = "ANY",     # matrix of subpopulations
    npatsub  = "numeric", # count of each subpopulation
    medianz  = "numeric", # median of V for each subpopulation
    minz     = "numeric", # minimum of V for each subpopulation, round to 4 digits
    maxz     = "numeric", # maximum of V for each subpopulation, round to 4 digits
    minc     = "numeric", # minimum of V for each subpopulation, actual
    maxc     = "numeric", # maximum of V for each subpopulation, actual
    neventsubTrt0 = "null_or_numeric",
    neventsubTrt1 = "null_or_numeric",
    init     = "logical"  # initialized
    )
)

setMethod("initialize", "stsubpop",
  function(.Object) {
    .Object@init <- FALSE
    return(.Object)
  }
)

setValidity("stsubpop", 
  function(object) {
    if (!is.numeric(object@colvar) || length(object@colvar) < object@win@r2) {
      print("Invalid argument.")
      return (FALSE)
    }
    return(TRUE)
  }
)

setGeneric("generate", function(.Object, win, covariate, coltype, coltrt, trts, minsubpops)
  standardGeneric("generate"))  

setMethod("generate", 
  signature="stsubpop",
  definition = function(.Object, win, covariate, coltype, coltrt, trts, minsubpops) {
    if (missing(minsubpops)) {
      minsubpops <- 2
    }
    .Object <- generate.all(.Object, win, covariate, coltype, coltrt, trts, minsubpops)
    return(.Object)
  }
)

setGeneric("merge", function(.Object, mergevector)
  standardGeneric("merge"))

setMethod("merge", 
  signature = "stsubpop",
  definition = function(.Object, mergevector) {
    subp.new <- new("stsubpop")

    if (.Object@nsubpop != length(mergevector)) stop("invalid merge vector length.")
    if (.Object@nsubpop == 1) {
      # if there is only one subpopulation left, there is nothing to merge
      subpop.new <- .Object
    } else {
      v1 <- c(mergevector[1], mergevector[-length(mergevector)])
      run <- c(1, which(v1 != mergevector), length(mergevector) + 1)
      merge <- mergevector[1]
    
      subpop.matrix <- NULL

      nsubpop.new <- 0
      for (i in 1:(length(run) - 1)) {
        beg <- run[i]
        end <- run[i + 1] - 1
        if (merge & (beg != end)) {
          subpop.matrix <- cbind(subpop.matrix, apply(.Object@subpop[, beg:end], 1, max))
          nsubpop.new <- nsubpop.new + 1
        } else {
          subpop.matrix <- cbind(subpop.matrix, .Object@subpop[, beg:end])
          nsubpop.new <- nsubpop.new + end - beg + 1
        }
        merge <- !merge
      }

      medianz <- rep(0, nsubpop.new)
      minz    <- rep(0, nsubpop.new)
      minc    <- rep(0, nsubpop.new)
      maxz    <- rep(0, nsubpop.new)
      maxc    <- rep(0, nsubpop.new)

      covariate <- .Object@colvar
      for (i in 1:nsubpop.new) {
        subpop.cov <- covariate[subpop.matrix[, i] == 1]
        medianz[i] <- round((median(subpop.cov)), digits = 2)
        minc[i]    <- min(subpop.cov)
        minz[i]    <- round(minc[i],digits=4)
        maxc[i]    <- max(subpop.cov)
        maxz[i]    <- round(maxc[i],digits=4)
      }

      subp.new@win     <- .Object@win
      subp.new@colvar  <- .Object@colvar
      subp.new@nsubpop <- nsubpop.new
      subp.new@npatsub <- apply(subpop.matrix,2,sum)
      subp.new@subpop  <- subpop.matrix
      subp.new@medianz <- medianz
      subp.new@minz    <- minz
      subp.new@minc    <- minc
      subp.new@maxz    <- maxz
      subp.new@maxc    <- maxc
      subp.new@neventsubTrt0 <- .Object@neventsubTrt0
      subp.new@neventsubTrt1 <- .Object@neventsubTrt1
      subp.new@init    <- TRUE
    }
    return(subp.new)
  }
)

setGeneric("edge", function(.Object, j, side)
  standardGeneric("edge"))  

setMethod("edge", 
  signature = "stsubpop",
  definition = function(.Object, j, side) {
    if (side=="L" & j ==1) stop("invalid arg j.")

    # sort the covariate value
    covariate <- .Object@colvar
    Z <- sort(covariate)
    ZU <- unique(Z)

    if (side == "L") {
      left.edge  <- .Object@minc[j-1]
      left.i   <- (which(ZU == left.edge))[1]
      rside  <- which(ZU == .Object@maxc[j-1])
      right.i    <- rside[length(rside)]
      n      <- min(right.i - left.i + 1, length(ZU)-right.i+1)
      start      <- left.i
    } else if (side == "R") {
      right.edge <- .Object@maxc[j+1]
      rside  <- which(ZU == right.edge)
      right.i  <- rside[length(rside)]
      left.i   <- (which(ZU == .Object@minc[j+1]))[1]
      n      <- min(right.i-left.i+1,left.i)
      start  <- left.i-n+1
    } else {
      stop("unknown side.")
    }

    # create a special stepp subpopulation just for this edge
    subpop.matrix <- matrix(rep(0, (length(covariate) * n)), ncol = n)
    medianz <- rep(0, n)
    minz    <- rep(0, n)
    minc    <- rep(0, n)
    maxz    <- rep(0, n)
    maxc    <- rep(0, n)

    # print(paste(side, "j=",j, "start=",start, "n=",n,ZU[start], ZU[start+n-1]))
    # print(ZU)
    
    for (i in 1:n) {
      # print(c(i, ZU[start + i - 1], ZU[start + i + n - 1]))
      subpop.matrix[,i ] <- (covariate >= ZU[start + i - 1]) * (covariate <= ZU[start + i+n - 2])

      subpop.cov <- covariate[subpop.matrix[, i] == 1]
      medianz[i] <- round((median(subpop.cov)), digits = 2)
      minc[i]    <- min(subpop.cov)
      minz[i]    <- round(minc[i],digits=4)
      maxc[i]    <- max(subpop.cov)
      maxz[i]    <- round(maxc[i],digits=4)
    }

    subp.new         <- new("stsubpop")
    subp.new@win     <- .Object@win
    subp.new@colvar  <- .Object@colvar
    subp.new@nsubpop <- n
    subp.new@subpop  <- subpop.matrix
    subp.new@npatsub <- apply(subpop.matrix,2,sum)
    subp.new@medianz <- medianz
    subp.new@minz    <- minz
    subp.new@minc    <- minc
    subp.new@maxz    <- maxz
    subp.new@maxc    <- maxc
    subp.new@neventsubTrt0 <- .Object@neventsubTrt0
    subp.new@neventsubTrt1 <- .Object@neventsubTrt1
    subp.new@init    <- TRUE
    
    return(subp.new)
  }
)

setMethod("summary", 
  signature = "stsubpop",
  definition = function(object) {
    summary(object@win)
    if (object@init) {
      if (object@win@type == "sliding_events") {
        write(paste("Number of subpopulations created:", object@nsubpop), file = "")
        cat("\n")
        write("Subpopulation summary information", file = "")
        nper <- apply(object@subpop, 2, sum)
        temp <- matrix(c(1:object@nsubpop, object@medianz, object@minz, object@maxz, nper,
          object@neventsubTrt0, object@neventsubTrt1), ncol = 7)
        write("                         Covariate Summary               Sample        Type 1 Events", file = "")
        write(" Subpopulation    Median      Minimum      Maximum         Size  Trt Group 1   Trt Group 2", file = "")
        for (i in 1:object@nsubpop) {
          write(paste(format(temp[i, 1], width = 8), format(temp[i, 2], width = 15, nsmall = 2),
            format(temp[i, 3], width = 12, nsmall = 4), format(temp[i, 4], width = 12, nsmall = 4),
            format(temp[i, 5], width = 12), format(temp[i, 6], width = 12),
            format(temp[i, 7], width = 13)), file = "")
        }
      } else {
        write(paste("Number of subpopulations created:", object@nsubpop), file = "")
        cat("\n")
        write("Subpopulation summary information (including all treatments)", file = "")
        nper <- apply(object@subpop, 2 ,sum)
        temp <- matrix(c(1:object@nsubpop, object@medianz, object@minz, object@maxz, nper), ncol = 5)
        write("                                  Covariate Summary                 Sample", file = "")
        write("     Subpopulation        Median       Minimum       Maximum          size", file = "")
        for (i in 1:object@nsubpop) {
          if (object@win@type == "tail-oriented" & i == length(object@win@r1) + 1 & length(object@win@r1) > 1) {
            write(paste(format(temp[i, 1], width = 12), format(temp[i, 2], width = 19, nsmall = 2),
              format(temp[i, 3], width = 13, nsmall = 4), format(temp[i, 4], width = 13, nsmall = 4),
              format(temp[i, 5], width = 13), "(entire cohort)"), file = "")
          } else if (object@win@type == "tail-oriented" & i == 1 & length(object@win@r2) > 1) {
            write(paste(format(temp[i, 1], width = 12), format(temp[i, 2], width = 19, nsmall = 2),
              format(temp[i, 3], width = 13, nsmall = 4), format(temp[i, 4], width = 13, nsmall = 4),
              format(temp[i, 5], width = 13), "(entire cohort)"), file = "")
          } else {
            write(paste(format(temp[i, 1], width = 12), format(temp[i, 2], width = 19, nsmall = 2),
              format(temp[i, 3], width = 13, nsmall = 4), format(temp[i, 4], width = 13, nsmall = 4),
              format(temp[i, 5], width = 13)), file = "")
          }
        }
      }
    } else {
      write("Subpopulations have not been generated yet.", file = "")
    }
  }
)

# constructor and accessor functions for stepp subpopulation
# 1. constructor
stepp.subpop <- function(swin, cov, coltype = NULL, coltrt = NULL, trts = NULL, minsubpops = NULL) {
  subp <- new("stsubpop")    
  subp <- generate(subp, win = swin, covariate = cov, coltype = coltype, coltrt = coltrt, trts = trts,
    minsubpops = minsubpops)
  return(subp)
}

# 2. merge the subpopulation according to the merge vector mv
stepp.merge.subpop <- function(subpop, mv) {
  subp.new <- merge(subpop, mv)
  return(subp.new)
}

#3. create the edge subpopulation for the "L" or "R" side of the jth subgroup
stepp.edge.subpop <- function(subpop, j, side) {
  subp.new <- edge(subpop, j, side)
  return(subp.new)
}
