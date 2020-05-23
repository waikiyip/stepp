#################################################################
#
# stsubpop.R
#
#######################
# stepp subpopulation #
#######################

# supporting routine to generate subpopulation based just on the coviarate and window
generate.all <- function(sp, win, covariate){
  r1 <- win@r1
  r2 <- win@r2

  if (win@type == "tail-oriented"){
    zvals   <- unique(sort(covariate))
    if (min(r1) < min(zvals) | max(r2) > max(zvals)) stop("Cannot create tail-oriented window")
    nsubpop <- length(r1)+length(r2)+1
    npatsub <- rep(0, nsubpop)
    I0      <- rep(1, nsubpop)
    I1      <- rep(length(zvals), nsubpop)
    if (length(r1) != 0){
      for (i in 1:length(r1)){
	  npatsub[i] <- sum(covariate <= r1[i])
	  sel 	 <- which(zvals <= r1[i])
	  I1[i]      <- sel[length(sel)]
	}
    }
    npatsub[length(r1)+1] <- length(covariate)
    if (length(r2) != 0){
      for (i in 1:length(r2)){
	  k <- i+length(r1)+1
	  npatsub[k] <- sum(covariate >= r2[i])
	  I0[k]	 <- which(zvals >= r2[i])[1]
	}
    }
  }
  else {
    # get the overlapping subgroups
    #
    zvals    <- rep(0, length(covariate))
    absfreq  <- rep(0, length(covariate))
    sortedz  <- sort(covariate)
    zvals[1] <- sortedz[1]
    j <- 1
    absfreq[1] <- 1
    for (i in 2:length(covariate)) {
      if (sortedz[i] != zvals[j]) {
        j <- j + 1
        zvals[j] <- sortedz[i]
        absfreq[j] <- 1
      }
      else {
        absfreq[j] <- absfreq[j] + 1
      }
    }
    zvals <- zvals[1:j]
    absfreq <- absfreq[1:j]
    cumfreq <- absfreq
    for (i in 2:length(cumfreq)) cumfreq[i] <- cumfreq[i] + cumfreq[i - 1]

    I0 <- rep(0, 1000)	# max size of I0 and I1 is 1000; these limit the 
    I1 <- rep(0, 1000)	# max no. of stepp subpopulation to 1000 also.
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
    I0	    <- I0[1:nsubpop]
    I1 	    <- I1[1:nsubpop]
  }
  medians <- rep(0,  nsubpop)
  minz    <- rep(NA, nsubpop)
  minc    <- rep(NA, nsubpop)
  maxz    <- rep(NA, nsubpop)
  maxc    <- rep(NA, nsubpop)

  npats   <- length(covariate)
  subpop  <- matrix(rep(0, (npats * nsubpop)), ncol = nsubpop)
  for (i in 1:nsubpop) {
    subpop[, i] <- (covariate >= zvals[I0[i]]) * (covariate <= zvals[I1[i]])
    medians[i]  <- round((median(covariate[subpop[, i] == 1])), digits = 2)
    minz[i]     <- round(zvals[I0[i]],digits=4)
    minc[i]     <- zvals[I0[i]]
    maxz[i]     <- round(zvals[I1[i]],digits=4)
    maxc[i]     <- zvals[I1[i]]
  }

  # update the object
  sp@win	 <- win
  sp@colvar  <- covariate
  sp@nsubpop <- nsubpop
  sp@subpop  <- subpop
  sp@npatsub <- npatsub
  sp@medianz <- medians
  sp@minz    <- minz
  sp@minc    <- minc
  sp@maxz    <- maxz
  sp@maxc    <- maxc
  sp@init    <- TRUE

  return(sp)

}


setClass("stsubpop",
	   representation(win      = "stwin",	# stepp window object
				colvar   = "numeric",   # vector of covariate of interest (V) 
				nsubpop  = "numeric",	# number of subpopulation generated
				subpop   = "ANY",   	# matrix of subpopulations
				npatsub  = "numeric",	# count of each subpopulation
				medianz  = "numeric",	# median of V for each subpopulation
				minz     = "numeric",	# minimum of V for each subpopulation, round to 4 digits
				maxz	   = "numeric",	# maximum of V for each subpopulation, round to 4 digits
				minc	   = "numeric",	# minimum of V for each subpopulation, actual
				maxc	   = "numeric",   # maximum of V for each subpopulation, actual
				init	   = "logical"	# initialized
				)
	   )

setMethod("initialize", "stsubpop",
	function(.Object){
		.Object@init    <- FALSE
		return(.Object)
	}
)


setValidity("stsubpop", 
	function(object){
	  if (!is.numeric(object@colvar) || length(object@colvar) < object@win@r2){
		print("Invalid argument")
		return (FALSE)
	  }
	  return(TRUE)
     }
)

setGeneric("generate", function(.Object, win, covariate)
		standardGeneric("generate"))	

setMethod("generate", 
	    signature="stsubpop",
	    definition=function(.Object, win, covariate){

		.Object <- generate.all (.Object, win, covariate)
	      return(.Object)
	  }
)


setGeneric("merge", function(.Object, mergevector)
		standardGeneric("merge"))	

setMethod("merge", 
	    signature="stsubpop",
	    definition=function(.Object, mergevector){
		subp.new    <- new("stsubpop")

		if (.Object@nsubpop != length(mergevector)) stop("invalid merge vector length")
		if (.Object@nsubpop == 1){
		  # if there is only one subpopulation left, there is nothing to merge
		  subpop.new <- .Object
		}
		else {

		  v1    <- c(mergevector[1], mergevector[-length(mergevector)])
		  run   <- c(1,which(v1 != mergevector),length(mergevector)+1)
		  merge <- mergevector[1]
		
		  subpop.matrix <- NULL

		  nsubpop.new <- 0
		  for (i in 1:(length(run)-1)){
		    beg <- run[i]
		    end <- run[i+1]-1
		    if (merge & (beg != end)){
			subpop.matrix <- cbind(subpop.matrix, apply(.Object@subpop[,beg:end],1,max))
			nsubpop.new <- nsubpop.new + 1
		    }
		    else {
			subpop.matrix <- cbind(subpop.matrix, .Object@subpop[,beg:end])
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

		  subp.new@win 	   <- .Object@win
		  subp.new@colvar    <- .Object@colvar
		  subp.new@nsubpop   <- nsubpop.new
	  	  subp.new@npatsub   <- apply(subpop.matrix,2,sum)
		  subp.new@subpop	   <- subpop.matrix
		  subp.new@medianz   <- medianz
		  subp.new@minz	   <- minz
		  subp.new@minc	   <- minc
		  subp.new@maxz	   <- maxz
		  subp.new@maxc	   <- maxc
		  subp.new@init	   <- TRUE
		}
	  	return(subp.new)
	    }
)


setGeneric("edge", function(.Object, j, side)
		standardGeneric("edge"))	

setMethod("edge", 
	    signature="stsubpop",
	    definition=function(.Object, j, side){

		if (side=="L" & j ==1) stop("invalid arg j")

		# sort the covariate value
		covariate <- .Object@colvar
		Z <- sort(covariate)
		ZU <- unique(Z)

		if (side == "L"){
		  left.edge  <- .Object@minc[j-1]
		  left.i	 <- (which(ZU == left.edge))[1]
		  rside	 <- which(ZU == .Object@maxc[j-1])
		  right.i    <- rside[length(rside)]
		  n 		 <- min(right.i - left.i + 1, length(ZU)-right.i+1)
		  start      <- left.i
		}
		else 
		if (side == "R"){
		  right.edge <- .Object@maxc[j+1]
		  rside	 <- which(ZU == right.edge)
		  right.i	 <- rside[length(rside)]
		  left.i	 <- (which(ZU == .Object@minc[j+1]))[1]
		  n 		 <- min(right.i-left.i+1,left.i)
		  start	 <- left.i-n+1
		}
		else
		stop("unknown side")

		# create a special stepp subpopulation just for this edge
	      subpop.matrix <- matrix(rep(0, (length(covariate) * n)), ncol = n)
	      medianz <- rep(0, n)
	      minz    <- rep(0, n)
	      minc    <- rep(0, n)
	      maxz    <- rep(0, n)
	      maxc    <- rep(0, n)

		#print(paste(side,"j=",j,"start=",start,"n=",n,ZU[start], ZU[start+n-1]))
		#print(ZU)
		
		for (i in 1:n){
		  #print(c(i, ZU[start+i-1], ZU[start+i+n-1]))
		  subpop.matrix[,i] <- (covariate >= ZU[start+i-1]) * (covariate <= ZU[start+i+n-2])

              subpop.cov <- covariate[subpop.matrix[, i] == 1]
              medianz[i] <- round((median(subpop.cov)), digits = 2)
		  minc[i]    <- min(subpop.cov)
              minz[i]    <- round(minc[i],digits=4)
		  maxc[i]    <- max(subpop.cov)
              maxz[i]    <- round(maxc[i],digits=4)

		}

		subp.new		   <- new("stsubpop")		
		subp.new@win 	   <- .Object@win
		subp.new@colvar      <- .Object@colvar
		subp.new@nsubpop     <- n
		subp.new@subpop	   <- subpop.matrix
	  	subp.new@npatsub     <- apply(subpop.matrix,2,sum)
		subp.new@medianz     <- medianz
		subp.new@minz	   <- minz
		subp.new@minc	   <- minc
		subp.new@maxz	   <- maxz
		subp.new@maxc	   <- maxc
		subp.new@init	   <- TRUE
		
	    return(subp.new)
	  }
)
	 

setMethod("summary", 
	    signature="stsubpop",
	    definition=function(object){
		summary(object@win)
		if (object@init){
		  write(paste("      Number Of Subpopulations Created :", object@nsubpop),file="")
    		  cat("\n")
    		  write("Subpopulation Summary Information (including all treatments)",file="")
    		  nper <- apply(object@subpop,2,sum)
    		  temp <- matrix(c(1:object@nsubpop,object@medianz,object@minz,object@maxz,nper),ncol=5)
    		  write("                                  Covariate Summary                  Sample",file="")
    		  write("     Subpopulation        Median       Minimum       Maximum          Size",file="")
    		  for (i in 1:object@nsubpop) {
		    if (object@win@type=="tail-oriented" & i == length(object@win@r1)+1){ 
       	      write(paste(format(temp[i,1],width=12),format(temp[i,2],width=19,nsmall=2),
          	     		  	format(temp[i,3],width=13,nsmall=4),format(temp[i,4],width=13,nsmall=4),
          	    		  	format(temp[i,5],width=13), "Entire Cohort"),file="")
		    } else {
       	      write(paste(format(temp[i,1],width=12),format(temp[i,2],width=19,nsmall=2),
          	    	  		format(temp[i,3],width=13,nsmall=4),format(temp[i,4],width=13,nsmall=4),
          	    			format(temp[i,5],width=13)),file="")
		    }
	        }
		}
		else
		  write("      Subpopulations have not been generated", file="")
	    }	
)

# constructor and accessor functions for stepp subpopulation
# 1. constructor
stepp.subpop <- function(swin, cov){
	subp    <- new("stsubpop")		
	subp    <- generate(subp, win=swin, covariate=cov)
	return(subp)
}

# 2. merge the subpopulation according to the merge vector mv
stepp.merge.subpop <- function(subpop, mv){
	subp.new <- merge(subpop, mv)
	return(subp.new)
}

#3. create the edge subpopulation for the "L" or "R" side of the jth subgroup
stepp.edge.subpop <- function(subpop, j, side){
	subp.new <- edge(subpop, j, side)
	return(subp.new)
}

