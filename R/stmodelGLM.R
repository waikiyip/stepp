#################################################################
#
# stmodelGLM.R
#
########################
# stepp model: GLM     #
########################
####################
# function: FiellerSE (ns, nt, r, mt, se.s, se.t)
#	returns the SE based on Fieller theorem for ratio of two independent normal RVs.
#	
#	ns - number of samples for the numerator
#	nt - number of samples for the denominator
#	r  - estimates of the ratio
#	mt - mean of the denominator samples
#	se.s - standard error for the numerator
#	se.t - standard error for the denominrator
#
FiellerSE <- function(ns, nt, r, mt, se.s, se.t){
  df     <- nt+ns-2
  t0.975 <- qt(0.975, df)
  sp     <- sqrt(((ns-1)*(se.s^2)+(nt-1)*(se.t^2))/df)
  g      <- (t0.975*sp)^2/(nt*(mt^2)) 
  if ((1-g)<0.0001) res <- NaN
  else {
    l1 <- (1/(1-g))*(r+t0.975*sp*sqrt((1-g)/ns+r^2/nt)/mt)
    l2 <- (1/(1-g))*(r-t0.975*sp*sqrt((1-g)/ns+r^2/nt)/mt)

    res <- (abs(l1-l2)/2)*t0.975
  }

  return(res) 

}


setClass("stmodelGLM",
  representation(
    coltrt   = "numeric",   # treatment
    colY     = "numeric",   # outcome 
    trts     = "numeric",   # trt encoding
    MM       = "ANY",       # model matrix
    glm      = "character", # glm type
    link     = "character"  # link function
  ),
  prototype(
   coltrt = 0,
   colY   = 0,
   trts   = 0,
   MM     = NULL,
   glm    = "",
   link   = ""
  )
)

setMethod("initialize", "stmodelGLM",
  function(.Object, coltrt, colY, trts, MM, glm, link, ...) {
    if (missing(colY)) {
       colY <- numeric()
    }
    if (missing(MM)) {
       MM <- NULL
    }
    if (missing(glm)) {
      glm <- "gaussian"
    }
    if (missing(link)) {
      if (glm == "gaussian") {
        link <- "identity"
      } else if (glm == "binomial") {
        link <- "logit"
      } else if (glm == "poisson") {
        link <- "log"
      }
    }
    if (missing(coltrt)) {
      coltrt <- rep(1, length(colY))
    }
    if (missing(trts)) {
      trts <- 1
    }
    .Object@coltrt <- coltrt
    .Object@colY <- colY
    .Object@trts <- trts
    .Object@MM <- MM
    .Object@glm <- glm
    .Object@link <- link
    if (!validObject(.Object)) stop("")
    callNextMethod(.Object, ...)
  }
)

setValidity("stmodelGLM", 
  function(object) {
    status <- TRUE

    # add conditions to verify

    return(status)
  }
)

setMethod("estimate", 
  signature="stmodelGLM",
  definition=function(.Object, sp, ...) {
    modtype   <- sp@win@type
    if (modtype == "sliding_events") {
      stop("Currently event-based sliding windows are available only for competing risks analyses.")
    }
    nsubpop   <- sp@nsubpop
    subpop    <- sp@subpop
    coltrt    <- .Object@coltrt
    trts      <- .Object@trts
    OuTcOmE   <- .Object@colY
    covariate <- .Object@MM
    glm.type  <- .Object@glm
    link      <- .Object@link

    if (glm.type!="gaussian"& glm.type!="binomial" & glm.type!="poisson") stop ("Unknown model!")
    ntrts     <- length(trts) 
    txassign  <- rep(NA, length(coltrt))

    for (j in 1:ntrts) txassign[which(coltrt == trts[j])] <- j

    Obs  <- list(sObs = rep(0, nsubpop),
                 sSE  = rep(0, nsubpop),
                 oObs = 0,
                 oSE  = 0)
    TrtEff <- array(list(Obs),ntrts)

    HRs    <- list(skmw     = 0,
                   logHR    = rep(0, nsubpop),
                   logHRSE  = rep(0, nsubpop),
                   ologHR   = 0,
                   ologHRSE = 0,
                   logHRw   = 0)
    Ratios <- array(list(HRs), ntrts - 1)

    ## Create a formula for the model:
    if (!is.null(covariate)) {
      has.intercept <- colnames(covariate)[1]=="(Intercept)"
      if (has.intercept) {
        cstart <- 2
        xnam <- c("txassign", colnames(covariate)[cstart:dim(covariate)[2]])
        fmla <- as.formula(paste("OuTcOmE ~ ", paste(xnam, collapse= "+")))
      } else {
        cstart <- 1
        xnam <- c("txassign", colnames(covariate)[cstart:dim(covariate)[2]])
        fmla <- as.formula(paste("OuTcOmE ~ 0+", paste(xnam, collapse= "+")))  
      }     
    } else {
      cstart <- 2 # default model has an intercept
      fmla <- as.formula("OuTcOmE ~ txassign")
    }

    if (ntrts > 1) {
      for (j in 2:ntrts) {
        DIFF       <- rep(0, nsubpop)
        DIFFSE     <- rep(0, nsubpop)

        for (i in 1:nsubpop) {

	  { # 1. Treatment Effect for the i subpopulation for 1st and jth treatment
	    #    teffect(subpop, glm.type, fmla, txassign, covariate, 1, j, i)
          seli <- (txassign == 1 | txassign == j) & (subpop[, i] == 1)
 
          if (glm.type == "gaussian") {
            m1 <- glm(fmla, subset=seli, family=gaussian(link="identity"))
          } else if (glm.type == "binomial") {
            m1 <- glm(fmla, subset=seli, family=binomial(link="logit"))
          } else if (glm.type == "poisson") {
            m1 <- glm(fmla, subset=seli, family=poisson(link="log"))
          }

          # use predict
          if (!is.null(covariate)) {
            mm1  <- covariate[which(seli & txassign==1), , drop = FALSE]
            mm1  <- cbind(txassign[which(seli & txassign==1)], mm1)
            cm   <- apply(mm1, 2, mean)
            if (has.intercept) cm <- cm[-2]	# WKY: remove intercept column, as we do use intercept in predict
            names(cm)[1] <- "txassign"
            cm[1]<- 1 # treatment indicator is 1
            p1   <- predict(m1, newdata=data.frame(t(cm)), se.fit=TRUE)

            mm1  <- covariate[which(seli & txassign==j), , drop = FALSE]
            mm1  <- cbind(txassign[which(seli & txassign==j)], mm1)
            cm   <- apply(mm1, 2, mean)
            if (has.intercept) cm <- cm[-2]	# WKY: remove intercept column, as we do use intercept in predict
            names(cm)[1] <- "txassign"
		cm[1]<- j # treatment indicator is j
            p2   <- predict(m1, newdata=data.frame(t(cm)), se.fit=TRUE)
          } else {
            p1   <- predict(m1, newdata=data.frame(txassign=1), se.fit=TRUE)
            p2   <- predict(m1, newdata=data.frame(txassign=j), se.fit=TRUE)
          }

	  } # treatment effect calculation

          if (glm.type == "gaussian") {
            TrtEff[[1]]$sObs[i]    <- p1$fit
            TrtEff[[1]]$sSE[i]     <- p1$se
            TrtEff[[j]]$sObs[i]    <- p2$fit
            TrtEff[[j]]$sSE[i]     <- p2$se

		DIFF[i]		     <- TrtEff[[1]]$sObs[i] - TrtEff[[j]]$sObs[i]
		DIFFSE[i]		     <- sqrt(TrtEff[[1]]$sSE[i]^2+TrtEff[[j]]$sSE[i]^2)
		
            Ratios[[j-1]]$logHR[i] <- TrtEff[[1]]$sObs[i] / TrtEff[[j]]$sObs[i]
            # use Fieller Theorem
            Ratios[[j-1]]$logHRSE[i] <- FiellerSE (length(which(seli & txassign==1)),	# ns
								   length(which(seli & txassign==j)),	# nt
								   Ratios[[j-1]]$logHR[i], 			# r
								   TrtEff[[j]]$sObs[i],				# mt
								   TrtEff[[1]]$sSE[i],				# se.s
								   TrtEff[[j]]$sSE[i])				# se.t
		
          } else if (glm.type == "binomial") {
            exp1          	     <- exp(p1$fit)
            TrtEff[[1]]$sObs[i]    <- exp1/(1+exp1)
            TrtEff[[1]]$sSE[i]     <- p1$se*(exp1/((1+exp1)^2))  # delta method
            exp2                   <- exp(p2$fit)
            TrtEff[[j]]$sObs[i]    <- exp2/(1+exp2)
            TrtEff[[j]]$sSE[i]     <- p2$se*(exp2/((1+exp2)^2))  # delta method

            DIFF[i]                <- TrtEff[[1]]$sObs[i] - TrtEff[[j]]$sObs[i]
            DIFFSE[i]              <- sqrt(TrtEff[[1]]$sSE[i]^2 + TrtEff[[j]]$sSE[i]^2)
            Ratios[[j-1]]$logHR[i] <- log(TrtEff[[1]]$sObs[i]) - log(TrtEff[[j]]$sObs[i])
            # use delta method
            Ratios[[j-1]]$logHRSE[i] <- sqrt((TrtEff[[1]]$sSE[i]^2)/(TrtEff[[1]]$sObs[i]^2) +
                                             (TrtEff[[j]]$sSE[i]^2)/(TrtEff[[j]]$sObs[i]^2))
          } else if (glm.type == "poisson") {
            TrtEff[[1]]$sObs[i]    <- exp(p1$fit)
            TrtEff[[1]]$sSE[i]     <- p1$se*exp(p1$fit)  # delta method

            TrtEff[[j]]$sObs[i]    <- exp(p2$fit)
            TrtEff[[j]]$sSE[i]     <- p2$se*exp(p2$fit)  # delta method

            DIFF[i]                <- TrtEff[[1]]$sObs[i] - TrtEff[[j]]$sObs[i]
            DIFFSE[i]              <- sqrt(TrtEff[[1]]$sSE[i]^2 + TrtEff[[j]]$sSE[i]^2)
            Ratios[[j-1]]$logHR[i] <- log(TrtEff[[1]]$sObs[i]) - log(TrtEff[[j]]$sObs[i])
            # use delta method
            Ratios[[j-1]]$logHRSE[i] <- sqrt((TrtEff[[1]]$sSE[i]^2)/(TrtEff[[1]]$sObs[i]^2) +
                                             (TrtEff[[j]]$sSE[i]^2)/(TrtEff[[j]]$sObs[i]^2))
          }
        }  # for each subpopulation

        Ratios[[j-1]]$skmw   <- sum(DIFF/DIFFSE)
        Ratios[[j-1]]$logHRw <- sum(Ratios[[j-1]]$logHR/Ratios[[j-1]]$logHRSE)

    { # 2. Overall Treatment Effect comparing 1st and jth treatment
        selall <- (txassign == 1 | txassign == j)

      if (glm.type == "gaussian") {
        m1all <- glm(fmla, subset=selall, family=gaussian(link="identity"))
      } else if (glm.type == "binomial") {
        m1all <- glm(fmla, subset=selall, family=binomial(link="logit"))
      } else if (glm.type == "poisson") {
        m1all <- glm(fmla, subset=selall, family=poisson(link="log"))
      }

      # use predict
      if (!is.null(covariate)) {
        mm1  <- covariate[which(selall & txassign==1), , drop = FALSE]
        mm1  <- cbind(txassign[which(selall & txassign==1)], mm1)
        cm   <- apply(mm1, 2, mean)
        if (has.intercept) cm <- cm[-2]	# WKY: remove intercept column, as we do use intercept in predict
        names(cm)[1] <- "txassign"
        cm[1]<- 1 # treatment indicator is 1
        p1   <- predict(m1all, newdata=data.frame(t(cm)), se.fit=TRUE)

        mm1  <- covariate[which(selall & txassign==j), , drop = FALSE]
        mm1  <- cbind(txassign[which(selall & txassign==j)], mm1)
        cm   <- apply(mm1, 2, mean)
        if (has.intercept) cm <- cm[-2]	# WKY: remove intercept column, as we do use intercept in predict
        names(cm)[1] <- "txassign"
        cm[1]<- j # treatment indicator is j
        p2   <- predict(m1all, newdata=data.frame(t(cm)), se.fit=TRUE)
      } else {
        p1   <- predict(m1all, newdata=data.frame(txassign=1), se.fit=TRUE)
        p2   <- predict(m1all, newdata=data.frame(txassign=j), se.fit=TRUE)
      }
    } # end of overall treatment effect

        if (glm.type == "gaussian") {
          TrtEff[[1]]$oObs  <- p1$fit
          TrtEff[[1]]$oSE   <- p1$se
          TrtEff[[j]]$oObs  <- p2$fit
          TrtEff[[j]]$oSE   <- p2$se

          Ratios[[j-1]]$ologHR   <- TrtEff[[1]]$oObs / TrtEff[[j]]$oObs

	    # Fieller Theorem
          Ratios[[j-1]]$ologHRSE <- FiellerSE (length(which(selall & txassign==1)),	# ns
							     length(which(selall & txassign==j)),	# nt
							     Ratios[[j-1]]$ologHR, 			# r
							     TrtEff[[j]]$oObs,				# mt
							     TrtEff[[1]]$oSE,				# se.s
							     TrtEff[[j]]$oSE)				# se.t

          # delta method
          # Ratios[[j-1]]$ologHRSE <- abs(Ratios[[j-1]]$ologHR)*
					# 	sqrt((TrtEff[[1]]$oSE^2)/(TrtEff[[1]]$oObs^2)
          # 						+(TrtEff[[j]]$oSE^2)/(TrtEff[[j]]$oObs^2))  

	    #print("Debug - Overall")
	    #print(c(p1$se,p2$se))
	    #print((TrtEff[[1]]$oSE^2)/(TrtEff[[1]]$oObs^2)
          #						+(TrtEff[[j]]$oSE^2)/(TrtEff[[j]]$oObs^2))  
	    #print("-Debug")
        } else if (glm.type == "binomial") {
          oexp1             <- exp(p1$fit)
          TrtEff[[1]]$oObs  <- oexp1/(1+oexp1)
          TrtEff[[1]]$oSE   <- p1$se*(oexp1/((1+oexp1)^2)) # delta method
          oexp2             <- exp(p2$fit)
          TrtEff[[j]]$oObs  <- oexp2/(1+oexp2)
          TrtEff[[j]]$oSE   <- p2$se*(oexp2/((1+oexp2)^2)) # delta method

          Ratios[[j-1]]$ologHR   <- log(TrtEff[[1]]$oObs) - log(TrtEff[[j]]$oObs)
          # delta method
          Ratios[[j-1]]$ologHRSE <- sqrt((TrtEff[[1]]$oSE^2)/(TrtEff[[1]]$oObs^2)
          						+(TrtEff[[j]]$oSE^2)/(TrtEff[[j]]$oObs^2))  
        } else if (glm.type == "poisson") {
          TrtEff[[1]]$oObs  <- exp(p1$fit)
          TrtEff[[1]]$oSE   <- p1$se*exp(p1$fit)  # delta method
          TrtEff[[j]]$oObs  <- exp(p2$fit)
          TrtEff[[j]]$oSE   <- p2$se*exp(p2$fit)  # delta method

          Ratios[[j-1]]$ologHR   <- log(TrtEff[[1]]$oObs) - log(TrtEff[[j]]$oObs)
          # delta method
          Ratios[[j-1]]$ologHRSE <- sqrt((TrtEff[[1]]$oSE^2)/(TrtEff[[1]]$oObs^2)
          						+(TrtEff[[j]]$oSE^2)/(TrtEff[[j]]$oObs^2))  
        }

        if (sum(is.na(TrtEff[[1]]$sObs)) != 0 | sum(is.na(TrtEff[[j]]$sObs)) != 0 | 
          is.na(TrtEff[[1]]$oObs) != FALSE | is.na(TrtEff[[j]]$oObs) != FALSE) {
          cat("\n")
          print(paste("Unable to estimate the effect because there are too few events within one or more subpopulation(s)."))
          print(paste("The problem may be avoided by constructing larger subpopulations."))
          stop()
        }
      } # for each trt j
    }
      
    if (ntrts == 1) {
      if (!is.null(covariate)) {
        has.intercept <- colnames(covariate)[1]=="(Intercept)"
        if (has.intercept) {
          cstart <- 2
          xnam <- colnames(covariate)[cstart:dim(covariate)[2]]
          fmla <- as.formula(paste("OuTcOmE ~ ", paste(xnam, collapse= "+")))
        } else {
          cstart <- 1
          xnam <- colnames(covariate)[cstart:dim(covariate)[2]]
          fmla <- as.formula(paste("OuTcOmE ~ 0+", paste(xnam, collapse= "+")))  
        }     
      } else {
        fmla <- as.formula("OuTcOmE ~ 1")
      }

      for (i in 1:nsubpop) {
        seli    <- (txassign==1) & (subpop[,i]==1)

        if (glm.type == "gaussian") {
          m1 <- glm(fmla, subset=seli, family=gaussian(link="identity"))
        } else if (glm.type == "binomial") {
          m1 <- glm(fmla, subset=seli, family=binomial(link="logit"))
        } else if (glm.type == "poisson") {
          m1 <- glm(fmla, subset=seli, family=poisson(link="log"))
        } else {
          stop("Unknown Model !")
        }
        # use predict 
        if (!is.null(covariate)) {
          mm1  <- covariate[which(seli),,drop=FALSE]
          cm   <- apply(mm1, 2, mean)
          if (has.intercept) cm <- cm[-1]
          p1   <- predict(m1, newdata=data.frame(t(cm)), se.fit=TRUE)
        } else {
          p1   <- predict(m1, newdata=data.frame(1), se.fit=TRUE)
        }

        if (glm.type == "gaussian") {
          TrtEff[[1]]$sObs[i]   <- p1$fit
          TrtEff[[1]]$sSE[i]    <- p1$se
        } else if (glm.type == "binomial") {
          exp1          <- exp(p1$fit)
          TrtEff[[1]]$sObs[i]   <- exp1/(1+exp1)
          TrtEff[[1]]$sSE[i]    <- p1$se*(exp1/((1+exp1)^2))  # delta method
        } else if (glm.type == "poisson") {
          TrtEff[[1]]$sObs[i]   <- exp(p1$fit)
          TrtEff[[1]]$sSE[i]    <- p1$se*exp(p1$fit)  # delta method
        }
      }  # for each subpopulation

      if (glm.type == "gaussian") {
        m1all     <- glm(fmla, family=gaussian(link="identity"))
      } else if (glm.type == "binomial") {
        m1all     <- glm(fmla, family=binomial(link="logit"))
      } else if (glm.type == "poisson") {
        m1all     <- glm(fmla, family=poisson(link="log"))
      }

      if (!is.null(covariate)) {
        om1    <- covariate
        onv1   <- apply(om1,2,mean)
        if (has.intercept) onv1 <- onv1[-2]
        op1    <- predict(m1all, newdata=data.frame(t(onv1)), se.fit=TRUE)
      } else {
        op1    <- predict(m1all, newdata=data.frame(1), se.fit=TRUE)
      }

      if (glm.type == "gaussian") {
        TrtEff[[1]]$oObs  <- op1$fit
        TrtEff[[1]]$oSE   <- op1$se
      } else if (glm.type == "binomial") {
        oexp1         <- exp(op1$fit)
        TrtEff[[1]]$oObs  <- oexp1/(1+oexp1)
        TrtEff[[1]]$oSE   <- op1$se*(oexp1/((1+oexp1)^2)) # delta method
      } else if (glm.type == "poisson") {
        TrtEff[[1]]$oObs  <- exp(op1$fit)
        TrtEff[[1]]$oSE   <- op1$se*exp(op1$fit)  # delta method
      }

      if (sum(is.na(TrtEff[[1]]$sObs)) != 0 | is.na(TrtEff[[1]]$oObs) != FALSE) {
        cat("\n")
        print(paste("Unable to estimate the effect because there are too few events within one or more subpopulation(s)."))
        print(paste("The problem may be avoided by constructing larger subpopulations."))
        stop()
      }
    }

    if (glm.type == "gaussian") {
      model = "GLMGe"
    } else if (glm.type == "binomial") {
      model = "GLMBe"
    } else if (glm.type == "poisson") {
      model = "GLMPe"
    }

    estimate <- list(model   = model,
                     ntrts   = ntrts,
                     TrtEff  = TrtEff,
                     Ratios  = Ratios)      
    return(estimate)
  }
)

setMethod("test",
  signature="stmodelGLM",
  definition=function(.Object, nperm, sp, effect, showstatus=TRUE){

  test <- NULL
  if (nperm > 0 & length(.Object@trts) > 1) {
    win         <- sp@win
    nsubpop     <- sp@nsubpop
    osubpop     <- sp@subpop
    coltrt      <- .Object@coltrt
    trts        <- .Object@trts
    OuTcOmE     <- .Object@colY
    covariate   <- .Object@MM
    glm.type    <- .Object@glm
    link        <- .Object@link

    ntrts       <- length(trts)
    txassign    <- rep(-1, length(coltrt)) # do not use NA

    for (j in 1:ntrts) 
      txassign[which(coltrt == trts[j])] <- j

    tsigma    <- matrix(rep(0, (nperm * nsubpop)), ncol = nsubpop)    
    PermTest  <- list(sigma      = tsigma,
                      HRsigma    = tsigma,
                      pvalue     = 0,
                      chi2pvalue = 0,
                      HRpvalue   = 0)
    Res       <- array(list(PermTest),ntrts-1)

    for (j in 2:ntrts){
    #
    #   do the permutations
    # compare trt1 with the rest
    #
      pvalue      <- NA
      differences <- matrix(rep(0, (nperm * nsubpop)), ncol = nsubpop)
      logratios   <- matrix(rep(0, (nperm * nsubpop)), ncol = nsubpop)

      tPerm       <- rep(0, nperm)
      no          <- 0
      p           <- 0
      terminate   <- 0
      Ntemp       <- nrow(osubpop)
      IndexSet1   <- (1:Ntemp)[txassign == 1]
      IndexSet2   <- (1:Ntemp)[txassign == j]

      ## Create a formula for the model:
      ## Cannot handle a no intercept model here
      if (!is.null(covariate)){
        has.intercept <- colnames(covariate)[1]=="(Intercept)"
        if (has.intercept){
          cstart <- 2
          xnam <- c("txassign", colnames(covariate)[cstart:dim(covariate)[2]])
          fmla <- as.formula(paste("OuTcOmE ~ ", paste(xnam, collapse= "+")))
        } else {
          cstart <- 1
          xnam <- c("txassign", colnames(covariate)[cstart:dim(covariate)[2]])
          fmla <- as.formula(paste("OuTcOmE ~ 0+", paste(xnam, collapse= "+")))  
        } 
      } else {
        cstart <- 2 # default model has an intercept
        fmla <- as.formula("OuTcOmE ~ txassign")
      }
      if (showstatus){
        title <- paste("\nComputing the p-value with", nperm)
        title <- paste(title, "permutations comparing trt", trts[j], "with trt", trts[1], "\n")
        cat(title)
        pb <- txtProgressBar(min=0, max=nperm-1, style=3)
      }

	Subpop1 <- as.matrix(osubpop)
      while (no < nperm) {
      ## PERMUTE THE VECTOR OF SUBPOPULATIONS WITHIN TREATMENTS ##
        if (showstatus) setTxtProgressBar(pb, no)

        #Subpop        <- as.matrix(osubpop)
        permuteSubpop <- matrix(0, nrow = Ntemp, ncol = nsubpop)
        permuteSubpop[txassign == 1] <- Subpop1[sample(IndexSet1),]
        permuteSubpop[txassign == j] <- Subpop1[sample(IndexSet2),]
        subpop        <- permuteSubpop
        sglm1         <- rep(0, nsubpop)
        sglm2         <- rep(0, nsubpop)
        sglmSE1       <- rep(0, nsubpop)
        sglmSE2       <- rep(0, nsubpop)
        slogratios    <- rep(0, nsubpop)
        slogratioSEs  <- rep(0, nsubpop)

        for (i in 1:nsubpop) {

	  { # 3. Treatment Effect for the i subpopulation for 1st and jth treatment
	    #    teffect(subpop, glm.type, fmla, txassign, covariate, 1, j, i)
          seli <- (txassign == 1 | txassign == j) & (subpop[, i] == 1)

          if (glm.type == "gaussian") {
            m1 <- glm(fmla, subset=seli, family=gaussian(link="identity"))
          } else if (glm.type == "binomial") {
            m1 <- glm(fmla, subset=seli, family=binomial(link="logit"))
          } else if (glm.type == "poisson") {
            m1 <- glm(fmla, subset=seli, family=poisson(link="log"))
          }
     
          # use predict
          if (!is.null(covariate)) {
            mm1  <- covariate[which(seli & txassign == 1), , drop = FALSE]
            mm1  <- cbind(txassign[which(seli & txassign == 1)], mm1)
            cm   <- apply(mm1, 2, mean)
            if (has.intercept) cm <- cm[-2]	# WKY: remove intercept column, as we do use intercept in predict
            names(cm)[1] <- "txassign"
            cm[1]<- 1 # treatment indicator is 1
            p1   <- predict(m1, newdata=data.frame(t(cm)), se.fit=TRUE)

            mm1  <- covariate[which(seli & txassign == j), , drop = FALSE]
            mm1  <- cbind(txassign[which(seli & txassign == j)], mm1)
            cm   <- apply(mm1, 2, mean)
            if (has.intercept) cm <- cm[-2]	# WKY: remove intercept column, as we do use intercept in predict
            names(cm)[1] <- "txassign"
            cm[1] <- j
            p2   <- predict(m1, newdata=data.frame(t(cm)), se.fit=TRUE)
          } else {
            p1   <- predict(m1, newdata=data.frame(txassign=1), se.fit=TRUE)
            p2   <- predict(m1, newdata=data.frame(txassign=j), se.fit=TRUE)
          }

	  } # treatment effect calculation

	    if (glm.type == "gaussian"){
            sglm1[i]   <- p1$fit
            sglmSE1[i] <- p1$se         
            sglm2[i]   <- p2$fit
            sglmSE2[i] <- p2$se

            slogratios[i]   <- sglm1[i]/sglm2[i]
		# Fieller Theorem
            slogratioSEs[i] <- FiellerSE (length(which(seli & txassign==1)),	# ns
							length(which(seli & txassign==j)),	# nt
							slogratios[i], 				# r
							sglm2[i],					# mt
							sglmSE1[i],					# se.s
							sglmSE2[i])					# se.t

	    } else
          if (glm.type == "binomial"){
            exp1       <- exp(p1$fit)
            sglm1[i]   <- exp1/(1+exp1)
            sglmSE1[i] <- p1$se*(exp1/((1+exp1)^2)) # delta method
            exp2       <- exp(p2$fit)
            sglm2[i]   <- exp2/(1+exp2)
            sglmSE2[i] <- p2$se*(exp2/((1+exp2)^2)) # delta method

            slogratios[i]   <- log(sglm1[i])-log(sglm2[i])
            slogratioSEs[i] <- sqrt((sglmSE1[i]^2)/(sglmSE1[i]^2) +
              				(sglmSE2[i]^2)/(sglmSE2[i]^2))
          } else
          if (glm.type == "poisson"){
		exp1       <- exp(p1$fit)
            sglm1[i]   <- exp1
            sglmSE1[i] <- p1$se*exp1 # delta method
		exp2       <- exp(p2$fit)
            sglm2[i]   <- exp2
            sglmSE2[i] <- p2$se*exp2 # delta method

            slogratios[i]   <- log(sglm1[i])-log(sglm2[i])
            slogratioSEs[i] <- sqrt((sglmSE1[i]^2)/(sglmSE1[i]^2) +
              				(sglmSE2[i]^2)/(sglmSE2[i]^2))
          } else
          stop ("Unsupported GLM !")

        } # for each subpop i

    { # 4. Overall Treatment Effect comparing 1st and jth treatment
        selall <- (txassign == 1 | txassign == j)

      if (glm.type == "gaussian") {
        m1all <- glm(fmla, subset=selall, family=gaussian(link="identity"))
      } else if (glm.type == "binomial") {
        m1all <- glm(fmla, subset=selall, family=binomial(link="logit"))
      } else if (glm.type == "poisson") {
        m1all <- glm(fmla, subset=selall, family=poisson(link="log"))
      }

      # use predict
      if (!is.null(covariate)) {
        mm1  <- covariate[which(selall & txassign == 1), , drop = FALSE]
        mm1  <- cbind(txassign[which(selall & txassign == 1)], mm1)
        cm   <- apply(mm1, 2, mean)
        if (has.intercept) cm <- cm[-2]	# WKY: remove intercept column, as we do use intercept in predict
        names(cm)[1] <- "txassign"
        cm[1]<- 1 # treatment indicator is 1
        p1   <- predict(m1, newdata=data.frame(t(cm)), se.fit=TRUE)

        mm1  <- covariate[which(selall & txassign == j), , drop = FALSE]
        mm1  <- cbind(txassign[which(selall & txassign == j)], mm1)
        cm   <- apply(mm1, 2, mean)
        if (has.intercept) cm <- cm[-2]	# WKY: remove intercept column, as we do use intercept in predict
        names(cm)[1] <- "txassign"
        cm[1] <- j
        p2   <- predict(m1all, newdata=data.frame(t(cm)), se.fit=TRUE)
      } else {
        p1   <- predict(m1all, newdata=data.frame(txassign=1), se.fit=TRUE)
        p2   <- predict(m1all, newdata=data.frame(txassign=j), se.fit=TRUE)
      }
    } # end of overall treatment effect

        if (glm.type == "gaussian"){
          overallSglm1    <- p1$fit       
          overallSglm2    <- p2$fit
          overalllogRatio <- overallSglm1/overallSglm2
        } else
        if (glm.type == "binomial"){
          exp1            <- exp(p1$fit)
          overallSglm1    <- exp1/(1+exp1)        
          exp2            <- exp(p2$fit)
          overallSglm2    <- exp2/(1+exp2)
          overalllogRatio <- log(overallSglm1)-log(overallSglm2)
        } else
        if (glm.type == "poisson"){
          overallSglm1    <- exp(p1$fit)
          overallSglm2    <- exp(p2$fit)
          overalllogRatio <- log(overallSglm1)-log(overallSglm2)
        }

        if (sum(is.na(sglm1)) == 0 & sum(is.na(sglm2)) == 0 & is.na(overallSglm1) == FALSE & is.na(overallSglm2) == FALSE) {
          no <- no + 1
          p  <- p + 1
          for (s in 1:nsubpop) {
            differences[p, s] <- (sglm1[s] - sglm2[s]) - (overallSglm1 - overallSglm2)
            logratios[p,s]    <- slogratios[s] - overalllogRatio
          }
        }

        terminate <- terminate + 1
        if (terminate >= nperm + 10000) {
          print(paste("After permuting ", nperm, "plus 10000, or ", 
                nperm + 10000, " times, the program is unable to generate the permutation distribution based on ", 
                nperm, "permutations of the data"))
          print(paste("Consider creating larger subpopulations or selecting a different timepoint for estimation"))
          stop()
        }
      }
      if (showstatus) close(pb)

      # generating the sigmas and p-values
      sObs   <- effect$TrtEff[[1]]$sObs-effect$TrtEff[[j]]$sObs
      oObs   <- effect$TrtEff[[1]]$oObs-effect$TrtEff[[j]]$oObs
      logHR  <- effect$Ratios[[j-1]]$logHR
      ologHR <- effect$Ratios[[j-1]]$ologHR

      mname  <- paste0("SP",seq(1,dim(differences)[2]),"-Overall")
      if (win@type == "tail-oriented"){
        # remove the trivial case of the full cohort
        rm <- length(win@r1)+1
        differences <- differences[,-rm]
        mname <- mname[-rm]
        sObs <- sObs[-rm]
        oObs <- oObs[-rm]
        logratios <- logratios[,-rm]
        logHR <- logHR[-rm]
        ologHR <- ologHR[-rm]
      }

      sigmas          <- ssigma(differences)
      rownames(sigmas$sigma)  <- mname
      colnames(sigmas$sigma)  <- mname
      Res[[j-1]]$sigma        <- sigmas$sigma
      Res[[j-1]]$chi2pvalue   <- ppv (differences, sigmas$sigmainv, sObs, oObs,nperm)
      Res[[j-1]]$pvalue       <- ppv2(differences, sObs, oObs, nperm)

      sigmas          <- ssigma(logratios)
      rownames(sigmas$sigma)  <- mname
      colnames(sigmas$sigma)  <- mname

      Res[[j-1]]$HRsigma      <- sigmas$sigma
      Res[[j-1]]$HRpvalue     <- ppv2(logratios, logHR, ologHR, nperm)

      } # for each trt j

      test = list(model = "GLMt",
                  ntrts = ntrts,
                  Res   = Res)
    } # nperm == 0 
    return(test)
  } # end of method

)

setGeneric("subgroup", function(.Object, subsample)
    standardGeneric("subgroup"))  

setMethod("subgroup", 
      signature="stmodelGLM",
      definition=function(.Object, subsample){

      #if (model.type == "stmodelKM" | mode.type == "stmodelCOX") {
      #  model.i <- stepp.KM(coltrt    = .Object@model@coltrt[bmatrix[,i]],
    #       survTime  = .Object@model@survTime[bmatrix[,i]],
    #       censor    = .Object@model@censor[bmatrix[,i]],
    #       trts      = .Object@model@trts, 
    #       timePoint = .Object@model@timePoint)
    #  }
    #  else
    #  if (model.type == "stmodelCI") {
    #    model.i <- stepp.CI(coltrt    = .Object@model@coltrt[bmatrix[,i]],
    #       coltime   = .Object@model@coltimeime[bmatrix[,i]],
    #       coltype   = .Object@model@coltype[bmatrix[,i]],
    #       trts      = .Object@model@trts, 
    #       timePoint = .Object@model@timePoint)
    #  }

      MM.new <- NULL
      if (!is.null(.Object@MM)) MM.new <- .Object@MM[subsample,]
      
        model.new  <- stepp.GLM(coltrt  = .Object@coltrt[subsample], 
              trts    = .Object@trts, 
              colY    = .Object@colY[subsample], 
              MM      = MM.new,
              glm     = .Object@glm)
      return(model.new)
  }
)

print.estimate.GLM <- function(x, family, trts){
  model <- x@model
  sp  <- x@subpop

  overall_lbl <- -1
  if (x@subpop@win@type == "tail-oriented") {
    if (length(x@subpop@win@r1) == 1) {
      overall_lbl <- 1
    } else if (length(x@subpop@win@r2) == 1) {
      overall_lbl <- length(x@subpop@win@r1) + 1
    }
  }

  for (j in 1:x@effect$ntrts) {
    # cat("\n")
    #
    # print out the glm ressults
    #
    est <- x@effect$TrtEff[[j]]

    if (family == "gaussian") {
      # cat("\n")
      if (x@effect$ntrts == 1) {
        write("Effect estimates", file = "")
      } else {
        write(paste("Effect estimates for treatment group", model@trts[j]), file = "")
      }
      temp <- matrix(c(1:sp@nsubpop, round(est$sObs, digits = 4), round(est$sSE, digits = 4)), ncol = 3)
      write("     Subpopulation     Effect       Std. Err.", file = "")
      for (i in 1:sp@nsubpop) {
        if (sp@win@type == "tail-oriented" & i == overall_lbl) {
                  write(paste(format(temp[i, 1], width = 12),
              format(temp[i, 2], width = 19, nsmall = 4),
              format(temp[i, 3], width = 15, nsmall = 4), "(entire cohort)"),
              file = "")
        } else {
                  write(paste(format(temp[i, 1], width = 12),
              format(temp[i, 2], width = 19, nsmall = 4),
              format(temp[i, 3], width = 15, nsmall = 4)),
              file = "")
        }
      }
      if (sp@win@type != "tail-oriented") {
        write(paste("        Overall",
          format(round(est$oObs, digits = 4), nsmall = 4, width = 16),
          format(round(est$oSE,  digits = 4), nsmall = 4, width = 15)),
          file = "")
      }
      cat("\n")
    } else if (family == "binomial") {
      # cat("\n")
      if (x@effect$ntrts == 1) {
        write("Risk estimates", file = "")
      } else {
        write(paste("Risk estimates for treatment group", model@trts[j]), file = "")
      }
      temp <- matrix(c(1:sp@nsubpop, round(est$sObs, digits = 4), round(est$sSE, digits = 4)), ncol = 3)
      write("     Subpopulation         Risk          Std. Err.", file = "")
      for (i in 1:sp@nsubpop) {
        if (sp@win@type == "tail-oriented" & i == overall_lbl) {
                  write(paste(format(temp[i, 1], width = 12),
              format(temp[i, 2], width = 19, nsmall = 4),
              format(temp[i, 3], width = 15, nsmall = 4), "(entire cohort)"),
              file = "")
        } else {
                  write(paste(format(temp[i, 1], width = 12),
              format(temp[i, 2], width = 19, nsmall = 4),
              format(temp[i, 3], width = 15, nsmall = 4)),
              file = "")
        }
      }
      if (sp@win@type != "tail-oriented") {
        write(paste("        Overall", 
          format(round(est$oObs, digits = 4), nsmall = 4, width = 16),
          format(round(est$oSE,  digits = 4), nsmall = 4, width = 15)),
          file = "")
      }
      cat("\n")
    } else if (family == "poisson") {
      # cat("\n")
      if (x@effect$ntrts == 1) {
        write("    Effect estimates", file = "")
      } else {
        write(paste("    Effect estimates for treatment group", model@trts[j]), file = "")
      }
      temp <- matrix(c(1:sp@nsubpop, round(est$sObs, digits = 4), round(est$sSE, digits = 4)), ncol = 3)
      write("     Subpopulation        Effect     Std. Err.", file = "")
      for (i in 1:sp@nsubpop) {
        if (sp@win@type == "tail-oriented" & i == overall_lbl) {
          write(paste(format(temp[i, 1], width = 12),
            format(temp[i, 2], width = 19, nsmall = 4),
            format(temp[i, 3], width = 15, nsmall = 4), "(entire cohort)"),
            file = "")
        } else {
          write(paste(format(temp[i, 1], width = 12), 
            format(temp[i, 2], width = 19, nsmall = 4), 
            format(temp[i, 3], width = 15, nsmall = 4)), 
            file = "")
        }
      }
      if (sp@win@type != "tail-oriented") {
        write(paste("        Overall", 
          format(round(est$oObs, digits = 4), nsmall = 4, width = 16),
          format(round(est$oSE,  digits = 4), nsmall = 4, width = 15)), 
          file = "")
      }
      cat("\n")
    }
  }

  if (x@effect$ntrts > 1) {
    # cat("\n")
    write("Effect differences and ratio estimates", file = "")
    est1 <- x@effect$TrtEff[[1]]

    for (j in 2:x@effect$ntrts) {
      cat("\n")
      write(paste("trt", trts[1], "vs. trt", trts[j]), file = "")

      est   <- x@effect$TrtEff[[j]]
      ratio <- x@effect$Ratios[[j-1]]

      if (family == "gaussian") {
        cat("\n")

        write(paste("Effect differences"), file = "")
        temp <- matrix(c(1:sp@nsubpop,
                       round(est1$sObs - est$sObs, digits = 4),
                       round(sqrt(est1$sSE^2 + est$sSE^2), digits = 4)), ncol = 3)
        #write("                          Effect", file = "")
        write("     Subpopulation      Difference      Std. Err.", file = "")
        for (i in 1:sp@nsubpop) {
          if (sp@win@type == "tail-oriented" & i == overall_lbl){
            write(paste(format(temp[i, 1], width = 12),
              format(temp[i, 2], width = 19, nsmall = 4),
              format(temp[i, 3], width = 15, nsmall = 4), "(entire cohort)"),
              file = "")
          } else {
              write(paste(format(temp[i, 1], width = 12),
            format(temp[i, 2], width = 19, nsmall = 4),
            format(temp[i, 3], width = 15, nsmall = 4)), 
            file = "")
          }
        }
        if (sp@win@type != "tail-oriented") {
          write(paste("        Overall",
            format(round(est1$oObs - est$oObs, digits = 4), nsmall = 4, width = 16), 
            format(round(sqrt(est1$oSE^2 + est$oSE^2), digits = 4), nsmall = 4, width = 15)),
            file = "")
        }
        cat("\n")

        write(paste("Effect ratios"), file = "")
        temp <- matrix(c(1:sp@nsubpop,
                       round(ratio$logHR, digits = 4),
                       round(ratio$logHRSE, digits = 4)), ncol = 3)
        write("                          Effect", file = "")
        write("     Subpopulation      Effect Ratio  Std. Err.  ", file = "")
        for (i in 1:sp@nsubpop) {
          if (sp@win@type == "tail-oriented" & i == overall_lbl) {
            write(paste(format(temp[i, 1], width = 12, ),
              format(round(temp[i, 2], digits = 4), width = 19, nsmall = 4),
              format(round(temp[i, 3], digits = 4), width = 15, nsmall = 4),
              "(entire cohort)"),
              file = "")
          } else {
            write(paste(format(temp[i, 1], width = 12, ),
              format(round(temp[i, 2], digits = 4), width = 19, nsmall = 4),
              format(round(temp[i, 3], digits = 4), width = 15, nsmall = 4)),
              file = "")
          }
        }
        if (sp@win@type != "tail-oriented") {
          write(paste("        Overall",
            format(round(ratio$ologHR,      digits = 4),  nsmall = 4, width = 16), 
            format(round(ratio$ologHRSE,    digits = 4),  nsmall = 4, width = 15)),
            file = "")
        }
        cat("\n")
      } else if (family == "binomial") {
        cat("\n")
        write(paste("Risk differences"), file = "")
        temp <- matrix(c(1:sp@nsubpop,
                       round(est1$sObs - est$sObs, digits = 4), 
                       round(sqrt(est1$sSE^2 + est$sSE^2), digits = 4)), ncol = 3)
        write("                           Risk", file = "")
        write("     Subpopulation      Difference       Std. Err.", 
          file = "")
        for (i in 1:sp@nsubpop) {
          if (sp@win@type == "tail-oriented" & i == overall_lbl) {
            write(paste(format(temp[i, 1], width = 12),
              format(temp[i, 2], width = 19, nsmall = 4),
              format(temp[i, 3], width = 15, nsmall = 4), "(entire cohort)"),
              file = "")
          } else {
            write(paste(format(temp[i, 1], width = 12),
              format(temp[i, 2], width = 19, nsmall = 4),
              format(temp[i, 3], width = 15, nsmall = 4)),
              file = "")
          }
        }
        if (sp@win@type != "tail-oriented") {
          write(paste("        Overall",
            format(round(est1$oObs - est$oObs, digits = 4), nsmall = 4, width = 16), 
            format(round(sqrt(est1$oSE^2 + est$oSE^2), digits = 4), nsmall = 4, width = 15)),
            file = "")
        }
        cat("\n")

        write(paste("Odds ratios"), file = "")
        temp <- matrix(c(1:sp@nsubpop,
                       round(ratio$logHR,      digits = 4),
                       round(ratio$logHRSE,    digits = 4),
                       round(exp(ratio$logHR), digits = 4)),
                      ncol = 4)
        write("                        log Odds", file = "")
        write("     Subpopulation        Ratio          Std. Err.          Odds Ratio", file = "")
        for (i in 1:sp@nsubpop) {
          if (sp@win@type == "tail-oriented" & i == overall_lbl) {
            write(paste(format(temp[i, 1], width = 12), 
              format(round(temp[i, 2], digits = 4), width = 19, nsmall = 4), 
              format(round(temp[i, 3], digits = 4), width = 15, nsmall = 4),
              format(round(temp[i, 4], digits = 4), width = 19, nsmall = 4), "(entire cohort)"),
              file = "")
          } else {
            write(paste(format(temp[i, 1], width = 12), 
              format(round(temp[i, 2], digits = 4), width = 19, nsmall = 4), 
              format(round(temp[i, 3], digits = 4), width = 15, nsmall = 4),
              format(round(temp[i, 4], digits = 4), width = 19, nsmall = 4)), 
              file = "")
          }
        }
        if (sp@win@type != "tail-oriented") {
          write(paste("        Overall",
            format(round(ratio$ologHR,     digits = 4), nsmall = 4, width = 16), 
            format(round(ratio$ologHRSE,   digits = 4), nsmall = 4, width = 15),
            format(round(exp(ratio$ologHR),digits = 4), nsmall = 4, width = 19)), 
            file = "")
        }

        cat("\n")
      } else if (family == "poisson") {
        cat("\n")
        write(paste("Effect differences"), file = "")
        temp <- matrix(c(1:sp@nsubpop, 
                       round(est1$sObs - est$sObs, digits = 4), 
                       round(sqrt(est1$sSE^2 + est$sSE^2), digits = 4)), ncol = 3)
        write("                          Effect ", file = "")
        write("     Subpopulation      Difference      Std. Err.", 
          file = "")
        for (i in 1:sp@nsubpop) {
          if (sp@win@type == "tail-oriented" & i == overall_lbl) {
            write(paste(format(temp[i, 1], width = 12), 
              format(temp[i, 2], width = 19, nsmall = 4), 
              format(temp[i, 3], width = 15, nsmall = 4), "(entire cohort)"),
              file = "")
          } else {
            write(paste(format(temp[i, 1], width = 12), 
              format(temp[i, 2], width = 19, nsmall = 4), 
              format(temp[i, 3], width = 15, nsmall = 4)), 
              file = "")
          }
        }
        if (sp@win@type != "tail-oriented") {
          write(paste("        Overall",
            format(round(est1$oObs - est$oObs, digits = 4), nsmall = 4, width = 16),
            format(round(sqrt(est1$oSE^2 + est$oSE^2), digits = 4), nsmall = 4, width = 15)),
            file = "")
        }
        cat("\n")

        write(paste("Relative    Effect "), file = "")
        temp <- matrix(c(1:sp@nsubpop,
                       round(ratio$logHR,      digits = 4),
                       round(ratio$logHRSE,    digits = 4),
                       round(exp(ratio$logHR), digits = 4)),
                       ncol = 4)
        write("                           log", file = "")
        write("     Subpopulation        Effect         Std. Err.        Relative Effect", 
          file = "")
        for (i in 1:sp@nsubpop) {
          if (sp@win@type == "tail-oriented" & i == overall_lbl){
            write(paste(format(temp[i, 1], width = 12), 
              format(temp[i, 2], width = 19, nsmall = 4), 
              format(temp[i, 3], width = 15, nsmall = 4),
              format(temp[i, 4], width = 19, nsmall = 4), "(entire cohort)"),
              file = "")
          } else {
            write(paste(format(temp[i, 1], width = 12), 
              format(temp[i, 2], width = 19, nsmall = 4), 
              format(temp[i, 3], width = 15, nsmall = 4),
              format(temp[i, 4], width = 19, nsmall = 4)),
              file = "")
          }
        }
        if (sp@win@type != "tail-oriented") {
          write(paste("        Overall", 
            format(round(ratio$ologHR,      digits = 4), nsmall = 4, width = 16), 
            format(round(ratio$ologHRSE,    digits = 4), nsmall = 4, width = 15),
            format(round(exp(ratio$ologHR), digits = 4), nsmall = 4, width = 19)), 
            file = "")
        }
        cat("\n")
      }
    }
  }
}

print.cov.GLM <- function(stobj, trts) {
  if (!is.null(stobj@testresults)) {
    for (j in 1:(stobj@testresults$ntrts - 1)) {
      ns <- stobj@subpop@nsubpop
      if (stobj@subpop@win@type == "tail-oriented") ns <- ns - 1
      # cat("\n")
      write(paste("The covariance matrix of the effect differences estimates for the",
        ns, "subpopulations is:"), file = "")
      write(paste("trt ", trts[1], "vs. trt", trts[j + 1]), file = "")
      print(stobj@testresults$Res[[j]]$sigma)

      cat("\n")
      write(paste("The covariance matrix of the effect ratios for the", 
        ns, "subpopulations is:"), file = "")
      print(stobj@testresults$Res[[j]]$HRsigma)
      cat("\n")
    }
  }
}

print.stat.GLM <- function(stobj, trts) {
  if (!is.null(stobj@testresults)) {
    for (j in 1:(stobj@testresults$ntrts-1)) {

      t <- stobj@testresults$Res[[j]]
      # cat("\n")
      write(paste("Supremum test results"), file = "")
      write(paste("trt", trts[1], "vs. trt", trts[j + 1]), file = "")

      write(paste("Interaction p-value based on effect estimate differences:", t$pvalue), file = "")
      cat("\n")
      write(paste("Chi-square interaction p-value based on effect estimate differences:", t$chi2pvalue), file = "")

      cat("\n")
    }
  }
}

setMethod("print",
  signature = "stmodelGLM",
  definition = function(x, stobj, estimate = TRUE, cov = TRUE, test = TRUE, ...){
    ntrts <- length(x@trts)

    #
    #  1. estimates
    #
    if (estimate) {
      print.estimate.GLM(stobj, x@glm, x@trts)
    }

    if (!is.null(stobj@testresults)) {
      #
      #   2. covariance matrices
      #
      if (cov & ntrts > 1) {
        print.cov.GLM(stobj, x@trts)
      }
  
      #
      #   3. Supremum test and Chi-square test results
      #
      if (test & ntrts > 1) {
        print.stat.GLM(stobj, x@trts)
      }
    }
 }
)

# constructor function for stmodelGLM
stepp.GLM <- function(coltrt, trts, colY, MM = NULL, glm, link){
  model <- new("stmodelGLM", coltrt = coltrt, trts = trts, colY = colY, MM = MM, glm = glm, link = link)
  return(model)
}
