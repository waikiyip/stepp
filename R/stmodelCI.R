#################################################################
#
# stmodelCI.R
#
#######################
# stepp model: CI     #
#######################
#
cuminc.HR <- function (ftime, fstatus, group, strata, rho = 0, cencode = 0, 
    subset, na.action = na.omit) {
#
    d <- data.frame(time = ftime, cause = fstatus, 
        group = as.factor(
            if (missing(group)) rep(1, length(ftime)) else group),
        strata = as.factor(if (missing(strata)) rep(1, length(ftime)) else strata))
    if (!missing(subset)) d <- d[subset, ]
    tmp <- nrow(d)
    d <- na.action(d)
    if (nrow(d) != tmp) 
      cat(format(tmp - nrow(d)), "cases omitted due to missing values\n")
    no <- nrow(d)
    cg <- "  "
    nst <- length(levels(d$strata))
    d <- d[order(d$time), ]
    ugg <- table(d$group)
    d$group <- factor(d$group, names(ugg)[ugg > 0])
    ugg <- levels(d$group)
    censind <- ifelse(d$cause == cencode, 0, 1)
    uc <- table(d$cause[censind == 1])
    if (is.factor(d$cause)) 
        uclab <- names(uc)[uc > 0]
    else uclab <- as.numeric(names(uc)[uc > 0])
    nc <- length(uclab)
    ng <- length(ugg)
    if (ng > 1) {
        ng1 <- ng - 1
        ng2 <- ng * ng1/2
        v <- matrix(0, nrow = ng1, ncol = ng1)
        storage.mode(v) <- "double"
        vt <- double(ng2)
        s <- double(ng1)
    }
    pf <- vector("list", ng * nc)
    stat <- double(nc)
    l <- 0
    for (ii in 1:nc) {
        causeind <- ifelse(d$cause == uclab[ii], 1, 0)
        for (jj in 1:length(ugg)) {
            cg <- c(cg, paste(ugg[jj], uclab[ii]))
            l <- l + 1
            cgind <- d$group == ugg[jj]
            ncg <- length(cgind[cgind])
            n2 <- length(unique(d$time[cgind & causeind == 1]))
            n2 <- 2 * n2 + 2
            tmp <- double(n2)
            z <- .Fortran("cinc", as.double(d$time[cgind]), as.integer(censind[cgind]), 
                as.integer(causeind[cgind]), as.integer(ncg), 
                x = tmp, f = tmp, v = tmp)
#                x = tmp, f = tmp, v = tmp, PACKAGE = "cmprsk")
            pf[[l]] <- list(time = z$x, est = z$f, var = z$v)
        }
        if (ng > 1) {
            causeind <- 2 * censind - causeind
            z2 <- .Fortran("crstm", as.double(d$time), as.integer(causeind), 
                as.integer(d$group), as.integer(d$strata), as.integer(no), 
                as.double(rho), as.integer(nst), as.integer(ng), 
                s, v, as.double(d$time), as.integer(causeind), 
                as.integer(d$group), vt, s, vt, double((4 + 3 * 
                  ng) * ng), integer(4 * ng))
#                  ng) * ng), integer(4 * ng), PACKAGE = "cmprsk")
            stat[ii] <- -1
            a <- qr(z2[[10]])
            if (a$rank == ncol(a$qr)) {
                b <- diag(dim(a$qr)[1])
                stat[ii] <- z2[[9]] %*% qr.coef(a, b) %*% z2[[9]]
            }
#
            if (ii == 1) {
               ome <- (-1)*z2[[9]]
               omevar <- z2[[10]]
            }
        }
    }
    names(pf) <- cg[2:length(cg)]
    if (ng > 1) {
        names(stat) <- uclab
        stat <- list(Tests = cbind(stat = stat, pv = 1 - pchisq(stat, 
            ng - 1), df = rep(ng - 1, length(stat))))
        pf <- c(pf, stat)
#
        omeres <- list(ome=ome,omevar=omevar)
        pf <- c(pf,omeres)
    }
    pf
}

setClass("stmodelCI",
  representation(
    coltrt    = "numeric",  # treatment
    coltime   = "numeric",  # time to event
    coltype   = "numeric",  # competing risk type
    trts      = "numeric",  # trt encoding
    timePoint = "numeric"   # evaluated time
  )
)

setMethod("initialize", "stmodelCI",
  function(.Object, coltrt, coltime, coltype, trts, timePoint, ...) {
    if (missing(coltime)) {
      coltime <- numeric()
    }
    if (missing(coltype)) {
      coltype <- numeric()
    }
    if (missing(timePoint)) {
      timePoint <- 1
    }
    if (missing(coltrt)) {
      coltrt <- rep(1, length(coltime))
    }
    if (missing(trts)) {
      trts <- 1
    }
    .Object@coltrt <- coltrt
    .Object@coltime <- coltime
    .Object@coltype <- coltype
    .Object@trts <- trts
    .Object@timePoint <- timePoint
    if (!validObject(.Object)) stop("")
    callNextMethod(.Object, ...)
  }
)

setValidity("stmodelCI", 
  function(object) {
    status <- TRUE

    # add conditions to verify

    return(status)
  }
)

setMethod("estimate",
  signature = "stmodelCI",
  definition = function(.Object, sp, ...) {
    nsubpop   <- sp@nsubpop
    subpop    <- sp@subpop
    coltrt    <- .Object@coltrt
    survTime  <- .Object@coltime
    trts      <- .Object@trts
    timePoint <- .Object@timePoint
    type      <- .Object@coltype

    ntrts     <- length(trts)

    # set up the return structure   
    Obs <- list(sObs = rep(0, nsubpop),
      sSE  = rep(0, nsubpop),
      oObs = 0,
      oSE  = 0)
    TrtEff <- array(list(Obs), ntrts)
    
    HRs <- list(skmw = 0,
      logHR = rep(0, nsubpop),
      logHRSE = rep(0, nsubpop),
      ologHR  = 0,
      ologHRSE = 0,
      logHRw = 0)
    Ratios <- array(list(HRs), ntrts - 1)

    # estimate the abs and rel trt effects comparing trt 1 with trt j
    if (ntrts > 1) {
      for (j in 2:ntrts) {
        txassign <- rep(-1, length(coltrt))
        txassign[which(coltrt == trts[1])] <- 1
        txassign[which(coltrt == trts[j])] <- 0

        sel <- (txassign == 1 | txassign == 0)

        # for each subgroup i
        for (i in 1:nsubpop) {
          result <- cuminc.HR(survTime[subpop[,i]==1& !(txassign==-1)],
          type[subpop[,i]==1& !(txassign==-1)],txassign[subpop[,i]==1 & !(txassign==-1)])
          if (max(result$"1 1"$time) >= timePoint) {
            index           <- sum(result$"1 1"$time <= timePoint)
            TrtEff[[1]]$sObs[i]  <- result$"1 1"$est[index]
            TrtEff[[1]]$sSE[i]   <- sqrt(result$"1 1"$var[index])
          }
          if (max(result$"0 1"$time) >= timePoint) {
            index           <- sum(result$"0 1"$time <= timePoint)
            TrtEff[[j]]$sObs[i]  <- result$"0 1"$est[index]
            TrtEff[[j]]$sSE[i]   <- sqrt(result$"0 1"$var[index])
          }
          Ratios[[j - 1]]$logHR[i]   <- result$ome/result$omevar
          Ratios[[j - 1]]$logHRSE[i] <- sqrt(1/result$omevar)
        }
        Ratios[[j - 1]]$logHRw <- sum(Ratios[[j - 1]]$logHR/Ratios[[j - 1]]$logHRSE)
        Ratios[[j - 1]]$skmw   <- sum((TrtEff[[1]]$sObs-TrtEff[[j]]$sObs)/sqrt(TrtEff[[1]]$sSE^2 + TrtEff[[j]]$sSE^2))

        # estimate the overall effect
        result <- cuminc.HR(survTime[sel],type[sel],txassign[sel])
        if (max(result$"1 1"$time) >= timePoint) {
          index <- length(result$"1 1"$time[result$"1 1"$time <= timePoint])
          TrtEff[[1]]$oObs <- result$"1 1"$est[index]
          TrtEff[[1]]$oSE  <- sqrt(result$"1 1"$var[index])
        }
        if (max(result$"0 1"$time) >= timePoint) {
          index <- sum(result$"0 1"$time <= timePoint)
          TrtEff[[j]]$oObs <- result$"0 1"$est[index]
          TrtEff[[j]]$oSE  <- sqrt(result$"0 1"$var[index])
        }
        Ratios[[j - 1]]$oLogHR   <- result$ome/result$omevar
        Ratios[[j - 1]]$oLogHRSE <- sqrt(1/result$omevar)
        if (sum(is.na(TrtEff[[1]]$sObs)) != 0 | sum(is.na(TrtEff[[j]]$sObs)) != 0 | 
              is.na(TrtEff[[1]]$oObs) != FALSE | is.na(TrtEff[[j]]$oObs) != FALSE) {
          cat("\n")
          print(paste("Unable to estimate survival time at ", timePoint, 
          " time-unit(s) because there are too few events within one or more subpopulation(s)."))
          print(paste("The problem may be avoided by constructing larger subpopulations and/or by selecting a different timepoint for estimation."))
          stop()
        }
      }
    }

    if (ntrts == 1) {
      txassign <- rep(-1, length(coltrt))
      txassign[which(coltrt == trts[1])] <- 1

      sel <- (txassign == 1 | txassign == 0)

      # for each subgroup i
      for (i in 1:nsubpop) {
        result <- cuminc.HR(survTime[subpop[,i]==1& !(txassign==-1)],
        type[subpop[,i]==1& !(txassign==-1)],txassign[subpop[,i]==1 & !(txassign==-1)])
        if (max(result$"1 1"$time) >= timePoint) {
          index           <- sum(result$"1 1"$time <= timePoint)
          TrtEff[[1]]$sObs[i]  <- result$"1 1"$est[index]
          TrtEff[[1]]$sSE[i]   <- sqrt(result$"1 1"$var[index])
        }
      }

      # estimate the overall effect
      result <- cuminc.HR(survTime[sel],type[sel],txassign[sel])
      if (max(result$"1 1"$time) >= timePoint) {
        index <- length(result$"1 1"$time[result$"1 1"$time <= timePoint])
        TrtEff[[1]]$oObs <- result$"1 1"$est[index]
        TrtEff[[1]]$oSE  <- sqrt(result$"1 1"$var[index])
      }
      if (sum(is.na(TrtEff[[1]]$sObs)) != 0 | is.na(TrtEff[[1]]$oObs) != FALSE) {
        cat("\n")
        print(paste("Unable to estimate survival time at ", timePoint, 
        " time-unit(s) because there are too few events within one or more subpopulation(s)."))
        print(paste("The problem may be avoided by constructing larger subpopulations and/or by selecting a different timepoint for estimation."))
        stop()
      }
    }
    #
    # create the estimates (Cumulative Incidence) object - to be saved
    #
    estimate <- list(
      model   = "CIe",
      ntrts   = ntrts,
      TrtEff  = TrtEff,
      Ratios  = Ratios)
    return(estimate)
  }
)

setMethod("test",
  signature = "stmodelCI",
  definition = function(.Object, nperm, sp, effect, showstatus = TRUE, ...) {
    test <- NULL

    # return immediately if nperm is 0
    if (nperm > 0 & length(.Object@trts) > 1) {
      win   <- sp@win
      nsubpop <- sp@nsubpop
      osubpop <- sp@subpop
      coltrt  <- .Object@coltrt
      survTime  <- .Object@coltime
      trts  <- .Object@trts
      timePoint <- .Object@timePoint
      type      <- .Object@coltype
      ntrts <- length(trts)

      # set up return structure
      tsigma    <- matrix(rep(0, (nperm * nsubpop)), ncol = nsubpop)    
      PermTest  <- list(sigma = tsigma,
        HRsigma    = tsigma,
        pvalue     = 0,
        chi2pvalue = 0,
        HRpvalue   = 0)
      Res <- array(list(PermTest), ntrts - 1)

      for (j in 2:ntrts) {
      #
      # do the permutations
      # compare trt1 with the rest
      #

      txassign <- rep(-1, length(coltrt))
      txassign[which(coltrt == trts[1])] <- 1
      txassign[which(coltrt == trts[j])] <- 0
      #
      #   do the permutations
      #

      # set up the intial values
      pvalue      <- NA
      differences <- matrix(rep(0, (nperm * nsubpop)), ncol = nsubpop)
      #diffha      <- matrix(rep(0, (nperm * nsubpop)), ncol = nsubpop)
      logHRs      <- matrix(rep(0, (nperm * nsubpop)), ncol = nsubpop)
      #logHRha     <- matrix(rep(0, (nperm * nsubpop)), ncol = nsubpop)
      tPerm       <- rep(0, nperm)
      no          <- 0
      p           <- 0
      terminate   <- 0
      Ntemp       <- nrow(osubpop)
      IndexSet1   <- (1:Ntemp)[txassign == 1]
      IndexSet2   <- (1:Ntemp)[txassign == 0]

      if (showstatus) {
        title <- paste("Computing the p-value with", nperm)
        title <- paste(title, "number of permutations comparing trt", trts[j], "with trt", trts[1], "\n")
        cat(title)
        pb <- txtProgressBar(min = 0, max = nperm - 1, style = 3)
      }

      while (no < nperm) {
        ## PERMUTE THE VECTOR OF SUBPOPULATIONS WITHIN TREATMENTS ##
        if (showstatus) setTxtProgressBar(pb, no)

        sel    <- (txassign == 1 | txassign == 0)

        Subpop <- as.matrix(osubpop)
        permuteSubpop <- matrix(0, nrow = Ntemp, ncol = nsubpop)
        permuteSubpop[txassign == 1] <- Subpop[sample(IndexSet1),]
        permuteSubpop[txassign == 0] <- Subpop[sample(IndexSet2),]
        subpop      <- permuteSubpop
        sCI1        <- rep(NA,nsubpop)
        sCI2        <- rep(NA,nsubpop)
        sCISE1      <- rep(NA,nsubpop)
        sCISE2      <- rep(NA,nsubpop)
        overallsCI1 <- NA
        overallsCI2 <- NA
        slogHRs     <- rep(0, nsubpop)
        slogHRSE    <- rep(0, nsubpop)

        for (i in 1:nsubpop) {
          result <- cuminc.HR(survTime[subpop[,i]==1 & !(txassign==-1)],
            type[subpop[,i]==1 & !(txassign==-1)], txassign[subpop[,i]==1 & !(txassign==-1)])
          if (max(result$"1 1"$time) >= timePoint) {
            index     <- sum(result$"1 1"$time <= timePoint)
            sCI1[i]   <- result$"1 1"$est[index]
            sCISE1[i] <- sqrt(result$"1 1"$var[index])
          }
          if (max(result$"0 1"$time) >= timePoint) {
            index     <- sum(result$"0 1"$time <= timePoint)
            sCI2[i]   <- result$"0 1"$est[index]
            sCISE2[i] <- sqrt(result$"0 1"$var[index])
          }
          slogHRs[i] <- result$ome/result$omevar
          slogHRSE[i] <- sqrt(1/result$omevar)
        }
        slogHRw <- sum(slogHRs/slogHRSE)
        #skmwha  <- sum((sCI1 - sCI2)/sqrt(sCISE1^2 + sCISE2^2))

        result <- cuminc.HR(survTime[sel], type[sel], txassign[sel])
        if (max(result$"1 1"$time) >= timePoint) {
          index <- length(result$"1 1"$time[result$"1 1"$time <= timePoint])
          overallsCI1 <- result$"1 1"$est[index]
        }
        if (max(result$"0 1"$time) >= timePoint) {
          index <- sum(result$"0 1"$time <= timePoint)
          overallsCI2 <- result$"0 1"$est[index]
        }
        overallSlogHR <- result$ome/result$omevar

        if (sum(is.na(sCI1)) == 0 & sum(is.na(sCI2)) == 0 & is.na(overallsCI1) == FALSE &
          is.na(overallsCI2) == FALSE) {
          no <- no + 1
          p <- p + 1
          for (s in 1:nsubpop) {
            differences[p, s] <- (sCI1[s] - sCI2[s]) - (overallsCI1 - overallsCI2)
            # diffha[p, s] <- (sCI1[s] - sCI2[s]) - skmwha
            logHRs[p, s] <- slogHRs[s] - overallSlogHR
            # logHRha[p, s] <- slogHRs[s] - slogHRw
          }
        }
        terminate <- terminate + 1
        if (terminate >= nperm + 10000) {
          print(paste("After permuting ", nperm, "plus 10000, or ", 
            nperm + 10000, " times, the program is unable to generate the permutation distribution based on ", 
            nperm, "permutations of the data."))
          print("Consider creating larger subpopulations or selecting a different timepoint for estimation.")
          stop()
        }
      }
      if (showstatus) close(pb)

      # generating the sigmas and p-values
      sObs   <- effect$TrtEff[[1]]$sObs-effect$TrtEff[[j]]$sObs
      oObs   <- effect$TrtEff[[1]]$oObs-effect$TrtEff[[j]]$oObs
      logHR  <- effect$Ratios[[j - 1]]$logHR
      ologHR <- effect$Ratios[[j - 1]]$ologHR

      mname  <- paste0("SP", seq(1, dim(differences)[2]), "-Overall")
      if (win@type == "tail-oriented") {
        # remove the trivial case of the full cohort
        rm <- length(win@r1)+1
        differences <- differences[,-rm]
        mname <- mname[-rm]
        sObs <- sObs[-rm]
        oObs <- oObs[-rm]
        logHRs <- logHRs[,-rm]
        logHR <- logHR[-rm]
        ologHR <- ologHR[-rm]
      }

      sigmas                <- ssigma(differences)
      rownames(sigmas$sigma)<- mname
      colnames(sigmas$sigma)<- mname
      Res[[j - 1]]$sigma      <- sigmas$sigma
      Res[[j - 1]]$chi2pvalue <- ppv (differences, sigmas$sigmainv, sObs, oObs, nperm)
      Res[[j - 1]]$pvalue     <- ppv2(differences,                  sObs, oObs, nperm)

      sigmas                <- ssigma(logHRs)
      rownames(sigmas$sigma)<- mname
      colnames(sigmas$sigma)<- mname
      Res[[j - 1]]$HRsigma    <- sigmas$sigma
      Res[[j - 1]]$HRpvalue   <- ppv2(logHRs,                       logHR, ologHR, nperm) 

      }
      test = list(
        model  = "CIt",
        ntrts  = ntrts,
        Res    = Res
      )
    }
    return(test)
  }
)

#
# printing support functions for cumulative incidence model
#
print.estimate.CI <- function(x, timePoint, trts) {
  for (j in 1:x@effect$ntrts) {
    cat("\n")
    if (x@effect$ntrts == 1) {
      write(paste("Cumulative incidence estimates ",
        "at time point", timePoint), file="")
    } else {
      write(paste("Cumulative incidence estimates for treatment group", trts[j],
        "at time point", timePoint), file="")
    }

    overall_lbl <- -1
    if (x@subpop@win@type == "tail-oriented") {
      if (length(x@subpop@win@r1) == 1) {
        overall_lbl <- 1
      } else if (length(x@subpop@win@r2) == 1) {
        overall_lbl <- length(x@subpop@win@r1) + 1
      }
    }

    temp <- matrix(c(1:x@subpop@nsubpop, round(x@effect$TrtEff[[j]]$sObs, digits=4),
      round(x@effect$TrtEff[[j]]$sSE, digits=4)), ncol=3)
    write("                        Cumulative",file="")
    write("     Subpopulation      Incidence        Std. Err.",file="")
    for (i in 1:x@subpop@nsubpop) {
      if (x@subpop@win@type == "tail-oriented" & i == overall_lbl) {
        write(paste(format(temp[i,1], width=12), format(temp[i,2], width=19, nsmall=4),
          format(temp[i,3], width=15, nsmall=4), "(entire cohort)"), file="")
      } else {
        write(paste(format(temp[i,1], width=12), format(temp[i,2], width=19, nsmall=4),
          format(temp[i,3], width=15, nsmall=4)), file="")
      }
    }
    if (x@subpop@win@type != "tail-oriented") {
      write(paste("        Overall", 
        format(round(x@effect$TrtEff[[j]]$oObs, digits=4), nsmall=4, width=16),
        format(round(x@effect$TrtEff[[j]]$oSE,  digits=4), nsmall=4, width=15)),
        file="")
    }
    cat("\n")
  }

  if (x@effect$ntrts > 1) {
    cat("\n")
    write(paste("Cumulative incidence differences at time point",timePoint), file="")

    for (j in 2:x@effect$ntrts) {
      cat("\n")
      write(paste("trt", trts[1], "vs. trt", trts[j]), file = "")

      temp <- matrix(c(1:x@subpop@nsubpop,
        round(x@effect$TrtEff[[1]]$sObs - x@effect$TrtEff[[j]]$sObs, digits = 4),
        round(sqrt(x@effect$TrtEff[[1]]$sSE^2+x@effect$TrtEff[[j]]$sSE^2), digits = 4)), ncol = 3)
      write("                        Cumulative",file="")
      write("                        Incidence",file="")
      write("     Subpopulation      Difference       Std. Err.",file="")
      for (i in 1:x@subpop@nsubpop) {
        if (x@subpop@win@type == "tail-oriented" & i == overall_lbl) {
          write(paste(format(temp[i,1],width=12),format(temp[i,2],width=19,nsmall=4),
            format(temp[i,3],width=15,nsmall=4),"(entire cohort)"),file="")
        } else {
          write(paste(format(temp[i,1],width=12),format(temp[i,2],width=19,nsmall=4),
            format(temp[i,3],width=15,nsmall=4)),file="")
        }
      }
      if (x@subpop@win@type != "tail-oriented") {
        write(paste("        Overall",
          format(round(x@effect$TrtEff[[1]]$oObs-x@effect$TrtEff[[j]]$oObs, digits = 4), nsmall = 4, width = 16),
          format(round(sqrt(x@effect$TrtEff[[1]]$oSE^2+x@effect$TrtEff[[j]]$oSE^2), digits = 4), nsmall = 4, width = 15)), file="")
      }
      cat("\n")
      write("Hazard ratio estimates",file="")
      temp <- matrix(c(1:x@subpop@nsubpop,
        round(x@effect$Ratios[[j - 1]]$logHR,digits=6),
        round(x@effect$Ratios[[j - 1]]$logHRSE,digits=6),
        round(exp(x@effect$Ratios[[j - 1]]$logHR),digits=2)),ncol=4)
      write("     Subpopulation        Log HR       Std. Err.       Hazard Ratio",file="")
      for (i in 1:x@subpop@nsubpop) {
        if (x@subpop@win@type == "tail-oriented" & i == overall_lbl){
          write(paste(format(temp[i,1],width=12),format(temp[i,2],width=19,nsmall=6),
            format(temp[i,3],width=14,nsmall=6),format(temp[i,4],width=15,nsmall=2),"(entire cohort)"),file="")
        } else {
          write(paste(format(temp[i,1],width=12),format(temp[i,2],width=19,nsmall=6),
            format(temp[i,3],width=14,nsmall=6),format(temp[i,4],width=15,nsmall=2)),file="")
        }
      }
      if (x@subpop@win@type != "tail-oriented") {
        write(paste("        Overall",
          format(round(x@effect$Ratios[[j - 1]]$oLogHR,      digits=6), nsmall=6, width=16),
          format(round(x@effect$Ratios[[j - 1]]$oLogHRSE,    digits=6), nsmall=6, width=14),
          format(round(exp(x@effect$Ratios[[j - 1]]$oLogHR), digits=2), nsmall=2, width=15)),
          file="")
      }
      cat("\n")
    }
  }
}

print.cov.CI <- function(stobj, timePoint, trts) {
  if (!is.null(stobj@result)) {
    for (j in 1:(stobj@result$ntrts - 1)) {
      ns <- stobj@subpop@nsubpop
      if (stobj@subpop@win@type == "tail-oriented") ns <- ns-1
      cat("\n")
      write(paste("The covariance matrix of the cumulative incidence differences at", 
        timePoint, "time units for the", ns, "subpopulations is:"), file = "")
      write(paste("trt ", trts[1], "vs. trt", trts[j + 1]), file = "")
      print(stobj@result$Res[[j]]$sigma)

      cat("\n")
      write(paste("The covariance matrix of the log hazard ratios for the", ns, "subpopulations is:"), file = "")
      print(stobj@result$Res[[j]]$HRsigma)
      cat("\n")
    }
  }
}

print.stat.CI <- function(stobj, trts) {
  if (!is.null(stobj@result)) {
    for (j in 1:(stobj@result$ntrts - 1)) {
      t <- stobj@result$Res[[j]]
      cat("\n")
      write(paste("Supremum test results"), file = "")
      write(paste("trt", trts[1], "vs. trt", trts[j + 1]), file = "")
      write(paste("Interaction p-value based on cumulative incidence estimates:", t$pvalue), file = "")
      write(paste("Interaction p-value based on hazard ratio estimates:", t$HRpvalue), file = "")

      cat("\n")
      write(paste("Chi-square test results"), file = "")
      write(paste("Interaction p-value based on cumulative incidence estimates:", t$chi2pvalue), file = "")

      #cat("\n")
      #write(paste("Homogeneous association test results"), file = "")
      #write(paste("Interaction p-value based on cumulative incidence estimates:", 
      #      hapvalue), file = "")
      #write(paste("Interaction p-value based on hazard ratio estimates:", 
      #      haHRpvalue), file = "")

      cat("\n")
    }
  }
}

setMethod("print",
  signature = "stmodelCI",
  definition = function(x, stobj, estimate = TRUE, cov = TRUE, test = TRUE, ...) {
    ntrts <- length(x@trts)

    #
    # 1. estimates
    #
    if (estimate) {
      print.estimate.CI(stobj, x@timePoint, x@trts)
    }

    #
    # 2. covariance matrices
    #
    if (cov & ntrts > 1) {
      print.cov.CI(stobj, x@timePoint, x@trts)
    }
  
    #
    # 3. Supremum test and Chi-square test results
    #
    if (test & ntrts > 1) {
      print.stat.CI(stobj, x@trts)
    }
  }
)

# constructor function for stmodelKM
stepp.CI <- function(coltrt, coltime, coltype, trts, timePoint) {
  model <- new("stmodelCI", coltrt = coltrt, coltime = coltime, coltype = coltype,
    trts = trts, timePoint = timePoint)
  return(model)
}
