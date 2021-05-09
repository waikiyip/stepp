#################################################################
#
# stmodelKM.R
#
#######################
# stepp model: KM     #
#######################
#
kmest1 <- function(y,m,n,ndf,t,s,v,ntpt,tpt,nrr,ndd) {
  # initialization 
  f<-1.0
  kr<-n
  nrr[1]<-n
  ndd[1]<-0

  for (i in 2:ntpt){
    nrr[i]<-0
    ndd[i]<-0
  }

  ltp<-1
  var<-0.0
  l<-1
  t[1]<-0
  s[1]<-1
  v[1]<-0
  i<-1

  # main loop 
  while (i<= n) {
    k<-i+1
    k2<-0

    while (k2<=0) {
      if (k > n) k2<-1
      else {
  if (y[k] != y[i]) k2<-1 else k<-k+1
      }
  
    } # end while 

    k<-k-1
    nd<-0

    for (j in i:k) nd<-nd+m[j]

    while (ltp<=ntpt && y[i]>tpt[ltp+1]) {
      ltp<-ltp+1
      nrr[ltp]<-kr
    }  # end while 

    ndd[ltp]<-ndd[ltp]+nd
    if (nd>0) {
      t1<- nd / kr  
      f<-f*(1-t1)
      if (nd<kr) var<-var+t1/(kr-nd)
      t[l+1]<-y[i]
      s[l+1]<-s[l]
      v[l+1]<-v[l]
      l<-l+2
      t[l]<-y[i]
      s[l]<-f
      v[l]<-var*f*f
    } # end if 

    i<-k+1
    kr<-n-k     # kr<-n-k-1
    k<-i
  } # end while - main loop 

  l<-l+1
  t[l]<-y[n]
  s[l]<-s[l-1]
  v[l]<-v[l-1]

  # return a list of all the arguments
  kmtest1<-list(y=y,m=m,n=n,ndf=ndf,t=t,s=s,v=v,ntpt=ntpt,tpt=tpt,nrr=as.integer(nrr),ndd=as.integer(ndd))
  kmtest1
  
}


kmest <- function (time,status,group,tpt,pv=TRUE,pv.strat,pv.sub,rho=0,subset,na.action=na.omit) {
  d <- data.frame(time=time,status=status,
    group=as.factor(if (missing(group)) rep(1,length(time)) else group),
    pv.strat=as.factor(if (missing(pv.strat)) rep(1,length(time)) else pv.strat),
    pv.sub = as.factor(if (missing(pv.sub)) rep(TRUE, length(time)) else pv.sub))
  if (!missing(subset)) d <- d[subset,]
  tmp <- nrow(d)
  d <- na.action(d)
  if (nrow(d) != tmp) cat(format(tmp-nrow(d)),'cases omitted due to missing values\n')
  no <- nrow(d)
  if (any(d$status != 0 & d$status != 1)) stop("invalid status values")
  T1 <- table(d$group)
  subgl <- T1>0
  T8 <- d$group[d$status==1]
  if (length(T8) > 0) {
    T8 <- cbind(T1,table(d$group[d$status==1]))
  } else {
    T8 <- cbind(T1,rep(0,length(T8)))
  }
  T8 <- T8[subgl,,drop=FALSE]
  T1 <- names(T1)[subgl]
  time <- ifelse(d$time<=0,.00001,d$time)
  if (missing(tpt)) {
    tpt <- pretty(time)}
  else {
    ym <- round(max(time),2)
    tpt <- c(0,tpt,ym)
  }
  ntpt <- length(tpt)-1
  lev <- vector("character",ntpt)
  for (i in 1:ntpt) lev[i] <- paste(format(tpt[i]),format(tpt[i+1]),sep="-")
  nrd <- rep(0, ntpt)
  o <- order(time)
  Tl <- length(T1)
  z <- as.list(1:Tl)
  for(i in 1:Tl) {
    Ty <- (time[o])[d$group[o]==T1[i]]
    Tm <- (d$status[o])[d$group[o]==T1[i]]
    ndf <- length(unique(Ty[Tm==1]))
    t <- double(2*ndf+2)

    a <- kmest1(Ty, Tm, length(Ty), ndf, t,t,t,ntpt, tpt, nrd, nrd)

    tt <- paste(a[[11]],a[[10]],sep="/")
    names(tt) <- lev
    z[[i]] <- list(time=a[[5]],est=a[[6]],var=a[[7]],tint=tt,nnd=T8[i,])
    }
  names(z) <- T1
  if (pv & Tl > 1) {
    stop("internal error")
    # pv <- logrank(time=d$time,status=d$status,group=d$group,strata=d$pv.strat,rho=rho,subset=(d$pv.sub=='TRUE'))$pv
    # attr(z,"pv") <- pv
  }

  class(z) <- "kmest"
  z
}

#
tpest1 <- function(x,n,ind,tp5,ntp) {
  l <- ntp
  for (i in rev(1:ntp)){
    if (x[n] >= tp5[i]) break;
    ind[l]<-0
    l <- l - 1
  }

  if (l <=0) return (ind);

  if (x[n] == tp5[l]){
    ind[l] <- n
    l <- l - 1
  }

  # assuming unique values sorted in ascending order
  k <- n-1
  loop <- TRUE
  while(loop){
    if (l <= 0) return (ind);

   loop <- FALSE
    for (i in rev(1:k)) {
      if (x[k] <= tp5[l]){
        ind[l] <- k+1
        l <- l - 1
        loop <- TRUE
        break #out of the for loop
      } else 
      k <- k-1
    } # end for loop
  } # end while loop

  # error in the following loop corrected 9-28-04
  for (i in 1:l) ind[i] <- 0
  return (ind); 
}


tpest <- function(w,times) {
  if (!is.null(w$Tests)) w <- w[names(w) != 'Tests']
  ng <- length(w)
  times <- sort(unique(times))
  nt <- length(times)
  storage.mode(times) <- "double"
  storage.mode(nt) <- "integer"
  ind <- matrix(0,ncol=nt,nrow=ng)
  oute <- matrix(NA,ncol=nt,nrow=ng)
  outv <- oute
  storage.mode(ind) <- "integer"
  slct <- rep(TRUE,ng)
  for (i in 1:ng) {
    if (is.null((w[[i]])$est)) {
      slct[i] <- FALSE
    } else {
      z1 <- as.integer(tpest1(w[[i]][[i]], length(w[[i]][[i]]), ind[i,], times, nt))
      ind[i,] <- z1
      oute[i, ind[i,]>0] <- w[[i]][[2]][z1]
      if (length(w[[i]])>2) outv[i,ind[i,]>0] <- w[[i]][[3]][z1]
    }
  }
  
  dimnames(oute) <- list(names(w)[1:ng],as.character(times))
  dimnames(outv) <- dimnames(oute)
  list(est=oute[slct,,drop=FALSE],var=outv[slct,,drop=FALSE])
}

setClass("stmodelKM",
  representation(
    coltrt    = "numeric",  # treatment
    survTime  = "numeric",  # time to event
    censor    = "numeric",  # time to censor
    trts      = "numeric",  # trt encoding
    timePoint = "numeric"   # evaluated time
  )
)

setMethod("initialize", "stmodelKM",
  function(.Object, coltrt, survTime, censor, trts, timePoint, ...) {
    if (missing(survTime)) {
      survTime <- numeric()
    }
    if (missing(censor)) {
      censor <- numeric()
    }
    if (missing(timePoint)) {
      timePoint <- 1
    }
    if (missing(coltrt)) {
      coltrt <- rep(1, length(survTime))
    }
    if (missing(trts)) {
      trts <- 1
    }
    .Object@coltrt <- coltrt
    .Object@survTime <- survTime
    .Object@censor <- censor
    .Object@trts <- trts
    .Object@timePoint <- timePoint
    if (!validObject(.Object)) stop("")
    callNextMethod(.Object, ...)
  }
)

setValidity("stmodelKM", 
  function(object) {
    status <- TRUE

    # add conditions to verify

    return(status)
  }
)

setMethod("estimate", 
  signature = "stmodelKM",
  definition = function(.Object, sp, ...) {
    modtype   <- sp@win@type
    if (modtype == "sliding_events") {
      stop("Currently event-based sliding windows are available only for competing risks analyses.")
    }
    nsubpop   <- sp@nsubpop
    subpop    <- sp@subpop
    coltrt    <- .Object@coltrt
    trts      <- .Object@trts
    timePoint <- .Object@timePoint
    survTime  <- .Object@survTime
    censor    <- .Object@censor
 
    # set up internal assignment
    ntrts     <- length(trts)
    txassign  <- rep(NA, length(coltrt))

    for (j in 1:ntrts) 
      txassign[which(coltrt == trts[j])] <- j
    
    # set up return structure
    Obs <- list(sObs = rep(0,nsubpop),
          sSE  = rep(0,nsubpop),
          oObs = 0,
          oSE  = 0)
    TrtEff <- array(list(Obs),ntrts)
    
    HRs <- list(skmw     = 0,
        logHR    = rep(0,nsubpop),
        logHRSE  = rep(0,nsubpop),
        ologHR   = 0,
        ologHRSE = 0,
        logHRw   = 0)
    Ratios <- array(list(HRs),ntrts-1)

    # Estimate treatment effect for each treatment j
    for (j in 1:ntrts) {
      skmObs <- rep(0,nsubpop)
      skmSE  <- rep(0,nsubpop)
          for (i in 1:nsubpop) {
        trteff    <- tpest(kmest(survTime[txassign==j & subpop[,i]==1],
                               censor[txassign==j & subpop[,i]==1]),timePoint)
                skmObs[i] <- max(trteff$est,0)
              skmSE[i]  <- sqrt(trteff$var)
      }
      overalltrteff <- tpest(kmest(survTime[txassign==j],censor[txassign==j]),timePoint)
          overallSkmObs <- max(overalltrteff$est, 0)
          overallSkmSE  <- sqrt(overalltrteff$var)

      # check to make sure that subpopulation estimates are OK.
          if (sum(is.na(skmObs)) != 0 | is.na(overallSkmObs) != FALSE) {
                  cat("\n")
                print(paste("Unable to estimate survival time at ", timePoint, " time-unit(s) for trt ",
      j, "because there are too few events within one or more subpopulation(s)."))
                print(paste("The problem may be avoided by constructing larger subpopulations and/or by selecting a different timepoint for estimation."))
                stop()
          }

      TrtEff[[j]]=list(sObs=skmObs, sSE=skmSE, oObs=overallSkmObs, oSE=overallSkmSE)
    }

    # Estimate the relative treatment effect comparing each one with trt 1
    if (ntrts > 1) {
      for (j in 2:ntrts) {
        txassign <- rep(-1, length(coltrt))
        txassign[which(coltrt == trts[1])] <- 1
        txassign[which(coltrt == trts[j])] <- 0

        logHR    <- rep(0,nsubpop)
        logHRSE  <- rep(0,nsubpop)
        sel      <- txassign==1 | txassign==0 #j
        for (i in 1:nsubpop) {
          seli       <- sel & subpop[,i]==1
          LogRank    <- survdiff(Surv(survTime[seli], censor[seli]) ~ txassign[seli])
          logHR[i]   <- -(LogRank$obs[1] - LogRank$exp[1])/LogRank$var[1, 1]
          logHRSE[i] <- sqrt(1/LogRank$var[1, 1])
        }
        LogRank        <- survdiff(Surv(survTime[sel], censor[sel]) ~ txassign[sel])
        overallLogHR   <- -(LogRank$obs[1] - LogRank$exp[1])/LogRank$var[1, 1]
        overallLogHRSE <- sqrt(1/LogRank$var[1, 1])
        logHRw <- sum(logHR/logHRSE)
        skmw   <- sum((TrtEff[[1]]$sObs - TrtEff[[j]]$sObs)/sqrt(TrtEff[[1]]$sSE^2 + TrtEff[[j]]$sSE^2))

        Ratios[[j-1]] = list(skmw     = skmw,
                             logHR    = logHR,
                             logHRSE  = logHRSE,
                             ologHR   = overallLogHR,
                             ologHRSE = overallLogHRSE,
                             logHRw   = logHRw)

      }
    }

    # return the estimate as a list
    estimate <- list(model    = "KMe",
                     ntrts    = ntrts,
                     TrtEff   = TrtEff,
                     Ratios   = Ratios)
    return(estimate)
  }
)

setMethod("test",
  signature = "stmodelKM",
  definition = function(.Object, nperm = 2500, sp, effect, showstatus = TRUE, Cox = FALSE, MM = NULL, ...) {

    test <- NULL

    # no permutation test is done if nperm is 0; return immediately
    if (nperm > 0 & length(.Object@trts) > 1) {
      win   <- sp@win
      nsubpop <- sp@nsubpop
      osubpop <- sp@subpop
      coltrt  <- .Object@coltrt
      trts  <- .Object@trts
      survTime  <- .Object@survTime
      timePoint <- .Object@timePoint
      censor    <- .Object@censor

      # set up internal assignment
      ntrts <- length(trts)
      txassign  <- rep(-1, length(coltrt))  # do not use NA here

      for (j in 1:ntrts) {
        txassign[which(coltrt == trts[j])] <- j
      }

      # set up return structure
      tsigma    <- matrix(rep(0, (nperm * nsubpop)), ncol = nsubpop)    
      PermTest  <- list(sigma      = tsigma,
                        HRsigma    = tsigma,
                        pvalue     = 0,
                        chi2pvalue = 0,
                        HRpvalue   = 0)
      Res       <- array(list(PermTest),ntrts-1)

      for (j in 2:ntrts) {
      #
      #   do the permutations
      # compare trt1 with the rest
      #
          differences <- matrix(rep(0, (nperm * nsubpop)), ncol = nsubpop)
            logHRs      <- matrix(rep(0, (nperm * nsubpop)), ncol = nsubpop)
            no          <- 0
            p           <- 0
            terminate   <- 0
            Ntemp       <- nrow(osubpop)
            IndexSet1   <- (1:Ntemp)[txassign == 1]
            IndexSet2   <- (1:Ntemp)[txassign == j]

        if (Cox) {
          if (is.null(MM)) {
            fmla <- as.formula(paste("Surv(survTime,censor) ~ txassign"))
          } else {
            xnam <- colnames(MM)[2:dim(MM)[2]]
            fmla <- as.formula(paste("Surv(survTime,censor) ~ txassign + ", paste(xnam, collapse= "+")))
          }
        }

        if (showstatus) {
          title <- paste("\nComputing the p-value with", nperm)
          title <- paste(title, "number of permutations comparing trt", trts[j], "with trt", trts[1], "\n")
          cat(title)
            pb    <- txtProgressBar(min=0, max=nperm-1, style=3)
        }

        while (no < nperm) {
          ## PERMUTE THE VECTOR OF SUBPOPULATIONS WITHIN TREATMENTS ##
          if (showstatus) setTxtProgressBar(pb, no)

          Subpop        <- as.matrix(osubpop)
          permuteSubpop <- matrix(0, nrow = Ntemp, ncol = nsubpop)
          permuteSubpop[txassign == 1] <- Subpop[sample(IndexSet1),]
          permuteSubpop[txassign == j] <- Subpop[sample(IndexSet2),]
          subpop        <- permuteSubpop
          skm1          <- rep(0, nsubpop)
          skm2          <- rep(0, nsubpop)
          skmSE1        <- rep(0, nsubpop)
          skmSE2        <- rep(0, nsubpop)
          slogHRs       <- rep(0, nsubpop)
          slogHRSEs     <- rep(0, nsubpop)

          for (i in 1:nsubpop) {
            seff1   <- tpest(kmest(survTime[txassign == 1 & subpop[, i] == 1], censor[txassign == 1 & subpop[, i] == 1]), timePoint)
            seff2   <- tpest(kmest(survTime[txassign == j & subpop[, i] == 1], censor[txassign == j & subpop[, i] == 1]), timePoint) 
            skm1[i]     <- max(seff1$est, 0)
            skm2[i]     <- max(seff2$est, 0)
            skmSE1[i]   <- sqrt(seff1$var)
            skmSE2[i]   <- sqrt(seff2$var)
            sel         <- (txassign == 1 | txassign == j) & subpop[,i]==1
            if (Cox) {
              CoxM        <- coxph(fmla, subset=(sel))
              if (is.null(MM)) {
                slogHRs[i]    <- CoxM$coefficient
                slogHRSEs[i]  <- sqrt(CoxM$var)
              } else {
                slogHRs[i]    <- CoxM$coefficient[1]
                slogHRSEs[i]  <- sqrt(CoxM$var[1,1])
              }
            } else {
              LogRank         <- survdiff(Surv(survTime[sel],censor[sel]) ~ txassign[sel])
              slogHRs[i]      <- -(LogRank$obs[1]-LogRank$exp[1])/LogRank$var[1,1]
              slogHRSEs[i]    <- sqrt(1/LogRank$var[1,1])
            }
          }
          selo <- (txassign == 1 | txassign == j)

          slogHRw <- sum(slogHRs/slogHRSEs)

          overallSkm1 <- max(tpest(kmest(survTime[txassign == 1], censor[txassign == 1]), timePoint)$est, 0)
          overallSkm2 <- max(tpest(kmest(survTime[txassign == j], censor[txassign == j]), timePoint)$est, 0)

          if (Cox) {
            CoxM <- coxph(fmla, subset=(selo))
            if (is.null(MM)) {
              overallSLogHR    <- CoxM$coefficient
              overallLogHRSE   <- sqrt(CoxM$var)
            } else {
              overallSLogHR    <- CoxM$coefficient[1]
              overallSLogHRSE  <- sqrt(CoxM$var[1,1])
            }
          } else {
            LogRank            <- survdiff(Surv(survTime[selo],censor[selo]) ~ txassign[selo])
            overallSLogHR      <- -(LogRank$obs[1]-LogRank$exp[1])/LogRank$var[1,1]
            overallSLogHRSE    <- sqrt(1/LogRank$var[1,1])
          }
          if (sum(is.na(skm1)) == 0 & sum(is.na(skm2)) == 0 & is.na(overallSkm1) == FALSE & is.na(overallSkm2) == FALSE) {
            no <- no + 1
            p <- p + 1
            for (s in 1:nsubpop) {
              differences[p, s] <- (skm1[s] - skm2[s]) - (overallSkm1 - overallSkm2)
              logHRs[p,s]       <- slogHRs[s] - overallSLogHR
            }
          }
          terminate <- terminate + 1
          if (terminate >= nperm + 10000) {
            print(paste("After permuting", nperm, "plus 10000, or", 
                  nperm + 10000, "times, the program is unable to generate the permutation distribution based on", 
                  nperm, "permutations of the data"))
            print(paste("Consider creating larger subpopulations or selecting a different timepoint for estimation"))
            print(paste("when comparing trt", trts[1], "vs. trt", trts[j]))
            stop()
          }
        }
        # generating the sigmas and p-values
        sObs   <- effect$TrtEff[[1]]$sObs - effect$TrtEff[[j]]$sObs
        oObs   <- effect$TrtEff[[1]]$oObs - effect$TrtEff[[j]]$oObs
        logHR  <- effect$Ratios[[j-1]]$logHR
        ologHR <- effect$Ratios[[j-1]]$ologHR

        mname  <- paste0("SP", seq(1, dim(differences)[2]), "-Overall")
        if (win@type == "tail-oriented") {
          # remove the trivial case of the full cohort
          rm <- length(win@r1) + 1
          differences <- differences[, -rm]
          mname <- mname[-rm]
          sObs <- sObs[-rm]
          oObs <- oObs[-rm]
          logHRs <- logHRs[,-rm]
          logHR <- logHR[-rm]
          ologHR <- ologHR[-rm]
        }

        sigmas                  <- ssigma(differences)
        rownames(sigmas$sigma)  <- mname
        colnames(sigmas$sigma)  <- mname
        Res[[j-1]]$sigma        <- sigmas$sigma
        Res[[j-1]]$chi2pvalue   <- ppv (differences, sigmas$sigmainv, sObs, oObs, nperm)
        Res[[j-1]]$pvalue       <- ppv2(differences, sObs, oObs, nperm)
        sigmas                  <- ssigma(logHRs)
        rownames(sigmas$sigma)  <- mname
        colnames(sigmas$sigma)  <- mname
        Res[[j-1]]$HRsigma      <- sigmas$sigma
        Res[[j-1]]$HRpvalue     <- ppv2(logHRs, logHR, ologHR, nperm)

        if (showstatus) close(pb)
      }
      test = list(model  = "KMt",
                  ntrts  = ntrts,
                  Res    = Res)
    }
    return(test)
  }
)

#
# printing support functions for KM model
#
print.estimate.KM <- function(x, timePoint, trts) {
  for (j in 1:x@effect$ntrts) {
    # cat("\n")
    if (x@effect$ntrts == 1) {
      write(paste("Survival estimates",
        "at time point", timePoint), file = "")
    } else {
      write(paste("Survival estimates for treatment group", trts[j], 
        "at time point", timePoint), file = "")
    }

    overall_lbl <- -1
    if (x@subpop@win@type == "tail-oriented") {
      if (length(x@subpop@win@r1) == 1) {
        overall_lbl <- 1
      } else if (length(x@subpop@win@r2) == 1) {
        overall_lbl <- length(x@subpop@win@r1) + 1
      }
    }

    temp <- matrix(c(1:x@subpop@nsubpop, round(x@effect$TrtEff[[j]]$sObs, digits = 4),
      round(x@effect$TrtEff[[j]]$sSE, digits = 4)), ncol = 3)
    write("                         Survival", file = "")
    write("     Subpopulation     Probability      Std. Err.", file = "")
    for (i in 1:x@subpop@nsubpop) {
      if (x@subpop@win@type == "tail-oriented" & i == overall_lbl){
        write(paste(format(temp[i, 1], width = 12), format(temp[i,2], width = 19, nsmall = 4),
          format(temp[i, 3], width = 15, nsmall = 4), "(entire cohort)"), file = "")
        } else {
        write(paste(format(temp[i, 1], width = 12), format(temp[i,2], width = 19, nsmall = 4),
          format(temp[i, 3], width = 15, nsmall = 4)), file = "")
      }
    }
    if (x@subpop@win@type != "tail-oriented"){
      write(paste("        Overall", 
        format(round(x@effect$TrtEff[[j]]$oObs, digits = 4), nsmall = 4, width = 16), 
        format(round(x@effect$TrtEff[[j]]$oSE, digits = 4), nsmall = 4, width = 15)), 
        file = "")
      }
    cat("\n")
  }

  if (x@effect$ntrts > 1) {
    # cat("\n")
    write("Survival differences at time point and hazard ratio estimates", file = "")

    for (j in 2:x@effect$ntrts) {
      cat("\n")
      write(paste("trt", trts[1], "vs. trt", trts[j]), file = "")

      temp <- matrix(c(1:x@subpop@nsubpop, 
      round(x@effect$TrtEff[[1]]$sObs - x@effect$TrtEff[[j]]$sObs, digits = 4), 
        round(sqrt(x@effect$TrtEff[[1]]$sSE^2 + x@effect$TrtEff[[j]]$sSE^2), digits = 4)), ncol = 3)
      cat("\n")
      write(paste("Survival differences at time point", timePoint), file = "")

      write(paste("Comparing trt", trts[1], "vs. trt", trts[j]), file = "")
      cat("\n")
      write("                         Survival", file = "")
      write("     Subpopulation      Difference      Std. Err.", 
        file = "")
      for (i in 1:x@subpop@nsubpop) {
        if (x@subpop@win@type == "tail-oriented" & i == overall_lbl){
          write(paste(format(temp[i, 1], width = 12), format(temp[i, 2], width = 19, nsmall = 4), 
            format(temp[i, 3], width = 15, nsmall = 4), "(entire cohort)"), file = "")
        } else {
          write(paste(format(temp[i, 1], width = 12), format(temp[i, 2], width = 19, nsmall = 4), 
            format(temp[i, 3], width = 15, nsmall = 4)), file = "")
        }
      }
      if (x@subpop@win@type != "tail-oriented") {
        write(paste("        Overall", 
          format(round(x@effect$TrtEff[[1]]$oObs - x@effect$TrtEff[[j]]$oObs, digits = 4), nsmall = 4, width = 16), 
          format(round(sqrt(x@effect$TrtEff[[1]]$oSE^2 + x@effect$TrtEff[[j]]$oSE^2), digits = 4), nsmall = 4,
            width = 15)), file = "")
      }
      cat("\n")
      write("Hazard ratio estimates", file = "")
      temp <- matrix(c(1:x@subpop@nsubpop, 
        round(x@effect$Ratios[[j-1]]$logHR, digits = 6), 
        round(x@effect$Ratios[[j-1]]$logHRSE, digits = 6), 
        round(exp(x@effect$Ratios[[j-1]]$logHR), digits = 2)), ncol = 4)
      write("     Subpopulation        Log HR       Std. Err.       Hazard Ratio", file = "")
      for (i in 1:x@subpop@nsubpop) {
        if (x@subpop@win@type == "tail-oriented" & i == overall_lbl){
          write(paste(format(temp[i, 1], width = 12),
            format(temp[i, 2], width = 19, nsmall = 6),
            format(temp[i, 3], width = 14, nsmall = 6),
            format(temp[i, 4], width = 15, nsmall = 2), "(entire cohort)"),
            file = "")
        } else {
          write(paste(format(temp[i, 1], width = 12),
            format(temp[i, 2], width = 19, nsmall = 6),
            format(temp[i, 3], width = 14, nsmall = 6),
            format(temp[i, 4], width = 15, nsmall = 2)),
            file = "")
        }
      }
      if (x@subpop@win@type != "tail-oriented"){
        write(paste("        Overall",
          format(round(x@effect$Ratios[[j-1]]$ologHR,      digits = 6), nsmall = 6, width = 16),
          format(round(x@effect$Ratios[[j-1]]$ologHRSE,    digits = 6), nsmall = 6, width = 14),
          format(round(exp(x@effect$Ratios[[j-1]]$ologHR), digits = 2), nsmall = 2, width = 15)),
          file = "")
      }
      cat("\n")
    }
  }
}

print.cov.KM <- function(stobj, timePoint, trts) {
  if (!is.null(stobj@result)) {
    for (j in 1:(stobj@result$ntrts - 1)) {
      ns <- stobj@subpop@nsubpop
      if (stobj@subpop@win@type == "tail-oriented") ns <- ns - 1
      # cat("\n")
      write(paste("The covariance matrix of the Kaplan-Meier differences at", 
                  timePoint, "time units for the", ns, "subpopulations is:"), 
                  file = "")
      write(paste("trt ", trts[1], "vs. trt", trts[j + 1]), file = "")
      print(stobj@result$Res[[j]]$sigma)
        
      cat("\n")
      write(paste("The covariance matrix of the log hazard ratios for the", 
                  ns, "subpopulations is:"), file = "")
      print(stobj@result$Res[[j]]$HRsigma)
      cat("\n")

      # write(paste("The covariance matrix (based on homogeneous association) of the Kaplan-Meier differences at", 
      #             timePoint, "time units for the", stobj@subpop@nsubpop, "subpopulations is:"), file = "")
      # print(stobj@result$hasigma)
        
      # cat("\n")
      # write(paste("The covariance matrix (based on homogeneous association) of the log hazard ratios for the", 
      #             stobj@subpop@nsubpop, "subpopulations is:"), file = "")
      # print(stobj@result$haHRsigma)
      # cat("\n")
    }
  }
}

print.stat.KM <- function(stobj, trts) {
  if (!is.null(stobj@result)) {
    for (j in 1:(stobj@result$ntrts - 1)) {
      t <- stobj@result$Res[[j]]
      # cat("\n")
      write(paste("Supremum test results"), file = "")
      write(paste("trt", trts[1], "vs. trt", trts[j + 1]), file = "")
      write(paste("Interaction p-value based on Kaplan-Meier estimates:", t$pvalue), file = "")
      write(paste("Interaction p-value based on hazard ratio estimates:", t$HRpvalue), file = "")

      cat("\n")
      write(paste("Chi-square test results"), file = "")
      write(paste("Interaction p-value based on Kaplan-Meier estimates:", 
                  t$chi2pvalue), file = "")

      # cat("\n")
      # write(paste("Homogeneous association test results"), file = "")
      # write(paste("Interaction p-value based on Kaplan-Meier estimates:", hapvalue), file = "")
      # write(paste("Interaction p-value based on hazard ratio estimates:", haHRpvalue), file = "")

      cat("\n")
    }
  }
}

setMethod("print",
  signature = "stmodelKM",
  definition = function(x, stobj, estimate = TRUE, cov = TRUE, test = TRUE, ...) {
    ntrts <- length(x@trts)

    #
    #  1. estimates
    #
    if (estimate) {
      print.estimate.KM(stobj, x@timePoint, x@trts)
    }

    #
    #   2. covariance matrices
    #
    if (cov & ntrts > 1) {
      print.cov.KM(stobj, x@timePoint, x@trts)
    }
  
    #
    #   3. Supremum test and Chi-square test results
    #
    if (test & ntrts > 1) {
      print.stat.KM(stobj, x@trts)
    }
   }
)

# constructor function for stmodelKM
stepp.KM <- function(coltrt, survTime, censor, trts, timePoint) {
  model <- new("stmodelKM", coltrt = coltrt, survTime = survTime, censor = censor,
      trts = trts, timePoint = timePoint)
  return(model)
}
