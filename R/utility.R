balance_patients <- function(range_r1, range_r2, maxnsubpops, covar, verbose = FALSE,
  plot = FALSE, contour = FALSE, nlevels = 5, border = FALSE) {
  unbalance <- function(r1, r2, maxnsubp, frq, allfrq) {
    # We work with indices, not values
    # Note: indices of values are from 1 to k
    k <- nrow(frq)
    subpops <- matrix(rep(NA, 5*maxnsubp), ncol = 5)
    colnames(subpops) <- c("IndLow", "IndHigh", "Size", "CutoffLow", "CutoffHigh")
    subpops[1, 1] <- 1
    # collect at least r2 subjects
    lasttemp <- sum(allfrq[1, ] < r2) + 1
    subpops[1, 2] <- lasttemp
    subpops[1, 3] <- allfrq[1, lasttemp]
    # drop at least (r2 - r1) subjects
    for (i in 2:maxnsubp) {
      preva <- subpops[i - 1, 1]
      prevb <- subpops[i - 1, 2]
      # First we find the potential new a 
      a <- prevb - sum(allfrq[preva:prevb, prevb] <= r1) + 1 # Careful: this could be >= prevb (CHECK)
      if (prevb < k) {
        if (allfrq[a, k] >= r2) { # if there is enough to build a new subpopulation
          subpops[i, 1] <- a
          b <- a + sum(allfrq[a,a:k] < r2) 
          subpops[i, 2] <- b
          subpops[i, 3] <- allfrq[a, b]
        } else {
          # last subpopulation
          subpops[i, 1] <- a
          b <- k
          subpops[i, 2] <- b
          subpops[i, 3] <- allfrq[a, b]
          break
        }
      } else {
        i <- i - 1
        break
      }
    }

    nsubpops <- i
    subpopsfin <- subpops[1:nsubpops, ]
    # Compute the unbalance measure
    res_var <- var(subpopsfin[, 3])
    # Add the covariate values that define the subpopulations
    subpopsfin[, 4] <- frq[subpopsfin[, 1], 1]
    subpopsfin[, 5] <- frq[subpopsfin[, 2], 1]

    return(list(var = res_var, nsubpops = nsubpops, subpops = subpopsfin))
  }

  if (length(range_r1) != 2) {
    stop("range_r1 must contain two elements.")
  }
  if (length(range_r2) != 2) {
    stop("range_r2 must contain two elements.")
  }
  if (diff(range_r1) < 0) {
    stop("the elements of range_r1 must be in ascending order.")
  }
  if (diff(range_r2) < 0) {
    stop("the elements of range_r2 must be in ascending order.")
  }
  if (range_r1[1] < 2) {
    stop("the first element of range_r1 must be at least 2.")
  }
  if (range_r2[1] < 2) {
    stop("the first element of range_r2 must be at least 2.")
  }
  range_r1 <- as.integer(range_r1)
  range_r2 <- as.integer(range_r2)
  if (maxnsubpops < 2) {
    stop("maxnsubpops must be at least 2.")
  }
  maxnsubpops <- as.integer(maxnsubpops)

  freqdist <- data.frame(table(covar))
  freqdist[, 1] <- as.numeric(as.character(freqdist[, 1]))
  freqdist$CumFreq <- cumsum(freqdist[, 2])

  # freq [a,b] for all combinations - do only once!
  k <- nrow(freqdist)
  allfreqs <- matrix(rep(NA, k*k), nrow = k)
  allfreqs[1,] <- freqdist$CumFreq
  for (i in 2:k) {
    allfreqs[i, i:k] <- freqdist$CumFreq[i:k] - freqdist$CumFreq[i - 1]
  }

  bestr1 <- bestr2 <- rep(NA, maxnsubpops)
  resmat <- matrix(rep(NA,
    5*(range_r1[2] - range_r1[1] + 1)*(range_r2[2] - range_r2[1] + 1)), ncol = 5)
  minvars <- rep(1e+99, maxnsubpops)
  nextind <- 1
  indbest <- 0
  varbest <- 1e+99
  for (i in range_r1[1]:range_r1[2]) {
    for (j in range_r2[1]:range_r2[2]) {
      resunb <- unbalance(i, j, maxnsubpops, freqdist, allfreqs)
      newvar <- resunb[[1]]
      resmat[nextind, 1] <- i
      resmat[nextind, 2] <- j
      resmat[nextind, 3] <- newvar
      resmat[nextind, 4] <- resunb[[2]]
      if (newvar < minvars[resunb[[2]]]) {
        bestr1[resunb[[2]]] <- i
        bestr2[resunb[[2]]] <- j
        minvars[resunb[[2]]] <- newvar
      }
      if (newvar < varbest) {
        r1best <- i
        r2best <- j
        indbest <- resunb[[2]]
        varbest <- newvar
      }
      nextind <- nextind + 1
    }
  }

  if (verbose) {
    cat("Balanced subpopulations determination (number of patients)\n")
    cat("\n")
    cat(paste("Range for number of patients in common in consecutive subpopulations (r1):", range_r1[1], "<-->", range_r1[2], "\n"))
    cat(paste("Range for number of patients per subpopulation (r2):", range_r2[1], "<-->", range_r2[2]), "\n")
    cat("\n")

    for (i in 1:maxnsubpops) {
      if (!is.na(bestr1[i])) {
        best <- unbalance(bestr1[i], bestr2[i], maxnsubpops, freqdist, allfreqs)
        cat(paste("Best result for", i, "subpopulations\n"))
        cat("  * Optimal r1 value:", bestr1[i], "\n")
        cat("  * Optimal r2 value:", bestr2[i], "\n")
        cat("  * Minimum variance of subpopulation sizes achieved:", minvars[i], "\n")
        cat("  * Subpopulation summary information\n")
        cat("                                 Covariate Summary                 Sample\n")
        cat("    Subpopulation        Median       Minimum       Maximum          size\n")
        temp <- matrix(NA, nrow = i, ncol = 5)
        temp[, 1] <- 1:i
        for (d in 1:i) {
          covar_tmp <- covar[(covar >= best$subpops[d, 4]) & (covar <= best$subpops[d, 5])]
          temp[d, 2] <- round(median(covar_tmp), digits = 2)
        }
        temp[, 3] <- round(best$subpops[, 4], digits = 4)
        temp[, 4] <- round(best$subpops[, 5], digits = 4)
        temp[, 5] <- best$subpops[, 3]
        for (d in 1:i) {
          cat(paste(format(temp[d, 1], width = 11),
                    format(temp[d, 2], width = 19, nsmall = 2),
                    format(temp[d, 3], width = 13, nsmall = 4),
                    format(temp[d, 4], width = 13, nsmall = 4),
                    format(temp[d, 5], width = 13), "\n"))
        }
        cat("\n")
      }
    }

    cat(paste("Overall best result\n"))
    cat(paste("  * Number of subpopulations:", indbest, "\n"))
    cat(paste("  * Best r1 value:", r1best, "\n"))
    cat(paste("  * Best r2 value:", r2best, "\n"))
    cat(paste("  * Minimum variance of subpopulation sizes achieved:", round(varbest, digits = 4), "\n"))
  }

  if (plot) {
    # Add the shade of color to the last column: darker means smaller variance.
    # It works on the log scale to highlight best results only (darker).
    # There are 8 (i.e. nlogvar - 1) shades of color from 0/8 to 7/8 (8/8 is white).
    nlogvar <- 9
    logvar <- log(resmat[, 3])
    cutoffs <- seq(min(logvar), max(logvar), length.out = nlogvar)
    resmat[, 5] <- sapply(as.list(logvar),
      function(x) (sum(cutoffs <= x) - 1)/(nlogvar - 1))

    if (diff(range(range_r1)) >= 10 && diff(range(range_r2)) >= 10) {
      def.par <- par(no.readonly = TRUE)
      nf <- graphics::layout(c(1, 2), heights = c(10, 2))
      nsubpops <- sort(unique(resmat[, 4]))
      clrs <- scales::viridis_pal(option = "viridis", direction = -1)(length(nsubpops))
      if (contour) {
        var_arr <- array(NA, dim = c(max(range_r1), max(range_r2), maxnsubpops),
          dimnames = list(r1 = 1:max(range_r1), r2 = 1:max(range_r2),
            nsubpops = 1:maxnsubpops))
      }
      if (border) {
        r2_arr <- array(NA, dim = c(max(range_r1), max(range_r2), maxnsubpops),
          dimnames = list(r1 = 1:max(range_r1), r2 = 1:max(range_r2),
            nsubpops = 1:maxnsubpops))
      }
      
      # plot heatmap with subpopulation variances
      par(mar = c(5, 5, 4, 1))
      plot(range_r1, range_r2, type = "n",
           xlab = expression(italic(r)[1]), ylab = expression(italic(r)[2]),
           main = "Subpopulations (colors, see the legend below)\nand variances (shades of color, darker is smaller)")
      for (i in 1:length(nsubpops)) {
        resmat_sub <- resmat[resmat[, 4] == nsubpops[i], ]
        r1_unique <- sort(unique(resmat_sub[, 1]))
        rect(resmat_sub[, 1] - .5, resmat_sub[, 2] - .5,
             resmat_sub[, 1] + .5, resmat_sub[, 2] + .5,
             col = scales::alpha(clrs[i], alpha = 1 - resmat_sub[, 5]),
             border = NA)
        if (contour) {
          for (j in 1:length(r1_unique)) {
            r2_unique <- sort(unique(resmat_sub[resmat_sub[, 1] == r1_unique[j], 2]))
            var_arr[r1_unique[j], r2_unique, nsubpops[i]] <- 
              resmat_sub[resmat_sub[, 1] == r1_unique[j], 3]
          }
        }
        if (border) {
          for (j in 1:length(r1_unique)) {
            r2_unique <- sort(unique(resmat_sub[resmat_sub[, 1] == r1_unique[j], 2]))
            r2_arr[r1_unique[j], r2_unique, nsubpops[i]] <- r2_unique
          }
        }
      }
      if (contour) {
        var_arr <- var_arr[-setdiff(1:max(range_r1), min(range_r1):max(range_r1)),
          -setdiff(1:max(range_r2), min(range_r2):max(range_r2)),
          -setdiff(1:maxnsubpops, nsubpops), drop = FALSE]
        for (i in 1:length(nsubpops)) {
          contour(min(range_r1):max(range_r1), min(range_r2):max(range_r2),
            var_arr[, , i], add = TRUE, nlevels = nlevels, col = gray(.4, alpha = .6))
        }
      }
      if (border) {
        r2_arr <- r2_arr[-setdiff(1:max(range_r1), min(range_r1):max(range_r1)),
          -setdiff(1:max(range_r2), min(range_r2):max(range_r2)),
          -setdiff(1:maxnsubpops, nsubpops), drop = FALSE]
        r2_min <- apply(r2_arr, c(1, 3),
          function(x) {
            tmp <- x[!is.na(x)]
            if (length(tmp)) {
              res <- min(tmp, na.rm = TRUE)
              ifelse(is.infinite(res), NA, res)
            } else {
              NA
            }
          }
        )
        for (i in 1:length(nsubpops)) {
          brd_dat <- na.omit(cbind(min(range_r1):max(range_r1), r2_min[, i]))
          brd_x <- brd_dat[, 1] - .5
          brd_x[1] <- brd_x[1] + .5
          brd_x[nrow(brd_dat)] <- brd_x[nrow(brd_dat)] + .5
          brd_x <- c(brd_x, brd_x[nrow(brd_dat)])
          brd_y <- brd_dat[, 2] - .5
          brd_y <- c(brd_y, max(range_r2))
          lines(brd_x, brd_y, type = "s", lwd = 2, col = gray(.5))
        }
      }
      points(bestr1, bestr2, pch = 1, col = "darkblue", lwd = 2)
      points(bestr1[indbest], bestr2[indbest], pch = 4, col = "red", lwd = 3, cex = 1.1)
      
      # plot legend
      par(mar = c(1.5, 5, 1, 1))
      plot(nsubpops, rep(1, length(nsubpops)), type = "n", xaxt = "n", yaxt = "n",
        bty = "n", xlim = c(min(nsubpops) - .5, max(nsubpops) + .5), ylim = c(.5, 3),
        xlab = "", ylab = "")
      title(main = "Number of subpopulations", cex.main = .75, line = 0)
      rect(nsubpops - .5, .5, nsubpops + .5, 1.5, col = scales::alpha(clrs, alpha = .6),
        border = "black")
      axis(side = 1, at = nsubpops, tick = FALSE, line = -.75, cex.axis = .6)
      
      par(def.par)
    } else {
      warning("the plot is produced only when at least 10 values for both r1 and r2 are provided.")
    }
  }
  
  return(list(r1_best = r1best, r2_best = r2best, var_best = varbest,
    nsubpops_best = indbest, all_res = resmat))
}
