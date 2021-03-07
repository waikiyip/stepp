\name{simdataKM}
\alias{simdataKM}
\docType{data}
\title{
  Simulated data for Kaplan-Meier STEPP analysis.
}
\description{
  Simulated data to use in a Kaplan-Meier STEPP analysis.
}
\usage{data(simdataKM)}
\format{
{
  A \code{\link{data.frame}} object with columns:
  \describe{
    \item{time}{time-to-event data}
    \item{censor}{censoring indicator}
    \item{trt}{treatment indicator}
    \item{covar}{covariate values}
  }
}
}
\examples{
### the following code generates the data in the object ###

set.seed(101)
n <- 1000

# generate the treatment indicator
Txassign <- sample(c(1, 2), n, replace = TRUE, prob = rep(.5, 2))
n1 <- sum(Txassign == 1)
n2 <- n - n1

# generate the covariate values
mean_cov <- 55
sd_cov <- 7
covariate <- rnorm(n, mean = mean_cov, sd = sd_cov)

# generate the survival times and censoring indicator assuming a treatment-covariate interaction
Entry <- sort(runif(n, 0, 5))
SurvT1 <- .5
beta0 <- -65/75
beta1 <- 2/75
Surv <- rep(0, n)
lambda1 <- -log(SurvT1)/4
Surv[Txassign == 1] <- rexp(n1, lambda1)
Surv[Txassign == 2] <- rexp(n2, (lambda1*(beta0 + beta1*covariate[Txassign == 2])))
EventTimes <- rep(0, n)
EventTimes <- Entry + Surv
censor <- time <- rep(0, n)
for (i in 1:n) {
  censor[i] <- ifelse(EventTimes[i] <= 7, 1, 0)
  time[i] <- ifelse(EventTimes[i] < 7, Surv[i], 7 - Entry[i])
}
simdataKM <- data.frame(time = time, censor = censor, trt = Txassign, covar = covariate)

# overall survival analysis
simdataKM_surv <- Surv(simdataKM$time, simdataKM$censor)
simdataKM_fit <- survfit(simdataKM_surv ~ trt, data = simdataKM)
plot(simdataKM_fit, lty = 2:3, lwd = rep(2, 2), col = 2:3)
legend("topright", c("Treatment 1", "Treatment 2"), lty = 2:3, lwd = rep(2, 2), col = 2:3)
title("Kaplan-Meier Curves\nfor simulated data")
survdiff(simdataKM_surv ~ trt, data = simdataKM, rho = 0)
}
\keyword{datasets}
