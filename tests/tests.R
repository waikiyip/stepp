library(stepp)

# library(steppevent)
# source("/Users/Sergio/Dropbox (Personal)/stepp/packages/lazar examples/STEPPpermCI2_event algorithm.R")

### Example 1 ###

data(bigKM)
rxgroup <- bigKM$trt
time    <- bigKM$time
evt     <- bigKM$event
cov     <- bigKM$ki67

## the following code performs the calculations using the steppevent package
# set.seed(101)
# noperm <- 100
# res_old <- STEPPpermCI2(coltrt = rxgroup, coltime = time, coltype = evt,
#   colvar = cov, trts = c(1, 2), eventpops = 20, mineventpops = 10,
#   timest = 4, noperm = noperm, minRequiredSubpops = 5, legendy = 30,
#   pline = -2.5, color = c("red", "black"), ylabel = "4 year DFS",
#   xlabel = "Subpopulations by Median Ki-67", ncex = 0.7,
#   tlegend = c("Tamoxifen", "Letrozole"), nlas = 3, alpha = 0.05,
#   pointwise = FALSE)

# generate event-based windows
swin_e <- new("stwin", type = "sliding_events", e1 = 10, e2 = 20)
subp_e <- new("stsubpop")
subp_e <- generate(subp_e, win = swin_e, covariate = cov, coltype = evt,
                   coltrt = rxgroup, trts = c(1, 2), minsubpops = 5)
summary(subp_e)

# estimate and test the CI model
steppes_e <- new("steppes")
modCI_e <- new("stmodelCI", coltrt = rxgroup, trts = c(1, 2), coltime = time,
  coltype = evt, timePoint = 4)
resCI_e <- estimate(steppes_e, subp_e, modCI_e)
print(resCI_e)
set.seed(101)
noperm <- 100
resCI_e <- test(resCI_e, noperm)
print(resCI_e)
plot(resCI_e, legendy = 50,
     pline = -2.5, color = c("red", "black"), ylabel = "4 year DFS",
     xlabel = "Subpopulations by Median Ki-67", ncex = 0.7,
     tlegend = c("Tamoxifen", "Letrozole"), nlas = 3, alpha = 0.05,
     pointwise = FALSE)

# estimate and test the KM model
steppes_e <- new("steppes")
modKM_e <- new("stmodelKM", coltrt = rxgroup, trts = c(1, 2), survTime = time,
  censor = evt, timePoint = 4)
resKM_e <- estimate(steppes_e, subp_e, modKM_e)
print(resKM_e)
set.seed(101)
noperm <- 100
resKM_e <- test(resKM_e, noperm)
print(resKM_e)
plot(resKM_e, ylabel = "4 year DFS",
     xlabel = "Subpopulations by Median Ki-67", ncex = 0.7,
     tlegend = c("Tamoxifen", "Letrozole"), nlas = 3, alpha = 0.05)

### Example 2 ###

# GENERATE THE DATA 
n <- 1000  # set the sample size 
mu <- 0   # set the mean and sd of the covariate 
sigma <- 1 
 
beta0 <- log(-log(0.5)) # set the intercept for the log hazard 
beta1 <- -0.2  # set the slope on the covariate 
beta2 <- 0.5  # set the slope on the treatment indicator 
beta3 <- 0.7  # set the slope on the interaction 

prob2 <- 0.2  # set the proportion type 2 events
cprob <- 0.3  # set the proportion censored

set.seed(7775432)  # set the random number seed
covariate <- rnorm(n,mean=mu,sd=sigma) # generate the covariate values
Txassign <- rbinom(n,1,0.5)  # generate the treatment indicator
x3 <- covariate*Txassign  # compute interaction term
lambda1 <- exp(beta0+beta1*covariate+beta2*Txassign+beta3*x3) # compute the hazard for type 1 event
lambda2 <- prob2*lambda1/(1-prob2) # compute the hazard for the type 2 event
# compute the hazard for censoring time
lambda0 <- cprob*(lambda1+lambda2)/(1-cprob) 
t1 <- rexp(n,rate=lambda1)  # generate the survival time for type 1 event 
t2 <- rexp(n,rate=lambda2)  # generate the survival time for type 2 event 
t0 <- rexp(n,rate=lambda0)  # generate the censoring time 
time <- pmin(t0,t1,t2)   # compute the observed survival time 
type <- rep(0,n) 
type[(t1 < t0)&(t1 < t2)] <- 1 
type[(t2 < t0)&(t2 < t1)] <- 2 
 
## the following code performs the calculations using the steppevent package
# STEPPpermCI2(coltrt = Txassign, coltime = time, coltype = type, colvar = covariate, trts = c(0, 1),
#   eventpops = 20, mineventpops = 10, timest = 1, noperm = 250, minRequiredSubpops = 5, legendy = 30,
#   pline = -2.5, color = c("red", "black"), ylabel = "PFS",
#   xlabel = "Subpopulations by Median Continuous Variable", ncex = 0.7,
#   tlegend = c("Trt A", "Trt B"), nlas = 3, alpha = 0.05, pointwise = FALSE) 

# generate event-based windows
swin_e <- new("stwin", type = "sliding_events", e1 = 10, e2 = 20)
subp_e <- new("stsubpop")
subp_e <- generate(subp_e, win = swin_e, covariate = covariate, coltype = type,
  coltrt = Txassign, trts = c(0, 1), minsubpops = 5)
summary(subp_e)

# estimate and test the CI model
steppes_e <- new("steppes")
modCI_e <- new("stmodelCI", coltrt = Txassign, trts = c(0, 1), coltime = time,
  coltype = type, timePoint = 1)
resCI_e <- estimate(steppes_e, subp_e, modCI_e)
print(resCI_e)
set.seed(101)
noperm <- 250
resCI_e <- test(resCI_e, noperm)
print(resCI_e)
plot(resCI_e, legendy = 30,
     pline = -2.5, color = c("red", "black"), ylabel = "PFS",
     xlabel = "Subpopulations by Median Continuous Variable", ncex = 0.7,
     tlegend = c("Trt A", "Trt B"), nlas = 3, alpha = 0.05,
     pointwise = FALSE)
