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
#   tlegend = c(Letrozole", "Tamoxifen"), nlas = 3, alpha = 0.05,
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
     tlegend = c("Letrozole", "Tamoxifen"), nlas = 3, alpha = 0.05,
     pointwise = FALSE)

# # estimate and test the KM model
# steppes_e <- new("steppes")
# modKM_e <- new("stmodelKM", coltrt = rxgroup, trts = c(1, 2), survTime = time,
#   censor = evt, timePoint = 4)
# resKM_e <- estimate(steppes_e, subp_e, modKM_e)
# print(resKM_e)
# set.seed(101)
# noperm <- 100
# resKM_e <- test(resKM_e, noperm)
# print(resKM_e)
# plot(resKM_e, ylabel = "4 year DFS",
#      xlabel = "Subpopulations by Median Ki-67", ncex = 0.7,
#      tlegend = c("Letrozole", "Tamoxifen"), nlas = 3, alpha = 0.05)

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

### Example 3 ###

data(bigKM)

ranger2 <- c(100, 450)
ranger1 <- c(50, 300)
maxnsubpops <- 30

res_bal <- balance_patients(ranger1, ranger2, maxnsubpops, bigKM$ki67,
                            plot = TRUE, verbose = TRUE, contour = TRUE,
                            nlevels = 6)

# ### Examples for 'tail-oriented' windows ###
# 
# library(stepp)
# 
# set.seed(101)
# Y <- rnorm(100)
# summary(Y)
# 
# nsubpop <- 10
# tt <- gen.tailwin(Y, nsub = nsubpop, dir = "LE")
# ss <- stepp.win(type = "tail-oriented", r1 = tt$v, r2 = rep(min(Y), nsubpop))
# 
# # create and generate the stepp subpopulation
# sp <- new("stsubpop")
# sp <- generate(sp, win = ss, cov = Y)
# summary(sp)
# 
# # ---
# 
# nsubpop <- 10
# tt <- gen.tailwin(Y, nsub = nsubpop, dir = "GE")
# ss <- stepp.win(type = "tail-oriented", r1 = rep(max(Y), nsubpop), r2 = tt$v)
# 
# # create and generate the stepp subpopulation
# sp <- new("stsubpop")
# sp <- generate(sp, win = ss, cov = Y)
# summary(sp)
# 
# # ---
# 
# win1 <- stepp.win(type="sliding", r1=5,r2=99)
# 
# # create and generate the stepp subpopulation
# sp <- new("stsubpop")
# sp <- generate(sp, win = win1, cov = Y)
# summary(sp)
# 
# ###
# 
# # debugonce(generate, signature = "stsubpop")
# # debugonce(generate.all)
# # debug(gen.tailwin)

###

library(stepp)

data(bigKM)

# bigKM$ki67[1:80] <- NA
# cov_na <- which(is.na(bigKM$ki67))
# bigKM <- bigKM[-cov_na, ]

rxgroup <- bigKM$trt
time    <- bigKM$time
evt     <- bigKM$event
cov     <- bigKM$ki67

# analyze using Cumulative Incidence method with
# sliding window size of 150 patients and a maximum of 50 patients in common
#
nsubpop_tmp <- 10
win_tmp <- gen.tailwin(cov, nsub = nsubpop_tmp, dir = "GE")
nsubpop <- length(win_tmp$v)
swin <- new("stwin", type = "tail-oriented", r1 = rep(max(cov), nsubpop), r2 = win_tmp$v) # create a tail-oriented window
subp <- new("stsubpop")                             # create subpopulation object
subp <- generate(subp, win = swin, covariate = cov) # generate the subpopulations
summary(subp)					                    # summary of the subpopulations

# create a stepp model using Kaplan Meier Method to analyze the data
#
smodel  <- new("stmodelKM", coltrt=rxgroup, trts=c(1,2), survTime=time, censor=evt, timePoint=4)

statKM  <- new("steppes")		  # create a test object based on subpopulation and window
statKM  <- estimate(statKM, subp, smodel) # estimate the subpopulation results
# Warning: IT IS RECOMMEND TO USE AT LEAST 2500 PERMUTATIONS TO  PROVIDE STABLE RESULTS.
statKM <- test(statKM, nperm = 10)       # permutation test with 10 iterations

print(statKM)				  # print the estimates and test statistics
plot(statKM, ncex=0.65, legendy=30, pline=-15.5, color=c("blue","gold"),
     pointwise=FALSE, 
     xlabel="Median Ki-67 LI in Subpopulation (% immunoreactivity)",
     ylabel="4-year Disease Free Survival", 
     tlegend=c("Letrozole", "Tamoxifen"), nlas=3)

### Example for single group analysis ###
library(stepp)

data(bigKM)

# bigKM$ki67[1:80] <- NA
# cov_na <- which(is.na(bigKM$ki67))
# bigKM <- bigKM[-cov_na, ]

rxgroup <- bigKM$trt
time    <- bigKM$time
evt     <- bigKM$event
cov     <- bigKM$ki67

# analyze using Cumulative Incidence method with
# sliding window size of 150 patients and a maximum of 50 patients in common
#
nsubpop_tmp <- 10
win_tmp <- gen.tailwin(cov, nsub = nsubpop_tmp, dir = "GE")
nsubpop <- length(win_tmp$v)
swin <- new("stwin", type = "tail-oriented", r1 = rep(max(cov), nsubpop), r2 = win_tmp$v) # create a tail-oriented window
subp <- new("stsubpop")                             # create subpopulation object
subp <- generate(subp, win = swin, covariate = cov) # generate the subpopulations
summary(subp)					                    # summary of the subpopulations

# create a stepp model using Kaplan Meier Method to analyze the data
#
smodel  <- new("stmodelKM", survTime=time, censor=evt, timePoint=4)

statKM  <- new("steppes")		  # create a test object based on subpopulation and window
statKM  <- estimate(statKM, subp, smodel) # estimate the subpopulation results
statKM <- test(statKM, nperm = 10)       # permutation test with 10 iterations

print(statKM)				  # print the estimates and test statistics
plot(statKM, ncex=0.65, pline=-15.5,
     xlabel="Median Ki-67 LI in Subpopulation (% immunoreactivity)",
     ylabel="Survival Estimates", 
     nlas=3, subplot = TRUE)

###

data(bigCI)

rxgroup <- bigCI$trt
time    <- bigCI$time
evt     <- bigCI$event
cov     <- bigCI$ki67

# analyze using Cumulative Incidence method with
# sliding window size of 150 patients and a maximum of 50 patients in common
#
swin    <- new("stwin", type="sliding", r1=50, r2=150) # create a sliding window
subp    <- new("stsubpop")                             # create subpopulation object
subp    <- generate(subp, win=swin, covariate=cov) # generate the subpopulations
summary(subp)					   # summary of the subpopulations

# create a stepp model using Cumulative Incidences to analyze the data
#
smodel  <- new("stmodelCI", coltime=time, coltype=evt, timePoint=4)

statCI  <- new("steppes")		  # create a test object based on subpopulation and window
statCI  <- estimate(statCI, subp, smodel) # estimate the subpo10ulation results
# Warning: In this example, the permutations have been set to 0 to allow the function
# to finish in a short amount of time.  IT IS RECOMMEND TO USE AT LEAST 2500 PERMUTATIONS TO 
# PROVIDE STABLE RESULTS.
statCI  <- test(statCI, nperm=10)       # permutation test with 0 iterations

print(statCI)				  # print the estimates and test statistics
plot(statCI, ncex=0.65, pline=-15.5,
     xlabel="Median Ki-67 LI in Subpopulation (% immunoreactivity)",
     ylabel="4-year Cumulative Incidence", 
     nlas=3, subplot = TRUE)

###

data(aspirin)

# remove cases with missing data
aspirinc <- aspirin[complete.cases(aspirin), ]

# make a subset of patients with placebo and 81 mg
attach(aspirinc)
subset1  <- DOSE == 0 | DOSE == 81
aspirin1 <- aspirinc[subset1,]
detach(aspirinc)

# set up treatment assignment
trtA     <- rep(0, dim(aspirin1)[1])
trtA[aspirin1[,"DOSE"] == 81] <- 1

# STEPP analysis A: placebo vs 81 mg aspirin
inc_win     <- stepp.win(type="sliding", r1=30, r2=100)
inc_sp      <- stepp.subpop(swin=inc_win, cov=aspirin1$AGE)

ADorLE      <- as.numeric(aspirin1$AD==1 | aspirin1$AL==1)
modelA      <- stepp.GLM(colY = ADorLE, glm = "binomial", link = "logit")
# Warning: In this example, the permutations have been set to 50 to allow the function
# to finish in a short amount of time.  IT IS RECOMMEND TO USE AT LEAST 2500 PERMUTATIONS TO 
# PROVIDE STABLE RESULTS.
statGLM  <- new("steppes")		  # create a test object based on subpopulation and window
statGLM  <- estimate(statGLM, inc_sp, modelA) # estimate the subpopulation results
# Warning: In this example, the permutations have been set to 0 to allow the function
# to finish in a short amount of time.  IT IS RECOMMEND TO USE AT LEAST 2500 PERMUTATIONS TO 
# PROVIDE STABLE RESULTS.
statGLM  <- test(statGLM, nperm=10)       # permutation test with 0 iterations

print(statGLM)				  # print the estimates and test statistics
plot(statGLM, ncex=0.70, legendy=30, pline=-4.5,
     xlabel="Subpopulations by Median Age", ylabel="Risk",
     nlas=3, noyscale=TRUE, subplot=TRUE)

### Example for balanced subpopulations ###
library(stepp)

data(bigKM)

rxgroup <- bigKM$trt
time    <- bigKM$time
evt     <- bigKM$event
cov     <- bigKM$ki67

res <- balance_patients(range_r1 = c(30, 100), range_r2 = c(50, 300),
                        maxnsubpops = 100, covar = cov, verbose = TRUE,
                        plot = TRUE, contour = FALSE)
