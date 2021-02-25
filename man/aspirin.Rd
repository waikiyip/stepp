\name{aspirin}
\alias{aspirin}
\docType{data}
\title{
  The aspirin data set.
}

\description{
  The "aspirin" data set was obtained with permission from the Polyp Prevention Study Group. The data are from a randomized clinical trial
  comparing two dosages of daily aspirin (81 mg or 325 mg/day) versus placebo for the prevention of colorectal adenomas among people having 
  a prior history of these lesions.  The study randomized a total of 1,121 participants to the three treatment groups.  After an approximately 
  three-year treatment period, 1,084 participants underwent follow-up colonoscopy. The primary endpoint of the trial was the occurrence of
  one or more large-bowel adenomas.  Secondary endpoints included the occurrence of any advanced lesion (i.e., adenomas having high-risk 
  characteristics or invasive cancer).
}

\usage{data(aspirin)}
\format{
  The data set has the following columns:
    DOSE - Randomized aspirin dose (0=placebo);
    AGE  - Age of the participant in years at baseline;
    G    - Gender of the participant (M or F)
    AD   - Occurrence of any adenoma at follow-up (0=no, 1=yes, .=missing);
    AL   - Occurrence of any advanced lesion at follow-up(0=no, 1=yes, .=missing)
}

\source{
  The Polyp Prevention Study Group is acknowledged for granting permission to use the data from 
  the Aspirin/Folate Polyp Prevention study.
}

\references{
  Baron J, Cole BF, Sandler RS, et al (2003), "A Randomized Trial of Aspirin to prevent Colorectal Adenomas." n engl j med 2003;234:891-9.
}

\examples{
data(aspirin)

# remove cases with missing data
aspirinc <- aspirin[complete.cases(aspirin),]

# make a subset of patients with placebo and 81 mg
attach(aspirinc)
subset1  <- DOSE == 0 | DOSE == 81
aspirin1 <- aspirinc[subset1,]
detach(aspirinc)

# set up treatment assignment
trtA     <- rep(0, dim(aspirin1)[1])
trtA[aspirin1[,"DOSE"] == 81] <- 1

# STEPP analysis A: placebo vs 81 mg aspirin
attach(aspirin1)
inc_win     <- stepp.win(type="sliding", r1=30, r2=100)
inc_sp      <- stepp.subpop(swin=inc_win, cov=AGE)

ADorLE      <- as.numeric(AD==1 | AL==1)
modelA      <- stepp.GLM(coltrt=trtA, trts=c(0,1), colY=ADorLE, glm="binomial")
# Warning: In this example, the permutations have been set to 50 to allow the function
# to finish in a short amount of time.  IT IS RECOMMEND TO USE AT LEAST 2500 PERMUTATIONS TO 
# PROVIDE STABLE RESULTS.
steppGLMA   <- stepp.test(inc_sp, modelA, nperm=50)
summary(steppGLMA)
print(steppGLMA)
plot(steppGLMA, ncex=0.70, legendy=30,
	pline=-4.5, color=c("red","black"),
	xlabel="Subpopulations by Median Age", ylabel="Risk",
	tlegend=c("Placebo", "81 mg aspirin"), nlas=3, pointwise=FALSE, noyscale=TRUE, rug=FALSE)

detach(aspirin1)
}

\keyword{datasets}

