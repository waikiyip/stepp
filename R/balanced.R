# This looks for the combination of r1 and r2 values that
# produces the most balanced subpopulations. By most balanced 
# we mean that the sizes of the subpopulations are the most
# similar as measured bu their variance.
# The input data is the vector of covariate values, and the
# two ranges of values for r1 and r2.
# Marco Bonetti May 2020
# r2 > r1

ranger2 <- c(600,900)
ranger1 <- c(100,400)
maxnsubpops <- 50

comprisk <- scan(file="~/Desktop/RESEARCH/STEPP BALANCED SUBPOPULATIONS/comprisk.aphinity.csv")

summary(comprisk)
hist(comprisk)
n <- length(comprisk)
n

freqdist <- data.frame(table(comprisk))
freqdist[,1] <- as.numeric(as.character(freqdist[,1]))
freqdist[,3] <- NA
names(freqdist)[3] <- "CumFreq"
k <- length(freqdist[,2])
freqdist[1,3] <- freqdist[1,2]
for(i in 2:k)
freqdist[i,3] <- freqdist[i-1,3]+ freqdist[i,2]

# freq [a,b] for all combinations - do only once!
allfreqs <- matrix(rep(NA,k*k),nrow=k)
allfreqs[1,] <- freqdist$CumFreq
for(i in 2:k)
allfreqs[i,i:k] <- freqdist$CumFreq[i:k]-freqdist$CumFreq[i-1]


unbalance <- function(pars,printind)
{
  r1 <- pars[1]
  r2 <- pars[2]
   # We work with indices, not values
  # Note: indices of values are from 1 to k
  subpops <- matrix(rep(NA,5*maxnsubpops),ncol=5)
  colnames(subpops) <- c("IndLow","IndHigh","Size","CutoffLow","CutoffHigh")
  subpops[1,1] <- 1
  # collect at least r2 subjects
  lasttemp <- sum(allfreqs[1,] < r2)+1
  subpops[1,2] <- lasttemp
  subpops[1,3] <- allfreqs[1,lasttemp]
  # drop at least (r2-r1) subjects
  for(i in 2:maxnsubpops)
  {
    preva <- subpops[i-1,1]
    prevb <- subpops[i-1,2]
    # First we find the potential new a 
    a <- prevb - sum(allfreqs[preva:prevb,prevb]<=r1)+1   # Careful: this could be >= prevb (CHECK)
    if(prevb<k)
      {
      if(allfreqs[a,k] >=r2)  # if there is enough to build a new subpopulation
      {
      subpops[i,1] <- a
      b <- a + sum(allfreqs[a,a:k]<r2) 
      subpops[i,2] <- b
      subpops[i,3] <- allfreqs[a,b]
      }  else {
      # last subpopulation
      subpops[i,1] <- a
      b <- k
      subpops[i,2] <- b
      subpops[i,3] <- allfreqs[a,b]
      break
      }
      } else {i <- i-1; break}
  }
nsubpops <- i
subpopsfin <- subpops[1:nsubpops,]
# Compute the unbalance measure
res <- var(subpopsfin[,3])
# Add the covariate values that define the subpopulations
subpopsfin[,4] <- freqdist[subpopsfin[,1],1]
subpopsfin[,5] <- freqdist[subpopsfin[,2],1]
if(printind == 1) {print(subpopsfin);
print("Variance of subpopulations sizes and number of subpopulations:")}
return(c(res,nsubpops))
}

resunb <- unbalance(c(250,800),1)
resunb


bestr1 <- ranger1[1]
bestr2 <- ranger2[1]  
resmat <- matrix(rep(NA,5*(ranger1[2]-ranger1[1]+1)*(ranger2[2]-ranger2[1]+1)),ncol=5)
print(paste("Range for r1 =",ranger1[1],ranger1[2]))
print(paste("Range for r2 =",ranger2[1],ranger2[2]))
minvar <- 10000000  
nextind <- 1
for(i in ranger1[1]:ranger1[2])
  for(j in ranger2[1]:ranger2[2])
  {
   resunb <- unbalance(c(i,j),0)
   newvar <- resunb[1]
   #   indx <- ranger1[1]-1+i
#   indy <- ranger2[1]-1+j
   resmat[nextind,1] <- i
   resmat[nextind,2] <- j
   resmat[nextind,3] <- newvar
   resmat[nextind,4] <- resunb[2]
#   print(c(i,j,newvar))
   if(newvar < minvar)
     {
     bestr1 <- i
     bestr2 <- j
     minvar <- newvar
     }
nextind <- nextind + 1
     }

# Add the cutoff values for the subpopulations
soln <- c(bestr1,bestr2,minvar)
names(soln) <- c("r1-Opt","r2-Opt","Variance")
print(soln)
unbalance(c(bestr1,bestr2),1)

# Check aphinity stepp
#unbalance(c(250,800),1)

# Add the shade of gray to the last column: darker means smaller variance
# work on log scael to highlight best results only (darker)
# qs <- quantile(logresmat3,c((0:8)/8))
logresmat3 <- log(resmat[,3])
cutoffs <- min(logresmat3) + ((0:8)/8)*(max(logresmat3)-min(logresmat3))
for(i in 1:length(logresmat3))
resmat[i,5] <- (sum(cutoffs <= logresmat3[i]) -1)/8

# there are 8 shades of gray from 0/8 to 7/8 (8/8 is white)
pdf(file="allcombs.pdf", paper="a4r")
plot(ranger1,ranger2,type="n",xlab="r1",ylab="r2",
     main="Subpopulations (labels) and variances (shades of gray, dark is small)")
cexval <- 0.20
text(resmat[,1],resmat[,2],labels=as.character(resmat[,4]),cex=cexval,col=gray(resmat[,5]))
points(bestr1,bestr2,pch=1,col="red")
dev.off()



