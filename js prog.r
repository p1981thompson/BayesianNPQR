  #
# Functions to simulate from a distribution
# comprising a normal body and a Pareto tail
#
# Function originally written by Paul Thompson,
# then modified by Julian Stander
#
#--------------------------------------------------------------
#
# These functions are used in the simulation function
#

#
# Function to work out the area below u for a normal distribution
# truncated at zero
# with mean gamma and standard deviation alpha.
# A value beta is subtracted for use in get.mean.below
#
area.below <- function(gamma, alpha, u , beta){
(pnorm(u, gamma, alpha) - pnorm(0, gamma, alpha)) / (1 - pnorm(0, gamma, alpha)) - beta
}

#
# Given beta, alpha and u, the mean gamma of the above truncated normal
# is determined.  This function calculates it.
#
get.mean.below <- function(beta, alpha, u){
# beta is proportion below
# alpha is standard deviation below
# u is threshold
uniroot(area.below,c(0,u), beta = beta, alpha = alpha, u = u)$root
}

#
# Example of its use
#
get.mean.below(0.9, 0.7, 2.9)


# Function to simulated from a normal body Pareto tail distribution
# with continuous density function (see Thompson, Cai, Reeve and Stander)

# q is the proportion of the density that is below the threshold
# Above q is beta and threshold is u
# sigma.below is the standard deviation of the truncated normal body
# Above sigma.below is gamma
#
# n is the number of realizations required

r.norm.pareto <- function(n=10000,  q = 0.9, threshold, sigma.below = 0.7, xi = 0.2, do.plot = FALSE)
{
#
# Work out mean.below (alpha above)
#
mean.below <- get.mean.below(q, sigma.below, threshold)

sigma.gpd <- (1-q) *  pnorm(0, mean.below, sigma.below, lower.tail = FALSE) / dnorm(threshold, mean.below, sigma.below)

# Vector for realizations

newdist<- vector(mode = "numeric", length = n)
#
# Now fill up the vector according to the algorithm in the paper
#
for(i in 1:n)
{
  newdist[i]<-rnorm(1,mean.below,sigma.below)
  if(newdist[i]<=0)
      {
      while(newdist[i]<=0)
        {
        newdist[i]<-rnorm(1,mean.below,sigma.below)
          if(newdist[i]>threshold)
            {
            newdist[i]<-rgpd(1,xi=xi,mu=threshold,beta=sigma.gpd)
            }
        }
      }
  if(newdist[i]>threshold)
      {
       newdist[i]<-rgpd(1,xi=xi,mu=threshold,beta=sigma.gpd)
      }
}
#
# Plot
#
if(do.plot){
par(mfrow = c(2,1))
hist(newdist,nclass = 100, xlab = "Realizations", main = "Histogram of realizations from a truncated normal body and a Pareto tail")
#
title(sub = substitute(paste(u == U, ", ", beta == Beta, ", ", gamma == Gamma, ", ", alpha == Alpha, ", ", xi == Xi, ", ", sigma == Sigma),list(U = signif(threshold,3), Beta = signif(q,3), Gamma = signif(mean.below,3), Alpha = signif(sigma.below,3),  Xi = signif(xi,3), Sigma = signif(sigma.gpd,3))))
#
# Rug
#
rug(newdist)
#
# Show the threshold
#
abline(v = threshold, col = "red", lwd = 2)
#
plot(newdist)
abline(h = threshold, col = "red", lwd = 2)

}
#
#
if(do.plot){
print(paste("Proportion below threshold", length(newdist[newdist < threshold])/length(newdist)))
}
#
# Return
#
return(invisible(list(data = newdist, mean.below = mean.below, sigma.gpd = sigma.gpd)))
}
#--------------------------------------------------------------------------
#
#
#--------------------------------------------------------------------------
#
# Function to estimate a 95% return level confidence interval from a GPD fit
#
gpd.ci <- function (z,  return.year = 100)
{
a <- z$mle
u <- z$threshold
la <- z$rate
n <- z$n
npy <- z$npy
mat <- z$cov
#
# a is the mle of sigma and xi
# u is the threshold
# la is the  mle estimate of the exceedance probability
# n is the number of data points
# npy is the number of observations per year, set by default to 365 in gpd.fit
# mat is the variance-covariance matrix for sigma and xi
# return.year is the return period in years
#
# m is the number of observations, i.e. the return period in years, multiplied by
# the number of observations per year
#
    m <- npy * return.year
#
    a <- c(la, a)
    eps <- 1e-06
    a1 <- a
    a2 <- a
    a3 <- a
    a1[1] <- a[1] + eps
    a2[2] <- a[2] + eps
    a3[3] <- a[3] + eps
    q <- gpdq2(a[2:3], u, la, m)
    d1 <- (gpdq2(a1[2:3], u, la, m) - q)/eps
    d2 <- (gpdq2(a2[2:3], u, la, m) - q)/eps
    d3 <- (gpdq2(a3[2:3], u, la, m) - q)/eps
    d <- cbind(d1, d2, d3)
    mat <- matrix(c((la * (1 - la))/n, 0, 0, 0, mat[1, 1], mat[1,
        2], 0, mat[2, 1], mat[2, 2]), nc = 3)
    v <- apply(d, 1, q.form, m = mat)
    return(list(return.level = q, lower.limit = q - 1.96 * sqrt(v), upper.limit = q + 1.96 * sqrt(v)))
}
#--------------------------------------------------------------------------

  ###THRESHOLD SELECTION FUNCTION###

     ppfitrange_idea5<-function (data, umin=0, umax, nint = 100,show=F)
{
  library(nortest)
lowfindthresh<-function (data, ne)
{
    data <- (sort(as.numeric(data)))
    thresholds <- unique(data)
    indices <- match(data[ne], thresholds)
    indices <- pmin(indices + 1, length(thresholds))
    thresholds[indices]
}

umax<-findthresh(data,50)
#print(umax)
umin<-lowfindthresh(data,50)
#print(umin)

    	m <- s <- paulvar <- up <- ul <- matrix(NA, nrow = nint, ncol = 2)
    	my.m <- m

  u <- seq(umin, umax, length = nint)
  #print(u)
  for (i in 1:nint) {
        z <- paulgpd5(data, u[i])
	#MAX LIKELIHOOD EST FOR EACH MODEL.
        m[i, ] <- z$mle
	  paulvar[i, ] <- diag(z$analytic_covar)

	#REPARAMETERIZATION OF GPD SCALE PARAMETER.
     #   m[i, 1] <- m[i, 1] - m[i, 2] * u[i]
        d <- matrix(c(1, -u[i]), ncol = 1)
      #VARIANCE AND STD ERRORS FOR EACH MODEL.
	 # v <- t(d) %*% z$analytic_covar %*% d
#        s[i, ] <- z$se2
#	  s[i, 1] <- sqrt(v)
	#95% CONF INTERVAL FOR MLE.
   #     up[i, ] <- m[i, ] + 1.96 * s[i, ]
   #     ul[i, ] <- m[i, ] - 1.96 * s[i, ]


diff<-matrix(NA,nrow = nint-1, ncol = 2)
diff2<-matrix(NA,nrow = nint-1, ncol = 2)

test<-vector(mode="numeric", length = (nint-1))


   names <- c("Differenced Modified Scale", "Differenced Shape")
   # oldpar <- par(mfrow = c(2, 1))

        		um <- max(up[, 2])
     			ud <- min(ul[, 2])
			my.m[,2] <- m[,2]
			#par(oldpar)
			invisible()
			}

for(g in 1:(length(m[,1])-1))
		{
		diff[g,]<-m[(g+1),]-m[g,]
		diff2[g,]<-paulvar[(g+1),]-paulvar[g,]
		}
p0old<-diff[,1]
uold<-u[2:length(u)]
#m<-na.omit(m)
#paulvar<-na.omit(paulvar)
diff.na<-is.na(diff[,1])
diff.na2<-is.na(diff2[,1])

## diff ##
#p0 <- diff[!diff.na,2]
#if(any(is.na(p0)))
#{print(p0)}
#u <- u[!diff.na]

#if(is.na(m[1,]) || is.na(paulvar[1,]))
#{na.omit(m[1,])
#na.omit(paulvar[1,])}

#m <- m[-1,]
#paulvar <- paulvar[-1,]

#m <- m[!diff.na,2]
#paulvar <- paulvar[!diff.na,2]

####diff2####

p0 <- diff[!diff.na2,1]
#if(any(is.na(p0)))
#{print(p0)}
u <- u[!diff.na2]

#if(is.na(m[1,]) || is.na(paulvar[1,]))
#{na.omit(m[1,])
#na.omit(paulvar[1,])}

m <- m[-1,]
paulvar <- paulvar[-1,]

m <- m[!diff.na2,1]
paulvar <- paulvar[!diff.na2,1]

####

#cat("m dimension",dim(m),"\n")
#cat("paulvar dimension",dim(paulvar),"\n")
#print(paulvar[,2])

pointset<- length(p0)
#print(pointset)
 quantilemin<-quantile(data, prob=0.98)
 newquant<-paulgpd5(data,quantilemin)
  #print(newquant$num_exceed)
 testex <- 0.25*newquant$num_exceed

loopy1<-function(pointset=pointset,u=u, m=m,paulvar=paulvar,umin=umin,umax=umax)
{
a<-list()
#while(testex <= 100)
# {

    testcase <- function(pointset=pointset,u=u, m=m,paulvar=paulvar,umin=umin,umax=umax)
        {
          z<-list()
          test<-vector(mode="numeric", length = (pointset))
          my.test1<-vector(mode="numeric", length = (pointset))
          my.test2<-vector(mode="numeric", length = (pointset))
          my.test3<-vector(mode="numeric", length = (pointset))
          my.test4<-vector(mode="numeric", length = (pointset))
          my.test5<-vector(mode="numeric", length = (pointset))
          testend<-length(test)-7

        for(ppit in 1:testend)
				{ #cat("ppit, iteration=",ppit,"\n")
          z$thres<-NA
        # print(paste("mean part=",m[ppit]))
        #  print(paste("variance part=",paulvar[ppit]*((u[ppit+1]-u[ppit])^2)))


         # my.test1[ppit] <- pearson.test(p0[ppit:length(test)])$p.value
          my.test2[ppit] <- sf.test(p0[ppit:length(test)])$p.value
          my.test3[ppit] <- ad.test(p0[ppit:length(test)])$p.value
          my.test4[ppit] <- lillie.test(p0[ppit:length(test)])$p.value
          #my.test5[ppit] <- cvm.test(p0[ppit:length(test)])$p.value

          #test[ppit]<-mypearson2(p0[1:ppit],meanpear = (m[ppit]*(umax-umin))/length(m), sdpear = paulvar[ppit]/(length(m)^2))
          test[ppit]<-mypearson2(p0[ppit:length(test)],meanpear = m[ppit], sdpear = paulvar[ppit]*((u[ppit+1]-u[ppit])^2))
          #print(paste("my test value1=",my.test1[ppit]))
         ##print(paste("my test value2=",my.test2[ppit]))
          ##print(paste("my test value3=",my.test3[ppit]))
          #print(paste("my test value4=",my.test4[ppit]))
          #print(paste("my test value5=",my.test5[ppit]))
          #cat("test value=",test[ppit],"\n")
          ##cat("threshold",u[ppit],"\n")
          #cat(" ","\n")

         if(my.test4[ppit] > 0.05)
					     {
                #if(my.test4[ppit-1] > 0.05)
					      #  {
                 #   if(my.test4[ppit-2] > 0.05)
	 				        #   {
	                     z$thres<-u[ppit]
                       z$b<-ppit
                   #  }
                  #}
					     }
					if(is.numeric(z$thres))
				  {
          #hist( p0[ppit:length(test)],10)
          #meanmine<-mean(p0[ppit:length(test)])
          #print(m[ppit])
         #abline(v=ppit,lwd=2,col="blue",lty=2)
          break
          }

        }
       # hist(my.test4,20)
        #plot(seq(c(1:length(my.test4))),my.test4,type="l",lty=1,xlab="normality test number",ylab="P values for normality tests")
        #lines(seq(c(1:length(my.test4))),my.test3,lty=2)
        #lines(seq(c(1:length(my.test4))),my.test2,lty=3)
        #abline(v=89,lwd=2,col="blue",lty=3)

  #      leg3.txt<-c("Shapiro-Francia","Lilliefors","Anderson-Darling")
  # legend(x=10,y=0.1,leg3.txt,cex=0.7,lty=c(3,2,1),lwd=1.5,bg="white")


         #print(z$thres)
	      z$thres<-z$thres
				z$b<-z$b
				invisible(z)
       }

       mytest<-testcase(pointset,u,m,paulvar,umin,umax)
    #   cat("b value=",mytest$b,"\n")
   #    cat("test thres=",mytest$thres,"\n")
      newtest<-paulgpd5(data,mytest$thres)
      testex<-newtest$num_exceed
  #     cat("test exceed=",testex,"\n")
      if(testex <= max(newquant$num_exceed))
        {pointset<-mytest$b}
 #     cat("pointset=",pointset,"\n")
    #rnd1<-round(mytest$thres,digits=3)
    a$choice<-mytest$thres
# }

 a$choice<-a$choice
 invisible(a)
}
selecta<-loopy1(pointset,u,m,paulvar,umin,umax)
    ulgth<-length(u)
    u<-u[2:ulgth]
#plot(uold, p0old,xlab = expression(paste("Threshold ", u[j-1])), ylab = expression(hat(tau)[u[j]] - hat(tau)[u[j-1]]), type = "b",ylim=c(-0.2,0.2),cex.lab=1.3, cex.axis=1.3)
#abline(v=1.672323,col="red",lwd=2)
#	abline(v=selecta$choice,col=10,lwd=2)
#  abline(v=umin,lty=3)
#  abline(v=umax,lty=3)
 # print(testcase$thres)
#cat("threshold choice=",selecta$choice)
#z<-list()
threshold<-selecta$choice
 #print(z$threshold)
#if (show) {
#           print(z[1])
#          }
invisible(z)
 #return(list(threshold=threshold))
 return(threshold)
 }
 
 #--------------------------------------------------------------------------
 
# We need the library evir for generating from the Pareto tail

library(evir)
library(ismev)

# n, Proportion, threshold, sigma.below, shape

# Paper's notation

beta <- 0.9
u <- 2.9
alpha <- 0.7
xi <- 0.2

sim.data <- r.norm.pareto(365*27, beta, u,  alpha, xi, do.plot = TRUE)

# r.norm.pareto also calculate the paper's sigma (and gamma)
 print(dim(sim.data))
sigma <- sim.data$sigma.gpd
#
ny <- 365
#
# zeta is probability of being in the tail
#
zeta <- 1 - beta
#
# *** 100 year return level using formula in Coles ***
#
N <- 100
true.return.level <- u + (sigma / xi)*((N*ny*zeta)^xi - 1)
true.return.level

#
#
# Check coverage of confidence intervals
#
# Set count to zero, generate n.rep data sets and work out and check confidence interval for each one
#
i.cover <- 0
n.rep <- 1
#
for(i in 1:n.rep){

if(i %% 10 == 0){print(paste("Iteration", i))}

#
sim.data <- r.norm.pareto(365*27, beta, u,  alpha, xi, do.plot = FALSE)
#
#
#********

 x<-sim.data$data
  bboot<-vector(mode="numeric",length=100)
 for(i in 1:100)
 {
  bboot[i]<-ppfitrange_idea5(sample(x,9855,replace=T))
 # boot[i]<-boot1$threshold
 }
  sortboot3<-sort(bboot)

   ###EXTRACT CONFIDENCE INTERVALS AND MEDIAN###

   boot.thres <- median(sortboot3)

   boot.thres.up <- quantile(sortboot3,0.975)

   boot.thres.low <- quantile(sortboot3,0.025)

   #hist(boot,main="Bootstrap Interval for Threshold Uncertainty",xlab="Threshold (m)")

   #abline(v=boot.thres)
   #abline(v=boot.thres.up,col="blue", lty=2)
  # abline(v=boot.thres.low,col="blue", lty=2)
   #abline(v=1.607,lty=3)
#********
## ***** PAUL:  THESE LINES OF CODE COULD BE REPLACED WITH YOUR BOOTSTRAP CONFIDENCE INTERVAL GENERATING CODE
##
## Fit model
##
#my.fit <- gpd.fit(sim.data$data, u, show = FALSE)
##
## Estimate the confidence interval
##
#ci.hat <- gpd.ci(my.fit, N)
##
## ***** PAUL:  TO HERE
#

ci.hat<-list(lower.limit=boot.thres.low,upper.limit=boot.thres.up)
print(ci.hat$lower.limit)
print(ci.hat$upper.limit)
#
# Check whether it contains the true value
#
i.cover <- i.cover + ifelse(ci.hat$lower.limit < true.return.level & true.return.level < ci.hat$upper.limit, 1, 0)
}
#
# Proportion covering
#
print(paste("Proportion of intervals covering", i.cover / n.rep))


