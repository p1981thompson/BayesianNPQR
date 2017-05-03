
        ##########################################################
        #                                                        #
        # 3D RETURN LEVEL PLOTS FOR QUANTILE REGRESSION SPLINES  #
        #                                                        #
        ##########################################################



      ##############################################
  ####LOAD AND SORT DATA FOR CORRECT QUANTILE CURVES####
      ##############################################

#new.dat<-read.table("F:/University work/R work/Rana QR/mysampsmlalt.txt",header=T)

#new.dat<-read.table("d://Local Data//p1thompson//University work//R work/blank//azdatgui.txt",header=T)
#
# Sort out the data
#
new.dat<-serum.dat
new.data<-data.frame(new.dat[,2],new.dat[,4])
#
# new.data<-new.data[ do.call(order, new.data) ,]
#
new.data <- new.data[order(new.data[,4]), ]

new.data<-as.matrix(new.data)

#
# Check by plotting
#
x <- new.data[,1]
y <- new.data[,2]


#x <- rnorm(5000)
#y <- rnorm(5000)
#y <- y[order(x)]
#x <- sort(x)

   #### SCATTER PLOT OF DATA #####

plot(x, y, xlab = "Age (years)", ylab = "IgG (grams/litre)",cex=0.5,pch=20,col="skyblue", cex.lab=1.3, cex.axis=1.3)

   #### LOAD REQUIRED LIBRARIES FOR PROGRAM####

library(rgl)
library(quantreg)
library(splines)
library(MCMCpack)
      ############################################
  #### MAIN QUANTILE REGESSION SPLINES  ALGORITHM ####
      ############################################
# --------------------------------------------------
	prob <- function(x, y, x.grid, g, p, K, lambda)
	{
    spline.fun.g <- interpSpline(x.grid, g)
       spline.pred<-predict(spline.fun.g,x)
   	u <- (y - spline.pred$y)

		ind <- ifelse(u < 0, 1, 0)
		z <-  -0.5*lambda*(t(g)%*%K%*%g) - sum(u * (p - ind))
		return(z)
	}
# --------------------------------------------------

# --------------------------------------------------
	prob2 <- function(x, y, x.grid, g, p, K, lambda, beta, alfa)
	{
	#print(beta)
#	print(alfa)
	#print(alfa*beta)
	#print(alfa*beta^2)
    n<-length(x.grid)
		z <-  ((n-2)/2+alfa)*log(lambda)-lambda*(0.5*(t(g)%*%K%*%g)+1/beta)
		return(z)
	}
# --------------------------------------------------



# --------------------------------------------------
	accept <- function(x, y, x.grid, gold, gnew, p,   K, lambda, denom)
	{

		#denom <- prob(y, x, p, thold)
		num <- prob(x, y, x.grid, gnew, p, K, lambda)
		mh <- (num - denom)
		mh <- exp(mh)
		mh <- min(mh, 1)
		q <- runif(1, 0, 1)
		if(mh > q) {
			th <- gnew
			den <- num
			i.accept <- 1
		}
		else {
			th <- gold
			den <- denom
			i.accept <- 0
		}

		return(list(th = th, den = den, i.accept = i.accept))
	}
# --------------------------------------------------

# --------------------------------------------------
	accept2 <- function(x, y, x.grid, gnew, p,   K, lambdaold,lambdanew, denom,beta,alfa)
	{

		#denom <- prob(y, x, p, thold)
		num <- prob2(x, y, x.grid, gnew, p, K, lambdanew,beta,alfa)
		mh <- (num - denom)
		mh <- exp(mh)
		mh <- min(mh, 1)
		q <- runif(1, 0, 1)
		if(mh > q) {
			lam <- lambdanew
			den2 <- num
			i.accept2 <- 1
		}
		else {
			lam <- lambdaold
			den2 <- denom
			i.accept2 <- 0
		}

		return(list(i.accept2 = i.accept2, lam = lam, den2 = den2))
	}
# --------------------------------------------------

# --------------------------------------------------
make.RQ <- function(unq.x)
{

Hlgth <- length(unq.x)

H <- diff(unq.x)
#print(H)

  R <- matrix( 0,Hlgth,Hlgth)
  Q <- matrix( 0,Hlgth,Hlgth)

         for(j in 2:(Hlgth-1))
                  {

                  Q[j-1,j] <- H[j-1]^(-1)
                  Q[j,j] <- -H[j-1]^(-1) - H[j]^(-1)
                  Q[j+1,j] <- H[j]^(-1)
                  R[j,j] <- (1/3)*(H[j-1] + H[j])

                  }

  #print(Q)

                  for(d in 2:(Hlgth-2))
                        {
                        R[d,d+1] <- (1/6)*H[d]
                        R[(d+1),d] <- (1/6)*H[d]
                        }



    Q <- Q[,-(Hlgth)]

    Q <- Q[,-1]
#print(Q)

   R <- R[-(Hlgth),]

   R <- R[,-(Hlgth)]

   R <- R[-1,-1]

   return(list(R = R, Q = Q))

}
# --------------------------------------------------

  ### GRID LENGTH SET AS DEFAULT=20 BUT CAN BE ALTERED MANUALLY IN FUNCTIONS BELOW ###

ky3alt<-function(x,y,p,  r,sigma, grid.length = 20,sigma2)
{
#
# This function requires the x's to be monotonically increasing
#


    unq.x <- unique(x)
    if(grid.length > length(unq.x)){
       warning("Finer grid than data")
    }

#    length.unq.x <- length(unq.x)
#    int.interval <- length.unq.x %/% grid.length
#    if(int.interval == 0) int.interval <- 1

#    ind.uni.x.grid <- seq(from = 1, to = length.unq.x, by = int.interval)
     #print(ind.uni.x.grid)
#    x.grid <- unique(c(unq.x[ind.uni.x.grid], max(unq.x)))
     #print(x.grid)
 #   grid.length <- length(x.grid)

 x.grid <- seq(from = min(x), to = max(x), length = grid.length )
    #print(x.grid)


    Hlgth<-length(x.grid)



    my.RQ <- make.RQ(x.grid)
    R <- my.RQ$R
    Q <- my.RQ$Q

    #print(R)
   # print(Q)

Rinv <- solve(R)

K <- Q%*%Rinv%*%t(Q)
#
# Singular value decomposition of K and generalized inverse
#
s.K <- svd(K)

s.K$d <- c(s.K$d[1:(grid.length - 2)], 0, 0)
d.ginv <- ifelse(s.K$d > 0 , 1 / s.K$d, 0)
K.ginv.half <- s.K$u  %*% sqrt(diag(d.ginv)) %*% t(s.K$v)

  #print("Rinv=")
  #print(Rinv)

  #print("K=")
  #print(K)

	G <- matrix(0,r,grid.length)
	Lambda <- vector(mode="numeric",length=r)
 # X<- model.matrix(y~bs(x,df=10))

 ####INITIAL INPUT INTO THE MCMC ROUTINE, df CAN BE ALTERED TO ALLOW FASTER CONVERGENCE####
 beta<-(0.1/smooth.spline(x,y)$lambda)
 alfa<-(smooth.spline(x,y)$lambda/beta)

qmodel <- rq(y ~ bs(x,df=15), tau = p)
y.hat <- predict(qmodel, newdata = data.frame(x = x.grid))

#qmodel <- rq(y ~ x + I(x^2) + I(x^3), tau = p)
#y.hat <- predict(qmodel, newdata = data.frame(x = x.grid))
#lines(x.grid, y.hat, col="red", lwd = 2)

 #y.hat <-rep(0.5,length(x.grid))
#y.hat <- rnorm(length(x.grid))

     gold <- y.hat

#     lines(x.grid, gold, col="green", lwd = 2)

    # print(gold)

     G[1,] <- gold
     # G[1,]<-rep(1,length(x.grid))
     # lambdaold <- smooth.spline(x,y)$lambda
       lambdaold <- 0.003
      Lambda[1] <- lambdaold
	 #  denom <- prob(x, y, x.grid, gold, p, K, lambda)

    denom <- prob(x, y, x.grid, gold, p, K, lambdaold)


	    no.acc <- 0
	    no.acc2 <- 0

	for(i in 2:r) {

	if(i %% 100 == 0){ print(paste("Iteration",i))}



#          	gnew <- gold + rnorm(grid.length, 0, sigma)

Z <- rnorm(grid.length)
V <- sigma*(K.ginv.half / sqrt(lambdaold))
eta <- as.vector(V %*% Z)


gnew <- gold + eta


            acc <- accept(x, y, x.grid, gold, gnew, p, K, lambdaold, denom)

            gold <- acc$th
            denom2 <-  prob2(x, y, x.grid, gold, p, K, lambdaold,beta,alfa)

            lambdanew <- rnorm(1,log(lambdaold),sigma2)
            lambdanew <- exp(lambdanew)

            acc2 <- accept2(x, y, x.grid, gold, p, K, lambdaold, lambdanew, denom2,beta,alfa)
            lambdaold<-acc2$lam
            denom <-  prob(x, y, x.grid, gold, p, K, lambdaold)
            G[i,]<-gold
            Lambda[i]<-lambdaold
            no.acc <- no.acc + acc$i.accept
            no.acc2 <- no.acc2 + acc2$i.accept2
                 }


                 print(paste("Proportion accepted", no.acc /(r - 1)))
                 print(paste("Proportion lambda accepted", no.acc2 /(r - 1)))
 return(list(G = G, x.grid = x.grid,K = K,R=R,Q=Q,Lambda = Lambda))
 }


    z1 <- ky3alt(x, y, 0.9 ,  250000, 0.35, grid.length = 20,1.4)

   thin.seq<-seq(from=1,to=250000,by=250)

  #par(mfrow=c(3,2))

   newG<-matrix(0,length(thin.seq),20)
  for (i in 1:20)
   {
      newG[,i]<-z1$G[thin.seq,i]
   }
   
	z1new <- newG[81:1000,]
  xr <- c(-1, 1)
	yr <- c(0, 3)
	r <- 920
	h <- 0.025
	hl <- floor(h * r)
	hu <- ceiling((1 - h) * r)
	ww1 <- apply(z1new, 2, sort)
	dim(ww)
	xl1 <- ww1[hl,  ]
	xu1 <- ww1[hu,  ]
	xm1 <- apply(ww1, 2, mean)
  xp1<-sort(z1$x.grid)

#  Lam1new <- z1$Lambda[5001:50000]
#  all<-as.matrix(data.frame(Lam1new,z1new))
# 	ww1 <- apply(all, 1,sort)
#  xl1 <- ww1[hl[2:length(x.grid)],  ]
#	xu1 <- ww1[hu[2:length(x.grid)],  ]
#	xm1 <- apply(ww1, 1, mean)
#	xmed1<-apply(ww1, 1, median)
#  xp1<-sort(z1$x.grid)

  ##INTERPOLATING SPLINES FITTED THROUGH MCMC POINTS TO SMOOTH##

  ispl1<-interpSpline(xp1,xm1)
  ispl1new<-spline(xp1,xm1,n=50,method="natural")
 # ispl1newmed<-spline(xp1,xmed1,n=1000,method="natural")
  ispl1up<-spline(xp1,xu1,method="natural")
  ispl1low<- spline(xp1,xl1,method="natural")

  lines(ispl1new, col="blue",lty=1,lwd=1.5)
 # lines(ispl1newmed, col="orange",lty=1)
   lines(ispl1up, col="blue",lty=2,lwd=1.5)
   lines(ispl1low, col="blue",lty=2,lwd=1.5)




   par(mfrow=c(3,2))
  for (i in 1:50)
   {
    ts.plot(z1$G[,i])
   }

   par(mfrow=c(1,1))

    ts.plot(z1$Lambda)
    abline(h = smooth.spline(x,y)$lambda,col="red")


    thin.seq<-seq(from=1,to=250000,by=100)

  par(mfrow=c(3,2))

   newG<-matrix(0,length(thin.seq),50)
  for (i in 1:50)
   {
      newG[,i]<-z1$G[thin.seq,i]
   }

   par(mfrow=c(3,2))
  for (i in 1:50)
   {
    ts.plot(newG[,i],ylab="MCMC sub-chain of g",xlab="Iteration")
   }



      newlambda<-z1$Lambda[thin.seq]
ts.plot(newlambda,ylab=expression(paste("MCMC chain of ",lambda)))

   library(MCMCpack)

 pal<-vector(mode="numeric",length=50)
  for (i in 1:50)
   {
    print(paste("chain",i))
    pal[i]<-gelman.diag(mcmc.list(as.mcmc(newG[10001:20000,i]),as.mcmc(newG[20001:30000,i]),as.mcmc(newG[30001:40000,i]),as.mcmc(newG[40001:50000,i])))
    print(pal[i])
   }


    pallambda<-gelman.diag(mcmc.list(as.mcmc(newlambda[10001:20000]),as.mcmc(newlambda[20001:30000]),as.mcmc(newlambda[30001:40000]),as.mcmc(newlambda[40001:50000])))
    print(paste("Gelman-Rubin Lambda",pallambda))


 ####Q-Q plot: to show goodness of fit of model####

   yemp<-vector(mode="numeric",length=50)
   for (i in 1:50)
   {
    y_emp<-subset(new.data,new.data[,1] <= (ispl1new$x[i]+0.25) & new.data[,1] >= (ispl1new$x[i]-0.25))
    yemp[i]<-quantile(y_emp[,2],0.9)
   }

 plot(yemp,ispl1new$y,main="",xlab="Empirical Quantiles",ylab="Fitted Quantiles",xlim=c(0,2.5),ylim=c(0,2.5),cex.lab=1.3, cex.axis=1.3)
    lines(c(2.5,0), c(2.5,0))


    my.mse<-function(x,y)
    {
      y_er_sqr<-vector(mode="numeric",length=length(x))
      y_er<-vector(mode="numeric",length=length(x))

     for(i in 1:length(x))
      {
       y_er[i]<-(x[i]-y[i])
       y_er_sqr[i]<-(x[i]-y[i])^2

      }
       #hist(y_er,xlab="residuals",main="Histogram of Residual from Bayesian quantile spline")
       mser<-mean(y_er_sqr)
       return(mser)
    }
    my.mse(yemp,ispl1new$y)
    
    
rana.leg<-c("0.9 Parametric","0.9 Nonparametric","Parametric 95\% credible interval","Nonparametric 95\% credible interval")
legend("topleft",rana.leg,cex=1,lty=c(1,1,4,2),lwd=1.5,col=c("red","blue","red","blue"),bg="white")
