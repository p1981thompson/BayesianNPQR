
        ##########################################################
        #                                                        #
        # 3D RETURN LEVEL PLOTS FOR QUANTILE REGRESSION SPLINES  #
        #                                                        #
        ##########################################################



      ##############################################
  ####LOAD AND SORT DATA FOR CORRECT QUANTILE CURVES####
      ##############################################

#new.dat<-read.table("F:/University work/R work/Rana QR/mysampsmlalt.txt",header=T)
new.dat<-read.table("d://Local Data//p1thompson//University work//R work/Rana QR//my.samp.txt",header=T)
#new.dat<-read.table("d://Local Data//p1thompson//University work//R work/blank//azdatgui.txt",header=T)
#
# Sort out the data
#
new.data<-data.frame(new.dat[,5],new.dat[,8])
#
# new.data<-new.data[ do.call(order, new.data) ,]
#
new.data <- new.data[order(new.data[,2]), ]

new.data<-as.matrix(new.data)

#
# Check by plotting
#
y <- new.data[,1]
x <- new.data[,2]

   #### SCATTER PLOT OF DATA #####

plot(x, y, xlab = "Cos(Wave direction)", ylab = "Wave Height",cex=0.5,pch=20,col="skyblue", cex.lab=1.3, cex.axis=1.3)

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

 # print(Q)

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
# --------------------------------------------------

get.foe <- function(G){
Glgth<-length(G[,1])-1
diffvec<-vector(mode = "numeric", length = Glgth)
for(i in 1:Glgth)
    {
     diffvec[i]<-sqrt(sum((G[i,]-G[(i+1),])^2))
    }
   diff2vec<-(diffvec)^2
   foe<-mean(diff2vec)

return(foe)
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
   # print(x.grid)


    Hlgth<-length(x.grid)



    my.RQ <- make.RQ(x.grid)
    R <- my.RQ$R
    Q <- my.RQ$Q

    #print(R)
   # print(Q)

Rinv <- solve(R)

K <- Q%*%Rinv%*%t(Q)

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

#qmodel <- rq(y ~ bs(x,df=50), tau = p)
#y.hat <- predict(qmodel, newdata = data.frame(x = x.grid))

qmodel <- rq(y ~ x + I(x^2) + I(x^3), tau = p)
y.hat <- predict(qmodel, newdata = data.frame(x = x.grid))
#y.hat <-rep(0.5,length(x.grid))
     gold <- y.hat

     lines(x.grid, gold, col="green", lwd = 2)

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

	#if(i %% 100 == 0){print(paste("Iteration",i))}



          	gnew <- gold + rnorm(grid.length, 0, sigma)
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
              FOE<-get.foe(G)
    alpha0<- no.acc /(r - 1)
                # cat(" ","\n")
                # print(paste("Proportion accepted", no.acc /(r - 1)))
                # print(paste("Proportion lambda accepted", no.acc2 /(r - 1)))
                # print(paste("FOE",FOE))

 plot.dat<-c(alpha0,FOE)
return(plot.dat)


# return(list(G = G, x.grid = x.grid,K = K,R=R,Q=Q,Lambda = Lambda))
 }


varoptim<- function(noint=50)
 {
  dat.new <- matrix(0, noint, 2)

  for (w in 1:noint)
    {
      print(paste("iteration number",w))
      upperlim0<-0.1
      lowerlim0<-0.0005
      sep0<-(upperlim0-lowerlim0)/49
      seq0<-seq(from=lowerlim0, to= upperlim0, by=sep0)
      dat.new[w,]<-ky3alt(x,y,0.9,5000,sigma=seq0[w],grid.length = 15,0.3)
   }

    seq0<-seq(from=0.0005, to= 0.1, by=(0.1-0.0005)/49)
    return(data.frame(dat.new,seq0))
 }

Paul<-varoptim()

    #### ASSOCIATED PLOTS ####

par(mfrow=c(1,1))

plot(Paul[,3],Paul[,2],main="Factor of Effiency vs Variance",xlab="Variance",ylab="FOE",pch=20)

plot(Paul[,1],Paul[,2],main="Factor of Efficiency vs Acceptance",xlab="Acceptance Rate",ylab="FOE",pch=20)

plot(Paul[,3],Paul[,1],main="Acceptance vs Variance", xlab="Variance", ylab="Acceptance Rate", pch=20)

plot(x, y, xlab = "Cos(Wave direction)", ylab = "Wave Height",cex=0.5,pch=20,col="skyblue", cex.lab=1.3, cex.axis=1.3)





    z1 <- ky3alt(x, y, 0.9 ,  100000, 0.5, grid.length = 40,0.3)


	z1new <- z1$G[101:2000,]
  xr <- c(-1, 1)
	yr <- c(0, 3)
	r <- 1900
	h <- 0.025
	hl <- floor(h * r)
	hu <- ceiling((1 - h) * r)
	ww1 <- apply(z1new, 2, sort)
	dim(ww)
	xl1 <- ww1[hl,  ]
	xu1 <- ww1[hu,  ]
	xm1 <- apply(ww1, 2, mean)
  xp1<-sort(z1$x.grid)

  ##INTERPOLATING SPLINES FITTED THROUGH MCMC POINTS TO SMOOTH##

  ispl1<-interpSpline(xp1,xm1)
  ispl1new<-spline(xp1,xm1,n=1000,method="natural")
  #ispl1up<-spline(xp1,xu1,method="natural")
  #ispl1low<- spline(xp1,xl1,method="natural")

  lines(ispl1new, col="blue",lty=2)

   par(mfrow=c(3,2))
  for (i in 1:64)
   {
    ts.plot(z1$G[,i])
   }

   par(mfrow=c(1,1))

    ts.plot(z1$Lambda)