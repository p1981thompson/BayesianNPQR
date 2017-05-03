
        ##########################################################
        #                                                        #
        # Rana suggetsion for spline work                        #
        #                                                        #
        ##########################################################



      ##############################################
  ####LOAD AND SORT DATA FOR CORRECT QUANTILE CURVES####
      ##############################################

#new.dat<-read.table("F:/University work/R work/Rana QR/mysampsmlalt.txt",header=T)
serum.dat<-load("d://Local Data//p1thompson//University work//R work//serum.txt")
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

   #### SCATTER PLOT OF DATA #####
   par(mfrow=c(1,1))
plot(x, y, xlab = "Age(Year)", ylab = "IgG (grams/litre)",cex=0.5,pch=20, cex.lab=1.3, cex.axis=1.3)

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
       spline.pred<-predict(spline.fun.g,x.grid)
   	u <- (y - spline.pred$y)
    print(paste("u values",u))
		ind <- ifelse(u < 0, 1, 0)
		z <-  -0.5*lambda*(t(g)%*%K%*%g) - sum(u * (p - ind))
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
make.RQ <- function(unq.x)
{

Hlgth <- length(unq.x)

H <- diff(unq.x)

  R <- matrix( 0,Hlgth,Hlgth)
  Q <- matrix( 0,Hlgth,Hlgth)

         for(j in 2:(Hlgth-1))
                  {

                  Q[j-1,j] <- H[j-1]^(-1)
                  Q[j,j] <- -H[j-1]^(-1) - H[j]^(-1)
                  Q[j+1,j] <- H[j]^(-1)
                  R[j,j] <- (1/3)*(H[j-1] + H[j])

                  }



                  for(d in 2:(Hlgth-2))
                        {
                        R[d,d+1] <- (1/6)*H[d]
                        R[(d+1),d] <- (1/6)*H[d]
                        }



    Q <- Q[,-(Hlgth)]

    Q <- Q[,-1]

   R <- R[-(Hlgth),]

   R <- R[,-(Hlgth)]

   R <- R[-1,-1]

   return(list(R = R, Q = Q))

}
# --------------------------------------------------

  ### GRID LENGTH SET AS DEFAULT=20 BUT CAN BE ALTERED MANUALLY IN FUNCTIONS BELOW ###

ky3alt<-function(x,y,p,  r,sigma, lambda, grid.length = 20)
{
#
# This function requires the x's to be monotonically increasing
#


    unq.x <- unique(x)
    x.grid<-unq.x
    grid.length <- length(x.grid)

    #print(x.grid)


    Hlgth<-length(x.grid)



    my.RQ <- make.RQ(x.grid)
    R <- my.RQ$R
    Q <- my.RQ$Q

    print(R)
   # print(Q)

Rinv <- solve(R)

K <- Q%*%Rinv%*%t(Q)

  #print("Rinv=")
  #print(Rinv)

  #print("K=")
  #print(K)

	G <- matrix(0,r,grid.length)
 # X<- model.matrix(y~bs(x,df=10))

 ####INITIAL INPUT INTO THE MCMC ROUTINE, df CAN BE ALTERED TO ALLOW FASTER CONVERGENCE####

#qmodel <- rq(y ~ bs(x,df=6), tau = p)
#y.hat <- predict(qmodel, newdata = data.frame(x = x.grid))

qmodel <- rq(y ~ x + I(x^2) + I(x^3), tau = p)
y.hat <- predict(qmodel, newdata = data.frame(x = x.grid))

     gold <- y.hat

     #lines(x.grid, gold, col="orange")

    # print(gold)

     G[1,] <- gold
    #  G[1,]<-rep(1,length(x.grid))
	   denom <- prob(x, y, x.grid, gold, p, K, lambda)



	    no.acc <- 0

	for(i in 2:r) {

	#if(i %% 100 == 0){ print(paste("Iteration",i))
	#}

          	gnew <- gold + rnorm(grid.length, 0, sigma)
            acc <- accept(x, y, x.grid, gold, gnew, p, K, lambda, denom)


            gold <- acc$th
            denom <- acc$den
            G[i,]<-gold
            no.acc <- no.acc + acc$i.accept
                 }


                 print(paste("Proportion accepted", no.acc /(r - 1)))
 return(list(G = G, x.grid = x.grid,K=K))
 }


        ##############################################
        #                                            #
        #   2D PLOTS OF QUANTILE REGRESSION SPLINES  #
        #                                            #
        ##############################################

        ###FUNCTIONS FOR DIFFERENT QUANTILE CURVES###

   ######90% QUANTILE REGRESSION CURVE (lambda 1) ######
   ######                                                            ######                                      

    z1 <- ky3alt(x, y, 0.9 ,  10, 0.061, lambda =  smooth.spline(x,y)$lambda)

	z1new <- z1$G[2:10,]
  xr <- c(-1, 1)
	yr <- c(0, 3)
	r <- 8
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

  lines(ispl1new, col="purple",lty=2)
	#points(ispl1up, col="red",lty=2)          
	#points(ispl1low, col="red",lty=2)

	##REQUIRED FOR PERSP PLOT##
  Rmat1<-as.matrix(cbind(xp1,xm1,rep(6,length(xp1))))
  Rmat1up<-as.matrix(cbind(xp1,xu1,rep(6,length(xp1))))
  Rmat1low<-as.matrix(cbind(xp1,xl1,rep(6,length(xp1))))

   ######90% QUANTILE REGRESSION CURVE (lambda 2) ######
   ######                                                          ######

    z2 <- ky3alt(x, y, 0.9 ,  100000, 0.045, lambda =  2*smooth.spline(x,y)$lambda)

	z2new <- z2$G[25001:100000,]
  xr <- c(-1, 1)
	yr <- c(0, 3)
	r <- 75000
	ww2 <- apply(z2new, 2, sort)
	dim(ww)
 	xl2 <- ww2[hl,  ]
	xu2 <- ww2[hu,  ]
	xm2 <- apply(ww2, 2, mean)
  xp2<-sort(z2$x.grid)

  ##INTERPOLATING SPLINES FITTED THROUGH MCMC POINTS TO SMOOTH##

 ispl2<-spline(xp2,xm2,n=1000,method="natural")
  #ispl2up<-spline(xp2,xu2,method="natural")
  #ispl2low<- spline(xp2,xl2,method="natural")

  lines(ispl2, col="red")
	#lines(ispl2up, col="red",lty=2)
	#lines(ispl2low, col="red",lty=2)

	##REQUIRED FOR PERSP PLOT##
  Rmat2<-as.matrix(cbind(xp2,xm2,rep(5,length(xp2))))
  Rmat2up<-as.matrix(cbind(xp2,xu2,rep(5,length(xp2))))
  Rmat2low<-as.matrix(cbind(xp2,xl2,rep(5,length(xp2))))

   ######99% QUANTILE REGRESSION CURVE (lambda 3) ######
   ######                                                       ######

    z3 <- ky3alt(x, y, 0.9 ,  100000, 0.0072, lambda =  10*smooth.spline(x,y)$lambda)

	z3new <- z3$G[2501:100000,]
  xr <- c(-1, 1)
	yr <- c(0, 3)
	r <- 75000
	ww3 <- apply(z3new, 2, sort)
	dim(ww)
 	xl3 <- ww3[hl,  ]
	xu3 <- ww3[hu,  ]
	xm3 <- apply(ww3, 2, mean)
  xp3<-sort(z3$x.grid)

   ##INTERPOLATING SPLINES FITTED THROUGH MCMC POINTS TO SMOOTH##

  ispl3<-spline(xp3,xm3,n=1000,method="natural")
  #ispl3up<-spline(xp3,xu3,method="natural")
  #ispl3low<- spline(xp3,xl3,method="natural")

  lines(ispl3, col="green")
	#lines(ispl3up, col="red",lty=2)
	#lines(ispl3low, col="red",lty=2)

	##REQUIRED FOR PERSP PLOT##
  Rmat3<-as.matrix(cbind(xp3,xm3,rep(4,length(xp3))))
  Rmat3up<-as.matrix(cbind(xp3,xu3,rep(4,length(xp3))))
  Rmat3low<-as.matrix(cbind(xp3,xl3,rep(4,length(xp3))))

   ######90% QUANTILE REGRESSION CURVE (lambda 4) ######
   ######                                                     ######
  z4 <- ky3alt(x, y, 0.9 ,  100000, 0.395, lambda =  0.01*smooth.spline(x,y)$lambda)

	z4new <- z4$G[25001:100000,]
  xr <- c(-1, 1)
	yr <- c(0, 3)
	r <- 75000
	ww4 <- apply(z4new, 2, sort)
	dim(ww)
 	xl4 <- ww4[hl,  ]
	xu4 <- ww4[hu,  ]
	xm4 <- apply(ww4, 2, mean)
  xp4<-sort(z4$x.grid)

   ##INTERPOLATING SPLINES FITTED THROUGH MCMC POINTS TO SMOOTH##

  ispl4<-spline(xp4,xm4,n=1000,method="natural")
  #ispl4up<-spline(xp4,xu4,method="natural")
  #ispl4low<- spline(xp4,xl4,method="natural")

  lines(ispl4, col="orange")
	#lines(ispl4up, col="red",lty=2)
	#lines(ispl4low, col="red",lty=2)

	##REQUIRED FOR PERSP PLOT##
  Rmat4<-as.matrix(cbind(xp4,xm4,rep(3,length(xp4))))
  Rmat4up<-as.matrix(cbind(xp4,xu4,rep(3,length(xp4))))
  Rmat4low<-as.matrix(cbind(xp4,xl4,rep(3,length(xp4))))

   ######90% QUANTILE REGRESSION CURVE (lambda 5) ######
   ######                                                      ######

  z5 <- ky3alt(x, y, 0.9 , 100000, 0.48, lambda =  0.0001*smooth.spline(x,y)$lambda, grid.length = 80)

	z5new <- z5$G[25001:100000,]
  xr <- c(-1, 1)
	yr <- c(0, 3)
	r <- 75000
	h <- 0.025
	hl <- floor(h * r)
	hu <- ceiling((1 - h) * r)
	ww5 <- apply(z5new, 2, sort)
	dim(ww5)
 	xl5 <- ww5[hl,  ]
	xu5 <- ww5[hu,  ]
	xm5 <- apply(ww5, 2, mean)
  xp5<-sort(z5$x.grid)

  ##INTERPOLATING SPLINES FITTED THROUGH MCMC POINTS TO SMOOTH##

  ispl5<-spline(xp5,xm5,n=1000,method="natural")
 # ispl5up<-spline(xp5,xu5,n=1000,method="natural")
  #ispl5low<- spline(xp5,xl5,n=100,method="natural")

  lines(ispl5, col="purple")
	#lines(ispl5up, col="red",lty=2,lwd=2)
#	lines(ispl5low, col="red",lty=2,lwd=2)

  thin.seq<-seq(from=1,to=1000000,by=100)

  par(mfrow=c(3,2))
  for (i in 1:64)
   {
    ts.plot(z1$G[,i])
   }

   newG<-matrix(0,10000,42)
  for (i in 1:22)
   {
      newG[,i]<-z5$G[thin.seq,i]
   }

   pal<-vector(mode="numeric",length=22)

   for (i in 1:22)
   {
    print(paste("chain",i))
    pal[i]<-gelman.diag(mcmc.list(as.mcmc(newG[2001:4000,i]),as.mcmc(newG[4001:6000,i]),as.mcmc(newG[6001:8000,i]),as.mcmc(newG[8001:10000,i])))
    print(pal[i])
   }

   ##GELMAN RUBIN PLOTS (CONVERGENCE CHECKS)##
  par(mfrow=c(1,1))
  gelman.plot(mcmc.list(as.mcmc(z5$G[100001:200000,5]),as.mcmc(z5$G[200001:300000,5]),as.mcmc(z5$G[300001:400000,5]),as.mcmc(z5$G[400001:500000,5])),main="chain5")
  gelman.plot(mcmc.list(as.mcmc(z5$G[100001:200000,10]),as.mcmc(z5$G[200001:300000,10]),as.mcmc(z5$G[300001:400000,10]),as.mcmc(z5$G[400001:500000,10])),main="chain10")
  gelman.plot(mcmc.list(as.mcmc(z5$G[100001:200000,15]),as.mcmc(z5$G[200001:300000,15]),as.mcmc(z5$G[300001:400000,15]),as.mcmc(z5$G[400001:500000,15])),main="chain15")
  gelman.plot(mcmc.list(as.mcmc(z5$G[100001:200000,21]),as.mcmc(z5$G[200001:300000,21]),as.mcmc(z5$G[300001:400000,21]),as.mcmc(z5$G[400001:500000,21])),main="chain21")
  gelman.plot(mcmc.list(as.mcmc(z5$G[100001:200000,1]),as.mcmc(z5$G[200001:300000,1]),as.mcmc(z5$G[300001:400000,1]),as.mcmc(z5$G[400001:500000,1])),main="chain1")

  gelman.diag(mcmc.list(as.mcmc(z5$G[100001:200000,5]),as.mcmc(z5$G[200001:300000,5]),as.mcmc(z5$G[300001:400000,5]),as.mcmc(z5$G[400001:500000,5])))
  gelman.diag(mcmc.list(as.mcmc(z5$G[100001:200000,10]),as.mcmc(z5$G[200001:300000,10]),as.mcmc(z5$G[300001:400000,10]),as.mcmc(z5$G[400001:500000,10])))
  gelman.diag(mcmc.list(as.mcmc(z5$G[100001:200000,15]),as.mcmc(z5$G[200001:300000,15]),as.mcmc(z5$G[300001:400000,15]),as.mcmc(z5$G[400001:500000,15])))
  gelman.diag(mcmc.list(as.mcmc(z5$G[100001:200000,21]),as.mcmc(z5$G[200001:300000,21]),as.mcmc(z5$G[300001:400000,21]),as.mcmc(z5$G[400001:500000,21])))
  gelman.diag(mcmc.list(as.mcmc(z5$G[100001:200000,1]),as.mcmc(z5$G[200001:300000,1]),as.mcmc(z5$G[300001:400000,1]),as.mcmc(z5$G[400001:500000,1])))

   par(mfrow=c(3,2))
  for (i in 1:22)
   {
    gelman.plot(mcmc.list(as.mcmc(z5$G[100001:200000,i]),as.mcmc(z5$G[200001:300000,i]),as.mcmc(z5$G[300001:400000,i]),as.mcmc(z5$G[400001:500000,i])))
   }
   pal<-vector(mode="numeric",length=22)
  for (i in 1:22)
   {
    print(paste("chain",i))
    pal[i]<-gelman.diag(mcmc.list(as.mcmc(z5$G[100001:200000,i]),as.mcmc(z5$G[200001:300000,i]),as.mcmc(z5$G[300001:400000,i]),as.mcmc(z5$G[400001:500000,i])))
    print(pal[i])
   }

	##REQUIRED FOR PERSP PLOT##
  Rmat5<-as.matrix(cbind(xp5,xm5,rep(2,length(xp5))))
  Rmat5up<-as.matrix(cbind(xp5,xu5,rep(2,length(xp5))))
  Rmat5low<-as.matrix(cbind(xp5,xl5,rep(2,length(xp5))))


   ######PLOT LEGEND SETUP######

   leg2.txt<-c("BQR using interpolating spline at 0.9,lambda1","BQR using interpolating spline at 0.9,lambda 2","BQR using interpolating spline at 0.9, lambda 3","BQR using interpolating spline at 0.9, lambda 4", "BQR using interpolating spline at 0.9, lambda 5")
   legend(x=1,y=12,leg2.txt,col=c("blue","red","green","orange","purple"),cex=0.7,lty=c(1,1,1,1,1),lwd=1.5,bg="white")



   ####Q-Q plot: to show goodness of fit of model####

   yemp<-vector(mode="numeric",length=100)
      for (i in 1:100)
   {
    y_emp<-subset(new.data,new.data[,2] == ispl1$x[i])
    print(y_emp)
    yemp[i]<-quantile(y_emp[,1],0.9)
    print(paste("yemp[i]",yemp[i]))
   }

 plot(yemp,ispl1$y,main="Bayesian Quantile Regression Spline",xlab="Empirical Quantiles",ylab="Fitted Quantiles",cex.lab=1.3, cex.axis=1.3)
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
    my.mse(yemp,ispl1$y)

       unq.x1<-as.vector(unique(x))

     yemp<-vector(mode="numeric",length=unique(x))
      for (i in 1:unique(x))
   {
    y_emp<-subset(new.data,new.data[,2] == unq.x1[i])
    print(y_emp)
    yemp[i]<-quantile(y_emp[,1],0.9)
    print(paste("yemp[i]",yemp[i]))
   }



   #####################################
   bin.seq<-seq(min(x),max(x),length=33)
   
  yemp<-vector(mode="numeric",length=32)
      for (i in 1:32)
   {
    y_emp<-subset(new.data,new.data[,1] >=bin.seq[i] & new.data[,1] < bin.seq[i+1])
   
    yemp[i]<-quantile(y_emp[,2],0.9)
    print(paste("yemp",i))
    print(yemp[i])
   }
    difer<-diff(bin.seq)
    xpos<-vector(mode="numeric",length=32)
   for (i in 1:32)
   {
   xpos[i]<-(difer[i]/2)+bin.seq[i]
   }
   
   #####################################
   
    emp.y<-predict(ispl1,xpos)
    
    plot(yemp,emp.y$y,main="Bayesian Quantile Regression Spline",xlab="Empirical Quantiles",ylab="Fitted Quantiles",cex.lab=1.3, cex.axis=1.3)
    lines(c(12,0), c(12,0))
    
    
    
 plot(x, y, xlab = "Age(Year)", ylab = "IgG (grams/litre)",cex=0.5,pch=20, cex.lab=1.3, cex.axis=1.3)
 abline(v=bin.seq)
 lines(xpos,emp.y$y,col="blue")
 lines(xpos,yemp,col="red")