ky3<-function(p, s, r)
{
	prob <- function(y, x, p, beta)
	{
		v <- x %*% beta
		u <- (y - v)
		ind <- ifelse(u < 0., 1., 0.)
		z <-  - u * (p - ind)
		return(sum(z))
	}
	accept <- function(y, x, p, denom, thold, thnew, beta)
	{
		prob <- function(y, x, p, beta)
		{
			v <- x %*% beta
			u <- (y - v)
			ind <- ifelse(u < 0., 1., 0.)
			z <-  - u * (p - ind)
			return(sum(z))
		}
		#denom <- prob(y, x, p, thold)
		num <- prob(y, x, p, beta)
		mh <- (num - denom)
		mh <- exp(mh)
		mh <- min(mh, 1.)
		q <- runif(1., 0., 1.)
		if(mh > q) {
			th <- thnew
			den <- num
		}
		else {
			th <- thold
			den <- denom
		}
		list(th = th, den = den)
	}
	x <- serum.dat[, 1.:3.]
	y <- serum.dat[, 4.]
	nr <- dim(serum.dat)[1.]
	av1 <- mean(x[, 2.])
	x[, 2.] <- x[, 2.] - av1
	x[, 3.] <- x[, 2.]^2.
	b0 <- vector(mode = "numeric", length = r)
	b1 <- vector(mode = "numeric", length = r)
	b2 <- vector(mode = "numeric", length = r)
	b0old <- 0.
	b1old <- 0.
	b2old <- 0.
	sig0 <- 1.3
	sig1 <- 1.3
	sig2 <- 0.1
	beta <- c(b0old, b1old, b2old)
	denom <- prob(y, x, p, beta)
	b0[1.] <- b0old
	b1[1.] <- b1old
	b2[1.] <- b2old
	#par(mfrow = c(2, 2))
	#plot(1:r, b0, xlab = "Iteration number", ylab = "beta_0", type = "l")
	#plot(1:r, b1, xlab = "Iteration number", ylab = "beta_1", type = "l")
	#plot(1:r, b2, xlab = "Iteration number", ylab = "beta_2", type = "l")
	#stat0 <- c(summary(b0[(s + 1):r]), stdev(b0[(s + 1):r]))
	#stat1 <- c(summary(b1[(s + 1):r]), stdev(b1[(s + 1):r]))
	#stat2 <- c(summary(b2[(s + 1):r]), stdev(b2[(s + 1):r]))
	#hist(b0[(s + 1):r])
	#hist(b1[(s + 1):r])
	#hist(b2[(s + 1):r])
	#par(mfrow = c(1, 1))
	#stat <- rbind(stat0, stat1, stat2)
	#return(stat)
	for(i in 2.:r) {
		thnew <- rnorm(1., b0old, sig0)
		beta <- c(thnew, b1old, b2old)
		acc <- accept(y, x, p, denom, b0old, thnew, beta)
		b0old <- acc$th
		denom <- acc$den
		thnew <- rnorm(1., b1old, sig1)
		beta <- c(b0old, thnew, b2old)
		acc <- accept(y, x, p, denom, b1old, thnew, beta)
		b1old <- acc$th
		denom <- acc$den
		thnew <- rnorm(1., b2old, sig2)
		beta <- c(b0old, b1old, thnew)
		acc <- accept(y, x, p, denom, b2old, thnew, beta)
		b2old <- acc$th
		denom <- acc$den
		b2[i] <- b2old
		b1[i] <- b1old - (2. * b2old * av1)
		b0[i] <- b0old - b1old * av1 + b2old * av1 * av1
	}
#	par(mfrow = c(2, 2))
#	plot(1:r, b0, xlab = "Iteration number", ylab = "beta_0", type = "l")
#	plot(1:r, b1, xlab = "Iteration number", ylab = "beta_1", type = "l")
#	plot(1:r, b2, xlab = "Iteration number", ylab = "beta_2", type = "l")
	b012 <- cbind(b0, b1, b2)
	return(b012)
}



	z1 <- ky3(0.05, 1000., 6000.)
	z2 <- ky3(0.25, 1000., 6000.)
	z3 <- ky3(0.5, 1000., 6000.)
	z4 <- ky3(0.75, 1000., 6000.)
	z5 <- ky3(0.9, 1000., 6000.)
	ky3out <- cbind(z1, z2, z3, z4, z5)



	#zm <- apply(z, 2, mean)
	z <- ky3out[1001.:6000., 13.:15.]
	z<-z5
	y <- serum.dat[, 4.]
	x <- serum.dat[, 1.:3.]
	xr <- c(0., 6.5)
	yr <- c(0., 15.)
plot(x[, 2.], y, xlim = xr, ylim = yr, xlab = "Age (year)",
	ylab = "IgG (grams/litre)")
	xp.cub <- seq(0.5, 6., len = 50.)
		xp <- seq(0.5, 6., len = 50.)
	w <- cbind(rep(1., 50.), xp, xp^2.)
	z5 <- w %*% t(z)
	r <- 5000.
	h <- 0.025
	hl <- floor(h * r)
	hu <- ceiling((1. - h) * r)
	ww <- apply(t(z5), 2., sort)
	xl <- ww[hl,  ]
	xu <- ww[hu,  ]
xm.cub <- apply(ww, 2., mean)
	xm<-xm.cub
	lines(xp, xm,col="red",lwd=1.5)
	lines(xp, xl, lty = 4.,col="red",lwd=1.5)
	lines(xp, xu, lty = 4.,col="red",lwd=1.5)
	#text(6.5, xm[100.], "0.9")
#	z <- ky3out[1001.:6000., 1.:3.]
#	z1 <- w %*% t(z)
#	ww <- apply(t(z1), 2., sort)
#	xl <- ww[hl,  ]
#	xu <- ww[hu,  ]
#	xm<- apply(ww, 2., mean)
	#lines(xp, xm)
	#lines(xp, xl, lty = 4.)
	#lines(xp, xu, lty = 4.)
#	text(6.5, xm[100.], "0.05")
   print(xm)
    ####Q-Q plot: to show goodness of fit of model####
   
     yempcub<-vector(mode="numeric",length=50)
   for (i in 1:50)
   {
    y_empcub<-subset(serum.dat,serum.dat[, 2.] <= (xp[i]+0.25) & serum.dat[, 2.] >= (xp[i]-0.25))
    yempcub[i]<-quantile(y_empcub[,4],0.9)
   }
   
 plot(yemp,xm,main="Bayesian Quantile Cubic Regression",xlab="Empirical Quantiles",ylab="Fitted Quantiles",xlim=c(0,2.5),ylim=c(0,2.5),cex.lab=1.3, cex.axis=1.3)
    lines(c(2.5,0), c(2.5,0))
    
    my.mse<-function(x,y)
    {
      y_er_sqr<-vector(mode="numeric",length=length(x))
      y_er<-vector(mode="numeric",length=length(x))
     for(i in 1:length(x))
      {
       y_er[i]<- x[i]-y[i]
       y_er_sqr[i]<-(x[i]-y[i])^2
      
      }
     # hist(y_er,xlab="residuals",main="Histogram of Residual from BQR for polynomial")
       mser<-mean(y_er_sqr)
       return(mser)
    }
    my.mse(yemp,xm)