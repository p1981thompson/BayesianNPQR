
  library(plotrix)
  new.dat<-read.table("d://Local Data//p1thompson//University work//R work/Rana QR//my.samp.txt",header=T)

myradial.plot<-function (lengths, radial.pos, labels, label.pos, rp.type = "r",
    label.prop = 1.1, main = "", xlab = "", ylab = "", line.col = par("fg"),
    mar = c(2, 2, 3, 2), show.grid = TRUE, show.radial.grid = TRUE,
    grid.col = "gray", grid.bg = "transparent", point.symbols = NULL,
    point.col = NULL, show.centroid = FALSE, radial.lim = NULL,
    ...)
{
    length.dim <- dim(lengths)
    if (is.null(radial.lim))
        radial.lim <- range(lengths)
    if (is.null(length.dim)) {
        npoints <- length(lengths)
        nsets <- 1
        lengths <- matrix(lengths, nrow = 1)
    }
    else {
        npoints <- length.dim[2]
        nsets <- length.dim[1]
    }
    if (missing(radial.pos))
        radial.pos <- seq(0, pi * (2 - 2/npoints), length = npoints)
    radial.pos.dim <- dim(radial.pos)
    if (is.null(radial.pos.dim))
        radial.pos <- matrix(rep(radial.pos, nsets), nrow = nsets,
            byrow = TRUE)
    if (show.grid) {
        grid.pos <- pretty(radial.lim)
        if (grid.pos[1] <= 0)
            grid.pos <- grid.pos[-1]
        maxlength <- max(grid.pos)
        angles <- seq(0, 1.96 * pi, by = 0.04 * pi)
    }
    else {
        grid.pos <- NA
        maxlength <- max(radial.lim)
    }
    oldpar <- par(no.readonly = TRUE)
    par(mar = mar, pty = "s")
    plot(c(-maxlength, maxlength), c(-maxlength, maxlength),
        type = "n", axes = FALSE, main = main, xlab = xlab, ylab = ylab,
        ...)
    par(xpd = TRUE)
    if (length(line.col) < nsets)
        line.col <- 1:nsets
    rp.type <- unlist(strsplit(rp.type, ""))
    if (match("s", rp.type, 0)) {
        if (is.null(point.symbols))
            point.symbols <- 1:nsets
        if (length(point.symbols) < nsets)
            point.symbols <- rep(point.symbols, length.out = nsets)
        if (is.null(point.col))
            point.col <- 1:nsets
        if (length(point.col) < nsets)
            point.col <- rep(point.col, length.out = nsets)
    }
    for (i in 1:nsets) {
        xpos <- cos(radial.pos[i, ]) * lengths[i, ]
        ypos <- sin(radial.pos[i, ]) * lengths[i, ]
        if (match("r", rp.type, 0))
            segments(0, 0, xpos, ypos, col = line.col[i], ...)
        if (match("p", rp.type, 0))
            polygon(xpos, ypos, border = line.col[i], col = NA,
                ...)
        if (match("s", rp.type, 0))
            points(xpos, ypos, pch = point.symbols[i], col = point.col[i],
                ...)
        if (show.centroid)
            points(mean(xpos), mean(ypos), col = point.col[i],
                pch = point.symbols[i], cex = 2, ...)
    }
    #print(length(xpos))
    #print(length(ypos))
    my.dat<-data.frame(xpos,ypos)
    if (missing(labels)) {
        if (length(radial.pos) <= 20) {
            labels <- as.character(round(radial.pos, 2))
            label.pos <- radial.pos
        }
        else {
            label.pos <- seq(0, 1.8 * pi, length = 9)
            labels <- as.character(round(label.pos, 2))
        }
    }
    if (missing(label.pos))
        label.pos <- seq(0, pi * (2 - 2/npoints), length = npoints)
    xpos <- cos(label.pos) * maxlength
    ypos <- sin(label.pos) * maxlength
    if (show.radial.grid)
        segments(0, 0, xpos, ypos, col = grid.col)
    xpos <- cos(label.pos) * maxlength * label.prop
    ypos <- sin(label.pos) * maxlength * label.prop
    boxed.labels(xpos, ypos, labels, ypad = 0.7, border = FALSE)
    if (show.grid) {
        for (i in seq(length(grid.pos), 1, by = -1)) {
            xpos <- cos(angles) * grid.pos[i]
            ypos <- sin(angles) * grid.pos[i]
            polygon(xpos, ypos, border = grid.col, col = grid.bg)
        }
        ypos <- rep(-maxlength/15, length(grid.pos))
        boxed.labels(grid.pos, ypos, as.character(grid.pos),
            border = FALSE)
    }
    par(oldpar)

    return(my.dat)
}

 paulradialplot<-myradial.plot(my.samp$Hs,my.samp$wavdir*pi/180,rp.type="s")

 library(rgl)
 open3d()

 points3d(paulradialplot$xpos,paulradialplot$ypos,my.samp$Tz,size=2)
     #npoints3d <- length(my.samp$Hs)
     npoints3d <- 16
    radial.lim3d <- range(my.samp$Hs)
    grid.pos3d <- pretty(radial.lim3d)
        if (grid.pos3d[1] <= 0)
            grid.pos3d <- grid.pos3d[-1]
        maxlength3d <- max(grid.pos3d)

    label.pos3d <- seq(0, pi * (2 - 2/npoints3d), length = npoints3d)
    xposnew <- cos(label.pos3d) * maxlength3d
    yposnew <- sin(label.pos3d) * maxlength3d
     for(g in 1:length(xposnew))
        {
        lines3d(x=c(0,xposnew[g]), y=c(0,yposnew[g]), c(0,0),col="grey")
        }
        labels <- as.character(round(label.pos3d, 2))
     texts3d(xposnew, yposnew,rep(0,length(xposnew)), labels)
     axis3d('z')

     for( b in 2:length(grid.pos3d) )
     {
      xposnew <- cos(label.pos3d) * grid.pos3d[b]
      yposnew <- sin(label.pos3d) * grid.pos3d[b]
      lines3d(x=c(0,xposnew), y=c(0,yposnew), c(0,0),col="grey")
     }

     yposlab <- rep(-maxlength3d/15, length(grid.pos3d))
        texts3d(grid.pos3d, yposlab,rep(0,length(grid.pos3d)), as.character(grid.pos3d))
      title3d(main = '3D bivariate scatter plot including direction covariate', zlab = 'Wave Period')



   rgl.snapshot("d:/Local Data/p1thompson/University work/R work/program plots/myradialplot1.png",fmt="png")


   ###########add a surface ############


 y<-paulradialplot$ypos
 x<-paulradialplot$xpos
 z<-my.samp$Tz
 surfacedat<-as.data.frame(cbind(x,y,z))


  surface.loess<-loess(z~x*y,surfacedat)
  surface.mar<-list(x=seq(from=min(x),to=max(x),by=(max(x)-min(x))/100),y=seq(from=min(y),to=max(y),by=(max(y)-min(y))/100))
  surface.lo<-predict(surface.loess,expand.grid(surface.mar))

 surface3d(surface.mar$x,surface.mar$y,surface.lo,col=heat.colors(100))

  rgl.snapshot("d:/Local Data/p1thompson/University work/R work/program plots/myradialplotsurface.png",fmt="png")
  
      ###############################
 ##### Varying threshold spline code ######
      ###############################

   #new.dat<-read.table("F:/University work/R work/Rana QR/mysampsmlalt.txt",header=T)
new.dat<-read.table("d://Local Data//p1thompson//University work//R work/Rana QR//my.samp.txt",header=T)
#new.dat<-read.table("d://Local Data//p1thompson//University work//R work/blank//azdatgui.txt",header=T)
#
# Sort out the data
#
new.data<-data.frame(new.dat[,6],new.dat[,7])
#
# new.data<-new.data[ do.call(order, new.data) ,]
#
new.data <- new.data[order(new.data[,2]), ]

new.data<-as.matrix(new.data)

#
# Check by plotting
#
y <- new.data[,1]
x <- new.data[,2]*pi/180

   #### SCATTER PLOT OF DATA #####

plot(x, y, xlab = "Cos(Wave direction)", ylab = "Wave Height",cex=0.5,pch=20,col="skyblue", cex.lab=1.3, cex.axis=1.3)

   #### LOAD REQUIRED LIBRARIES FOR PROGRAM####

library(rgl)
library(quantreg)
library(splines)

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
    if(grid.length > length(unq.x)){
       warning("Finer grid than data")
    }

    length.unq.x <- length(unq.x)
    int.interval <- length.unq.x %/% grid.length
    if(int.interval == 0) int.interval <- 1

    ind.uni.x.grid <- seq(from = 1, to = length(unq.x), by = int.interval)

    x.grid <- unique(c(unq.x[ind.uni.x.grid], max(unq.x)))

    grid.length <- length(x.grid)

    #print(x.grid)


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
 # X<- model.matrix(y~bs(x,df=10))

 ####INITIAL INPUT INTO THE MCMC ROUTINE, df CAN BE ALTERED TO ALLOW FASTER CONVERGENCE####

qmodel <- rq(y ~ bs(x,df=7), tau = p)
y.hat <- predict(qmodel, newdata = data.frame(x = x.grid))

     gold <- y.hat

    # lines(x.grid, gold, col="yellow")

    # print(gold)

     G[1,] <- gold

	   denom <- prob(x, y, x.grid, gold, p, K, lambda)



	    no.acc <- 0

	for(i in 2:r) {

	if(i %% 100 == 0){ print(paste("Iteration",i))
	}

          	gnew <- gold + rnorm(grid.length, 0, sigma)
            acc <- accept(x, y, x.grid, gold, gnew, p, K, lambda, denom)


            gold <- acc$th
            denom <- acc$den
            G[i,]<-gold
            no.acc <- no.acc + acc$i.accept
                 }


                 print(paste("Proportion accepted", no.acc /(r - 1)))
 return(list(G = G, x.grid = x.grid))
 }


   ######90% QUANTILE REGRESSION CURVE (10 YEAR RETURN PERIOD) ######
   ######                                                      ######
    z5 <- ky3alt(x, y, 0.9 ,  10000, 0.0155, lambda =  smooth.spline(x,y)$lambda, grid.length = 20)

	z5new <- z5$G[2501:10000,]
  xr <- c(-1, 1)
	yr <- c(0, 3)
	r <- 7500
	ww5 <- apply(z5new, 2, sort)
	dim(ww)
 	xl5 <- ww5[hl,  ]
	xu5 <- ww5[hu,  ]
	xm5tz <- apply(ww5, 2, mean)
  xp5<-sort(z5$x.grid)


      ###################
  #### varying threshold #####
      ###################

    Hs.dim <- length(xm5)
    npointsv <- Hs.dim
    radial.posv <- xp5
  

        xposv <- cos(radial.posv) * xm5
        yposv <- sin(radial.posv) * xm5

    #lines3d(xposv, yposv,rep(0,length(xposv)))
    lines3d(xposv, yposv,xm5tz,col="red")
    ###########################
    
    rgl.snapshot("d:/Local Data/p1thompson/University work/R work/program plots/myradialplotvariablethres2.png",fmt="png")