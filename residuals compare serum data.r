
##you need to run rana code.r and serumsplineprogram.r##


    yemp<-vector(mode="numeric",length=50)
   for (i in 1:50)
   {
    y_emp<-subset(serum.dat,serum.dat[, 2.] <= (ispl1new$x[i]+0.25) & serum.dat[, 2.] >= (ispl1new$x[i]-0.25))
    yemp[i]<-quantile(y_emp[,4],0.9)
   }


 res1<-abs(yemp-round(ispl1new$y,2))
res2<-abs(yemp-round(xm.cub,2))
par(mfrow=c(1,1))


plot(xp.cub,res2,main="",xlab="Cos(Wave direction)",ylab="|Residuals|",cex.lab=1.3, cex.axis=1.3,pch=4)
points(xp.cub,res1,pch=20)
lines(loess.smooth(xp.cub,res2),lty=2,lwd=2)
lines(loess.smooth(xp.cub,res1),lwd=2)
#lines(xp.cub,yemp)
legend("topleft", legend = c("Residuals from spline", "Residuals from cubic"), pch = c(20,4),lty=1:2,lwd=2,cex=1.2)