affine2distance<-function(affVals) {
	
	affVals<-as.numeric(affVals)
	
	dx<-affVals[10]
	dy<-affVals[11]
	dz<-affVals[12]

	rotx<-asin(affVals[7])
	roty<-atan2(affVals[8]/cos(rotx),affVals[9]/cos(rotx))
	rotz<-atan2(affVals[4]/cos(rotx),affVals[1]/cos(rotx))

	rotx<-rotx*360/(2*pi)
	roty<-roty*360/(2*pi)
	rotz<-rotz*360/(2*pi)

	return(c(dx, dy, dz, rotx, roty, rotz))
	
}

getMotionParams<-function(m) {

	if (length(as.numeric(as.matrix(m)))==16) {
		affVals<-as.numeric(m)[c(1:3,5:7,9:11,13:15)]
		p<-affine2distance(affVals)
		return(p)
	} else if (length(as.numeric(as.matrix(m)))==12) {
		affVals<-as.numeric(m)
		p<-affine2distance(affVals)
		return(p)
	} else if (dim(m)[2]==14) {
		m<-m[,3:14]
		pmat<-t(apply(m, 1, FUN=affine2distance))
		p.df<-as.data.frame(pmat)
		names(p.df)<-c("x","y","z","pitch","roll","yaw")
		write.csv(file="motion_6p.csv",p.df,quote=FALSE,row.names=FALSE)
		return(pmat)
	} else if (dim(m)[2]==12) {
		pmat<-t(apply(m, 1, FUN=affine2distance))
		p.df<-as.data.frame(pmat)
		names(p.df)<-c("x","y","z","pitch","roll","yaw")
		write.csv(file="motion_6p.csv",p.df,quote=FALSE,row.names=FALSE)
		return(pmat)
	}

}

plotMotion<-function(params) {
# Creates a simple line plot of 6 motion parameters.
# Input is a matrix or data frame with 6 columns:
# x-translation, y-translation, z-translation, pitch, roll, yaw

params<-as.matrix(params)
cmap<-c("orange","green2","steelblue2", "tomato3","forestgreen","navy")
plot(x=1:dim(params)[1],params[,1],type="l",col=cmap[1],ylim=c(-3,3),lwd=1,xlab="TR",ylab="translation (mm) or rotation (degrees)", main="Distance from Average BOLD")
for (i in 2:6) {
	lines(x=1:dim(params)[1],params[,i],type="l",col=cmap[i],ylim=c(-3,3),lwd=1)
}

legend(x=0.01*dim(params)[1],y=3,legend=c("x","y","z","pitch","roll","yaw"),fill=cmap, ncol=2)

}

