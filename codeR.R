########################################
## 
## Purpose:
## 
##
## Author:
##   Dorian Giulliano Talledo Torres
##
## Version:
##   1.0    
##
## Date:
##   30/04/2020
########################################
###############
## Libraries ##
library(mgcv)
###############


readTheData = function(filepath){
########################################
#
# Purpose: Read the csv file
#   
# Inputs:  filepath
#   
# Outputs: dataset
#   
########################################

  data = read.csv(filepath, header = TRUE, sep = ",")

  return(data)
}

grid = function(x){
########################################
#
# Purpose: create a real line
#   
# Inputs:  x = variable of interest
#   
# Outputs: sequence
#   
########################################
	
	n = length(x)
	grid = seq((min(x)), (max(x)), length.out = n)

	return(grid)

}


K = function(u,d){
########################################
#
# Purpose: define the Kernel function 
#   
# Inputs:  u = ((z - x[i])/h), d = distribution
#   
# Outputs: K(.)
#   
########################################
    
    Uniform = 1/2
    Triangle = (1 - abs(u))
    Epanechnikov = (3/4)*(1-(u^2))
    Quartic = (15/16)*(1-(u)^2)^2
    Triweight = (35/32)*(1-(u)^2)^3
	Gaussian = (1/sqrt(2*pi))*exp(-0.5*(u)^2)
	Cosine = (pi/4)*cos((pi/2)*u)


	if (d =='Uniform'){
		return(Uniform)
	}
	else if (d =='Triangle'){
		return(Triangle)
	}
	else if (d =='Epanechnikov'){
		return(Epanechnikov)
	}
	else if (d =='Quartic'){
		return(Quartic)
	}
	else if (d =='Triweight'){
		return(Triweight)
	}
	else if (d =='Gaussian'){
		return(Gaussian)
	}
	else if (d =='Cosine'){
		return(Cosine)
	}
}

fh = function(x,h,d){
########################################
#
# Purpose: create Kernel density estimator
#   
# Inputs: x = variable of interest, h = bandwidth parameter, d = distribution
#   
# Outputs: Kernel density estimator
#   
########################################

	n = length(x)
	z = grid(x)

	fh = 0
	for(i in 1:n){
		u = ((z - x[i])/h)
		fh = fh + K(u,d)*(1/(n*h))
	}
	return(fh)

}

SilvermanRule = function(x){
########################################
#
# Purpose: Obtain the bandwidth parameter via Silverman's rule of thumb
#   
# Inputs: x = variable of interest
#   
# Outputs: h0 (optimal h)
#   
########################################

	n = length(x)
	h =  sd(x)*(4/(3*n))^(1/5)

	return(h)
}

mDynamic = function(z,p){
########################################
#
# Purpose: Create a dynamic matrix (n x p)
#   
# Inputs: z = variable of interest, p = degree of the polynomial
#   
# Outputs: Dynamic matrix (n x p)
#   
########################################

    
    mDynamic = c()

	for (i in seq(1, by=length(z), length=p+1)){
		mDynamic[seq(i, length=length(z))] = z
	}

	dim(mDynamic) = c(length(z),p+1)


	if (p==0){

    	mDynamic = z^p
    	dim(mDynamic) = c(length(z),p+1)

    }else if (p>0) {
    	for (i in 1:p+1){
    		mDynamic[,i] = mDynamic[,i]^(i-1)
    	}
    	mDynamic[,1]=mDynamic[,1]^0

    	dim(mDynamic) = c(length(z),p+1)


    }
	
	return(mDynamic)

}

WLS = function(x,y,d,h,p){
########################################
#
# Purpose: create a WLS estimator
#   
# Inputs: x,y = variables of interest, d = distribution, h = bandwidth parameter p = degree of the polynomial
#   
# Outputs: WLS estimator llest
#   
########################################

	n = length(x)
	xval = seq(min(x),max(x),length.out = n)

	llest = rep(0, n)

	for(i in 1:n){
		z = xval[i]-x
		X = mDynamic(z,p)
		Y = y
		W = K((z/h),d)
		B = (solve((t(X)%*%diag(W)%*%X))%*%(t(X))%*%diag(W)%*%Y)
		llest[i] = B[1]
	}

	return(list(llest,xval))

}

NWe = function(x,y,hCV,d){
########################################
#
# Purpose: define Nadaraya-Watson estimator
#   
# Inputs: x,y = variables of interest, d = distribution, h = CV bandwidth parameter
#   
# Outputs: Nadaraya-Watson estimator
#   
########################################

	z = grid(x)
	NWe = ksmooth(x,y, kernel=d, bandwidth=hCV, x.points=z)

	return(NWe)
}

CV = function(x,y,d,h){
########################################
#
# Purpose: define the CV function (-i) to be optimized to obtain hCV
#   
# Inputs: x,y = variables of interest, d = distribution, h = initial bandwidth parameter
#   
# Outputs: CV
#   
########################################

	n = length(x)
	mcv = rep(0,n)

	for (i in 1:n){
		mcv[i] = ksmooth(x[-i],y[-i], kernel=d,bandwidth=h, x.points=x[i])$y
	}

	CV = mean((y-mcv)^2)
	return(CV)

}

CI = function(mu,q,s,n){
########################################
#
# Purpose: define confidence intervals for Nadaraya-Watson estimator
#   
# Inputs: mu = prediction of the NWe, alpha lv, s = s2*K2, n = n*hCV*fx$x
#   
# Outputs: confidence intervals
#   
########################################

	CIub = mu$y + qnorm(q)*sqrt(s)/sqrt(n)
	CIlb = mu$y - qnorm(q)*sqrt(s)/sqrt(n)

	CI = list(CIub,CIlb)

	return(CI)

}

standardize = function(x){
########################################
#
# Purpose: create a function that standardizes variables
#   
# Inputs: x = variable of interest
#   
# Outputs: z = standardized variable
#   
########################################

	n = length(x)
	mu = mean(x)
	sigma = sd(x)
	rate = sqrt(n)

	z = (x - mu)/(sigma)

	return(z)

}


ex1 = function(tritium,n){

	h = SilvermanRule(tritium)
	fx = fh(tritium,h,'Gaussian')
	d = density(tritium, bw=h, kernel = "gaussian")
	x = grid(tritium)


	plot(x, fx, type="l", xlab="Tritium", ylab="Kernel density", main="Kernel density for tritium using a Gaussian kernel K(.)")
	points(tritium,rep(0,n), pch = 3)
	legend("topright", c("Kernel Density"), lty=c(1), col=c("black"),bty = "n")

	plot(d$x, d$y, type="l", xlab="Tritium", ylab="Kernel density", main="Kernel density for tritium using a Gaussian kernel K(.)")
	points(tritium,rep(0,n), pch = 3)
	legend("topright", c("Kernel Density"), lty=c(1), col=c("black"),bty = "n")

	print('Question 1')
	print('h chosen via Silverman s rule of thumb')
	print(h)
}


ex2 = function(pressure,tritium,n){

	h0 = optim(par=0.2, fn=function(h) CV(pressure,tritium,'normal',h), method="BFGS")
    hCV = h0$par

    NWe = NWe(pressure,tritium,hCV,'normal')

	plot(pressure,tritium, col="gray40", ylab="Tritium", xlab="Pressure", main="Regression of tritium on pressure using Nadaraya-Watson estimates")
	points(NWe$x,NWe$y,type="l",col=2,lwd=2)
	legend("topright", c("NW estimator"), lty=c( 1), col=c("red"),bty = "n")

	print('Question 2')
	print('hCV')
	print(hCV)

}


ex3 = function(pressure,tritium,n){

	h0 = optim(par=0.2, fn=function(h) CV(pressure,tritium,'normal',h), method="BFGS")
    hCV = h0$par

	NWe1 = ksmooth(pressure,tritium, kernel="normal", bandwidth=hCV, x.points=pressure)
	fx = density(pressure,from=min(pressure),to=max(pressure), n=n)

	s2 = mean((tritium[order(pressure)]-NWe1$y)^2)
	K2 = 2*qnorm(0.75)/sqrt(pi)
	NWe2 = ksmooth(pressure,tritium, kernel='normal', bandwidth=hCV, x.points=fx$x)
	s = s2*K2       #numerator
	n1 = n*hCV*fx$y #denominator

	alpha = 0.01
	q = 1-alpha/2
	

	CIub = CI(NWe2,q,s,n1)[[1]] # Upper Bound
	CIlb = CI(NWe2,q,s,n1)[[2]] # Lower Bound

	plot(pressure,tritium,col="gray50", ylab="Tritium", xlab="Pressure", main="Tritium on pressure using Nadaraya-Watson estimates")
	points(NWe2$x,NWe2$y,type = "l",col=2,lwd=2)
	points(NWe2$x,CIub,type = "l",col="blue",lwd=2,lty=4)
	points(NWe2$x,CIlb,type = "l",col="blue",lwd=2,lty=4)
	legend("topright", c("NW estimator",expression("CI"[99])), lty=c(1,4), col=c("red","blue"),bty = "n")
	
}

ex4 = function(pressure,tritium,n){

	hCV = 100000 # high h

	NWe1 = ksmooth(pressure,tritium, kernel="normal", bandwidth=hCV, x.points=pressure)
	fx = density(pressure,from=min(pressure),to=max(pressure), n=n)

	s2 = mean((tritium[order(pressure)]-NWe1$y)^2)
	K2 = 2*qnorm(0.75)/sqrt(pi)
	NWe2 = ksmooth(pressure,tritium, kernel='normal', bandwidth=hCV, x.points=fx$x)
	s = s2*K2       #numerator
	n1 = n*hCV*fx$y #denominator

	alpha = 0.01
	q = 1-alpha/2
	

	CIub = CI(NWe2,q,s,n1)[[1]] # Upper Bound
	CIlb = CI(NWe2,q,s,n1)[[2]] # Lower Bound

	plot(pressure,tritium,col="gray50", ylab="Tritium", xlab="Pressure", main="Tritium on pressure using Nadaraya-Watson estimates")
	points(NWe2$x,NWe2$y,type = "l",col=2,lwd=2)
	points(NWe2$x,CIub,type = "l",col="blue",lwd=2,lty=4)
	points(NWe2$x,CIlb,type = "l",col="blue",lwd=2,lty=4)
	legend("topright", c("NW estimator",expression("CI"[99])), lty=c(1,4), col=c("red","blue"),bty = "n")


}

ex5 = function(pressure,tritium,n){

    hCV = 0.26
    
    llest = WLS(pressure,tritium,'Gaussian',hCV,2)[[1]]
	x = WLS(pressure,tritium,'Gaussian',hCV,2)[[2]]

	h0 = optim(par=0.2, fn=function(h) CV(pressure,tritium,'normal',h), method="BFGS")
    hCV1 = h0$par
    NWe = NWe(pressure,tritium,hCV1,'normal')

	plot(pressure,tritium,col="gray50", ylab="Tritium", xlab="Pressure", main="Regression of tritium on pressure using Local polynomial estimates")
	points(x,llest,type = "l",col="red",lwd=2)
	points(NWe$x,NWe$y,type="l",col="blue",lwd=2)
    legend("topright", c("Local polynomial estimator (p=2)","NW estimator"), lty=c(1,1), col=c("red","blue"),bty = "n")

}

ex6 = function(tritium,longitude,latitude,n){

	y = tritium
	x1 = longitude
	x2 = latitude

	# Fit an a thin plate splines #

	Est = gam(y~s(x1,x2))
	x1val = grid(x1)
	x2val = grid(x2)

	xval = expand.grid(x1=x1val, x2=x2val)
	yhat = predict(Est, newdata = xval)

	ymat = matrix(yhat,n,n)

	contour(x = x1val, y = x2val, z = ymat, xlab ="longitude", ylab = "latitude")
	points(x1,x2,col=2,pch=20,cex=0.5)

	print('Question 6')

}

ex7 = function(tritium,longitude,latitude,n){

	print('Question 7 is not a programming exercise')

}

ex8 = function(tritium,longitude,latitude,salinity,pressure,n){

	y = tritium
	x1 = longitude
	x2 = latitude
	x3 = salinity
	x4 = pressure

    print('Question 8')

	# Fit an additive model #
	AdditiveEst = gam(y~s(x1)+s(x2)+s(x3)+s(x4))

	par(mfrow=c(2,2))
	plot(AdditiveEst)

}

main = function(){
    
    # filepath #
	filepath ='tritium.csv'
	data = readTheData(filepath)
    
    # Variables #
	tritium = data$tritium
	pressure = data$pressure
	longitude = data$longitude
	latitude = data$latitude
	salinity = data$salinity
    
    # vector length #
	n = length(tritium)
    
    # Run the exercises #
    # Part 1a #
	ex1(tritium,n)
	ex2(pressure,tritium,n)
	ex3(pressure,tritium,n)
	ex4(pressure,tritium,n)
	ex5(pressure,tritium,n)

	# Part 1b #
	ex6(tritium,longitude,latitude,n)
	ex7(tritium,longitude,latitude,n)
	ex8(tritium,longitude,latitude,salinity,pressure,n)

}

main()
