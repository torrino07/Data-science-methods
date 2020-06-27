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
##   20/05/2020
########################################
###############
## Libraries ##
library(mgcv)
library(glmnet)
library(nnet)
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
    m = length(x[1,])
    z = x

	for(i in 1:m){
		z[,i] = (x[,i]-mean(x[,i]))/sd(x[,i])
	}

	return(z)
}

logistic = function(x){
########################################
#
# Purpose: create the logistic function
#   
# Inputs: x = variable of interest
#   
# Outputs: logistic function 
#   
########################################

	logistic = exp(x)/(1+exp(x))

	return(logistic)

}

CI = function(z,q,SE){
########################################
#
# Purpose: define confidence intervals for Nadaraya-Watson estimator
#   
# Inputs: mu = prediction, alpha lv, SE = standard error
#   
# Outputs: confidence intervals
#   
########################################

	CIub = z + qnorm(q)*SE
	CIlb = z - qnorm(q)*SE

	CI = list(CIub,CIlb)

	return(CI)

}

Hitratio = function(Y,pred,treshold){
########################################
#
# Purpose: count the percentage of correct predictions
#   
# Inputs: Y= actual values , pred = latent prediction variable, treshhold
#   
# Outputs: hit-ratio
#   
########################################
	n = length(pred)

	for (j in 1:n){
		if (pred[j] > treshold){
			pred[j]=1
		}
		else{
			pred[j]=0

		}
	}

	yhat = pred
    
	correctNrPred = 0
	for (k in 1:n){
		if(yhat[k]==Y[k]){
			correctNrPred = correctNrPred + 1

		}
	}

	Hitratio = correctNrPred/n

	return(Hitratio)

}


# Part 1a #
Ex1a = function(quality,alcohol){
    
	y = quality
	x = alcohol
	n = length(x)

	# Logistic parametric regression estimation
	est = gam(y~x,family = binomial)
	xval = grid(x)
	yhat = predict(est, data.frame(x=xval), se.fit = TRUE)
	z = yhat$fit
	logistic = logistic(z)

	# Plot
	plot(x,y,ylab = "Quality",xlab = "Alcohol",main="Parametric Logistic Regression - Quality on Alcohol")
	legend(13.6,0.15,c('LPR'), lty=c(7),col=c('red'), bty = "n")
	points(xval,logistic,type = "l",col=2,lwd=2)

}

Ex2a = function(quality,alcohol){

	y = quality
	x = alcohol
	n = length(x)

	# Logistic non-parametric regression estimation
	est = gam(y~s(x),family = binomial)
	xval = grid(x)
	yhat = predict(est, data.frame(x=xval), se.fit = TRUE)
	SE = yhat$se.fit
	z = yhat$fit
	logistic = logistic(z)

	# 99 % CI

	# Confidence Bands

	alpha = 0.01
	q = 1-alpha/2
	

	CIub = CI(z,q,SE)[[1]] # Upper Bound
	CIlb = CI(z,q,SE)[[2]] # Lower Bound

	UB = logistic(CIub) #logistic(.)
	LB = logistic(CIlb) #logistic(.)

	# Plot
	plot(x,y,ylab = "Quality",xlab = "Alcohol",main="Non - Parametric Logistic Regression - Quality on Alcohol")
	legend(8.6, 0.17,c('LPR','99% CI'), lty=c(7,2), col=c('red',2), bty = "n")
	points(xval,logistic,type = "l",col=2,lwd=2)
	points(xval,UB,type = "l",col=2,lwd=2,lty=2)
	points(xval,LB,type = "l",col=2,lwd=2,lty=2)
}

# Part 2a #

Ex1b = function(data){

	m = as.matrix(data)
	Xs = m[,-12]

	y = m[1:1000,12]
	X = Xs[1:1000,]

	# Neural Network
	est = nnet(X,y, size=2, decay=0, linout=F, maxit=2000)

	summary(est)

}

Ex2b = function(data,treshold){

	m = as.matrix(data)
	Xs = standardize(m[,-12])

	y = m[1:1000,12]
	X = Xs[1:1000,]

	# Training sample
	y_train = y[1:700]
	X_train = X[1:700,]

    # Validation sample (to choose λ)
	y_sel = y[701:1000]
	X_sel = X[701:1000,]


	a = seq(0.01,0.5, length.out = 20)
	out = rep(0,20) # Criterion function

	for(i in 1:20){
		est = nnet(X_train,y_train, size=50, decay=a[i], linout=F, maxit=2000)
		yhat = predict(est, newdata = X_sel)
		out[i] = mean(y_sel*log(yhat)+(1-y_sel)*log(1-yhat))
	}

	LambdaOpt = a[which(out==max(out))]
	print(sprintf("Optimal Decay Parameter λ  = %g ",LambdaOpt))

    
	# Probability of correct classification #
	# Neural Network
	est = nnet(X,y, size=2, decay=LambdaOpt, linout=F, maxit=2000)

	# Prediction #

	X = Xs[1001:1599,]
	Y = m[1001:1599,12]

	pred = predict(est,X)

	Hitratio = Hitratio(Y,pred,treshold)

	print(sprintf("Hit-ratio NN = %g ", Hitratio))

}


Ex3b = function(data,treshold){

	m = as.matrix(data)
	Xs = m[,-12]

	y = m[1:1000,12]
	X = Xs[1:1000,]

	x1 = X[,1]
	x2 = X[,2]
	x3 = X[,3]
	x4 = X[,4]
	x5 = X[,5]
	x6 = X[,6]
	x7 = X[,7]
	x8 = X[,8]
	x9 = X[,9]
	x10 = X[,10]
	x11 = X[,11]

	## Non-Parametric logistic additive model ##

	estnpa = gam(y~s(x1)+s(x2)+s(x3)+s(x4)+s(x5)+s(x6)+s(x7)+s(x8)+s(x9)+s(x10)+s(x11), family=binomial)
	estpa = gam(y~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11, family=binomial)


	X = Xs[1001:1599,]
	Y = m[1001:1599,12]
	

	theta = predict(estnpa, newdata = data.frame(x1=X[,1],x2=X[,2],x3=X[,3],x3=X[,3],x4=X[,4],x5=X[,5],x6=X[,6],x7=X[,7],x8=X[,8],x9=X[,9],x10=X[,10],x11=X[,11]))
	ypredNPA = logistic(theta)

	HitratioNPA = Hitratio(Y,ypredNPA,treshold)

	print(sprintf("Hit-ratio Non-Parametric model = %g ",HitratioNPA))


	## Parametric logistic additive model ##
	beta = predict(estpa, newdata = data.frame(x1=X[,1],x2=X[,2],x3=X[,3],x3=X[,3],x4=X[,4],x5=X[,5],x6=X[,6],x7=X[,7],x8=X[,8],x9=X[,9],x10=X[,10],x11=X[,11]))
	ypredPA = logistic(beta)

	HitratioPA = Hitratio(Y,ypredPA,treshold)

	print(sprintf("Hit-ratio Parametric model = %g", HitratioPA))


}

Ex4b = function(data){
	

}

main = function(){

	# Magic number/data #
	treshold = 0.5

	filepath = 'red_wine.csv'
	data = readTheData(filepath)

	# Initialization #
	quality = data$quality
	alcohol = data$alcohol

    # Estimation #
	Ex1a(quality,alcohol)
	Ex2a(quality,alcohol)

	Ex1b(data)
	Ex2b(data,treshold)
	Ex3b(data,treshold)
	Ex4b(data)

}

main()