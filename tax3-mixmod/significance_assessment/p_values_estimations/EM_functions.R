#########################################
#					
#    EM ALGORITHM : BETA MIXTURE	
#					
#########################################


# Required libraries
# as root: install.packages(gtools)
# The library gtools allows the generation of random mixing proportions.
library(gtools)

# LOG
verbose=TRUE

# TODO : option to start with mirrored random values


#======================
# Internal functions :
#======================

# MOMENT :
# This function estimates the parameters alpha and beta
# for a beta mixture by using the method of moments. 
# The method of moments provides an efficient estimation 
# of the parameters alpha and beta.
# x = vector of resampling
# t = vector of conditional probability
# sum_t = sum of conditional probability vector
#
moment<-function( x , t , sum_t)
{
	w<-t/sum_t # weight
	x_bar<-sum(w*x, na.rm=TRUE) # mean
	s2=sum( w*(x-x_bar)^2, na.rm=TRUE ) # variance
	if(s2==0) return(NaN)
	z<-(x_bar*(1-x_bar)/s2-1)

	# return alpha0 & beta0
	return(c(x_bar*z,(1-x_bar)*z))
}

#===================================================================================
# GRADIENT :
# This function calculates the gradient of the log-likelihood
# for 1 Beta component. 
# x = vector of resampling
# t = vector of conditional probability
# sum_t = sum of vector of conditional probability
#
g<-function( x, t, sum_t, alpha, beta )
{
	# get constants
	z<-digamma(alpha+beta)

	return( c( ( sum_t*( z - digamma(alpha) ) + sum( t*log(x),   na.rm=TRUE ) ),
	           ( sum_t*( z - digamma(beta)  ) + sum( t*log(1-x), na.rm=TRUE ) ) ) )
}

#===================================================================================
# JACOBIAN :
# This function returns the Jacobian matrix of the gradient
# for a given component k.
# sum_t = sum of vector of conditional probability
#
Jacob<-function( sum_t, alpha, beta )
{
	A<-matrix(NA,2,2)

	# get constants
	z<-trigamma(alpha+beta)

	A[1,1]<-sum_t*( z - trigamma(alpha) ) 
	A[2,2]<-sum_t*( z - trigamma(beta)  ) 
	A[1,2]<-A[2,1]<-sum_t*z
	return(A)
}

#===================================================================================
# NEWTON :
# This function returns the solution of the numerical 
# Newton algorithm to search the zeros of the gradient g
# for a given component k.
# x = vector of resampling
# t = vector of conditional probability
#
# Global Parameters for the newton-raphson algorithm :
# Tolerance was 1e-10 : Christophe Biernacki okayed lowering to 1e-2 or 1e-3, as not critical
# 10-2 gives wrong results (one beta)
tolerance <- 1e-3 
iteration_max <- 1000 

newton<-function(x , t, alpha, beta )
{
	# number of iterations
	count=0

	# get sum of vector t
	sum_t<-sum(t, na.rm=TRUE)

	# start values
	theta<-moment( x, t, sum_t )
	if( is.nan(sum(theta)) ) {return( NULL) }

	b<-g( x, t, sum_t, theta[1], theta[2] )

	# Loop : until the error is negligible
	while( (max(abs(b)) > tolerance) && (count < iteration_max) )
	{
		# Jacobian matrix
		A<--Jacob( sum_t, theta[1], theta[2] )

		# diagonal precondition number
		D_inv<-solve( diag(diag(A)) )
	
		# equation solver 
      # 1e-7 is default
		sol<-qr.solve( D_inv%*%A, D_inv%*%b, tol=1e-7 ) 
		
		# update the solution	
		theta<-sol+theta

		# update the error
		b<-g( x, t, sum_t, theta[1], theta[2] )
		if( is.infinite(sum(b)) | is.nan(sum(b)) ) {cat("\nb=",b,"\n") ;return( NULL) }

		# increase counter
		count=count+1
	}

   # OD
   if (verbose) {
      print(paste("     Newton:  count=",count,"  error=",max(abs(b)), sep=""),quote = FALSE)
   }

	return(theta)
}

#===================================================================================
# LOG MIXTURE DENSITY:
# x_i = quantiles i
# p = vector of proportions
# alpha = vector of parameters
# beta = vector of parameters
#
log_mixture<-function(x_i, p, alpha, beta)
{ return(log(p)+dbeta(x_i,alpha,beta,log=TRUE)) }

#===================================================================================
# MIXTURE DENSITY:
# x_i = quantiles i
# p = vector of proportions
# alpha = vector of parameters
# beta = vector of parameters
#
mixture<-function(x_i, p, alpha, beta)
{ return(p*dbeta(x_i,alpha,beta)) }


####################################################################################
####################################################################################

# SHORT_EM FUNCTION
#===================================================================================
# Run EM_algorithm with random position until epsilon=0.01 
# Return parameters and Likelihood if success
# If failed, then return -Inf for Likelihood and 0 for parameters
SHORT_EM<-function(x, K, epsilon=1e-02, nbiteration=200, model="Beta_free_proportions")
{
   # log
   if (verbose) {
   	cat("Short EM algorithm\n" )
   	cat("------------------\n" )
      cat("  - K       = ",K,"\n")
      cat("  - epsilon = ",epsilon,"\n\n")
   }

	# Dimension of the dataset
	n<-length(x) 

	# ---------------------------------------------------
	# find suitable starting values, with no empty class
	flag <- TRUE
	counter <- 0
	ClassSize <- rep(0,K)
	alpha <- rep(1,K)
	beta <- rep(1,K)
	while ( ( sum( ClassSize < 2 ) != 0 ) && (counter < nbiteration) ) {
		counter=counter+1
		# if all classes empty (first time or after reset)
		# get new mixing proportions
		if ( sum(ClassSize == 0) == K ) { 
			if (model == "Beta_free_proportions")	{ 
				proportions<-rdirichlet(K,rep(1,K))[1,] 
			}
			if (model == "Beta_equal_proportions")	{ 
				proportions<-rep((1/K),K) 
			}
		}
		# get new parameters for empty classes
		for(i in 1:K){
			if (ClassSize[i] < 2) {
				# Beta function parameters where most of the time alpha<beta
				alpha[i] <-runif(1,1,1000)
				if ( runif(1,0,1) > 0.1 ) {
					beta[i] <- runif(1, alpha[i] , 1000)		
				} else {
					beta[i] <- runif(1, 1 , 1000)		
				}			
			}
		}
		
		# get density & condition probabilities for ClassSize		
		if(K==1) mixture_density<-as.matrix(apply(as.matrix(x),1,mixture,p=proportions,alpha=alpha,beta=beta))
		else     mixture_density<-t(apply(as.matrix(x),1,mixture,p=proportions,alpha=alpha,beta=beta))
		conditional_probability<-mixture_density/apply(mixture_density,1,sum,na.rm=TRUE)
		ClassSize <- apply(conditional_probability,2,sum,na.rm=TRUE)
		
		# mixing proportions reset criteria: one class full
		if ( ( K != 1) && ( sum(ClassSize == n) != 0 ) ) {
			ClassSize <- rep(0,K)
		}
		
	}
	
	# log
	if (verbose) {
		cat("Init parameters\n" )
		cat("  - tries        : ",counter,"\n" )
		cat("  - ClassSize    : ",formatC(floor(ClassSize),width=8),"\n" )
		cat("  - alpha        : ",formatC(alpha,format="f",digits=1,width=8),"\n" )
		cat("  - beta         : ",formatC(beta,format="f",digits=1,width=8),"\n" )
		cat("  - proportions  : ",formatC(proportions,format="f",digits=3,width=8),"\n" )
	}
	
	# ---------------------------------------------------
	# EM ALGORITHM
   #
	
	# Initial mixture density
	if(K==1)	{
		mixture_density<-as.matrix(apply(as.matrix(x),1,mixture,p=proportions,alpha=alpha,beta=beta))
	}	else	{
		mixture_density<-t(apply(as.matrix(x),1,mixture,p=proportions,alpha=alpha,beta=beta))
	}

	# starting parameters
	parameters=1:(3*K)
	index=1:(3*K)
	loglikelihood=-Inf
	abs_difference=1   
	counter=1

	# main loop
	while ( (abs_difference > epsilon) && (counter < iteration_max) )
	{
      # ------------
		# E step 

		# Conditional probabilities
		conditional_probability<-mixture_density/apply(mixture_density,1,sum,na.rm=TRUE)

		# Stopping conditions : stop EM algorithm if one class is too small: '<2' threshold suggested by C.Biernacki
		if ( sum(apply(conditional_probability,2,sum,na.rm=TRUE)<2) ) 	{
			return(c(-Inf, rep(0,(3*K)) ))
		}
	
      # ------------
		# M step 
		if (model == "Beta_free_proportions")	{ 
			proportions<-apply(conditional_probability,2,sum,na.rm=TRUE)/n 
		}

		# Maximization of the observed Log-Likelihood
		for (k in 1:K)		{
			theta<-newton(x, conditional_probability[,k], alpha[k], beta[k])
			if(is.null(theta)) {
				return(c(-Inf, rep(0,(3*K)) ))
			}
			alpha[k]=theta[1]
			beta[k]=theta[2]
		}
	
		# Update mixture density
		if(K==1)	{
			mixture_density<-as.matrix(apply(as.matrix(x),1,mixture,p=proportions,alpha=alpha,beta=beta))
		}	else	{
			mixture_density<-t(apply(as.matrix(x),1,mixture,p=proportions,alpha=alpha,beta=beta))
		}

		# Observed Log-likelihood
		old_loglikelihood<-loglikelihood
      loglikelihood<-sum( log( apply( mixture_density, 1, sum, na.rm=TRUE ) ) )

      # ------------
		# Stopping conditions	
		#
		# if loglikelihood is not finite, stop EM algorithm
		if( is.infinite(loglikelihood) | is.nan(loglikelihood))	{
			return(c(-Inf, rep(0,(3*K))))
		}
      # ------------

		# Absolute difference between 2 consecutive values of the loglikelihood.
		abs_difference<-abs(loglikelihood - old_loglikelihood)

      # log
      if (verbose) {
        cat("  counter=",counter,"   loglikelihood=",loglikelihood,"\n");
      }

		# Increasing
		counter=counter+1
	}
	# END OF ALGORITHM

   # log
   if (verbose) {
     cat("Short_EM completed.\n\n");
   }

	# EM is done : get parameters
	parameters[which(index%%3 == 1)]=proportions
	parameters[which(index%%3 == 2)]=alpha
	parameters[which(index%%3 == 0)]=beta
	return(c(loglikelihood, parameters))
}



# SMALL_EM FUNCTION
#===================================================================================
# Run EM_algorithm with 10 initial random positions. 
# Keep the best one according to the observed Log-likelihood of datas
# Return parameters of this best results
# If the 10 positions failed, then return random parameters
SMALL_EM<-function(x, K, epsilon=1e-02, nbiteration=200, model="Beta_free_proportions", n.positions=10)
{
	# Dimension of the dataset
	n<-length(x) 

	all_parameters=matrix(0, ncol=n.positions, nrow=(3*K))
	loglikelihood=1:n.positions
	index=1:(3*K)

	# loop on the positions : get all parameters and all Likelihoods
	for (index.position in 1:n.positions)
	{
		short_em=SHORT_EM(x, K, epsilon=1e-02, nbiteration=200, model="Beta_free_proportions")
		loglikelihood[index.position]=short_em[1]
		all_parameters[, index.position]=short_em[2:length(short_em)]
	}

	# get the best results of this 10 short_EM :
	# check before if at least 1 of them succeed
	res=sum(loglikelihood == rep(-Inf, n.positions))
	if ( res != n.positions)	{	
		index_best_results=which(loglikelihood==max(loglikelihood))
		best_parameters=all_parameters[ ,index_best_results]
	}	else	{
		# SMALL_EM failed 
		best_parameters <- NULL
	}
	return(best_parameters)
}


# EM ALGORITHM FOR MIXTURE DENSITY :
#===================================================================================
# _ This function estimates a pdf by using Beta Mixture. 
# _ This is an iterative algorithm. 
# _ The threshold epsilon must be defined by the user. 
# _ x is the dataset and K is the number of components. 
# _ model indicates if the proportions are free or not.
EM_algorithm<-function(x, K, epsilon=1e-06, nbiteration=200, model="Beta_free_proportions", output.path="./", init.file="none", SMALL_EM=TRUE)
{
	# Dimension of the dataset
	n<-length(x) 

	index=1:(3*K)	

	#
	# Starting parameters 
	#
	if (init.file != "none")	{
	
		# ---------------------------------------------------
		# user defined starting parameters
		
		cat("\nEM with init.file is running...\n")
		
		init.file.parameters=as.matrix(read.table(init.file), header=FALSE)
		proportions=init.file.parameters[which(index%%3 == 1),1]
		alpha=init.file.parameters[which(index%%3 == 2),1]
		beta=init.file.parameters[which(index%%3 == 0),1]
		
	} else {

		init.parameters <- NULL

		if(SMALL_EM == TRUE)	{

			# ---------------------------------------------------
			# starting parameters derived from SMALL_EM run
			# (10 positions of short_runs_EM )

			cat("\nSMALL_EM is running...\n\n")
			init.parameters=SMALL_EM(x, K, epsilon=1e-02, model="Beta_free_proportions", n.positions=10)

		}


		if ( ! is.null(init.parameters) ) {

			# ---------------------------------------------------
			# SMALL_EM sucess : get parameters from Small_EM output
			
			init.parameters <- as.matrix(init.parameters)
			proportions		 <- init.parameters[which(index%%3 == 1),1]
			alpha 			 <- init.parameters[which(index%%3 == 2),1]
			beta 				 <- init.parameters[which(index%%3 == 0),1]
			
		}	else	{

			# ---------------------------------------------------
			# either SMALL_EM failed or USER ask for random parameters
			
			cat("\nRandom start values...\n")		
			
			# ---------------------------------------------------
			# find suitable starting values, with no empty class
			# (cette partie pourrait être remplacée par la méthode des moments)
			
			counter = 0
			ClassSize <- rep(0,K)
			alpha <- rep(1,K)
			beta <- rep(1,K)
			
			while ( ( sum( ClassSize < 2 ) != 0 ) && (counter < nbiteration) ) {
				counter=counter+1
			
				# if all classes empty (first time or after a reset)
				# get new mixing proportions
				if ( sum(ClassSize == 0) == K ) { 
					if (model == "Beta_free_proportions")	{ 
						proportions<-rdirichlet(K,rep(1,K))[1,] 
					}
					if (model == "Beta_equal_proportions")	{ 
						proportions<-rep((1/K),K) 
					}
				}
			
				# get new parameters for empty classes
				for(i in 1:K){
					if (ClassSize[i] < 2) {
						# Beta function parameters where most of the time alpha<beta
						alpha[i] <-runif(1,1,1000)
						if ( runif(1,0,1) > 0.1 ) {
							beta[i] <- runif(1, alpha[i] , 1000)		
						} else {
							beta[i] <- runif(1, 1 , 1000)		
						}			
					}
				}
	
				# get density & condition probabilities for ClassSize		
				if(K==1) mixture_density<-as.matrix(apply( as.matrix(x),1,mixture,p=proportions,alpha=alpha,beta=beta ))
				else     mixture_density<-t(apply( as.matrix(x),1,mixture,p=proportions,alpha=alpha,beta=beta ))
				conditional_probability<-mixture_density/apply(mixture_density,1,sum,na.rm=TRUE)
				ClassSize <- apply(conditional_probability,2,sum,na.rm=TRUE)
	
				# mixing proportions reset criteria: one class full
				if ( ( K != 1) && ( sum(ClassSize == n) != 0 ) ) {
					ClassSize <- rep(0,K)
				}		
			} # end while
			
			# log
			if (verbose) {
				cat("Random init parameters\n" )
				cat("  - tries        : ",counter,"\n" )
				cat("  - ClassSize    : ",formatC(floor(ClassSize),width=8),"\n" )
				cat("  - alpha        : ",formatC(alpha,format="f",digits=1,width=8),"\n" )
				cat("  - beta         : ",formatC(beta,format="f",digits=1,width=8),"\n" )
				cat("  - proportions  : ",formatC(proportions,format="f",digits=3,width=8),"\n\n" )
			}
		}
	}

	# Print out starting parameters
	cat("----------------------------\n")
	cat("Starting parameters:\n")
	cat("K           =",format(seq(1:K),digit=2,nsmall=2,justify="centre",width=8),"\n")
	cat("proportions =",format(proportions,digit=2,nsmall=2,justify="centre",width=8),"\n")
	cat("alpha       =",format(alpha,digit=2,nsmall=2,justify="centre",width=8),"\n")
	cat("beta        =",format(beta,digit=2,nsmall=2,justify="centre",width=8),"\n")

	# ---------------------------------------------------
	# EM ALGORITHM
   #

	# Initialization of algorithm parameters
	loglikelihood=-Inf		# initial LogL
	abs_difference=+Inf		# LogL difference between iterations
	counter=1					# number of iterations

	# Initial mixture density
	# 'apply’ returns an array of dimension ‘c(K, dim(x)[MARGIN])’ if ‘n > 1’.
	# If ‘K’ equals ‘1’, ‘apply’ returns a vector
	if(K==1) mixture_density<-as.matrix(apply( as.matrix(x),1,mixture,p=proportions,alpha=alpha,beta=beta ))
	else     mixture_density<-t(apply( as.matrix(x),1,mixture,p=proportions,alpha=alpha,beta=beta ))

	# main loop
	while ( (abs_difference > epsilon) && (counter < iteration_max) )
	{
		# ----------
		# E step 
		#

		# Conditional probabilities
		conditional_probability<-mixture_density/apply(mixture_density,1,sum,na.rm=TRUE)

		# Stopping conditions : stop EM algorithm if one class is too small: '<2' threshold suggested by C.Biernacki
      ClassSize <- apply(conditional_probability,2,sum,na.rm=TRUE)
		if (verbose) {
			cat("ClassSize=",floor(ClassSize),"\n" ) 
		}
		if ( sum(ClassSize<2) ) {
			return(cat("\nError: class too small.\n"))
		}

		# ----------
		# M step 
		#

		# Update proportions (only if they are free)
		if (model == "Beta_free_proportions") { 
            proportions<-apply(conditional_probability,2,sum,na.rm=TRUE)/n 
      }
	
		# Maximization of the observed Log-Likelihood
		for (k in 1:K)	{
			theta<-newton(x, conditional_probability[,k], alpha[k], beta[k])
			if(is.null(theta)) {
				return(cat("\nError: Newton failed\n") )
			}
			alpha[k]=theta[1]
			beta[k]=theta[2]
		}
		
		# Updated mixture density with new parameters
		if(K==1) mixture_density<-as.matrix(apply( as.matrix(x),1,mixture,p=proportions,alpha=alpha,beta=beta ))
		else     mixture_density<-t(apply( as.matrix(x),1,mixture,p=proportions,alpha=alpha,beta=beta ))

		# Observed Log-likelihood
		old_loglikelihood<-loglikelihood
      loglikelihood<-sum( log( apply( mixture_density, 1, sum, na.rm=TRUE ) ) )

		# if loglikelihood is not finite, stop EM algorithm
		if( is.infinite(loglikelihood) | is.nan(loglikelihood) ) return( cat("\nError: loglikelihood is not finite.\n") )

		# Absolute difference between 2 consecutive values of the loglikelihood.
		abs_difference<-abs( loglikelihood - old_loglikelihood )

		# log
		if (verbose) {
			print(paste("  counter=",counter,"   loglikelihood=",loglikelihood,sep=""),quote = FALSE);
		}

		# Increment counter
		counter=counter+1
	}
	# END OF ALGORITHM
	
	# BIC criterion (corrected as per C.Biernacki's email)
	if (model == "Beta_free_proportions")  BIC=-2*loglikelihood + (3*K-1)*log(n)
	if (model == "Beta_equal_proportions") BIC=-2*loglikelihood + 2*K*log(n)
	
	# Parameters
	parameters=matrix(data=c(0), nrow=(3*K), ncol=1)
	index=1:(3*K)
	parameters[which(index%%3 == 1),1]=proportions
	parameters[which(index%%3 == 2),1]=alpha
	parameters[which(index%%3 == 0),1]=beta

	# RESULTS
	write.table(BIC, file=paste(output.path,"BIC.txt",sep=""), col.names=FALSE, row.names=FALSE, append=FALSE, quote=FALSE)
	write.table(parameters, file=paste(output.path,"parameters.txt",sep=""), col.names=FALSE, row.names=FALSE, append=FALSE, eol="\n", quote=FALSE)
	cat("----------------------------")
	cat("\nBIC=",BIC,"\n")
	cat("----------------------------")
	cat("\nFinal parameters:\n")
	cat("K           =",format(seq(1:K),digit=2,nsmall=2,justify="centre",width=8),"\n")
	cat("proportions =",format(proportions,digit=2,nsmall=2,justify="centre",width=8),"\n")
	cat("alpha       =",format(alpha,digit=2,nsmall=2,justify="centre",width=8),"\n")
	cat("beta        =",format(beta,digit=2,nsmall=2,justify="centre",width=8),"\n")
	cat("----------------------------\n")
	return( cat("\nEM algorithm done!\nResults into '",output.path,"' directory.\n") )
}


#
# Plotting functions
#
#===================================================================================
# PLOT MIXTURE :
# This function plot the Mixture density 
# according to the parameters given by the EM
# algorithm and the number of components.
plot_mixture<-function(path,x,K,parameters)
{
	n=length(x)
	abscisse=seq(0,1,1/(n-1))
	graph=array(data=c(0),dim=length(abscisse))
	graph=as.vector(graph)

	index=1:(3*K)
	proportions=parameters[which(index%%3 == 1),1]
	alpha=parameters[which(index%%3 == 2),1]
	beta=parameters[which(index%%3 == 0),1]

	for (k in 1:K)
	{
		graph = graph + proportions[k]*dbeta(abscisse,alpha[k],beta[k])
	}

	png(path)
	leg.txt=c("EM beta mixture")
	hist(x, freq=FALSE, main=sprintf("BETA MIXTURE WITH THE E.M. ALGORITHM : %d COMPONENTS",K),
	     xlab="", ylab="Mixture Density", xlim=c(0,1), col="green", border="darkgreen")
	lines(abscisse, graph, col="darkblue", lwd=2, type="l", lty=1)
	legend("topright", legend=leg.txt, col="darkblue", lwd=2, cex=0.8, lty=1)
	dev.off()
}

#===================================================================================
# PLOT INCOMPLETE LOG-LIKELIHOOD :
# This function plot the uncomplete log likelihood
# to verify that it is an increasing function.
plot_likelihood<-function(path,loglikelihood)
{
	n=length(loglikelihood)
	abscisse=1:n
	png(path)
	leg.txt=c("Uncomplete Log-likelihood")
	plot(abscisse, loglikelihood, col="blue", lwd=2, type="l", lty=1, main="OBSERVED COMPLETE LOG-LIKELIHOOD")
	legend("bottomright", legend=leg.txt, col="blue", lwd=2, cex=1.0, lty=1)
	dev.off()
}

