#### function used


##### SIR model - this code was slightly modified and adapted from the script of Sherry Towers (informations and copyright below)
##################################################################################
# An R script to solve ODE's of an SIR model with age structure
# Further info at http://sherrytowers.com/2012/12/11/sir-model-with-age-classes/
#
# Author: Sherry Towers
#         smtowers@asu.edu
# Created: Dec 2nd, 2012
#
# Copyright Sherry Towers, 2012
#
# This script is not guaranteed to be free of bugs and/or errors.
#
# This script can be freely used and shared as long as the author and 
# copyright information in this header remain intact.
##################################################################################
##################################################################################
# this is an age structured SIR model
# the parameters in the vparameters list are:
#    the recovery period, gamma
#    the probability of transmission on contact, beta
#    the contact matrix, C, that is the # contacts per day among age groups
#
# Note that x is a vector of length (#model compartment types)*(#age classes)
# For the SIR model, there are 3 model compartment types (S, I, and R)
# The code at the beginning of the function fills the age classes for each
# model compartment type in turn.
# Thus, S, I and R are vectors, all of length nage
##################################################################################
calculate_derivatives <- function(t, x, vparameters){
	ncompartment <- 3
	nage <- length(x)/ncompartment
	S    <- as.matrix(x[1:nage])
	I    <- as.matrix(x[(nage+1):(2*nage)])
	R    <- as.matrix(x[(2*nage+1):(3*nage)])
	
	I[I<0] = 0
	with(vparameters,{
		# note that because S, I and R are all vectors of length nage, so will N,
		# and dS, dI, and dR
		N = S+I+R
		dS = -as.matrix(S*beta)*(as.matrix(C)%*%as.matrix(I/N))
		dI = -dS - gamma*as.matrix(I)
		dR = +gamma*as.matrix(I)
		# remember that you have to have the output in the same order as the model
		# compartments are at the beginning of the function
		out=c(dS,dI,dR)
		list(out)
	})
}


##### Function that calculate initial parameters for the model
##################################################################################
# needs country name (as in the "Table Name" column of the Country.xlsx) as argument
# needs (beta OR R0) AND gamma as arguments
# if beta is known, R0 is calculated
# if R0 is known, beta is calculated 
parameters_SIR_COVID <- function(countryname, beta = 0.01270518, R0 = NULL, gamma = 1/14, ages=age_cl){

	# search and load the contact matrix
	if (countryname %in% data_C_1 ){
		C <- CI_act_df <- read_excel("MUestimates_all_locations_1.xlsx", sheet = countryname) %>% as.matrix()
	} else if (countryname %in% data_C_2){
		C <- CI_act_df <- read_excel("MUestimates_all_locations_2.xlsx", sheet = countryname, col_names = FALSE) %>% as.matrix()
	}
	row.names(C) <- colnames(C) <- ages
	# load population data
	data_pop <- data_pop_all %>%	filter(`Country Name`==countryname)
	
	# select usefull data (population per age classes)
	pop <- data_pop  %>%
		filter(str_detect(`Indicator Code`,"SP\\.POP\\.[[:digit:]]{2}[[:alnum:]]{2}\\.(FE|MA)$")) %>%
		separate(`Indicator Code`, c("ind1", "ind2","Age", "Sexe")) %>%
		group_by(Age) %>%
		summarise(pop = sum(`2018`)) %>%
		ungroup()
	
	# sum age classes 75-79 and 80+ into one classe 75+
	new_line <- pop %>% mutate(sel = c(rep(1,15),0,0)) %>%
		group_by(sel) %>%
		summarise(sum=sum(pop)) %>%
		filter(sel==0)
	pop <- pop %>% filter(!(Age %in% c("7579","80UP"))) %>% rbind(c("75UP",as.numeric(new_line[1,2])))
	pop$pop <- as.numeric(pop$pop)
	
	
	# 
	f <- pop$pop / sum(pop$pop) # frequency per age classe
	N <- pop$pop # Number of individuals per age classe
	nage = length(f) # Number of age classes
	
	# Calculate beta or R0 according to Towers et al. 2012
	ex <- expand.grid(1:nage,1:nage) 
	M <-map2(ex$Var1, ex$Var2, .f =fun <- function(i,j){C[i,j]*f[i]/f[j]}) %>% unlist() %>% matrix(nrow = nage)
	eig <- eigen(M)
	
	
	if (is.null(beta) & !(is.null(R0)))	{
		
		beta = R0*gamma/max(Re(eig$values)) 
		
	} else if (is.null(R0) & !(is.null(beta)))	{
		
		R0 <- beta / gamma * max(Re(eig$values))
		
	} else {
		print("one argument among beta and R0 should be filled, the other should be set to NULL")
	}
	
	# initial population (N individuals in each compartment and age classe)
	I_0    = rep(0,nage) # vector of Infected
	I_0[9] <- 1				 # put one infected person in 60-64 age classe
	S_0    = N-I_0				# vector of Susceptible
	R_0    = rep(0,nage) # vector of Recovered
	
	# return parameters
	return(paramSIR = list(vparameters = list(gamma=gamma,beta=beta,C=C), inits = c(S=S_0,I=I_0,R=R_0), R0 = R0))
}


####
lsoda2 <- function(param, vt){
	model <- lsoda(param$inits, vt, calculate_derivatives, param$vparameters) %>% as.data.frame()
	return(model)
}

##### Function that calculate R0 from a vector of daily new cases (from the first to the last case)
##################################################################################
# the fonction used the estimate.R function of the R0 library
# take as argument a vector of daily new cases (x) and a Generation Time repartition function (GT, see help of the estimate.R function) 
# the function remove leading and ending zeros in the vector of daily new cases
# return NA in case when function estimate.R fails.
R0x <- function(x){
	x <- x[cumsum(x) & rev(cumsum(rev(x)))] # remove leading and ending zeros
	tryCatch(R0 <- estimate.R(epid = x, GT = GT, begin = as.integer("1"), end = as.integer(length(x)), methods=c("ML")), error = function(e) NULL)
	if(is.null(return(R0$estimates[[1]]$R))){
		return(NA) } 
	else {
		return(R0$estimates[[1]]$R)
	}
}
