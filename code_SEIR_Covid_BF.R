#####
# N MOIROUX : nicolas.moiroux@ird.fr

# load required libraries
require("deSolve")
require(tidyverse)
require(readxl)
require(tmap)
require(sf)


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
calculate_derivatives=function(t, x, vparameters){
	ncompartment = 3
	nage = length(x)/ncompartment
	S    = as.matrix(x[1:nage])
	I    = as.matrix(x[(nage+1):(2*nage)])
	R    = as.matrix(x[(2*nage+1):(3*nage)])
	
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

##################################################################################
##################################################################################
data_pop_all <- read_excel("HNP_StatsEXCEL.xlsx") %>% as.data.frame() # data from World Bank
data_C_1 <- excel_sheets("MUestimates_all_locations_1.xlsx") # data from https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005697
data_C_2 <- excel_sheets("MUestimates_all_locations_2.xlsx")

#### function that calculate initial parameters for the model, needs country name (as in the "Table Name" column of the Country.xlsx)
## needs beta OR R0 AND gamma
## if beta is known, R0 is calculated
## if R0 is known, beta is calculated 
parameters_SIR_COVID <- function(countryname, beta = 0.01270518, R0 = NULL, gamma = 1/14){
	
	# search and load the contact matrix
	if (countryname %in% data_C_1 ){
		C <- CI_act_df <- read_excel("MUestimates_all_locations_1.xlsx", sheet = countryname) %>% as.matrix()
	} else if (countryname %in% data_C_2){
		C <- CI_act_df <- read_excel("MUestimates_all_locations_2.xlsx", sheet = countryname, col_names = FALSE) %>% as.matrix()
	}
	
	# load population data
	data_pop <- data_pop_all %>%	filter(`Country Name`==countryname)
	
	# select usefull data (population per age classes)
	pop <- data_pop  %>%
		separate(`Indicator Code`, c("ind1", "ind2","Age", "Sexe", "Sy")) %>% 
		filter(ind1 == "SP", ind2 == "POP", Sexe %in% c("MA","FE"), is.na(Sy)) %>%
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
	
	# Calculate beta / R0 according to Towers et al. 2012
	ex <- expand.grid(1:nage,1:nage) 
	M <-map2(ex$Var1, ex$Var2, .f =fun <- function(i,j){C[i,j]*f[i]/f[j]}) %>% unlist() %>% matrix(nrow = nage)
	eig <- eigen(M)

	 
	if (is.null(beta) & !(is.null(R0)))	{
		
		beta = R0*gamma/max(Re(eig$values)) 
		
	} else if (is.null(R0) & !(is.null(beta)))	{
		
		R0 <- beta / gamma * max(Re(eig$values))
		
	} else {
		print("one parameter beta or R0 should have a value, the other should be NULL")
	}
	
	# initial population (N individuals in each compartment and age classe)
	I_0    = rep(0,nage) # vector of Infected
	I_0[13] <- 1				 # put one infected person in 60-64 age classes
	S_0    = N-I_0				# vector of Susceptible
	R_0    = rep(0,nage) # vector of recovered
	
	# return parameters
	return(paramSIR = list(vparameters = list(gamma=gamma,beta=beta,C=C), inits = c(S=S_0,I=I_0,R=R_0), R0 = R0))
}


#### how much is beta (probability of transmission on contact) for france giving a nR0 of 2,5 (as calculated by ETE team)
## and a gamma (recovery rate) of 1/14 days-1
paramFR <- parameters_SIR_COVID("France", beta = NULL, R0 = 2.5, gamma = 1/14)
# ignore warnings
beta <- paramFR$vparameters[2]
beta


#### Calculate R0 for all countries based on the same probability of transmission on contact (beta)

#### dataframe of countries with population data available for 2018
liste_pays <- data_pop_all %>%  
	separate(`Indicator Code`, c("ind1", "ind2","Age", "Sexe", "Sy")) %>% 
	filter(ind1 == "SP", ind2 == "POP", Sexe %in% c("MA","FE"), is.na(Sy)) %>%
	group_by(`Country Name`) %>%
	summarise(n = sum(!is.na(`2018`))) %>%
	filter(n == 34) %>%
	ungroup() %>%
	rename(country=`Country Name`)

#### dataframe of countries with available contact matrix
liste_mat <- c(data_C_1, data_C_2)
liste_mat <- data.frame(country = liste_mat)

#### join the two dataframe
liste_pays <- left_join(liste_mat,liste_pays, by="country") %>%
	filter(n == 34) %>%
	mutate(R=NA) 

#### calculate R0 for all countries (146/152 having contact matrix data)
for (i in 1:nrow(liste_pays)){
	countryname <- liste_pays[i,1] %>% as.character()
	par <- parameters_SIR_COVID(countryname=countryname, beta = beta, R0 = NULL, gamma = 1/14)
	liste_pays$R[i] <- par$R0
}


#### Map R0
# load "World" data from package tmap
data(World) 

# join World data to the list of countries with calculated R0
df_countries <- read_excel("Country.xlsx")[,1:3] %>%
	rename(iso_a3 = `Country Code`, country = `Table Name`) %>%
	left_join(World) %>%
	left_join(liste_pays, by=c("country")) 

# plot the map
st_geometry(df_countries) <- df_countries$geometry
tm_shape(df_countries) +
	tm_polygons("R")



##################################################################################
# solve the model for selected countries
##################################################################################
paramBF <- parameters_SIR_COVID(countryname="Burkina Faso")
paramFR <- parameters_SIR_COVID("France")
paramCH <- parameters_SIR_COVID("China")

##################################################################################
# determine the values of S,I and R at times in vt
##################################################################################
vt = seq(0,350,1)  

# solve models and store reults in a dataframe
mymodel_results_BF <- lsoda(paramBF$inits, vt, calculate_derivatives, paramBF$vparameters) %>% as.data.frame()
mymodel_results_FR <- lsoda(paramFR$inits, vt, calculate_derivatives, paramFR$vparameters) %>% as.data.frame()
mymodel_results_CH <- lsoda(paramCH$inits, vt, calculate_derivatives, paramCH$vparameters) %>% as.data.frame()


##### extract interesting indicators
# calculate total Number of habitants
N_BF <- mymodel_results_BF[1,2:17] %>% sum()
N_FR <- mymodel_results_FR[1,2:17] %>% sum()
N_CH <- mymodel_results_CH[1,2:17] %>% sum()

# calculate total Number of Infected per day
I_BF <- mymodel_results_BF[,18:33] %>% rowSums()
I_FR <- mymodel_results_FR[,18:33] %>% rowSums()
I_CH <- mymodel_results_CH[,18:33] %>% rowSums()

# calculate total Number of Susceptible per day
S_BF <- mymodel_results_BF[,2:17] %>% rowSums()
S_FR <- mymodel_results_FR[,2:17] %>% rowSums()
S_CH <- mymodel_results_CH[,2:17] %>% rowSums()

# calculate total Number of Susceptible per day
R_BF <- mymodel_results_BF[,34:49] %>% rowSums()
R_FR <- mymodel_results_FR[,34:49] %>% rowSums()
R_CH <- mymodel_results_CH[,34:49] %>% rowSums()

# prevalence per day
Res_BF <- I_BF / N_BF
Res_FR <- I_FR / N_FR
Res_CH <- I_CH / N_CH

# plot prevalence of infection (time from the first case)
prev <- data.frame(BF = Res_BF, FR = Res_FR, CH = Res_CH) %>% mutate(days= vt)

prev %>% gather(1:3, key=Country, value=prev) %>%
ggplot(aes(x = days, y = prev,  group = Country, color = Country)) +
	geom_line()

# plot prevalence of immunised people (recovered)
recov <- data.frame(BF = R_BF / N_BF, FR = R_FR / N_FR, CH = R_CH / N_CH) %>% mutate(days= vt)
recov %>% gather(1:3, key=Country, value=recov) %>%
	ggplot(aes(x = days, y = recov,  group = Country, color = Country)) +
	geom_line()


# R0
paramFR$R0
paramBF$R0 
paramCH$R0 


####################################################
# Estimate R0 from case data in all countries (using the method described in report 1 of ETE team)
####################################################

library(R0)
library(lubridate)


##### calculate parameters of the GT distribution (needed to evaluate R0)
interval <- read.csv("interval.txt")
int <- mdy(interval$InfecteeOnset) - mdy(interval$InfectorOnset) 
GT <- est.GT(serial.interval = as.integer(int))

### function that calculate R0 from a vector of daily new cases (from the first to the last case)
R0x <- function(x){
	x <- x[cumsum(x) & rev(cumsum(rev(x)))] # remove leading and ending zeros
	tryCatch(R0 <- estimate.R(epid = x, GT = GT, begin = as.integer("1"), end = as.integer(length(x)), methods=c("ML")), error = function(e) NA)
	tryCatch(return(R0$estimates[[1]]$R), error = function(e) NA)

}

### load dataframe of cases
case <- read.csv("cases_covid_27032020.csv")
case$date <- dmy(case$date)

### test: calculate R0 in France before lockdown (as in report 1 of ETE team)
case %>% filter(location=="France" & date < dmy("16/03/2020")) %>%
	dplyr::select(new_cases) %>%
	R0x()

## R0 in France before 27/03/2020
case %>% filter(location=="France") %>%
	dplyr::select(new_cases) %>%
	R0x()

## R0 in France before 27/03/2020
case %>% filter(location=="Burkina Faso") %>%
	dplyr::select(new_cases) %>%
	R0x()

## R0 for all countries
R0_world <- case %>% 	
	group_by(location) %>%
	summarise(R0 = R0x(new_cases)) 


### plot calculated R0:
# join World data to the list of countries with calculated R0
df_countries_2 <- read_excel("Country.xlsx")[,1:3] %>%
	rename(iso_a3 = `Country Code`, country = `Table Name`) %>%
	left_join(World) %>%
	left_join(R0_world , by=c("country"="location")) 

# plot the map
st_geometry(df_countries_2) <- df_countries_2$geometry
tm_shape(df_countries_2) +
	tm_polygons("R0", breaks = c(0,1,2,3,5,7,10)) 


# compare R0 from case data to R0 from the age-structured model 
R0_comp <- liste_pays %>% left_join(R0_world, by=c("country"="location")) 
R0_comp %>%
	lm(R~R0, data=.)%>%
	summary()
plot(R0_comp$R~R0_comp$R0)
abline(lm(R~R0, data=R0_comp ))
