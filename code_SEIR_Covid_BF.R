#####
# N MOIROUX : nicolas.moiroux@ird.fr

# load required libraries
require(deSolve)
require(R0)
require(lubridate)
require(tidyverse)
require(readxl)
require(tmap)
require(sf)
require(plotly)
require(multipanelfigure)
source("Function_COVID_SIR.R")

#### load data ----
### dataframe with age structure (World bank)
data_pop_all <- read_excel("HNP_StatsEXCEL.xlsx") %>% as.data.frame() # data from World Bank

### contact matrices data from https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005697
data_C_1 <- excel_sheets("MUestimates_all_locations_1.xlsx") 
data_C_2 <- excel_sheets("MUestimates_all_locations_2.xlsx")

### age classes (from matrix), used to rename rox- and col-names of the matrices
age_cl <- c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75+")


#### how much is beta (probability of transmission on contact) for France giving an R0 of 2,5 (as calculated by ETE team) and a gamma (recovery rate) of 11 days-1 ----
#
gamma <- 1/11
paramFR <- parameters_SIR_COVID("France", beta = NULL, R0 = 2.5, gamma = gamma, age_cl)
# ignore warnings
beta <- paramFR$vparameters[2] %>% unlist()


#### Calculate R0 for all countries based on the same probability of transmission on contact (beta)----

## dataframe of countries with population data available for 2018
liste_pays <-data_pop_all %>%  
	filter(str_detect(`Indicator Code`,"SP\\.POP\\.[[:digit:]]{2}[[:alnum:]]{2}\\.(FE|MA)$")) %>%
	group_by(`Country Name`) %>%
	summarise(n = sum(!is.na(`2018`))) %>%
	filter(n == 34) %>%
	ungroup() %>%
	rename(country=`Country Name`)

## dataframe of countries with available contact matrix
liste_mat <- c(data_C_1, data_C_2)
liste_mat <- data.frame(country = liste_mat)

## join the two dataframes
liste_pays <- inner_join(liste_mat,liste_pays, by="country") %>%
	mutate(R=NA) 

## calculate R0 for all countries (146/152 having contact matrix data)
for (i in 1:nrow(liste_pays)){
	countryname <- liste_pays[i,1] %>% as.character()
	par <- parameters_SIR_COVID(countryname=countryname, beta = beta, R0 = NULL, gamma = gamma, age_cl)
	liste_pays$R[i] <- par$R0
}


#### Map R0 ----
# load "World" data from package tmap
data(World) 

# join World data to the list of countries with calculated R0
df_countries <- read_excel("Country.xlsx")[,1:3] %>%
	rename(iso_a3 = `Country Code`, country = `Table Name`) %>%
	left_join(World) %>%
	left_join(liste_pays, by=c("country")) 

# plot the map
tmap_mode("view")
st_geometry(df_countries) <- df_countries$geometry
map1 <- tm_shape(df_countries) +
	tm_polygons("R")
map1



##################################################################################
# solve the model for selected countries ----
##################################################################################
countries <- c("Niger", "France", "China", "Germany")
param_list <- map(countries, parameters_SIR_COVID)

##################################################################################
# determine the values of S,I and R at times in vt
##################################################################################
vt = seq(0,600,1)  

# solve models and store results in a list

mymodel_results <- map(param_list,lsoda2, vt)

# extract prevalence of infection
prev.df <- map(mymodel_results, function(x){rowSums(x[,18:33]/sum(x[1,2:17]))}) %>% 
	unlist() %>% 
	as.data.frame() %>% 
	mutate(days= rep(vt, length(countries)), 
				 Country = rep(countries, each=length(vt)))

colnames(prev.df)[1]<- "prev"

# plot prevalence of infection (time from the first case)
p <- ggplot(prev.df,aes(x = days, y = prev,  group = Country, color = Country)) +
	geom_line()

ggplotly(p)

# plot prevalence of immunised people (recovered)
recov.df <- map(mymodel_results, function(x){rowSums(x[,34:49]/sum(x[1,2:17]))}) %>% unlist() %>% as.data.frame() %>% mutate(days= rep(vt, length(countries)), Country = rep(countries, each=length(vt)))
colnames(recov.df)[1]<- "recov"

p2 <-	ggplot(recov.df,aes(x = days, y = recov,  group = Country, color = Country)) +
	geom_line()
ggplotly(p2)


# R0
for (i in 1:length(countries)){
	cat(countries[i], "R0=", param_list[[i]]$R0, "\n")
}

#### Figure of contact matrices and age strucure for selected countries ----
# list of matrix plots (1 per selected country)
fig_mat <- map(param_list, function(x){
	C <- x$vparameters$C %>% as.table() %>% as.data.frame()
	g <- ggplot(C, aes(Var2, Var1, fill= Freq)) + 
		geom_tile() +
		scale_fill_gradient(low="white",high="darkblue", limits=c(0,25),guide = FALSE) +
		labs(y="Age of contact", x="Age of individual") + 
		theme(axis.text.x = element_text(angle = 90, vjust=0.5))
		

	return(g)
})

# list of age distribution plots
fig_age <- map(param_list, function(x){
	S <- x$inits[1:16]
	P <- data.frame(age=age_cl, freq=S / sum(S))
	P$age <- factor(P$age, levels = P$age) # lock levels order
	bp <- ggplot(P, aes(x = age, y=freq))+
		geom_bar(stat = "identity") +
		ylab("frequency") + xlab("")+
		theme(axis.text.x = element_blank())
	return(bp)
})

# daily number of contact per age classes
fig_ctc <- map(param_list, function(x){
	C <- x$vparameters$C %>% as.table() %>% as.data.frame()
	S <- C %>% group_by(Var1) %>% 
		summarise(sum = sum(Freq)) %>%
		mutate(age=Var1)
	S$age <- factor(S$age, levels = S$age) # lock levels order
	bp <- ggplot(S, aes(x = age, y=sum))+
		geom_bar(stat = "identity") +
		ylim(c(0,45)) +
		ylab("daily no. of contacts") + xlab("") +
		theme(axis.text.x = element_blank())
	return(bp)
})

# list of list of figures
list_l_fig <- list(fig_age, fig_ctc, fig_mat)

# expl, multipanel fig

figure_pop <- multi_panel_figure(columns = length(countries), rows = length(list_l_fig))   # create multipanel figure

for (l in 1:length(list_l_fig)){
	for (i in 1:length(countries)){
		figure_pop %<>%	fill_panel(list_l_fig[[l]][[i]], col=i, row=l)
	}	
}

for (i in 1:length(countries)){
	figure_pop %<>%	fill_panel(fig_age[[i]], col=i, row=1)
}	

for (i in 1:length(countries)){
	figure_pop %<>%	fill_panel(fig_mat[[i]], col=i, row=2)
}	
figure_pop


fill_panel(Fig2A, column = 1, row = 1:2) %<>%
	fill_panel(Fig2B, column = 2, row = 1) %<>%
	fill_panel(Fig2C, column = 3, row = 1) %<>%
	fill_panel(Fig2D, column = 2, row = 2)

figure2



####################################################
# Estimate R0 from case data in all countries (using the method described in report 1 of ETE team) ----
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
	tryCatch(R0 <- estimate.R(epid = x, GT = GT, begin = as.integer("1"), end = as.integer(length(x)), methods=c("ML")), error = function(e) NULL)
	if(is.null(return(R0$estimates[[1]]$R))){
		return(NA) } 
		else {
			return(R0$estimates[[1]]$R)
		}
	}


### load dataframe of cases
case <- read.csv("cases_covid_27032020.csv")
case$date <- dmy(case$date)

### test: calculate R0 in France before lockdown (as in report 1 of ETE team)
case %>% filter(location=="France" & date < dmy("16/03/2020")) %>%
	dplyr::select(new_cases) %>%
	R0x()

## R0 in germany 
x <- case %>% filter(location=="Germany") %>%
	dplyr::select(new_cases) 
x <- x[,1]
i<- 2
for (i in 1:70){
	y <- x[i:length(x)]
	y <- y[cumsum(y) & rev(cumsum(rev(y)))] # remove leading and ending zeros
	try(print(estimate.R(y,GT, begin=as.integer("1"), end=as.integer(length(y)),methods=c("ML"))))
}
GT
R0x(x)
## R0 in France before 27/03/2020
case %>% filter(location=="Russia") %>%
	dplyr::select(new_cases) %>% as.vector() %>%
	R0x()

## R0 for all countries
R0_world <- case %>% 	
	group_by(location) %>%
	summarise(R0 = R0x(new_cases)) 

R0_by_date <- function(date1, country){
	tryCatch(x <- case %>% filter(location == country & date < date1 & date > (date1-15)) %>%
		dplyr::select(new_cases) %>% R0x, error = function(e) NA)
	return(x)
}


caseFR <- case %>% filter(location=="France")
R0_time <- map2(caseFR$date, caseFR$location, R0_by_date)
date1 <- ymd("2020-03-27")
date1-15

country <- "France"
R0_by_date(date1, country)

caseFR$R0 <- R0_time

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


#### data french departements
hosp <- read.csv("donnees-hospitalieres-covid19-2020-04-01-19h00.csv", sep=";")

# function that gives the daily number of new hospitalization based on daily number of people hospitalized
firstdiff <- function(x) {
	shifted <- c(NA,x[1:(length(x)-1)])
	return(x-shifted)
}

# calculate R0 since the bigining the first hospitalization
R0_dpt <- hosp %>% filter(sexe ==0 & dep!="") %>%
	group_by(dep) %>%
	mutate(new_cases = firstdiff(hosp)) %>%
	filter(!(is.na(new_cases))) %>%
	summarise(R0 = R0x(new_cases)) 

R02_dpt <- hosp %>% filter(sexe ==0 & dep!="") %>%
	group_by(dep) %>%
	mutate(new_rea = firstdiff(rea)) %>%
	filter(!(is.na(new_rea))) %>%
	summarise(R02 = R0x(new_rea)) 

R0 <- cbind(R0_dpt, R02_dpt)
	
chif_cle <- read.csv("chiffres-cles.csv", sep=",")
	
	
Ro_71 <- hosp %>% filter(sexe ==0 & dep=="71") %>%
	mutate(new_cases = firstdiff(hosp+rea)) 