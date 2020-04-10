#####
# N MOIROUX : nicolas.moiroux@ird.fr

# load required libraries
require(deSolve)
require(lubridate)
require(tidyverse)
require(readxl)
require(tmap)
require(sf)
require(plotly)
require(multipanelfigure)
source("R/Function_COVID_SIR.R")

#### load data ----
### dataframe with age structure (World bank)
data_pop_all <- read_excel("data/HNP_StatsEXCEL.xlsx") %>% as.data.frame() # data from World Bank

### contact matrices data from https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005697 
# (some sheet names (country) of the oricinal data have been modified to fit with data used for mapping)
data_C_1 <- excel_sheets("data/MUestimates_all_locations_1.xlsx") 
data_C_2 <- excel_sheets("data/MUestimates_all_locations_2.xlsx")

### age classes (from matrix), used to rename row- and col-names of the matrices
age_cl <- c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75+")


#### how much is beta (probability of transmission on contact) for France giving an R0 of 2,5 (as calculated by ETE team) and a gamma (recovery rate) of 11 days-1 ----
#
gamma <- 1/7
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


#### Map R0
# load "World" data from package tmap
data(World) 

# join World data to the list of countries with calculated R0
df_countries <- read_excel("data/Country.xlsx")[,1:3] %>%
	rename(iso_a3 = `Country Code`, country = `Table Name`) %>%
	left_join(World) %>%
	left_join(liste_pays, by=c("country")) 

# plot the map
tmap_mode("view")
st_geometry(df_countries) <- df_countries$geometry
map1 <- tm_shape(df_countries) +
	tm_polygons("R")

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

#### Figure of contact matrices and age strucure for selected countries ----
# list of matrix plots (1 per selected country)
fig_mat <- map(param_list, function(x){
	C <- x$vparameters$C %>% as.table() %>% as.data.frame()
	g <- ggplot(C, aes(Var2, Var1, fill= Freq)) + 
		geom_tile() +
		scale_fill_gradient(low="white",high="darkblue", limits=c(0,25),guide = FALSE, trans="sqrt") +
		labs(y="Age of contact", x="Age of individual") + 
		theme(axis.text.x = element_text(angle = 90, vjust=0.5, color=rep(c(1,0),times=8)))+
		theme(axis.text.y = element_text(color=rep(c(0,1),times=8)))
	
	
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
	S <- x$inits[1:16]
	sC <- C %>% group_by(Var1) %>% 
		summarise(sum = sum(Freq)) %>%
		mutate(age=Var1, pop = S) %>%
		mutate(tot_ctc = sum * pop)
	mean <- sum(sC$tot_ctc) / sum(sC$pop)# calculate mean no. of contact per individual
	sC$age <- factor(sC$age, levels = sC$age) # lock levels order
	bp <- ggplot(sC, aes(x = age, y=sum))+
		geom_bar(stat = "identity") +
		geom_hline(yintercept=mean, linetype="dashed", color = "red") +
		ylim(c(0,45)) +
		ylab("daily no. of contacts") + xlab("") +
		theme(axis.text.x = element_blank())
	return(bp)
})

# list of list of figures
list_l_fig <- list(fig_age, fig_ctc, fig_mat)

# plot a multipanel figure

figure_pop <- multi_panel_figure(columns = length(countries), rows = length(list_l_fig))   # create multipanel figure

for (l in 1:length(list_l_fig)){
	for (i in 1:length(countries)){
		figure_pop %<>%	fill_panel(list_l_fig[[l]][[i]], col=i, row=l)
	}	
}

# mean number of daily contacts in each conties
av_d_ctc <- map(param_list, function(x){
	C <- x$vparameters$C %>% as.table() %>% as.data.frame()
	S <- x$inits[1:16]
	sC <- C %>% group_by(Var1) %>% 
		summarise(sum = sum(Freq)) %>%
		mutate(age=Var1, pop = S) %>%
		mutate(tot_ctc = sum * pop)
	mean <- sum(sC$tot_ctc) / sum(sC$pop)
	return(mean)})

# extract prevalence of infection
prev.df <- map(mymodel_results, function(x){rowSums(x[,18:33]/sum(x[1,2:17]))}) %>% unlist() %>% as.data.frame() %>% mutate(days= rep(vt, length(countries)), Country = rep(countries, each=length(vt)))

colnames(prev.df)[1]<- "prev"

# plot prevalence of infection (time from the first case)
p <- ggplot(prev.df,aes(x = days, y = prev,  group = Country, linetype = Country, color = Country)) +
	geom_line() +
	ylab("Prevalence of infection (% of individuals in stage I)")+
	xlab("Time (days since first case)")

# plot prevalence of immunised people (recovered)
recov.df <- map(mymodel_results, function(x){rowSums(x[,34:49]/sum(x[1,2:17]))}) %>% unlist() %>% as.data.frame() %>% mutate(days= rep(vt, length(countries)), Country = rep(countries, each=length(vt)))
colnames(recov.df)[1]<- "recov"

p2 <-	ggplot(recov.df,aes(x = days, y = recov,  group = Country, color = Country, linetype = Country)) +
	geom_line()+
	ylab("Prevalence of immunization (% of individuals in stage R)")+
	xlab("Time (days since first case)")

#### Severe cases (needing hospitalization) prediction ----
# data needed
inf_sev <- c(0, 0.00408, 0.0104, 0.0343, 0.0425, 0.0816, 0.118, 0.166, 0.184) #Proportion of infected individuals hospitalised (Verity et al. 2020)
age_cl_fat <- c("0_9", "10_19", "20_29", "30_39", "40_49", "50_59", "60_69", "70_79", "80+") # age classes
time_fat <- 600 # time in days at which cumulated nb of severe case are predicted (need to be in vt)

# plot fatality by age classes
list_figure_sev <- map2(mymodel_results,countries, function(x,y){
	v_recov <- x[time_fat,34:49] %>% t() %>% as.vector()
	v_pop1 <- x[1,2:17] %>% t() %>% as.vector()
	
	data_pop <- data_pop_all %>%	filter(`Country Name`==y)
	v_pop2 <-	data_pop  %>%
		filter(str_detect(`Indicator Code`,"SP\\.POP\\.[[:digit:]]{2}[[:alnum:]]{2}\\.(FE|MA)$")) %>%
		separate(`Indicator Code`, c("ind1", "ind2","Age", "Sexe")) %>%
		group_by(Age) %>%
		summarise(pop = sum(`2018`)) %>%
		ungroup() %>% 
		select(pop) %>%
		t() %>%
		as.vector()
	
	# calculate recovery rate for age classes 75-79 and 80+ (we assume same infection ratios)
	pr <- v_recov[16]/v_pop1[16]
	v_recov[16] <- v_pop2[16]*pr
	v_recov[17] <- v_pop2[17]*pr
	
	
	
	col_age_cl_fat <- rep(age_cl_fat, each=2)
	col_age_cl_fat <- col_age_cl_fat[-length(col_age_cl_fat)]
	
	df_fat <- data.frame(recov = v_recov, pop = v_pop2, age = col_age_cl_fat)
	
	P <- df_fat %>% 
		group_by(age) %>%
		summarise(sum_r = sum(recov), sum_p = sum(pop)) %>%
		mutate(r_sev = inf_sev) %>%
		mutate(sev = r_sev*sum_r) %>%
		mutate(p_sev = sev/sum(sev))
	
	p_sev_t <- sum(P$sev) / sum(P$sum_p)
	P$age <- factor(P$age, levels = P$age) # lock levels order
	bp <- ggplot(P, aes(x = age, y=p_sev))+
		geom_bar(stat = "identity") +
		labs(y="frequency", x="Age") + 
		theme(axis.text.x = element_text(angle = 90, vjust=0.5, color=rep(c(1,0),times=8)))
	return(list(bp, p_sev_t ))
	
})


# plot a multipanel figure

figure_sev <- multi_panel_figure(columns = length(countries), rows = 1)   # create multipanel figure

for (i in 1:length(countries)){
	figure_sev %<>%	fill_panel(list_figure_sev[[i]][[1]], col=i)
}	
