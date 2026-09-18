#!/usr/bin/R
library(tidyverse)
library(data.table)
library(abcrf)

model_comparison=function(observation, statistics, param, liste_models, ntree, ncores){
	nTot = 10000
	
	tmp = tibble(liste_models)
	groups = tmp %>% group_by_at(param)
	
	groups = groups %>% summarise(nModels=n())
	nGroups = nrow(groups)	
	
	colnames(groups)[1] = 'ID'

	stat_tmp = NULL
	modindex_tmp = NULL
	for(i in groups$ID){
		nModels = groups$nModels[which(groups$ID==i)]
		nSimulation_per_model_to_get = floor(nTot/nModels)
		
		models_tmp = tmp %>% dplyr::filter_at(param, all_vars(.==i)) %>% dplyr::select(models)
		
		for(model_tmp in models_tmp$models){
			sampled = sample(x=1:nrow(statistics[[model_tmp]]), size=nSimulation_per_model_to_get, replace=FALSE)
			stat_tmp = rbind(stat_tmp, statistics[[model_tmp]][sampled,])
			
			modindex_tmp = c(modindex_tmp, rep(i, nSimulation_per_model_to_get))
		}
	}

	# remove unvariant stats
	toRemove = c(1)
	for(i in 1:ncol(stat_tmp)){
		if(sd(stat_tmp[,i])<0.0001){
			toRemove = c(toRemove, i)
		}
	}
	toRemove = unique(toRemove)
	
	# inferences
	modindex_tmp = as.factor(modindex_tmp)
	data1 = data.frame(modindex_tmp, stat_tmp[, -toRemove])
	model.rf1 = abcrf(modindex_tmp~., data = data1, ntree=ntree, paral=TRUE, ncores=ncores)
	print(model.rf1)
	prediction = predict(model.rf1, observation[-toRemove], data1, ntree=ntree, paral=TRUE, ncores=ncores)

	print(groups$ID)
	return(prediction)
}


parameter_estimates=function(obs, sumsta, r, ntree, ncores){
	data2 = data.frame(r, sumsta)
	model.rf.r = regAbcrf(r~., data2, ntree=ntree, ncores=ncores, paral=TRUE)
	predicted_param = predict(model.rf.r, obs, data2)
	return(predicted_param)
}


## get data 
#datapath='/home/croux/Documents/zoe/DILS/4pop_v1'
#nIteration=10
ntree=1000
ncores=5

for(arg in commandArgs()){
	tmp=strsplit(arg, '=')
	if(tmp[[1]][1]=='nIteration'){
		nIteration = as.numeric(tmp[[1]][2])
	}
	if(tmp[[1]][1]=='datapath'){
		datapath = tmp[[1]][2]
	}
	if(tmp[[1]][1]=='subdir'){ # param / opti1 / opti2
		subdir = tmp[[1]][2]
	}
}

observation = read.table(paste(datapath, '/ABCstat.txt', sep=''), h=T)

statistics = NULL
parameters = NULL
for(iteration_tmp in 0:(nIteration-1)){
	statistics_tmp = read.table(paste(datapath, '/', subdir, '/best_', iteration_tmp, '/ABCstat.txt', sep=''), sep='\t', h=T)
	parameters_tmp = read.table(paste(datapath, '/', subdir, '/best_', iteration_tmp, '/priorfile.txt', sep=''), sep='\t', h=T)
	
	statistics = rbind(statistics, statistics_tmp)
	parameters = rbind(parameters, parameters_tmp)
}

## parameters
toRemove = c(1)
for(i in 1:ncol(statistics)){
	if(sd(statistics[,i])<0.0001){
		toRemove = c(toRemove, i)
	}
}
toRemove = unique(toRemove)
 
sumsta=statistics[, -toRemove]
obs=observation[-toRemove]
 
liste_params = c()
expectation = c()
median = c()
variance = c()
variance_cdf = c()
quantiles = c()
NMAEmean = c()
minPrior = c()
medianPrior = c()
maxPrior = c()
 
nParams = ncol(parameters)
for(i in 1:nParams){
	param_Name = colnames(parameters)[i]
	r=parameters[, i]

#	if(param_Name%in%c("N1","N2","N3","N4","Na_34","Na_24","Na","Tsplit_34","Tsplit_24","Tsplit","Tsc_23","Tsc_24")){
#		r=r/parameters[[selection]][,7]
#	}

	predicted_param = parameter_estimates(obs=obs, sumsta=sumsta, r=r, ncores=ncores, ntree=ntree)
	
	liste_params = c(liste_params, param_Name)
	expectation = c(expectation, predicted_param$expectation)
	median = c(median, predicted_param$med)
	variance = c(variance, predicted_param$variance)
	variance_cdf = c(variance_cdf, predicted_param$variance.cdf)
 	quantiles = c(quantiles, predicted_param$quantiles)
	NMAEmean = c(NMAEmean, predicted_param$post.NMAE.mean)
	minPrior = c(minPrior, min(r))
	maxPrior = c(maxPrior, max(r))
	medianPrior = c(medianPrior, median(r))
}
 
res = tibble(parameters=liste_params, expectation=expectation, median=median, variance=variance, variance_cdf=variance_cdf, post_NMAE_mean=NMAEmean, quantile_0_025=matrix(quantiles, ncol=2, byrow=T)[,1], quantile_0_975=matrix(quantiles, ncol=2, byrow=T)[,2], minPrior=minPrior, medianPrior=medianPrior, maxPrior=maxPrior)
 
write.table(x=res, file=paste(datapath, '/', subdir, '/posterior.txt', sep=''), quote=F, sep="\t", col.names=T, row.names=F)

