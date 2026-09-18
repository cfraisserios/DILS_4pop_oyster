#!/usr/bin/R
library(tidyverse)
library(FactoMineR)
library(viridis)
library(ggpubr)

## get data 
#datapath='/home/croux/Documents/zoe/DILS/4pop_v1'
#nIteration=10

for(arg in commandArgs()){
	tmp=strsplit(arg, '=')
	if(tmp[[1]][1]=='nIteration'){
		nIteration = as.numeric(tmp[[1]][2])
	}
	if(tmp[[1]][1]=='datapath'){
		datapath = tmp[[1]][2]
	}
}

observation = read.table(paste(datapath, '/ABCstat.txt', sep=''), h=T)

prior = posterior = opti5 = NULL
for(iteration_tmp in 0:(nIteration-1)){
	statistics_tmp = read.table(paste(datapath, '/param/best_', iteration_tmp, '/ABCstat.txt', sep=''), sep='\t', h=T)
	prior = rbind(prior, statistics_tmp)
	
	statistics_tmp = read.table(paste(datapath, '/gof_param/best_', iteration_tmp, '/ABCstat.txt', sep=''), sep='\t', h=T)
	posterior = rbind(posterior, statistics_tmp)
	
	statistics_tmp = read.table(paste(datapath, '/gof_opti5/best_', iteration_tmp, '/ABCstat.txt', sep=''), sep='\t', h=T)
	opti5 = rbind(opti5, statistics_tmp)
}

## parameters
dataset = rbind(observation, prior, posterior, opti5)
toRemove = c(1)
for(i in 1:ncol(dataset)){
	if(sd(dataset[,i])<0.0001){
		toRemove = c(toRemove, i)
	}
}
toRemove = unique(toRemove)
dataset = dataset[, -toRemove]

origin = c('observation', rep('prior', nrow(prior)), rep('posterior', nrow(posterior)), rep('optimization', nrow(opti5)))

dataset = cbind(dataset, origin)

PCA = PCA(dataset[,-ncol(dataset)], graph=F)

dataset = cbind(dataset, PC1=PCA$ind$coord[,1], PC2=PCA$ind$coord[,2], PC3=PCA$ind$coord[,3], PC4=PCA$ind$coord[,4], PC5=PCA$ind$coord[,5])

PCA_tbl = tibble(PC1=PCA$ind$coord[,1], PC2=PCA$ind$coord[,2], PC3=PCA$ind$coord[,3], PC4=PCA$ind$coord[,4], PC5=PCA$ind$coord[,5], origin=dataset$origin)
PCA_tbl = PCA_tbl %>% dplyr::mutate(order=0)
PCA_tbl = PCA_tbl %>% dplyr::mutate(sizePoint=1)

PCA_tbl = PCA_tbl %>% dplyr::mutate(order=ifelse(origin=='prior', 1, order))
PCA_tbl = PCA_tbl %>% dplyr::mutate(order=ifelse(origin=='posterior', 2, order))
PCA_tbl = PCA_tbl %>% dplyr::mutate(order=ifelse(origin=='optimization', 3, order))
PCA_tbl = PCA_tbl %>% dplyr::mutate(order=ifelse(origin=='observation', 4, order))
PCA_tbl = PCA_tbl %>% dplyr::mutate(sizePoint=ifelse(origin=='observation', 2, sizePoint))
PCA_tbl = PCA_tbl %>%  arrange(order) %>%  dplyr::mutate(order=factor(order), sizePoint=factor(sizePoint))

P1 = PCA_tbl %>% ggplot(aes(x=PC1, y=PC2, col=order, size=sizePoint)) + geom_point() +
	scale_color_manual(values=c('grey', viridis(n=3, option='D')), labels=c('prior', 'posterior', 'optimization', 'observation')) +
	theme_bw(base_size=25) +
	labs(color=NULL) +
	guides(size='none', colour = guide_legend(override.aes = list(size=5))) +
	xlab(paste('PC1 (', round(PCA$eig[1,2], 2), '%)', sep='')) +
	ylab(paste('PC2 (', round(PCA$eig[2,2], 2), '%)', sep=''))

P2 = PCA_tbl %>% ggplot(aes(x=PC1, y=PC3, col=order, size=sizePoint)) + geom_point() +
	scale_color_manual(values=c('grey', viridis(n=3, option='D')), labels=c('prior', 'posterior', 'optimization', 'observation')) +
	theme_bw(base_size=25) +
	labs(color=NULL) +
	guides(size='none', colour = guide_legend(override.aes = list(size=5))) +
	xlab(paste('PC1 (', round(PCA$eig[1,2], 2), '%)', sep='')) +
	ylab(paste('PC3 (', round(PCA$eig[3,2], 2), '%)', sep=''))

P3 = PCA_tbl %>% ggplot(aes(x=PC3, y=PC2, col=order, size=sizePoint)) + geom_point() +
	scale_color_manual(values=c('grey', viridis(n=3, option='D')), labels=c('prior', 'posterior', 'optimization', 'observation')) +
	theme_bw(base_size=25) +
	labs(color=NULL) +
	guides(size='none', colour = guide_legend(override.aes = list(size=5))) +
	xlab(paste('PC3 (', round(PCA$eig[3,2], 2), '%)', sep='')) +
	ylab(paste('PC2 (', round(PCA$eig[2,2], 2), '%)', sep=''))

final_plot = ggarrange(P1, P3, P2, labels = c('A', 'B', 'C'), align='hv', common.legend = T, font.label = list(size = 30))

ggsave(filename = paste(datapath, '/gof.pdf', sep=''), plot=final_plot, units='in', dpi=300, width=10, heigh=10)

write.table(x=dataset, file=paste(datapath, '/gof.txt', sep=''), quote=F, sep='\t', col.names=T, row.names=F)

