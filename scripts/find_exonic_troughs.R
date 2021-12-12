#! /usr/bin/env Rscript
library("changepoint")


############################################################################################
# Setup variables
############################################################################################

args <- commandArgs()
input_filename<-args[6]
output_filename<-args[7]
cpu<-as.numeric(args[8])

all_data<-read.table(input_filename,sep=",",header=FALSE)
colnames(all_data)<-c("type","condition","chromosome","start","end","exon_loc","coveragewithq")
all_data$type<-as.character(all_data$type)
all_data$condition<-as.character(all_data$condition)
all_data$exon_loc<-as.character(all_data$exon_loc)
all_data$chromosome<-as.character(all_data$chromosome)
all_data$coveragewithq<-as.character(all_data$coveragewithq)
all_data$coveragewithq<-lapply(strsplit(all_data$coveragewithq,";"),as.numeric)

############################################################################################
# For testing purposes
############################################################################################
#all_data<-head(all_data,100)

computeCPD <- function(datawithq)
{
	q<-datawithq[[1]]
	data<-datawithq[2:length(datawithq)]
	#print(data)
	#print(q)
	result=cpt.meanvar(data,test.stat="Exponential",method="BinSeg",Q=q)
	cpts.ts(result)
}

############################################################################################
# Unparallelized
############################################################################################
all_results<-lapply(all_data$coveragewithq,computeCPD)

############################################################################################
# Parallelized
############################################################################################
#all_results<-mclapply(all_data$coveragewithq,computeCPD,mc.cores = cpu-1)

results_with_coverage<-data.frame(type=character(),
                 condition=character(),
                 chromosome=character(),
                 start=integer(),
                 end=integer(),
                 k_bkps=integer(),
                 exon_loc=character(),
                 changepoints=character(),
                 #coverage=character(),
                 stringsAsFactors=FALSE) 

for(i in 1:length(all_results))
{
	results_with_coverage[i,]$type <- all_data[i,]$type
	results_with_coverage[i,]$condition<-all_data[i,]$condition
	results_with_coverage[i,]$chromosome <- all_data[i,]$chromosome
	results_with_coverage[i,]$start <- as.numeric(all_data[i,]$start)
	results_with_coverage[i,]$end <- as.numeric(all_data[i,]$end)
	results_with_coverage[i,]$k_bkps <- all_data[i,]$coveragewithq[[1]][1]
	results_with_coverage[i,]$exon_loc <- all_data[i,]$exon_loc
	results_with_coverage[i,]$changepoints <- paste(all_results[[i]],collapse=";")
	
	#results_with_coverage[i,]$coverage <- paste(all_data$coveragewithq[2:length(all_data$coveragewithq)],collapse=";")
}	

write.table(results_with_coverage,file=output_filename,sep=",",row.names=FALSE,quote = FALSE,col.names=FALSE)
