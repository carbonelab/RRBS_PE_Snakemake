#!/usr/bin/env Rscript

source("scripts/helperFunctions.R")
source("scripts/methylkit_diff_functions.R")

library(methylKit)
library(ggplot2)
library(umap)
library(tidyr)
library(scales)
library(gridExtra)
library(impute)    

#Validate script inputs

args = commandArgs(trailingOnly=TRUE)
samples <- checkSamples()
config <- checkConfiguration()
checkInputs(args, config)

#Extract data and store human-readable variables

outputDirectory <- args[2]
inputDirectory <- args[1]
sampleFiles <- getSampleFiles(samples[,c(1)], inputDirectory)
treatment <- getTreatmentVector(samples[,c(2)])
sampleNames <- getList(samples[,c(1)])

#Create objects for data analysis 

myObj <- getObject(treatment, sampleFiles, sampleNames, config)
meth = getMergedRegions(myObj, config)

#move to comparison directory
dir.create(args[3])
setwd(args[3])


#Wrapper script generates the information below:

#reorganize() helper function
	#function to create object of sample names
		#need group info
	#function to create treatment vector 
		#need group sizes

#vector of the two group IDs to be compared
groups <- c(1,2)
reorganizedSamples <- NULL

#loop over all samples in the sample file
for(sample in samples) {
	#add sample to vector for reorganization 
	if(sample[2] == groups[0]) {
		reorganizedSamples <- append(reorganizedSamples, sample)
	}
	else if(sample[2] == groups[1]) {
		reorganizedSamples <- append(reorganizedSamples, sample)
	}
}


c1 <- reorganize(meth,
		sample.ids=getList(reorganizedSamples[,c(1)]),
		treatment=getTreatmentVector(reorganizedSamples[,c(2)])
) 

#pass in comparison name 

res1 <- calc.DMRs(c1,
		covariate=NULL,
		overdispersion="MN",
		test="Chisq",
		comparison=args[3],
		meth.diff=10,
		qval=0.1,
		type="DMR",
		mc=8
)

#error handling here to prevent job failure

makeBED(res1,args[3],"DMR")