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

c1 <- reorganize(meth,
		sample.ids=c("UMN_1_1","UMN_3_1","UMN_5_1","UMN_7_1","UMN_11_1","UMN_1_10","UMN_5_10","UMN_3_10","UMN_7_10","UMN_11_12"),
		treatment=c(0,0,0,0,0,1,1,1,1,1)
) 

covariates1 <- data.frame(individual = c("1","3","5","7","11","1","5","3","7","11"))

res1 <- calc.DMRs(c1,
		covariate=covariates1,
		overdispersion="MN",
		test="Chisq",
		comparison="PatientPre_v_PatientPost",
		meth.diff=10,
		qval=0.1,
		type="DMR",
		mc=8
)

#error handling here to prevent job failure

makeBED(res1,"PatientPre_v_PatientPost","DMR")