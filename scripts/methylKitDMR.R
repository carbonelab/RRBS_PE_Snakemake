#!/usr/bin/env Rscript

source("scripts/helperFunctions.R")

#Validate script inputs
args = commandArgs(trailingOnly=TRUE)
samples <- checkSamples()
config <- checkConfiguration()

outputDirectory <- args[2]
inputDirectory <- args[1]
sampleFiles <- getSampleFiles(samples[,c(1)], inputDirectory)
treatment <- getTreatmentVector(samples[,c(2)])
sampleNames <- getList(samples[,c(1)])


myObj <- getObject(treatment, sampleFiles, sampleNames, config)
meth = getMergedRegions(myObj, config)

#Merge/Filter

#DMR 

#Bed files