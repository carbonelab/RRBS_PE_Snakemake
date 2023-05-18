#!/usr/bin/env Rscript

source("scripts/helperFunctions.R")

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

#Output IDE files

getPlots(myObj, outputDirectory, meth)

getStats(myObj, outputDirectory, meth)
