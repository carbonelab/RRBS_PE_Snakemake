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

if (config$plots || config$all) {
    getMethPlots(myObj, outputDirectory) 
    getCovPlots(myObj, outputDirectory)
    getMergedPlots(meth, outputDirectory)
}
if (config$tables || config$all) {
    getMethStats(myObj, outputDirectory)
    getCovStats(myObj, outputDirectory)
    getMergedTables(meth, outputDirectory)
}