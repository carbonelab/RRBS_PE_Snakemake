#!/usr/bin/env Rscript

source("scripts/helperFunctions.R")

args = commandArgs(trailingOnly=TRUE)
#checkInputs(args)

#get sample information
samples <- read.csv("data/samples.info", header=TRUE)
executionConfiguration <- yaml::read_yaml("proj_config.yaml")

sampleFiles <- getSampleFiles(samples[,c(1)], args[1])
sampleNames <- lapply(samples[,c(1)], function(x) x)
treatment <- getTreatmentVector(samples[,c(2)])

myObj <- getObject(treatment, sampleFiles, sampleNames, executionConfiguration$minimum_coverage)

if (executionConfiguration$plots || executionConfiguration$all) {
    getMethPlots(myObj, args[2]) 
    getCovPlots(myObj, args[2])
}
if (executionConfiguration$tables || executionConfiguration$all) {
    getMethStats(myObj, args[2])
    getCovStats(myObj, args[2])
}

meth = getMergedRegions(myObj, executionConfiguration)

if (executionConfiguration$tables || executionConfiguration$all) {
    getMergedTables(meth, args[2])
}
if (executionConfiguration$plots || executionConfiguration$all) {
    getMergedPlots(meth, args[2])
}
