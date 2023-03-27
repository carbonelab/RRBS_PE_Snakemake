#!/usr/bin/env Rscript

library(methylKit)
source("helperFunctions.R")
 
#executionConfiguration_parser = executionConfigurationionParser(executionConfigurationion_list=getexecutionConfigurationionsList());
#executionConfiguration = parse_args(executionConfiguration_parser);

args = commandArgs(trailingOnly=TRUE)
checkInputs(args)

#get sample information
samples <- read.csv("data/samples.info", header=TRUE)

#if all files found, return list in same order as sample.info file
sampleFiles <- getSampleFiles(samples, args[1])

treatment <- getTreatmentVector(samples)

#get configuration settings
executionConfiguration <- read_yaml("proj_config.yaml")

# sampleFiles <- getSampleFiles(args[1])
# sampleNames <- getSampleNames(args[1])

myObj <- getObject(treatement, sampleFiles, samples$Sample,executionConfiguration$minimum_coverage)

if (executionConfiguration$plots or executionConfiguration$all) {
    getMethPlots(myObj, args[2]) #this function must be made
    getCovPlots(myObj, args[2])
}
if (executionConfiguration$tables or executionConfiguration$all) {
    getMethStats(myObj, args[2])
    getCovStats(myObj, args[2])
}

meth = getMergedRegions(myObj, executionConfiguration)

if (executionConfiguration$tables or executionConfiguration$all) {
    getCorrelation(meth, plot=FALSE)
    clusterSamples(meth, dist="correlation", method="ward", plot=FALSE)
}
if (executionConfiguration$plots or executionConfiguration$all) {
    clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
    PCASamples(meth, screeplot=TRUE)
    PCASamples(meth)
}
