#!/usr/bin/env Rscript

library(methylKit)
source("helperFunctions.R")
 
opt_parser = OptionParser(option_list=getOptionsList());
opt = parse_args(opt_parser);

checkInputOptions(opt)

sampleFiles = getSampleFiles(opt$input)
sampleNames = getSampleNames(opt$input)

myObj = getObject(opt, sampleFiles, sampleNames)

if (opt$plots or opt$all) {
    getMethPlots(myObj, opt$output) #this function must be made
    getCovPlots(myObj, opt$output)
}
if (opt$tables or opt$all) {
    getMethStats(myObj, opt$output)
    getCovStats(myObj, opt$output)
}

meth = getMergedRegions(myObj, opt)

if (opt$tables or opt$all) {
    getCorrelation(meth, plot=FALSE)
    clusterSamples(meth, dist="correlation", method="ward", plot=FALSE)
}
if (opt$plots or opt$all) {
    clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
    PCASamples(meth, screeplot=TRUE)
    PCASamples(meth)
}
