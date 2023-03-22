#!/usr/bin/env Rscript

library(methylKit)
source("helperFunctions.R")
 
opt_parser = OptionParser(option_list=getOptionsList());
opt = parse_args(opt_parser);

checkInputOptions(opt)

sampleFiles <- getSampleFiles(opt$input)
sampleNames <- getSampleNames(opt$input)

myObj <- getObject(opt, sampleFiles, sampleNames)

if (opt$plots or opt$all) {
    makeMethPlots(myObj) #this function must be made
    makeCovPlots(myObj)
}

if (opt$tables or opt$all) {
    makeMethStats(myObj)
    makeCovStats(myObj)
}

#merge samples

#myObj = filterByCoverage(myObj, lo.count=5, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)

meth = unite(myObj, destrand=FALSE,c.cores=12)

if (opt$tables or opt$all) {
    getCorrelation(meth, plot=FALSE)
    clusterSamples(meth, dist="correlation", method="ward", plot=FALSE)
}

if (opt$plots or opt$all) {
    clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
    PCASamples(meth, screeplot=TRUE)
    PCASamples(meth)
}
