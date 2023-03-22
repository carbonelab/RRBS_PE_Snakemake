#!/usr/bin/env Rscript

library(methylKit)
source("helperFunctions.R")
 
opt_parser = OptionParser(option_list=getOptionsList());
opt = parse_args(opt_parser);

checkInputOptions(opt)

sampleFiles <- getSampleFiles(opt$input)
sampleNames <- getSampleNames(opt$input)

myObj = methRead(sampleFiles,
		   	sample.id = sampleNames,
			header = opt$header,
			assembly = "genome",
			treatment = opt$treatment,
			pipeline = "bismarkCoverage",
			mincov = opt$coverage) 

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
    #file name
    getCorrelation(meth, plot=FALSE)
    #close file

    #file name
    clusterSamples(meth, dist="correlation", method="ward", plot=FALSE)
    #close file
}


if (opt$plots or opt$all) {
    #file name
    clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
    #close file

    #file name
    PCASamples(meth, screeplot=TRUE)
    #close file

    #file name
    PCASamples(meth)
    #close file
}
