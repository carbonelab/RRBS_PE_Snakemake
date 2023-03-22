library(methylKit)
#!/usr/bin/env Rscript
#library("optparse")

source("helperFunctions.R")

#command line inputs
#
# -p [0,1]          : output plots ?
# -t [0,1]          : output tables ?
# -a [0,1]          : all plots/tables or just merged plots/tables
# -i [directory]    : input directory
# -o [directory]    : output directory
# -s "c(1,1,0,0)"   : treatement??? #is this needed??
# -c [1-100]        : minimum read coverage
# -h [0,1]          : header line present on files?
 
opt_parser = OptionParser(option_list=getOptionsList());
opt = parse_args(opt_parser);

checkInputOptions(opt)

#check input directory exists
#Enumerate list of input files

#inputDirectory <- paste(getwd(), "/cov_files", sep="")
#if(!file.exists(inputDirectory)) {

inputFiles <- list.files(opt$input, "*.cov.gz", full=T)
if(!length(inputFiles) > 1) {
	print("No input files found")
	quit()
}

#check if input files are present
#check if input files are accessible?
#store input file names

#create list of sample names from file names

inputNames <- gsub("_val.*", "", inputFiles)
inputNames <- gsub(paste0(opt$input,"/"), "", inputNames)

sampleFiles <- lapply(inputFiles, function(x) x)
sampleNames <- lapply(inputNames, function(x) x)


#check if output directory exists
#create if not

#outputDirectory <- paste(getwd(), "/output/", sep="")
#if(!file.exists(outputDirectory)) {

#create data object
#   include all input parameters as potential command line inputs

myObj = methRead(sampleFiles,
		   	sample.id = sampleNames,
			header = opt$header,
			assembly = "genome",
			treatment = opt$treatment,
			pipeline = "bismarkCoverage",
			mincov = opt$coverage) 

if (opt$plots or opt$all) {
    for(i in myObj) {
        pdf(paste0(getSampleID(i),"_methStats.pdf"))
        getMethylationStats(i,plot=TRUE,both.strands=FALSE)
        dev.off()

        pdf(paste0(getSampleID(i),"_covStats.pdf"))
        getCoverageStats(i,plot=TRUE,both.strands=FALSE)
        dev.off()

        # file.rename(from = "./methylationStats.pdf", to = paste0(outputDirectory, "methylationStats.pdf"))
        # file.rename(from = "./coverageStats.pdf", to = paste0(outputDirectory, "coverageStats.pdf"))
    }
}

if (opt$tables or opt$all) {
    for(i in myObj) {
        pdf(paste0(getSampleID(i),"_methStats.txt"))
        getMethylationStats(i,plot=FALSE,both.strands=FALSE)
        dev.off()

        pdf(paste0(getSampleID(i),"_covStats.txt"))
        getCoverageStats(i,plot=FALSE,both.strands=FALSE)
        dev.off()

        # file.rename(from = "./coverageStats.pdf", to = paste0(outputDirectory, "coverageStats.pdf"))
        # file.rename(from = "./methylationStats.pdf", to = paste0(outputDirectory, "methylationStats.pdf"))
    }
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
