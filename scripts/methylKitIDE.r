library(methylKit)
#!/usr/bin/env Rscript
library("optparse")

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

option_list = list(
    make_option(c("-p", "--plots"), action="store_true", default="store_false", dest="plots", help="create plot output files"),
    make_option(c("-t", "--tables"), action="store_true", default="store_false", dest="tables", help="create table output files"),
    make_option(c("-a", "--all"), action="store_true", default="store_false", dest="all", help="create plot and table output files"),
    make_option(c("-i", "--input"), type="character", default=NULL, help="input file directory", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL, help="output file directory", metavar="character"),
    make_option(c("-s", "--treatment"), type="character", default=NULL, help="treatement vector for input files", metavar="character"),
    make_option(c("-c", "--coverage"), type="integer", default=NULL, help="minimum coverage amount", metavar="integer"),
    make_option(c("-h", "--header"), action="store_true", default="store_true", dest="header", help="designate the whether the input files contain a header line")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#check input directory exists
#Enumerate list of input files

#inputDirectory <- paste(getwd(), "/cov_files", sep="")
#if(!file.exists(inputDirectory)) {
if(!file.exists(opt$input)) {
	print("input directory not found")
	quit()
}
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
inputNames <- gsub(paste0(inputDirectory,"/"), "", inputNames)

sampleFiles <- lapply(inputFiles, function(x) x)
sampleNames <- lapply(inputNames, function(x) x)


#check if output directory exists
#create if not

#outputDirectory <- paste(getwd(), "/output/", sep="")
#if(!file.exists(outputDirectory)) {
if(!file.exists(opt$output)) {
	print("output directory not found, creating one now")
	dir.create(opt$output)
}

#create data object
#   include all input parameters as potential command line inputs

myObj = methRead(sampleFiles,
		   	sample.id = sampleNames,
			header = opt$header,
			assembly = "genome",
			treatment = opt$treatment,
			pipeline = "bismarkCoverage",
			mincov = opt$coverage) 

if (plots) {
    #Loop over samples
        #coverage plot
            #create
            #store
        #methylation plot
            #create
            #store
}

if (tables) {
    #Loop over samples
        #coverage table
            #create
            #store
        #methylation table
            #create
            #store
}

#merge samples

meth = unite(myObj, destrand=FALSE,c.cores=12)

if (tables) {
    #correlation table
    #cluster sample table
}


if (plots) {
    #cluster sample ploit
    #PCA scree
    #PCA scatter plot
}
