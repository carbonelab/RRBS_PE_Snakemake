# Includes custom functions that do not require methylKit
# Operates on general objects that may be shared across different frameworks

#------------------------------------------------------------------------------------#

# Functions:

# getMergedStats        Outputs merged data statistics
# getMergedPlots        Outputs merged data plots
# getMergedRegions      Returns object of merged dataset
# getMethPlots          Outputs methylation plots for all samples
# getCovPlots           Outputs coverage plots for all samples
# getMethStats          Outputs methylation stats for all samples
# getCovStats           Outputs coverage statistics for all samples
# getObject             Returns data object 
# getList               Returns list object
# getTreatmentVector    Returns treatment vector from sample info file
# getSampleFiles        Returns file names with full directory location 
# checkInputs           Validates input parameters
# checkConfifuration    Verifies project configuration file exists
# checkSamples          Verifies sample info file exists

#------------------------------------------------------------------------------------#

# Load Libraries
library(methylKit)

#------------------------------------------------------------------------------------#

getPlots <- function(myObj, outputDirectory, meth) {
    getMethPlots(myObj, outputDirectory) 
    getCovPlots(myObj, outputDirectory)
    getMergedPlots(meth, outputDirectory)
}

getStats <- function(myObj, outputDirectory, meth) {
    getMethStats(myObj, outputDirectory)
    getCovStats(myObj, outputDirectory)
    getMergedTables(meth, outputDirectory)
}

############------------------------------------#
# FUNCTION # outputs merged data statistics #
############------------------------------------#

getMergedTables <- function(meth, outputDirectory) {
  if(!file.exists(outputDirectory)) {
    dir.create(file.path(outputDirectory, "merged_stats"), showWarnings = FALSE)
  }
  setwd(file.path(outputDirectory, "merged_stats"))

  capture.output((getCorrelation(meth, plot=FALSE)),file="correlation.txt")

  #Doesn't output anything 
  #capture.output((clusterSamples(meth, dist="correlation", method="ward", plot=FALSE)),file="clustering.txt")

  setwd("../../../")
}

############------------------------------------#
# FUNCTION # outputs merged data plots #
############------------------------------------#

getMergedPlots <- function(meth, outputDirectory) {
  
  if(!file.exists(outputDirectory)) {
    dir.create(file.path(outputDirectory, "merged_stats_plots"), showWarnings = FALSE)
  }
  setwd(file.path(outputDirectory, "merged_stats_plots"))
  
  pdf("clusteringDendro.pdf")
  clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
  dev.off()

  pdf("pcaScree.pdf")
  PCASamples(meth, screeplot=TRUE)
  dev.off()

  pdf("pcaScatter.pdf")
  PCASamples(meth)
  dev.off()

  setwd("../../../")
}

############------------------------------------#
# FUNCTION # Returns object of merged dataset #
############------------------------------------#

getMergedRegions <- function(myObj, executionConfiguration) {

  myObj.filtered <- filterByCoverage(myObj, lo.count=10, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)

  myObj.filtered <- normalizeCoverage(myObj.filtered, method="median")
  
  my.tiles <- tileMethylCounts(myObj.filtered,
                                win.size=1000,
                                step.size=1000,
                                cov.bases=10,
                                mc.cores=executionConfiguration$processors
  )

  meth <- methylKit::unite(my.tiles, min.per.group=NULL, mc.cores=executionConfiguration$processors)

  return(meth)
}

############------------------------------------#
# FUNCTION # Outputs methylation plots for all samples #
############------------------------------------#

getMethPlots <- function(myObj, outputDirectory) {
  
  if(!file.exists(outputDirectory)) {
    dir.create(file.path(outputDirectory, "meth_stats_plots"), showWarnings = FALSE)
  }

  setwd(file.path(outputDirectory, "meth_stats_plots"))
 
  for (i in myObj) {
      pdf(paste0(getSampleID(i),"_methstats.pdf"))
      getMethylationStats(i, plot=TRUE, both.strands=FALSE)
      dev.off()
  }
  setwd("../../../")
}

############------------------------------------#
# FUNCTION # Outputs coverage plots for all samples #
############------------------------------------#

getCovPlots <- function(myObj, outputDirectory) {

  if(!file.exists(outputDirectory)) {
    dir.create(file.path(outputDirectory, "cov_stats_plots"), showWarnings = FALSE)
  }

  setwd(file.path(outputDirectory, "cov_stats_plots"))
  
  for (i in myObj) {
      pdf(paste0(getSampleID(i),"_covstats.pdf"))
      getCoverageStats(i, plot=TRUE, both.strands=FALSE)
      dev.off()
  }
  setwd("../../../")
}

############------------------------------------#
# FUNCTION # Outputs methylation stats for all samples #
############------------------------------------#

getMethStats <- function(myObj, outputDirectory) {

  if(!file.exists(outputDirectory)) {
    dir.create(file.path(outputDirectory, "meth_stats"), showWarnings = FALSE)
  }

  setwd(file.path(outputDirectory, "meth_stats"))

  for (i in myObj) {
      capture.output((getMethylationStats(i, plot=FALSE, both.strands=FALSE)),file=paste0(getSampleID(i),"_methstats.txt"))
  }
  setwd("../../../")
}

############------------------------------------#
# FUNCTION # Outputs coverage statistics for all samples #
############------------------------------------#

getCovStats <- function(myObj, outputDirectory) {

  if(!file.exists(outputDirectory)) {
    dir.create(file.path(outputDirectory, "cov_stats"), showWarnings = FALSE)
  }

  setwd(file.path(outputDirectory, "cov_stats"))

  for (i in myObj) {
      capture.output((getCoverageStats(i, plot=FALSE, both.strands=FALSE)),file=paste0(getSampleID(i),"_covstats.txt"))
  }
  setwd("../../../")
}

############------------------------------------#
# FUNCTION # Returns data object #
############------------------------------------#

getObject <- function(treatment, sampleFiles, sampleNames, configuration) {

  myObj <- methRead(sampleFiles,
                    sample.id=sampleNames,
                    header=configuration$header_present,
                    assembly="genome",
                    treatment=treatment,
                    pipeline="bismarkCoverage",
                    context="CpG",
                    mincov=configuration$minimum_coverage)
  return(myObj)
}

############------------------------------------#
# FUNCTION # returns list object #
############------------------------------------#

getList <- function(temp) {
  return(lapply(temp, function(x) x))
}

############------------------------------------#
# FUNCTION # returns treatment vector from sample.info file #
############------------------------------------#

getTreatmentVector <- function(samples) {
  treatment <- c()
  for (i in samples) {
    treatment <- append(treatment, i) 
  }
  return(treatment)
}

############------------------------------------#
# FUNCTION # verifies sample files exist and returns list of files #
############------------------------------------#

getSampleFiles <- function(samples, inputDirectory) {
  files = list()
  for(i in samples) {
    temp <- paste0(inputDirectory,i)
    temp <- paste0(temp,"_val_1_bismark_bt2_pe.bismark.cov.gz")
    if(!file.exists(temp)) {
        print("Sample file not found")
        quit()
    } else {
      files <- append(files,temp)
    }
  }
  return(files)
}

############------------------------------------#
# FUNCTION # validates input parameters #
############------------------------------------#

checkInputs <- function(args, executionConfiguration) {

  if(!file.exists(args[2])) {
	  print("output directory not found, creating one now")
	  dir.create(executionConfiguration$output)
  }

  temp <- list.files(args[1], "*.cov.gz", full=T)
  if(!length(temp) > 1) {
    print("No input files found")
    quit()
  }
}

############------------------------------------#
# FUNCTION # Verifies project configuration file exists #
############------------------------------------#

checkConfiguration <- function() {
  if(!file.exists("proj_config.yaml")) {
    print("Project configuration file not found")
    quit()
  } else {
      return(yaml::read_yaml("proj_config.yaml"))
  }
}

############------------------------------------#
# FUNCTION # Verifies sample info file exists #
############------------------------------------#

checkSamples <- function() {
  if(!file.exists("data/samples.info")) {
    print("Sample information file not found")
    quit()
  } else {
      return(suppressWarnings(read.csv("data/samples.info", header=TRUE)))
  }
}
