# Includes custom functions that do not require methylKit
# Operates on general objects that may be shared across different frameworks

#------------------------------------------------------------------------------------#

# Functions:

# getOptionsList        Returns list of command line options
# getMergedRegions      Returns object of merged dataset
# getCovStats           Outputs coverage statistics for all samples
# getCovPlots           Outputs coverage plots for all samples
# getMethPlots          Outputs methylation plots for all samples
# getMethStats          Outputs methylation stats for all samples
# getObject             Returns data object 
# getSampleFiles        Returns file names with full directory location 
# getSampleNames        Returns file names without file extentsions or locations
# checkInputOptions     Validates input parameters 
# describe_pca            
# plot_pca



#------------------------------------------------------------------------------------#

# Load Libraries

library(ggplot2)
library("optparse")

#------------------------------------------------------------------------------------#


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

getOptionsList <- function() {
  option_list = list(
    make_option(c("-p", "--plots"), action="store_true", default="store_false", dest="plots", help="create plot output files"),
    make_option(c("-t", "--tables"), action="store_true", default="store_false", dest="tables", help="create table output files"),
    make_option(c("-a", "--all"), action="store_true", default="store_false", dest="all", help="create plot and table output files"),
    make_option(c("-i", "--input"), type="character", default=NULL, help="input file directory", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL, help="output file directory", metavar="character"),
    make_option(c("-s", "--treatment"), type="character", default=NULL, help="treatement vector for input files", metavar="character"),
    make_option(c("-c", "--coverage"), type="integer", default=NULL, help="minimum coverage amount", metavar="integer"),
    make_option(c("-h", "--header"), action="store_true", default="store_true", dest="header", help="designate the whether the input files contain a header line"),
    make_option(c("-r", "--regional"), action="store_true", default="store_false", dest="regional", help="specify whether or not the data should be tiled"),
    make_option(c("-n", "--normalize"), action="store_true", default="store_false", dest="normalize", help="specify whether or not to normalize coverage after filtering"),
    make_option(c("-m", "--processors"), action="store_true", default="store_false", dest="processors", help="specofy the number of processors to utilize in available functions")
  ); 
  return(option_list)
}

getMergedRegions <- function(myObj, opt) {#myObj=myObj, min=10, max=99.9, normalize=TRUE, regional=TRUE, tile_size=1000, step_size=1000, min_cpg=10, mpg=NULL, destrand=FALSE
  
  myObj.filtered <- filterByCoverage(myObj, lo.count=10, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)

  if (opt$normalize) {
    myObj.filtered <- normalizeCoverage(myObj.filtered, method="median")
  }

  if (opt$regional) {
    my.tiles <- tileMethylCounts(myObj.filtered,
                                 win.size=1000,
                                 step.size=1000,
                                 cov.bases=10,
                                 mc.cores=opt$processors
    )
    for (i in 1:length(my.tiles)) {
      print(nrow(my.tiles[[i]]))
    }
    # merge the tiles
    print("Merging tiled regions across samples.")
    meth <- methylKit::unite(my.tiles, min.per.group=NULL, mc.cores=opt$processors)
    print(paste0("Returning ", nrow(meth), " regions."))

    return(meth)
  }

  # if not regional, just merge cpgs with no tiling
  else {
    print("merging CpGs")
    meth <- methylKit::unite(myObj.filtered,
                  min.per.group=NULL,
                  destrand=FALSE)
    print(paste0("Returning ", nrow(meth), " CpGs."))
    return(meth)
  }
}

getMethPlots <- function(myObj, outputDirectory) {

  dir.create(file.path(outputDirectory, "meth_stats_plots"), showWarnings = FALSE)
  setwd(file.path(outputDirectory, "meth_stats_plots"))

  for (i in myObj) {
      pdf(paste0(getSampleID(i),"_methstats.pdf"))
      getMethylationStats(i, plot=TRUE, both.strands=FALSE)
      dev.off()
  }
  setwd(outputDirectory)
}

getCovPlots <- function(myObj, outputDirectory) {

  dir.create(file.path(outputDirectory, "cov_stats_plots"), showWarnings = FALSE)
  setwd(file.path(outputDirectory, "cov_stats_plots"))

  for (i in myObj) {
      pdf(paste0(getSampleID(i),"_covstats.pdf"))
      getMethylationStats(i, plot=TRUE, both.strands=FALSE)
      dev.off()
  }
  setwd(outputDirectory)
}

getMethStats <- function(myObj, outputDirectory) {

  dir.create(file.path(outputDirectory, "meth_stats"), showWarnings = FALSE)
  setwd(file.path(outputDirectory, "meth_stats"))

  for (i in myObj) {
      pdf(paste0(getSampleID(i),"_methstats.txt"))
      getMethylationStats(i, plot=FALSE, both.strands=FALSE)
      dev.off()
  }
  setwd(outputDirectory)
}

getCovStats <- function(myObj, outputDirectory) {

  dir.create(file.path(outputDirectory, "cov_stats"), showWarnings = FALSE)
  setwd(file.path(outputDirectory, "cov_stats"))

  for (i in myObj) {
      pdf(paste0(getSampleID(i),"_covstats.txt"))
      getMethylationStats(i, plot=FALSE, both.strands=FALSE)
      dev.off()
  }
  setwd(outputDirectory)
}

getObject <- function(opt, sampleFiles, sampleNames) {
  myObj <- methRead(sampleFiles,
                    sample.id=sampleNames,
                    assembly="genome",
                    treatment=opt$treatement,
                    pipeline="bismarkCoverage",
                    context="CpG",
                    mincov=opt$coverage)
  return(myObj)
}

getSampleFiles <- function(inputDirectory) {
  inputFiles <- list.files(inputDirectory, "*.cov.gz", full=T)
  inputNames <- gsub("_val.*", "", inputFiles)
  inputNames <- gsub(paste0(inputDirectory,"/"), "", inputNames)
  return(lapply(inputFiles, function(x) x))
}

getSampleNames <- function(inputDirectory) {
  inputFiles <- list.files(inputDirectory, "*.cov.gz", full=T)
  inputNames <- gsub("_val.*", "", inputFiles)
  inputNames <- gsub(paste0(inputDirectory,"/"), "", inputNames)
  return(lapply(inputNames, function(x) x))
}

checkInputOptions <- function(opt) {
  if(!file.exists(opt$input)) {
	  print("input directory not found")
	  quit()
  }

  #if(opt$treatment)

  if(!file.exists(opt$output)) {
	  print("output directory not found, creating one now")
	  dir.create(opt$output)
  }

  if(opt$plots == FALSE AND opt$tables == FALSE AND opt$all == FALSE) {
    print("No data output specified")
    quit()
  }

  temp <- list.files(opt$input, "*.cov.gz", full=T)
  if(!length(temp) > 1) {
    print("No input files found")
    quit()
  }
}

############------------------------------------#
# FUNCTION # to describe a prcomp result object #
############------------------------------------#

describe_pca <- function(mypca) {
  # Determine how many PC's were returned
  print(paste0("Number of PC's returned: ", ncol(mypca$x)))

  # Obtain the eigenvalues (can get values proportional to eigenvalues by taking sd^2)
  eigs <- mypca$sdev^2

  # Determine number of PC's with eigenvalue > 1 (considered important)
  print(paste0("Number of PC's with eigenvalue > 1: ", sum(eigs > 1)))

  # Determine how much of total variance is explained by first PC
  print(paste0("Percent of total variance explained by first PC: ", round((eigs[1]/sum(eigs))*100, digits=2), "%"))

  # How many PC's are needed to explain at least 80% of total variance
  my_sum = 0
  num_pc = 0
  for (i in 1:ncol(mypca$x)) {
    my_sum = my_sum + (eigs[i] / sum(eigs))
    num_pc = num_pc + 1
    if (my_sum >= 0.8) {
      break
    }
  }
  print(paste0("Number of PC's required to explain at least 80% of the variance: ", num_pc))

  return(num_pc)
}

#----------------------------------------------------------------------------------------------------------#

############-------------------------------#
# FUNCTION # to create a biplot from a PCA #
############-------------------------------#
plot_pca <- function(df, pc_a="PC1", pc_b="PC2", color_var, shape_var, label_var, eigs=eigs, num_cpg, tiles=FALSE) {

  # get the percent variance explained by the two PC's
  pc_a_var <- round(((eigs[as.numeric(gsub("PC", "", pc_a))] / sum(eigs)) * 100), 1)
  pc_b_var <- round(((eigs[as.numeric(gsub("PC", "", pc_b))] / sum(eigs)) * 100), 1)

  # subtitle
  #subtitle <- paste0(pc_a, " by ", pc_b)

  # check to see if any variables are "NULL" and change them to NULL
  if (color_var == "NULL") {
    color_var <- NULL
  }
  if (shape_var == "NULL") {
    shape_var <- NULL
  }
  if (label_var == "NULL") {
    label_var <- NULL
  }

  # filename
  if (is.null(shape_var)) {
    filename <- paste0(pc_a, "_", pc_b, "_", "color_by_", color_var, ".pdf")
  } else {
    filename <- paste0(pc_a, "_", pc_b, "_", "color_", color_var, "_shape_", shape_var, ".pdf")
  }

  # General Format of plotting PCA (with snakemake and string variables)
  myplot <- ggplot(df, aes_string(pc_a, pc_b, color=color_var, shape=shape_var)) +
    geom_point(size=5) +
    xlab(paste0(pc_a, ": ", pc_a_var, "% variance")) +
    ylab(paste0(pc_b, ": ", pc_b_var, "% variance")) +
    #if (tiles==TRUE) {
    #  ggtitle(paste0("PCA: CpG Methylation: ", num_cpg, " tiles")) +
    #}
    #if (tiles==FALSE) {
    #  ggtitle(paste0("PCA: CpG Methylation: ", num_cpg, " CpGs")) +
    #}
    geom_text(aes_string(label=label_var),hjust=0.7, vjust=-1.1, show.legend = FALSE) +
    theme_classic() +
    theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) +
    if (tiles==TRUE) {
      ggtitle(paste0("PCA: CpG Methylation: ", num_cpg, " tiles"))
    }
    if (tiles==FALSE) {
      ggtitle(paste0("PCA: CpG Methylation: ", num_cpg, " CpGs"))
    }
  ggsave(filename=filename)

  # General Format of plotting PCA (with manual variables)
  #myplot <- ggplot(df, aes_string(pc_a, pc_b, color=color_var, shape=shape_var)) +
  #  geom_point(size=5) +
  #  xlab(paste0(pc_a, ": ", pc_a_var, "% variance")) +
  #  ylab(paste0(pc_b, ": ", pc_b_var, "% variance")) +
  #  ggtitle("PCA: CpG Methylation") +
  #  geom_text(aes_string(label=label_var),hjust=0.7, vjust=-1.1, show.legend = FALSE) +
  #  theme_classic() +
  #  theme(legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))
  #ggsave(filename=filename)
}

#------------------------------------------------------------------------------------------------------------#

