# Includes custom functions that do not require methylKit
# Operates on general objects that may be shared across different frameworks

#------------------------------------------------------------------------------------#

# Functions:

# getMergedRegions      Returns object of merged dataset
# getCovPlots           Outputs coverage plots for all samples
# getMethPlots          Outputs methylation plots for all samples
# getCovStats           Outputs coverage statistics for all samples
# getMethStats          Outputs methylation stats for all samples
# getObject             Returns data object 
# getSampleFiles        Returns file names with full directory location 
# getSampleNames        Returns file names without file extentsions or locations
# checkInputexecutionConfigurationions     Validates input parameters
# describe_pca            
# plot_pca



#------------------------------------------------------------------------------------#

# Load Libraries

library(ggplot2)
library(methylKit)

#------------------------------------------------------------------------------------#

############------------------------------------#
# FUNCTION # Returns object of merged dataset #
############------------------------------------#

getMergedRegions <- function(myObj, executionConfiguration) {

  myObj.filtered <- filterByCoverage(myObj, lo.count=10, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)

  #if (executionConfiguration$normalize) {
    myObj.filtered <- normalizeCoverage(myObj.filtered, method="median")
  #}

  #if (executionConfiguration$regional) {
    my.tiles <- tileMethylCounts(myObj.filtered,
                                 win.size=1000,
                                 step.size=1000,
                                 cov.bases=10,
                                 mc.cores=executionConfiguration$processors
    )
    for (i in 1:length(my.tiles)) {
      print(nrow(my.tiles[[i]]))
    }
    # merge the tiles
    print("Merging tiled regions across samples.")
    meth <- methylKit::unite(my.tiles, min.per.group=NULL, mc.cores=executionConfiguration$processors)
    print(paste0("Returning ", nrow(meth), " regions."))

    return(meth)
  #}

  # if not regional, just merge cpgs with no tiling
  # else {
  #   print("merging CpGs")
  #   meth <- methylKit::unite(myObj.filtered,
  #                 min.per.group=NULL,
  #                 destrand=FALSE)
  #   print(paste0("Returning ", nrow(meth), " CpGs."))
  #   return(meth)
  # }
}

############------------------------------------#
# FUNCTION # Outputs methylation plots for all samples #
############------------------------------------#

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

############------------------------------------#
# FUNCTION # Outputs coverage plots for all samples #
############------------------------------------#

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

############------------------------------------#
# FUNCTION # Outputs methylation stats for all samples #
############------------------------------------#

getMethStats <- function(myObj, outputDirectory) {

  dir.create(file.path(outputDirectory, "meth_stats"), showWarnings = FALSE)
  setwd(file.path(outputDirectory, "meth_stats"))

  for (i in myObj) {
      pdf(paste0(getSampleID(i),"_methstats.txt"))
      getMethylationStats(i, plot=FALSE, both.strands=FALSE)
      dev.off()
  }

  print("outputDirectory:")
  print(outputDirectory)
  
  setwd(outputDirectory)
}

############------------------------------------#
# FUNCTION # Outputs coverage statistics for all samples #
############------------------------------------#

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

############------------------------------------#
# FUNCTION # Returns data object #
############------------------------------------#

getObject <- function(treatment, sampleFiles, sampleNames, minimumCoverage) {

  myObj <- methRead(sampleFiles,
                    sample.id=sampleNames,
                    assembly="genome",
                    treatment=treatment,
                    pipeline="bismarkCoverage",
                    context="CpG",
                    mincov=minimumCoverage)
  return(myObj)
}

# ############------------------------------------#
# # FUNCTION # Returns file names with full directory location #
# ############------------------------------------#

# getSampleFiles <- function(inputDirectory) {
#   inputFiles <- list.files(inputDirectory, "*.cov.gz", full=T)
#   return(lapply(inputFiles, function(x) x))
# }

# ############------------------------------------#
# # FUNCTION # Returns file names without file extentsions or locations #
# ############------------------------------------#

# getSampleNames <- function(inputDirectory) {
#   inputFiles <- list.files(inputDirectory, "*.cov.gz", full=T)
#   inputNames <- gsub("_val.*", "", inputFiles)
#   inputNames <- gsub(paste0(inputDirectory,"/"), "", inputNames)
#   return(lapply(inputNames, function(x) x))
# }

getTreatmentVector <- function(samples) {
  treatment <- c()
  for (i in samples) {
    treatment <- append(treatment, i) 
  }
  return(treatment)
}

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

checkInputs <- function() {
  # if(!file.exists("")) {
	#   print("input directory not found")
	#   quit()
  # }

  #if(executionConfiguration$treatment)

  if(!file.exists(executionConfiguration$output)) {
	  print("output directory not found, creating one now")
	  dir.create(executionConfiguration$output)
  }

  if(executionConfiguration$plots == FALSE && executionConfiguration$tables == FALSE && executionConfiguration$all == FALSE) {
    print("No data output specified")
    quit()
  }

  temp <- list.files(executionConfiguration$input, "*.cov.gz", full=T)
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

