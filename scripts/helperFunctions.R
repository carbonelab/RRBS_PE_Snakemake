# Includes custom functions that do not require methylKit
# Operates on general objects that may be shared across different frameworks

#------------------------------------------------------------------------------------#

# Load Libraries

library(ggplot2)

#-------------------------------------------------------------------------------------------#

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

