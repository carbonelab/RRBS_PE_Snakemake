# Using methylKit framework/functions to analyze bisulfite sequencing data
# Uses output of Bismark as input: cytosine reports or coverage files
# Requires a metadata excel file with sample names and groups/covariates

#------------------------------------------------------------------------------#

# Scripts:
# 1.  create.myObj()              # read in coverage files and create methylKit object 
# 2.  plot_cov_dist()             # lineplot or boxplots of cpg coverage distribution per sample (all on one plot)
# 3.  number_cpg_at_cov_barplot() # Barplot of number of CpGs >= 10X cov (using cpg cov table as input)
# 4.  get_promoter_gr()           # use methylKit to create a GRanges object of promoter coordinates from a gtf file 
# 5.  summarize_prom_meth()       # use methylKit to summarize coverage and methylation of myObj over the promoters GRanges object
# 6.  promoter_cov_pie()          # promoter cov pie charts, one plot file per sample AND all samples on one plot file
# 7.  promoter_cov_bar()          # bar height represents percent of gene promoters with X coverage
# 8.  cpg_cov()                   # Produce CpG Coverage Table and filter on coverage
# 9.  makeCovPlots()              # Make a methylKit coverage histogram per sample
# 10. makeMethStats()             # Make a methylKit percent methylation histogram per sample
# 11. get_merged_regions()        # merging and more, cpg or regional

#------------------------------------------------------------------------------#
# Load Libraries

library(methylKit)
library(ggplot2)
library(umap)
library(tidyr)
library(scales)
library(gridExtra) # for arranging multi-sample plots on one file
library(impute)    # for imputing missing values

#------------------------------------------------------------------------------#

############------------------------------------------------#
# FUNCTION # to create the methyl object from the cov files #
############------------------------------------------------#

# cov_file_dir: path to canonical coverage files directory
# metadata: path to metadata txt file
# group_var: column name in metadata file that signifies group membership
# suffix: extension of canonical coverage files to remove, leaving sample name
create.myObj <- function(cov_file_dir,  metadata, group_var="group", suffix="_can_cov.gz") {
  # capture the cov.file paths
  my_cov_files <- list.files(cov_file_dir, full=T)
  cov_list <- lapply(my_cov_files, function(x) x)

  # capture the samp.names
  my_samp_names <- gsub(cov_file_dir, "", my_cov_files)
  my_samp_names <- gsub("/", "", my_samp_names)
  my_samp_names <- gsub(suffix, "", my_samp_names)
  samp_name_list <- lapply(my_samp_names, function(x) x)

  # read in metadata table
  meta <- read.delim(metadata, header=T)

  # create vector of treatments from matching sample with group in metadata file
  my_treat <- sapply(my_samp_names, function(x) meta[meta$sample == x, group_var])
  my_treat_numeric <- sapply(my_treat, function(x) as.numeric(x))

  # use methRead function from methylKit to read in the cov.files
  myObj <- methRead(cov_list,
                    sample.id=samp_name_list,
                    assembly="genome",
                    treatment=my_treat_numeric,
                    pipeline="bismarkCoverage",
                    context="CpG",
                    mincov=1)
  return(myObj)
}

#----------------------------------------------------------------------------------------#

############--------------------------------#
# FUNCTION # to plot coverage distributions #
############--------------------------------#
# Can produce lineplots or boxplots
# Includes all samples on one plot

plot_cov_dist <- function(myObj=myObj, filename="cov.dist.pdf", type="line", fill_col="lightskyblue1", width=7, height=7) {
  # get a list of vectors of coverages from each sample df in myObj
  coverage <- sapply(myObj, function(x) getData(x)$coverage)
  # combine into a single vector
  coverage_vec <- unlist(coverage)

  # get a list of vectors of sample.ids
  sample_list <- list()
  for (i in 1:length(coverage)) {
    sample_list[[i]] <- c(rep(getSampleID(myObj[[i]]), length(coverage[[i]])))
  }
  # combine into single vector
  sample_vec <- unlist(sample_list)

  # create df for plotting
  df <- data.frame(coverage=coverage_vec,
                   sample=sample_vec
  )

  if (type=="line") {
    # plot
    myplot <- ggplot(df, aes(x=coverage, color=sample)) +
              geom_freqpoly() +
              scale_y_log10() +
              theme_minimal()
    ggsave(filename=filename, width=width, height=height)
  }
  else if (type=="box") {
    # plot
    myplot <- ggplot(df, aes(x=sample, y=coverage)) +
              geom_boxplot(fill=fill_col) +
              scale_y_log10() +
              coord_flip() +
              theme_minimal()
    ggsave(filename=filename, width=width, height=height)
  }
}




############--------------------------------#
# FUNCTION # to plot coverage distributions #
############--------------------------------#
# Can produce lineplots or boxplots
# Includes all samples on one plot

#plot_cov_dist <- function(myObj=myObj, filename="cov.dist.pdf", type="line", fill_col="lightskyblue1") {
#  # get a list of vectors of coverages from each sample df in myObj
#  coverage <- sapply(myObj, function(x) getData(x)$coverage)
#  # combine into a single vector
#  coverage_vec <- unlist(coverage)
#
#  # get a list of vectors of sample.ids
#  sample_list <- list()
#  for (i in 1:length(coverage)) {
#    sample_list[[i]] <- c(rep(getSampleID(myObj[[i]]), length(coverage[[i]])))
#  }
#  # combine into single vector
#  sample_vec <- unlist(sample_list)
#
#  # create df for plotting
#  df <- data.frame(coverage=coverage_vec,
#                   sample=sample_vec
#  )
#
#  if (type=="line") {
#    # plot
#    myplot <- ggplot(df, aes(x=coverage, color=sample)) +
#              geom_freqpoly() +
#              scale_y_log10() +
#              theme_minimal()
#    ggsave(filename=filename)
#  }
#  else if (type=="box") {
#    # plot
#    myplot <- ggplot(df, aes(x=sample, y=coverage)) +
#              geom_boxplot(fill=fill_col) +
#              scale_y_log10() +
#              coord_flip() +
#              theme_minimal()
#    ggsave(filename=filename)
#  }
#}

#-----------------------------------------------------------------------------------------------#

############--------------------------------------------------------------#
# FUNCTION # to make barplot of number of CpGs with at least 10X coverage #
############--------------------------------------------------------------#
# use the CpG.coverage.table.txt file created from "cpg_cov" function

number_cpg_at_cov_barplot <- function(cpg_cov_table, filename="cpv_cov_barplot.png") {
  # read in the cpg cov table
  df <- read.delim(cpg_cov_table, header=T)

  myplot <- ggplot(df, aes(x=reorder(sample, -CpG.10.), y=CpG.10.)) +
    geom_bar(stat="identity") +
    ylab("Number of CpGs with >= 10X cov") +
    xlab("Sample") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) +
    scale_y_continuous(breaks=c(500000,1000000,1500000,2000000,2500000))+
    scale_y_continuous(labels = comma)
  ggsave(filename=filename)
}

#-----------------------------------------------------------------------------------------------#

#-------------------#
# Promoter Coverage #
#-------------------#

############------------------------------------------------#
# FUNCTION # to obtain GRanges object of promoters from gtf #
############------------------------------------------------#
# get_promoter_gr: use methylKit's "regionCounts" to summarize methylation and
#                  coverage over the promoters of a user-provided gtf file
#                  from a given methylkit object

# gtf: gtf file to use
# skip: number of lines to skip of gtf (num of lines in header)
# biotype: type of genes to include from gtf (protein_coding, miRNA, lncRNA, etc)
# tss_up: bp upstream from TSS to define promoter
# tss_dn: bp downstream from TSS to define promoter

get_promoter_gr <- function(gtf, skip=5, biotype="protein_coding", tss_up=3000, tss_dn=200) {
  # read in the gtf file
  annotation_gtf <- read.delim(gtf, header=F, skip=skip)
  # subset to genes
  anno_gtf <- annotation_gtf[annotation_gtf$V3 == "gene", ]
  # subset to protein coding genes
  # separate the 9th field on semicolon delimiter
  #   first remove all spaces
  anno_gtf$V9 <- gsub(" ", "", anno_gtf$V9)
  foo <- anno_gtf %>% separate(V9, c("GeneID","Version","gene.name","gene.source","biotype"), ";")
  # remove "gene_id" text
  foo$GeneID <- gsub("gene_id", "", foo$GeneID)
  # remove "gene_biotype" text
  foo$biotype <- gsub("gene_biotype", "", foo$biotype)
  # remove "gene_name" text
  foo$gene.name <- gsub("gene_name", "", foo$gene.name)
  # subset to only protein coding biotype genes
  print(paste0("Subsetting to genes from ", biotype))
  foo2 <- foo[foo$biotype %in% biotype, ]
  print(paste0("Returning ", nrow(foo2), " promoter regions."))
  # convert to GRanges object
  annot_gr <- makeGRangesFromDataFrame(foo2,
                                       seqnames.field="V1",
                                       start.field="V4",
                                       end.field="V5",
                                       strand.field="V7"
  )
  # add gene id column to granges obj
  annot_gr$gene_id <- foo2$GeneID
  # add gene name to granges obj
  annot_gr$gene_name <- foo2$gene.name

  # create new granges obj of the promoter region of each gene in original granges obj
  # promoter region = (TSS - 2000bp):(TSS + 200bp - 1)
  my_promoters <- promoters(annot_gr, upstream=tss_up, downstream=tss_dn)

  return(my_promoters)
}

#------------------------------------------------------------------------------------------------------#

############-------------------------------------------------#
# FUNCTION # to produce cov and meth_rate table of promoters #
############-------------------------------------------------#

# summarize_prom_meth: produce a summary table of coverage and meth rates per promoter
#                      from a user-provided GRanges promoter object and methylkit myObj

# myObj: methylkit object to use for cov/meth info
# prom_gr: GRanges object of promoter coordinates
# export: if TRUE, export the full table of

summarize_prom_meth <- function(myObj=myObj, prom_gr=prom_gr, export=TRUE) {
  # collect methylation and coverage info over each promoter from the filtered gene object
  # uses methylKit's "regionCounts" function
  # store in a methylRawList object that contains methylRawObjects
  prom.sum <- regionCounts(myObj, prom_gr)

  # Convert the "my_promoters" genomic ranges object to a data.frame
  my_promoters_df <- as.data.frame(prom_gr)
  # subset to only include columns: seqnames, start, gene_id, gene_name
  my_promoters_df <- my_promoters_df[ , c("seqnames","start","gene_id", "gene_name")]

  # extract the data from the methylRawList and store in a list of data frames
  prom.sum_dflist <- lapply(prom.sum, function(x) getData(x))

  # get a list of vectors of sample.ids
  sample_list <- list()
  for (i in 1:length(prom.sum_dflist)) {
    sample_list[[i]] <- c(getSampleID(myObj[[i]]))
  }
  # combine into single vector
  sample_vec <- unlist(sample_list)

  # add names to each element in this list of df's, use "sample_names" variable initiated above
  #names(prom.sum_dflist) <- sample_vec

  # create a list of the df's containing meth,coverage info and gene id info
  prom.sum_mergedlist <- lapply(prom.sum_dflist, function(x) merge(x,
                                                                   my_promoters_df,
                                                                   by.x = c("chr","start"),
                                                                   by.y = c("seqnames","start"))
  )

  # set names of each entry in the df list as sample name
  names(prom.sum_mergedlist) <- sample_vec

  # change the class of the coverage and numCs columns to numeric
  for (i in 1:length(prom.sum_mergedlist)) {prom.sum_mergedlist[[i]]$coverage <- as.numeric(prom.sum_mergedlist[[i]]$coverage)}
  for (i in 1:length(prom.sum_mergedlist)) {prom.sum_mergedlist[[i]]$numCs <- as.numeric(prom.sum_mergedlist[[i]]$numCs)}

  # add a meth_rate column to each df
  for (i in 1:length(prom.sum_mergedlist)) {prom.sum_mergedlist[[i]]$meth_rate <- round(100 * (prom.sum_mergedlist[[i]]$numCs / prom.sum_mergedlist[[i]]$coverage), 2)}

  # reorganize the columns in each df
  for (i in 1:length(prom.sum_mergedlist)) {prom.sum_mergedlist[[i]] <- prom.sum_mergedlist[[i]][ ,c(1,2,3,4,8,9,5,6,7,10)]}

  # Export? Probably not necessary, too many files
  #setwd("/u1/bd/Projects/ECP25/ide/promoter_coverage/1x_3200")
  #for (i in 1:length(prom.sum_mergedlist)) {write.table(prom.sum_mergedlist[[i]],
  #                                                      file=paste0(samp_names[i], ".1x3200.prom_stats.txt"),
  #                                                     sep="\t",
  #                                                     quote=FALSE,
  #                                                     row.names=FALSE)
  #}

  # prepare for merging all files
  # add the sample name to the "coverage","numCs", "numTs", and "meth_rate" columns in each df
  for (i in 1:length(prom.sum_mergedlist)) {colnames(prom.sum_mergedlist[[i]])[7] <- paste0("coverage","_",names(prom.sum_dflist)[i])}
  for (i in 1:length(prom.sum_mergedlist)) {colnames(prom.sum_mergedlist[[i]])[8] <- paste0("numCs","_",names(prom.sum_dflist)[i])}
  for (i in 1:length(prom.sum_mergedlist)) {colnames(prom.sum_mergedlist[[i]])[9] <- paste0("numTs","_",names(prom.sum_dflist)[i])}
  for (i in 1:length(prom.sum_mergedlist)) {colnames(prom.sum_mergedlist[[i]])[10] <- paste0("meth_rate","_",names(prom.sum_dflist)[i])}

  return(prom.sum_mergedlist)

  # merge all of these dfs on chr and start
  #final_table <- Reduce(function(df1, df2) merge(df1, df2, by = c("chr","start")), prom.sum_mergedlist)

  #if (export==TRUE) {
  #  # export the final table
  #  write.table(final_table, "promoter_stats_table.txt", sep="\t", col.names=T, row.names=F, quote=F)
  #}

  #return(final_table)
}

#-------------------------------------------------------------------------------------------------------------#

############--------------------------------------------#
# FUNCTION # to produce pie charts of promoter coverage #
############--------------------------------------------#
# Produce individual pie charts of each sample showing percent of promoters
#  with various coverage amount categories
# Option to produce individual plot files per sample or all plots on one file

# promoter_cov_pie: function to plot pie charts of promoter cpg coverage per sample
# total: total number of gene promoters considered
# outdir: desired name of output directory (to be created) to store pie charts

promoter_cov_pie <- function(my_df_list, total, outdir="cpg_cov_pie_charts") {
  # initialize a plotdf list to hold data frames for each sample to be used for plotting
  plotdf_list <- list()
  # Loop through each df in the list and create a new plotdf for each that contains the coverage at each coverage group
  # Also, change the factor levels order for proper plotting
  for (i in 1:length(my_df_list)){
    plotdf_list[[i]] <- data.frame(coverage=c("0","1-4","5-9","10-24","25-49","50-99",">99"),
                                   value=c(total - nrow(my_df_list[[i]]),
                                          nrow(my_df_list[[i]][my_df_list[[i]]$coverage > 0 & my_df_list[[i]]$coverage < 5, ]),
                                          nrow(my_df_list[[i]][my_df_list[[i]]$coverage > 4 & my_df_list[[i]]$coverage < 10, ]),
                                          nrow(my_df_list[[i]][my_df_list[[i]]$coverage > 9 & my_df_list[[i]]$coverage < 25, ]),
                                          nrow(my_df_list[[i]][my_df_list[[i]]$coverage > 24 & my_df_list[[i]]$coverage < 50, ]),
                                          nrow(my_df_list[[i]][my_df_list[[i]]$coverage > 49 & my_df_list[[i]]$coverage < 100, ]),
                                          nrow(my_df_list[[i]][my_df_list[[i]]$coverage > 99, ]))
    )

    plotdf_list[[i]]$coverage <- factor(plotdf_list[[i]]$coverage, levels = c("0","1-4","5-9","10-24","25-49","50-99",">99"))
  }
  names(plotdf_list) <- names(my_df_list)

  # PLOT #
  # create outdir
  dir.create(outdir)
  setwd(outdir)

  # create blank theme
  blank_theme <- theme_minimal()+
  theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.border = element_blank(),
  panel.grid=element_blank(),
  axis.ticks = element_blank(),
  plot.title=element_text(size=14, face="bold"))

  # loop through the plotdf_list and plot each sample
  for (i in 1:length(plotdf_list)){
    myplot <- ggplot(plotdf_list[[i]], aes(x="", y=value, fill=coverage)) +
              geom_bar(width = 1, stat = "identity") +
              coord_polar("y", start=0) +
              scale_fill_brewer(palette="Blues") +
              blank_theme +
              theme(axis.text.x=element_blank()) +
              labs(title=paste0(names(plotdf_list)[i]," Promoter Coverage: ", total, " promoters"))

    ggsave(myplot, filename=paste0(names(plotdf_list)[i],"_promoter_cov.png"))
  }

  #-----------------------------------#
  # put all sample plots on one page

  # order list of df's by ascending number of 0 coverage
  ### store names of df's (sample names)
  samp_names <- names(plotdf_list)
  ### loop through the list of df's
  for (i in 1:length(plotdf_list)) {
    # add sample name to a colname for storage
    colnames(plotdf_list[[i]])[2] <- paste0(names(plotdf_list)[i], "_value")
    # rename each df to be the 0 coverage value
    names(plotdf_list)[i] <- plotdf_list[[i]][1,2]
  }
  # store a vector of names (now, representing 0 coverage value) in ascending order
  foo <- as.character(sort(as.numeric(names(plotdf_list))))
  # reorder the list of df's
  newlist <- plotdf_list[foo]
  # rename each df to be sample name and remove sample name from colnames
  for (i in 1:length(newlist)) {
    names(newlist)[i] <- gsub("_value", "", colnames(newlist[[i]])[2])
    colnames(newlist[[i]])[2] <- "value"
  }

  # Define plotting function of pie chart
  pie_fxn <- function(df, samp_name, legend=TRUE) {
    if(legend==TRUE) {
      myplot <- ggplot(df, aes(x="", y=value, fill=coverage)) +
        geom_bar(width = 1, stat = "identity") +
        coord_polar("y", start=0) +
        scale_fill_brewer(palette="Blues") +
        blank_theme +
        theme(axis.text.x=element_blank()) +
        #labs(title=paste0(names(plotdf_list)[i]," Promoter Coverage: ", total, " promoters"))
        labs(title=samp_name)
    }
    if(legend==FALSE) {
      myplot <- ggplot(df, aes(x="", y=value, fill=coverage)) +
        geom_bar(width = 1, stat = "identity") +
        coord_polar("y", start=0) +
        scale_fill_brewer(palette="Blues") +
        blank_theme +
        theme(axis.text.x=element_blank()) +
        #labs(title=paste0(names(plotdf_list)[i]," Promoter Coverage: ", total, " promoters"))
        labs(title=samp_name) +
        theme(legend.position="none")
    }
    return(myplot)
  }

  # create one plot to get legend
  #p1 <- pie_fxn(plotdf_list[[1]])

  # extract legend
  #g_legend <- function(a.gplot){
  #  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  #  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  #  legend <- tmp$grobs[[leg]]
  #  return(legend)
  #}
  #mylegend <- g_legend(p1)

  plotdf_list <- newlist
  # re-run plot function, this time on list of df's and without legends
  myplots <- lapply(seq_along(plotdf_list), function(i) {
    pie_fxn(df=plotdf_list[[i]], samp_name=names(plotdf_list)[i], legend=FALSE)
  })

  n <- length(myplots)
  nCol <- floor(sqrt(n))
  #do.call("grid.arrange", c(myplots, ncol=nCol))
  grid.arrange(grobs=myplots)
  ggsave(file="promoter_cov_all_piecharts.png", arrangeGrob(grobs=myplots, ncol=nCol), width=16, height=16)
  #g <- arrangeGrob(myplots)
  #ggsave(filename="promoter_cov_all_piecharts.png", g)

  setwd("..")
  return(myplots)
}











#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX#
############--------------------------------------------#
# FUNCTION # to produce pie charts of promoter coverage #
############--------------------------------------------#
# Produce individual pie charts of each sample showing percent of promoters
#  with various coverage amount categories
# Option to produce individual plot files per sample or all plots on one file

# promoter_cov_pie: function to plot pie charts of promoter cpg coverage per sample
# total: total number of gene promoters considered
# outdir: desired name of output directory (to be created) to store pie charts

#promoter_cov_pie <- function(my_df_list, total, outdir="cpg_cov_pie_charts") {
#  # initialize a plotdf list to hold data frames for each sample to be used for plotting
#  plotdf_list <- list()
#  # Loop through each df in the list and create a new plotdf for each that contains the coverage at each coverage group
#  # Also, change the factor levels order for proper plotting
#  for (i in 1:length(my_df_list)){
#    plotdf_list[[i]] <- data.frame(coverage=c("0","1-4","5-9","10-24","25-49","50-99",">99"),
#                                   value=c(total - nrow(my_df_list[[i]]),
#                                          nrow(my_df_list[[i]][my_df_list[[i]]$coverage > 0 & my_df_list[[i]]$coverage < 5, ]),
#                                          nrow(my_df_list[[i]][my_df_list[[i]]$coverage > 4 & my_df_list[[i]]$coverage < 10, ]),
#                                          nrow(my_df_list[[i]][my_df_list[[i]]$coverage > 9 & my_df_list[[i]]$coverage < 25, ]),
#                                          nrow(my_df_list[[i]][my_df_list[[i]]$coverage > 24 & my_df_list[[i]]$coverage < 50, ]),
#                                          nrow(my_df_list[[i]][my_df_list[[i]]$coverage > 49 & my_df_list[[i]]$coverage < 100, ]),
#                                          nrow(my_df_list[[i]][my_df_list[[i]]$coverage > 99, ]))
#    )
#
#    plotdf_list[[i]]$coverage <- factor(plotdf_list[[i]]$coverage, levels = c("0","1-4","5-9","10-24","25-49","50-99",">99"))
#  }
#  names(plotdf_list) <- names(my_df_list)
#
#  # PLOT #
#  # create outdir
#  dir.create(outdir)
#  setwd(outdir)
#  # create blank theme
#  blank_theme <- theme_minimal()+
#  theme(
#  axis.title.x = element_blank(),
#  axis.title.y = element_blank(),
#  panel.border = element_blank(),
#  panel.grid=element_blank(),
#  axis.ticks = element_blank(),
#  plot.title=element_text(size=14, face="bold"))
#
#  # loop through the plotdf_list and plot each sample
#  for (i in 1:length(plotdf_list)){
#    myplot <- ggplot(plotdf_list[[i]], aes(x="", y=value, fill=coverage)) +
#              geom_bar(width = 1, stat = "identity") +
#              coord_polar("y", start=0) +
#              scale_fill_brewer(palette="Blues") +
#              blank_theme +
#              theme(axis.text.x=element_blank()) +
#              labs(title=paste0(names(plotdf_list)[i]," Promoter Coverage: ", total, " promoters"))
#
#    ggsave(myplot, filename=paste0(names(plotdf_list)[i],"_promoter_cov.png"))
#  }
#  setwd("..")
#}

#-----------------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------------#








############------------------------------------------------------------------------------------------------------------#
# FUNCTION # to produce promoter coverage barplot where bar height represents percent of gene promoters with X coverage #
############------------------------------------------------------------------------------------------------------------#
# plot barchart of all samples in one plot,
# bar height represents percent of gene promoters with X coverage

promoter_cov_bar <- function(my_df_list, cov=5, prefix="") {
  # loop through each sample df to get number of promoters with at least X coverage
  num_proms_list <- list()
  for (i in 1:length(my_df_list)) {
    num_proms_list[[i]] <- nrow(my_df_list[[i]][my_df_list[[i]]$coverage >= cov, ])
  }

  # create new df for plotting
  plotdf <- data.frame(sample=names(my_df_list),
                       number=unlist(num_proms_list)
  )

  # PLOT #
  myplot <- ggplot(plotdf, aes(x=sample, y=number)) +
                   geom_bar(stat="identity", fill="royalblue") +
                   ylab("Number of Promoters") +
                   theme_minimal() +
                   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(myplot, filename=paste0(prefix, "promoter_cov_", cov, "X_bar.png"))

}

#-------------------------------------------------------------------------------------------------------------------------#

############-------------------------------#
# FUNCTION # to produce CpG Coverage Table #
############-------------------------------#
# This function can primarily be used to filter the CpGs to a defined coverage cutoff
#  default is cov4 and perc4 (set to 10X and 99.9% respectively)
# if table is set to true (default), a CpG coverage table will be written to stdout

cpg_cov <- function(myObj=myObj, cov1=1, perc1=NULL, cov2=10, perc2=NULL, cov3=40, perc3=NULL, cov4=10, perc4=99.9, table=TRUE) {

  # obtain sample.ids from myObj
  samp_names <- lapply(myObj, function(x) getSampleID(x))

  # cov1
  myObj.filtered=filterByCoverage(myObj,lo.count=cov1,lo.perc=NULL,hi.count=NULL,hi.perc=perc1)
  my.cov1 <- sapply(myObj.filtered, function(x) nrow(x))

  # cov2
  myObj.filtered=filterByCoverage(myObj,lo.count=cov2,lo.perc=NULL,hi.count=NULL,hi.perc=perc2)
  my.cov2 <- sapply(myObj.filtered, function(x) nrow(x))

  # cov3
  myObj.filtered=filterByCoverage(myObj,lo.count=cov3,lo.perc=NULL,hi.count=NULL,hi.perc=perc3)
  my.cov3 <- sapply(myObj.filtered, function(x) nrow(x))

  # cov4
  myObj.filtered=filterByCoverage(myObj,lo.count=cov4,lo.perc=NULL,hi.count=NULL,hi.perc=perc4)
  my.cov4 <- sapply(myObj.filtered, function(x) nrow(x))

  if (table == TRUE) {
    # Create data frame of CpG coverages
    cpg_covs <- data.frame(sample=unlist(samp_names, use.names=FALSE),
                           V1=my.cov1,
                           V2=my.cov2,
                           V3=my.cov3,
                           V4=my.cov4
    )

    colnames(cpg_covs) <- c('sample',
                            paste('CpG', cov1, perc1, sep='.'),
                            paste('CpG', cov2, perc2, sep='.'),
                            paste('CpG', cov3, perc3, sep='.'),
                            paste('CpG', cov4, perc4, sep='.')
    )

    # Export CpG Coverage Table
    write.table(cpg_covs, "CpG.coverage.table.cov_files.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
  }

  return(myObj.filtered)
}

#---------------------------------------------------------------------------------------------------------------------------------------#

############-----------------------------------------#
# FUNCTION # to create coverage plots from methylKit #
############-----------------------------------------#
makeCovPlots <- function(myObj=myObj, basedir=getwd(), outdir="coverage_plots") {

  dir.create(file.path(basedir, outdir), showWarnings = FALSE)
  setwd(file.path(basedir, outdir))

  # obtain sample.ids from myObj
  samp_names <- lapply(myObj, function(x) getSampleID(x))

  for (i in 1:length(myObj)) {
      pdf(paste0(samp_names[[i]],"_covplot.pdf"))
      getCoverageStats(myObj[[i]], plot=TRUE, both.strands=FALSE)
      dev.off()
  }
  setwd(basedir)
}

#--------------------------------------------------------------------------------------------------------------------#

############--------------------------------------------------#
# FUNCTION # to create methylation stats plots from methylKit #
############--------------------------------------------------#

makeMethStats <- function(myObj=myObj, basedir=getwd(), outdir="meth_stats_plots") {

  dir.create(file.path(basedir, outdir), showWarnings = FALSE)
  setwd(file.path(basedir, outdir))

  # obtain sample.ids from myObj
  samp_names <- lapply(myObj, function(x) getSampleID(x))

  for (i in 1:length(myObj)) {
      pdf(paste0(samp_names[[i]],"_methstats.pdf"))
      getMethylationStats(myObj[[i]], plot=TRUE, both.strands=FALSE)
      dev.off()
  }
  setwd(basedir)
}

#--------------------------------------------------------------------------------------------------------------------#

############------------------------------------------------------------------------------------#
# FUNCTION # to cov filter, optional normalize by cov, merge on cpgs or tile and merge on tiles #
############------------------------------------------------------------------------------------#

# myObj
# min: minimum cpg cov number (default 10)
# max: maximum cpg cov percent (default 99.9)
# normalize: whether or not to normalize by coverage after filtering (default TRUE)
# tile_size: bp size of genomic tile region (default 1000)
# step_size: bp size of step between tiles (default 1000)
# min_cpg: remove regions with less than this number of covered cytosines (default 10)
# mpg: minimin per group, number of samples per group the tile must appear in to be kept after merging (default NULL)
# returns the meth object of merged tiled regions
get_merged_regions <- function(myObj=myObj, min=10, max=99.9, normalize=TRUE, regional=TRUE, tile_size=1000, step_size=1000, min_cpg=10, mpg=NULL, destrand=FALSE) {

  # filter cpgs
  print("Filtering CpGs with low or high coverage.")
  myObj.filtered <- filterByCoverage(myObj, lo.count=min, lo.perc=NULL, hi.count=NULL, hi.perc=max)

  # normalize if parameter is TRUE
  if (normalize) {
    print("Normalizing coverage.")
    myObj.filtered <- normalizeCoverage(myObj.filtered, method="median")
  }

  # if regional=TRUE, tile into regions and summarize
  if (regional) {
    # Summarize CpG methylation in tiled regions
    print("Summarizing CpG methylation per tiled region.")
    my.tiles <- tileMethylCounts(myObj.filtered,
                                 win.size=tile_size,
                                 step.size=step_size,
                                 cov.bases=min_cpg,
                                 mc.cores=8
    )
    for (i in 1:length(my.tiles)) {
      print(nrow(my.tiles[[i]]))
    }
    # merge the tiles
    print("Merging tiled regions across samples.")
    meth <- methylKit::unite(my.tiles, min.per.group=mpg, mc.cores=8)
    print(paste0("Returning ", nrow(meth), " regions."))

    return(meth)
  }

  # if not regional, just merge cpgs with no tiling
  else {
    print("merging CpGs")
    meth <- methylKit::unite(myObj.filtered,
                  min.per.group=mpg,
                  destrand=destrand)
    print(paste0("Returning ", nrow(meth), " CpGs."))
    return(meth)
  }
}

#----------------------------------------------------------------------------------------------------------------------------------#
