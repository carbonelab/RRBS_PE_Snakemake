

# Functions:
# 1. calc.DMRs
# 2. makeBED
# 3. IN PROGRESS (maybe unnecessary?): remove zero variance rows from methylbase object

############--------------------------------------#
# FUNCTION # to calculate DMRs and export results #
############--------------------------------------#
# will utilize methylKit's "MN" overdispersion and "Chisq" test
calc.DMRs <- function(my.meth, covariate=NULL, overdispersion="MN", test="Chisq", comparison, meth.diff=10, qval=0.1, type="DMR", mc=8) {
  options(scipen=999)
  
  print("Calculating DMRs.")
  # Calculate DMRs: Overdispersion:YES, Test:Chisq
  myDiff <- calculateDiffMeth(my.meth,
                              covariates=covariate,
                              overdispersion=overdispersion,
                              test=test,
                              mc.cores=mc
  )
  # convert results to a data frame
  myDiff.df <- getData(myDiff)

  print(paste0("nrow of unfiltered diff results is: ", nrow(myDiff.df)))

  # Plot p-value distribution
  png(paste0(comparison, ".pval_dist.png"))
  hist(myDiff$pvalue)
  dev.off()

  # Subset to significant DMRs
  print(paste0("Number of sig ", type, "s: ", nrow(getMethylDiff(myDiff, difference=meth.diff, qvalue=qval))))
  myDiff.sig <- getMethylDiff(myDiff, difference=meth.diff, qvalue=qval)
  # convert to data frame
  myDiff.sig.df <- getData(myDiff.sig)
  print(paste0("nrow of myDiff.sig.df is: ", nrow(myDiff.sig.df)))

  if (nrow(myDiff.sig.df) > 0) {
    #----------------------------#
    # HYPO and HYPER SEPARATIONS #
    #----------------------------#

    # separate hyper and hypo sig DMRs
    myDiff.hyper <- myDiff.sig.df[myDiff.sig.df$meth.diff > 0, ]
    myDiff.hypo <- myDiff.sig.df[myDiff.sig.df$meth.diff < 0, ]

    print(paste0("number of hyper ", type, "s: ", nrow(myDiff.hyper)))
    print(paste0("number of hypo ", type, "s: ", nrow(myDiff.hypo)))

    # df for pie chart of hyper and hypo sig DMRs
    pie.df <- data.frame(class=c("hyper", "hypo"),
                         percent=c(nrow(myDiff.hyper) / nrow(myDiff.sig.df),
                                   nrow(myDiff.hypo) / nrow(myDiff.sig.df))
    )
    pie.df$label=paste0(round(as.numeric(pie.df$percent)*100, digits=0), '%')

    # define blank theme for prettier pie chart
    blank_theme <- theme_minimal()+
    theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
    )

    # plot ratio of hyper vs hypo sig DMRs
    myplot <- ggplot(pie.df, aes(x="", y=percent, fill=class)) +
               geom_bar(width = 1, stat = "identity") +
               coord_polar("y", start=0) +
               blank_theme +
               theme(axis.text.x=element_blank()) +
               #geom_text(aes(x=1, y=pos, label=label), size = 6) +
               geom_text(aes(x=1, y=cumsum(percent/sum(percent)) - percent/sum(percent)/2,
                 label=label[c(2,1)]), size=6) +  #label = paste(round(percent/sum(percent)*100),"%")), size = 6)
               ggtitle(paste0("Ratio of Hyper and Hypo Methylated sig ", type, "s\ntotal: ", nrow(myDiff.sig.df)))
    ggsave(filename=paste0(comparison, ".hyper.hypo.sig.", type, ".pieChart.png"))

    # -------------------------------#
    # continue with both hyper and hypo export

    # add DMR_ID column
    myDiff.sig.df$unq_ID <- paste(type, seq(1:nrow(myDiff.sig.df)), sep='_')

    # reorder columns
    myDiff.sig.df <- myDiff.sig.df[ ,c(8,1,2,3,5,6,7)]

    # export sig results
    write.table(myDiff.sig.df,
                paste0(comparison, ".sig", type, "s.txt"),
                sep="\t",
                col.names=TRUE,
                row.names=FALSE,
                quote=FALSE
    )

    #--------------#
    # sig dmr beds #
    #--------------#

    # BED file containing all sig DMRs and ensembl contig format
    # subset sig DMR results to a bed file for annotation
    res <- myDiff.sig.df[ ,c(2,3,4,1)]
    res$start <- res$start - 1

    # export sig DMR bed for annotation
    write.table(res,
                paste0(comparison, ".sig", type, "s.bed"),
                sep="\t",
                col.names=FALSE,
                row.names=FALSE,
                quote=FALSE
    )

    # 3 BED files for GREAT
    # UCSC contig format, 1 bed file of hyper only, 1 of hypo only, 1 of both
    res$chr <- paste0("chr", res$chr)
    write.table(res,
                paste0(comparison, ".sig", type, "s.both.GREAT.bed"),
                sep="\t",
                col.names=FALSE,
                row.names=FALSE,
                quote=FALSE
    )

    my.hyper <- myDiff.sig.df[myDiff.sig.df$meth.diff > 0, c(2,3,4,1)]
    # check if there are 0 sig hyper DMRs
    if (nrow(my.hyper) > 0) {
      my.hyper$start <- my.hyper$start - 1
      my.hyper$chr <- paste0("chr", my.hyper$chr)
      write.table(my.hyper,
                  paste0(comparison, ".sig", type, "s.hyper.GREAT.bed"),
                  sep="\t",
                  col.names=FALSE,
                  row.names=FALSE,
                  quote=FALSE
      )
    }
    else if (nrow(my.hyper) < 1) {
      print("Zero sig hyper DMRs for GREAT.")
    }

    my.hypo <- myDiff.sig.df[myDiff.sig.df$meth.diff < 0, c(2,3,4,1)]
    # check if there are 0 sig hypo DMRs
    if (nrow(my.hypo) > 0) {
      my.hypo$start <- my.hypo$start - 1
      my.hypo$chr <- paste0("chr", my.hypo$chr)
      write.table(my.hypo,
                  paste0(comparison, ".sig", type, "s.hypo.GREAT.bed"),
                  sep="\t",
                  col.names=FALSE,
                  row.names=FALSE,
                  quote=FALSE
      )
    }
    else if (nrow(my.hypo) < 1) {
      print("Zero sig hypo DMRs for GREAT.")
    }
  }
  else {
    print(paste0("There are zero ", type, "s passing signficance thresholds."))
  }
  return(myDiff.sig.df)
}


############---------------------------------------------------------#
# FUNCTION # to make a sig DMR bed file track from results dataframe #
############---------------------------------------------------------#

makeBED <- function(res, comparison, type="DMR") {
  res2 <- res[ ,c(2,3,4,1,7)]
  res2$start <- res2$start - 1
  res2$chr <- paste0("chr", res2$chr)

  res2[ ,6] <- "."
  res2[ ,7] <- res2$start
  res2[ ,8] <- res2$end
  res2[ ,9] <- ifelse(res2[ ,5] > 0, '255,0,0', ifelse(res2[ ,5] < 0, '0,0,255', '0,0,0'))

  # add track line to exported file
  cat(paste0("track type=bed name=", comparison, " ", "itemRgb=On"), file=paste0(comparison, ".sig", type, "track.bed"))
  cat("\n", file=paste0(comparison, ".sig", type, "track.bed"), append=TRUE)

  # export bed info
  write.table(res2,
              paste0(comparison, ".sig", type, "track.bed"),
              sep="\t",
              col.names=FALSE,
              row.names=FALSE,
              quote=FALSE,
              append=TRUE
  )
}

############---------------------------------------------------------#
# FUNCTION # to remove zero variance rows from methylbase object     #
############---------------------------------------------------------#

# Useful to lower the number of tests performed during diff analysis

rm_zero_var <- function(methylbase_obj) {
  # extract methylbase info to df
  df <- getData(methylbase_obj)

  # define number of samples based on number of columns
  n_samp <- (ncol(df) - 4) / 3

  # calculate percent meth from numCs and coverage cols per sample
  # ignore first 4 cols
  perc_meth_df <- df[, -c(1,2,3,4)]
  # create vector of numCs columns
  numc_vec <- seq(2, ncol(perc_meth_df)-1, by=3)
  # create vector of coverage columns
  cov_vec <- seq(1, ncol(perc_meth_df)-2, by=3)
  # loop through each sample and take numc over coverage to add percMeth column
  for (i in 1:n_samp) {
    new_col_name <- paste0("percMeth", i)
    perc_meth_df[[new_col_name]] <- perc_meth_df[, numc_vec[i]] / perc_meth_df[, cov_vec[i]]
  }

  # determine zero var rows
  perc_meth_df2 <- perc_meth_df[apply(perc_meth_df, 1, var) != 0, ]
  #perc_meth_df2 <- perc_meth_df[apply(perc_meth_df, 1, var) >=0.3, ]
  # on my current rrbs dataset of 160,118 tiles, neither of these removed any rows
}