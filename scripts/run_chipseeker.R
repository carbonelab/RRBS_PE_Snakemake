

# conda activate chipseeker

library(GenomicFeatures)
library(ChIPseeker)
library(writexl)

#-----------------------------#

# INPUT VARIABLES

# input a gtf file
gtf_file <- paste0("/home/groups/hoolock2/u0/genomes/ensembl/",
  "homo_sapiens/primary_assembly/annotation/Homo_sapiens.GRCh38.103.gtf")

# gene info file
gene_info_file <- paste0("/home/groups/hoolock2/u0/genomes/ensembl/",
  "homo_sapiens/primary_assembly/annotation/GRCh38.103.gene_info.txt")

# sig dmr results txt file
sig_dmrs_txt <- args[2]

#paste0("/home/groups/hoolock2/u0/jtw/newell/data/",
 # "dmr/comparison4/PatientPost_v_ControlPost.sigDMRs.txt")

# comparison name (for naming output files)
comp_name <- args[1]

#"PatientPost_v_ControlPost"

setwd(paste0("data/dmr/",args[1]))

#-----------------------------#

# make a TxDB object from the annotation gtf file
txdb <- makeTxDbFromGFF(gtf_file,
  format="gtf", "GRCh38.103", "Homo sapiens")

# read in DMR results to annotate
mydat <- makeGRangesFromDataFrame(read.delim(sig_dmrs_txt, header=T),
  keep.extra.columns=TRUE, seqnames.field="chr", start.field="start",
  end.field="end", ignore.strand=T)

# annotate via Chipseeker
res <- annotatePeak(mydat,
         tssRegion = c(-3000, 200),
         TxDb = txdb,
  	 level = "transcript",
  	 assignGenomicAnnotation = TRUE,
  	 overlap="TSS"
)
# convert to df
foo <- as.data.frame(res)

# combine with bioinfo data to get gene info
ref <- read.delim(gene_info_file, header=T)
colnames(ref)[1] <- "geneId"
foo2 <- merge(foo, ref, by="geneId")
#  clean cols
foo3 <- foo2[ ,c(7,2,3,4,5,6,8,9,10,11,12,13,14,15,16,1,17,18,19,20,21)]

# export txt and excel file
# define outfile names from input results txt file name
write.table(foo3,
  gsub(".sigDMRs.txt", ".sigDMRs.annot.txt", sig_dmrs_txt),
  sep="\t",
  col.names=T,
  row.names=F,
  quote=F)

# export excel file
write_xlsx(foo3,
  gsub(".sigDMRs.txt", ".sigDMRs.annot.xlsx", sig_dmrs_txt))

#---------------#
# Plots

# anno pie
pdf(paste0(comp_name, ".anno_pie_chart.pdf"))
plotAnnoPie(res)
dev.off()

# vennpie
pdf(paste0(comp_name, ".vennPie.pdf"))
#par(mar=c(0.5,2,1,3))
vennpie(res)
dev.off()
