#!/usr/bin/env Rscript

source("scripts/helperFunctions.R")

#Validate script inputs
args = commandArgs(trailingOnly=TRUE)
samples <- checkSamples()
config <- checkConfiguration()

outputDirectory <- args[2]
inputDirectory <- args[1]
sampleFiles <- getSampleFiles(samples[,c(1)], inputDirectory)
treatment <- getTreatmentVector(samples[,c(2)])
sampleNames <- getList(samples[,c(1)])


myObj <- getObject(treatment, sampleFiles, sampleNames, config)
meth = getMergedRegions(myObj, config)

#DMR 
res <- calc.DMRs(meth,
                 covariate=NULL,
                 overdispersion="MN",
                 test="Chisq",
                 comparison="TS_BAV_v_TS",
                 meth.diff=10,
                 qval=0.1,
                 type="DMR",
		 mc=1
)


#Bed files
makeBED(res, "TS_BAV_v_TS", "DMR")