#!/usr/bin/env Rscript

source("/home/groups/hoolock2/u0/jtw/roselli/scripts/helperFunctions.R")
source("/home/groups/hoolock2/u0/jtw/roselli/scripts/methylkit_diff_functions.R")

library(methylKit)
library(ggplot2)
library(umap)
library(tidyr)
library(scales)
library(gridExtra)
library(impute)    

#Validate script inputs

args = commandArgs(trailingOnly=TRUE)
samples <- checkSamples()
config <- checkConfiguration()
checkInputs(args, config)

#Extract data and store human-readable variables

outputDirectory <- args[2]
inputDirectory <- args[1]
sampleFiles <- getSampleFiles(samples[,c(1)], inputDirectory)
treatment <- getTreatmentVector(samples[,c(2)])
sampleNames <- getList(samples[,c(1)])

#Create objects for data analysis 

myObj <- getObject(treatment, sampleFiles, sampleNames, config)
meth = getMergedRegions(myObj, config)

#group 1 vs group 2
c1 <- reorganize(meth,
		sample.ids=c("UMN_1_1","UMN_3_1","UMN_5_1","UMN_7_1","UMN_11_1","UMN_1_10","UMN_5_10","UMN_3_10","UMN_7_10","UMN_11_12"),
		treatment=c(0,0,0,0,0,1,1,1,1,1)
) 

covariates1 <- data.frame(individual = c("1","3","5","7","11","1","5","3","7","11"))

res1 <- calc.DMRs(c1,
		covariate=covariates1,
		overdispersion="MN",
		test="Chisq",
		comparison="PatientPre_v_PatientPost",
		meth.diff=10,
		qval=0.1,
		type="DMR",
		mc=8
)

makeBED(res1,"PatientPre_v_PatientPost","DMR")

#group 1 vs group 3
c2 <- reorganize(meth,
		sample.ids=c("UMN_1_1","UMN_3_1","UMN_5_1","UMN_7_1","UMN_11_1","UMN_2_1","UMN_6_1","UMN_4_1","UMN_8_1","UMN_12_1"),
		treatment=c(0,0,0,0,0,1,1,1,1,1)
)

covariates2 <- data.frame(individual = c("1","3","5","7","11","2","6","4","8","12"))

res2 <- calc.DMRs(c2,
                covariate=covariates2,
                overdispersion="MN",
                test="Chisq",
                comparison="PatientPre_v_ControlPre",
                meth.diff=10,
                qval=0.1,
                type="DMR",
                mc=8
)

makeBED(res2,"PatientPre_v_ControlPre","DMR")

#group 3 vs group 4
c3 <- reorganize(meth,
		sample.ids=c("UMN_2_1","UMN_6_1","UMN_4_1","UMN_8_1","UMN_12_1","UMN_2_10","UMN_6_10","UMN_4_10","UMN_8_10","UMN_12_12"),
		treatment=c(0,0,0,0,0,1,1,1,1,1)
)

covariates3 <- data.frame(individual = c("2","6","4","8","12","2","6","4","8","12"))

res3 <- calc.DMRs(c3,
                covariate=covariates3,
                overdispersion="MN",
                test="Chisq",
                comparison="ControlPre_v_ControlPost",
                meth.diff=10,
                qval=0.1,
                type="DMR",
                mc=8
)

makeBED(res3,"ControlPre_v_ControlPost","DMR")

#group 2 vs group 4F
c4 <- reorganize(meth,
		sample.ids=c("UMN_1_10","UMN_5_10","UMN_3_10","UMN_7_10","UMN_11_12","UMN_2_10","UMN_6_10","UMN_4_10","UMN_8_10","UMN_12_12"),
		treatment=c(0,0,0,0,0,1,1,1,1,1)
)

covariates4 <- data.frame(individual = c("1","5","3","7","11","2","6","4","8","12"))

res4 <- calc.DMRs(c4,
                covariate=covariates4,
                overdispersion="MN",
                test="Chisq",
                comparison="PatientPost_v_ControlPost",
                meth.diff=10,
                qval=0.1,
                type="DMR",
                mc=8
)

makeBED(res4,"PatientPost_v_ControlPost","DMR")

