library(methylKit)

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

#check input directory exists
#Enumerate list of input files

#check if output directory exists
#create if not

#check if input files are present
#check if input files are accessible?
#store input file names

#create list of sample names from file names

#create data object
#   include all input parameters as potential command line inputs

#if (plots)
    #Loop over samples
        #coverage plot
            #create
            #store
        #methylation plot
            #create
            #store
#if (tables)
    #Loop over samples
        #coverage table
            #create
            #store
        #methylation table
            #create
            #store

#merge samples

#if (tables)
    #correlation table
    #cluster sample table

#if (Plots)
    #cluster sample ploit
    #PCA scree
    #PCA scatter plot
