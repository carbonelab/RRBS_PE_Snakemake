#---------------------------------------------------------------------------#
# Purpose: to parse relevant info from a directory of bismark log files     #
#          and produce an outfile containing a table of stats per sample    #
#                                                                           #
# Useage: ./parse.bismark.logs.py -d /path/to/reports -o path/to/outfile    #
#---------------------------------------------------------------------------#

# Notes:
#       - Assumes the bismark log files end with "report.txt"
#       - Written for PE files (looks for "Sequence pairs analysed in total" line)
#       - Don't need path to trimmed files anymore

# TO DO: - Make it PE or SE agnostic? Grab for something that is in both reports
#       - OR make a PE version and a SE version


# Import libraries
import optparse
import pandas as pd
import re
import os


# Initialize command line options for input directory and output file
p=optparse.OptionParser()
p.add_option("-d", action = "store", dest = "directory")
p.add_option("-o", action = "store", dest = "outfile")

opts,args=p.parse_args()
in_dir=opts.directory
outfile=opts.outfile

# open outfile for writing
fhw = open(outfile, "w+")
# write colnames to outfile
fhw.write('Sample_Name' + '\t' + 'Input_Reads' + '\t' + 'Unique_Alignments' + '\t' + 'Unq_Alignment_Rate' + '\t' + 'No_Alignments' + '\t' + 'Multimapped' + '\t' + 'Discarded' + '\t' + 'Total_Cs_Analyzed' + '\t' + 'CpG_Context' + '\t' + 'CHG_Context' + '\t' + 'CHH_Context' + '\t' + 'Unknown_Context' + '\n')

# Loop through files in input directory, and work on files ending in *report.txt
#   to grab lines of interest and record the metric to an output file
for file in sorted(os.listdir(in_dir)):
    filename = os.fsdecode(file)
    if filename.endswith("report.txt"):
        path_file=in_dir + '/' + filename
        theFile = open(path_file,'r')
        FILE = theFile.readlines()
        theFile.close()
        printList = []
        for line in FILE:
            # Obtain fastq file name (trimmed fq name, r1 if pe)
            if ('Bismark report for' in line):
                intname=line.split()[3]
                intname2=intname.split('/')[-1]
                intname3=re.sub("(_val_1.*)$", "", intname2)
                printList.append(intname3 + '\t')
            # Obtain Number of Seq Pairs Analyzed
            if ('Sequence pairs analysed in total' in line):
                intname=line.split('\t')[1].rstrip()
                printList.append(intname + '\t')
            # Obtain number of unique best hit alignments
            if ('unique best hit' in line):
                intname=line.split('\t')[1].rstrip()
                printList.append(intname + '\t')
            # Obtain mapping percentage
            if ('Mapping efficiency' in line) or ('Paleontologist' in line):
                intname=line.split('\t')[1].rstrip()
                printList.append(intname + '\t')
            # Get no aln
            if ('no alignments under any condition' in line):
                intname=line.split('\t')[1].rstrip()
                printList.append(intname + '\t')
            # Get multimap
            if ('not map uniquely:' in line):
                intname=line.split('\t')[1].rstrip()
                printList.append(intname + '\t')
            # Get discarded
            if ('genomic sequence could not be extracted' in line):
                intname=line.split('\t')[1].rstrip()
                printList.append(intname + '\t')
            # Get Total C's analyzed
            if ('Total number of C' in line):
                intname=line.split('\t')[1].rstrip()
                printList.append(intname + '\t')
            # Get CpG Meth
            if ('C methylated in CpG context' in line):
                intname=line.split('\t')[1].rstrip()
                printList.append(intname + '\t')
            # Get CHG Meth
            if ('C methylated in CHG context' in line):
                intname=line.split('\t')[1].rstrip()
                printList.append(intname + '\t')
            # Get CHH Meth
            if ('C methylated in CHH context' in line):
                intname=line.split('\t')[1].rstrip()
                printList.append(intname + '\t')
            # Get Unknown (CN or CHN) Meth
            if ('C methylated in unknown context' in line):
                intname=line.split('\t')[1]
                printList.append(intname)

        for item in printList:
            fhw.write(item)

    else:
        continue

fhw.close()
