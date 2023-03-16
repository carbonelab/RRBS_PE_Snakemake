#!/bin/bash

# Purpose:

# Input variables:
#   - fqc_path: path to folder containing FastQC results
#   - outdir: desired output directory of the summary results file
#             defaults to fqc_path

fqc_path=$1

#-----------------------#
# Get Stats from FastQC #
#-----------------------#
cd $fqc_path
for f in `ls --color=none *.zip`; do unzip $f; done
# Grab Metrics
grep "Total Sequences" */fastqc_data.txt > tmp_reads.txt
grep "Sequence length" */fastqc_data.txt > tmp_readLength.txt
grep "Total Deduplicated" */fastqc_data.txt > tmp_duplicated.txt
grep "^%GC" */fastqc_data.txt > tmp_GC.txt
# Clean up Metrics
### gather the R2 entries into a separate file
awk 'NR%2==1' tmp_reads.txt > tmp_readsR1.txt
awk 'NR%2==0' tmp_reads.txt > tmp_readsR2.txt
awk 'NR%2==1' tmp_readLength.txt > tmp_readLengthR1.txt
awk 'NR%2==0' tmp_readLength.txt > tmp_readLengthR2.txt
awk 'NR%2==1' tmp_GC.txt > tmp_GCR1.txt
awk 'NR%2==0' tmp_GC.txt > tmp_GCR2.txt
### transform the deduplicated percentage to duplicated percentage (by subtracting from 100)
awk -F'\t' '{print $1,(100-$2)}' tmp_duplicated.txt > tmp_duplicated2.txt
### gather the R2 entries into a separate file
awk 'NR%2==1' tmp_duplicated2.txt > tmp_duplicatedR1.txt
awk 'NR%2==0' tmp_duplicated2.txt > tmp_duplicatedR2.txt
### Keep only the fastq prefix name in the first field of the tmp_readsR1.txt file
sed -i 's/\_R1\_fastqc\/fastqc_data\.txt\:Total Sequences/''/g' tmp_readsR1.txt
### Keep only the last field for remaining tmp R1,R2 files
cut -f2 tmp_readsR2.txt > tmp2_readsR2.txt
cut -f2 tmp_readLengthR1.txt > tmp2_readLengthR1.txt
cut -f2 tmp_readLengthR2.txt > tmp2_readLengthR2.txt
cut -f2 tmp_GCR1.txt > tmp2_GCR1.txt
cut -f2 tmp_GCR2.txt > tmp2_GCR2.txt
cut -d ' ' -f4 tmp_duplicatedR1.txt > tmp2_duplicatedR1.txt
cut -d ' ' -f4 tmp_duplicatedR2.txt > tmp2_duplicatedR2.txt
# Combine FQC stats into one table
paste -d "\t" tmp_readsR1.txt tmp2_readsR2.txt tmp2_readLengthR1.txt tmp2_readLengthR2.txt tmp2_GCR1.txt tmp2_GCR2.txt tmp2_duplicatedR1.txt tmp2_duplicatedR2.txt > tmp_fqc_stats.tmp.txt
sort tmp_fqc_stats.tmp.txt > fqc_stats.table.txt

# add a header
sed -i '1s/^/Sample_Name\tR1_raw_reads\tR2_raw_reads\tR1_read_length\tR2_read_length\tR1_GC%\tR2_GC%\tR1_dup%\tR2_dup%\n/' fqc_stats.table.txt
rm tmp*.txt

# Move output file and change directories
#cd $dir
#cd ..
#mv $fqc_path/fqc_stats.table.txt .
