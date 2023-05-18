#!/bin/bash

cd $1

join fastqc/raw/fqc_stats.table.txt trimming/trimgalore_stats.txt > tempMetrics.txt

join tempMetrics.txt bismark_aln/bismark_stats.txt > metrics_summary.txt

sed -e 's/\s/,/g' metrics_summary.txt > metrics_summary.csv

unix2dos metrics_summary.csv

ssconvert metrics_summary.csv metrics_summary.xlsx

rm -rf tempMetrics.txt
rm -rf metrics_summary.txt
rm -rf metrics_summary.csv
