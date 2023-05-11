#!/bin/bash

join data/fastqc/raw/fqc_stats.table.txt data/trimming/trimgalore_stats.txt > data/tempMetrics.txt

join data/tempMetrics.txt data/bismark_aln/bismark_stats.txt > data/metrics_summary.txt

sed -e 's/\s/,/g' data/metrics_summary.txt > data/metrics_summary.csv

unix2dos data/metrics_summary.csv

ssconvert data/metrics_summary.csv data/metrics_summary.xlsx

rm -rf data/tempMetrics.txt
rm -rf data/metrics_summary.txt
rm -rf data/metrics_summary.csv
