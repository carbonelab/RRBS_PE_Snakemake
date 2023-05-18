#!/bin/bash

#move into data directory
cd $1

mkdir homer
cd homer
mkdir inputs
cp ../dmr/*/*.sigDMRs.bed inputs/.

dir="./inputs"
genome_path="/home/groups/hoolock2/u0/genomes/ensembl/homo_sapiens/primary_assembly/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"

for file in $dir/*.bed; do

    name=${file##*/}
    dirname=${name##*/}
    outdir="${dirname%.*}"

    mkdir $outdir

    findMotifsGenome.pl $file $genome_path $outdir -cpg -size given -len 8,10,12 -p 8;
done
