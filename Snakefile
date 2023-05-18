# Snakefile to analyze RRBS PE data
# 

configfile:"proj_config.yaml"
#project_id = config["project_id"]


SAMPLES, = glob_wildcards("data/fastq/{sample}_R1.fastq.gz")

localrules: collect_fqc_metrics, collect_trimgalore_metrics, collect_bismark_metrics

rule all:
    input:
        expand("data/fastqc/raw/{sample}_{dir}_fastqc.zip", sample = SAMPLES, dir = ["R1", "R2"]),
        expand("data/trimming/{sample}_val_{dir}.fq.gz", sample = SAMPLES, dir = ["1", "2"]),
        expand("data/fastqc/trim/{sample}_val_{dir}_fastqc.zip", sample = SAMPLES, dir = ["1", "2"]),
        expand("data/bismark_aln/{sample}_val_1_bismark_bt2_pe.bam", sample = SAMPLES),
        "data/fastqc/raw/fqc_stats.table.txt",
        "data/trimming/trimgalore_stats.txt",
        "data/bismark_aln/bismark_stats.txt",
        "data/metrics_summary.xlsx",
        expand("data/meth_extract/{sample}_val_1_bismark_bt2_pe.CpG_report.txt.gz", sample = SAMPLES),
        expand("data/ide/meth_stats_plots/{sample}_methstats.pdf", sample = SAMPLES),
        expand("data/ide/cov_stats_plots/{sample}_covstats.pdf", sample = SAMPLES),
        expand("data/ide/meth_stats/{sample}_methstats.txt", sample = SAMPLES),
        expand("data/ide/cov_stats/{sample}_covstats.txt", sample = SAMPLES),
        "data/ide/merged_stats/correlation.txt",
        "data/ide/merged_stats_plots/clusteringDendro.pdf",
        "data/ide/merged_stats_plots/pcaScree.pdf",
        "data/ide/merged_stats_plots/pcaScatter.pdf",
        expand("data/dmr/{comparison}/{comparison}.sigDMRs.bed", comparison = config["comparisons"]),
        expand("data/homer/{comparison}.sigDMRs/knownResults.html", comparison = config["comparisons"])

rule fastqc_raw:
    input:
        fwd = "data/fastq/{sample}_R1.fastq.gz",
        rev = "data/fastq/{sample}_R2.fastq.gz"
    output:
        fwd = "data/fastqc/raw/{sample}_R1_fastqc.zip",
        rev = "data/fastqc/raw/{sample}_R2_fastqc.zip"
    conda:
        "envs/fastqc.yaml"
    params:
        outdir = "data/fastqc/raw"
    shell:
        "fastqc -o {params.outdir} {input.fwd} {input.rev}"

rule trim_galore:
    input:
        fwd = "data/fastq/{sample}_R1.fastq.gz",
        rev = "data/fastq/{sample}_R2.fastq.gz"
    output:
        "data/trimming/{sample}_R2.fastq.gz_trimming_report.txt",
        fwd = "data/trimming/{sample}_val_1.fq.gz",
        rev = "data/trimming/{sample}_val_2.fq.gz"
    conda:
        "envs/trimgalore.yaml"
    params:
        basename = "{sample}",
        outdir = "data/trimming",
    shell:
        "trim_galore --paired --basename {params.basename} -o {params.outdir} --rrbs {input.fwd} {input.rev}"

rule fastqc_trim:
    input:
        fwd = "data/trimming/{sample}_val_1.fq.gz",
        rev = "data/trimming/{sample}_val_2.fq.gz"
    output:
        fwd = "data/fastqc/trim/{sample}_val_1_fastqc.zip",
        rev = "data/fastqc/trim/{sample}_val_2_fastqc.zip"
    conda:
        "envs/fastqc.yaml"
    params:
        outdir = "data/fastqc/trim"
    shell:
        "fastqc -o {params.outdir} {input.fwd} {input.rev}"

rule bismark_aln:
    input:
        fwd = "data/trimming/{sample}_val_1.fq.gz",
        rev = "data/trimming/{sample}_val_2.fq.gz"
    output:
        "data/bismark_aln/{sample}_val_1_bismark_bt2_PE_report.txt",
        bam = "data/bismark_aln/{sample}_val_1_bismark_bt2_pe.bam"
    conda:
        "envs/bismark.yaml"
    params:
        genome_dir = config["bismark_ref_genome"],
	    outdir = "data/bismark_aln"
    shell:
        "bismark -p 4 {params.genome_dir} -1 {input.fwd} -2 {input.rev} -o {params.outdir} --bam"

rule collect_fqc_metrics:
    input:
        expand("data/fastqc/raw/{sample}_{dir}_fastqc.zip", sample = SAMPLES, dir = ["R1", "R2"])
    output:
        "data/fastqc/raw/fqc_stats.table.txt"
    params:
        inpath = "data/fastqc/raw"
    shell:
        "scripts/collect_fastqc_metrics_PE.sh {params.inpath}"

rule collect_trimgalore_metrics:
    input:
        expand("data/trimming/{sample}_R2.fastq.gz_trimming_report.txt", sample = SAMPLES)
    output:
        "data/trimming/trimgalore_stats.txt"
    conda:
        "envs/python3_general.yaml"
    params:
        inpath = "data/trimming",
        outfile = "data/trimming/trimgalore_stats.txt"
    shell:
        "python scripts/parse.trimgalore.rrbs.pe.logs.py -d {params.inpath} -o {params.outfile}"

rule collect_bismark_metrics:
    input:
        expand("data/bismark_aln/{sample}_val_1_bismark_bt2_PE_report.txt", sample = SAMPLES)
    output:
        "data/bismark_aln/bismark_stats.txt"
    conda:
        "envs/python3_general.yaml"
    params:
        inpath = "data/bismark_aln",
        outfile = "data/bismark_aln/bismark_stats.txt"
    shell:
        "python scripts/parse.bismark.pe.logs.py -d {params.inpath} -o {params.outfile}"

rule summarize_metrics:
    input:
        "data/bismark_aln/bismark_stats.txt",
        "data/trimming/trimgalore_stats.txt",
        "data/fastqc/raw/fqc_stats.table.txt"
    output:
        "data/metrics_summary.xlsx"
    params:
        inpath = "data"
    shell:
        "scripts/combine_metrics.sh {params.inpath}"

rule meth_extract:
    input:
        bam = "data/bismark_aln/{sample}_val_1_bismark_bt2_pe.bam"
    output:
        "data/meth_extract/{sample}_val_1_bismark_bt2_pe.CpG_report.txt.gz",
        "data/meth_extract/{sample}_val_1_bismark_bt2_pe.bismark.cov.gz",
        "data/meth_extract/{sample}_val_1_bismark_bt2_pe.bedGraph.gz",
        "data/meth_extract/{sample}_val_1_bismark_bt2_pe.M-bias.txt"
    conda:
        "envs/bismark.yaml"
    params:
        genome_dir = config["bismark_ref_genome"],
        outdir = "data/meth_extract"
    shell:
        "bismark_methylation_extractor -p --comprehensive --merge_non_CpG --bedGraph --cytosine_report --gzip --genome_folder {params.genome_dir} -o {params.outdir} {input.bam}" 

rule ide:
    input:
        cov = expand("data/meth_extract/{sample}_val_1_bismark_bt2_pe.bismark.cov.gz", sample = SAMPLES)
    output:
        expand("data/ide/meth_stats_plots/{sample}_methstats.pdf", sample = SAMPLES),
        expand("data/ide/cov_stats_plots/{sample}_covstats.pdf", sample = SAMPLES),
        expand("data/ide/meth_stats/{sample}_methstats.txt", sample = SAMPLES),
        expand("data/ide/cov_stats/{sample}_covstats.txt", sample = SAMPLES),
        "data/ide/merged_stats/correlation.txt",
        "data/ide/merged_stats_plots/clusteringDendro.pdf",
        "data/ide/merged_stats_plots/pcaScree.pdf",
        "data/ide/merged_stats_plots/pcaScatter.pdf"
    conda:
        "envs/methylKit.yaml"
    params:
        outdir = "data/ide",
        inpath = "data/meth_extract/"
    shell:
        "Rscript scripts/methylKitIDE.R {params.inpath} {params.outdir}"

rule dmr:
    input:
        cov = expand("data/meth_extract/{sample}_val_1_bismark_bt2_pe.bismark.cov.gz", sample = SAMPLES)
    output:
        expand("data/dmr/{comparison}/{comparison}.sigDMRs.bed", comparison = config["comparisons"])
    conda:
        "envs/methylKit.yaml" #create this yaml
    params:
        outdir = "data/dmr", 
        inpath = "data/meth_extract/", 
        comparison = expand("data/dmr/{comparison}", comparison = config["comparisons"])
    shell:
        "Rscript scripts/methylKitDMR.R {params.inpath} {params.outdir} {params.comparison}" #needs to autogenerate comparison reorganization

rule chipseeker:
    input:
        dmrs = expand("data/dmr/{comparison}/{comparison}.sigDMRs.txt", comparison = config["comparisons"])
    output:
        expand("data/dmr/{comparison}/{comparison}.sigDMRs.annot.xlsx", comparison = config["comparisons"])
    conda:
        "envs/chipseeker.yaml" #create this yaml
    params:
        comparisons = expand("{comparison}", comparison = config["comparisons"])
    shell:
        "Rscript run_chipseeker.R {params.comparisons} {input.dmrs}"

rule homer:
    input:
        bed = expand("data/dmr/{comparison}/{comparison}.sigDMRs.bed", comparison = config["comparisons"])
    output:
        expand("data/homer/{comparison}.sigDMRs/knownResults.html", comparison = config["comparisons"])
    conda:
        "envs/homer.yaml" #create this yaml
    params:
        inpath = "data"
    shell:
        "scrips/run_homer_findMotifsGenome.sh {params.inpath}"


