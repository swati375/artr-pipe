 ARTR-seq Processing Workflow

This document describes how to run the preprocessing, alignment, and MACS3 peak-calling workflows included in this repository.

The workflow supports both single-end (SE) and paired-end (PE) sequencing data where indicated.

## Overview

The processing workflow consists of:

```text
FASTQ files
    ↓
FastQC
    ↓
Adapter trimming and UMI extraction
    ↓
rRNA filtering
    ↓
STAR alignment
    ↓
Normalization for IGV
    ↓
MACS3 peak calling
```

# 1. Preprocessing and alignment

All preprocessing and alignment scripts are located in:

```text
preprocess-and-align/
```

## Input directory structure

Before starting, arrange the raw sequencing files should be placed in a separate folder inside the project directory, as follows:

```text
project_directory/
└── fastq-files/
    ├── sample1.fastq.gz
    ├── sample2.fastq.gz
    └── ...
```
For paired-end data, both R1 and R2 files must be present.

## Step 1: Checking Quality of squencing files using FASTQ

Run:

```bash
bash preprocess-and-align/fastqc.sh
```

The script will ask for:

1. Project Base directory
2. File containing the FASTQ filenames to analyze (example preprocess-and-align/examples/list-fastq.txt)
3. Output directory name for FastQC results
4. Directory containing the input FASTQ files (should contain all FASTQ files from file 2.)

The output directory contains the FastQC results and the MultiQC summary.

---

## Step 2: Sequence Preprocessing

Run:

```bash
bash preprocess-and-align/trim-filter-reads.sh
```

The script will ask whether the sequencing data are paired-end: 

```text
TRUE
```

or single-end:

```text
FALSE
```

It will then ask for the Project base directory.

The script performs: Adapter trimming, UMI extraction 

For paired-end data, the script automatically identifies the matching R2 file from the R1 filename and checks that the R2 file exists.

This script then automatically runs rRNA filtering, genome alignment (using STAR) and read deduplication on your data.

The resulting files are written to:

```text
project_directory/trimmed-files/
project_directory/staralign-bam-files/
```

### Important

The STAR genome index must be configured correctly in:

```text
preprocess-and-align/align-dedup.sh
```

Please also note that the script is currently configured for the human genome and would need to be modified if a different organism is being used.


## Step 3: Generate BigWig files for viaulization

After preprocessing and alignment are complete, run:

```bash
bash preprocess-and-align/bamtobed.sh
```

The script asks for the project base directory and generates CPM-normalized BigWig coverage files in:

```text
project_directory/crosslink_nt/
```

One BigWig file is generated for each sample.


# 2. MACS3 peak calling

For Peak-calling, MACS3 peak caller has been used and options chosen to run with RNA reads. Check MACS3 documentation (https://macs3-project.github.io/MACS/) for details and different parameters)

## Sample ID file

Peak calling uses a sample identification file to associate IP samples with their corresponding input controls.

An example file is provided at:

```text
peak-calling-macs3/examples/sample-id.txt
```

For example:

```text
WT:48966_WT
Input_WT:48967_WT_input
protAKO:49197_protAKO
Input_protAKO:49200_protAKO_input
```

The input sample naming is important.

Each IP sample:

```text
WT 
protAKO
DKO
```

will have its corresponding input named:

```text
Input_WT
Input_protAKO
Input_DKO
```

MACS3 can also find peaks if no Input/ control is provided by user, however I would highly recommend it for getting clear signal enrichment over background.

The paired-end (PE) data can either be used as PE or as single-end (SE). Script for both are available. To use same parameters as in the ARTR-seq paper https://pubmed.ncbi.nlm.nih.gov/38200227/,
use:

---

## peak calling

Run the single-end workflow using:

```bash
python peak-calling-macs3/macs3_calling_se.py
```

Run the paired-end workflow using:

```bash
python peak-calling-macs3/macs3_calling_pe.py
```

When prompted, provide the path to the sample ID file.

For example:

```text
peak-calling-macs3/examples/sample-id.txt
```

The script will ask whether to run peak calling for all proteins or for a selected protein from your sample ID file.

---

# Configuration requirements

Before running the workflow on a new system, check the reference index paths in:

These currently require the appropriate:

1. Bowtie2 rRNA index

This can be created one time for human genome using:

```bash
mkdir rrna-files
esearch -db nucleotide -query "NR_003285.3 OR NR_003286.4 OR NR_003287.4 OR NR_023363.1" | efetch -format fasta > human_rRNA_refs.fasta
bowtie2-build human_rRNA_refs.fasta human_rRNA_index
```

2. STAR genome index (got GRCh38 are provided in genome_files folder, any other can be downloaded from gencode/ USCS browser)
to be available on the system. 

The paths should be configured for the reference genome being analyzed. They are currently working for human GRCh38. The paths can be edited in script align-dedup.sh

## Software environment

The recommended Conda environment is provided in:

```text
environment.yml
```

Create and activate it using:

```bash
conda env create -f environment.yml
conda activate seq-py312
```