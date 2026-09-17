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

# Configuration requirements

Before running the workflow on a new system, check the reference index paths in:

These currently require the appropriate:

1. Genome files
All genome fasta and annotation must be downloaded as required. These will be required for indexing, alignment and peak annotation. 

Download the files below and place everything in the genome_files folder under correct genome. All files are required for the genome of choice. As an example a dummy empty folder for GRCh38 has been created where you can place all downloaded files below. Please also change paths with file names downloaded in scripts. 
These can be downloaded from the GENCODE database
a. GRCh38 (hg38 human genome): https://www.gencodegenes.org/human/release_39.html

i. Download the comprehensive gene annotation (CHR) GTF and GFF3 files

ii. Genome sequence, primary assembly (GRCh38)	(PRI) fasta file

b. Repeats: https://genome.ucsc.edu/cgi-bin/hgTables for hg38

Choose Group: Repeats, Track: RepeatMasker, Table:rmsk

Region: Genome

Add output name and output field separator csv(for excel)

2. Bowtie2 rRNA index

This can be created one time for human genome using:

```bash
mkdir rrna-files
esearch -db nucleotide -query "NR_003285.3 OR NR_003286.4 OR NR_003287.4 OR NR_023363.1" | efetch -format fasta > human_rRNA_refs.fasta
bowtie2-build human_rRNA_refs.fasta human_rRNA_index
```

3. STAR genome index 

GRCh38 files downloaded in 1. above should be available on the system as placed in the directory as instructed.  (Any other genome versions can be downloaded from gencode/ USCS browser ina  similar way)

The paths should be configured for the reference genome being analyzed. The paths can be edited in script align-dedup.sh

4. Input data files

a. Create a base directory for the project.

b. Create a folder inside eg. fastq-files that stores all fastq files. Place fastq files in .gz format

c. create list-fastq.txt and sample-id.txt files. You can follow example files in preprocess-and-align and peak-calling-macs3 folder.

d. Download all genome files and create index files as instructed above. Change paths and create directories as shown for GRCh38 above if using another version of genome.

