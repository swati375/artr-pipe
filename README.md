# artr-seq-processing
A command line workflow to process Paired and single end ARTR-seq reads from fastq files to peak calling and downstream analysis, including:

- FASTQ quality control
- Adapter trimming and UMI extraction
- rRNA filtering
- Genome alignment
- UMI-based deduplication
- Strand-specific BAM generation
- BigWig coverage generation
- MACS3 peak calling
- Basic downstream statistics

## Installation

```bash
git clone <REPOSITORY_URL>
cd artr-seq-processing
conda env create -f environment.yml
conda activate seq-py312
```
