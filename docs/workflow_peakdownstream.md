```bash
Rscript peaks-analysis/cmd_peak_analysis.R -p "/home/user/Desktop/data/ptbp1/peak-files-se-ext30" -f "peak.WT_PTBP1.narrowPeak.bed" -P "PTBP1"
```
It takes as input output from macs3 peak calling. Details about required inputs can be seen using:

```bash
Rscript peaks-analysis/cmd_peak_analysis.R -h
```



Motif Search

A. Using MEME-CHIP suite. You can go to website https://meme-suite.org/meme/doc/meme-chip.html

For motif dicovery, sleect Motif Discovery-> XTREME from the left menu. Fill up form and upload the generated fasta sequence file for selected sequences. It will generate a zip file of dicovered motifs. You can also get it emailed (fill details before submitting job)

B. using homer

Use a conda environment of this. 
```bash
conda activate environment_motif
```