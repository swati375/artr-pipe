## Process narrowPeak output files 

The program creats granges object to easily process and store for R. It also does basic processing, assigns gene names, go analysis, and peak regions etc.
The program can also create fasta files for all peaks or peaks from selected region eg. all 3UTR peaks. An example to call the program is shown below:

You can call the program in two ways. 

A. The first one below can directly be used to run the entire script.
It takes as input output from macs3 peak calling. Details about required inputs can be seen using:

```bash
conda activate rbase_43
Rscript peaks-analysis/cmd_peak_analysis.R -h
```

```bash
Rscript peaks-analysis/cmd_peak_analysis.R -p "/home/user/Desktop/data/ptbp1/peak-files-se-ext30" -f "peak.WT_PTBP1.narrowPeak.bed" -P "PTBP1"
```

B. The second way is to run command wise in Rstudio. You can use peak-analysis/Rstudio_peak_analysis.R Just open the script in Rstudio. Make sure to activate the rbase_43 environment.

Both of the scripts do the same job but one can be run directly on command line with providing all required inputs, while the other can prompt an input with a separate command.


## stats_scripts

Other scripts for performing basic statistics on peak files, motifs and for plots are also in the stats-scripts folder.

### Differential analysis
For comparing peaks between conditions, it is first required to find consensus locations on the genome and count reads for each consensus location across all conditions being compared. This is done by deg-consensus-generation.py script

Use a conda environment
```bash
conda activate seq-py312
python stats-scripts/deg-consensusfiles-generation.py
```
Follow the instructions on the screen and enter input files as required (eg. the sample-id.txt file we created earlier), conditions to be compared etc.

Then, in Rstudio you can use the 'Rstudio-deg-conditionwise.R' script by processing it line by line. It reads the consensus files and does DEG analysis on the files generated above. In the same script you can also plot scatter plots. 

To derive fasta sequences for selected peaks for motif analysis, in Rstudio you can use 'deg_motifanalysis.R'


### Motif Search

A. Using MEME-CHIP suite. You can go to website https://meme-suite.org/meme/doc/meme-chip.html

For motif discovery, select Motif Discovery-> XTREME from the left menu. Fill up form and upload the generated fasta sequence file for selected sequences. It will generate a zip file of discovered motifs. You can also get it emailed (fill details before submitting job)

B. using homer

We already create the fasta files from sequences after peak calling R scripts, based on region of peaks and again in deg R scripts for selected 
Use a conda environment of this. 
```bash
conda activate motif
bash stats-scripts/getmotif_fromfastaseq.sh
```

Now, follow the instructions on the screen. the programs asks if you want to do motif discoovery or search, what motif you want to match it to etc. For motif match, example files to use as input are in the stats-scripts folder (nsun2.motif, nsun6.motif)


### read and peak statistics

A. peakfile_stats.py: The program calculates and plots average peak length for all peaks for an input narrowPeak macs3 output file

B. bam_sizestats.py: This program plots histogram for tlen/ insert size from bam file

C. filter_peaks.py: This program filters macs3 called peaks by percentile. It removes 15th percentile (can be chnaged) using values in column 7 of narrowpeaks file. It also draws a distribution to show which peaks were removed.

D. dist_motif2peaksummit.R: This program finds distance of given motifs from peak summit. It can find both all motif occurrences and the closest motif occurrence to peak summit and plots the same.

E. get_motif.bash: This program finds motifs using homer. It directly takes the narrowpeak file and extracts +-20 nt from peak start and end to run homer.