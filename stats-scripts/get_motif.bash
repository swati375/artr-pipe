#!/usr/bin/env bash

date
echo "Getting the motif using Homer2"

#####
# dir="/home/swati/Desktop/ARTR-seq/artr-origpaper-PTBP1/peak-files"
#oripeak="/home/swati/Desktop/ARTR-seq/artr-origpaper-PTBP1/peak-files/peak.PTBP1.narrowPeak.bed"
read -p "enter peaks file with complete path: " oripeak 
printf "Your narrowpeak.bed file -- %s\n" $oripeak

read -p "enter the directory for homer output with full path. " motifdir
mkdir $motifdir

 printf "Your directory for the homer2 output files" $motifdir


# set this variable if your data is not from human
genomeHomer2="hg38"
printf "Your homer2 genome tag -- %s\n" ${genomeHomer2}

############################################################
annoforHomer2="/home/swati/Desktop/ARTR-seq/common_genome_files/GRCh38/file-homer/gencode.v39.annotation-bghomer.bed12"
printf "Your annotation bed12 file for Homer2 %s \n" ${annoforHomer2}

nc=2
printf "Core number for homer2 -- %s\n" ${nc}
######

extpeak="$(awk -F/ '{print $NF}' <<< ${oripeak})"
odir=${extpeak}

extpeak="${extpeak/peak./ext.}"


cd $motifdir
pwd

if [ ! -d "${odir}" ]; then
	echo "NO, mkdir" ${odir}
	mkdir ${odir};
else
	echo "YES"
fi

awk -v wid=20 'BEGIN{FS="\t"; OFS="\t"}; {$2=$2-wid; $3=$3+wid; print $0}' $oripeak > $motifdir/$extpeak
odir=${odir/.narrowpeak*/}

findMotifsGenome.pl $extpeak ${genomeHomer2} $odir/ -p ${nc} -rna -S 10 -len 4,5,6 \
	-bg ${annoforHomer2}

# findMotifsGenome.pl $extpeak ${genomeHomer2} $odir/ -p ${nc} -rna -S 10 -len 5,6,7,8,9 \
# 	-bg ${annoforHomer2}
#####
date
echo "Finish!"
