#!/usr/bin/env bash
date

echo "Step2. mapping reads to a genome, and deduplicate with umi-tools, and spliting bams by strand"

base_dir="$2"
pair="$1"
# # echo $base_dir
# # echo $pair

indir="$base_dir/trimmed-files"

stardir="$base_dir/staralign-bam-files"
mkdir $stardir

staridx=/scratch/swatig/ARTR-seq/artr-origpaper-YTH/index-files
mkdir $staridx
printf "Your STAR index directory -- %s\n" $staridx

## generate star index
##download refrence files gft annotation(CHR) and fasta files (PRI)
#index-files already available
# STAR --runThreadN 8 --runMode genomeGenerate --genomeDir $staridx \
# --genomeFastaFiles $base/GRCh38.primary_assembly.genome.fa \
# --sjdbGTFfile $base/gencode.v39.annotation.gtf --sjdbOverhang 100

nc=5
# printf "Your core number used for STAR -- %s\n" $nc
##

startt="$(date +%s)"

for fi in $(ls $indir/norRNA.*)
do
	
	filename=$(basename "$fi" )

	if [[ "$pair" =~ ^([Tt][Rr][Uu][Ee])$ && "$filename" == norRNA.2* ]]; then
	    echo "Skipping mate 2 file: $filename"
	    continue
	fi

	{
	echo "3. mapping reads to a genome with STAR"
	
	# echo "The file is $filename"
	startt1="$(date +%s)"
	inrawR1=$fi

	odir=${filename/norRNA./}
	odir=${stardir}/${odir}
   	# echo $odir

	if [ "$pair" == FALSE ]; then
		inf=${inrawR1}
	else
		inrawR2=${inrawR1/.1./.2.}
		inf="${inrawR1} ${inrawR2}"
		echo $inrawR1 $inrawR2
		odir=${odir/1./}
	fi
	# echo $odir

	if [ ! -d $odir ];then
		echo "NO, mkdir" $odir
		mkdir $odir
	else
		echo "YES"
	fi

	STAR --runMode alignReads --runThreadN ${nc} \
		--readFilesCommand zcat \
		--genomeDir $staridx \
		--alignEndsType EndToEnd \
		--genomeLoad NoSharedMemory \
		--quantMode TranscriptomeSAM \
		--alignMatesGapMax 15000 \
		--readFilesIn ${inf} \
		--outFileNamePrefix $odir/ \
		--outFilterMultimapNmax 1 \
		--outSAMattributes All \
		--outSAMtype BAM SortedByCoordinate \
		--outFilterType BySJout \
		--outReadsUnmapped Fastx \
		--outFilterScoreMin 10 \
		--outFilterMatchNmin 24

	samtools index -@ ${nc} $odir/Aligned.sortedByCoord.out.bam
	endt="$(date +%s)"
	printf "Mapping reads to a genome for %s done, elapsed time -- %.0f s\n\n" $inrawR1 "$((( $endt - $startt1)))"
	printf "\n"

	startt2="$(date +%s)"
	echo "4. deduplicating reads with umi-tools"
	
	cd $odir ## so temporary files are not mixed due to same input file names although in different folders
	umi_tools dedup --method unique\
		-I $odir/Aligned.sortedByCoord.out.bam \
		--output-stats=$odir/dedup \
		-L $odir/umitools.log \
		-S $odir/dedup.bam

	samtools index -@ ${nc} $odir/dedup.bam
	endt="$(date +%s)"
	printf "Deduplicating for %s done, elapsed time -- %.0f s\n\n" $filename "$((( $endt - $startt2)))"

	# ##
	startt3="$(date +%s)"
	echo "5. splitting reads by strand"
	ofwd=$odir/dedup.fwd.bam
	orev=${ofwd/.fwd.bam/.rev.bam}
	curinbam=$odir/dedup.bam
	echo ${curinbam}
	echo $ofwd

	samtools view -@ ${nc} -F 16 ${curinbam} -b -o $ofwd &&
		samtools index -@ ${nc} $ofwd;

	echo $orev
	samtools view -@ ${nc} -f 16 ${curinbam} -b -o $orev &&
		samtools index -@ ${nc} $orev
	endt="$(date +%s)"
	printf "Splitting reads for %s done, elapsed time -- %.0f s\n\n" $curinbam "$((( $endt - $startt3)))"
	} &
done
wait
##
date
endt="$(date +%s)"
printf "Total elapsed time -- %.2f min \n\n" "$((($endt - $startt) / 60))"
echo "Finish!"
######
















