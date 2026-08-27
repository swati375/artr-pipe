#!/usr/bin/env bash

date
echo "Step 1.3 Removing reads mapping to rRNA"

base_dir="$2"
pair="$1"
if [[ "$pair" != "TRUE" && "$pair" != "FALSE" ]]; then
    echo "ERROR: First argument must be TRUE (paired-end) or FALSE (single-end)."
    exit 1
fi
outfdir=$base_dir/trimmed-files

nc=4

bt2idx=rrna-files/human_rRNA_index
#################################### we assume that index files already created once, so below lines are commented #############################
# inside the bt2idx folder
# esearch -db nucleotide -query "NR_003285.3 OR NR_003286.4 OR NR_003287.4 OR NR_023363.1" | efetch -format fasta > human_rRNA_refs.fasta
# bowtie2-build human_rRNA_refs.fasta human_rRNA_index

printf "Your bowtie2 index directory -- %s\n" $bt2idx

for fi in $(ls $outfdir/trimmed*.fastq.gz)
do
	(filename=$(basename "$fi" | sed 's/^trimmed\.//')
	echo "The file is $filename"
	inrawR1=$filename
	echo "Discarding reads mapped to rRNA with bowtie2"

	startt2="$(date +%s)"
	if [ "$pair" = FALSE ]; then
		echo "Single-end"
		#inrawR2=""
		#printf "Your Read 1 fq.gz -- %s, Read 2 fq.gz -- %s\n" $inrawR1 $inrawR2
		printf "Your Read fq.gz -- %s\n" $inrawR1
		t2=${outfdir}/trimmed.${inrawR1}

		bw2o=${t2/trimmed./norRNA.}
		filetag=$t2
		
		# echo $t2 $bw2o
		printf "Your input files fq.gz -- %s\n" $t2
		bowtie2 --threads ${nc} --seedlen=15 -x $bt2idx \
			-U $t2 --un-gz $bw2o > /dev/null
		

	else
		if [[ $inrawR1 == *"R1"* ]]; then
	    		echo "Paired-end";
	    		inrawR2=${inrawR1/R1/R2}
	    		if [[ ! -f "${outfdir}/trimmed.${inrawR2}" ]]; then
			    echo "ERROR: Matching R2 file not found for $inrawR1"
			    echo "Expected: ${outfdir}/trimmed.${inrawR2}"
			    exit 1
			fi
	    		printf "Your Read 1 fq.gz -- %s, Read 2 fq.gz -- %s\n" $inrawR1 $inrawR2
			t21_trim3=${outfdir}/trimmed.${inrawR1}
			t22_trim3=${outfdir}/trimmed.${inrawR2}

			r1base=$(basename "$t21_trim3")
		    	core=${r1base#trimmed.}
		    	core=${core%_R1.fastq.gz}

		    	bw2o="${outfdir}/norRNA.${core}"

		    	echo "Bowtie2 un-conc prefix: $bw2o"

			filetag=$t21_trim3

			# echo $bw2o
			bowtie2 --threads ${nc} --seedlen=15 -x $bt2idx \
				-1 $t21_trim3 -2 $t22_trim3 --un-conc-gz $bw2o > /dev/null
		fi

	fi
	endt="$(date +%s)"
	printf "Filtering rRNA for %s done, elapsed time -- %.2f min \n\n" $filetag "$((($endt - $startt2) / 60))"
	
	) &
done
wait

