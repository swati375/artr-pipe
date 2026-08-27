#!/usr/bin/env bash

date
echo "Step 1 trimming adapters"

###
##
# note of paired-end reads
# TRUE -- paired-end reads
# FALSE -- single-end read
echo -n "Please enter if fastq reads/ sequencing is pairend end (TRUE) or not (FALSE): "
read bool_var
pair=$bool_var

# Check that paired-end option is valid
if [[ "$pair" != "TRUE" && "$pair" != "FALSE" ]]; then
    echo "ERROR: Please enter TRUE for paired-end reads or FALSE for single-end reads."
    exit 1
fi

echo -n "Please enter base directory: "
# eg. /home/desktop/seqdir/projA
read base_dir

rawfiledir=$base_dir/fastq-files/
printf "Your directory containing the raw fq.gz -- %s\n" $rawfiledir

outfdir=$base_dir/trimmed-files
mkdir $outfdir
printf "Your directory containing the output files fq.gz -- %s\n" $outfdir

###
startt="$(date +%s)"

###
for fi in $(ls $rawfiledir/*.fastq.gz)
do
	
	(filename=$(basename "$fi" )
	echo "The file is $filename"
	inrawR1=$filename

	echo "1 Trimming adapters with cutadapt"
	if [ "${pair}" == FALSE ]; then
		echo "Single-end"
		inrawR2=""
		printf "Your Read 1 fq.gz -- %s, Read 2 fq.gz -- %s\n" $inrawR1 $inrawR2

		t1=${outfdir}/trim1.${inrawR1}
		# log1=${t1/R1/} #Uses bash string substitution to remove the substring "R1" from t1 and assigns the result to log1.
		log1=${t1/.fastq.gz/.log}
		# echo $log1
		
		t2=${t1/trim1/trim}
		log=${log1/trim1/trim}
		
		echo "1.1 trimming the adapter"
		# echo $t1 $log1
		
		# trimming the 3'-adapter
		## removed -h ${nc} caz parallel not supported in python2 -m 32 retains reads with length atleast 32 so that after barcode and
		# 4 nt removal still 20 nt remain. Matches authors command to retain 20 nt reads
		# after trimming
		cutadapt --nextseq-trim=20 --action=trim \
			-a AGATCGGAAGAGCACACGTCTGAACTCCAG -m 32\
			-o $t1 $rawfiledir/$inrawR1 >$log1

		echo "1.2 extracting the umi"
		# echo $t2 $log

		umi_tools extract --bc-pattern=NNNNNNNN \
	  		--stdin=$t1 --stdout=$t2 \
	  		--log=$log


	  	echo "1.3 removing last 4 nt from 3 prime end"

	  	t2_trim3=${outfdir}/trimmed.${inrawR1}
		log_trim3=${log1/trim1/trimmed}
			
		cutadapt -q 20 -m 20 --action=trim -u -4 -o $t2_trim3 $t2 > $log_trim3

	else
		if [[ $inrawR1 == *"R1"* ]]; then
			echo "in Paired-end loop"
			inrawR2=${inrawR1/R1/R2}
			
			# Check that the matching R2 file exists
			if [[ ! -f "$rawfiledir/$inrawR2" ]]; then
			    echo "ERROR: Matching R2 file not found for $inrawR1"
			    echo "Expected: $rawfiledir/$inrawR2"
			    exit 1
			fi
			printf "Your Read 1 fq.gz -- %s, Read 2 fq.gz -- %s\n" $inrawR1 $inrawR2
			
			t11=${outfdir}/trim1.${inrawR1}
			t12=${t11/R1/R2}
			
			# log1=${t11/-R1/}
			log1=${t11/R1.fastq.gz/.log}

			t21=${t11/trim1/trim}
			t22=${t12/trim1/trim}
			log=${log1/trim1/trim}
			
			# echo $inrawR1 $inrawR2

			echo "1.1 trimming the adapter"
			# # echo $t11 $t12 $log1

			cutadapt --nextseq-trim=20 -m 28 --action=trim \
				-a AGATCGGAAGAGCACACGTCTGAACTCCAG \
				-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAG \
				-o $t11 -p $t12 $rawfiledir/$inrawR1 $rawfiledir/$inrawR2 > $log1

		                
			echo "1.2 extracting the umi"
			# # echo $t21 $t22 $log

			umi_tools extract --bc-pattern2=NNNNNNNN \
		  		--stdin=$t11 --stdout=$t21 \
		  		--read2-in=$t12 --read2-out=$t22 \
		  		--log=$log

		  	t21_trim3=${outfdir}/trimmed.${inrawR1}
		  	t22_trim3=${outfdir}/trimmed.${inrawR2}
			log_trim3=${log1/trim1/trimmed}

			cutadapt -q 20 -m 20 --action=trim \
				-u -4 -U -4 --pair-filter=any \
				-o $t21_trim3 -p $t22_trim3 $t21 $t22 > $log_trim3
		fi
		
	fi
) &

done
wait

endt="$(date +%s)"
printf "Trimming adapters for %s done, elapsed time -- %.0f s\n\n" $inrawR1 "$((( $endt - $startt)))"

bash filter-rrna-new.sh $pair $base_dir
bash align-dedup-new.sh $pair $base_dir
#####


### check trimmed file pairs have same no. of reads (for paired end data)
## zcat file_R1.fastq.gz | echo $((`wc -l`/4))
# zcat file_R2.fastq.gz | echo $((`wc -l`/4))
