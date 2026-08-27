#!/usr/bin/env bash

date
echo "Step 0. check fastq file quality"


echo -n "Please enter base directory: "
# eg. /home/desktop/seqdir/projA
read base_dir
echo "All output files will be created inside this folder or the input files should be inside this folder"

echo -n "Please enter file with fastq files for quality testing: \n The file should be placed inside the base dirextory provided before"
# eg. list-fastq.txt
read file
file_fastq=$base_dir/$file

echo -n "Please enter output directory name for fastqc results: \n"
read out_dir
echo $base_dir/$out_dir

echo -n "Please enter directory of fastq files: "
#eg. fastqc
read dir
fastq_dir=$base_dir/$dir


if [ ! -d $base_dir/$out_dir ];then
  # echo "Output directory doesnt exist, making" $out_dir
  mkdir $base_dir/$out_dir
  else
    echo "YES exists" $base_dir/$out_dir
  fi

# Run FastQC on all files from the list
dir_qc="$base_dir/$out_dir"
mkdir -p "$dir_qc"

fastqc -t 8 -o "$dir_qc" $(sed "s|^|$fastq_dir/|" "$file_fastq")

# Aggregate with MultiQC
multiqc -o "$dir_qc" "$dir_qc"
