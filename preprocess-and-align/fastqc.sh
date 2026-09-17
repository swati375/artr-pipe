#!/usr/bin/env bash

date
echo "Step 0. check fastq file quality"


echo -n "Please enter base directory (Please provide fill path): "
# eg. /home/desktop/seqdir/projA
read base_dir
echo "All input files should be placed inside and corresponding folders for output will be created inside this folder"

echo -n "Please enter file with fastq files for quality testing (eg. list-fastq.txt; see example file): "
# eg. list-fastq.txt
read file
file_fastq=$base_dir/$file

echo -n "Please enter output directory name for fastqc results (eg. fastqc; it will be created inside the base directory):"
read out_dir
echo $base_dir/$out_dir

echo -n "Please enter directory of fastq files (eg. fastq-files; all project fastq files should be placed inside the folder): "
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
