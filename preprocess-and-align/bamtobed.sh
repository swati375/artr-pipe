#!/bin/bash
# Move to the parent directory

## in conda env chipseq_env

echo -n "Please enter base directory: "
read base_dir
echo $base_dir

bam_dir="$base_dir/staralign-bam-files"
echo $bam_dir
crosslink_dir="$base_dir/crosslink_nt"
mkdir $crosslink_dir

cd $bam_dir

# ## for piranha
# for bam in */Aligned.sorted*.bam; do
#     # Extract a clean sample name (folder name)
#     (sample=$(basename "$(dirname "$bam")")
#     echo "Processing $sample ..."
#     echo $bam
#     $out=${bam/.out.bam/.bed}
#     bedtools bamtobed -bedpe -i "$bam" > $out

#     $out2=${out/.sortedByCoord.bed/.sorted.bed}
#     sort -k1,1 -k2,2n $out > $out2

# done


## CPM coverage files for IGV
for bam in */dedup.bam; do
  sample=$(basename "$(dirname "$bam")")
  echo "Processing $sample ..."
  echo $bam
  bamCoverage -b "$bam" -o "$crosslink_dir/${sample}_CPM.bw" \
  --normalizeUsing CPM
done


 
