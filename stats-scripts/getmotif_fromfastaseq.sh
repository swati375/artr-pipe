#!/usr/bin/env bash

## activate env motif

date
echo "Running HOMER motif analysis with FASTA file of selected peaks"

#####

# Input FASTA (peaks)
read -p "Enter FASTA (.fa) file (peaks): " fasta_file
printf "Peaks FASTA file-- %s\n" "$fasta_file"

# Background FASTA (optional)
read -p "Enter background FASTA (.fa) file (press Enter to skip): " bg_file

if [ -n "$bg_file" ]; then
    printf "Background FASTA -- %s\n" "$bg_file"
else
    echo "No background FASTA provided. Running without background."
fi

# Number of cores
nc=2
printf "Cores -- %s\n" "$nc"

#####

# Extract filename for subfolder
fname="$(basename "$fasta_file")"
odir="${fname%.fa}"

printf "Output subdirectory -- %s\n" "$odir"
printf "FASTA file -- %s\n" "$fasta_file"

#####

# Output directory
read -p "Enter output directory for HOMER results with path: " motifdir
mkdir -p "$motifdir"
printf "Output directory -- %s\n" "$motifdir"

cd "$motifdir" || exit 1
pwd

# Create output subdirectory
if [ ! -d "$odir" ]; then
    echo "Creating directory: $odir"
    mkdir -p "$odir"
else
    echo "Directory already exists: $odir"
fi

#####

# Ask user whether to perform motif discovery or defined motif search

echo ""
echo "Choose analysis type:"
echo "1. Motif discovery"
echo "2. Search for defined motif"
read -p "Enter choice (1 or 2): " analysis_choice

#####

if [ "$analysis_choice" -eq 1 ]; then

    echo ""
    echo "Running HOMER motif discovery"

    if [ -n "$bg_file" ]; then

        echo "Using background FASTA..."

        findMotifs.pl "$fasta_file" fasta "$odir/" \
            -bg "$bg_file" \
            -p "$nc" -rna \
            -len 4,5,6 \
            -S 10

    else

        echo "Running without background FASTA..."

        findMotifs.pl "$fasta_file" fasta "$odir/" \
            -p "$nc" -rna \
            -len 4,5,6 \
            -S 10

    fi

#####

elif [ "$analysis_choice" -eq 2 ]; then

    echo ""
    echo "Searching for defined motif(s)..."

    # Output motif file
    read -p "Enter motif file to match with path: eg at stats-scripts/nsun2.motif" motif_file


    #####

    # Define output filename

    fasta_base="${fname%.fa}"
    read -p "Enter suffix for output file with target motif search hits with path" suffix

    if [ -n "$bg_file" ]; then

        hitsfilename="${fasta_base}_${suffix}_hits_withbg.txt"

        echo "Searching defined motifs with background FASTA..."

        findMotifs.pl "$fasta_file" fasta "$odir/" \
            -find "$motif_file" \
            -bg "$bg_file" \
            -rna \
            > "$hitsfilename"

    else

        motiffilename="${fasta_base}_${suffix}_hits.txt"

        echo "Searching defined motifs without background FASTA..."

        findMotifs.pl "$fasta_file" fasta "$odir/" \
            -find "$motif_file" \
            -rna \
            > "$motiffilename"

    fi

    echo "Motif search results saved to: $motiffilename"

#####

else

    echo "Invalid choice. Please enter 1 or 2."
    exit 1

fi

#####

date
echo "Finish!"

#/home/swati/Desktop/ARTR-seq/scripts_server/edited/stats-scripts/nsun2.motif
#/home/swati/Desktop/ARTR-seq/ezgi-data/2026-08-tia/consensus_water/nsun2ko_lost_water.fa