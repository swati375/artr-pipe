import pysam
import matplotlib.pyplot as plt
import os,sys

#### plot histogram for tlen/ insert size from bam file

def hist_bam_tlen(genome_bam):
    insert_sizes = []

    with pysam.AlignmentFile(genome_bam, "rb") as bam:
        for read in bam:
            if read.is_proper_pair and read.is_read1 and read.template_length > 0:
                insert_sizes.append(read.template_length)

    plt.figure()
    counts,bins,bars=plt.hist(insert_sizes, bins=[0,100,200,300,400,500,600,1000])
    plt.title("Insert Size Distribution (Genome BAM)")
    plt.xlabel("Insert Size (bp)")
    plt.ylabel("Frequency")
    plt.show()

    print("Fragments:", len(insert_sizes))
    print("Mean:", sum(insert_sizes)/len(insert_sizes))
    print(counts,bins)


def hist_bam_readlength(genome_bam):
    read_lengths = []

    with pysam.AlignmentFile(genome_bam, "rb") as bam:
        for read in bam:
            if not read.is_unmapped:
                read_lengths.append(read.query_length)

    plt.figure()
    counts, bins, bars = plt.hist(read_lengths, bins=[0,30,50,70,100,250,500])
    plt.title("Read Length Distribution")
    plt.xlabel("Read Length (bp)")
    plt.ylabel("Frequency")
    plt.show()

    print("Reads:", len(read_lengths))
    print("Mean read length:", sum(read_lengths)/len(read_lengths))
    print(counts,bins)

#### extract reads mapping to introns- have N in cigar

##In the case of STAR, any gap less than alignIntronMin (21 bases last I looked) is considered a deletion.

def intronicreads_frombam(genome_bam,bam_dir,outfile):
    command='samtools view -h '+genome_bam+" |awk '$1 ~ /^@/ || $6 ~ /N/' | samtools view -b >"+bam_dir+outfile
    os.system(command)


if __name__=="__main__":
    dir_bam=input('\nEnter bam dir for analysis:')
    # bam_file=input('\nEnter bam file:')
    bam_file='dedup.bam'#'Aligned.sortedByCoord.out.bam'
    genome_bam=dir_bam+bam_file
    #/scratch/swatig/ARTR-seq/2026-01-26-tia/staralign-bam-files/48966_ET1_WT_TIA1_E6440_A3/Aligned.sortedByCoord.out.bam
    hist_bam_readlength(genome_bam)
    # hist_bam_tlen(genome_bam)
    # # outfile=input('\nENter output file :')
    outfile='intronic_reads.bam'
    intronicreads_frombam(genome_bam,dir_bam,outfile)
    # hist_bam_tlen(dir_file+outfile)