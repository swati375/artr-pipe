import pysam
import matplotlib.pyplot as plt
import os,sys

#### plot histogram for tlen/ insert size from bam file

def hist_bam_tlen_readlength(genome_bam,dir_bam):
    fig,(ax0,ax1) = plt.subplots(nrows=2,ncols=1)
    insert_sizes = []
    with pysam.AlignmentFile(genome_bam, "rb") as bam:
        for read in bam:
            if read.is_proper_pair and read.is_read1 and read.template_length > 0:
                insert_sizes.append(read.template_length)

    # plt.figure()
    bins=[0,100,200,300,400,500,600,1000]
    counts,bins,bars=ax0.hist(insert_sizes, bins=bins)
    ax0.set_xticks(bins)
    ax0.set_title("Insert Size Distribution (Genome BAM)")
    ax0.set_xlabel("Insert Size (bp)")
    ax0.set_ylabel("Frequency")
    print("Fragments:", len(insert_sizes))
    print("Mean:", sum(insert_sizes)/len(insert_sizes))
    print(counts,bins)

    read_lengths = []
    with pysam.AlignmentFile(genome_bam, "rb") as bam:
        for read in bam:
            if not read.is_unmapped:
                read_lengths.append(read.query_length)

    # plt.figure()
    bins=[0,30,50,70,100,250,500]
    counts, bins, bars = ax1.hist(read_lengths, bins=bins)
    ax1.set_xticks(bins)
    ax1.set_title("Read Length Distribution")
    ax1.set_xlabel("Read Length (bp)")
    ax1.set_ylabel("Frequency")

    print("Reads:", len(read_lengths))
    print("Mean read length:", sum(read_lengths)/len(read_lengths))
    print(counts,bins)
    fig.tight_layout()
    fig.savefig(dir_bam+'/bam_stats.png')
    print('plot saved in bam file directory')
    plt.show()


#### extract reads mapping to introns- have N in cigar

##In the case of STAR, any gap less than alignIntronMin (21 bases last I looked) is considered a deletion.

def intronicreads_frombam(genome_bam,outfile):
    command='samtools view -h '+genome_bam+" |awk '$1 ~ /^@/ || $6 ~ /N/' | samtools view -b >"+bam_dir+outfile
    os.system(command)
    print('\n Intronic reads filtered')



if __name__=="__main__":
    dir_bam=input('\nEnter directory with bam file:')
    #/home/swati/Desktop/ARTR-seq/ezgi-data/2025-09/staralign-bam-files/48184_ET6_E6442_G8/dedup.bam
    bam_file=input('\nEnter bam file for analysis:')
    genome_bam=dir_bam+'/'+bam_file
    print(genome_bam)
    hist_bam_tlen_readlength(genome_bam,dir_bam)

    filter_intronic = input(
    '\nDo you want to filter intronic reads? (y/n): ').strip().lower()

    if filter_intronic in ('y', 'yes'):
        outputfile = input('\nEnter output BAM filename: ').strip()

        # Add directory if only a filename was provided
        outfile_path = dir_bam + '/' + outputfile

        intronicreads_frombam(genome_bam, outfile_path)