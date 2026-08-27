import matplotlib.pyplot as plt
import os
import sys

def avg_peak_length(narrowpeak_file):
    lengths = []

    with open(narrowpeak_file) as f:
        for line in f:
            if line.strip() == "":
                continue
            fields = line.strip().split("\t")
            start = int(fields[1])
            end = int(fields[2])
            lengths.append(end - start)

    avg_len = sum(lengths) / len(lengths)

    print("Number of peaks:", len(lengths))
    print("Average peak length:", avg_len)

    return avg_len


def peak_length_stats(narrowpeak_file):
    lengths = []

    with open(narrowpeak_file) as f:
        for line in f:
            if line.strip() == "":
                continue
            fields = line.strip().split("\t")
            start = int(fields[1])
            end = int(fields[2])
            lengths.append(end - start)

    plt.figure()
    counts,bins,bars=plt.hist(lengths, bins=[0,50,100,150,200,300,400,500])
    plt.title("Peak Length Distribution")
    plt.xlabel("Peak Length (bp)")
    plt.ylabel("Frequency")
    plt.show()

    print("Number of peaks:", len(lengths))
    print("Mean:", sum(lengths)/len(lengths))
    print("Min:", min(lengths))
    print("Max:", max(lengths))
    print(counts,bins)

if __name__=='__main__':
    dir_file=input('\nEnter directory where narrowpeak bed file output from macs3 is located:')
    file_list=[item for item in os.listdir(dir_file) if '.narrowPeak.bed' in item]
    for item in file_list:
        print(item)
        avg_peak_length(dir_file+item)
        peak_length_stats(dir_file+item)
