import os
import sys
import itertools
import pandas as pd

def consensus_peaks(base_dir,sample,dir_name):

	## find a union of peaks across all datasets
	print("Step. build consensus peak set for selected samples\n")
	peakdir=base_dir+'/peak-files-se-ext30/'
	consensus_dir=base_dir+'/consensus_'+dir_name+'/'
	os.system('mkdir '+consensus_dir)
	os.system('mkdir '+consensus_dir+'consensus_peak/')
	file_list=[peakdir+'peak.'+item[0]+'.narrowPeak.bed' for item in sample.values()]
	# print(file_list)
	nm=" ".join(file_list)
	print(nm)

	## with strand info
	command='cat '+nm+' | cut -f1-6 |sort -k1,1 -k2,2n > '+\
			consensus_dir+'consensus_peak/allpeaks-withstrand.bed'
	os.system(command)
	command='bedtools merge -s -c 4,6 -o collapse,distinct -delim "|" -i '+\
			consensus_dir+'consensus_peak/allpeaks-withstrand.bed > '+\
			consensus_dir+'consensus_peak/consensus-withstrand_and_sample.bed'
	os.system(command)
	
##########################################################################################################################

###counting reads in each peak
def peak_read_counts(base_dir,samples,dir_name):
	print("Step. counting reads\n")
	consensus_dir=base_dir+'/consensus_'+dir_name+'/'
	consensus_count_dir=consensus_dir+'count_reads/'
	os.system('mkdir '+consensus_count_dir)
	staralign_dir=base_dir+'/staralign-bam-files/'
	# dir_list=os.listdir(staralign_dir)
	dir_list=samples.keys()
	for ele in dir_list:
		print(ele)
		command='bedtools coverage -a '+consensus_dir+'consensus_peak/consensus-withstrand_and_sample.bed -b '+\
				staralign_dir+ele+'/dedup.bam > '+consensus_count_dir+samples[ele][0]+'.counts'
		print(command)
		os.system(command)

		command='cut -f6 '+consensus_count_dir+samples[ele][0]+'.counts > '+\
				consensus_count_dir+samples[ele][0]+'.txt'
		os.system(command)

def read_sample_ids(filename):
    samples = {}

    with open(filename) as f:
        for line in f:
            line = line.strip()

            if not line or ":" not in line:
                continue

            sample_name, sample_id = line.split(":", 1)
            samples[sample_name] = sample_id

    return samples


if __name__=='__main__':
	base_dir=input('Enter Path of base directory:')
	samples={}
	print('Enter text file with sample ids created for peak calling before')
	filein=input('\nEnter file with Path:')# /scratch/swatig/ARTR-seq/2025-09-03_AAGK7MVM5/Marijke_P._A._Baltissen/sample-id.txt
	samples=read_sample_ids(filein)
	# Only show non-input samples
	sample_options = [
        name for name in samples
        if not name.startswith("Input_")  ]
	
	print("Available samples:")
	for i, sample_name in enumerate(sample_options, start=1):
		print(f"  {i}. {sample_name}")

	selection = input("\nSelect samples (e.g. 1,3): ").strip()

	selected_samples = [ sample_options[int(i.strip()) - 1] for i in selection.split(",")]

	selected_dict = {}
	selected_noinput_dict={}
	for sample_name in selected_samples:
		# Add selected sample
		selected_dict[samples[sample_name]] = [sample_name]
		selected_noinput_dict[samples[sample_name]] = [sample_name]
		# Find corresponding Input sample
		input_name = f"Input_{sample_name}"
		if input_name not in samples:
			raise ValueError(
				f"Input sample not found for {sample_name}: "
				f"expected '{input_name}'"
				)
		selected_dict[samples[input_name]] = [input_name]

	print("\nSelected samples:")
	print(selected_dict)
	print(selected_noinput_dict)

	dir_name=input('\nEnter name for this selection. It will be used to create a directory. Dont use spaces but "_" between words :')# /scratch/swatig/ARTR-seq/2025-09-03_AAGK7MVM5/Marijke_P._A._Baltissen/sample-id.txt
	consensus_peaks(base_dir,selected_noinput_dict,dir_name)
	peak_read_counts(base_dir,selected_dict,dir_name)
	
