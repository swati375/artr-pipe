import os
import sys

## modified for paired end files by using -f BAMPE parameter

print("Step 6. calling peaks with macs3\n")
base_dir=input('Enter  Base directory Path:')
peakdir=base_dir+'/peak-files-default/'#/peak-files_gap30_len100/'
command='mkdir '+peakdir
os.system(command)
# print("Your directory for the output peak file "+peakdir)

#change genomeForMacs3 if data is not from human. Follow the macs3 manual
genomeForMacs3="hs"
stardir=base_dir+"/staralign-bam-files/"

print(' Create input file with bam files for peak calling. Enter text file containg protein IP files as mentioned in readme, with control labelled as Input')
filein=input('\nEnter file with Path:')#eg. sample-id.txt
prot_dict={}
with open(filein, 'r') as file:
	for line in file:
		line = line.rstrip()
		var = line.split(':')
		if ',' in var[1]:
			varval = var[1].split(',')
			prot_dict[var[0]] = varval
		else:
			prot_dict[var[0]] = [var[1]]
print(prot_dict)

print("You have following proteins :"+ str(prot_dict.keys()))


ans=input("Run for all proteins? [y/n]:")
if ans=='y' or ans=='Y':
	prot_list=list(prot_dict.keys())
	prot_list=[item for item in prot_list if 'Input' not in item]
	if 'Input' in prot_list:
		prot_list.remove('Input')
else:
	prot=input("\nEnter protein of choice:")
	prot_list=[prot]

for prot in prot_list:
	if prot in prot_dict.keys():
		print(prot_dict[prot])
		fwdIPbams=[stardir+str(prot_dict[prot][i])+'/dedup.fwd.bam' for i in range(len(prot_dict[prot]))]
		revIPbams=[stardir+str(prot_dict[prot][i])+'/dedup.rev.bam' for i in range(len(prot_dict[prot]))]
		print(fwdIPbams)

	outdir=peakdir+prot
	print("Finding peaks for "+prot+" using:\n Protein IP id-- "+str(prot_dict[prot])+'\n')
	## forward
	# by default q<0.05 used
	##max-gap by default is read length
	##min-length is calculated d
	input_sample='Input_'+prot
	if input_sample in prot_dict.keys():
		print('Input exists forward')
		fwdCtrbams=[stardir+str(prot_dict[input_sample][i])+'/dedup.fwd.bam' for i in range(len(prot_dict[input_sample]))]
		command="macs3 callpeak --treatment "+" ".join(fwdIPbams)+" --control "+" ".join(fwdCtrbams)+" -f BAMPE -n "+prot+".fwd -g "+genomeForMacs3+ \
		" -B --keep-dup all -q 0.05 --outdir "+outdir+ " --tempdir "+outdir+" --nomodel --extsize 30"# --max-gap 30 --min-length 100"
	else:
		command="macs3 callpeak --treatment "+" ".join(fwdIPbams)+" -f BAMPE -n "+prot+".fwd -g "+genomeForMacs3+ \
		" -B --keep-dup all -q 0.05 --outdir "+outdir+ " --tempdir "+outdir+" --nomodel --extsize 30"# --max-gap 30 --min-length 100"

	os.system(command)

	##reverse
	if input_sample in prot_dict.keys():
		print('Input exists reverse')
		revCtrbams=[stardir+str(prot_dict[input_sample][i])+'/dedup.rev.bam' for i in range(len(prot_dict[input_sample]))]
		command="macs3 callpeak --treatment "+" ".join(revIPbams)+" --control "+" ".join(revCtrbams)+" -f BAMPE -n "+prot+".rev -g "+ \
				genomeForMacs3+" -B -q 0.05 --keep-dup all -q 0.05 --outdir "+outdir+ \
				" --tempdir "+outdir+" --nomodel --extsize 30"# --max-gap 30 --min-length 100"
	else:
		command="macs3 callpeak --treatment "+" ".join(revIPbams)+" -f BAMPE -n "+prot+".rev -g "+ \
				genomeForMacs3+" -B -q 0.05 --keep-dup all --outdir "+outdir+ \
				" --tempdir "+outdir+" --nomodel --extsize 30"# --max-gap 30 --min-length 100"

	os.system(command)

	outpeak=peakdir+'peak.'+prot+'.narrowPeak.bed'
	command= "cat "+outdir+"/*peaks.narrowPeak | awk 'BEGIN{FS=\"\\t\";OFS=\"\\t\"} {if ($4 ~ /fwd/) {$6 = \"+\"} "+ \
			"else {$6 = \"-\"}; print $0 }' > "+outpeak

	os.system(command)


