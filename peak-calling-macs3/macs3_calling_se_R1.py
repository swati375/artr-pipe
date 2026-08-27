import os
import sys
## selecting only R1 reads from bam file for macs3 call
print("Step 6. calling peaks with macs3\n")
base_dir=input('Enter Path:')
peakdir=base_dir+'/peak-files-se-ext30-R1/'
command='mkdir '+peakdir
os.system(command)
# print("Your directory for the output peak file "+peakdir)

#change genomeForMacs3 if data is not from human. Follow the macs3 manual
genomeForMacs3="hs"
# print("Your genome setting for macs3 -- "+genomeForMacs3)
stardir=base_dir+"/staralign-bam-files/"

print(' Create input file with bam files for peak calling. Enter text file containg protein IP files as mentioned in readme, with control labelled as Input')
filein=input('\nEnter file with Path:')# /scratch/swatig/ARTR-seq/2025-09-03_AAGK7MVM5/Marijke_P._A._Baltissen/sample-id.txt
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
else:
	prot=input("\nEnter protein of choice:")
	prot_list=[prot]

for prot in prot_list:
	if prot in prot_dict.keys():
		print(prot_dict[prot])
		# if len(prot_dict[prot])>1:
		fwdIPbams=[stardir+str(prot_dict[prot][i])+'/R1-bam/dedup.fwd.bam' for i in range(len(prot_dict[prot]))]
		revIPbams=[stardir+str(prot_dict[prot][i])+'/R1-bam/dedup.rev.bam' for i in range(len(prot_dict[prot]))]
		# else:
			# fwdIPbams=[stardir+str(prot_dict[prot])+'/dedup.fwd.bam']
			# revIPbams=[stardir+str(prot_dict[prot])+'/dedup.rev.bam']
		print(fwdIPbams)

	outdir=peakdir+prot
	print("Finding peaks for "+prot+"using:\n Protein IP id-- "+str(prot_dict[prot])+'\n')
	## forward
	# by default q<0.05 used
	input_sample='Input_'+prot
	if 'Input' in prot_dict.keys():
		print('Input exists forward')
		fwdCtrbams=[stardir+str(prot_dict[input_sample][i])+'/R1-bam/dedup.fwd.bam' for i in range(len(prot_dict[input_sample]))]
		command="macs3 callpeak --treatment "+" ".join(fwdIPbams)+" --control "+" ".join(fwdCtrbams)+" -f BAM -n "+prot+".fwd -g "+genomeForMacs3+ \
		" -B --keep-dup all -q 0.05 --outdir "+outdir+ " --tempdir "+outdir+" --nomodel --extsize 30"
	else:
		command="macs3 callpeak --treatment "+" ".join(fwdIPbams)+" -f BAM -n "+prot+".fwd -g "+genomeForMacs3+ \
		" -B --keep-dup all -q 0.05 --outdir "+outdir+ " --tempdir "+outdir+" --nomodel --extsize 30"

	# print(command)
	os.system(command)

	##reverse
	if 'Input' in prot_dict.keys():
		print('Input exists reverse')
		revCtrbams=[stardir+str(prot_dict[input_sample][i])+'/R1-bam/dedup.rev.bam' for i in range(len(prot_dict[input_sample]))]
		command="macs3 callpeak --treatment "+" ".join(revIPbams)+" --control "+" ".join(revCtrbams)+" -f BAM -n "+prot+".rev -g "+ \
				genomeForMacs3+" -B --keep-dup all -q 0.05 --outdir "+outdir+ \
				" --tempdir "+outdir+" --nomodel --extsize 30"
	else:
		command="macs3 callpeak --treatment "+" ".join(revIPbams)+" -f BAM -n "+prot+".rev -g "+ \
				genomeForMacs3+" -B -q 0.05 --keep-dup all --outdir "+outdir+ \
				" --tempdir "+outdir+" --nomodel --extsize 30"

	os.system(command)

	outpeak=peakdir+'peak.'+prot+'.narrowPeak.bed'
	command= "cat "+outdir+"/*peaks.narrowPeak | awk 'BEGIN{FS=\"\\t\";OFS=\"\\t\"} {if ($4 ~ /fwd/) {$6 = \"+\"} "+ \
			"else {$6 = \"-\"}; print $0 }' > "+outpeak
	# print(command)
	os.system(command)

