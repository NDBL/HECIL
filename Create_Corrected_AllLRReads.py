
#This program adds the uncorrected LR reads (not in Pileup file)
#to the output file of HECIL

import fileinput, os, sys

#USAGE: python Create_Corrected_AllLRReads.py Orig_LR.fa Num_OrigLR Corrected_LR.fa Num_CorrLR

ref_file=sys.argv[1]
num_LRread=int(sys.argv[2])
corrref_file=sys.argv[3]
num_corrLRread=int(sys.argv[4])

out=ref_file+".corrected"
fout=open(out,'w')
fref=open(ref_file,'r')

dict_ref={}

# Store all pre-corrected reads
for c in range(int(num_LRread)):
	# Check for header
	l1=fref.readline()
	l2=fref.readline()
	l1=l1.rstrip('\n')
	l2=l2.rstrip('\n')

	dict_ref[l1]=l2


fcorrref=open(corrref_file,'r')

dict_corrref={}

# Store all corrected reads (subset of original list)
for c1 in range(int(num_corrLRread)):
	lc1=fcorrref.readline()
	lc2=fcorrref.readline()
	lc1=lc1.rstrip('\n')
	lc2=lc2.rstrip('\n')

	dict_corrref[lc1]=lc2

for ref in dict_ref.keys():
	if (len(ref)>2):
		if ref in dict_corrref.keys():
			header=ref+'\n'
			read=dict_corrref[ref]+'\n'
		else:
			header=ref+'\n'
			read=dict_ref[ref]+'\n'

		fout.write(header)
		fout.write(read)


fout.close()
