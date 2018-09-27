
#This program adds the uncorrected LR reads (not in Pileup file)
#to the output file of HECIL

import fileinput, os, sys, re

#USAGE: python Create_Corrected_AllLRReads.py Orig_LR.fa Num_OrigLR LR_Corrected.fa Num_CorrLR

ref_file=sys.argv[1]
num_LRread=int(sys.argv[2])
corrref_file=sys.argv[3]
num_corrLRread=int(sys.argv[4])

#Split filename and extension
col=re.split('.fa|.fasta ',ref_file)
out=col[0]+'_Corrected.fasta'

fout=open(out,'w')

# Store all pre-corrected reads
dict_ref={}
header=''
read=''
for line in fileinput.input(ref_file):
	if (line[0]=='>'):
		dict_ref[header]=read
		header=line.rstrip('\n')
		read=''

	else:
		read=read+line.rstrip('\n')

#For the last read
dict_ref[header]=read


dict_corrref={}
header=''
read=''

for line in fileinput.input(corrref_file):
        if (line[0]=='>'):
                dict_corrref[header]=read
                header=line.rstrip('\n')
                read=''

        else:
                read=read+line.rstrip('\n')

#For the last read
dict_corrref[header]=read



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
