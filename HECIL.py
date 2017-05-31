
import time, argparse, os

# Parse command-line input
parser = argparse.ArgumentParser(description='HECIL : Hubrid Error Correction of Long Reads using Iterative Learning')
# Required arguments
requiredNamed = parser.add_argument_group('required arguments')
requiredNamed.add_argument('-l','--LR', help='Erroneous long reads',required=True)
requiredNamed.add_argument('-s','--SR', help='Accurate short reads',required=True)
requiredNamed.add_argument('-len','--length_SR', help='length of short reads', required=True)
requiredNamed.add_argument('-o','--output', help='Output file in FASTA format containing corrected long reads',required=True)
# Optional arguments
parser.add_argument('-lc','--low_confidence', help='Output file containing low-confidence corrections in pileup format', default='LowConf.txt')
parser.add_argument('-c','--cutoff_QC', help='Cutoff value (between 0 and 1) for quick correction', default='0.85')
parser.add_argument('-k','--confidence_threshold', help='Confidence threshold (between 0 ans 2) to categorize low and high confidence corrections', default='1.2')


args = vars(parser.parse_args())

LR=args['LR']
SR=args['SR']
len_SR=args['length_SR']
out=args['output']
lc=args['low_confidence']
c=args['cutoff_QC']
k=args['confidence_threshold']


cmd_align='./Align_Corr.sh '+LR+' '+SR+' '+len_SR+' '+out+' '+lc+' '+c+' '+k
os.system(cmd_align)



#cutoff_quickcorrect=0.85
#SR_readlength=202
#num_LRread=33360
#max_qualscore=60
#conf_threshold=1.2





