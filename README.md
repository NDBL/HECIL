# HECIL
Hybrid Error Correction of Long Reads using Iterative Learning

Dependencies:

BWA (tested on version 0.7.15)
Samtools (tested on version 0.1.18)
Python (tested on version 2.7.8)

Running HECIL:
usage: python HECIL.py [-h] -l LR -s SR -len LENGTH_SR -o OUTPUT [-lc LOW_CONFIDENCE]
                [-c CUTOFF_QC] [-k CONFIDENCE_THRESHOLD]
                
Example command: python HECIL.py -l LongRead.fa -s ShortRead.fq -len 202 -o Out



Contact: ochoudhu@nd.edu or semrich@nd.edu
