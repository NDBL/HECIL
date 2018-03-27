# HECIL
Hybrid Error Correction of Long Reads using Iterative Learning

Dependencies:

BWA (tested on version 0.7.15),
Samtools (tested on version 0.1.18),
Python (tested on version 2.7.8)

Running HECIL:

usage: python HECIL.py [-h] -l LR -s SR -len LENGTH_SR -o OUTPUT [-lc LOW_CONFIDENCE]
                [-c CUTOFF_QC] [-k CONFIDENCE_THRESHOLD]
                
Example command: python HECIL.py -l LongRead.fa -s ShortRead.fq -len 202 -o Out

The output containing all (corrected and uncorrected) long reads will be saved in Corrected_LongRead.fa


Troubleshooting:

1. If Align_Corr.sh, bwa, and samtools cannot execute after cloning the repository, run the following commands:
    chmod +x bwa
    chmod +x samtools
    chmod +x Align_Corr.sh




Contact: olivia.choudhury1@ibm.com
