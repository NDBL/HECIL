
#!/bin/bash

# Align short and long reads
echo 'Started running alignment'
./bwa index $1 2>stdout.txt
./bwa mem -t 12 $1 $2 > Out.sam 2>stdout.txt

# Create Pileup File
./samtools view -bS Out.sam | ./samtools sort - Out.sort 2>stdout.txt
./samtools mpileup -s -f $1 Out.sort.bam  > Pileup.txt 2>stdout.txt

# Remove intermediate file
rm $1.sa
rm $1.pac
rm $1.bwt
rm $1.ann
rm $1.amb
rm Out.sort.bam

echo -e 'Finished running alignment\n'

#----------------------------------------------

# HECIL Core Algorithm

#Get no. of long reads
d=2
num_refline="$(wc -l "$1" | awk '{print $1}')"
num_ref=$(($num_refline / $d))

#'./Align.sh '+LR+' '+SR+' '+len_SR+' '+out+' '+lc+' '+c+' '+k
#$1=LR, $2=SR, $3=len_SR, $4=Output, $5=Uncorr (LowConf Corr), $6=cutoff_QC, $7=conf_threshold

echo 'Started running correction'
python Correction.py Pileup.txt $1 $5 $num_ref $3 $6 $7 > $4

echo -e 'Finished running correction\n'

num_correfline="$(wc -l "$4" | awk '{print $1}')"
num_corref=$(($num_correfline / $d))

echo -e 'Generating output file with all (corrected and uncorrected) reads\n'

python Create_Corrected_AllLRReads.py $1 $num_ref $4 $num_corref

rm $1.fai
rm Out.sam

echo -e "Finished running HECIL\n"

