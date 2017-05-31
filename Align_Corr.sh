
#!/bin/bash

# Align short and long reads
echo 'Alignment started'
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

echo 'Alignment finished'

#----------------------------------------------

# HECIL Core Algorithm

#Get no. of long reads
d=2
num_refline="$(wc -l "$1" | awk '{print $1}')"
num_ref=$(($num_refline / $d))

#'./Align.sh '+LR+' '+SR+' '+len_SR+' '+out+' '+lc+' '+c+' '+k
#$1=LR, $2=SR, $3=len_SR, $4=Output, $5=Uncorr (LowConf Corr), $6=cutoff_QC, $7=conf_threshold


echo 'Correction started'
python Correction.py Pileup.txt $1 $5 $num_ref $3 > $4

echo 'Correction finished'


