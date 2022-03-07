#!/bin/sh

date

module purge
module load GCC/8.2.0-2.31.1 
module load Bowtie2/2.3.5.1-GCC-8.2.0-2.31.1

mkdir consensus_0
cd consensus_0

echo Col0
cat ../reference/Araport.fasta |\
 awk '{if (substr($1,1,1)==">") print $1 "_Col0"; else print $1;}' \
 > Col0.fasta

bowtie2-build Col0.fasta Col0

echo Tsu1
cat ../reference/Araport.fasta |\
   awk '{if (substr($1,1,1)==">") print $1 "_Tsu1"; else print $1;}' \
   > Tsu1.fasta

bowtie2-build Tsu1.fasta Tsu1

cd ..
date
