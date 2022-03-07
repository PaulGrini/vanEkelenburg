#!/bin/sh

# Command line parameter 1 gives the number for this round.
# Give a one for the initial round where we mapped to Araport11 i.e consensus_0.
# This will build a pilon consensus for Col0 and another for Tsu1.

if [ $# -eq 0 ]; then
    echo "Please provide the number of this round."
    exit 1
fi
ROUND=$1
let PREV=${ROUND}-1
WORK_DIR=consensus_${ROUND}
PREV_REF=consensus_${PREV}
PREV_MAP=map_to_consensus_${PREV}
echo WORK_DIR $WORK_DIR
echo PREV_REF $PREV_REF
echo PREV_MAP $PREV_MAP
echo

date

echo MODULE
module purge
module load GCC/8.2.0-2.31.1 
module load SAMtools/1.9-GCC-8.2.0-2.31.1   ## saga
module load Pilon/1.23-Java-11  ## saga
module list
echo

mkdir ${WORK_DIR}
cd ${WORK_DIR}
pwd
echo

echo "Col0 Early + Late"
ln -s ../${PREV_REF}/Col0.fasta prev.Col0.fasta
samtools view -b -H ../${PREV_MAP}/01-Col-0-Early-BR1.Reference_1/Sorted.bam > Col0.header.bam  
echo MERGE
samtools merge -h Col0.header.bam -@ 4 -O BAM Col0.bam \
../${PREV_MAP}/01-Col-0-Early-BR1.Reference_1/Sorted.bam \
../${PREV_MAP}/01-Col-0-Early-BR1.Reference_2/Sorted.bam \
../${PREV_MAP}/02-Col-0-Early-BR2.Reference_1/Sorted.bam \
../${PREV_MAP}/02-Col-0-Early-BR2.Reference_2/Sorted.bam \
../${PREV_MAP}/03-Col-0-Early-BR3.Reference_1/Sorted.bam \
../${PREV_MAP}/03-Col-0-Early-BR3.Reference_2/Sorted.bam \
../${PREV_MAP}/07-Col-0-Late-BR1_.Reference_1/Sorted.bam \
../${PREV_MAP}/07-Col-0-Late-BR1_.Reference_2/Sorted.bam \
../${PREV_MAP}/08-Col-0-Late-BR2_.Reference_1/Sorted.bam \
../${PREV_MAP}/08-Col-0-Late-BR2_.Reference_2/Sorted.bam \
../${PREV_MAP}/09-Col-0-Late-BR3_.Reference_1/Sorted.bam \
../${PREV_MAP}/09-Col-0-Late-BR3_.Reference_2/Sorted.bam

date
echo "Tsu1 Early + Late"
ln -s ../${PREV_REF}/Tsu1.fasta prev.Tsu1.fasta
samtools view -b -H ../${PREV_MAP}/04-Tsu-1-Early-BR1.Reference_1/Sorted.bam > Tsu1.header.bam
echo MERGE
samtools merge -h Col0.header.bam -@ 4 -O BAM Tsu1.bam \
../${PREV_MAP}/04-Tsu-1-Early-BR1.Reference_1/Sorted.bam \
../${PREV_MAP}/04-Tsu-1-Early-BR1.Reference_2/Sorted.bam \
../${PREV_MAP}/05-Tsu-1-Early-BR2.Reference_1/Sorted.bam \
../${PREV_MAP}/05-Tsu-1-Early-BR2.Reference_2/Sorted.bam \
../${PREV_MAP}/06-Tsu-1-Early-BR3.Reference_1/Sorted.bam \
../${PREV_MAP}/06-Tsu-1-Early-BR3.Reference_2/Sorted.bam \
../${PREV_MAP}/10-Tsu-1-Late-BR1_.Reference_1/Sorted.bam \
../${PREV_MAP}/10-Tsu-1-Late-BR1_.Reference_2/Sorted.bam \
../${PREV_MAP}/11-Tsu-1-Late-BR2_.Reference_1/Sorted.bam \
../${PREV_MAP}/11-Tsu-1-Late-BR2_.Reference_2/Sorted.bam \
../${PREV_MAP}/12-Tsu-1-Late-BR3_.Reference_1/Sorted.bam \
../${PREV_MAP}/12-Tsu-1-Late-BR3_.Reference_2/Sorted.bam 

cd ..
date
