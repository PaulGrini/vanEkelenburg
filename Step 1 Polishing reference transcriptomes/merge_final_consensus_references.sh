#!/bin/sh

# Saga setup
# Constrained to use GCC 8.2 in order to use samtools & bowtie

echo SETUP
module --force purge
module load StdEnv 
module load GCC/8.2.0-2.31.1 
module load Python/3.7.2-GCCcore-8.2.0
module load SAMtools/1.9-GCC-8.2.0-2.31.1
module load Bowtie2/2.3.5.1-GCC-8.2.0-2.31.1

echo DIRECTORY
mkdir final_reference
cd final_reference

CONSENSUS_DIR= # Location to the final consensus_ directory
PILON_SUFFIX="_pilon_pilon_pilon_pilon_pilon_pilon_pilon_pilon_pilon_pilon"

echo EXTRACT COL
cat ${CONSENSUS_DIR}/Col0.fasta | \
sed 's/_pilon_pilon_pilon_pilon_pilon_pilon_pilon_pilon_pilon_pilon//' \
> ColAndTsu.fasta

echo EXTRACT TSU
cat ${CONSENSUS_DIR}/Tsu1.fasta | \
sed 's/_pilon_pilon_pilon_pilon_pilon_pilon_pilon_pilon_pilon_pilon//' \
>> ColAndTsu.fasta

echo INDEX
# It would cause headaches later if index != fasta name.
bowtie2-build ColAndTsu.fasta ColAndTsu

cd ..
echo DONE


