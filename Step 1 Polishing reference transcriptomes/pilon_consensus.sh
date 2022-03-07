#!/bin/sh

#SBATCH --account=$ACCOUNT
#SBATCH --job-name=FACSP
#SBATCH --time=24:00:00 # pilon takes about 5 hours
#SBATCH --mem-per-cpu=18G  # pilon can require 18 GB RAM
#SBATCH --cpus-per-task=1  # pilone multi-threading is inefficient
# Jobs submit with 'arrayrun 1-10' receive a $TASK_ID.
# Jobs submit with 'sbatch --array=1-10' recieve a $SLURM_ARRAY_TASK_ID.

set -o errexit # exit on errors
# Pilon expected outputs
#savefile pilon.changes
#savefile pilon.fasta 

if [ $# -ne 2 ]; then
    echo "Please specify <Col0|Tsu1> and number of this round."
    exit 1
fi
GENOME=$1
echo GENOME $GENOME
ROUND=$2
echo ROUND $ROUND
let PREV=${ROUND}-1
PREV_REF=consensus_${PREV}
PREV_MAP=map_to_consensus_${PREV}
echo PREV_REF $PREV_REF
echo PREV_MAP $PREV_MAP
INITIALDIR=`pwd`
echo INITIALDIR ${INITIALDIR}
echo

echo MODULES
module purge
module load SAMtools/1.9-GCC-8.2.0-2.31.1   
module load Bowtie2/2.3.5.1-GCC-8.2.0-2.31.1
module load Pilon/1.23-Java-11  ## saga
# To execute Pilon run: java -Xmx8G -jar $EBROOTPILON/pilon.jar
# Pilon Requirement: fasta file, bam file, bai index
# Pilon Requirement: about 18 GB RAM
JARFILE=$EBROOTPILON/pilon.jar
echo PILON JAR; ls -l $JARFILE
echo HEADER; ls -l $HEADERFILE
echo "Java version is:"
JAVA=`which java`
${JAVA} -showversion |& head -n 4
PILON="${JAVA} -Xmx16G -jar ${JARFILE}"
echo "Pilon command is: ${PILON}"
echo

# Here are some interesting pilon options
OPTIONS="--diploid "  # assume diploid
OPTIONS="--vcf"       # output a VCF
OPTIONS="--variant"   # heuristic for variants not assembly
OPTIONS="--fix all"   # correct bases, indels, and local misassembly and write a FASTA file
OPTIONS="--changes"   # show changes to the FASTA
# Use this set of options to generate the homozygous genome informed by reads.
OPTIONS="--fix all --changes --output ${GENOME}"
echo "Pilon options set to: ${OPTIONS}"
echo

echo "SAMTOOLS MERGE"
if [ "${GENOME}" == "Col0" ] ; then
    echo "Col0 Early + Late"
    ln -s ../${PREV_REF}/Col0.fasta prev.Col0.fasta
    samtools view -b -H ../${PREV_MAP}/01-Col-0-Early-BR1.Reference_1/Sorted.bam > Col0.header.bam
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
fi
if [ "${GENOME}" == "Tsu1" ] ; then
    echo "Tsu1 Early + Late"
    ln -s ../${PREV_REF}/Tsu1.fasta prev.Tsu1.fasta
    samtools view -b -H ../${PREV_MAP}/04-Tsu-1-Early-BR1.Reference_1/Sorted.bam > Tsu1.header.bam
    samtools merge -h Tsu1.header.bam -@ 4 -O BAM Tsu1.bam \
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
fi
echo

date
echo SAMTOOLS INDEX
samtools index -@ 4 ${GENOME}.bam
echo

date
echo START PILON
${PILON} --genome prev.${GENOME}.fasta --frags ${GENOME}.bam ${OPTIONS}
echo -n $?; echo " exit status"
echo DONE PILON
date
echo

echo START BOWTIE BUILD
bowtie2-build ${GENOME}.fasta ${GENOME}
date
echo DONE
