#!/bin/sh

# Command line parameter 1 gives the number for this round.
# Give a zero for the initial round where we map to Araport11.
# This will assume the target is in the consensus_0 directory.
# The consensus_0 directory should contain Col0.fasta and Col0.*.bt2
# as well as Tsu1.fasta and Tsu1.*.bt2.

if [ $# -eq 0 ]; then
    echo "Please provide the number of this round."
    exit 1
fi
ROUND=$1
echo ROUND $ROUND

TRIM_BASE='/Trimmed' #Location of trimmed reads
CONS_BASE='/HomozygousConsensus' #Location where pilon polishing of reference transcriptomes takes place
CONS_DIR="consensus_${ROUND}"
MAP_DIR="map_to_consensus_${ROUND}"
mkdir ${MAP_DIR}
cd ${MAP_DIR}
pwd

DATA1[0]="01-Col-0-Early-BR1"
DATA1[1]="02-Col-0-Early-BR2"
DATA1[2]="03-Col-0-Early-BR3"
DATA1[3]="07-Col-0-Late-BR1_"
DATA1[4]="08-Col-0-Late-BR2_"
DATA1[5]="09-Col-0-Late-BR3_"
DATA1[6]="04-Tsu-1-Early-BR1"
DATA1[7]="05-Tsu-1-Early-BR2"
DATA1[8]="06-Tsu-1-Early-BR3"
DATA1[9]="10-Tsu-1-Late-BR1_"
DATA1[10]="11-Tsu-1-Late-BR2_"
DATA1[11]="12-Tsu-1-Late-BR3_"

for DD in `seq 0 11` ; do
    SAMPLE=${DATA1[${DD}]}
    DIR1="${SAMPLE}.Reference_1"
    rm -rf $DIR1
    mkdir $DIR1
    cd $DIR1
    if [ $DD -lt 6 ]; then
         ln -s ../../${CONS_DIR}/Col0.* .
    else
         ln -s ../../${CONS_DIR}/Tsu1.* .
    fi
    ln -s ${TRIM_BASE}/Reference_1/*${SAMPLE}*.fq.gz .
    cd ..
    # 
    DIR2="${SAMPLE}.Reference_2"
    rm -rf $DIR2
    mkdir $DIR2
    cd $DIR2
    if [ $DD -lt 6 ]; then
         ln -s ../../${CONS_DIR}/Col0* .
    else
         ln -s ../../${CONS_DIR}/Tsu1* .
    fi
    ln -s ${TRIM_BASE}/Reference_2/*${SAMPLE}*.fq.gz .
    cd ..
done

cd ..
echo DONE
