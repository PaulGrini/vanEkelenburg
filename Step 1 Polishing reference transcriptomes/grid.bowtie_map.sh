#!/bin/sh

# Command line parameter gives the number for this round.
# Give a zero for the initial round where we map to Araport11.
# Give a 1 to map to consensus_1, etc.

if [ $# -eq 0 ]; then
    echo "Please provide the number of this round."
    exit 1
fi
ROUND=$1
echo ROUND $ROUND
TRIM_BASE='/Trimmed' #Location of trimmed reads
CONS_DIR="consensus_${ROUND}"
MAP_DIR="map_to_consensus_${ROUND}"

echo "DIRECTORY SETUP"
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

echo
echo "GRID SUBMIT"

export ACCOUNT= #Account name/number to submit to the server
export MOLBAR_HOME= #Location of scripts

# Visit every subdirectory
for D in ??-???-*/;
do
    cd $D
    sbatch --account=${ACCOUNT} ../../bowtie_map.sh
    cd ..
done
echo DONE
cd ..
