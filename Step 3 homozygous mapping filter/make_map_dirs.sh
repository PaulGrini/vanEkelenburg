#!/bin/sh

# Reference_1 = 20190705
# Reference_2 = 20190914
# Prepare a folder '/reference' that contains the polished Col-0 and Tsu-' reference transcriptomes

DATA1[0]="01-Col-0-Early-BR1"
DATA1[1]="02-Col-0-Early-BR2"
DATA1[2]="03-Col-0-Early-BR3"
DATA1[3]="04-Tsu-1-Early-BR1"
DATA1[4]="05-Tsu-1-Early-BR2"
DATA1[5]="06-Tsu-1-Early-BR3"
DATA1[6]="07-Col-0-Late-BR1_"
DATA1[7]="08-Col-0-Late-BR2_"
DATA1[8]="09-Col-0-Late-BR3_"
DATA1[9]="10-Tsu-1-Late-BR1_"
DATA1[10]="11-Tsu-1-Late-BR2_"
DATA1[11]="12-Tsu-1-Late-BR3_"

for DD in `seq 0 11` ; do
    SAMPLE=${DATA1[${DD}]}
    DIR1="${SAMPLE}.20190705"
    rm -rf $DIR1
    mkdir $DIR1
    cd $DIR1
    ln -s ../reference/* .
    ln -s /Trimmed/Reference_1/*${SAMPLE}*.fq.gz . #Location of trimmed reads
    cd ..
    # 
    DIR2="${SAMPLE}.20190914"
    rm -rf $DIR2
    mkdir $DIR2
    cd $DIR2
    ln -s ../reference/* .
    ln -s /Trimmed/Reference_2/*${SAMPLE}*.fq.gz . #Location of trimmed reads
    cd ..
done

echo DONE
