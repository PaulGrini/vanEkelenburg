#!/bin/sh

DATA[1]="EE1pos"
DATA[2]="EE3pos"
DATA[3]="EE4pos"
DATA[4]="ESR2pos"
DATA[5]="ESR3pos"
DATA[6]="ESR4pos"
DATA[7]="TE1-2pos"
DATA[8]="TE1-3pos"
DATA[9]="TE1-4pos"

for DD in `seq 1 9`; do
    DATUM=${DATA[${DD}]}
    DIR="map_${DATUM}"
    echo $DIR
    rm -rf $DIR
    mkdir $DIR
    cd $DIR
    rm Sorted.bam
    ln -s /map_${DATUM}/Primary.bam . # Location for Primary.bam files after mapping of reads performed during step 4
    cd ..
done
