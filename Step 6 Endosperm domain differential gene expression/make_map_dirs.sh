#!/bin/sh

# Make a directory 'reference' that contains the ColAndTsu.fasta polished reference transcriptomes generated during step 1

DATA[1]="DAL1-5-2pos"
DATA[2]="DAL1-5-4pos"
DATA[3]="EE1pos"
DATA[4]="EE3pos"
DATA[5]="EE4pos"
DATA[6]="ESR2pos"
DATA[7]="ESR3pos"
DATA[8]="ESR4pos"
DATA[9]="TE1-2pos"
DATA[10]="TE1-3pos"
DATA[11]="TE1-4pos"

for DD in `seq 1 11`; do
    DATUM=${DATA[${DD}]}
    DIR="map_${DATUM}"
    echo $DIR
    rm -rf $DIR
    mkdir $DIR
    cd $DIR
    ln -s ../reference/* .
    ln -s /Trimmed/${DATUM}_*.fq.gz . # Location to trimmed marker reads
    cd ..
done
