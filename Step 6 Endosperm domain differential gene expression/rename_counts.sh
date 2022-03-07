#!/bin/sh

# Rename counts.tsv to ${DATUM}.tsv

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

# Rename counts.tsv to ${DATUM}.tsv

for DD in `seq 1 11`; do
    DATUM=${DATA[${DD}]}
    DIR="map_${DATUM}"
    echo $DIR
    cd $DIR
    mv counts.tsv count_by_gene_${DATUM}.tsv
    cd ..
done
