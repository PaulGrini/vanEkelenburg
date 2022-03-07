#!/bin/sh

echo modules
module purge
module load Python/3.6.6-foss-2018b
module load R-bundle-Bioconductor/3.8-foss-2018b-R-3.5.1

p=python3
s= # Location to scripts (src) at MOLBAR_HOME

# This script will generate csv files.
# There shouldn't be any csv files in this directory.
# If there are, they are assumed to be input and lead to redundant outputs.
# In case we rerun this script, we need to get rid of the previously generated csv files.
echo
echo clean up

${p} ${s}/normalize_counts.py --debug pos.EE1.counts.csv EE1pos.out.csv > EE1pos.normfactors.csv
${p} ${s}/normalize_counts.py --debug pos.ESR.counts.csv ESRpos.out.csv > ESRpos.normfactors.csv
${p} ${s}/normalize_counts.py --debug pos.TE1.counts.csv TE1pos.out.csv > TE1pos.normfactors.csv

echo
echo done
