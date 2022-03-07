#!/bin/sh

echo modules
module purge
module load Python/3.6.6-foss-2018b
module load R-bundle-Bioconductor/3.8-foss-2018b-R-3.5.1

p=python3
s= # Location to scripts (src) on MOLBAR_HOME

# This script will generate csv files.
# There shouldn't be any csv files in this directory.
# If there are, they are assumed to be input and lead to redundant outputs.
# In case we rerun this script, we need to get rid of the previously generated csv files.
echo
echo clean up
rm -v *.log
rm -v *.csv

#ColTsu.model.tsv contains normalization factors obtained from step 5
#Copy $.tsv files generated during step 4 to this directory
echo
echo prepare
# Collate three replicates into one row per gene like MAT,MAT,MAT,PAT,PAT,PAT.
# Generate several files like this one called pos.EE1.counts.csv
# AT4G38740.1,26,36,39,33,37,30
# AT1G10760.1,87,343,286,134,291,274
# ...
$p ${s}/prepare_heterozygous_counts.py ColTsu.model.tsv 
echo -n $?
echo " exit status"
echo
echo "Prepared these csv files:"
ls *.csv

# Place previously generated filter lists (step 2 and 3) in this directory 

echo 
echo filter
function filter() {
    echo "Filter $1 excluding $2"
    cat ${1}.counts.csv | grep -v '^ATMG' > ${1}.mito
    cat ${1}.mito | grep -v '^ATCG' > ${1}.noplastid
    ${p} ${s}/apply_filter_to_counts.py ${2} ${1}.noplastid > ${1}.homozygous
    ${p} ${s}/apply_filter_to_counts.py ${3} ${1}.homozygous > ${1}.filtered
    #echo "Normalize $1"
    #${p} ${s}/normalize_counts.py --debug ${1}.filtered ${1}.normalized
    echo "Statistics $1"
    ${p} ${s}/statistics_from_counts.py --debug ${s}/limma.foldchange.r ${1}.filtered
}
filter neg.EE1 Early.genes_pass_filter Early_ecotype_pass_filter
filter pos.EE1 Early.genes_pass_filter Early_ecotype_pass_filter
filter neg.ESR Late.genes_pass_filter Late_ecotype_pass_filter
filter neg.TE1 Late.genes_pass_filter Late_ecotype_pass_filter
filter pos.ESR Late.genes_pass_filter Late_ecotype_pass_filter
filter pos.TE1 Late.genes_pass_filter Late_ecotype_pass_filter

echo
echo clean up
rm -v *.json
# statistics_from_counts.py leaves behind timestamps in *json files
#echo JUST FOR DEBUGGING NO clean up

echo
echo "Look here for P-value statistics."
ls -l *.final.csv

echo
echo done
