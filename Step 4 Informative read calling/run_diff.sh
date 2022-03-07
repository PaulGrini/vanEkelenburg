#!/bin/sh

# Process every heterozygous bam file.
# Create text files of maps per read.
# The --debug flag generates stats on stderr.

#SBATCH --account=${ACCOUNT}
#SBATCH --job-name=Take2
#SBATCH --time=04:00:00   # Bowtie takes about 35 minutes
#SBATCH --mem-per-cpu=4G  # 16 GB total
#SBATCH --cpus-per-task=4  # 4 cpu is optimal for 4 threads

set -o errexit # exit on errors

module --force purge
module load StdEnv 
module load GCC/8.2.0-2.31.1 
module load Python/3.7.2-GCCcore-8.2.0
module load SAMtools/1.9-GCC-8.2.0-2.31.1
module load Bowtie2/2.3.5.1-GCC-8.2.0-2.31.1
module list

export SRC= # Location to scripts in MOLBAR_HOME

echo $1

date
echo differential_mapping.py
python3 ${SRC}/differential_mapping.py  map_${1}/Primary.bam --debug 1> ${1}.out 2> ${1}.err

date
echo count_diffmap_per_gene.py
python3 ${SRC}/count_diffmap_per_gene.py  ${1}.out > ${1}.tsv


echo "DONE"
date
