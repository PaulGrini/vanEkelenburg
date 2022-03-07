#!/bin/sh

# Process every heterozygous bam file.
# Create text files of maps per read.
# The --debug flag generates stats on stderr.

echo "EE"
sbatch --account=${ACCOUNT} run_diff.sh EE1pos
sbatch --account=${ACCOUNT} run_diff.sh EE3pos
sbatch --account=${ACCOUNT} run_diff.sh EE4pos
echo "ESR"
sbatch --account=${ACCOUNT} run_diff.sh ESR2pos
sbatch --account=${ACCOUNT} run_diff.sh ESR3pos
sbatch --account=${ACCOUNT} run_diff.sh ESR4pos
echo "TE1"
sbatch --account=${ACCOUNT} run_diff.sh TE1-2pos
sbatch --account=${ACCOUNT} run_diff.sh TE1-3pos
sbatch --account=${ACCOUNT} run_diff.sh TE1-4pos
echo "DONE"
