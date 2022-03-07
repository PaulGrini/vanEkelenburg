#!/bin/sh

# Process every heterozygous bam file.
# Create text files of maps per read.
# The --debug flag generates stats on stderr.

echo "Col Early"
sbatch --account=${ACCOUNT} run_diff.sh 01-Col-0-Early-BR1.20190705
sbatch --account=${ACCOUNT} run_diff.sh 01-Col-0-Early-BR1.20190914
sbatch --account=${ACCOUNT} run_diff.sh 02-Col-0-Early-BR2.20190705
sbatch --account=${ACCOUNT} run_diff.sh 02-Col-0-Early-BR2.20190914
sbatch --account=${ACCOUNT} run_diff.sh 03-Col-0-Early-BR3.20190705
sbatch --account=${ACCOUNT} run_diff.sh 03-Col-0-Early-BR3.20190914
echo "Tsu Early"
sbatch --account=${ACCOUNT} run_diff.sh 04-Tsu-1-Early-BR1.20190705
sbatch --account=${ACCOUNT} run_diff.sh 04-Tsu-1-Early-BR1.20190914
sbatch --account=${ACCOUNT} run_diff.sh 05-Tsu-1-Early-BR2.20190705
sbatch --account=${ACCOUNT} run_diff.sh 05-Tsu-1-Early-BR2.20190914
sbatch --account=${ACCOUNT} run_diff.sh 06-Tsu-1-Early-BR3.20190705
sbatch --account=${ACCOUNT} run_diff.sh 06-Tsu-1-Early-BR3.20190914
echo "Col Late"
sbatch --account=${ACCOUNT} run_diff.sh 07-Col-0-Late-BR1_.20190705
sbatch --account=${ACCOUNT} run_diff.sh 07-Col-0-Late-BR1_.20190914
sbatch --account=${ACCOUNT} run_diff.sh 08-Col-0-Late-BR2_.20190705
sbatch --account=${ACCOUNT} run_diff.sh 08-Col-0-Late-BR2_.20190914
sbatch --account=${ACCOUNT} run_diff.sh 09-Col-0-Late-BR3_.20190705
sbatch --account=${ACCOUNT} run_diff.sh 09-Col-0-Late-BR3_.20190914
echo "Tsu Late"
sbatch --account=${ACCOUNT} run_diff.sh 10-Tsu-1-Late-BR1_.20190705
sbatch --account=${ACCOUNT} run_diff.sh 10-Tsu-1-Late-BR1_.20190914
sbatch --account=${ACCOUNT} run_diff.sh 11-Tsu-1-Late-BR2_.20190705
sbatch --account=${ACCOUNT} run_diff.sh 11-Tsu-1-Late-BR2_.20190914
sbatch --account=${ACCOUNT} run_diff.sh 12-Tsu-1-Late-BR3_.20190705
sbatch --account=${ACCOUNT} run_diff.sh 12-Tsu-1-Late-BR3_.20190914
echo "DONE"
