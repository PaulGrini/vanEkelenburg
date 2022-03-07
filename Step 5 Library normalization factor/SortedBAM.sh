#!/bin/sh 
#SBATCH --account=nn9525k
#SBATCH --job-name=Sort
#SBATCH --time=04:00:00   
#SBATCH --mem-per-cpu=4G  # 16 GB total
#SBATCH --cpus-per-task=4  # 4 cpu is optimal for 4 threads

# Obtain Sorted bam files

THREADS=4
# Visit every subdirectory
for D in map_*/;
do
    cd $D
    samtools sort --threads $THREADS -T tmp --output-fmt BAM -o Sorted.bam Primary.bam
    cd ..
done
echo DONE
