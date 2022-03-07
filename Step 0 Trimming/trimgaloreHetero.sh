#!/bin/bash
#SBATCH --job-name=trimgalore
#SBATCH --account=
#SBATCH --output=slurm-%j.base
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=6
##SBATCH --nodes=1
#SBATCH --mem-per-cpu=8G
##SBATCH --partition=normal

module purge
module load Trim_Galore/0.6.5-GCCcore-8.3.0-Java-11-Python-3.7.4
/TrimGalore-0.6.5Dev/TrimGalore_Dev/trim_galore $1 $2 --paired -a "file:/adapters.fa" -a2 "file:/adapters.fa" -j 6
