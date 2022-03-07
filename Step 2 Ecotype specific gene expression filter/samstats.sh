#!/bin/sh 
#SBATCH --account=${ACCOUNT}
#SBATCH --job-name=Het3a
#SBATCH --time=04:00:00   # Probably need 1 hour
#SBATCH --mem-per-cpu=4G  # total = per-cpu * cpus
#SBATCH --cpus-per-task=1  # samtools view in 1 thread
## source /cluster/bin/jobsetup   ## abel only
set -o errexit # exit on errors
# Our python will generate smaller Aligned.bam which we keep.
savefile counts.tsv

module --force purge
module load StdEnv 
module load GCC/8.2.0-2.31.1 
module load Python/3.7.2-GCCcore-8.2.0
module load SAMtools/1.9-GCC-8.2.0-2.31.1
module load Bowtie2/2.3.5.1-GCC-8.2.0-2.31.1
module list
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/Python/3.7.2-GCCcore-8.2.0/lib/
echo LD_LIBRARY_PATH $LD_LIBRARY_PATH
PYTHON_VENV= #Location to virtual environment
echo PYTHON_VENV ${PYTHON_VENV}
source ${PYTHON_VENV}/activate
echo WHICH PYTHON3
which python3
python3 --version
MOLBAR_HOME= #Location to scripts
echo MOLBAR_HOME ${MOLBAR_HOME}

date
echo COPY INPUTS TO GRID 
echo RUN THIS FROM THE DIRECTORY CONTAINING R1 AND R2.fastq.gz
INITIALDIR=`pwd`
echo INITIALDIR ${INITIALDIR}
echo SCRATCH ${SCRATCH}
cd ${SCRATCH}
cp -vHR ${INITIALDIR}/Sorted.bam .
date
echo run counts
python3 ${MOLBAR_HOME}/src/samstats.py --debug --file Sorted.bam > counts.tsv
echo -n $?
echo " exit status"
date
ls -l
echo DONE
