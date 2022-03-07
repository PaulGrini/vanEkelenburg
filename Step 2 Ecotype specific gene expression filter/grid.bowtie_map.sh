#!/bin/sh

# Run this from a directory such as map_bowtie_0.
# Expect scripts in parent directory.
# Expect several subdirectories created by make_map_dir.sh

export ACCOUNT= #Account
export PYTHON_VENV= #Location of virtual environment
export MOLBAR_HOME= #Location of scripts
# Visit every subdirector
for D in ??-???-*/;
do
    cd $D
    sbatch --account=${ACCOUNT} ../bowtie_map.sh
    cd ..
done
echo DONE
