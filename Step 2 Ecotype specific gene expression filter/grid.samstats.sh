#!/bin/sh

# Run this from a directory such as Heterozygous.
# Expect scripts in parent directory.
# Expect several subdirectories created by make_map_dir.sh

export ACCOUNT= #Account
export PYTHON_VENV= #Location to virtual environment
export MOLBAR_HOME= #Location to scripts
# Visit every subdirectory
for D in ??-???-*/;
do
    cd $D
    sbatch --account=${ACCOUNT} ../samstats.sh
    cd ..
done
echo DONE
