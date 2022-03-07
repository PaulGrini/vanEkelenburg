#!/bin/sh

# Launch pilon on the grid.
# Run this after make_consensus_dir.sh

# Command line parameter 1 gives the number for this round.
# Give a one for the initial round where pilon runs for the first time.

# After last pilon round, merge the polished reference transcriptomes to generate one ColAndTsu.fasta file

if [ $# -eq 0 ]; then
    echo "Please provide the number of this round."
    exit 1
fi
ROUND=$1
WORKDIR=consensus_${ROUND}
echo WORKDIR $WORKDIR
mkdir "$WORKDIR"
if [ $? -ne 0 ] ; then
    echo "FAIL: Directory already exists."
    exit 2
fi

export ACCOUNT= #Account name/number to submit to the server
# It is critical to launch the script in the subdirectory.
# The launch determines where the pilon.fasta result gets written.
# Launching two jobs in one directory leads to overrittedn files.

cd $WORKDIR
sbatch --account=${ACCOUNT} ../pilon_consensus.sh Col0 $ROUND
sbatch --account=${ACCOUNT} ../pilon_consensus.sh Tsu1 $ROUND
cd ..
