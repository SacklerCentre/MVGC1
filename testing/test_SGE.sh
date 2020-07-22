#!/bin/sh

echo "Starting job 'test_SGE.sh' at `date`"

#$ -j y
#$ -q inf.q
#$ -cwd
#$ -N test_SGE
#$ -o test_SGE.out
#$ -S /bin/bash

export rundir=`pwd`

cd $LOCALREPO/matlab/ncomp/mvgc

matlab -nojvm -nodisplay -r "cd(getenv('rundir')); pwd"
