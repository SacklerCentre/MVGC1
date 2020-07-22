#!/bin/sh

echo "Starting job 'accu_test_n' at `date`"

#$ -j y
#$ -cwd
#$ -N accu_test_n
#$ -o accu_test_n.out
#$ -S /bin/bash
#$ -q serial_lowmem.q
#$ -pe openmp 8

cd $LOCALREPO/matlab/ncomp/mvgc

matlab -nojvm -nodisplay -r "cd(fullfile(getenv('DATADIR'),'accu_test')); maxNumCompThreads(8); accu_test_n"

echo "Finished job 'accu_test_n' at `date`"
