#!/bin/tcsh

#$ -V
#$ -cwd
#$ -N snsphtools
#$ -j y
#$ -o tools.out
#$ -pe 1way 12
#$ -q serial
#$ -l h_rt=03:30:00
./addabundance2v3 $SCRATCH/runwvt/wvtfinal.sdf_float $SCRATCH/runsnsph/sn87a/sn87a_10M2.sdf 1 >& outfile.dat
