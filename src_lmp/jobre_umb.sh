#!/bin/bash -l
#PBS -l walltime=12:00:00,nodes=1:ppn=24,pmem=2580mb
#PBS -m abe
#PBS -M vsethura@umn.edu
cd ${PBS_O_WORKDIR}
echo job_start
mpirun -np 24 ./lmp_mesabi -in in.umbre -e screen
