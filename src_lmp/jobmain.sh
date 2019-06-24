#!/bin/bash -l
#PBS -l walltime=24:00:00,nodes=2:ppn=24,pmem=2580mb
#PBS -m abe
#PBS -M vsethura@umn.edu
cd ${PBS_O_WORKDIR}
echo job_start
mpirun -np 48 ./lmp_mesabi -in in.init -e screen
wait
mpirun -np 48 ./lmp_mesabi -in in.run1 -e screen
wait
mpirun -np 48 ./lmp_mesabi -in in.longrun -e screen
