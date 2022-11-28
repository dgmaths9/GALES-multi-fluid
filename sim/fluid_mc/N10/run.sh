#!/bin/bash


#SBATCH --job-name=sp
#SBATCH --ntasks=192
#SBATCH -p CLUSTER
#SBATCH --time=01:00:00


mpirun --mca btl vader,self,openib  ./executable 2> log.txt
