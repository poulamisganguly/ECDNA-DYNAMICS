#!/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=1:0:0
#$ -l h_vmem=1G
#$ -N ecDNA_evolution
#$ -t 1-2
#$ -j y

INPUT_ARGS=$(sed -n "${SGE_TASK_ID}p" list_of_args.txt)
module load boost
module load gcc
make -C /data/home/ahw899/ECDNA-DYNAMICS/src clean
make -C /data/home/ahw899/ECDNA-DYNAMICS/src
/data/home/ahw899/ECDNA-DYNAMICS/src/ecDNA_evolution_2types.exe $INPUT_ARGS
