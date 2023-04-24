#!/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=1:0:0
#$ -l h_vmem=1G
#$ -N ecDNA_evolution

cd ../src
make clean
make 
./ecDNA_evolution_2types.exe
