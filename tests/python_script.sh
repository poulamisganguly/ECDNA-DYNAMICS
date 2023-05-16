#!/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=1:0:0
#$ -l h_vmem=5G
#$ -N ecDNA_plots
#$ -j y

module load anaconda3
conda activate biasedDoubling
python /data/home/ahw899/ECDNA-DYNAMICS/python/make_plots.py 
