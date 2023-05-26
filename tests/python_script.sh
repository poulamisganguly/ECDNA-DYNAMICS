#!/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l h_rt=1:0:0
#$ -l h_vmem=5G
#$ -N ecDNA_plots
#$ -j y
#$ -t 1-1

INPUT_ARGS=$(sed -n "${SGE_TASK_ID}p" python_list_of_args.txt)
module load anaconda3
conda activate biasedDoubling
python /data/home/ahw899/ECDNA-DYNAMICS/python/plotting_funcs.py $INPUT_ARGS
