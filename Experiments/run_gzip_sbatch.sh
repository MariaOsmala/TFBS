#!/bin/bash
#SBATCH --job-name=gzip
#SBATCH --output=outs/gzip_%A_%a.out
#SBATCH --error=errs/gzip_%A_%a.out
#SBATCH --account=project_2006203
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=2G #max 382G
#SBATCH --cpus-per-task=1
#SBATCH --array=0


#DONE 0-4
#running 5, 35-39


srun gzip_better.sh ${SLURM_ARRAY_TASK_ID}




seff $SLURM_JOBID
