#!/bin/bash
#SBATCH --job-name=rows
#SBATCH --output=outs/rows_%A_%a.out
#SBATCH --error=errs/rows_%A_%a.out
#SBATCH --account=project_2006203
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=1G #max 382G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-39





srun rows.sh ${SLURM_ARRAY_TASK_ID}




seff $SLURM_JOBID
