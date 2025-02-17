#!/bin/bash -l
#SBATCH --job-name=SSTAT
#SBATCH --output=outs/SSTAT_%A_%a.txt
#SBATCH --error=errs/SSTAT_%A_%a.txt
#SBATCH --account=project_2006472
#SBATCH --partition=small
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --gres=nvme:20 
#SBATCH --array=0-1000

mkdir ../../Results/MOSTA-SSTAT

srun run_SSTAT_version2.2.sh ${SLURM_ARRAY_TASK_ID} 

seff $SLURM_JOBID


