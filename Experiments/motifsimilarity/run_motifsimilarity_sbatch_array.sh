#!/bin/bash -l
#SBATCH --job-name=motifsimilarity
#SBATCH --output=outs/motifsimilarity_%A_%a.txt
#SBATCH --error=errs/motifsimilarity_%A_%a.txt
#SBATCH --account=project_2006472
#SBATCH --partition=small
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5
#SBATCH --gres=nvme:20
#SBATCH --array=701-1000


srun run_motifsimilarity_final_version2.2.sh ${SLURM_ARRAY_TASK_ID}

seff $SLURM_JOBID



