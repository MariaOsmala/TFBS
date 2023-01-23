#!/bin/bash
#SBATCH --job-name=MOODS
#SBATCH --output=outs/MOODS_%A_%a.out
#SBATCH --error=errs/MOODS_%A_%a.err
#SBATCH --account=project_2006203
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=15:00
#SBATCH --mem-per-cpu=5G #max 382G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-399# 0-399 #399



#srun MOODS_better.sh ${SLURM_ARRAY_TASK_ID}

#srun MOODS_Vierstra.sh ${SLURM_ARRAY_TASK_ID}
srun MOODS_Teemu.sh ${SLURM_ARRAY_TASK_ID}


seff $SLURM_JOBID
