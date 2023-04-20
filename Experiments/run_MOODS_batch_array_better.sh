#!/bin/bash
#SBATCH --job-name=MOODS
#SBATCH --output=outs/MOODS_%A_%a.out
#SBATCH --error=errs/MOODS_%A_%a.err
#SBATCH --account=project_2006203
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=5G #
#SBATCH --cpus-per-task=1
#SBATCH --array=0-57 #1-102# 0-399 #399 0-102 0-76 0-57


#files=($(awk -F '"' '{print $4}' rerunMOODS_lower_threshold.txt)) #1024 ->102
#files=($(awk -F '"' '{print $4}' rerunMOODS_lower_threshold_4.txt)) #761 ->76
#files=($(awk -F '"' '{print $4}' rerunMOODS_lower_threshold_3.txt)) #578 -> 57


#srun MOODS_better.sh ${SLURM_ARRAY_TASK_ID}

#srun MOODS_Vierstra.sh ${SLURM_ARRAY_TASK_ID}
#srun MOODS_Teemu.sh ${SLURM_ARRAY_TASK_ID}
srun MOODS_Teemu_lower_threshold.sh ${SLURM_ARRAY_TASK_ID}


seff $SLURM_JOBID
