#!/bin/bash -l
#SBATCH --job-name=MOODS_artificial
#SBATCH --output=outs/MOODS_artificial_%A_%a.out
#SBATCH --error=errs/MOODS_artificial_%A_%a.err
#SBATCH --account=project_2006203
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=200G
#SBATCH --gres=nvme:10
#SBATCH --array=0  #0:786


#srun MOODS_artificial_motifs.sh ${SLURM_ARRAY_TASK_ID}
srun MOODS_artificial_motifs_artificialHTSELEX_halfsites.sh ${SLURM_ARRAY_TASK_ID}

echo "success"

seff $SLURM_JOBID
