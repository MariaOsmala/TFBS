#!/bin/bash -l
#SBATCH --job-name=true_vs_artificial
#SBATCH --output=outs/true_vs_artificial_%A_%a.txt
#SBATCH --error=errs/true_vs_artificial_%A_%a.txt
#SBATCH --account=project_2007567
#SBATCH --partition=small
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5
#SBATCH --gres=nvme:2 
#SBATCH --array=393 #0-393 

#This is to compute similarities between the true and the corresponding artificial motifs
  
srun run_motifsimilarity_true_vs_corresponding_artificial.sh ${SLURM_ARRAY_TASK_ID} #Old

seff $SLURM_JOBID



