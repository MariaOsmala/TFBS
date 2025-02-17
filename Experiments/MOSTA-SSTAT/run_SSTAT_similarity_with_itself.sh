#!/bin/bash -l
#SBATCH --job-name=SSTAT_itself
#SBATCH --output=outs/SSTAT_itself.txt
#SBATCH --error=errs/SSTAT_itself.txt
#SBATCH --account=project_2006203
#SBATCH --partition=small
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G


srun run_SSTAT_with_itself_version2.2.sh 

seff $SLURM_JOBID


