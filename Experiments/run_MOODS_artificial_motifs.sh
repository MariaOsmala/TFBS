#!/bin/bash
#SBATCH --job-name=MOODS_artificial
#SBATCH --output=outs/MOODS_artificial_%A_%a.out
#SBATCH --error=errs/MOODS_artificial_%A_%a.err
#SBATCH --account=project_2007567
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=100G #
#SBATCH --cpus-per-task=1
#SBATCH --gres=nvme:50
#SBATCH --array=644

#Request local storage using the --gres flag in the job submission:
#--gres=nvme:<local_storage_space_per_node>
#The amount of space is given in GB (check maximum sizes from the list above). 
#For example, to request 100 GB of storage, use option --gres=nvme:100. 
#The local storage reservation is on a per node basis.

#Use the environment variable $LOCAL_SCRATCH in your batch job scripts to access the local storage on each node.

srun MOODS_artificial_motifs.sh ${SLURM_ARRAY_TASK_ID}

seff $SLURM_JOBID
