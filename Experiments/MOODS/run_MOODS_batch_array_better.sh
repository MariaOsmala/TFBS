#!/bin/bash
#SBATCH --job-name=MOODS
#SBATCH --output=outs/MOODS_%A_%a.out
#SBATCH --error=errs/MOODS_%A_%a.err
#SBATCH --account=project_2006472
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=80G #
#SBATCH --cpus-per-task=1
#SBATCH --gres=nvme:50
#SBATCH --array=1-393 # #0-393  


#Request local storage using the --gres flag in the job submission:
#--gres=nvme:<local_storage_space_per_node>
#The amount of space is given in GB (check maximum sizes from the list above). 
#For example, to request 100 GB of storage, use option --gres=nvme:100. 
#The local storage reservation is on a per node basis.

#Use the environment variable $LOCAL_SCRATCH in your batch job scripts to access the local storage on each node.

srun MOODS_final.sh ${SLURM_ARRAY_TASK_ID} #Final version2.2
#srun MOODS_final_artificialHTSelex.sh ${SLURM_ARRAY_TASK_ID} #Not an array job

seff $SLURM_JOBID
