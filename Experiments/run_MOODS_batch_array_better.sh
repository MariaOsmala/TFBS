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
#SBATCH --array=301-393 # #0-393  version2.2 #version1 0-329 

#0-100 21186319
#101-200 21186322
#201-300 21186323
#301-393 21186325

#sacct --format JobID%-20,State -j 21137238

#Request local storage using the --gres flag in the job submission:
#--gres=nvme:<local_storage_space_per_node>
#The amount of space is given in GB (check maximum sizes from the list above). 
#For example, to request 100 GB of storage, use option --gres=nvme:100. 
#The local storage reservation is on a per node basis.

#Use the environment variable $LOCAL_SCRATCH in your batch job scripts to access the local storage on each node.

#srun MOODS_Teemu.sh ${SLURM_ARRAY_TASK_ID} #USE THIS
#srun MOODS_mouse.sh ${SLURM_ARRAY_TASK_ID} #USE THIS
#srun MOODS_final.sh ${SLURM_ARRAY_TASK_ID} #Manuscript version2.2
#srun MOODS_final_version2.2.sh ${SLURM_ARRAY_TASK_ID} #Version2 motifs which are not in union
srun MOODS_final_version2.2_correct.sh ${SLURM_ARRAY_TASK_ID} #Version2 motifs which are not in union
seff $SLURM_JOBID
