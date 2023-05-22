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
#SBATCH --gres=nvme:50
#SBATCH --array=1-100 #0-57 #1-102# 0-399 #399 0-102 0-76 0-57


#Request local storage using the --gres flag in the job submission:
#--gres=nvme:<local_storage_space_per_node>
#The amount of space is given in GB (check maximum sizes from the list above). 
#For example, to request 100 GB of storage, use option --gres=nvme:100. 
#The local storage reservation is on a per node basis.

#Use the environment variable $LOCAL_SCRATCH in your batch job scripts to access the local storage on each node.


#files=($(awk -F '"' '{print $4}' rerunMOODS_lower_threshold.txt)) #1024 ->102
#files=($(awk -F '"' '{print $4}' rerunMOODS_lower_threshold_4.txt)) #761 ->76
#files=($(awk -F '"' '{print $4}' rerunMOODS_lower_threshold_3.txt)) #578 -> 57


#srun MOODS_better.sh ${SLURM_ARRAY_TASK_ID}

#srun MOODS_Vierstra.sh ${SLURM_ARRAY_TASK_ID}
#srun MOODS_Teemu.sh ${SLURM_ARRAY_TASK_ID}
#srun MOODS_Teemu_lower_threshold.sh ${SLURM_ARRAY_TASK_ID}
#srun MOODS_Teemu_Morgunova2015_bHLH_HOMEOBOXYimeng.sh ${SLURM_ARRAY_TASK_ID}
#MOODS_Teemu_Morgunova2015_bHLH_HOMEOBOXYimeng.sh
#srun MOODS_mouse.sh ${SLURM_ARRAY_TASK_ID}
srun MOODS_mouse_lower_threshold.sh ${SLURM_ARRAY_TASK_ID} #0-100
seff $SLURM_JOBID
