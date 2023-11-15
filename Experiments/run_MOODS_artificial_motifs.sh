#!/bin/bash -l
#SBATCH --job-name=MOODS_artificial
#SBATCH --output=outs/MOODS_artificial_%A_%a.out
#SBATCH --error=errs/MOODS_artificial_%A_%a.err
#SBATCH --account=project_2007567
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=08:00:00
#SBATCH --mem-per-cpu=150G 
#SBATCH --cpus-per-task=1
#SBATCH --gres=nvme:50
#SBATCH --array=336

#   19014924_[336]

#serer 300,307,308 #-643,645-765, 767-1030 #0-1030 #644, 766, 1-100, 101-400, 401-600 901-1030 


# 18903955_[1-67] #OK

# 18919054_[68-100] #OK

# 18921005_[101-300] #OK

# 18921006_[301-400] #OK

#18930251_[401-600] #OK

#18930254_[601-700] #OK

#18943905_[601-800] #OK

#18943916_[801-900] #OK

#18979629_[901-1000 #OK

#1001-1030 OK

#sacct --format JobID%-20,Submit,Start,End,Cluster,State,ExitCode,User,Group,JobName,QOS,AllocCPUS,NNodes,NTasks,TotalCPU,REQMEM,MaxRSS,TIMELIMIT,Elapsed -j 18274045 | grep CANCELLED

#sacct --format JobID%-20,State -j 18903955 | grep CANCELLED

#Request local storage using the --gres flag in the job submission:
#--gres=nvme:<local_storage_space_per_node>
#The amount of space is given in GB (check maximum sizes from the list above). 
#For example, to request 100 GB of storage, use option --gres=nvme:100. 
#The local storage reservation is on a per node basis.

#Use the environment variable $LOCAL_SCRATCH in your batch job scripts to access the local storage on each node.

#cd $LOCAL_SCRATCH
#cp MOODS_"${lines[$array]}".csv.gz $outfolder 

srun MOODS_artificial_motifs.sh ${SLURM_ARRAY_TASK_ID}

echo "success"

seff $SLURM_JOBID
