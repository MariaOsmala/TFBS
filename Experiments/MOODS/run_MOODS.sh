#!/bin/bash
#SBATCH --job-name=MOODS
#SBATCH --output=outs/MOODS_%A_%a.out
#SBATCH --error=errs/MOODS_%A_%a.err
#SBATCH --account=project_2006203
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=03:00:00
#SBATCH --mem-per-cpu=80G #
#SBATCH --cpus-per-task=1
#SBATCH --gres=nvme:50
#SBATCH --array=0-393 # #0-393  

# 27928718-0-393 #T2T
# Max memory (KB):9,741,227
# Max elapsed time (seconds): 9757 2.7h

# sacct --format=JobID,MaxRSS,Elapsed -j  27928718 --parsable2 | \
# awk -F'|' '
#   function to_sec(t) {
#     split(t, a, ":")
#     return a[1]*3600 + a[2]*60 + a[3]
#   }
# 
#   $2 ~ /^[0-9]+K$/ {
#     mem = substr($2, 1, length($2)-1)
#     if (mem > max_mem) max_mem = mem
#   }
# 
#   $3 ~ /^[0-9]{2}:[0-9]{2}:[0-9]{2}$/ {
#     sec = to_sec($3)
#     if (sec > max_time) max_time = sec
#   }
# 
#   END {
#     printf "Max memory (KB): %d\n", max_mem
#     print "Max elapsed time (seconds):", max_time
#   }
#   '



# Max memory (KB): 76 477 467 
# Max elapsed time (seconds): 8227 = 2.3h



#which runs fail 
#sacct --format JobID%-20,Submit,Cluster,State,ExitCode,User,Group,JobName,QOS,AllocCPUS,NNodes,NTasks,TotalCPU,REQMEM,MaxRSS,TIMELIMIT,Elapsed -j 27682712 | grep CANCELLED

#sacct --format JobID%-20,State -j 27682471 | grep CANCELLED

#Request local storage using the --gres flag in the job submission:
#--gres=nvme:<local_storage_space_per_node>
#The amount of space is given in GB (check maximum sizes from the list above). 
#For example, to request 100 GB of storage, use option --gres=nvme:100. 
#The local storage reservation is on a per node basis.

#Use the environment variable $LOCAL_SCRATCH in your batch job scripts to access the local storage on each node.

#T2T
genome=/projappl/project_2007567/Genomes/T2T-CHM13/chm13v2.0_maskedY_rCRS_shorter_chr_names.fa
outfolder=/scratch/project_2006203/TFBS/Results/MOODS_T2T-CHM13/ 
pwm_filenames="../../Data/SELEX-motif-collection/filenames.csv"

#hg38 analysis set
#genome=/projappl/project_2007567/Genomes/GRCh38.p14_31102023/hg38.analysisSet.fa
#outfolder=/scratch/project_2006203/TFBS/Results/MOODS_hg38.analysisSset/
#pwm_filenames="../../Data/SELEX-motif-collection/filenames.csv"

srun MOODS.sh ${SLURM_ARRAY_TASK_ID} $genome $outfolder $pwm_filenames 


seff $SLURM_JOBID
