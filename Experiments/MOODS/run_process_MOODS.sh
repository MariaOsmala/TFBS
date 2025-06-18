#!/bin/bash
#SBATCH --job-name=process_MOODS
#SBATCH --account=project_2006203
#SBATCH --output=outs/process_MOODS_%A_%a.txt
#SBATCH --error=errs/process_MOODS_%A_%a.txt
#SBATCH --partition=small
#SBATCH --time=3-00:00:00 #08:00:00
#SBATCH --mem-per-cpu=250G #300G
#SBATCH --gres=nvme:50
#SBATCH --array=0-39 #0:39

#sacct --format JobID%-20,State -j 27787360 | grep CANCELLED

# hg38: 0 ok 
# 27936410_[1-39] 
#Max memory (KB): 97,282,418
#Max elapsed time (seconds): 83526=23 h

#T2T: 27950920_[0-39] OK
#Max memory (KB): 95765430
#Max elapsed time (seconds): 78830

# sacct --format=JobID,MaxRSS,Elapsed -j 27950920 --parsable2 | \
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


#for true motifs #0:39 1 day, 150G

#for artificial motifs: 0-786 #

# Load r-env
module load r-env/430

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2006203/tmp///" >> ~/.Renviron

results_path=/scratch/project_2006203/TFBS/Results/MOODS_T2T-CHM13_processed/ 
MOODS_path=/scratch/project_2006203/TFBS/Results/MOODS_T2T-CHM13/ 

#results_path=/scratch/project_2006203/TFBS/Results/MOODS_hg38.analysisSset_processed/ 
#MOODS_path=/scratch/project_2006203/TFBS/Results/MOODS_hg38.analysisSset/

echo $MOODS_path
echo $results_path

#To process true motif matches
#srun apptainer_wrapper exec Rscript --no-save ../../code/MOODS/MOODS_results_to_database.R $SLURM_ARRAY_TASK_ID $MOODS_path $results_path #repeat-masted
srun apptainer_wrapper exec Rscript --no-save ../../code/MOODS/MOODS_results_to_database_T2T-CHM13.R $SLURM_ARRAY_TASK_ID $MOODS_path $results_path 
#srun apptainer_wrapper exec Rscript --no-save ../../code/MOODS/MOODS_results_to_database_hg38.analysisSet.R $SLURM_ARRAY_TASK_ID $MOODS_path $results_path 

#To process artificial half-site motif matches
#srun apptainer_wrapper exec Rscript --no-save ../RProjects/TFBS/MOODS_artificialHTSelex_results_to_database.R $SLURM_ARRAY_TASK_ID

#To process scrambled motif matches
#srun apptainer_wrapper exec Rscript --no-save ../RProjects/TFBS/MOODS_artificial_results_to_database.R $SLURM_ARRAY_TASK_ID

#To process scrambled artificial half-site motif matches
#srun apptainer_wrapper exec Rscript --no-save ../RProjects/TFBS/MOODS_artificial_artificialHTSelex_results_to_database.R $SLURM_ARRAY_TASK_ID




seff $SLURM_JOBID
