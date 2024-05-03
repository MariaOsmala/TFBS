#!/bin/bash -l
#SBATCH --job-name=process_MOODS
#SBATCH --account=project_2006203
#SBATCH --output=outs/process_MOODS_%A_%a.txt
#SBATCH --error=errs/process_MOODS_%A_%a.txt
#SBATCH --partition=small
#SBATCH --time=06:00:00 #08:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=200G #300G
#SBATCH --gres=nvme:50
#SBATCH --array=786

# 21395129_[786]

#for true motifs #0:39 1 day, 150G
#2001678
#for artificial motifs: 0-786 #

#21346873_0,100,158

#  21347748 1-99
# rerun 21350283 100G
# 1,2,5,6,11,12,13,17,20,21,25,27-29,32,33,35,39-46,50,51,54,56-60,65,67-69,71,77-99 

#still memory runs out 21352368 OK
#25,32,39,56,69,82,83,84,87,89,92,95,97 

# 101-300 21352457 OK

#21352374 301-600 OK

#21347770 601-786

#rerun 21352434 OK
#601-620,622,625-629,631-633,636,642,644-647,649,651,653-655,660,666-671,673-680,682-684,687,688,690,692,693,696,697,699,703,704,706,707,709,711,712,714,715,719-722,724-728,730-733,736,739,740,742,743,746,748,750,755,758,759,760,762-765,769,771,777,784,785 




#sacct --format JobID%-20,State -j 21347755
# MOODS_results_to_database.R $SLURM_ARRAY_TASK_ID=0-40
# srun apptainer_wrapper exec Rscript --no-save ../RProjects/TFBS/MOODS_artificial_results_to_database.R $SLURM_ARRAY_TASK_ID=0-1030



##SBATCH --gres=nvme:50
#with 1 CPU required 1-04:33:39
#Memory Utilized: 55.83 GB 55.83% of 100.00 GB

#a-put tar command, you can see what is inside the compressed file a-find, where is your tar file

#sacct -o jobname,jobid,reqmem,maxrss,timelimit,elapsed,state --units=[KMGTP] -j <jobid> make some alias on this
#Add the current working directory to $PATH:
#export PATH=$PWD:$PATH

#To add paths automatically, you can add the export command to your $HOME/.bashrc file. Instead of $PWD, use the full path:
#export PATH=/projappl/<project>/$USER/gcta-1.94.1-linux-kernel-3-x86_64:$PATH   # replace <project> with your CSC project, e.g. project_2001234

# Load r-env
module load r-env/421

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2006203/tmp///" >> ~/.Renviron

# Run the R script
#srun apptainer_wrapper exec Rscript --no-save ../RProjects/TFBS/collect_MOODS_results.R $SLURM_ARRAY_TASK_ID
#srun apptainer_wrapper exec Rscript --no-save ../RProjects/TFBS/MOODS_results_to_database.R $SLURM_ARRAY_TASK_ID
srun apptainer_wrapper exec Rscript --no-save ../RProjects/TFBS/MOODS_artificial_results_to_database.R $SLURM_ARRAY_TASK_ID
#srun apptainer_wrapper exec Rscript --no-save ../RProjects/TFBS/MOODS_results_to_database_mouse.R $SLURM_ARRAY_TASK_ID

seff $SLURM_JOBID
