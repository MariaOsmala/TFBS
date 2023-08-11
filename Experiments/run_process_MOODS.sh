#!/bin/bash -l
#SBATCH --job-name=process_MOODS
#SBATCH --account=project_2006203
#SBATCH --output=outs/process_MOODS_%j_%a.txt
#SBATCH --error=errs/process_MOODS_%j_%a.txt
#SBATCH --partition=small
#SBATCH --time=3-00:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=150G #3G with database
#SBATCH --array=1-32 ##0-32

#3,4,5,6,11
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
srun apptainer_wrapper exec Rscript --no-save ../RProjects/TFBS/MOODS_results_to_database.R $SLURM_ARRAY_TASK_ID
#srun apptainer_wrapper exec Rscript --no-save ../RProjects/TFBS/MOODS_results_to_database_mouse.R $SLURM_ARRAY_TASK_ID

seff $SLURM_JOBID
