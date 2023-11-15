#!/bin/bash -l
#SBATCH --job-name=process_MOODS
#SBATCH --account=project_2007567
#SBATCH --output=outs/process_MOODS_%A_%a.txt
#SBATCH --error=errs/process_MOODS_%A_%a.txt
#SBATCH --partition=small
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=300G #3G with database
#SBATCH --gres=nvme:100
#SBATCH --array=616,684,726 




# 19184582_[1-300] #OK

#301-500: OK

#501-600: OK

#19209909_[601-800]

#FAILED: 616,684,726 #19386445

#OK 642,685,



#19235452_[901-1000 OK!

#1001-1031: 19322725_[0-30] OK!

#0-32 human and mouse true motifs

#0-1030 artificial representative motifs

#0-15: 18913194 #DONE

#16-32: 18919036 Memory runs out, rerun 17-19,21-24,26-27,29,31

#rerun 18930020

#Figure out which artificial results are still missing and their job array numbers




#0-10
#3,4,5,6,11 #SBATCH --gres=nvme:50
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
