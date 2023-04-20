#!/bin/bash -l
#SBATCH --job-name=process_MOODS
#SBATCH --account=project_2006203
#SBATCH --output=outs/process_MOODS_%j_%a.txt
#SBATCH --error=errs/process_MOODS_%j_%a.txt
#SBATCH --partition=small
#SBATCH --time=02-00:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=10G #3G with database
#SBATCH --array=0 # #0-39

# Load r-env
module load r-env

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2006203/tmp///" >> ~/.Renviron

# Run the R script
#srun apptainer_wrapper exec Rscript --no-save ../RProjects/TFBS/collect_MOODS_results.R $SLURM_ARRAY_TASK_ID
srun apptainer_wrapper exec Rscript --no-save ../RProjects/TFBS/MOODS_results_to_database.R $SLURM_ARRAY_TASK_ID

seff $SLURM_JOBID
