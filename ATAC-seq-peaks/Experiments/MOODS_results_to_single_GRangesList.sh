#!/bin/bash -l
#SBATCH --job-name=GRangesList_MOODS
#SBATCH --account=project_2006203
#SBATCH --output=GRangesList_MOODS.out
#SBATCH --error=GRangesList_MOODS.err
#SBATCH --partition=small
#SBATCH --time=02-12:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=200G


# Load r-env
module load r-env

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2006203/tmp///" >> ~/.Renviron

# Run the R script
srun apptainer_wrapper exec Rscript --no-save ../code/MOODS_results_to_single_GRangesList.R

seff $SLURM_JOBID
