#!/bin/bash -l
#SBATCH --job-name=motif_logos
#SBATCH --account=project_2006203
#SBATCH --output=outs/motif_logos.out
#SBATCH --error=errs/motif_logos.err
#SBATCH --partition=small
#SBATCH --time=01-00:00:00 #2.5days?
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G  #200G

# Load r-env
module load r-env

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2006203/tmp///" >> ~/.Renviron

# Run the R script
srun apptainer_wrapper exec Rscript --no-save ../code/motif_logos_better.R

seff $SLURM_JOBID
