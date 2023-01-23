#!/bin/bash -l
#SBATCH --job-name=r_multicore
#SBATCH --account=project_2006203
#SBATCH --output=out_%j.txt
#SBATCH --error=err_%j.txt
#SBATCH --partition=test
#SBATCH --time=00:05:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=1000

# Load r-env
module load r-env

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/<project>" >> ~/.Renviron

# Run the R script
srun apptainer_wrapper exec Rscript --no-save analysis_SQLite_CRMs.R

seff $SLURM_JOBID
