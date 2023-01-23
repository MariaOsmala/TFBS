#!/bin/bash -l
#SBATCH --job-name=analysis_RDS
#SBATCH --account=project_2006203
#SBATCH --output=outs/analysis_RDS.out
#SBATCH --error=errs/analysis_RDS.err
#SBATCH --partition=small
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G


# Load r-env
module load r-env

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2006203/tmp///" >> ~/.Renviron

# Run the R script
#srun apptainer_wrapper exec Rscript --no-save ../code/analysis_SQLite.R
srun apptainer_wrapper exec Rscript --no-save ../code/analysis_Vierstra_RDS.R

seff $SLURM_JOBID
