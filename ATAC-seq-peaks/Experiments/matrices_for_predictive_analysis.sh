#!/bin/bash -l
#SBATCH --job-name=pred_matrix
#SBATCH --account=project_2006203
#SBATCH --output=outs/pred_matrix.out
#SBATCH --error=errs/pred_matrix.err
#SBATCH --partition=longrun
#SBATCH --time=5-00:00:00 #1-12:46:56
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G #54.02 GB

#    21298579  


# Load r-env
module load r-env/421

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2006203/tmp///" >> ~/.Renviron

# Run the R script

srun apptainer_wrapper exec Rscript --no-save ../code/matrices_for_predictive_analysis.R

seff $SLURM_JOBID
