#!/bin/bash -l
#SBATCH --job-name=logistic_reg
#SBATCH --account=project_2006203
#SBATCH --output=outs/logistic_reg_%j_%a.txt
#SBATCH --error=errs/logistic_reg_%j_%a.txt
#SBATCH --partition=small
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G #3G with database
#SBATCH --array=2-111 #1-111


# Load r-env
module load r-env/421

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2006203/tmp///" >> ~/.Renviron


srun apptainer_wrapper exec Rscript --no-save ../code/predictive_analysis_logistic_regression.R $SLURM_ARRAY_TASK_ID


seff $SLURM_JOBID
