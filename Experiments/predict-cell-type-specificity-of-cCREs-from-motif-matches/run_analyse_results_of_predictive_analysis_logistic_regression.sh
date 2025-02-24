#!/bin/bash -l
#SBATCH --job-name=figures_logistic_reg
#SBATCH --account=project_2001678
#SBATCH --output=outs/figures_logistic_reg_%j_%a.txt
#SBATCH --error=errs/figures_logistic_reg_%j_%a.txt
#SBATCH --partition=small
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G #3G with database
#SBATCH --array=1 #-15 #1-15


#21339926_[1-15]

# Load r-env
module load r-env/421

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2006203/tmp///" >> ~/.Renviron

#srun apptainer_wrapper exec Rscript --no-save ../code/collect_logistic_regression_results.R
srun apptainer_wrapper exec Rscript --no-save ../code/analyse_results_of_predictive_analysis_logistic_regression.R $SLURM_ARRAY_TASK_ID



seff $SLURM_JOBID
