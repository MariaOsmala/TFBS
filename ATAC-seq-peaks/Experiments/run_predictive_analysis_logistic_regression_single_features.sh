#!/bin/bash -l
#SBATCH --job-name=logistic_reg
#SBATCH --account=project_2001678 
#SBATCH --output=outs/logistic_reg_%j_%a.txt
#SBATCH --error=errs/logistic_reg_%j_%a.txt
#SBATCH --partition=small
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=200G #3G with database
#SBATCH --array=1-2 #1-111


#21450805_[1-2] 

# Load r-env
module load r-env/430

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2006203/tmp///" >> ~/.Renviron


srun apptainer_wrapper exec Rscript --no-save ../code/predictive_analysis_logistic_regression_single_features.R $SLURM_ARRAY_TASK_ID


seff $SLURM_JOBID
