#!/bin/bash -l
#SBATCH --job-name=process_logistic_reg
#SBATCH --account=project_2001678
#SBATCH --output=outs/process_logistic_reg.txt
#SBATCH --error=errs/process_logistic_reg.txt
#SBATCH --partition=small
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=200G #3G with database



#21339926_[1-15]

# Load r-env
module load r-env/430

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2006203/tmp///" >> ~/.Renviron

srun apptainer_wrapper exec Rscript --no-save ../code/collect_logistic_regression_results.R
