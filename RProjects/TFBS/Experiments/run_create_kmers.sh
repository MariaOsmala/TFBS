#!/bin/bash -l
#SBATCH --job-name=create_kmers
#SBATCH --account=project_2006203
#SBATCH --output=outs/create_kmers_%A_%a.txt
#SBATCH --error=errs/create_kmers_%A_%a.txt
#SBATCH --partition=small
#SBATCH --time=06:00:00 #2.5days?
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=350G  #200G
#SBATCH --array=15


#24031139_14-17
#24031351_15-17
#24033438_16-17
# Load r-env

module load r-env/430


# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2006203/tmp///" >> ~/.Renviron

# Run the R script

srun apptainer_wrapper exec Rscript --no-save ../code/create_kmers.R ${SLURM_ARRAY_TASK_ID}

seff $SLURM_JOBID
