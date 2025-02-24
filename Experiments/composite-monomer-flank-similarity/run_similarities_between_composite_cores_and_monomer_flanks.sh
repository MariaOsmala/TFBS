#!/bin/bash -l
#SBATCH --job-name=core_flanks
#SBATCH --account=project_2007567
#SBATCH --output=outs/core_flanks_%A_%a.txt
#SBATCH --error=errs/core_flanks_%A_%a.txt
#SBATCH --partition=hugemem
#SBATCH --time=1-00:00:00 #2.5days?
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=400G  #50G mostly enough
#SBATCH --array=892,869


# Load r-env
module load r-env/430


# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2006203/tmp///" >> ~/.Renviron

# Run the R script

srun apptainer_wrapper exec Rscript --no-save ../code/similarities_between_composite_core_and_overlapping_monomer_flanks.R ${SLURM_ARRAY_TASK_ID}


seff $SLURM_JOBID
