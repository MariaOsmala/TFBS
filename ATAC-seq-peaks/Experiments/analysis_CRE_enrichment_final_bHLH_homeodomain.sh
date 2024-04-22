#!/bin/bash -l
#SBATCH --job-name=bHLH_homeo
#SBATCH --account=project_2007567
#SBATCH --output=outs/bHLH_homeo.out
#SBATCH --error=errs/bHLH_homeo.out
#SBATCH --partition=small
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=150G


# 21291268
#21291311

# Load r-env
#module load r-env/421
module load r-env/430

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2006203/tmp///" >> ~/.Renviron

# Run the R script

#srun apptainer_wrapper exec Rscript --no-save ../code/enrichment_at_CREs_final.R #only representatives
srun apptainer_wrapper exec Rscript --no-save ../code/enrichment_at_CREs_final_bHLH_homeodomain_only.R 

seff $SLURM_JOBID
