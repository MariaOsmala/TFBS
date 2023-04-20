#!/bin/bash -l
#SBATCH --job-name=CRE_enrichment_Teemu
#SBATCH --account=project_2006203
#SBATCH --output=outs/CRE_enrichment_Teemu.out
#SBATCH --error=errs/CRE_enrichment_Teemu.err
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

srun apptainer_wrapper exec Rscript --no-save ../code/enrichment_at_CREs_Teemu.R

seff $SLURM_JOBID
