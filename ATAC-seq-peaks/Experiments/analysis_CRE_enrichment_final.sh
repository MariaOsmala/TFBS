#!/bin/bash -l
#SBATCH --job-name=CRE_enrichment
#SBATCH --account=project_2006203
#SBATCH --output=outs/CRE_enrichment.out
#SBATCH --error=errs/CRE_enrichment.err
#SBATCH --partition=small
#SBATCH --time=3-00:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=200G


# Load r-env
#module load r-env/421
module load r-env #/430

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2006203/tmp///" >> ~/.Renviron

# Run the R script

#srun apptainer_wrapper exec Rscript --no-save ../code/enrichment_at_CREs_final.R #only representatives
srun apptainer_wrapper exec Rscript --no-save ../code/enrichment_at_CREs_final_all_motifs.R #all 3294 motifs

seff $SLURM_JOBID
