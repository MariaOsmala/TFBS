#!/bin/bash -l
#SBATCH --job-name=parts_CRE_enrichment
#SBATCH --account=project_2006203
#SBATCH --output=outs/parts_CRE_enrichment_%A_%a.txt
#SBATCH --error=errs/parts_CRE_enrichment_%A_%a.txt
#SBATCH --partition=small
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G
#SBATCH --array=301-394 #1:394

#1 OK
#21360459 2-100 #OK
#21360462 101-200 #OK
#21360464 201-300 #OK
# 21360466 301-394 #OK

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
srun apptainer_wrapper exec Rscript --no-save ../code/enrichment_at_CREs_final_all_motifs_parts.R $SLURM_ARRAY_TASK_ID #all 3294 motifs

seff $SLURM_JOBID
