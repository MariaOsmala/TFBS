#!/bin/bash
#SBATCH --job-name=GRangesList_MOODS
#SBATCH --account=project_2007567
#SBATCH --output=outs/GRangesList_MOODS_%A.out
#SBATCH --error=errs/GRangesList_MOODS_%A.err
#SBATCH --partition=small
#SBATCH --time=3-00:00:00 #2.5days?
#SBATCH --mem-per-cpu=250G 

#

#28078258 hg38.analysisSet

# Load r-env
module load r-env/430

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2006203/tmp///" >> ~/.Renviron

#Old
#input=Results/MOODS_human_final_version2.2_correct_processed/MOODS_RDS
#output=/ATAC-seq-peaks/RData/top_motif_matches_human_final_version2.2.Rds

#hg38.analysisSet

input=Results/MOODS_hg38.analysisSset_processed/MOODS_RDS/ 
output=Results/RData/hg38.analysisSset.Rds

#input=Results/MOODS_T2T-CHM13_processed/MOODS_RDS/ 
#output=Results/RData/T2T-CHM13.Rds

echo "test"

# Run the R script
srun apptainer_wrapper exec Rscript --no-save ../../code/MOODS/MOODS_results_to_single_GRangesList.R $input $output
#srun apptainer_wrapper exec Rscript --no-save ../code/MOODS_results_to_single_GRangesList_from_bedFiles.R

seff $SLURM_JOBID
