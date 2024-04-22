#!/bin/bash -l
#SBATCH --job-name=motifs_process
#SBATCH --account=project_2006203
#SBATCH --output=outs/motifs_process.out
#SBATCH --error=errs/motifs_process.err
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G

#21237593

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


srun apptainer_wrapper exec Rscript --no-save ../code/process_and_analyse_top_MOODS_hits.R 

seff $SLURM_JOBID
