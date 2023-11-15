#!/bin/bash -l
#SBATCH --job-name=Kvon
#SBATCH --account=project_2006203
#SBATCH --output=outs/Kvon.out
#SBATCH --error=errs/Kvon.err
#SBATCH --partition=small
#SBATCH --time=03-00:00:00 #2.5days?
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G  #200G

# Load r-env
module load r-env/430

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2006203/tmp///" >> ~/.Renviron

# Run the R script
srun apptainer_wrapper exec Rscript --no-save ../motif_matches_at_Kvon_limb_enhancers_mouse.R

seff $SLURM_JOBID
