#!/bin/bash -l
#SBATCH --job-name=process_MOODS
#SBATCH --account=project_2007567
#SBATCH --output=outs/process_MOODS_%A_%a.txt
#SBATCH --error=errs/process_MOODS_%A_%a.txt
#SBATCH --partition=small
#SBATCH --time=06:00:00 #08:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=150G #300G
#SBATCH --gres=nvme:50
#SBATCH --array=1-39 #0:39


#for true motifs #0:39 1 day, 150G

#for artificial motifs: 0-786 #

# Load r-env
module load r-env/421

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2006203/tmp///" >> ~/.Renviron

code_path=

#To process true motif matches
srun apptainer_wrapper exec Rscript --no-save ../RProjects/TFBS/MOODS_results_to_database.R $SLURM_ARRAY_TASK_ID

#To process artificial half-site motif matches
#srun apptainer_wrapper exec Rscript --no-save ../RProjects/TFBS/MOODS_artificialHTSelex_results_to_database.R $SLURM_ARRAY_TASK_ID

#To process scrambled motif matches
#srun apptainer_wrapper exec Rscript --no-save ../RProjects/TFBS/MOODS_artificial_results_to_database.R $SLURM_ARRAY_TASK_ID

#To process scrambled artificial half-site motif matches
#srun apptainer_wrapper exec Rscript --no-save ../RProjects/TFBS/MOODS_artificial_artificialHTSelex_results_to_database.R $SLURM_ARRAY_TASK_ID




seff $SLURM_JOBID
