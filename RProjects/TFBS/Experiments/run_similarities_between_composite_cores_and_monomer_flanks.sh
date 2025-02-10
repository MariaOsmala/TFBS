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

#Another way
#24120893_1 #OK
#24120907_1-300
#24121110_301-400
#24121117_401-700
#24121564_701-1000
#24121564_869.0       CANCELLED+ 
#24121564_892.0       CANCELLED+ 
#rerun 24123721_869 TIMEOUT
#rerun 24124209_892 TIMEOUT

#24136692_869,892

#24123635_1-23 #1001-1023 #OK
 
#24041919_1-300 #OK
#24042124-301-400 #OK
#24042489_401-600 #OK
#24042695_601-700 #OK
#24042701-701-1000 #OK

#24044315_1-23 #1001-1023 #OK
 

# Load r-env
module load r-env/430


# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2006203/tmp///" >> ~/.Renviron

# Run the R script

#srun apptainer_wrapper exec Rscript --no-save ../code/similarities_between_composite_core_and_overlapping_monomer_flanks.R ${SLURM_ARRAY_TASK_ID}
srun apptainer_wrapper exec Rscript --no-save ../code/similarities_between_composite_core_and_overlapping_monomer_flanks_another_way.R ${SLURM_ARRAY_TASK_ID}

seff $SLURM_JOBID
