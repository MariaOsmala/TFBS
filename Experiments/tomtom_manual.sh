#!/bin/bash
#SBATCH --job-name=tomtom
#SBATCH --output=outs/tomtom.out
#SBATCH --error=errs/tomtom.err
#SBATCH --account=project_2006203
#SBATCH --time=2-00:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --gres=nvme:20 
#SBATCH --partition=small
##SBATCH --mail-type=BEGIN #uncomment to enable mail


##SBATCH --gres=nvme:50 This is in GB


export PATH="/projappl/project_2006203/softwares/conda_envs/motif-clustering-Vierstra/bin:$PATH"

cp ../PWMs_final/all.meme $LOCAL_SCRATCH

tomtom -dist kullback -motif-pseudo 0.1 -text -min-overlap 1 $LOCAL_SCRATCH/all.meme $LOCAL_SCRATCH/all.meme > $LOCAL_SCRATCH/tomtom.all.txt

cd $LOCAL_SCRATCH
cp tomtom.all.txt /scratch/project_2006203/Results/tomtom_final/tomtom/

