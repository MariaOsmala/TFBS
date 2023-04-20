#!/bin/bash
#SBATCH --job-name=MOODS
#SBATCH --output=outs/MOODS_%A_%a.out
#SBATCH --error=errs/MOODS_%A_%a.err
#SBATCH --account=project_2006203
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=3-00:00:00
#SBATCH --mem-per-cpu=300G #max 382G
#SBATCH --cpus-per-task=1
#SBATCH --array=7


case $SLURM_ARRAY_TASK_ID in
1) ARGS="Nitta2015 Homo_sapiens";; #10 #34 mins 5.2GB
2) ARGS="Jolma2013 Homo_sapiens";; #708
3) ARGS="Jolma2013 Mus_musculus";; #134 X 71 GB, 7 hours
4) ARGS="Jolma2015 Homo_sapiens";; #658
5) ARGS="Yin2017 Homo_sapiens";; #1794
6) ARGS="fromYimeng pfm_composite_new";; #503
7) ARGS="fromYimeng pfm_spacing_new";; #202 #12 hours wasn't enough, needs 100GB
esac



srun MOODS.sh $ARGS

seff $SLURM_JOBID
