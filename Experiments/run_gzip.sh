#!/bin/bash
#SBATCH --job-name=gzip
#SBATCH --output=outs/gzip.out
#SBATCH --error=errs/gzip.err
#SBATCH --account=project_2006203
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=10G #max 382G
#SBATCH --cpus-per-task=1


#DONE 0-4
#running 5, 35-39


cd /scratch/project_2006203/TFBS/Results/MOODS

gzip MOODS_36.csv
gzip MOODS_37.csv


cd /scratch/project_2006472/MOODS

gzip MOODS_0.csv
gzip MOODS_1.csv
gzip MOODS_2.csv
gzip MOODS_3.csv
gzip MOODS_4.csv
