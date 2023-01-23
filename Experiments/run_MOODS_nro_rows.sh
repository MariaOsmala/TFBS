#!/bin/bash
#SBATCH --job-name=rows
#SBATCH --output=rows.out
#SBATCH --error=rows.err
#SBATCH --account=project_2006203
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --partition=small
##SBATCH --mail-type=BEGIN #uncomment to enable mail



outfile=../Results/MOODS/nro_rows.txt

unpigz -c MOODS_0.csv.gz | wc -l
unpigz -c 1.csv.gz | wc -l
unpigz -c 10.csv.gz | wc -l
unpigz -c 11.csv.gz | wc -l
unpigz -c 12.csv.gz | wc -l
unpigz -c 13.csv.gz | wc -l
unpigz -c 14.csv.gz | wc -l
unpigz -c 15.csv.gz | wc -l
unpigz -c 16.csv.gz | wc -l
unpigz -c 17.csv.gz | wc -l
unpigz -c 18.csv.gz | wc -l
unpigz -c 19.csv.gz | wc -l
unpigz -c 2.csv.gz | wc -l
unpigz -c 20.csv.gz | wc -l
unpigz -c 21.csv.gz | wc -l
unpigz -c 22.csv.gz | wc -l
unpigz -c 23.csv.gz | wc -l
unpigz -c 24.csv.gz | wc -l
unpigz -c 25.csv.gz | wc -l
unpigz -c 26.csv.gz | wc -l
unpigz -c 27.csv.gz | wc -l
unpigz -c 28.csv.gz | wc -l
unpigz -c 29.csv.gz | wc -l
unpigz -c 3.csv.gz | wc -l
unpigz -c 30.csv.gz | wc -l
unpigz -c 31.csv.gz | wc -l
unpigz -c 32.csv.gz | wc -l
unpigz -c 33.csv.gz | wc -l
unpigz -c 34.csv.gz | wc -l
unpigz -c 35.csv.gz | wc -l
unpigz -c 36.csv.gz | wc -l
unpigz -c 37.csv.gz | wc -l
unpigz -c 38.csv.gz | wc -l
unpigz -c 39.csv.gz | wc -l
unpigz -c 4.csv.gz | wc -l
unpigz -c 5.csv.gz | wc -l
unpigz -c 6.csv.gz | wc -l
unpigz -c 7.csv.gz | wc -l
unpigz -c 8.csv.gz | wc -l
unpigz -c 9.csv.gz | wc -l




seff $SLURM_JOBID
