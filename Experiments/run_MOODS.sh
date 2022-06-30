#!/bin/bash
#SBATCH --job-name=MOODS
#SBATCH --output=MOODS.out
#SBATCH --error=MOODS.err
#SBATCH --account=project_2006203
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --partition=small
##SBATCH --mail-type=BEGIN #uncomment to enable mail



export PATH="/projappl/project_2006203/softwares/conda_envs/MOODS/bin:$PATH"



S=/projappl/project_2006203/Genomes/Homo_sapiens/chr_sequences.fa

n=300000 #the number of best hits

outfile=../Results/MOODS/TCF7.csv

srun moods-dna.py -S ../PWMs/Yin2017/pwms/Homo_sapiens/TCF7_eDBD_HMG_HT-SELEX_TCGGAC40NCCA_KS_ASATCAAAS_1_4_NO.pfm \
../PWMs/Yin2017/pwms/Homo_sapiens/TCF7_eDBD_HMG_Methyl-HT-SELEX_TTGTCT40NGCT_KS_ASATCAAAS_1_4_NO.pfm \
../PWMs/Yin2017/pwms/Homo_sapiens/TCF7_eDBD_HMG_HT-SELEX_TCGGAC40NCCA_KS_WCATCGRGRCGCTGW_2_4_NO.pfm \
../PWMs/Yin2017/pwms/Homo_sapiens/TCF7_eDBD_HMG_Methyl-HT-SELEX_TTGTCT40NGCT_KS_WCATCGRGRCGCTGW_2_4_NO.pfm \
 -s $S -B $n -o $outfile 
 
seff $SLURM_JOBID
