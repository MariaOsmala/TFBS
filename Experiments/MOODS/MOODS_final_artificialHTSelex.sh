#!/bin/bash

#!/bin/bash
array=$1 #When 66 motifs varies between 0 and 6


export PATH="/projappl/project_2006203/softwares/conda_envs/MOODS/bin:$PATH"


S=/projappl/project_2006203/Genomes/Homo_sapiens/chr_sequences.fa #This is repeat masked, contains only chromsomes chr1-22, X, Y

n=300000 #the number of best hits


outfolder=/scratch/project_2006203/TFBS/Results/MOODS_human_final_artificialHTSelex/ 
mkdir $outfolder
#outfile=../Results/MOODS/MOODS_"$array".csv"

directory=/scratch/project_2006203/TFBS/PWMs_final_version2.2/fromYimeng/half-site-motifs-20240426/pwms_space

pwms=("$directory"/*)

export PATH="/projappl/project_2006203/softwares/conda_envs/MOODS/bin:$PATH"

moods-dna.py -m ${pwms[@]:0:7} --threshold 2 -s $S | gzip > $outfolder/MOODS.csv.gz


