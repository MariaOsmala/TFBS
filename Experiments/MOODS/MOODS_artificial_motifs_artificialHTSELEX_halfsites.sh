#!/bin/bash
array=$1 #

export PATH="/projappl/project_2006203/softwares/conda_envs/MOODS/bin:$PATH"


S=/projappl/project_2006203/Genomes/Homo_sapiens/chr_sequences.fa #This is repeat masked, contains only chromsomes chr1-22, X, Y

n=300000 #the number of best hits

outfolder=/scratch/project_2006203/TFBS/Results/MOODS_artificial_human_artificialHTSELEX/

mkdir $outfolder

#Read in motif IDs
readarray -t lines < Data/SELEX-motif-collection/artificial-half-site-motifs/motifnames.csv #7 Need to create this!

start_ind=0
end_ind=7

pwms_all=()

for (( index=$start_ind; index<$end_ind; index++ )); do


  pwms=($(ls Data/SELEX-motif-collection/artificial-half-site-motifs/pfms_space/${lines[$index]}*.pfm))
  echo ${lines[$index]}
  echo ${#pwms[@]}
  
  true_motif=${lines[$index]}

  pwms_all+=( "${pwms[@]}"  )

done



moods-dna.py -m ${pwms_all[@]} --threshold 2 -s $S | gzip > $LOCAL_SCRATCH"/MOODS_"$array".csv.gz" 


cd $LOCAL_SCRATCH
cp MOODS_"$array".csv.gz $outfolder 

