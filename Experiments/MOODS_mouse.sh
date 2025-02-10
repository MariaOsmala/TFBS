#!/bin/bash

#!/bin/bash
array=$1 #this varies between 0 and 393



export PATH="/projappl/project_2006203/softwares/conda_envs/MOODS/bin:$PATH"


S=/projappl/project_2006203/Genomes/Mus_musculus/GRCm39_mm39/chr_sequences.fa #This is repeat masked (hard), contains only chromsomes chr1-22, X, Y
#S=/projappl/project_2007567/Genomes/mm39/mm39.fa non-repeat masked

#Old #outfolder=/scratch/project_2006203/TFBS/Results/MOODS_mouse_mm39_final/
#outfolder=/scratch/project_2006203/TFBS/Results/MOODS_mouse_mm39_final_whole_genome/ #non-repeat masked
outfolder=/scratch/project_2006203/TFBS/Results/MOODS_mouse_mm39_version2.2/
mkdir $outfolder


n=300000 #the number of best hits

pwms=( $(cut -d ',' -f1 ../PWMs_final_version2.2/filenames.csv ) )


#echo "${pwms[3]##*/}"

#echo "${#pwms[@]}"
#echo "${pwms[@]}"

nro_pwms=${#pwms[@]}


start_ind=$(($array*10)) #100
end_ind=$((($array+1)*10 -1)) #100
length=10 #100

if [[ $end_ind -gt $(($nro_pwms-1)) ]]
then
     echo $end_ind is greater than $(($nro_pwms-1))
     length=$(($nro_pwms-$start_ind+1))
fi



moods-dna.py -m ${pwms[@]:$start_ind:$length} --threshold 2 -s $S | gzip > $LOCAL_SCRATCH"/MOODS_"$array".csv.gz"

cd $LOCAL_SCRATCH
cp MOODS_"$array".csv.gz $outfolder







