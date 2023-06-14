#!/bin/bash

#!/bin/bash
array=$1 #this varies between 0 and 399

files=($(awk -F '"' '{print $4}' rerunMOODS_mouse_lower_threshold_3.txt)) #1001 736

#echo ${files[0]}
#ELK1_TBX21_TACTCT40NAAC_YYI_NAGGTGTTACTTCCGGYN_m2_c3b0u_short_composite
#${files[1000]}
#IRF4_HT-SELEX_TCAAGG20NCG_AD_NCGAAACCGAAACYN_1_3_YES


export PATH="/projappl/project_2006203/softwares/conda_envs/MOODS/bin:$PATH"


S=/projappl/project_2006203/Genomes/Mus_musculus/GRCm39_mm39/chr_sequences.fa #This is repeat masked (hard), contains only chromsomes chr1-22, X, Y
#mkdir /scratch/project_2006203/TFBS/Results/MOODS_mouse_mm39_4
#mkdir /scratch/project_2006203/TFBS/Results/MOODS_mouse_mm39_3
mkdir /scratch/project_2006203/TFBS/Results/MOODS_mouse_mm39_2

n=300000 #the number of best hits

pwms=($(ls ../PWMs/*/pwms_space/*/*.pfm)) #3982+15=3997

#which pwms are in files
filenames=(${pwms[@]##*/})
filenames2=(${filenames[@]%%.pfm})

#($(awk 'BEGIN{RS = FS} NR == FNR {files[$1] = 1; next} $1 in files' \
#    < (echo "${files[*]}") < (echo "${filenames2[*]}")))

declare -A indices
declare -A pwms_new

i=0
for file in "${files[@]}"
do
   #echo $i
   #echo $file
   #echo ${filenames2[@]/$file//} | cut -d/ -f1 | wc -w | tr -d ' '
   #tmp=$(echo ${filenames2[@]/$file//} | cut -d/ -f1 | wc -w | tr -d ' ')
   #echo $tmp
   indices[$i]=$(echo ${filenames2[@]/$file//} | cut -d/ -f1 | wc -w | tr -d ' ')
   pwms_new[$i]=${pwms[${indices[$i]}]}
   i=$((i+1))
done


nro_pwms=${#pwms_new[@]}


start_ind=$(($array*10)) #100
end_ind=$((($array+1)*10 -1)) #100
length=10 #100

if [[ $end_ind -gt $(($nro_pwms-1)) ]]
then
     echo $end_ind is greater than $(($nro_pwms-1))
     length=$(($nro_pwms-$start_ind+1))
fi

export PATH="/projappl/project_2006203/softwares/conda_envs/MOODS/bin:$PATH"

moods-dna.py -m ${pwms_new[@]:$start_ind:$length} --threshold 2 -s $S | gzip > $LOCAL_SCRATCH"/MOODS_"$array".csv.gz"

cd $LOCAL_SCRATCH
#cp MOODS_"$array".csv.gz /scratch/project_2006203/TFBS/Results/MOODS_mouse_mm39_4/
#cp MOODS_"$array".csv.gz /scratch/project_2006203/TFBS/Results/MOODS_mouse_mm39_3/
cp MOODS_"$array".csv.gz /scratch/project_2006203/TFBS/Results/MOODS_mouse_mm39_2/


