#!/bin/bash

#!/bin/bash
array=$1 


export PATH="/projappl/project_2006203/softwares/conda_envs/MOODS/bin:$PATH"


S=/projappl/project_2006203/Genomes/Homo_sapiens/chr_sequences.fa #This is repeat masked, contains only chromsomes chr1-22, X, Y

n=300000 #the number of best hits


outfolder=/scratch/project_2006203/TFBS/Results/MOODS_human_final_version2.2_correct/ 
mkdir $outfolder


pwms=( $(cut -d ',' -f1 Data/SELEX-motif-collection/filenames.csv ) )

nro_pwms=${#pwms[@]} #3933


start_ind=$(($array*10)) #100
end_ind=$((($array+1)*10 -1)) #100
length=10 #100


if [[ $end_ind -gt $(($nro_pwms-1)) ]] 
then
     echo $end_ind is greater than $(($nro_pwms-1))
     length=$(($nro_pwms-$start_ind+1))
fi

export PATH="/projappl/project_2006203/softwares/conda_envs/MOODS/bin:$PATH"

#threshold selection (exactly one required):
#  -p p, --p-value p     compute threshold from p-value
#  -t T, --threshold T   use specified absolute threshold
#  -B n, --best-hits n   return at least the specified amount of best matches


#-m to give count/frequency matrices, will be converted to PWMs


#--bg pA pC pG pT      background distribution for computing thresholds from
#                        p-value with --batch (default is 0.25 for all alleles)
#--ps p                total pseudocount added to each matrix column in log-
#                        odds conversion (default = 0.01)

#  --log-base x          logarithm base for log-odds conversion (default
#                       natural logarithm)

#  --lo-bg pA pC pG pT   background distribution for log-odds conversion
#                        (default is 0.25 for all alleles)

moods-dna.py -m ${pwms[@]:$start_ind:$length} --threshold 2 -s $S | gzip > $LOCAL_SCRATCH"/MOODS_"$array".csv.gz" 

cd $LOCAL_SCRATCH
cp MOODS_"$array".csv.gz $outfolder 

