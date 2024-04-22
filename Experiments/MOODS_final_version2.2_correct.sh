#!/bin/bash

#!/bin/bash
array=$1 #When 66 motifs varies between 0 and 6


export PATH="/projappl/project_2006203/softwares/conda_envs/MOODS/bin:$PATH"


S=/projappl/project_2006203/Genomes/Homo_sapiens/chr_sequences.fa #This is repeat masked, contains only chromsomes chr1-22, X, Y

n=300000 #the number of best hits

#outfolder=/scratch/project_2006203/TFBS/Results/MOODS_human_final/ #version1
outfolder=/scratch/project_2006203/TFBS/Results/MOODS_human_final_version2.2_correct/ #version2
mkdir $outfolder
#outfile=../Results/MOODS/MOODS_"$array".csv"

pwms=( $(cut -d ',' -f1 ../PWMs_final_version2.2/filenames.csv ) )

nro_pwms=${#pwms[@]} #3933


start_ind=$(($array*10)) #100
end_ind=$((($array+1)*10 -1)) #100
length=10 #100


if [[ $end_ind -gt $(($nro_pwms-1)) ]] 
then
     echo $end_ind is greater than $(($nro_pwms-1))
     length=$(($nro_pwms-$start_ind+1))
fi

#threshold selection (exactly one required):
#  -p p, --p-value p     compute threshold from p-value
#  -t T, --threshold T   use specified absolute threshold
#  -B n, --best-hits n   return at least the specified amount of best matches

# search and model behaviour (optional):
#  --no-snps             ignore IUPAC symbols coding multiple nucleotides
#  --threshold-precision x
#                        specify the precision used for computing the
#                        thresholds from p-values (default = 2000.0)

#Write this directly to a database

#--p-value 1e-4 --lo-bg 0.2977 0.2023 0.2023 0.2977

#moods-dna.py -m ${pwms[@]:$start_ind:$length}  -s $S -B $n -o $outfile #${array[@]:START:LENGTH}

#removed -B $n

#by default, MOODS assume that the threshold is given by a p-value x, 
#and the actual threshold T is chosen so that the probability that 
#the background distribution π generates 
#a sequence u of length m with score W_L(u)>=T is p.
# p-value 0.0001

#-m to give count/frequency matrices, will be converted to PWMs
export PATH="/projappl/project_2006203/softwares/conda_envs/MOODS/bin:$PATH"

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

