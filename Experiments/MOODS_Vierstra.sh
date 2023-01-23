#!/bin/bash

#!/bin/bash
array=$1 #this varies between 0 and 399



export PATH="/projappl/project_2006203/softwares/conda_envs/MOODS/bin:$PATH"


S=/projappl/project_2006203/Genomes/Homo_sapiens/chr_sequences.fa #What is this?, this is repeat masked, contains only chromsomes chr1-22, X, Y
#S=/projappl/project_2006203/Genomes/Homo_sapiens/hg38.fa.masked #Repeats (RepeatMasker& Tandem Repeats Finder) masked by capital N, contains also other genomes
#S=/projappl/project_2006203/Genomes/Homo_sapiens/hg38.fa

n=300000 #the number of best hits

outfile="/scratch/project_2006203/TFBS/Results/MOODS_Vierstra/MOODS_"$array".csv.gz"
#outfile=../Results/MOODS/MOODS_"$array".csv"

pwms=($(ls ../PWMs/*/pwms_space/*/*.pfm)) #3982
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


moods-dna.py -m ${pwms[@]:$start_ind:$length} --p-value 1e-4  --lo-bg 0.2977 0.2023 0.2023 0.2977 -s $S | gzip > $outfile #${array[@]:START:LENGTH}

#gzip $outfile

#outfile to parquet, where does the file go, requires a lot of space
#csvcli $outfile convert -to "parquet"


#FATAL:   container creation failed: hook function for tag layer returns error: failed to create /tmp/nvme/job_14238673 directory: mkdir /tmp/nvme/job_14238673: permission denied

#/projappl/project_2006203/softwares/conda_envs/MOODS/bin/moods-dna.py: line 22: 705506 Segmentation fault
#   /usr/bin/singularity --silent exec -B $DIR/../$SQFS_IMAGE:$INSTALLATION_PATH:image-src=/
# $DIR/../$CONTAINER_IMAGE bash -c "eval \"\$(/CSC_CONTAINER/miniconda/bin/conda shell.bash hook )\"
# && conda activate env1 &>/dev/null &&  exec -a $_O_SOURCE $DIR/moods-dna.py $(
#test $# -eq 0 || printf " %q" "$@" )"
