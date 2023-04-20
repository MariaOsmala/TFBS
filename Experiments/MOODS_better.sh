#!/bin/bash

#!/bin/bash
array=$1 #this varies between 0 and 399



export PATH="/projappl/project_2006203/softwares/conda_envs/MOODS/bin:$PATH"


S=/projappl/project_2006203/Genomes/Homo_sapiens/chr_sequences.fa #What is this?, this is repeat masked, contains only chromsomes chr1-22, X, Y
#S=/projappl/project_2006203/Genomes/Homo_sapiens/hg38.fa.masked #Repeats (RepeatMasker& Tandem Repeats Finder) masked by capital N, contains also other genomes
#S=/projappl/project_2006203/Genomes/Homo_sapiens/hg38.fa

n=300000 #the number of best hits

outfile="/scratch/project_2006472/MOODS/MOODS_"$array".csv.gz"
#outfile=../Results/MOODS/MOODS_"$array".csv"

pwms=($(ls ../PWMs/*/pwms_space/*/*.pfm)) #3982
echo "${#pwms[@]}"
echo "${pwms[@]}"

nro_pwms=${#pwms[@]}


start_ind=$(($array*10)) #100
end_ind=$((($array+1)*10 -1)) #100
length=10 #100

if [[ $end_ind -gt $(($nro_pwms-1)) ]]
then
     echo $end_ind is greater than $(($nro_pwms-1))
     length=$(($nro_pwms-$start_ind+1))
fi



#Write this directly to a database

#threshold selection (exactly one required):
#  -p p, --p-value p     compute threshold from p-value
#  -t T, --threshold T   use specified absolute threshold
#  -B n, --best-hits n   return at least the specified amount of best matches

#What Teemu has used: 
#olen hakenut löysällä affineettirajalla ja valinnut sitten kärjestä halutun määrän. 
#Suuria määriä hakiessa jokin minimiaffiniteettiraja (esim. 9) on kuitenkin ollut käytössä, 
#ettei todella pitkille motiiveille väkisin etsitä ihan järjettömiä osumia.

#Jos muistan oikein, tuo esimerkkinä mainitsemani raja on 9 on alunperin EEL:istä, 
#jonka affiniteetit ovat log2 kun taas moodsissa oletuskantaluku on e (--log-base). 

#--log-base x          logarithm base for log-odds conversion (default natural logarithm)

#Jolma2015: using the program MOODS54 with a loose cut-off (P value <10−3 with flat background distribution) to obtain a large excess of putative
binding sites for each motif. All found sites were merged into one list and 10,000
non-overlapping highest affinity sites selected for conservation analysis regardless
of the motif identity (heterodimer or control).


#moods-dna.py -m ${pwms[@]:$start_ind:$length}  -s $S -B $n -o $outfile #${array[@]:START:LENGTH}
moods-dna.py -m ${pwms[@]:$start_ind:$length}  -s $S -B $n | gzip > $outfile #${array[@]:START:LENGTH}

gzip $outfile

#outfile to parquet, where does the file go, requires a lot of space
#csvcli $outfile convert -to "parquet"


#FATAL:   container creation failed: hook function for tag layer returns error: failed to create /tmp/nvme/job_14238673 directory: mkdir /tmp/nvme/job_14238673: permission denied

#/projappl/project_2006203/softwares/conda_envs/MOODS/bin/moods-dna.py: line 22: 705506 Segmentation fault
#   /usr/bin/singularity --silent exec -B $DIR/../$SQFS_IMAGE:$INSTALLATION_PATH:image-src=/
# $DIR/../$CONTAINER_IMAGE bash -c "eval \"\$(/CSC_CONTAINER/miniconda/bin/conda shell.bash hook )\"
# && conda activate env1 &>/dev/null &&  exec -a $_O_SOURCE $DIR/moods-dna.py $(
#test $# -eq 0 || printf " %q" "$@" )"
