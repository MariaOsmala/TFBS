#!/bin/bash

#!/bin/bash
array=$1 #this varies between 0 and 399



export PATH="/projappl/project_2006203/softwares/conda_envs/MOODS/bin:$PATH"


S=/projappl/project_2006203/Genomes/Homo_sapiens/chr_sequences.fa #What is this?, this is repeat masked, contains only chromsomes chr1-22, X, Y
#S=/projappl/project_2006203/Genomes/Homo_sapiens/hg38.fa.masked #Repeats (RepeatMasker& Tandem Repeats Finder) masked by capital N, contains also other genomes
#S=/projappl/project_2006203/Genomes/Homo_sapiens/hg38.fa

n=300000 #the number of best hits

mkdir /scratch/project_2006203/TFBS/Results/MOODS_Teemu_rest
outfile="/scratch/project_2006203/TFBS/Results/MOODS_Teemu_rest/MOODS_"$array".csv.gz"
#outfile=../Results/MOODS/MOODS_"$array".csv"

#rename Yimengs motifs, add 

#filenames=($(ls ../PWMs/fromYimeng/pwms_space/pfm_composite_new/*.pfm))
#for filename in "${filenames[@]}"
#do
#   newname=${filename%.pfm}"_composite.pfm"
#   echo $newname
#   mv $filename $newname
   
#done

#filenames=($(ls ../PWMs/fromYimeng/pwms_space/pfm_spacing_new/*.pfm))
#for filename in "${filenames[@]}"
#do
#   newname=${filename%.pfm}"_spacing.pfm"
#   echo $newname
#   mv $filename $newname
   
#done



pwms1=($(ls ../PWMs/Morgunova2015/pwms_space/*/*.pfm)) #3982

pwms2=($(ls ../PWMs/fromYimeng/pwms_space/PWM_bHLH-Homeobox/*.pfm)) #3982

pwms=( ${pwms1[*]} ${pwms2[*]}  )


#echo "${pwms[3]##*/}"

#echo "${#pwms[@]}"
#echo "${pwms[@]}"

nro_pwms=${#pwms[@]}



export PATH="/projappl/project_2006203/softwares/conda_envs/MOODS/bin:$PATH"





moods-dna.py -m ${pwms[@]} --threshold 5 -s $S | gzip > $outfile #${array[@]:START:LENGTH}

#gzip $outfile

#outfile to parquet, where does the file go, requires a lot of space
#csvcli $outfile convert -to "parquet"


#FATAL:   container creation failed: hook function for tag layer returns error: failed to create /tmp/nvme/job_14238673 directory: mkdir /tmp/nvme/job_14238673: permission denied

#/projappl/project_2006203/softwares/conda_envs/MOODS/bin/moods-dna.py: line 22: 705506 Segmentation fault
#   /usr/bin/singularity --silent exec -B $DIR/../$SQFS_IMAGE:$INSTALLATION_PATH:image-src=/
# $DIR/../$CONTAINER_IMAGE bash -c "eval \"\$(/CSC_CONTAINER/miniconda/bin/conda shell.bash hook )\"
# && conda activate env1 &>/dev/null &&  exec -a $_O_SOURCE $DIR/moods-dna.py $(
#test $# -eq 0 || printf " %q" "$@" )"
