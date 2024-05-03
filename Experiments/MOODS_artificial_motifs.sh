#!/bin/bash
array=$1 #this varies between 0 and 786 (number of all motifs)

export PATH="/projappl/project_2006203/softwares/conda_envs/MOODS/bin:$PATH"


S=/projappl/project_2006203/Genomes/Homo_sapiens/chr_sequences.fa #What is this?, this is repeat masked, contains only chromsomes chr1-22, X, Y
#S=/projappl/project_2006203/Genomes/Homo_sapiens/hg38.fa.masked #Repeats (RepeatMasker& Tandem Repeats Finder) masked by capital N, 
#contains also other genomes
#S=/projappl/project_2006203/Genomes/Homo_sapiens/hg38.fa

n=300000 #the number of best hits

#outfolder=/scratch/project_2006203/TFBS/Results/MOODS_artificial_human_final/
#mkdir $outfolder

#outfolder=/scratch/project_2006203/TFBS/Results/MOODS_artificial_human_final_union_zstd/ Can not read in zstd
#mkdir $outfolder

#outfolder=/scratch/project_2006203/TFBS/Results/MOODS_artificial_human_final_union/ 
#outfolder=/scratch/project_2006203/TFBS/Results/MOODS_artificial_human_final_version2.2/ #1000 arrays
outfolder=/scratch/project_2006203/TFBS/Results/MOODS_artificial_human_final_version2.2_correct/ #2000 arrays 
mkdir $outfolder


#Read in motif IDs
readarray -t lines < /projappl/project_2006203/TFBS/PWMs_final_version2.2/motifnames.csv #4003

start_ind=$(($array*5)) 
end_ind=$((($array+1)*5 ))

if [ "$end_ind" -gt 3933 ]; then
   
   end_ind=3933
fi

pwms_all=()

for (( index=$start_ind; index<$end_ind; index++ )); do


  pwms=($(ls /scratch/project_2006203/TFBS/PWMs_final_version2.2/artificial_motifs_space/${lines[$index]}*.pfm))
  echo ${lines[$index]}
  echo ${#pwms[@]}
  
  true_motif=${lines[$index]}

  #If the true_motif does not contain v2, remove v2 containing motifs them from the artificial
  # Test if the string contains the substring
  if [[ "$true_motif" == *"v2"* ]]; then
    #echo "The string '$true_motif' contains the substring v2. Do nothing"
    :
  else
    #echo "The string '$true_motif' does not contain the substring v2. Remove possible v2-motifs from the artificial list"
    substring="v2"
    # Print each element of the array on a new line, filter with grep -v (invert match), and read into a new array
    mapfile -t filtered_artificial_files < <(printf "%s\n" "${pwms[@]}" | grep -v "$substring")
    # Print the filtered_array to verify
    #printf "%s\n" "${filtered_artificial_files[@]}"
    pwms=("${filtered_artificial_files[@]}")
  fi

  pwms_all+=( "${pwms[@]}"  )

done



moods-dna.py -m ${pwms_all[@]} --threshold 2 -s $S | gzip > $LOCAL_SCRATCH"/MOODS_"$array".csv.gz" 

#moods-dna.py -m ${pwms_all[@]} --threshold 2 -s $S | zstd > $LOCAL_SCRATCH"/MOODS_"$array".zst" 

cd $LOCAL_SCRATCH
cp MOODS_"$array".csv.gz $outfolder 
#cp MOODS_"$array".zst $outfolder 

#gzip $outfile

#outfile to parquet, where does the file go, requires a lot of space
#csvcli $outfile convert -to "parquet"


#FATAL:   container creation failed: hook function for tag layer returns error: failed to create /tmp/nvme/job_14238673 directory: mkdir /tmp/nvme/job_14238673: permission denied

#/projappl/project_2006203/softwares/conda_envs/MOODS/bin/moods-dna.py: line 22: 705506 Segmentation fault
#   /usr/bin/singularity --silent exec -B $DIR/../$SQFS_IMAGE:$INSTALLATION_PATH:image-src=/
# $DIR/../$CONTAINER_IMAGE bash -c "eval \"\$(/CSC_CONTAINER/miniconda/bin/conda shell.bash hook )\"
# && conda activate env1 &>/dev/null &&  exec -a $_O_SOURCE $DIR/moods-dna.py $(
#test $# -eq 0 || printf " %q" "$@" )"
