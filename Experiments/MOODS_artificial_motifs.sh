#!/bin/bash

#!/bin/bash
array=$1 #this varies between 0 and 1030 (number of representatives)
echo $array

#for array > 1000
#num2=1000
#array=$[array + num2]

echo $array

export PATH="/projappl/project_2006203/softwares/conda_envs/MOODS/bin:$PATH"


S=/projappl/project_2006203/Genomes/Homo_sapiens/chr_sequences.fa #What is this?, this is repeat masked, contains only chromsomes chr1-22, X, Y
#S=/projappl/project_2006203/Genomes/Homo_sapiens/hg38.fa.masked #Repeats (RepeatMasker& Tandem Repeats Finder) masked by capital N, 
#contains also other genomes
#S=/projappl/project_2006203/Genomes/Homo_sapiens/hg38.fa

n=300000 #the number of best hits

outfolder=/scratch/project_2006203/TFBS/Results/MOODS_artificial_human_final/
mkdir $outfolder

#outfolder=/scratch/project_2006203/TFBS/Results/MOODS_artificial_human_zstd/
#mkdir $outfolder



#pwms=($(ls /scratch/project_2006203/TFBS/artificial_motifs/artificial_motifs_space/PROX1_HOXA2_*.pfm))

#export PATH="/projappl/project_2006203/softwares/conda_envs/MOODS/bin:$PATH"

readarray -t lines < /projappl/project_2006203/TFBS/PWMs_final/representatives.csv

# substring="ERF_DLX2_CAP-SELEX_TCGAAA40NAAT_AAC_RSCGGAANNNNNYMATTA_1_3_v2" #336
# substring=MTF1_HT-SELEX_TGCGCA40NCAA_KP_NYTTTGCACACGRYNYN_2_3 #1001
# substring=ZBTB12_HT-SELEX_TACTGG40NGGA_KV_NGGCCTGCCGTCNT_1_4 #1020
# substring=ZBTB18_HT-SELEX_TGCAGC40NATT_KR_NWTCCAGATGTKN_1_2 #1026
# substring=ZBTB22_HT-SELEX_TGCCTC40NCAA_KP_NKCACTANNNTAGTGMN_1_3 #1022
# substring=ZBTB2_HT-SELEX_TTAAAT40NAAT_KT_NTTTMCGGTWAN_1_4 #1030

#1,20,22,26,30

#substring="PROX1_HOXA2"
#indices=$(printf "%s\n" "${lines[@]}" | awk -v s="$substring" '$0 ~ s {print NR-1}') #644
#echo $indices
#pwms=($(ls /scratch/project_2006203/TFBS/artificial_motifs/artificial_motifs_space/${lines[644]}"_"*".pfm"))

#"/scratch/project_2006203/TFBS/artificial_motifs/artificial_motifs_space/"${lines[644]} *.pfm

pwms=($(ls /scratch/project_2006203/TFBS/artificial_motifs/artificial_motifs_space/${lines[$array]}*.pfm))

echo ${#pwms[@]}

true_motif=${lines[$array]}

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



export PATH="/projappl/project_2006203/softwares/conda_envs/MOODS/bin:$PATH"

moods-dna.py -m ${pwms[@]} --threshold 2 -s $S | gzip > $LOCAL_SCRATCH"/MOODS_"${lines[$array]}".csv.gz" 

#moods-dna.py -m ${pwms[@]} --threshold 2 -s $S | zstd > $LOCAL_SCRATCH"/MOODS_"${lines[$array]}".zst" 

cd $LOCAL_SCRATCH
cp MOODS_"${lines[$array]}".csv.gz $outfolder 
#cp MOODS_"${lines[$array]}".zst $outfolder 

#gzip $outfile

#outfile to parquet, where does the file go, requires a lot of space
#csvcli $outfile convert -to "parquet"


#FATAL:   container creation failed: hook function for tag layer returns error: failed to create /tmp/nvme/job_14238673 directory: mkdir /tmp/nvme/job_14238673: permission denied

#/projappl/project_2006203/softwares/conda_envs/MOODS/bin/moods-dna.py: line 22: 705506 Segmentation fault
#   /usr/bin/singularity --silent exec -B $DIR/../$SQFS_IMAGE:$INSTALLATION_PATH:image-src=/
# $DIR/../$CONTAINER_IMAGE bash -c "eval \"\$(/CSC_CONTAINER/miniconda/bin/conda shell.bash hook )\"
# && conda activate env1 &>/dev/null &&  exec -a $_O_SOURCE $DIR/moods-dna.py $(
#test $# -eq 0 || printf " %q" "$@" )"
