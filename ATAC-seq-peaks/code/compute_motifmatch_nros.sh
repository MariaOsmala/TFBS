cd /scratch/project_2006203/TFBS/Results/MOODS_Teemu_processed/MOODS_bigbed



top=($(ls *_top.bed)) #3982

tmp=(${top[@]%%_top.bed})
all=("${tmp[@]/%/.bed}")

#count the motif matches in bed files
#awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' ZSCAN5A_Methyl-HT-SELEX_TCTTGC40NTTC_KR_NYGTCCCYCCCCAAANMN_2_2_NO_top.bed

# get length of an array
length=${#top[@]}

echo "motif top all" > "../../Teemu_5_motif_match_numbers.csv"
 
# use C style for loop syntax to read all values and indexes
for (( j=0; j<length; j++ ));
do
  a=($(wc -l ${top[$j]}))
 
  top_nro=${a[0]}
  a=($(wc -l ${all[$j]}))
  all_nro=${a[0]}
  motif=(${top[j]%%_top.bed})
  printf "%s %d %d \n" "${motif}" $top_nro $all_nro >> "../../Teemu_5_motif_match_numbers.csv"
done

#/scratch/project_2006203/TFBS/Results/Teemu_5_motif_match_numbers.csv

# for i in "${top[@]}"
# do
# 	echo "$i"
# done


# i=0
# for file in "${top[@]}"
# do
#    #echo $i
#    #echo $file
#    #echo ${filenames2[@]/$file//} | cut -d/ -f1 | wc -w | tr -d ' '
#    tmp=$(echo ${filenames2[@]/$file//} | cut -d/ -f1 | wc -w | tr -d ' ')
#    all_file=
#    
#    #echo $tmp
#    indices[$i]=$(echo ${filenames2[@]/$file//} | cut -d/ -f1 | wc -w | tr -d ' ')
#    pwms_new[$i]=${pwms[${indices[$i]}]}
#    i=$((i+1))
# done