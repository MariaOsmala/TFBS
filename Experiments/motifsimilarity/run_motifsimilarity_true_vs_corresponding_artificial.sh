#!/bin/bash
array=$1 #this varies between 0 and 393


cd Experiments/motifsimilarity


results_path=Results/motifsimilarity/results_final_true_vs_corresponding_artificial
mkdir - p $results_path

readarray -t motifs < Data/SELEX-motif-collection/motifnames.csv


filenames=( $(cut -d ',' -f1 Data/SELEX-motif-collection/filenames.csv ) ) #These are the full filenames of all motifs

start_ind=$(($array*10)) 
end_ind=$((($array+1)*10 ))

if [ "$end_ind" -gt 3933 ]; then
   
   end_ind=3933
fi




TF_PATH=TFBS/

write_new_file=1

for (( index=$start_ind; index<$end_ind; index++ )); do
 

  #True motif
  true_motif=${motifs[$index]}".pfm"
  
  true_motif_file=$(printf "%s\n" "${filenames[@]}" | grep "$true_motif")
  
  
  search_dir="Data/PWMs/scrambled_pfms_tab"
  

  # Empty array to hold the resulting files
  artificial_files=()

  true_motif=${motifs[$index]}

  # For each prefix/lines use glob to get files and add them to the result array
  for file in "$search_dir/$true_motif"*; do
      # Check if file exists (glob might return pattern if no files are found)
      if [[ -f $file ]]; then
          artificial_files+=("$file")
      fi
  done
  
  #If the true_motif does not contain v2, remove v2 containing motifs them from the artificial
  # Test if the string contains the substring
  if [[ "$true_motif" == *"v2"* ]]; then
    #echo "The string '$true_motif' contains the substring v2. Do nothing"
    :
  else
    #echo "The string '$true_motif' does not contain the substring v2. Remove possible v2-motifs from the artificial list"
    substring="v2"
    # Print each element of the array on a new line, filter with grep -v (invert match), and read into a new array
    mapfile -t filtered_artificial_files < <(printf "%s\n" "${artificial_files[@]}" | grep -v "$substring")
    # Print the filtered_array to verify
    #printf "%s\n" "${filtered_artificial_files[@]}"
    artificial_files=("${filtered_artificial_files[@]}")
  fi
  
  
  
  PWM1=$true_motif

  for TF in "${artificial_files[@]}"; do #Loop over artificial
  
    PWM2=$(basename "$TF")
    PWM2="${PWM2%.*}"
    if [ $write_new_file -eq 1 ];
    then
       #echo $index is equal to start
       echo $PWM1"-"$PWM2 > $LOCAL_SCRATCH"/result_"$array".out"
       ../motifsimilarity-parallel $TF_PATH$true_motif_file $TF gapped 10 >> $LOCAL_SCRATCH"/result_"$array".out"
       
       write_new_file=0
    else
       #echo $index is greater then start
       echo $PWM1"-"$PWM2 >> $LOCAL_SCRATCH"/result_"$array".out"
       ../motifsimilarity-parallel $TF_PATH$true_motif_file $TF gapped 10 >> $LOCAL_SCRATCH"/result_"$array".out"
    fi
  done #Loop over artificial
  
done #Loop over representatives

cd $LOCAL_SCRATCH
cp "result_"$array".out" $results_path"/"



