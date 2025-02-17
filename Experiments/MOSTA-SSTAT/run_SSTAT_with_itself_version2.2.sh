#!/bin/bash



cd /projappl/project_2006203/MOSTA-SSTAT/Experiments
# cd csc_scratch/motifsimilarity-private/experiments/

results_path=/scratch/project_2006203/MOSTA-SSTAT



start_ind=0
end_ind=3993 #3933


for (( index=$start_ind; index<$end_ind; index++ )); do

i=$index
j=$index

#echo "index: "$index
echo "i: "$i
echo "j: "$j

#From i j back to index
#index=j+i*(i-1)/2



filenames=( $(cut -d ',' -f1 ../../TFBS/PWMs_final_version2.2/filenames.csv ) )
length=${#filenames[@]}
TF_PATH=../../TFBS/
TF1=${filenames[$i]}
TF2=${filenames[$j]}

#convert "pwms" to "transfac"
TF1="${TF1/pwms/transfac}"
TF2="${TF2/pwms/transfac}"

echo $TF1
echo $TF2

motifnames=( $(cut -d ',' -f1 ../../TFBS/PWMs_final_version2.2/motifnames.csv ) )

PWM1=${motifnames[$i]}
PWM2=${motifnames[$j]}



echo "../../TFBS/"$TF1 > $results_path/data/$PWM1"_"$PWM2".list"
echo "../../TFBS/"$TF2 >> $results_path/data/$PWM1"_"$PWM2".list"

if [[ $index -eq $start_ind ]]
then
     #echo $index is equal to start

     ../sstat .5 list:$results_path/data/$PWM1"_"$PWM2".list" typeI 0.01 > $results_path"/results_final_version2.2/result_with_itself.out"



else
     #echo $index is greater then start
   
     ../sstat .5 list:$results_path/data/$PWM1"_"$PWM2".list" typeI 0.01 >> $results_path"/results_final_version2.2/result_with_itself.out"


fi

rm $results_path/data/$PWM1"_"$PWM2".list"


done
