#!/bin/bash
array=$1 #this varies between 0 and 1000


cd /projappl/project_2006203/motifsimilarity-private/experiments
# cd csc_scratch/motifsimilarity-private/experiments/

results_path=/scratch/project_2006203/motifsimilarity-private/results_final_version2.2
mkdir $results_path



start_ind=$(($array*7730))
end_ind=$((($array+1)*7730 ))



for (( index=$start_ind; index<$end_ind; index++ )); do

tmp1=$((1+8*$index))
tmp=`echo sqrt "($tmp1)" | bc -l`
tmp1=`echo "($tmp)+1" | bc -l`
i=`echo "($tmp1)/2" | bc -l`
i=${i%.*}
j=$(($index-$i*($i-1)/2))

#echo "index: "$index
echo "i: "$i
echo "j: "$j


filenames=( $(cut -d ',' -f1 ../../TFBS/PWMs_final_version2.2/filenames_similarity_computations.csv ) ) #These are the full filenames

length=${#filenames[@]}
TF_PATH=../../TFBS/
TF1=${filenames[$i]}
TF2=${filenames[$j]}



echo $TF1
echo $TF2
#head -n 1 $TF_PATH$TF1 | awk '{print NF}'
#10
#head -n 1 $TF_PATH$TF2 | awk '{print NF}'


motifnames=( $(cut -d ',' -f1 ../../TFBS/PWMs_final_version2.2/motifnames.csv ) )
PWM1=${motifnames[$i]}
PWM2=${motifnames[$j]}


if [[ $index -eq $start_ind ]]
then
     #echo $index is equal to start
     echo $PWM1"-"$PWM2 > $LOCAL_SCRATCH"/result_"$array".out"
     ../motifsimilarity-parallel $TF_PATH$TF1 $TF_PATH$TF2 gapped 10 >> $LOCAL_SCRATCH"/result_"$array".out"


else
     #echo $index is greater then start
     echo $PWM1"-"$PWM2 >> $LOCAL_SCRATCH"/result_"$array".out"
     ../motifsimilarity-parallel $TF_PATH$TF1 $TF_PATH$TF2 gapped 10 >> $LOCAL_SCRATCH"/result_"$array".out"

fi



done

cd $LOCAL_SCRATCH
cp "result_"$array".out" $results_path"/"



