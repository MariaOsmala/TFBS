#!/bin/bash -l
#SBATCH --job-name=copy_MOODS
#SBATCH --account=project_2006472
#SBATCH --output=outs/copy_MOODS.txt
#SBATCH --error=errs/copy_MOODS.txt
#SBATCH --partition=small
#SBATCH --time=12:00:00 
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G #300G


#21093609

outfolder=/scratch/project_2006203/TFBS/Results/MOODS_human_final_version2.2_processed/MOODS_bigbed
mkdir -p $outfolder

filenames=( $(cut -d ',' -f1 ../../TFBS/PWMs_final_version2.2/MOODS_union.txt ) )
length=${#filenames[@]}

MOODS_union=/scratch/project_2006203/TFBS/Results/MOODS_human_final_union_processed/MOODS_bigbed/

for (( index=0; index<$length; index++ )); do
  echo $index
  TF1=${filenames[$index]}

  cp $MOODS_union$TF1 $outfolder

done

