#!/bin/bash
#SBATCH --job-name=sort_bed
#SBATCH --output=outs/sort_bed.out
#SBATCH --error=errs/sort_bed.err
#SBATCH --account=project_2006203
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=50G
#SBATCH --cpus-per-task=1




export PATH="/projappl/project_2007567/softwares/ucsc-tools:$PATH" 


cd /scratch/project_2006203/TFBS/Results/MOODS_mouse_mm39_version2.2_processed/

fetchChromSizes mm39 > mm39.chrom.sizes

mkdir MOODS_bigbed_sorted

readarray -t lines < /projappl/project_2006203/TFBS/PWMs_final_version2.2/representatives.csv #1232

# Loop through the array elements
for TF in "${lines[@]}"
do
    echo $TF
    bedSort MOODS_bigbed/$TF"_top.bed" MOODS_bigbed_sorted/$TF"_top.bed"
    
done




#bedToBigBed in.bed chrom.sizes out.bb

#Where in.bed is in one of the ascii bed formats, but not including track lines
#and chrom.sizes is a two-column file/URL: <chromosome name> <size in bases>
#and out.bb is the output indexed big bed file.



