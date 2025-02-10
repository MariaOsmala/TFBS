#!/bin/bash
#SBATCH --job-name=liftOver_TFBS
#SBATCH --output=outs/liftOver_TFBS.out
#SBATCH --error=errs/liftOver_TFBS.err
#SBATCH --account=project_2006203
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=50G
#SBATCH --cpus-per-task=1

#22768136

export PATH="/projappl/project_2007567/softwares/ucsc-tools:$PATH" 

liftOver_chainfile_path=/projappl/project_2006203/liftOver/


#fetchChromSizes hg19 > hg19.chrom.sizes

module load biokit
module load bedops



#mouse mm39 -> mm10
cd /scratch/project_2006203/TFBS/Results/MOODS_mouse_mm39_version2.2_processed/


mkdir MOODS_bigbed_mm10/
mkdir MOODS_bigbed_mm10_sorted/
mkdir MOODS_bigbed_mm10_representatives_sorted/

pwms=( $(cut -d ',' -f1 /projappl/project_2006203/TFBS/PWMs_final_version2.2/representatives.csv ) )
for element in "${pwms[@]}"
do
    # There are 3 extra columns in addition to chr start end
    liftOver -bedPlus=3 MOODS_bigbed/$element"_top.bed"  $liftOver_chainfile_path"mm39ToMm10.over.chain" MOODS_bigbed_mm10/$element"_top.bed" unMapped.$element"_top.bed"
    #echo "bedSort"
    rm unMapped.$element"_top.bed"
    bedSort MOODS_bigbed_mm10/$element"_top.bed" MOODS_bigbed_mm10_representatives_sorted/$element"_top.bed"
done

#human hg39 -> 19, DONE!
# cd /scratch/project_2006203/TFBS/Results/MOODS_human_final_version2.2_correct_processed 
# 
# mkdir MOODS_bigbed_hg19/
# mkdir MOODS_bigbed_hg19_sorted/
# mkdir MOODS_bigbed_hg19_representatives_sorted/
# 
# for file in MOODS_bigbed/*; do
#     filename=$(basename "$file")
#     # There are 3 extra columns in addition to chr start end
#     liftOver -bedPlus=3 MOODS_bigbed/$filename $liftOver_chainfile_path"hg38ToHg19.over.chain" MOODS_bigbed_hg19/$filename unMapped.$filename
#     #echo "bedSort"
#     bedSort MOODS_bigbed_hg19/$filename MOODS_bigbed_hg19_sorted/$filename
#     
# done
# 
# #Move representatives to separate folder
# 
# pwms=( $(cut -d ',' -f1 /projappl/project_2006203/TFBS/PWMs_final_version2.2/representatives.csv ) )
# for element in "${pwms[@]}"
# do
#     echo "Processing $element"
#     # You can perform other operations with $element here
#     mv MOODS_bigbed_hg19_sorted/$element"_top.bed" MOODS_bigbed_hg19_representatives_sorted/
# done


