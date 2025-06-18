export PATH="/projappl/project_2007567/softwares/ucsc-tools:$PATH" 


cd /scratch/project_2006203/TFBS/Results/MOODS_hg38.analysisSset_processed
cd /scratch/project_2006203/TFBS/Results/MOODS_T2T-CHM13_processed

module load biokit #for tabix

fetchChromSizes hg38 > hg38.chrom.sizes #Is this needed?

mkdir MOODS_bed_sorted

# metadata=read.table("/projappl/project_2006203/TFBS/Data/SELEX-motif-collection/metadata_final.tsv", header=TRUE, sep="\t")
# 
# library("tidyverse")
# 
# reps=metadata %>% filter(representative=="YES") %>% select(ID) %>% pull(ID)
# 
# write.table(reps, file="/projappl/project_2006203/TFBS/Data/SELEX-motif-collection/representatives.csv", sep=",", 
#           col.names=FALSE, row.names = FALSE,quote=FALSE)



readarray -t lines < /projappl/project_2006203/TFBS/Data/SELEX-motif-collection/representatives.csv #1232

mkdir -p MOODS_bed_sorted/bgzip
mkdir -p MOODS_bed_sorted/tabix

# Loop through the array elements
for TF in "${lines[@]}"
do
    echo $TF
    bedSort MOODS_bigbed/$TF"_top.bed" MOODS_bed_sorted/bgzip/$TF"_top.bed"
    
    bgzip MOODS_bed_sorted/bgzip/$TF"_top.bed"
    tabix -p bed MOODS_bed_sorted/bgzip/$TF"_top.bed.gz"
    mv MOODS_bed_sorted/bgzip/$TF"_top.bed.gz.tbi" MOODS_bed_sorted/tabix/
    
done


#tar -cjvf MOODS_matches_human_hg38analysisSet_representative.tar.bz MOODS_bed_sorted/
tar -cjvf MOODS_matches_human_T2T-CHM13_representative.tar.bz MOODS_bed_sorted/

#Move to allas

#Creates a new bucket with public access and uploads the data to the bucket. 

#Command a-publish creates the bucket and uploads the selected files into it. 
#Parameter -b is used to define the name for the bucket, in this case TFBS-project-public.

module load allas
allas-conf project_2006203
a-publish -b TFBS-project-public MOODS_matches_human_hg38analysisSet_representative.tar.bz MOODS_matches_human_hg38analysisSet_representative.tar.bz
a-publish -b TFBS-project-public MOODS_matches_human_T2T-CHM13_representative.tar.bz MOODS_matches_human_T2T-CHM13_representative.tar.bz

#bedToBigBed in.bed chrom.sizes out.bb

#Where in.bed is in one of the ascii bed formats, but not including track lines
#and chrom.sizes is a two-column file/URL: <chromosome name> <size in bases>
#and out.bb is the output indexed big bed file.



