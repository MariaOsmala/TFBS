export PATH="/projappl/project_2007567/softwares/ucsc-tools:$PATH" 

#cd /scratch/project_2006203/TFBS/Results/MOODS_human_final_processed/ old
cd /scratch/project_2006203/TFBS/Results/MOODS_human_final_version2.2_correct_processed

module load biokit #for tabix

fetchChromSizes hg38 > hg38.chrom.sizes

mkdir MOODS_bigbed_sorted
#mkdir MOODS_bigbed_representatives

# metadata=read.table("/projappl/project_2006203/TFBS/PWMs_final_version2.2/metadata_representatives.tsv", header=TRUE, sep="\t")
# 
# library("tidyverse")
# 
# reps=metadata %>% filter(new_representative=="YES") %>% select(ID) %>% pull(ID)
# 
# write.table(reps, file="/projappl/project_2006203/TFBS/PWMs_final_version2.2/representatives.csv", sep=",", 
#            col.names=FALSE, row.names = FALSE,quote=FALSE)



readarray -t lines < /projappl/project_2006203/TFBS/PWMs_final_version2.2/representatives.csv #1232

mkdir -p MOODS_bigbed_sorted/bgzip
mkdir -p MOODS_bigbed_sorted/tabix

# Loop through the array elements
for TF in "${lines[@]}"
do
    echo $TF
    bedSort MOODS_bigbed/$TF"_top.bed" MOODS_bigbed_sorted/bgzip/$TF"_top.bed"
    
    bgzip MOODS_bigbed_sorted/bgzip/$TF"_top.bed"
    tabix -p bed MOODS_bigbed_sorted/bgzip/$TF"_top.bed.gz"
    mv MOODS_bigbed_sorted/bgzip/$TF"_top.bed.gz.tbi" MOODS_bigbed_sorted/tabix/
    
    #bedtobedGraph
    #awk '{ print $1"\t"$2"\t"$3"\t"$5 }' MOODS_bigbed_sorted/$TF"_top.bed" > MOODS_bedGraph_representatives/$TF"_top.bedGraph";
    
    #bedGraph to bigWig
    
    
    #bedToBigBed MOODS_bigbed_sorted/$TF"_top.bed" hg38.chrom.sizes MOODS_bigbed_representatives/$TF"_top.bb"
done


tar -cjvf MOODS_matches_human_hg38_representative.tar.bz MOODS_bigbed_sorted/

#Move to allas

#Creates a new bucket with public access and uploads the data to the bucket. 

#Command a-publish creates the bucket and uploads the selected files into it. 
#Parameter -b is used to define the name for the bucket, in this case TFBS-project-public.

module load allas
allas-conf project_2006203
a-publish -b TFBS-project-public MOODS_matches_human_hg38_representative.tar.bz MOODS_matches_human_hg38_representative.tar.bz


#bedToBigBed in.bed chrom.sizes out.bb

#Where in.bed is in one of the ascii bed formats, but not including track lines
#and chrom.sizes is a two-column file/URL: <chromosome name> <size in bases>
#and out.bb is the output indexed big bed file.



