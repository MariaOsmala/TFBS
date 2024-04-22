#!/bin/bash
#SBATCH --job-name=tomtom_figures
#SBATCH --output=outs/tomtom_figures.out
#SBATCH --error=errs/tomtom_figures.err
#SBATCH --account=project_2006203
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --gres=nvme:20 
#SBATCH --partition=small



#/scratch/project_2006203/motif-clustering-Viestra-private/motif-clustering-Vierstra.yml (this has been likely changed)
#Before this you need to run code/run_scpd2meme.sh
export PATH="/projappl/project_2006203/softwares/conda_envs/motif-clustering-Vierstra/bin:$PATH"
#meme -version 5.4.1
#Draw the the figures for the alignment of the TF-TF-pair motifs and the individual TFs

#!/bin/bash

data_file=/projappl/project_2006203/TFBS/PWMs_final_union/metadata_google_drive_supplement_for_tomtom_figures.tsv

# Define empty arrays
dimer_array=()
monomer1_array=()
monomer2_array=()

counter=0

# Read file
while IFS=$'\t' read -r col1 col2 col3
do
    ((counter++))
    
    # Skip the first line
    if [ $counter -eq 1 ]; then
        continue
    fi

    # echo $col1
    dimer_array+=("$col1")
    # echo ${dimer_array[0]}
    monomer1_array+=("$col2")
    monomer2_array+=("$col3")
done < $data_file

# To demonstrate that the data was read into the arrays, print them
echo "First Column:"
printf "%s\n" "${#dimer_array[@]}"

echo "Third Column:"
printf "%s\n" "${#monomer2_array[@]}"



#dimer="ALX4_EOMES_CAP-SELEX_TGCACG40NTTG_AAD_NGYGYTAAYNNNNNNTNACACNN_1_3" 
#monomer1="ALX4_HT-SELEX_TCTATT40NCAT_KW_CYAATTAN_1_3"
#monomer2="EOMES_HT-SELEX_TGTAAA30NAAG_AI_NAGGTGTGAAAWN_1_3"

meme_file=../PWMs_final_union/all.meme

#results_path=/scratch/project_2006203/Results/tomtom_final/tomtom_relaxed/
    results_path=/scratch/project_2006203/Results/tomtom_union_final/tomtom_relaxed/figures

mkdir -p $results_path




# Iterate from 0 to 1417
for ((i = 0; i < 1418; i++)); do
    echo "Index: $i"
    
    dimer=${dimer_array[i]}
    monomer1=${monomer1_array[i]}
    monomer2=${monomer2_array[i]}
    
    grep -A 4 $dimer ../PWMs_final_union/all.scpd > $LOCAL_SCRATCH"/"$dimer".scpd"
    grep -A 4 $monomer1 ../PWMs_final_union/all.scpd > $LOCAL_SCRATCH"/"$monomer1".scpd"
    grep -A 4 $monomer2 ../PWMs_final_union/all.scpd > $LOCAL_SCRATCH"/"$monomer2".scpd"
    #convert to meme

    scpd2meme $LOCAL_SCRATCH"/"$dimer".scpd" -pseudo 1 > $LOCAL_SCRATCH"/"$dimer".meme"
    scpd2meme $LOCAL_SCRATCH"/"$monomer1".scpd" -pseudo 1 > $LOCAL_SCRATCH"/"$monomer1".meme"
    scpd2meme $LOCAL_SCRATCH"/"$monomer2".scpd" -pseudo 1 > $LOCAL_SCRATCH"/"$monomer2".meme"

    meme_dimer=$LOCAL_SCRATCH"/"$dimer".meme" 
    meme_monomer1=$LOCAL_SCRATCH"/"$monomer1".meme"
    meme_monomer2=$LOCAL_SCRATCH"/"$monomer2".meme"

    tomtom -dist kullback -motif-pseudo 0.1 -png -eps -no-ssc -thresh 1 -min-overlap 1 $meme_dimer $meme_monomer1 -oc $results_path
    tomtom -dist kullback -motif-pseudo 0.1 -png -eps -no-ssc -thresh 1 -min-overlap 1 $meme_dimer $meme_monomer2 -oc $results_path
    
done







#-motif-pseudo <pseudo count>
#                   Apply the pseudocount to the query and target motifs;
#                    default: apply a pseudocount of 0.1
#-min-overlap <int>
#                   Minimum overlap between query and target;
#                    default: 1
#Might be good to try: 
#-thresh <float>  Significance threshold; default: 0.5
#-incomplete-scores Ignore unaligned columns in computing scores default: use complete set of columns
# -png             Create PNG logos; default: don't create PNG logos
#  -eps             Create EPS logos; default: don't create EPS logos
#  -no-ssc          Don't apply small-sample correction to logos;
#                   default: use small-sample correction


