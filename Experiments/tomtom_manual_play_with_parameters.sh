#!/bin/bash
#SBATCH --job-name=tomtom_relaxed
#SBATCH --output=outs/tomtom_relaxed.out
#SBATCH --error=errs/tomtom_relaxed.err
#SBATCH --account=project_2006203
#SBATCH --time=20:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --partition=small
##SBATCH --mail-type=BEGIN #uncomment to enable mail


#/scratch/project_2006203/motif-clustering-Viestra-private/motif-clustering-Vierstra.yml (this has been likely changed)
#Before this you need to run code/run_scpd2meme.sh
export PATH="/projappl/project_2006203/softwares/conda_envs/motif-clustering-Vierstra/bin:$PATH"
#meme -version 5.4.1

tomtom -h

#-png -eps -no-ssc Draw motif alignment figures, produces millions of files, not sensible

#Draw the the figures for the alignment of the TF-TF-pair motifs and the individual TFs

#meme_file=../PWMs_final/all.meme
#meme_file=../PWMs_final_union/all.meme
meme_file=../PWMs_final_version2.2/all.meme

#results_path=/scratch/project_2006203/Results/tomtom_final/tomtom_relaxed/
results_path=/scratch/project_2006203/Results/tomtom_final_version2.2/tomtom_relaxed/

mkdir -p $results_path

tomtom -dist kullback -motif-pseudo 0.1 -thresh 1 -min-overlap 1 $meme_file $meme_file -oc $results_path

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


