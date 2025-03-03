#!/bin/bash
#SBATCH --job-name=tomtom
#SBATCH --output=outs/tomtom.out
#SBATCH --error=errs/tomtom.err
#SBATCH --account=project_2006203
#SBATCH --time=2-00:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --gres=nvme:20 
#SBATCH --partition=small
##SBATCH --mail-type=BEGIN #uncomment to enable mail
##SBATCH --gres=nvme:50 This is in GB


#/scratch/project_2006203/motif-clustering-Viestra-private/motif-clustering-Vierstra.yml (this has been likely changed)
#Before this you need to run code/run_scpd2meme.sh
export PATH="/projappl/project_2006203/softwares/conda_envs/motif-clustering-Vierstra/bin:$PATH"
#meme -version 5.4.1

cp ../PWMs_final/all.meme $LOCAL_SCRATCH
cp ../PWMs_final_union/all.meme $LOCAL_SCRATCH

#results_path=/scratch/project_2006203/Results/tomtom_final/tomtom/ Earlier results
#results_path=/scratch/project_2006203/Results/tomtom_final/tomtom_test/ just testing
results_path=/scratch/project_2006203/Results/tomtom_union_final/tomtom/ #final set of motifs

mkdir -p $results_path

tomtom -dist kullback -motif-pseudo 0.1 -text -min-overlap 1 $LOCAL_SCRATCH/all.meme $LOCAL_SCRATCH/all.meme > $LOCAL_SCRATCH/tomtom.all.txt

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

cd $LOCAL_SCRATCH
cp tomtom.all.txt $results_path
#cp tomtom.all.txt /scratch/project_2006203/Results/tomtom_final/tomtom/
