#!/bin/bash
#SBATCH --job-name=csvcli
#SBATCH --output=csvcli.out
#SBATCH --error=csvcli.err
#SBATCH --account=project_2006203
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=30G
#SBATCH --partition=small
##SBATCH --mail-type=BEGIN #uncomment to enable mail



#export PATH="/projappl/project_2006203/softwares/conda_envs/MOODS/bin:$PATH"
export PATH="/projappl/project_2006203/softwares/conda_envs/MOODS-csvcli/bin":$PATH
#export PATH="/scratch/project_2006203/TFBS/Experiments/moods/bin:$PATH"


#S=/projappl/project_2006203/Genomes/Homo_sapiens/chr_sequences.fa #What is this?, this is repeat masked, contains only chromsomes chr1-22, X, Y
#S=/projappl/project_2006203/Genomes/Homo_sapiens/hg38.fa.masked #Repeats (RepeatMasker& Tandem Repeats Finder) masked by capital N, contains also other genomes
#S=/projappl/project_2006203/Genomes/Homo_sapiens/hg38.fa

#n=300000 #the number of best hits

#outfile=../Results/MOODS/TCF7.gz
outfile=../Results/MOODS/TCF7.csv

#srun moods-dna.py -m ../PWMs/Yin2017/pwms_space/Homo_sapiens/TCF7_HT-SELEX_TCGGAC40NCCA_KS_ASATCAAAS_1_4_NO.pfm \
#../PWMs/Yin2017/pwms_space/Homo_sapiens/TCF7_HT-SELEX_TCGGAC40NCCA_KS_WCATCGRGRCGCTGW_2_4_NO.pfm \
#../PWMs/Yin2017/pwms_space/Homo_sapiens/TCF7_Methyl-HT-SELEX_TTGTCT40NGCT_KS_ASATCAAAS_1_4_NO.pfm \
#../PWMs/Yin2017/pwms_space/Homo_sapiens/TCF7_Methyl-HT-SELEX_TTGTCT40NGCT_KS_WCATCGRGRCGCTGW_2_4_NO.pfm \
#-s $S -B $n | gzip > $outfile

#srun moods-dna.py -m ../PWMs/Yin2017/pwms_space/Homo_sapiens/TCF7_HT-SELEX_TCGGAC40NCCA_KS_ASATCAAAS_1_4_NO.pfm \
#../PWMs/Yin2017/pwms_space/Homo_sapiens/TCF7_HT-SELEX_TCGGAC40NCCA_KS_WCATCGRGRCGCTGW_2_4_NO.pfm \
#../PWMs/Yin2017/pwms_space/Homo_sapiens/TCF7_Methyl-HT-SELEX_TTGTCT40NGCT_KS_ASATCAAAS_1_4_NO.pfm \
#../PWMs/Yin2017/pwms_space/Homo_sapiens/TCF7_Methyl-HT-SELEX_TTGTCT40NGCT_KS_WCATCGRGRCGCTGW_2_4_NO.pfm \
#-s $S -B $n -o $outfile


csvcli $outfile convert -to parquet

seff $SLURM_JOBID

#FATAL:   container creation failed: hook function for tag layer returns error: failed to create /tmp/nvme/job_14238673 directory: mkdir /tmp/nvme/job_14238673: permission denied

#/projappl/project_2006203/softwares/conda_envs/MOODS/bin/moods-dna.py: line 22: 705506 Segmentation fault
#   /usr/bin/singularity --silent exec -B $DIR/../$SQFS_IMAGE:$INSTALLATION_PATH:image-src=/
# $DIR/../$CONTAINER_IMAGE bash -c "eval \"\$(/CSC_CONTAINER/miniconda/bin/conda shell.bash hook )\"
# && conda activate env1 &>/dev/null &&  exec -a $_O_SOURCE $DIR/moods-dna.py $(
#test $# -eq 0 || printf " %q" "$@" )"
