cd /home/osmalama/Dropbox/Taipale-lab/TFBS/Experiments

source activate /home/osmalama/softwares/conda_envs/MOODS



export PATH=/home/osmalama/softwares/MOODS-python-1.9.4.1/scripts:$PATH


../PWMs/Yin2017/pwms/Homo_sapiens/TCF7_eDBD_HMG_HT-SELEX_TCGGAC40NCCA_KS_ASATCAAAS_1_4_NO.pfm
../PWMs/Yin2017/pwms/Homo_sapiens/TCF7_eDBD_HMG_Methyl-HT-SELEX_TTGTCT40NGCT_KS_ASATCAAAS_1_4_NO.pfm
../PWMs/Yin2017/pwms/Homo_sapiens/TCF7_eDBD_HMG_HT-SELEX_TCGGAC40NCCA_KS_WCATCGRGRCGCTGW_2_4_NO.pfm
../PWMs/Yin2017/pwms/Homo_sapiens/TCF7_eDBD_HMG_Methyl-HT-SELEX_TTGTCT40NGCT_KS_WCATCGRGRCGCTGW_2_4_NO.pfm

#moods-dna.py [-h] [-v]
#[-m M [M ...]] pfm's
#[-S M [M ...]] pwm's
#[-s S [S ...]] sequence
#                   [-p p] [-t T] [-B n] [-o outfile] [--sep S] [-R]
#                   [--no-snps] [--batch] [--bg pA pC pG pT] [--ps p]
#                   [--log-base x] [--lo-bg pA pC pG pT]
#                   [--threshold-precision x]


#M are PWM matrices of size 4 x length, space separated

#S is the sequence
#S=/home/osmalama/Genomes/Homo_sapiens/hg38.fa.masked
S=/home/osmalama/Genomes/Homo_sapiens/chr_sequences.fa #keep only chr1-22 chrX-Y

n=300000 #the number of best hits

outfile=../Results/MOODS/TCF7.csv

moods-dna.py -S ../PWMs/Yin2017/pwms/Homo_sapiens/TCF7_eDBD_HMG_HT-SELEX_TCGGAC40NCCA_KS_ASATCAAAS_1_4_NO.pfm \
../PWMs/Yin2017/pwms/Homo_sapiens/TCF7_eDBD_HMG_Methyl-HT-SELEX_TTGTCT40NGCT_KS_ASATCAAAS_1_4_NO.pfm \
../PWMs/Yin2017/pwms/Homo_sapiens/TCF7_eDBD_HMG_HT-SELEX_TCGGAC40NCCA_KS_WCATCGRGRCGCTGW_2_4_NO.pfm \
../PWMs/Yin2017/pwms/Homo_sapiens/TCF7_eDBD_HMG_Methyl-HT-SELEX_TTGTCT40NGCT_KS_WCATCGRGRCGCTGW_2_4_NO.pfm \
 -s $S -B $n -o $outfile \
 /
