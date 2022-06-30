

cd /projappl/project_2006203/softwares

module load tykky

conda-containerize new --mamba --prefix conda_envs/MOODS MOODS.yml

export PATH="/projappl/project_2006203/softwares/conda_envs/MOODS/bin:$PATH"



#when moods installed somewhere else
#export PATH=/projappl/project_2006203/softwares/MOODS-python-1.9.4.1/scripts:$PATH

#python setup.py build_ext --inplace
#This builds the MOODS extensions and puts everything necessary under the MOODS/ directory. However, this does not install MOODS to Python library path (see below), so Python will not find the library unless MOODS/ directory is in the script path or you explicitly add it to sys.path. To use the example scripts under scripts/ directory, for instance, you can symlink the library there:

#cd scripts/
#ln -s ../MOODS/  


cd /scratch/project_2006203/TFBS/Experiments


S=/projappl/project_2006203/Genomes/Homo_sapiens/chr_sequences.fa

n=300000 #the number of best hits

outfile=../Results/MOODS/TCF7.csv

moods-dna.py -S ../PWMs/Yin2017/pwms/Homo_sapiens/TCF7_eDBD_HMG_HT-SELEX_TCGGAC40NCCA_KS_ASATCAAAS_1_4_NO.pfm \
../PWMs/Yin2017/pwms/Homo_sapiens/TCF7_eDBD_HMG_Methyl-HT-SELEX_TTGTCT40NGCT_KS_ASATCAAAS_1_4_NO.pfm \
../PWMs/Yin2017/pwms/Homo_sapiens/TCF7_eDBD_HMG_HT-SELEX_TCGGAC40NCCA_KS_WCATCGRGRCGCTGW_2_4_NO.pfm \
../PWMs/Yin2017/pwms/Homo_sapiens/TCF7_eDBD_HMG_Methyl-HT-SELEX_TTGTCT40NGCT_KS_WCATCGRGRCGCTGW_2_4_NO.pfm \
 -s $S -B $n -o $outfile \
 /
