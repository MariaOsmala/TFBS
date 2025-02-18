# "Soft-masked" assembly sequence in one file.
#     Repeats from RepeatMasker and Tandem Repeats Finder (with period of 12 or
#     less) are shown in lower case; non-repeating sequence is shown in upper
#     case. (again, the most current version of this file is latest/hg38.fa.gz)


cd ~/projects/TFBS/Genomes/hg38

rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz .

gunzip hg38.fa.gz

# "Hard-masked" assembly sequence in one file.
#     Repeats are masked by capital Ns; non-repeating sequence is shown in
#     upper case.

rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.masked.gz .

rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/README.txt .

rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes .

rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromAlias.txt .

gunzip hg38.fa.masked.gz


 
# mamba create -n bioinfo -c bioconda -c conda-forge samtools
# mamba activate bioinfo
# samtools --version
# samtools 1.21
# Using htslib 1.21

cd ~/projects/TFBS/Genomes/hg38

samtools faidx hg38.fa.masked

for chr in {1..22} X Y; do
    samtools faidx hg38.fa.masked "chr${chr}" >> chr_sequences.fa
done
