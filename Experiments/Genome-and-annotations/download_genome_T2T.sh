#!/bin/bash -l
#SBATCH --job-name=download
#SBATCH --account=project_2006203
#SBATCH --output=outs/download.txt
#SBATCH --error=errs/download.txt
#SBATCH --partition=small
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G 


#See https://github.com/marbl/CHM13

#chm13v2.0_maskedY.rCRS.fa.gz: PARs on chrY hard masked to "N" and mitochondrion replaced with rCRS (AC:NC_012920.1)
 
cd /projappl/project_2007567/Genomes/T2T-CHM13

#aws s3 --no-sign-request cp s3://human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY_rCRS.fa.gz . #DONE


#gunzip chm13v2.0_maskedY_rCRS.fa.gz #DONE
#grep ">" chm13v2.0_maskedY_rCRS.fa | wc -l #DONE

# >chr1 CP068277.2 Homo sapiens isolate CHM13 chromosome 1
# >chr2 CP068276.2 Homo sapiens isolate CHM13 chromosome 2
# >chr3 CP068275.2 Homo sapiens isolate CHM13 chromosome 3
# >chr4 CP068274.2 Homo sapiens isolate CHM13 chromosome 4
# >chr5 CP068273.2 Homo sapiens isolate CHM13 chromosome 5
# >chr6 CP068272.2 Homo sapiens isolate CHM13 chromosome 6
# >chr7 CP068271.2 Homo sapiens isolate CHM13 chromosome 7
# >chr8 CP068270.2 Homo sapiens isolate CHM13 chromosome 8
# >chr9 CP068269.2 Homo sapiens isolate CHM13 chromosome 9
# >chr10 CP068268.2 Homo sapiens isolate CHM13 chromosome 10
# >chr11 CP068267.2 Homo sapiens isolate CHM13 chromosome 11
# >chr12 CP068266.2 Homo sapiens isolate CHM13 chromosome 12
# >chr13 CP068265.2 Homo sapiens isolate CHM13 chromosome 13
# >chr14 CP068264.2 Homo sapiens isolate CHM13 chromosome 14
# >chr15 CP068263.2 Homo sapiens isolate CHM13 chromosome 15
# >chr16 CP068262.2 Homo sapiens isolate CHM13 chromosome 16
# >chr17 CP068261.2 Homo sapiens isolate CHM13 chromosome 17
# >chr18 CP068260.2 Homo sapiens isolate CHM13 chromosome 18
# >chr19 CP068259.2 Homo sapiens isolate CHM13 chromosome 19
# >chr20 CP068258.2 Homo sapiens isolate CHM13 chromosome 20
# >chr21 CP068257.2 Homo sapiens isolate CHM13 chromosome 21
# >chr22 CP068256.2 Homo sapiens isolate CHM13 chromosome 22
# >chrX CP068255.2 Homo sapiens isolate CHM13 chromosome X
# >chrY CP086569.2 Homo sapiens isolate NA24385 chromosome Y
# >chrM NC_012920.1 Homo sapiens mitochondrion, complete genome

sed -E 's/^>.*chromosome ([0-9XY]+).*/>chr\1/; s/^>.*mitochondrion.*/>chrM/' chm13v2.0_maskedY_rCRS.fa > chm13v2.0_maskedY_rCRS_shorter_chr_names.fa
                  
module load biokit      
samtools faidx chm13v2.0_maskedY_rCRS_shorter_chr_names.fa                  
                    
                    
#https://projects.ensembl.org/hprc/

#mkdir -p /projappl/project_2007567/Genomes/T2T-CHM13/Annotations

#cd /projappl/project_2007567/Genomes/T2T-CHM13/Annotations

#wget https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/geneset/2022_07/Homo_sapiens-GCA_009914755.4-2022_07-genes.gff3.gz #DONE
#wget https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/geneset/2022_07/Homo_sapiens-GCA_009914755.4-2022_07-genes.gtf.gz #DONE


# Feature               | GFF3                                       | GTF
# Version               | GFF3 (most common)                         | Based on GFF2
# Attribute syntax      | ID=gene1;Name=XYZ                          | gene_id "gene1"; Name "XYZ";
# Relationship tracking | Explicit with ID/Parent                    | Inferred via gene/transcript ID
# Use cases             | Genome annotation pipelines, Ensembl, NCBI | UCSC Genome Browser, RNA-seq tools    
#      
     
mkdir -p /projappl/project_2007567/Genomes/T2T-CHM13/Cytobands
cd /projappl/project_2007567/Genomes/T2T-CHM13/Cytobands

# Cytobands (short for cytogenetic bands) are distinct regions on a chromosome 
# that appear as light and dark bands when stained and viewed under a microscope,
# especially using a technique called Giemsa staining (hence, G-banding).

# Chromosome regions that are distinguishable by differences in how they absorb stain.
# Dark bands (G-positive): Rich in adenine-thymine (AT) base pairs, gene-poor, more condensed chromatin.
# Light bands (G-negative): Rich in guanine-cytosine (GC) base pairs, gene-rich, less condensed chromatin.

# Help identify structural features of chromosomes.
# Used to pinpoint gene locations (e.g., the BRCA1 gene is located at 17q21.31, 
# where "17" is the chromosome, "q" is the long arm, and "21.31" is the cytoband).

#Important in karyotyping, identifying chromosomal abnormalities (like deletions, duplications, translocations).

#Each cytoband is labeled based on:
#Chromosome number (1–22, X, Y),
#Arm: p (short arm) or q (long arm),
#Region and band number, counted outward from the centromere.
#For example: 2p16.3
#Chromosome 2
#Short arm (p)
#Region 1, band 6, sub-band 3

wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_cytobands_allchrs.bed

#Segmental Duplications, v2022-03-11 in simple and full bed format

mkdir -p /projappl/project_2007567/Genomes/T2T-CHM13/Segmental-Duplications
cd /projappl/project_2007567/Genomes/T2T-CHM13/Segmental-Duplications

wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_SD.bed
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_SD.full.bed

#A more comprehensive centromere/satellite repeat annotation. (Re colored to be consistent with the primates Cen/Sat tracks)
https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_censat_v2.1.bed

RepeatMasker v4.1.2p1.2022Apr14 in bed or native out. 


https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed
https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.out

Here is a great resource for building a custom RepeatMasker library with new repeat models from the T2T genomes and a walk through for running RepeatMasker.
https://github.com/jessicaStorer88/RepeatMasker_library_CHM13

Composite repeats
https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_composite-repeats_2022DEC.bed

New satellites
https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_new-satellites_2022DEC.bed

chrXY sequence class, v1
https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_chrXY_sequence_class_v1.bed

TElomere
https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_telomere.bed


Y specific annotation

Palindromes and Inverted Repeats, v1 https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0Y_inverted_repeats_v1.bed
Amplicons v1 https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0Y_amplicons_v1.bed
AZFa, AZFb, AZFc and DYZ v1: https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0Y_AZF_DYZ_v1.bed


# Variant calls
# 1000 Genomes Project, recalled on T2T-CHM13v2.0. Now available for all chromosomes, for the entire 3,202 samples or the unrelated 2504 samples. Reference sets, bam, and vcf files are also available on AnVIL_T2T_CHRY.
# 1000 Genomes Project - Allele Frequency by Population, of the unrelated samples, further excluding 14 individuals discovered as first and second degree relatives (more details here).
# 1000 Genomes Project - Phased with SHAPEIT5, using the above variant calls.
# Simons Genome Diversity Project, recalled on T2T-CHM13v2.0. Reference sets, bam, and vcf files are also available on AnVIL_T2T_CHRY.
# gnomAD v3.1.2 from FTP: This is a lifted over version from GRCh38, annotated with predicted molecular consequence and transcript-specific variant deleteriousness scores from PolyPhen-2 and SIFT using Ensembl Variant Effect Predictor.
# Short-Read Accessibility Mask, with the three masks used to make the combined_mask are available here. See description
# ClinVar 20220313, lifted over from GRCh38. See description
# GWAS v1.0, lifted over from GRCh38. See description
# dbSNP build 155, lifted over from GRCh38. See description
# Variants disappearing in GRCh38-Y coordinates, v0.005 when using T2T-Y as a reference, more details are here.

#Liftover resources
#1:1 Liftover GRCh38 <-> T2T-CHM13v2.0, see description https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/v1_nflo_description.html
#GRCh38/hg38 -> T2T-CHM13v2.0: grch38-chm13v2.chain 
  https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/grch38-chm13v2.chain
#GRCh38/hg38 <- T2T-CHM13v2.0: chm13v2-grch38.chain 
https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/chm13v2-grch38.chain
# Alignment grch38-chm13v2.paf
# 1:1 Liftover hg19 <-> T2T-CHM13v2.0
# GRCh37/hg19 -> T2T-CHM13v2.0: hg19-chm13v2.chain
# GRCh37/hg19 <- T2T-CHM13v2.0: chm13v2-hg19.chain
# Alignment hg19-chm13v2.paf
