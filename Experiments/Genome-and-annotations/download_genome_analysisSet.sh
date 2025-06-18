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

cd /projappl/project_2007567/Genomes/GRCh38.p14_31102023


# The Dec. 2013  assembly of the human genome (GRCh38 Genome Reference Consortium
# Human Reference 38), is called hg38 at UCSC. This directory contains the genome
# as released by UCSC, selected annotation files and updates. The directory
# "genes/" contains GTF/GFF files for the main gene transcript sets.

# The sequences of the main chromosomes are identical to the genome files distributed
# by NCBI and the EBI, but the sequence names are different. For example, the
# name of chromosome 1 is called "chr1" at UCSC, "NC_000001.11" at NCBI, and "1"
# at the EBI.  Also, the lowercasing in the files is not exactly identical, as
# UCSC, NCBI and EBI run Repeatmasker with slightly different settings.
# 
# The NCBI accession of the UCSC hg38 genome is GCA_000001405.15. The version 
# that includes the updates for patch release 14 GRCh38.p14 has the NCBI
# accession GCA_000001405.29.

# The GRCh38 assembly contains more than just the chromosome sequences, 

#but also 
# a mitochondrial genome, #chrM	16569

#unplaced sequences, 

#centromeric sequences

# and alternates. 


# Sequence names
# ^^^^^^^^^^^^^^
# 
# For historical reasons, what UCSC calls "chr1", Ensembl calls "1" and NCBI
# calls "NC_000067.6". The sequences are identical though. To map between UCSC,
# Ensembl and NCBI names, use our table "chromAlias", available via our Table
# Browser or as file:
# https://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/chromAlias.txt.gz We
# also provide a Python command line tool to convert sequence names in the most
# common genomics file formats:
# http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/chromToUcsc
# 
# During genome assembly, reads are assembled into "contigs" (a few kbp long),
# which are then joined into longer "scaffolds" of a few hundred kbp. These are 
# finally placed, often manually e.g. with FISH assays, onto chromosomes.
# As a result, the hg38 genome sequence files contains different types of sequences:
# 
# Chromosomes: 
# - made from scaffolds placed onto chromosome locations, 95% of the genome file
# - format: chr{chromosome number or name} 
# - e.g. chr1 or chrX, chrM for the mitochondrial genome.
# 
# Unlocalized scaffolds: KEEP THESE
# - a sequence found in an assembly that is associated with a specific 
# chromosome but cannot be ordered or oriented on that chromosome. 
# - format: chr{chromosome number or name}_{sequence_accession}v{sequence_version}_random
# - e.g. chr17_GL000205v2_random

# chr11_KI270721v1_random
# chr22_KI270737v1_random
# chr17_KI270730v1_random
# chr1_KI270708v1_random
# chr22_KI270731v1_random
# chr2_KI270716v1_random
# chr3_GL000221v1_random
# chr2_KI270715v1_random
# chr22_KI270734v1_random
# chr14_KI270725v1_random
# chr1_KI270706v1_random
# chr1_KI270712v1_random
# chr9_KI270719v1_random
# chr22_KI270733v1_random
# chr22_KI270736v1_random
# chr17_GL000205v2_random
# chr16_KI270728v1_random
# chr14_GL000194v1_random
# chr14_KI270722v1_random
# chr14_GL000009v2_random
# chr4_GL000008v2_random
# chr14_GL000225v1_random
# chr17_KI270729v1_random
# chr1_KI270707v1_random
# chrY_KI270740v1_random
# chr9_KI270718v1_random
# chr14_KI270723v1_random
# chr9_KI270720v1_random
# chr14_KI270724v1_random
# chr9_KI270717v1_random
# chr1_KI270710v1_random
# chr1_KI270713v1_random
# chr22_KI270732v1_random
# chr1_KI270714v1_random
# chr1_KI270711v1_random
# chr22_KI270735v1_random
# chr14_KI270726v1_random
# chr15_KI270727v1_random
# chr1_KI270709v1_random
# chr22_KI270739v1_random
# chr5_GL000208v1_random
# chr22_KI270738v1_random

# 
# Unplaced scaffolds: KEEP THESE
# - a sequence found in an assembly that is not associated with any chromosome.  
# - format: chrUn_{sequence_accession}v{sequence_version}
# - e.g. chrUn_GL000220v1

# chrUn_KI270442v1
# chrUn_KI270743v1
# chrUn_KI270747v1
# chrUn_KI270742v1
# chrUn_GL000195v1
# chrUn_GL000224v1
# chrUn_GL000219v1
# chrUn_GL000216v2
# chrUn_KI270744v1
# chrUn_GL000213v1
# chrUn_GL000220v1
# chrUn_GL000218v1
# chrUn_KI270749v1
# chrUn_KI270741v1
# chrUn_KI270751v1
# chrUn_KI270750v1
# chrUn_KI270519v1
# chrUn_GL000214v1
# chrUn_KI270438v1
# chrUn_KI270748v1
# chrUn_KI270435v1
# chrUn_KI270538v1
# chrUn_KI270756v1
# chrUn_KI270757v1
# chrUn_KI270746v1
# chrUn_KI270753v1
# chrUn_KI270589v1
# chrUn_KI270745v1
# chrUn_KI270754v1
# chrUn_KI270317v1
# chrUn_KI270755v1
# chrUn_KI270579v1
# chrUn_KI270752v1
# chrUn_KI270512v1
# chrUn_KI270322v1
# chrUn_GL000226v1
# chrUn_KI270311v1
# chrUn_KI270366v1
# chrUn_KI270511v1
# chrUn_KI270448v1
# chrUn_KI270521v1
# chrUn_KI270581v1
# chrUn_KI270582v1
# chrUn_KI270515v1
# chrUn_KI270588v1
# chrUn_KI270591v1
# chrUn_KI270522v1
# chrUn_KI270507v1
# chrUn_KI270590v1
# chrUn_KI270584v1
# chrUn_KI270320v1
# chrUn_KI270382v1
# chrUn_KI270468v1
# chrUn_KI270467v1
# chrUn_KI270362v1
# chrUn_KI270517v1
# chrUn_KI270593v1
# chrUn_KI270528v1
# chrUn_KI270587v1
# chrUn_KI270364v1
# chrUn_KI270371v1
# chrUn_KI270333v1
# chrUn_KI270374v1
# chrUn_KI270411v1
# chrUn_KI270414v1
# chrUn_KI270510v1
# chrUn_KI270390v1
# chrUn_KI270375v1
# chrUn_KI270420v1
# chrUn_KI270509v1
# chrUn_KI270315v1
# chrUn_KI270302v1
# chrUn_KI270518v1
# chrUn_KI270530v1
# chrUn_KI270304v1
# chrUn_KI270418v1
# chrUn_KI270424v1
# chrUn_KI270417v1
# chrUn_KI270508v1
# chrUn_KI270303v1
# chrUn_KI270381v1
# chrUn_KI270529v1
# chrUn_KI270425v1
# chrUn_KI270396v1
# chrUn_KI270363v1
# chrUn_KI270386v1
# chrUn_KI270465v1
# chrUn_KI270383v1
# chrUn_KI270384v1
# chrUn_KI270330v1
# chrUn_KI270372v1
# chrUn_KI270548v1
# chrUn_KI270580v1
# chrUn_KI270387v1
# chrUn_KI270391v1
# chrUn_KI270305v1
# chrUn_KI270373v1
# chrUn_KI270422v1
# chrUn_KI270316v1
# chrUn_KI270338v1
# chrUn_KI270340v1
# chrUn_KI270583v1
# chrUn_KI270334v1
# chrUn_KI270429v1
# chrUn_KI270393v1
# chrUn_KI270516v1
# chrUn_KI270389v1
# chrUn_KI270466v1
# chrUn_KI270388v1
# chrUn_KI270544v1
# chrUn_KI270310v1
# chrUn_KI270412v1
# chrUn_KI270395v1
# chrUn_KI270376v1
# chrUn_KI270337v1
# chrUn_KI270335v1
# chrUn_KI270378v1
# chrUn_KI270379v1
# chrUn_KI270329v1
# chrUn_KI270419v1
# chrUn_KI270336v1
# chrUn_KI270312v1
# chrUn_KI270539v1
# chrUn_KI270385v1
# chrUn_KI270423v1
# chrUn_KI270392v1
# chrUn_KI270394v1

# 
# Alternate loci scaffolds: DO NOT INCLUDE THESE
# - a scaffold that provides an alternate representation of a locus found
#   in the primary assembly. These sequences do not represent a complete
#   chromosome sequence although there is no hard limit on the size of the
#   alternate locus; currently these are less than 1 Mb. These could either 
#   be NOVEL patch sequences, added through patch releases, or present in the 
#   initial assembly release.
# - format: chr{chromosome number or name}_{sequence_accession}v{sequence_version}_alt
# - e.g. chr6_GL000250v2_alt

# chr19_GL949751v2_alt
# chr19_KV575257v1_alt
# chr5_GL383530v1_alt
# chr22_KI270877v1_alt
# chr8_KZ559107v1_alt
# chr13_KI270843v1_alt
# chr18_GL383568v1_alt
# chr19_GL949748v2_alt
# chr19_GL949750v2_alt
# chr19_KI270938v1_alt
# chr11_KI270902v1_alt
# chr17_KI270859v1_alt
# chr19_GL949749v2_alt
# chr3_KI270783v1_alt
# chr12_KQ090023v1_alt
# chr1_KI270760v1_alt
# chr2_KI270768v1_alt
# chr1_GL383519v1_alt
# chr2_KI270771v1_alt
# chr7_KI270803v1_alt
# chr18_KI270864v1_alt
# chr4_KI270787v1_alt
# chr3_KI270781v1_alt
# chr5_KI270897v1_alt
# chr21_GL383581v2_alt
# chr17_KZ559114v1_alt
# chr20_KI270869v1_alt
# chr7_GL383534v2_alt
# chr12_KI270834v1_alt
# chr4_KI270785v1_alt
# chr2_KI270769v1_alt
# chr12_GL383549v1_alt
# chr19_KI270918v1_alt
# chr13_KQ090025v1_alt
# chr2_GL383522v1_alt
# chr6_GL383533v1_alt
# chr5_KI270793v1_alt
# chr7_KI270807v1_alt
# chr20_GL383577v2_alt
# chr5_KI270898v1_alt
# chr5_KI270795v1_alt
# chr8_KI270815v1_alt
# chr2_KI270772v1_alt
# chr17_GL383564v2_alt
# chr8_KI270819v1_alt
# chr9_KQ090019v1_alt
# chr16_KI270854v1_alt
# chr14_KI270846v1_alt
# chr5_KZ208910v1_alt
# chr2_KI270770v1_alt
# chr8_KI270901v1_alt
# chr17_KI270907v1_alt
# chr2_KI270775v1_alt
# chr12_GL383552v1_alt
# chr5_KN196477v1_alt
# chr1_KZ208905v1_alt
# chr2_KZ208908v1_alt
# chr1_KQ458382v1_alt
# chr8_KI270814v1_alt
# chr17_KI270908v1_alt
# chr2_GL383521v1_alt
# chr21_KI270873v1_alt
# chrX_KI270881v1_alt
# chr22_KQ759761v1_alt
# chr8_KI270818v1_alt
# chr19_KV575260v1_alt
# chr14_KI270847v1_alt
# chr6_KI270799v1_alt
# chr12_GL383553v2_alt
# chr22_KN196486v1_alt
# chr12_KZ559112v1_alt
# chr11_GL383547v1_alt
# chr19_KI270888v1_alt
# chr19_GL383574v1_alt
# chr22_KQ458387v1_alt
# chr22_KN196485v1_alt
# chr19_KV575258v1_alt
# chr19_KI270884v1_alt
# chr17_KI270910v1_alt
# chr18_KI270911v1_alt
# chr7_KI270804v1_alt
# chr7_KI270806v1_alt
# chr4_KI270788v1_alt
# chr8_KI270817v1_alt
# chr19_KV575251v1_alt
# chr18_GL383572v1_alt
# chr19_KV575255v1_alt
# chr2_KI270893v1_alt
# chr5_GL339449v2_alt
# chr2_KI270767v1_alt
# chr1_KI270892v1_alt
# chr3_KI270782v1_alt
# chr22_GL383582v2_alt
# chr3_KI270895v1_alt
# chr9_GL383539v1_alt
# chr18_KZ559116v1_alt
# chr3_KI270934v1_alt
# chr4_KQ090014v1_alt
# chr9_KQ090018v1_alt
# chr19_KV575246v1_alt
# chr3_KZ559101v1_alt
# chr3_KI270936v1_alt
# chr4_GL383527v1_alt
# chr5_KI270794v1_alt
# chr18_GL383570v1_alt
# chr3_KI270937v1_alt
# chr1_KI270761v1_alt
# chr1_KZ208904v1_alt
# chr3_KI270924v1_alt
# chr19_KV575253v1_alt
# chr21_KI270874v1_alt
# chr12_GL877875v1_alt
# chr18_GL383569v1_alt
# chr18_KI270863v1_alt
# chr19_KV575248v1_alt
# chr13_KQ090024v1_alt
# chr13_KI270841v1_alt
# chr16_KQ031390v1_alt
# chr12_GL383550v2_alt
# chr19_KI270931v1_alt
# chr19_KV575247v1_alt
# chr19_GL383575v2_alt
# chr19_KI270883v1_alt
# chr19_KI270933v1_alt
# chr19_KI270915v1_alt
# chr19_KI270891v1_alt
# chr19_KI270889v1_alt
# chr19_KI270919v1_alt
# chr19_KI270885v1_alt
# chr19_KV575259v1_alt
# chr9_GL383541v1_alt
# chr14_KZ208919v1_alt
# chr7_KZ559106v1_alt
# chr5_KI270796v1_alt
# chr3_JH636055v2_alt
# chr5_GL383531v1_alt
# chr3_KI270777v1_alt
# chr18_KI270912v1_alt
# chr2_KI270776v1_alt
# chr22_KQ458388v1_alt
# chr12_KZ208918v1_alt
# chr6_KI270800v1_alt
# chr3_KZ208909v1_alt
# chr22_KI270928v1_alt
# chr11_KI270830v1_alt
# chr19_GL000209v2_alt
# chr19_KV575252v1_alt
# chr17_KI270860v1_alt
# chr5_KI270792v1_alt
# chr10_GL383545v1_alt
# chr13_KI270839v1_alt
# chr3_GL383526v1_alt
# chr14_KI270845v1_alt
# chr11_KZ559111v1_alt
# chr10_KI270824v1_alt
# chr2_KZ208907v1_alt
# chr17_GL000258v2_alt
# chr1_GL383518v1_alt
# chr20_KI270870v1_alt
# chr12_GL383551v1_alt
# chr3_KI270784v1_alt
# chr19_KI270890v1_alt
# chr19_KI270916v1_alt
# chr1_KI270765v1_alt
# chr10_KQ090020v1_alt
# chr6_KB021644v2_alt
# chr11_KI270826v1_alt
# chr19_KI270929v1_alt
# chr22_KI270878v1_alt
# chr19_KI270922v1_alt
# chrX_KV766199v1_alt
# chr19_GL383576v1_alt
# chr10_KI270825v1_alt
# chr19_KI270923v1_alt
# chr7_KI270899v1_alt
# chr19_KI270917v1_alt
# chr11_JH159137v1_alt
# chr13_KI270840v1_alt
# chr16_GL383556v1_alt
# chr3_KZ559105v1_alt
# chr5_KI270791v1_alt
# chr15_KI270906v1_alt
# chr17_KI270861v1_alt
# chr3_KI270935v1_alt
# chr6_KI270797v1_alt
# chr3_KZ559102v1_alt
# chr19_KI270920v1_alt
# chr18_GL383571v1_alt
# chr19_KI270930v1_alt
# chr11_JH159136v1_alt
# chr21_GL383579v2_alt
# chr11_KN538368v1_alt
# chr15_MU273375v1_alt
# chr11_KI270829v1_alt
# chr19_KI270886v1_alt
# chr18_KQ458385v1_alt
# chr19_KI270914v1_alt
# chr3_KI270779v1_alt
# chr4_KQ983258v1_alt
# chr4_KI270789v1_alt
# chr19_KI270887v1_alt
# chr7_KI270809v1_alt
# chr7_KI270805v1_alt
# chr11_KI270832v1_alt
# chr1_KQ458384v1_alt
# chr2_KI270894v1_alt
# chr11_KI270903v1_alt
# chr3_ML143343v1_alt
# chr19_KI270932v1_alt
# chr11_KI270927v1_alt
# chr4_KI270790v1_alt
# chr19_KV575256v1_alt
# chr2_KI270774v1_alt
# chr17_GL383565v1_alt
# chr3_KI270780v1_alt
# chr5_GL949742v1_alt
# chr8_KI270926v1_alt
# chr16_KI270855v1_alt
# chr19_KI270867v1_alt
# chr17_KI270858v1_alt
# chr4_KQ090015v1_alt
# chr15_KQ031389v1_alt
# chr12_KI270835v1_alt
# chr19_KV575250v1_alt
# chr4_KI270786v1_alt
# chr15_KI270849v1_alt
# chr17_KV766197v1_alt
# chr3_KI270778v1_alt
# chr19_KI270882v1_alt
# chr1_KI270766v1_alt
# chr22_KI270875v1_alt
# chr11_MU273368v1_alt
# chr15_KI270851v1_alt
# chr22_KI270876v1_alt
# chr14_ML143368v1_alt
# chr16_KI270853v1_alt
# chr16_KQ090027v1_alt
# chr7_KI270808v1_alt
# chr6_KI270798v1_alt
# chrX_KI270913v1_alt
# chr17_KV766198v1_alt
# chr17_JH159146v1_alt
# chr1_KQ983255v1_alt
# chr19_KI270921v1_alt
# chr8_KI270812v1_alt
# chrX_KI270880v1_alt
# chr2_MU273340v1_alt
# chr17_KI270857v1_alt
# chr18_GL383567v1_alt
# chr8_KI270811v1_alt
# chr19_KV575249v1_alt
# chrX_MU273396v1_alt
# chr15_GL383554v1_alt
# chr11_KI270831v1_alt
# chr8_KI270813v1_alt
# chr11_KZ559110v1_alt
# chr5_MU273356v1_alt
# chr3_KZ559103v1_alt
# chr22_KI270879v1_alt
# chr8_KI270816v1_alt
# chr13_KI270838v1_alt
# chr4_MU273349v1_alt
# chr10_GL383546v1_alt
# chr8_KI270900v1_alt
# chr14_KI270844v1_alt
# chr17_KI270909v1_alt
# chr15_KI270848v1_alt
# chrX_MU273397v1_alt
# chr1_MU273332v1_alt
# chr1_KQ458383v1_alt
# chr1_KI270762v1_alt
# chr5_KV575243v1_alt
# chr8_KI270820v1_alt
# chr1_GL383520v2_alt
# chr17_MU273378v1_alt
# chr13_KI270842v1_alt
# chr8_KI270810v1_alt
# chr17_GL383563v3_alt
# chr4_GL383528v1_alt
# chr4_KI270896v1_alt
# chr6_MU273357v1_alt
# chr19_GL383573v1_alt
# chr15_GL383555v2_alt
# chr17_KI270862v1_alt
# chr12_KI270837v1_alt
# chr12_GL877876v1_alt
# chr4_KV766193v1_alt
# chr1_KI270759v1_alt
# chr15_KI270850v1_alt
# chr19_KI270866v1_alt
# chr2_MU273337v1_alt
# chr9_KI270823v1_alt
# chr6_GL000252v2_alt
# chr6_GL000255v2_alt
# chr7_MU273358v1_alt
# chr6_GL000250v2_alt
# chr6_GL000253v2_alt
# chr15_KI270852v1_alt
# chr6_GL000251v2_alt
# chr6_GL000254v2_alt
# chr6_GL000256v2_alt
# chr2_MU273339v1_alt
# chr1_KI270764v1_alt
# chr15_KI270905v1_alt
# chr1_MU273330v1_alt
# chr19_KI270865v1_alt
# chr2_KQ983256v1_alt
# chr2_MU273338v1_alt
# chr1_KV880763v1_alt
# chr4_KI270925v1_alt
# chr12_KI270836v1_alt
# chr12_KI270904v1_alt
# chr4_GL000257v2_alt
# chr20_KI270871v1_alt
# chr16_KQ090026v1_alt
# chr9_GL383542v1_alt
# chr19_KI270868v1_alt
# chrX_MU273395v1_alt
# chr8_KI270822v1_alt
# chr21_GL383578v2_alt
# chr16_KI270856v1_alt
# chr11_KI270827v1_alt
# chr7_KZ208913v1_alt
# chr17_JH159147v1_alt
# chr2_KI270773v1_alt
# chr9_GL383540v1_alt
# chr19_GL949747v2_alt
# chr22_KB663609v1_alt
# chr21_GL383580v2_alt
# chr6_KI270802v1_alt
# chr12_KI270833v1_alt
# chr6_KI270758v1_alt
# chr16_KZ208921v1_alt
# chr19_GL949753v2_alt
# chr6_KQ090017v1_alt
# chr21_KI270872v1_alt
# chr5_GL383532v1_alt
# chr1_MU273331v1_alt
# chr6_KI270801v1_alt
# chr17_JH159148v1_alt
# chr19_MU273387v1_alt
# chr16_GL383557v1_alt
# chr17_GL383566v1_alt
# chr4_KQ090013v1_alt
# chr1_KI270763v1_alt
# chr2_GL582966v2_alt
# chr22_GL383583v2_alt
# chr8_KI270821v1_alt
# chr19_GL949752v1_alt
# chr19_GL949746v1_alt
# chr19_KV575254v1_alt

# 
# Fix loci scaffolds: DO NOT INCLUDE THESE
# - a patch that corrects sequence or reduces an assembly gap in a given
#   major release. FIX patch sequences are meant to be incorporated into
#   the primary or existing alt-loci assembly units at the next major
#   release.
# - these sequences are not part of the files in the initial/ directory
# - format: chr{chromosome number or name}_{sequence_accession}v{sequence_version}_fix
# - e.g. chr2_KN538362v1_fix

# chr22_KQ759762v1_fix
# chr22_KQ759762v2_fix
# chrY_KN196487v1_fix
# chr21_MU273391v1_fix
# chr12_MU273372v1_fix
# chr12_KZ208916v1_fix
# chr3_KZ559104v1_fix
# chr8_MU273361v1_fix
# chr11_KN196481v1_fix
# chr4_MU273350v1_fix
# chr15_MU273374v1_fix
# chr2_MU273341v1_fix
# chr1_KN196474v1_fix
# chr11_MU273371v1_fix
# chr22_ML143379v1_fix
# chr4_ML143348v1_fix
# chr19_MU273385v1_fix
# chr17_ML143374v1_fix
# chr6_KV766194v1_fix
# chrX_MU273394v1_fix
# chr11_KV766195v1_fix
# chr7_KV880764v1_fix
# chr10_KN538365v1_fix
# chr17_MU273381v1_fix
# chr2_ML143341v1_fix
# chrX_ML143384v1_fix
# chr13_KN538373v1_fix
# chr8_MU273359v1_fix
# chr17_KV575245v1_fix
# chr8_KV880766v1_fix
# chr1_MU273333v1_fix
# chr13_ML143364v1_fix
# chr8_KZ208914v1_fix
# chr11_ML143357v1_fix
# chr3_KQ031386v1_fix
# chr1_KN196473v1_fix
# chr11_ML143360v1_fix
# chr17_MU273383v1_fix
# chrX_ML143385v1_fix
# chr2_MU273345v1_fix
# chr4_ML143347v1_fix
# chr7_KQ031388v1_fix
# chr11_KQ090022v1_fix
# chr1_KN196472v1_fix
# chr17_MU273382v1_fix
# chr21_MU273392v1_fix
# chr12_ML143362v1_fix
# chr16_KV880768v1_fix
# chr10_MU273367v1_fix
# chr11_KQ759759v1_fix
# chr11_KQ759759v2_fix
# chr4_MU273351v1_fix
# chr13_KN538371v1_fix
# chr8_MU273363v1_fix
# chr2_KN538362v1_fix
# chr5_MU273353v1_fix
# chrY_KZ208924v1_fix
# chr5_MU273354v1_fix
# chr1_MU273334v1_fix
# chr12_KN196482v1_fix
# chr1_MU273335v1_fix
# chr11_ML143359v1_fix
# chr19_MU273386v1_fix
# chr4_KQ983257v1_fix
# chr18_KZ559115v1_fix
# chr17_MU273379v1_fix
# chr4_ML143344v1_fix
# chr6_KZ208911v1_fix
# chr2_MU273344v1_fix
# chr6_KQ090016v1_fix
# chr1_MU273336v1_fix
# chr9_ML143353v1_fix
# chr7_ML143352v1_fix
# chr10_KQ090021v1_fix
# chr8_KV880767v1_fix
# chr6_KN196478v1_fix
# chr11_ML143358v1_fix
# chr16_ML143373v1_fix
# chr20_MU273388v1_fix
# chr4_ML143349v1_fix
# chr10_KN196480v1_fix
# chr11_KZ559109v1_fix
# chr17_KV766196v1_fix
# chrX_ML143382v1_fix
# chr10_ML143355v1_fix
# chr12_ML143361v1_fix
# chr3_MU273347v1_fix
# chr11_KZ559108v1_fix
# chr1_KN538361v1_fix
# chr3_KN196476v1_fix
# chr12_KQ759760v1_fix
# chr6_KQ031387v1_fix
# chr1_KZ208906v1_fix
# chr9_KN196479v1_fix
# chr19_MU273384v1_fix
# chr16_MU273377v1_fix
# chr21_MU273390v1_fix
# chr9_MU273364v1_fix
# chr4_ML143345v1_fix
# chr5_MU273352v1_fix
# chr11_MU273370v1_fix
# chr13_KN196483v1_fix
# chr20_MU273389v1_fix
# chr13_KN538372v1_fix
# chr2_KN538363v1_fix
# chr15_ML143370v1_fix
# chr19_KN196484v1_fix
# chr3_KQ031385v1_fix
# chr8_MU273360v1_fix
# chr15_ML143372v1_fix
# chr14_ML143367v1_fix
# chrX_ML143381v1_fix
# chr19_KQ458386v1_fix
# chr18_KQ090028v1_fix
# chr13_ML143366v1_fix
# chr3_KV766192v1_fix
# chr22_ML143380v1_fix
# chr3_KN538364v1_fix
# chr10_KN538367v1_fix
# chr8_MU273362v1_fix
# chr11_MU273369v1_fix
# chr1_KZ559100v1_fix
# chr3_KN196475v1_fix
# chr11_ML143356v1_fix
# chr10_ML143354v1_fix
# chr1_KN538360v1_fix
# chr22_ML143378v1_fix
# chr1_KQ031383v1_fix
# chr7_KV880765v1_fix
# chr3_MU273346v1_fix
# chr3_MU273348v1_fix
# chr16_KZ559113v1_fix
# chr2_KQ031384v1_fix
# chr9_MU273365v1_fix
# chrY_KZ208923v1_fix
# chr2_MU273343v1_fix
# chr19_ML143376v1_fix
# chr15_KN538374v1_fix
# chr5_MU273355v1_fix
# chr21_ML143377v1_fix
# chr4_ML143346v1_fix
# chr17_MU273380v1_fix
# chr12_KN538369v1_fix
# chr15_ML143371v1_fix
# chr17_ML143375v1_fix
# chr9_MU273366v1_fix
# chr7_KZ208912v1_fix
# chr8_KZ208915v1_fix
# chr12_KZ208917v1_fix
# chr13_ML143365v1_fix
# chr5_KV575244v1_fix
# chrX_ML143383v1_fix
# chrX_MU273393v1_fix
# chr14_KZ208920v1_fix
# chr14_MU273373v1_fix
# chr13_ML143363v1_fix
# chr6_ML143351v1_fix
# chr2_ML143342v1_fix
# chr10_KN538366v1_fix
# chr12_KN538370v1_fix
# chrY_MU273398v1_fix
# chr16_MU273376v1_fix
# chr5_ML143350v1_fix
# chr18_KZ208922v1_fix
# chr2_MU273342v1_fix
# chr15_ML143369v1_fix


#the UCSC hg38 file. hg38.fa.gz - "Soft-masked" assembly sequence in one file.
#Repeats from RepeatMasker and Tandem Repeats Finder (with period of 12 or
#less) are shown in lower case; non-repeating sequence is shown in upper case. 
#(again, the most current version of this file is latest/hg38.fa.gz)

#The "latest/" symbolic link points to the subdirectory for the most recent
#patch version. (patch release 14).

# To unpack the *.tar.gz files:
#     tar xvzf <file>.tar.gz
# To uncompress the fa.gz files:

#rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz . 
#gunzip hg38.fa.gz. results in hg38.fa
#grep ">" hg38.fa | wc -l #771

#rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.chromFa.tar.gz .
#tar xvzf hg38.chromFa.tar.gz #creates chroms/*.fa
#ls chroms/ | wc -l #771




#Analysis set, alternate sequences removed:
rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz .
gunzip hg38.analysisSet.fa.gz

#grep ">" hg38.analysisSet.fa | wc -l #195





# Analysis set
# ^^^^^^^^^^^^
# 
# The GRCh38 assembly contains more than just the chromosome sequences, 

#but also 
# a mitochondrial genome, 
#chrM	16569

#unplaced sequences, 

#centromeric sequences
# and alternates. 

#To better capture variation in the human genome across the world
# it contains more copies of some loci than hg19. Some of these additions, like
# the EBV genome, are mostly relevant for genomic analysis, i.e. alignment.
# For an overview of the different types and reasons for the additions see
# https://software.broadinstitute.org/gatk/documentation/article?id=11010
# 
# This means that if you want to use the genome sequence for alignment and
# especially for variant calling, you should use the optimal genome file for your
# aligner. The genome file can make a big difference, especially for variant
# calling.  In most cases, the authors of your alignment program will provide
# advice on which hg38 genome version to use and usually they recommend one of
# the files in our analysisSet/ directory, like the GATK link above. These 
# special genome files sometimes remove the alternate sequences, sometimes they
# add decoys or change single nucleotides towards the major allele, but they never
# insert or delete sequences, so the annotation coordinates remain the same.
# 
# - for BWA see also https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use
# - for Novoalign see its manual at http://www.novocraft.com/userfiles/file/Novocraft.pdf
# - For Bowtie, see the different versions of the human genome that the Bowtie authors
#   provide: http://bowtie-bio.sourceforge.net/index.shtml
# 
# Also see analysisSet/README.txt for further details
# 
# Patches
# ^^^^^^^

# Like hg19, hg38 has been updated with patches since its release in 2013. GRC
# patch releases do not change any previously existing sequences; they simply add
# small, new sequences for fix patches or alternate haplotypes that correspond to
# specific regions of the main chromosome sequences (see below). For most users,
# the patches are unlikely to make a difference and may complicate the analysis
# as they introduce more duplication. If you want a version of the genome
# without these complexities, look at the analysisSet/ subdirectory.
# 
# The initial/ subdirectory contains files for the initial release of GRCh38,
# which includes the original alternate sequences (261) and no fix sequences.
# 
# The p11/ subdirectory contains files for GRCh38.p11 (patch release 11).
# 
# The p12/ subdirectory contains files for GRCh38.p12 (patch release 12).
# 
# The p13/ subdirectory contains files for GRCh38.p13 (patch release 13).
# 
# The p14/ subdirectory contains files for GRCh38.p14 (patch release 14).

# Files
# ^^^^^
# 
# Files in this directory reflect the initial 2013 release of the genome, 
# the most current versions are in the "latest/" subdirectory:
# 
# 
# 
# hg38.2bit - contains the complete human/hg38 genome sequence
#     in the 2bit file format.  Repeats from RepeatMasker and Tandem Repeats
#     Finder (with period of 12 or less) are shown in lower case; non-repeating
#     sequence is shown in upper case.  The utility program, twoBitToFa (available
#     from the kent src tree), can be used to extract .fa file(s) from
#     this file.  A pre-compiled version of the command line tool can be
#     found at:
#         http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
#     See also:
#         http://genome.ucsc.edu/admin/git.html
# 	http://genome.ucsc.edu/admin/jk-install.html

# hg38.agp.gz - Description of how the assembly was generated from
#     fragments.



# hg38.chromFaMasked.tar.gz - The assembly sequence in one file per chromosome.
#     Repeats are masked by capital Ns; non-repeating sequence is shown in
#     upper case.

# hg38.fa.masked.gz - "Hard-masked" assembly sequence in one file.
#     Repeats are masked by capital Ns; non-repeating sequence is shown in
#     upper case.

# hg38.fa.out.gz - RepeatMasker .out file.  RepeatMasker was run with the
#     -s (sensitive) setting.
#     June 20 2013 (open-4-0-3) version of RepeatMasker
#     RepBase library: RELEASE 20130422

rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.out.gz .
gunzip hg38.fa.out.gz



# hg38.fa.align.gz - RepeatMasker .align file.  RepeatMasker was run with the
#     -s (sensitive) setting.
#     June 20 2013 (open-4-0-3) version of RepeatMasker
#     RepBase library: RELEASE 20130422

rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.align.gz .
gunzip hg38.fa.align.gz

# hg38.trf.bed.gz - Tandem Repeats Finder locations, filtered to keep repeats
#     with period less than or equal to 12, and translated into UCSC's BED
#     format.

rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.trf.bed.gz .
gunzip hg38.trf.bed.gz

# md5sum.txt - checksums of files in this directory

# mrna.fa.gz - Human mRNA from GenBank. This sequence data is updated
#     regularly via automatic GenBank updates.

# refMrna.fa.gz - RefSeq mRNA from the same species as the genome.
#     This sequence data is updated regularly via automatic GenBank
#     updates.

# upstream1000.fa.gz - Sequences 1000 bases upstream of annotated
#     transcription starts of RefSeq genes with annotated 5' UTRs.
#     This file is updated regularly. It might be slightly out of sync with
#     the RefSeq data shown on the browser, as is it updated daily for most assemblies.

# upstream2000.fa.gz - Same as upstream1000, but 2000 bases.
# 
# upstream5000.fa.gz - Same as upstream1000, but 5000 bases.
# 
# xenoMrna.fa.gz - GenBank mRNAs from species other than that of
#     the genome. This sequence data is updated regularly via
#     automatic GenBank updates.


# Two-column tab-separated text file containing assembly
#     sequence names and sizes.
    
#rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.chrom.sizes .

# hg38.gc5Base.wigVarStep.gz - ascii data wiggle variable step values used
#                            - to construct the GC Percent track
# hg38.gc5Base.bw - binary bigWig data for the gc5Base track.


 # sequence name alias file, one line
 #    for each sequence name.  First column is sequence name followed by
 #    tab separated alias names.
    
#rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.chromAlias.txt .

# hg38.chromAlias.bb - bigBed file for alias sequence names, one line
#     for each sequence name. The first three columns are the sequence in 
#     BED format, followed by tab-separated alias names.
#     The .bb file is used by bedToBigBed as a URL to avoid having to download
#     the entire chromAlias.txt file.  From the usage message:
#         -sizesIsChromAliasBb -- If set, then chrom.sizes file is assumed to be a
#         chromAlias bigBed file or a URL to a such a file (see above).

# More documentation is found here:
# https://genomewiki.ucsc.edu/index.php?title=Chrom_Alias

# Dropped in Genbank and Refseq official releases patch14 since these 2 old versions are obsolete and no longer needed.
# This patch contains their v2 replacements.
# chr11_KQ759759v1_fix
# chr22_KQ759762v1_fix

# Dropped in Refseq official release Patch14.
# These 3 are contamination or obsolete.
# chr10_KI270825v1_alt
# chr22_KI270734v1_random
# chr11_KI270721v1_random

# Because of the difficulty of removing the old chroms chr11_KQ759759v1_fix and chr22_KQ759762v1_fix from all of the database tables and bigData files,
# custom tracks, and hubs, we are not dropping them from the UCSC hg38 patch 14 .2bit and chromInfo.
# However, we have dropped them from chromAlias to accord with the Genbank and Refseq official releases for patch14.
# 
# 
# Dropped in Patch13 from Refseq
# chrUn_KI270752v1      HSCHRUN_RANDOM_CTG29    KI270752.1
# KI270752.1 is no longer part of the RefSeq assembly because it is hamster sequence
# derived from the human-hamster CHO cell line. 
# https://www.ncbi.nlm.nih.gov/grc/human/issues/HG-2587


# ------------------------------------------------------------------
# How to Download
# ^^^^^^^^^^^^^^^
# 
# If you plan to download a large file or multiple files from this
# directory, we recommend that you use ftp rather than downloading the
# files via our website. To do so, ftp to hgdownload.cse.ucsc.edu
# [username: anonymous, password: your email address], then cd to the
# directory goldenPath/hg38/bigZips. To download multiple files, use
# the "mget" command:
# 
#     mget <filename1> <filename2> ...
#     - or -
#     mget -a (to download all the files in the directory)
# 
# Alternate methods to ftp access.
# 
# Using an rsync command to download the entire directory:
#     rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/ .
# For a single file, e.g. chromFa.tar.gz
#     rsync -avzP 
#         rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/chromFa.tar.gz .
# 
# Or with wget, all files:
#     wget --timestamping 
#         'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/*'
# With wget, a single file:
#     wget --timestamping 
#         'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/chromFa.tar.gz' 
#         -O chromFa.tar.gz
# 
# To unpack the *.tar.gz files:
#     tar xvzf <file>.tar.gz
# To uncompress the fa.gz files:
#     gunzip <file>.fa.gz

