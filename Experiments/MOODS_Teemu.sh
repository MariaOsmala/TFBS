#!/bin/bash

#!/bin/bash
array=$1 #this varies between 0 and 399 when 3982 motifs
#When 3294 motifs varies between 0 and 329



export PATH="/projappl/project_2006203/softwares/conda_envs/MOODS/bin:$PATH"


#S=/projappl/project_2006203/Genomes/Homo_sapiens/chr_sequences.fa #What is this?, this is repeat masked, contains only chromsomes chr1-22, X, Y
#S=/projappl/project_2006203/Genomes/Homo_sapiens/hg38.fa.masked #Repeats (RepeatMasker& Tandem Repeats Finder) masked by capital N, contains also other genomes
#S=/projappl/project_2006203/Genomes/Homo_sapiens/hg38.fa

#USE THIS INSTEAD?
S=/projappl/project_2007567/Genomes/GRCh38.p14_31102023/hg38.analysisSet.fa

n=300000 #the number of best hits

outfolder=/scratch/project_2006203/TFBS/Results/MOODS_human_final_whole_genome/
mkdir $outfolder
#outfile=../Results/MOODS/MOODS_"$array".csv"

#rename Yimengs motifs, add 

#filenames=($(ls ../PWMs/fromYimeng/pwms_space/pfm_composite_new/*.pfm))
#for filename in "${filenames[@]}"
#do
#   newname=${filename%.pfm}"_composite.pfm"
#   echo $newname
#   mv $filename $newname
   
#done

#filenames=($(ls ../PWMs/fromYimeng/pwms_space/pfm_spacing_new/*.pfm))
#for filename in "${filenames[@]}"
#do
#   newname=${filename%.pfm}"_spacing.pfm"
#   echo $newname
#   mv $filename $newname
   
#done



pwms=($(ls ../PWMs_final/*/pwms_space/*/*.pfm)) #3294




#echo "${pwms[3]##*/}"

#echo "${#pwms[@]}"
#echo "${pwms[@]}"

nro_pwms=${#pwms[@]} #3294


start_ind=$(($array*10)) #100
end_ind=$((($array+1)*10 -1)) #100
length=10 #100

#index is from 0 to 3293

if [[ $end_ind -gt $(($nro_pwms-1)) ]] #3309 > 3293
then
     echo $end_ind is greater than $(($nro_pwms-1))
     length=$(($nro_pwms-$start_ind+1))
fi

#threshold selection (exactly one required):
#  -p p, --p-value p     compute threshold from p-value
#  -t T, --threshold T   use specified absolute threshold
#  -B n, --best-hits n   return at least the specified amount of best matches

# search and model behaviour (optional):
#  --no-snps             ignore IUPAC symbols coding multiple nucleotides
#  --threshold-precision x
#                        specify the precision used for computing the
#                        thresholds from p-values (default = 2000.0)

#Write this directly to a database

#--p-value 1e-4 --lo-bg 0.2977 0.2023 0.2023 0.2977

#moods-dna.py -m ${pwms[@]:$start_ind:$length}  -s $S -B $n -o $outfile #${array[@]:START:LENGTH}

#removed -B $n

#by default, MOODS assume that the threshold is given by a p-value x, 
#and the actual threshold T is chosen so that the probability that 
#the background distribution π generates 
#a sequence u of length m with score W_L(u)>=T is p.
# p-value 0.0001

#-m to give count/frequency matrices, will be converted to PWMs
export PATH="/projappl/project_2006203/softwares/conda_envs/MOODS/bin:$PATH"

#--bg pA pC pG pT      background distribution for computing thresholds from
#                        p-value with --batch (default is 0.25 for all alleles)
#--ps p                total pseudocount added to each matrix column in log-
#                        odds conversion (default = 0.01)

#  --log-base x          logarithm base for log-odds conversion (default
#                       natural logarithm)

#  --lo-bg pA pC pG pT   background distribution for log-odds conversion
#                        (default is 0.25 for all alleles)

#Teemu: Yleensä olen hakenut löysällä affineettirajalla ja 
#valinnut sitten kärjestä halutun määrän. Suuria määriä hakiessa 
#jokin minimiaffiniteettiraja (esim. 9) on kuitenkin ollut käytössä, 
#ettei todella pitkille motiiveille väkisin etsitä ihan järjettömiä osumia.


#raja on 9 on alunperin EEL:istä, jonka affiniteetit ovat log2 kun taas moodsissa oletuskantaluku on e (--log-base). 

#Jokaiselle motiiville etsitään 300 000 matchiä?

#keskeinen ongelma tässä on se, että mielivaltaisella matriisilla ei 
#välttämättä ole mitään pisterajaa, jolla se antaa täsmälleen 300 000 matchia, 
#tai edes lähelle sitä. Tämä on erityisesti ongelma matriiseilla j
#otka ovat melko tasaisesti jakautuneita joissain positioissa.

#MOODSin -B parametri muistaakseni yrittää etsiä jonkun rajan jolla tulee vähintään 
#k matchia ja enintään 2*k matchia, ja jos tämä ei ole mahdollista, 
#se valitsee pienimmän rajan jolla tulee yli k matchia. Valitettavasti tämä 
#pakolti tarkoittaa sitä, että joissain tapauksissa saat kasan matcheja 
#joista isolla osalla on sama score.

#Vastaavasti --p-value P antaa noin P*m matchia, kun sekvenssin pituus on m, 
#mutta vain tilastollisessa mielessä eli olettaen että sekvenssi on oikeasti 
#jakautunut taustajakauman mukaan ja on pitkä. 
#Käytännössä matchien määrä on tietysti aika vaihteleva.

#Eli käytännössä joillekin matriiseille voi löytyä esim. 1000 hyvää hittiä, 
#mutta sitten seuraavaksi parhailla 4 miljoonalla hitillä on kaikilla sama pistemäärä. 
#Näissä tapauksissa filtteröinti on sitten tehtävä jollain muulla kriteerillä. 


#Olikos niin, että jos ei aseta tuota taustajakaumaa (--lo-bg) tausta estimoidaan genomisekvenssistä?

#Joo, juurikin näin. Eli --bg asettaa p-arvojen laskennassa käytettävän taustajakauman, ja 
#--lo-bg asettaa log-odds-matriisien laskennassa käytettävän taustajakauman. 
#Jos jompaa kumpaa ei anneta, tausta kyseiseen tarkoitukseen estimoidaan aina jokaisen sekvenssin kohdalla erikseen kyseisestä sekvenssistä.

#Yleensä ajattelen, että --lo-bg on loogista antaa parametrina, 
#koska motiivien jakaumaa verrataan johonkin fiksattuun prosessiin 
#joka generoi taustajakauman (en tosin ole biologi, joten en tiedä kannattaako minuun luottaa tässä.)


moods-dna.py -m ${pwms[@]:$start_ind:$length} --threshold 2 -s $S | gzip > $LOCAL_SCRATCH"/MOODS_"$array".csv.gz" 

cd $LOCAL_SCRATCH
cp MOODS_"$array".csv.gz $outfolder 

#gzip $outfile

#outfile to parquet, where does the file go, requires a lot of space
#csvcli $outfile convert -to "parquet"


#FATAL:   container creation failed: hook function for tag layer returns error: failed to create /tmp/nvme/job_14238673 directory: mkdir /tmp/nvme/job_14238673: permission denied

#/projappl/project_2006203/softwares/conda_envs/MOODS/bin/moods-dna.py: line 22: 705506 Segmentation fault
#   /usr/bin/singularity --silent exec -B $DIR/../$SQFS_IMAGE:$INSTALLATION_PATH:image-src=/
# $DIR/../$CONTAINER_IMAGE bash -c "eval \"\$(/CSC_CONTAINER/miniconda/bin/conda shell.bash hook )\"
# && conda activate env1 &>/dev/null &&  exec -a $_O_SOURCE $DIR/moods-dna.py $(
#test $# -eq 0 || printf " %q" "$@" )"
