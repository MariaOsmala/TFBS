#!/bin/bash

#!/bin/bash

array=$1 #this varies between 0 and 102
#the last start 1020 and end 1029 (-> 1023)

#For how many motifs you need to run your analysis
#wc -l rerunMOODS_lower_threshold.txt #1024 #Are these representatives?

#files=($(awk -F '"' '{print $4}' rerunMOODS_lower_threshold.txt)) #1024
#files=($(awk -F '"' '{print $4}' rerunMOODS_lower_threshold_4.txt)) #761
files=($(awk -F '"' '{print $4}' rerunMOODS_lower_threshold_3.txt)) #578

#PBNET

export PATH="/projappl/project_2006203/softwares/conda_envs/MOODS/bin:$PATH"


S=/projappl/project_2006203/Genomes/Homo_sapiens/chr_sequences.fa #What is this?, this is repeat masked, contains only chromsomes chr1-22, X, Y
#S=/projappl/project_2006203/Genomes/Homo_sapiens/hg38.fa.masked #Repeats (RepeatMasker& Tandem Repeats Finder) masked by capital N, contains also other genomes
#S=/projappl/project_2006203/Genomes/Homo_sapiens/hg38.fa

n=300000 #the number of best hits
#mkdir /scratch/project_2006203/TFBS/Results/MOODS_Teemu_4/
#mkdir /scratch/project_2006203/TFBS/Results/MOODS_Teemu_3/
mkdir /scratch/project_2006203/TFBS/Results/MOODS_Teemu_2/
#outfile="/scratch/project_2006203/TFBS/Results/MOODS_Teemu_3/MOODS_"$array".csv.gz"
outfile="/scratch/project_2006203/TFBS/Results/MOODS_Teemu_2/MOODS_"$array".csv.gz"
#outfile=../Results/MOODS/MOODS_"$array".csv"

pwms=($(ls ../PWMs/*/pwms_space/*/*.pfm)) #3982
#echo "${#pwms[@]}"
#echo "${pwms[@]}"

#which pwms are in files
filenames=(${pwms[@]##*/})
filenames2=(${filenames[@]%%.pfm})

($(awk 'BEGIN{RS = FS} NR == FNR {files[$1] = 1; next} $1 in files' \
    <(echo "${files[*]}") <(echo "${filenames2[*]}")))

declare -A indices
declare -A pwms_new

i=0
for file in "${files[@]}"
do
   #echo $i
   #echo $file
   #echo ${filenames2[@]/$file//} | cut -d/ -f1 | wc -w | tr -d ' '
   tmp=$(echo ${filenames2[@]/$file//} | cut -d/ -f1 | wc -w | tr -d ' ')
   #echo $tmp
   indices[$i]=$(echo ${filenames2[@]/$file//} | cut -d/ -f1 | wc -w | tr -d ' ')
   pwms_new[$i]=${pwms[${indices[$i]}]}
   i=$((i+1))
done


#echo "${#indices[@]}"
#echo "${indices[2]}"

#echo "${#pwms_new[@]}"
#echo "${pwms_new[0]}"
#echo "${files[0]}" SAME


#str='test1@test2'
#echo "${str#*@}"
#The # character says Remove the smallest prefix of the expansion matching the pattern.
#The % character means Remove the smallest suffix of the expansion matching the pattern. (So you can do "${str%@*}" to get the "test1" part.)
#The / character means Remove the smallest and first substring of the expansion matching the following pattern. Bash has it, but it's not POSIX.
#If you double the pattern character it matches greedily.

## means Remove the largest prefix of the expansion matching the pattern.
#%% means Remove the largest suffix of the expansion matching the pattern.
#// means Remove all substrings of the expansion matching the pattern.    

    
#1024 -> 0-102
#762 -> 0-76
nro_pwms=${#pwms_new[@]} 


start_ind=$(($array*10)) #100
end_ind=$((($array+1)*10 -1)) #100
length=10 #100

if [[ $end_ind -gt $(($nro_pwms-1)) ]] #
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
#the background distribution ?? generates 
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

#Teemu: Yleens?? olen hakenut l??ys??ll?? affineettirajalla ja 
#valinnut sitten k??rjest?? halutun m????r??n. Suuria m????ri?? hakiessa 
#jokin minimiaffiniteettiraja (esim. 9) on kuitenkin ollut k??yt??ss??, 
#ettei todella pitkille motiiveille v??kisin etsit?? ihan j??rjett??mi?? osumia.


#raja on 9 on alunperin EEL:ist??, jonka affiniteetit ovat log2 kun taas moodsissa oletuskantaluku on e (--log-base). 

#Jokaiselle motiiville etsit????n 300 000 matchi???

#keskeinen ongelma t??ss?? on se, ett?? mielivaltaisella matriisilla ei 
#v??ltt??m??tt?? ole mit????n pisterajaa, jolla se antaa t??sm??lleen 300 000 matchia, 
#tai edes l??helle sit??. T??m?? on erityisesti ongelma matriiseilla j
#otka ovat melko tasaisesti jakautuneita joissain positioissa.

#MOODSin -B parametri muistaakseni yritt???? etsi?? jonkun rajan jolla tulee v??hint????n 
#k matchia ja enint????n 2*k matchia, ja jos t??m?? ei ole mahdollista, 
#se valitsee pienimm??n rajan jolla tulee yli k matchia. Valitettavasti t??m?? 
#pakolti tarkoittaa sit??, ett?? joissain tapauksissa saat kasan matcheja 
#joista isolla osalla on sama score.

#Vastaavasti --p-value P antaa noin P*m matchia, kun sekvenssin pituus on m, 
#mutta vain tilastollisessa mieless?? eli olettaen ett?? sekvenssi on oikeasti 
#jakautunut taustajakauman mukaan ja on pitk??. 
#K??yt??nn??ss?? matchien m????r?? on tietysti aika vaihteleva.

#Eli k??yt??nn??ss?? joillekin matriiseille voi l??yty?? esim. 1000 hyv???? hitti??, 
#mutta sitten seuraavaksi parhailla 4 miljoonalla hitill?? on kaikilla sama pistem????r??. 
#N??iss?? tapauksissa filtter??inti on sitten teht??v?? jollain muulla kriteerill??. 


#Olikos niin, ett?? jos ei aseta tuota taustajakaumaa (--lo-bg) tausta estimoidaan genomisekvenssist???

#Joo, juurikin n??in. Eli --bg asettaa p-arvojen laskennassa k??ytett??v??n taustajakauman, ja 
#--lo-bg asettaa log-odds-matriisien laskennassa k??ytett??v??n taustajakauman. 
#Jos jompaa kumpaa ei anneta, tausta kyseiseen tarkoitukseen estimoidaan aina jokaisen sekvenssin kohdalla erikseen kyseisest?? sekvenssist??.

#Yleens?? ajattelen, ett?? --lo-bg on loogista antaa parametrina, 
#koska motiivien jakaumaa verrataan johonkin fiksattuun prosessiin 
#joka generoi taustajakauman (en tosin ole biologi, joten en tied?? kannattaako minuun luottaa t??ss??.)


moods-dna.py -m ${pwms_new[@]:$start_ind:$length} --threshold 2 -s $S | gzip > $outfile #${array[@]:START:LENGTH}

#gzip $outfile

#outfile to parquet, where does the file go, requires a lot of space
#csvcli $outfile convert -to "parquet"


#FATAL:   container creation failed: hook function for tag layer returns error: failed to create /tmp/nvme/job_14238673 directory: mkdir /tmp/nvme/job_14238673: permission denied

#/projappl/project_2006203/softwares/conda_envs/MOODS/bin/moods-dna.py: line 22: 705506 Segmentation fault
#   /usr/bin/singularity --silent exec -B $DIR/../$SQFS_IMAGE:$INSTALLATION_PATH:image-src=/
# $DIR/../$CONTAINER_IMAGE bash -c "eval \"\$(/CSC_CONTAINER/miniconda/bin/conda shell.bash hook )\"
# && conda activate env1 &>/dev/null &&  exec -a $_O_SOURCE $DIR/moods-dna.py $(
#test $# -eq 0 || printf " %q" "$@" )"
