#!/bin/bash

cd /Users/osmalama/spacek

PATH=/Users/osmalama/spacek:$PATH

bash

IFS=$'\n' read -d '' -r -a IDs < ~/projects/TFBS/Data/SELEX-motif-collection/motifnames.csv
IFS=$'\n' read -d '' -r -a files < ~/projects/TFBS/Data/SELEX-motif-collection/filenames.csv #/Data/pfms_tab/<motif>.pfm


#Copy motifs to data folder. The data folder needs to be in the spacek folder

mkdir -p /Users/osmalama/spacek/Data/PWMs/pfms_tab

cd /Users/osmalama/projects/TFBS/



for ((i = 0; i < 3933; i++)); do
    echo "Index: $i"
    
    file=${files[i]}
    cp $file /Users/osmalama/spacek/Data/PWMs/pfms_tab
done

#Copy reverse complements motifs to data folder. The data folder needs to be in the spacek folder

mkdir -p /Users/osmalama/spacek/Data/PWMs-reverse-complement/pfms_tab

cd /Users/osmalama/projects/TFBS/
#How to replace PWMs with PWMs-reverse-complement in files, use bash script

for ((i = 0; i < 3933; i++)); do
    echo "Index: $i"
    
    file=${files[i]}
    
    new_file="${file/PWMs/PWMS-reverse-complement}"
    echo "$new_file"
    
    #cp $file /Users/osmalama/spacek/Data/PWMs/pfms_tab
done


mkdir -p /Users/osmalama/spacek/Logos/PWMs/svg
mkdir -p /Users/osmalama/spacek/Logos/PWMs/png
mkdir -p /Users/osmalama/spacek/Logos/PWMs-reverse-complement/svg
mkdir -p /Users/osmalama/spacek/Logos/PWMs-reverse-complement/png
mkdir -p /Users/osmalama/spacek/Logos/PWMs/barcode_svg
mkdir -p /Users/osmalama/spacek/Logos/PWMs/barcode_png
mkdir -p /Users/osmalama/spacek/Logos/PWMs-reverse-complement/barcode_svg
mkdir -p /Users/osmalama/spacek/Logos/PWMs-reverse-complement/barcode_png


#spacek40 --logo -barcodelogo -heightscaledbars -colorscaledbars -path -noname -o=/Users/osmalama/spacek/Figures/test ${monomer1_array[0]} 

#spacek40 --logo -noname ${monomer1_array[0]}

#convert command is deprecated, use magick or magick convert

cd /Users/osmalama/spacek/Data/PWMS/pfms_tab
for ((i = 0; i < 3933; i++)); do
    echo $i
    spacek40 --logo -heightscaledbars -colorscaledbars -path -noname ${IDs[$i]}".pfm"
    
    mv ${IDs[$i]}".pfm.svg" /Users/osmalama/spacek/Logos/PWMs/svg/${IDs[i]}".svg"
    mv ${IDs[$i]}".pfm.png" /Users/osmalama/spacek/Logos/PWMs/png/${IDs[i]}".png"
    
    
    spacek40 --logo -barcodelogo -heightscaledbars -colorscaledbars -path -noname ${IDs[$i]}".pfm"
    mv ${IDs[$i]}".pfm.svg" /Users/osmalama/spacek/Logos/PWMs/barcode_svg/${IDs[i]}".svg"
    mv ${IDs[$i]}".pfm.png" /Users/osmalama/spacek/Logos/PWMs/barcode_png/${IDs[i]}".png"
     
done

#Remove empty space from .svg files? There is no empty space?
#cd /Users/osmalama/spacek/Logos/PWMs/svg/
#for file in *.svg; do
#    echo "Processing $file"
#    sed -i '' -e 's/width="100%" height="100%"/width="100%" height="100%"/g' $file
#done


#Draw also reverse complement logos, correct paths!




cd /Users/osmalama/spacek/Data_version2.2_rev_comp
for ((i = 0; i < 3932; i++)); do
     echo $i
    spacek40 --logo -heightscaledbars -colorscaledbars -path -noname ${IDs[$i]}".pfm"
    
    mv ${IDs[$i]}".pfm.svg" /Users/osmalama/spacek/Figures_version2.2_rev_comp_Logos/svg/${IDs[i]}".svg"
    mv ${IDs[$i]}".pfm.png" /Users/osmalama/spacek/Figures_version2.2_rev_comp_Logos/png/${IDs[i]}".png"
    
    
    spacek40 --logo -barcodelogo -heightscaledbars -colorscaledbars -path -noname ${IDs[$i]}".pfm"
    mv ${IDs[$i]}".pfm.svg" /Users/osmalama/spacek/Figures_version2.2_rev_comp_Logos/barcode_svg/${IDs[i]}".svg"
    mv ${IDs[$i]}".pfm.png" /Users/osmalama/spacek/Figures_version2.2_rev_comp_Logos/barcode_png/${IDs[i]}".png"
     
done



#Artificial half sites

IFS=$'\n' read -d '' -r -a IDs < ~/projects/TFBS/Data/SELEX-motif-collection/artificial-half-site-motifs/motifnames.csv
IFS=$'\n' read -d '' -r -a files < ~/projects/TFBS/Data/SELEX-motif-collection/artificial-half-site-motifs/filenames.csv #/Data/pfms_tab/<motif>.pfm


#Copy files to data folder. The data folder needs to be in the spacek folder

target_dir=/Users/osmalama/spacek/Data/artificial-half-site-motifs/PWMs/pfms_tab
mkdir -p$target_dir

cd /Users/osmalama/projects/TFBS/

for ((i = 0; i < 7; i++)); do
    #echo "Index: $i"
    
    file=${files[i]}
    echo $file
    cp $file $target_dir 
done

mkdir -p /Users/osmalama/spacek/Logos/artificial-half-site-motifs/PWMs/svg
mkdir -p /Users/osmalama/spacek/Logos/artificial-half-site-motifs/PWMs/png
mkdir -p /Users/osmalama/spacek/Logos/artificial-half-site-motifs/PWMs-reverse-complement/svg
mkdir -p /Users/osmalama/spacek/Logos/artificial-half-site-motifs/PWMs-reverse-complement/png
mkdir -p /Users/osmalama/spacek/Logos/artificial-half-site-motifs/PWMs/barcode_svg
mkdir -p /Users/osmalama/spacek/Logos/artificial-half-site-motifs/PWMs/barcode_png
mkdir -p /Users/osmalama/spacek/Logos/artificial-half-site-motifs/PWMs-reverse-complement/barcode_svg
mkdir -p /Users/osmalama/spacek/Logos/artificial-half-site-motifs/PWMs-reverse-complement/barcode_png


#spacek40 --logo -barcodelogo -heightscaledbars -colorscaledbars -path -noname -o=/Users/osmalama/spacek/Figures/test ${monomer1_array[0]} 

#spacek40 --logo -noname ${monomer1_array[0]}

#convert command is deprecated, use magick or magick convert

cd $target_dir
for ((i = 0; i < 7; i++)); do
    echo $i
    spacek40 --logo -heightscaledbars -colorscaledbars -path -noname ${IDs[$i]}".pfm"
    
    mv ${IDs[$i]}".pfm.svg" /Users/osmalama/spacek/Logos/artificial-half-site-motifs/PWMs/svg/${IDs[i]}".svg"
    mv ${IDs[$i]}".pfm.png" /Users/osmalama/spacek/Logos/artificial-half-site-motifs/PWMs/png/${IDs[i]}".png"
    
    
    spacek40 --logo -barcodelogo -heightscaledbars -colorscaledbars -path -noname ${IDs[$i]}".pfm"
    mv ${IDs[$i]}".pfm.svg" /Users/osmalama/spacek/Logos/artificial-half-site-motifs/PWMs/barcode_svg/${IDs[i]}".svg"
    mv ${IDs[$i]}".pfm.png" /Users/osmalama/spacek/Logos/artificial-half-site-motifs/PWMs/barcode_png/${IDs[i]}".png"
     
done


#Draw also reverse complement logos, correct paths!

IFS=$'\n' read -d '' -r -a IDs < /Users/osmalama/projects/TFBS/PWMs_final_version2.2/motifnames.csv


cd /Users/osmalama/spacek/Data_version2.2_rev_comp
for ((i = 0; i < 3932; i++)); do
     echo $i
    spacek40 --logo -heightscaledbars -colorscaledbars -path -noname ${IDs[$i]}".pfm"
    
    mv ${IDs[$i]}".pfm.svg" /Users/osmalama/spacek/Figures_version2.2_rev_comp_Logos/svg/${IDs[i]}".svg"
    mv ${IDs[$i]}".pfm.png" /Users/osmalama/spacek/Figures_version2.2_rev_comp_Logos/png/${IDs[i]}".png"
    
    
    spacek40 --logo -barcodelogo -heightscaledbars -colorscaledbars -path -noname ${IDs[$i]}".pfm"
    mv ${IDs[$i]}".pfm.svg" /Users/osmalama/spacek/Figures_version2.2_rev_comp_Logos/barcode_svg/${IDs[i]}".svg"
    mv ${IDs[$i]}".pfm.png" /Users/osmalama/spacek/Figures_version2.2_rev_comp_Logos/barcode_png/${IDs[i]}".png"
     
done


