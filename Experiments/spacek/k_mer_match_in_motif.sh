cd /Users/osmalama/spacek

PATH=/Users/osmalama/spacek:$PATH



IFS=$'\n' read -d '' -r -a IDs < /Users/osmalama/projects/TFBS/PWMs_final_version2.2/motifnames.csv


cd /Users/osmalama/spacek/Data_version2.2
k_mer=CACGTG
spacek40 --pwm - E_box_kmer.txt MAX_HT-SELEX_TGACCT20NGA_Y_NNCACGTGNN_1_2.pfm  


"MAX_HT-SELEX_TGACCT20NGA_Y_NNCACGTGNN_1_2.pfm""


cd /Users/osmalama/spacek/Data_version2.2
for ((i = 0; i < 3932; i++)); do
     echo $i
    spacek40 --logo -heightscaledbars -colorscaledbars -path -noname ${IDs[$i]}".pfm"
    
    mv ${IDs[$i]}".pfm.svg" /Users/osmalama/spacek/Figures_version2.2_Logos/svg/${IDs[i]}".svg"
    mv ${IDs[$i]}".pfm.png" /Users/osmalama/spacek/Figures_version2.2_Logos/png/${IDs[i]}".png"
    
    
    spacek40 --logo -barcodelogo -heightscaledbars -colorscaledbars -path -noname ${IDs[$i]}".pfm"
    mv ${IDs[$i]}".pfm.svg" /Users/osmalama/spacek/Figures_version2.2_Logos/barcode_svg/${IDs[i]}".svg"
    mv ${IDs[$i]}".pfm.png" /Users/osmalama/spacek/Figures_version2.2_Logos/barcode_png/${IDs[i]}".png"
     
done

#Remove empty space from .svg files
cd /Users/osmalama/spacek/Figures_version2.2_Logos/svg/
for file in *.svg; do
    echo "Processing $file"
    sed -i '' -e 's/width="100%" height="100%"/width="100%" height="100%"/g' $file
done



