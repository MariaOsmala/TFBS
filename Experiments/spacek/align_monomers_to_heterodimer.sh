cd /Users/osmalama/spacek

PATH=/Users/osmalama/spacek:$PATH



#data_file=/Users/osmalama/projects/TFBS/PWMs_final_version2.2/metadata_for_tomtom_figures.tsv
data_file=/Users/osmalama/projects/TFBS/PWMs_final_version2.2/metadata_for_composite_spacek_figures.tsv
# Define empty arrays
dimer_ID=()
first_monomer_ID=()
second_monomer_ID=()

dimer_array=()
monomer1_array=()
monomer2_array=()

counter=0

mkdir /Users/osmalama/spacek/Data_version2.2_composites

# Read file
while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6 col7
do
    ((counter++))
    
    
    if [ $counter -eq 1 ]; then
        continue
    fi

    echo $col1
    dimer_ID+=("$col1")
    first_monomer_ID+=("$col2")
    second_monomer_ID+=("$col3")
    
    
    
    dimer_array+=("$col4")
    echo ${dimer_array[0]}
    monomer1_array+=("$col5")
    monomer2_array+=("$col6")
done < $data_file

#Copy files to data folder

cd /Users/osmalama/projects/TFBS/RProjects/TFBS

for ((i = 0; i < 1132; i++)); do #i< 1132(needs to be one larger than the number of composites(1131) as zero is empty)
    echo "Index: $i"

    dimer=${dimer_array[i]}
    monomer1=${monomer1_array[i]}
    monomer2=${monomer2_array[i]}

    #cp $dimer /Users/osmalama/spacek/Data_version2.2
    #cp $monomer1 /Users/osmalama/spacek/Data_version2.2
    #cp $monomer2 /Users/osmalama/spacek/Data_version2.2
    cp $dimer /Users/osmalama/spacek/Data_version2.2_composites
    cp $monomer1 /Users/osmalama/spacek/Data_version2.2_composites
    cp $monomer2 /Users/osmalama/spacek/Data_version2.2_composites
    
done



#spacek40 --logo -barcodelogo -heightscaledbars -colorscaledbars -path -noname -o=/Users/osmalama/spacek/Figures/test ${monomer1_array[0]} 

#spacek40 --logo -noname ${monomer1_array[0]}
mkdir /Users/osmalama/spacek/Figures_version2.2_composites/
mkdir /Users/osmalama/spacek/Alignment_info_version2.2_composites/


cd /Users/osmalama/spacek/Data_version2.2_composites
for ((i = 0; i < 1132; i++)); do #i< 1132(needs to be one larger than the number of composites as zero is empty)
    spacek40 --pwmalign ${first_monomer_ID[i]}".pfm" ${dimer_ID[i]}".pfm" ${second_monomer_ID[i]}".pfm" > ${first_monomer_ID[i]}"_"${dimer_ID[i]}"_"${second_monomer_ID[i]}".txt"
    mv ${first_monomer_ID[i]}".pfm_"${dimer_ID[i]}".pfm_"${second_monomer_ID[i]}".pfm.svg" /Users/osmalama/spacek/Figures_version2.2_composites/${dimer_ID[i]}".svg"
    mv ${first_monomer_ID[i]}"_"${dimer_ID[i]}"_"${second_monomer_ID[i]}".txt" /Users/osmalama/spacek/Alignment_info_version2.2_composites/
done

#Remove empty space from .svg files
cd /Users/osmalama/spacek/Figures_version2.2_composites/
for file in *.svg; do
    echo "Processing $file"
    sed -i '' -e 's/width="100%" height="100%"/width="100%" height="100%"/g' $file
done



