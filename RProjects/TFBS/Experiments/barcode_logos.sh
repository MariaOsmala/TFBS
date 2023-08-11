

cd /Users/osmalama/projects/TFBS/

pwms=($(ls PWMs_final/*/pwms/*/*.pfm)) 

for file in "${pwms[@]}"
do
     echo $file
    /Users/osmalama/softwares/codes_from_Nitta2015/spacek40 --logo -barcodelogo -heightscaledbars -colorscaledbars -path -noname $file

done

#move all png and svg files into barcode_logos

mkdir barcode_logos/png
mkdir barcode_logos/svg


pngs=($(ls PWMs_final/*/pwms/*/*.png)) 

svgs=($(ls PWMs_final/*/pwms/*/*.svg)) 


for file in "${pngs[@]}"
do 
  mv $file barcode_logos/png/
done

for file in "${svgs[@]}"
do 
  mv $file barcode_logos/svg/
done