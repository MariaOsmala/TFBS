
ln -s /Applications/Inkscape.app/Contents/MacOS/inkscape \
      /usr/local/bin/inkscape


#Trim the Whitespace and Convert to PDF/PNG
#1. Trim and convert SVG to PNG:

#You can use Inkscape to export the SVG to PNG 
#while trimming the extra whitespace. 
#Inkscape's --export-area-drawing option helps crop the content to the drawing area (no extra whitespace).

mkdir -p /Users/osmalama/spacek/Logos/PWMs/pdf

cd /Users/osmalama/spacek/Logos/PWMs/svg



for file in *.svg; do
  basename="${file%.*}"
  echo $basename
  #inkscape $file --export-area-drawing --export-filename=~/spacek/Figures_version2.2_composites_png/$basename.png
  inkscape $file --export-area-drawing --export-filename=/Users/osmalama/spacek/Logos/PWMs/pdf/$basename.pdf
done

