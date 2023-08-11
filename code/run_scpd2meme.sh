i=./PWMs_final/all.scpd

scpd2meme $i -pseudo 1 > ${i/scpd/meme}

grep "^MOTIF" ./PWMs_final/all.meme | awk '{print $2}' > code/motifs_in_meme.txt

find ./PWMs_final -name "*.scpd*"



paths=($(find ./PWMs_final -name "*.scpd"))

printf "%s\n" "${paths[@]}"

for i in "${paths[@]}"
do
   :
   echo ${i/scpd/meme}
   scpd2meme $i -pseudo 1 > ${i/scpd/meme}
   # do whatever on "$i" here
done

Converted 1787 motifs.
Skipped 0 motifs.
0 conversion errors.
./PWMs_final/Morgunova2015/all.meme
Converted 1 motifs.
Skipped 0 motifs.
0 conversion errors.
./PWMs_final/Nitta2015/all.meme
Converted 10 motifs.
Skipped 0 motifs.
0 conversion errors.
./PWMs_final/Jolma2013/all.meme
Converted 820 motifs.
Skipped 0 motifs.
0 conversion errors.
./PWMs_final/fromYimeng/Homo_sapiens_all.meme
Converted 709 motifs.
Skipped 0 motifs.
0 conversion errors.
./PWMs_final/Jolma2015/all.meme
Error: Motif position 11 summed to zero.
Error: Motif position 11 summed to zero.
Error: Motif position 11 summed to zero.
Error: Motif position 11 summed to zero.
Error: Motif position 9 summed to zero.
Error: Motif position 11 summed to zero.
Error: Motif position 10 summed to zero.
Error: Motif position 11 summed to zero.
Error: Motif position 11 summed to zero.
Error: Motif position 11 summed to zero.
Converted 583 motifs.
Skipped 0 motifs.
10 conversion errors.

wc -l ./PWMs_final/Yin2017/metadata.csv #1161
wc -l ./PWMs_final/Morgunova2015/metadata.csv #1
wc -l ./PWMs_final/Nitta2015/metadata.csv #10
wc -l ./PWMs_final/Jolma2013/metadata.csv #820
wc -l ./PWMs_final/fromYimeng/metadata.csv #719
wc -l./PWMs_final/Jolma2015/metadata.csv #593



wc -l ./PWMs/Jolma2013/metadata.csv #844
wc -l ./PWMs/Jolma2015/metadata.csv #663
wc -l ./PWMs/fromYimeng/metadata.csv #706
wc -l ./PWMs/Yin2017/metadata.csv #1795
wc -l ./PWMs/Nitta2015/metadata.csv #11


./PWMs/Jolma2013/Mus_musculus_all.meme
Converted 134 motifs.
Skipped 0 motifs.
0 conversion errors.

./PWMs/Jolma2013/Homo_sapiens_all.meme
Converted 709 motifs.
Skipped 0 motifs.
0 conversion errors.

./PWMs/Jolma2015/Homo_sapiens_all.meme
Error: Motif position 11 summed to zero.
Error: Motif position 11 summed to zero.
Error: Motif position 11 summed to zero.
Error: Motif position 11 summed to zero.
Error: Motif position 9 summed to zero.
Error: Motif position 11 summed to zero.
Error: Motif position 10 summed to zero.
Error: Motif position 11 summed to zero.
Error: Motif position 11 summed to zero.
Error: Motif position 11 summed to zero.
Converted 652 motifs.
Skipped 0 motifs.
10 conversion errors.

./PWMs/fromYimeng/Homo_sapiens_all.meme
Error: Expected count matrix but got probability matrix. Can't continue without site count.
Error: Expected count matrix but got probability matrix. Can't continue without site count.
Error: Expected count matrix but got probability matrix. Can't continue without site count.
Error: Expected count matrix but got probability matrix. Can't continue without site count.
Error: Expected count matrix but got probability matrix. Can't continue without site count.
Error: Expected count matrix but got probability matrix. Can't continue without site count.
Error: Expected count matrix but got probability matrix. Can't continue without site count.
Error: Expected count matrix but got probability matrix. Can't continue without site count.
Error: Expected count matrix but got probability matrix. Can't continue without site count.
Error: Expected count matrix but got probability matrix. Can't continue without site count.
Converted 695 motifs.
Skipped 0 motifs.
10 conversion errors.

./PWMs/Yin2017/Homo_sapiens_all.meme
Error: Expected count matrix but got probability matrix. Can't continue without site count.
Error: Empty matrix.
Error: Empty matrix.
Error: Empty matrix.
Error: Empty matrix.
Error: Expected count matrix but got probability matrix. Can't continue without site count.
Error: Expected count matrix but got probability matrix. Can't continue without site count.
Converted 1787 motifs.
Skipped 0 motifs.
7 conversion errors.

./PWMs/Nitta2015/Homo_sapiens_all.meme
Converted 10 motifs.
Skipped 0 motifs.
0 conversion errors.
