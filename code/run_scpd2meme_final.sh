export PATH="/projappl/project_2006203/softwares/conda_envs/motif-clustering-Vierstra/bin:$PATH"
cd /projappl/project_2006203/TFBS

i=./PWMs_final/all.scpd #Version 1 motifs

i=./PWMs_final_union/all.scpd #Union of Version 1 and Version 2 motifs

i=./PWMs_final_versio2.2/all.scpd #Final version2.2

scpd2meme PWMs_final_version2.2/all.scpd -pseudo 1 > PWMs_final_version2.2/all.meme

scpd2meme $i -pseudo 1 > ${i/scpd/meme} #-pseudo <total pseudocounts> add <total pseudocounts> times letter background to each frequency; default: 0

grep "^MOTIF" ./PWMs_final/all.meme | awk '{print $2}' > code/motifs_in_meme.txt
