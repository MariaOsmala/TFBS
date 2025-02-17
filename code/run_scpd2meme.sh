export PATH="/projappl/project_2006203/softwares/conda_envs/motif-clustering-Vierstra/bin:$PATH"

cd /projappl/project_2006203/TFBS

mkdir -p Data/PWMs/meme

#This works
scpd2meme Data/PWMs/pfms_scpd/all.scpd -pseudo 1 > Data/PWMs/meme/all.meme

#Artificial HT-SELEX motifs
scpd2meme Data/PWMs/pfms_scpd/half-site-motifs.scpd -pseudo 1 > Data/PWMs/meme/half-site-motifs.meme

#Check that there are all motifs
#grep "^MOTIF" Data/PWMs/meme/all.meme | awk '{print $2}' > motifs_in_meme.txt
