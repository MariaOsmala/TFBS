
source activate /home/local/osmalama/projects/motif-clustering-Viestra-private/conda_envs/motif-clustering-Vierstra
cd /home/osmalama/projects/TFBS/Experiments

Vierstra_path=/home/local/osmalama/projects/motif-clustering-Viestra-private

pfm_dir="../Results/tomtom/tomtom/pfms/"
mkdir $pfm_dir

meme_file=../Results/tomtom/all.dbs.meme
python $Vierstra_path/bin/meme2jaspar.py ${meme_file} ${pfm_dir}
