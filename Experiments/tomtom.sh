#!/bin/bash
#SBATCH --job-name=tomtom
#SBATCH --output=/logs/%J.out
#SBATCH --error=/logs/%J.err
#SBATCH --account=project_2006203
#SBATCH --time=2-00:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --partition=small
##SBATCH --mail-type=BEGIN #uncomment to enable mail

export PATH="/projappl/project_2006203/softwares/conda_envs/motif-clustering-Vierstra/bin:/projappl/project_2006203/softwares/conda_envs/motif-clustering-Vierstra/bin:/projappl/project_2006203/softwares/conda_envs/motif-clustering-Vierstra/bin:/appl/spack/v018/install-tree/gcc-8.5.0/git-2.35.2-z47wj6/bin:/appl/opt/csc-cli-utils/bin:/appl/spack/v018/install-tree/gcc-11.3.0/openmpi-4.1.4-w2aekq/bin:/appl/spack/v018/install-tree/gcc-11.3.0/intel-oneapi-mkl-2022.1.0-p2wpb4/mkl/2022.1.0/bin/intel64:/appl/spack/v018/install-tree/gcc-8.5.0/gcc-11.3.0-i44hho/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/appl/bin:/users/mariaosm/.local/bin:/users/mariaosm/bin"

results_dir="/scratch/project_2006203//Results/tomtom_final"

tomtom -dist kullback -motif-pseudo 0.1 -text -min-overlap 1 ../PWMs_final/all.meme ../PWMs_final/all.meme > /scratch/project_2006203/Results/tomtom_final/tomtom/tomtom.all.txt

