#!/bin/bash

#Before this you need to run code/run_scpd2meme.sh

export PATH="/projappl/project_2006203/softwares/conda_envs/motif-clustering-Vierstra/bin:$PATH"


results_dir="/scratch/project_2006203/Results/tomtom_final"


#meme2meme ../PWMs/*/*.meme > ${base_dir}/all.dbs.meme

#meme2meme \
#	databases/jaspar2022/*.meme \
#	databases/cisbp_v2.0/*.meme \
#	databases/Grand2021/*.meme \
#> ${base_dir}/all.dbs.meme

#rm -rf ${results_dir}/logs && mkdir -p ${results_dir}/logs
rm -rf ${results_dir}/tomtom && mkdir -p ${results_dir}/tomtom

cat <<__SCRIPT__ > tomtom.sh
#!/bin/bash
#SBATCH --job-name=tomtom
#SBATCH --output=/outs/%J.out
#SBATCH --error=/errs/%J.err
#SBATCH --account=project_2006203
#SBATCH --time=2-00:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --partition=small
#SBATCH --gres=nvme:20 
##SBATCH --mail-type=BEGIN #uncomment to enable mail

##SBATCH --gres=nvme:50 This is in GB

export PATH="/projappl/project_2006203/softwares/conda_envs/motif-clustering-Vierstra/bin:$PATH"

results_dir="/scratch/project_2006203/Results/tomtom_final"

tomtom \
	-dist kullback \
	-motif-pseudo 0.1 \
	-text \
	-min-overlap 1 \
	../PWMs_final/all.meme ../PWMs_final/all.meme \
> ${results_dir}/tomtom/tomtom.all.txt

__SCRIPT__

JOB0=$(sbatch --export=ALL --parsable \
	--job-name=tomtom.chunk \
	tomtom.sh
echo $JOB0

#-dist allr|ed|kullback|pearson|sandelin|blic1|blic5|llr1|llr5
#                   Distance metric for scoring alignments;
#                    default: ed
