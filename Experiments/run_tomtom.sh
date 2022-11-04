#!/bin/bash

export PATH="/projappl/project_2006203/softwares/conda_envs/motif-clustering-Vierstra/bin:$PATH"


base_dir="../Results/tomtom"


meme2meme ../PWMs/*/*.meme > ${base_dir}/all.dbs.meme

#meme2meme \
#	databases/jaspar2022/*.meme \
#	databases/cisbp_v2.0/*.meme \
#	databases/Grand2021/*.meme \
#> ${base_dir}/all.dbs.meme

rm -rf ${base_dir}/logs && mkdir -p ${base_dir}/logs
rm -rf ${base_dir}/tomtom && mkdir -p ${base_dir}/tomtom

cat <<__SCRIPT__ > ${base_dir}/slurm.tomtom
#!/bin/bash
#SBATCH --job-name=motif-clustering
#SBATCH --output=${base_dir}/logs/%J.out
#SBATCH --error=${base_dir}/logs/%J.err
#SBATCH --account=project_2006203
#SBATCH --time=2-00:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --partition=small
##SBATCH --mail-type=BEGIN #uncomment to enable mail

export PATH="/projappl/project_2006203/softwares/conda_envs/motif-clustering-Vierstra/bin:$PATH"

tomtom \
	-dist kullback \
	-motif-pseudo 0.1 \
	-text \
	-min-overlap 1 \
	${base_dir}/all.dbs.meme ${base_dir}/all.dbs.meme \
> ${base_dir}/tomtom/tomtom.all.txt

__SCRIPT__

JOB0=$(sbatch --export=ALL --parsable \
	--job-name=tomtom.chunk \
	${base_dir}/slurm.tomtom)
echo $JOB0
