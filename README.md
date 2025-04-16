# Transcription factor motif and binding site analysis

This repository contains the code, computing environment description, and software/package versions for the analysis of transcription factor motifs and binding sites in the human genome.

The repository also contains the Taipale Lab collection of 3933 motifs of which 1232 are representative. The motifs are in different formats (.pfm tab and space separated, transfac, and cspd) in the following folder: See also [here](Data/SELEX-motif-collection/README.md)

```         
Data/PWMs/
```

Motif metadata e.g. with info for representativeness of the motifs (representative=="YES") is in

```         
Data/SELEX-motif-collection/metadata_final.tsv
```

The .pfm and transfac files for only the representatives are in

```         
Data/forSharing/representatives
```

The codes

-   create the Taipale Lab SELEX motif collection including heterodimeric motifs from the new CAP-SELEX experiments create motif in different formats and visualises logos

-   create artificial half-site motifs of homodimers

-   create scrambled control motifs (artificial control motifs)

-   contain commands to perform various motif similarity computations

-   perform minimum dominating set analysis

-   analyse the difference between composite motif center and overlapping flanks of monomeric motifs aligned against the composite motif

-   contain commands to perform motif matching

-   perform motif enrichment analysis at cell-type-specific candidate cis-regulatory elements (cCREs)

-   perform logistic regression analysis to predict cell-type-specificity of a cCRE based on motif matches

## Cmputing environment and package/tool installation

```{=html}
<!-- This is a comment and will not be displayed
[devtools](https://devtools.r-lib.org):
-->
```

The codes have been run on [CSC Puhti](https://docs.csc.fi/computing/systems-puhti/) or Apple MacBook M2 Pro Sequoai 15.3.

Used [R](https://www.r-project.org)-environments in Puhti are 4.2.1 and 4.3.0 [info](https://docs.csc.fi/apps/r-env/) [Container recipies](https://github.com/CSCfi/singularity-recipes/tree/main/r-env-singularity)

The documentation lists [Tykky](https://docs.csc.fi/computing/containers/tykky/) environments

The code depends on software tools, R-packages, and Bioconductor packages (listed in documentation)

To install the required software, R packages etc., on your laptop, workstation, cluster, etc. we recommend using [anaconda or miniconda](https://www.anaconda.com/products/individual) and the provided `renv-4.2.1.yml`, `renv-4.3.0.yml`, and `renv-4.4.2.yml` files:

Installation on Mac M2: The current conda channels (see conda_info.txt) do not support R-4.2.1 & R-4.3.0 with bioconductor packages so the installation is done with R-4.4.2.

Conda does not install all cran and bioconductor packages automatically so need to install some manually.

Install mamba to the conda base environment.

```         
conda info > conda_info.txt
conda config --set auto_update_conda false
conda config --set subdir osx-arm64
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority flexible
conda config --show

conda install -n base -c conda-forge mamba
conda activate base

mamba env create -f renv-4.4.2.yml -n TFBS
mamba activate TFBS


Rscript install_bioconductor_packages.4.4.2.R #Remember to change here the path to the mamba env library
```

## Running the analysis in the paper

The steps to reproduce the analysis as described in the paper are [here](Experiments/Steps.qmd)

The steps are also described here: `/RProjects/TFBS/code/motif_analysis_workflow.R`

## Citation:

Xie et al. (2025). DNA-guided transcription factor interactions extend human gene regulatory code

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Contact

Maria Osmala University of Helsinki\
Faculty of Medicine, Applied Tumor Genomics Research Program Unit Email: [firstname.surname\@helsinki.fi](mailto:firstname.surname@helsinki.fi){.email}\
[homepage](https://www.helsinki.fi/en/about-us/people/people-finder/maria-osmala-9460935)

## Acknowledgments

-   Thank you for Teemu Kivioja providing support for various steps.
