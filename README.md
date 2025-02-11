# Transcription factor motif and binding site analysis

This repository contains the code, computing environment description, and 
software/package versions for the analysis of transcription factor motifs and binding sites in the human genome.

The codes 
  - create the Taipale Lab SELEX motif collection including heterodimeric motifs from the new CAP-SELEX experiments
  - contain commands to perform various motif similarity computations 
  - perform minimum dominating set analysis
  - contain commands to perform motif matching
  - perform motif enrichment analysis at cell-type-specific candidate cis-regulatory elements (cCREs)
  - perform logistic regression analysis to predict cell-type-specificity of a cCRE based on motif matches
  - analyse the difference between composite motif center and overlapping flanks of monomeric motifs aligned againt the composite motif

**More documentation and instructions on how to run the code are coming soon.**

## Installation

<!-- This is a comment and will not be displayed
[devtools](https://devtools.r-lib.org):
-->


You can install the package using 

```R
```

The code depends on software tools, R-packages, and Bioconductor packages which are listed in the following files.

  - [R](https://www.r-project.org) (version 3)

To install the required software, R packages etc., 
we recommend using [anaconda](https://www.anaconda.com/products/individual) and the provided `conda_environment.yml` file:

```
conda env create -f conda_environment.yml -n TFBS
source activate TFBS
```


## Running the analysis in the paper

The scripts to reproduce the analysis pipeline as described in the paper, are in the `name/name` folder.

to run the analysis by
first editing `name/config.yaml` to specify the folders in which to place the data.


## Citation:

Xie et al. (2025). DNA-guided transcription factor interactions extend human gene regulatory code

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Contact

Maria Osmala
University of Helsinki  
Faculty of Medicine, Applied Tumor Genomics Research Program Unit
Email: firstname.surname@helsinki.fi  
[homepage](https://www.helsinki.fi/en/about-us/people/people-finder/maria-osmala-9460935)

## Acknowledgments

* Thank you for Teemu Kivioja providing support for various steps.
