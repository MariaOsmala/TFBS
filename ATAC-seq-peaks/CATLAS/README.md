NOTE
====

1. Raw data can be downloaded at: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184462.

2. Additional data for generating the figures in the paper can be found here: https://dx.doi.org/10.17632/yv4fzv6cnm,
   and source codes can be found here: https://gitlab.com/kaizhang/hea2.

fragment
========

Fragment files.

Cell_by_cCRE
============

This directory stores raw TN5 insertion counts for each cCRE in each cell, containing three files:

1. matrix.mtx.gz: count matrix in matrix market format.
2. features.txt.gz: column names (cCREs) of the matrix.
3. barcodes.txt.gz: row names (cell barcodes) of the matrix.

Cell_by_gene
============

This directory stores raw TN5 insertion counts for each promoter region in each cell, containing three files:

1. matrix.mtx.gz: count matrix in matrix market format.
2. features.txt.gz: column names (genes) of the matrix.
3. barcodes.txt.gz: row names (cell barcodes) of the matrix.

cCRE_hg38.tsv.gz
================

The list of 1,154,611 cCREs. There are 7 columns in this file:

1. chromosome: name of the chromosome.
2. hg38_Start: 0-based starting location of the cCRE in the genome (hg38).
3. hg38_End: End of the cCRE in the genome (hg38).
4. class: Promoter (-200 to +200 of TSS), Promoter Proximal (less) or Distal
5. Present in fetal tissues: if this cCRE is detected in at least one fetal tissue
6. Present in adult tissues: if this cCRE is detected in at least one adult tissue
7. CRE module: The ID of CRE module that the cCRE belongs to.

Cell_metadata.tsv.gz
====================

This file contains metadata for each cell/barcode. It has 6 columns:

1. cellID: unique ID.
2. logUMI: log10 number of unique fragment passing QC.
3. tsse: TSS enrichment
4. tissue: The ID of tissue sample
5. cell type: The name of the cell type
6. Life stage: Fetal or Adult

UMAP_embedding
==============

This directory contains UMAP embeddings of cells:

1. adult.tsv.gz: UMAP embedding of cells from adult tissues
2. fetal+adult.tsv.gz: UMAP embedding of cells from fetal and adult tissues

Cell_ontology.tsv
=================

Matching between annotated cell type names and cell ontology terms.

cCRE_by_cell_type
=================

This directory stores cCRE by cell type, containing three files:

1. matrix.mtx.gz: matrix in matrix market format, each line represents a cCRE-cell type pair. For example, "2 22" means cCRE No. 2 is accessibile in cell type No. 22.
2. celltypes.txt.gz: column names (cell types) of the matrix.
3. cCREs.bed.gz: row names (cCREs) of the matrix.

(NOTE: the original files were wrong and the correct files were uploaded on Nov 26, 2021.)

Bigwig
======

This directory contains bigwig files for 222 fetal and adult cell types

Peaks
======

This directory contains peaks called for each of the 222 cell types, in NarrowPeak format: https://genome.ucsc.edu/FAQ/FAQformat.html#format12.

We first called peaks in each individual cell type using different FDR cutoffs
adapted to read depths (see paper for details). These peaks can be found in the
`RAW` subdirectory.

We then merge all peaks to obtain fix-width (400bp) non-overlapping peak list.
For each cell type, we intersected this peak list with originally called peaks
and further filter these peaks based on the fraction of accessible cells (see paper for details).

ABC_scores
==========

This directory contains cCRE to gene linkage predicted by the ABC model.


yv4fzv6cnm-4 contains data used to produce the figures, cCRE modules
