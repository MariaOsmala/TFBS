#!/bin/bash

# Die on unset variables
set -o nounset
# Die on errors
set -o errexit
# Die if any part of a pipe fails
set -o pipefail

# An example of creating PWM dominating set for our published SELEX matrices
# T. Kivioja, January 4, 2017

results_path=Results/SELEX-Dominating-Set-Analysis

## 1. Create ILP dat file from SSTAT output

#sstat_to_DSA.tsv is obtaoined using code/motif-similarity/extract_SSTAT.R

perl scripts/sstat2ilpdat.pl \
<../TFBS/PWMs_final_version2.2/sstat_to_DSA.tsv \
> $results_path""data/Selex_all_version2.2.dat

## 1. Create ILP dat file from SSTAT output, only composites

perl scripts/sstat2ilpdat.pl \
<../TFBS/PWMs_final_version2.2/sstat_to_DSA_new_composites.tsv \
> $results_path""data/Selex_new_composites_version2.2.dat



# 3. The actual ILP solver call
#GLPSOL: GLPK LP/MIP Solver, v4.65

#-m filename, --model filename
#read model section and optional data section from
#                    filename (same as --math)
#-d filename, --data filename
#read data section from filename (for --math only);
#if model file also has data section, it is ignored
#-o filename, --output filename
#write solution to filename in printable format

glpsol -m models/min_dom_set.mod -d $results_path""data/Selex_all_version2.2.dat -o $results_path""solutions/SELEX_all_version2.2.sol #fromYimeng_version2.2

glpsol -m models/min_dom_set.mod -d $results_path""data/Selex_new_composites_version2.2.dat -o $results_path""solutions/SELEX_new_composites_version2.2.sol #fromYimeng_version2.2 only composites


perl scripts/parse_dom_set.pl <$results_path""solutions/SELEX_all_version2.2.sol \
> $results_path""solutions/SELEX_all_version2.2_min_dom_set_list.txt

#new composites
perl scripts/parse_dom_set.pl <$results_path""solutions/SELEX_new_composites_version2.2.sol \
>$results_path""solutions/SELEX_new_composites_version2.2_min_dom_set_list.txt



