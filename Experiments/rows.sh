#!/bin/bash

array=$1 #this varies between 0 and 39

#cd /scratch/project_2006472/MOODS
cd /scratch/project_2006203/TFBS/Results/MOODS

unpigz -c MOODS_"$array".csv.gz | wc -l > "nro_rows_"$array".txt"
