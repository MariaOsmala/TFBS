#!/bin/bash

array=$1 #this varies between 0 and 39

cd /scratch/project_2006203/TFBS/Results/MOODS

gzip MOODS_"$array".csv
