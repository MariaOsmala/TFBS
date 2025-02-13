# Steps to reproduce the TFBS analysis

First edit `Experiments/config.yaml` to specify the folders in which to place the data.

These are described in `/code/motif_analysis_workflow.R`

cp PWMs/Jolma2013/mmc3.xls Data/
cp PWMs/Jolma2015/41586_2015_BFnature15518_MOESM33_ESM.xlsx Data/
cp PWMs/Nitta2015/elife-04837-supp1-v1.xlsx Data/
cp PWMs/Morgunova2015/E2F8_cycle4.adm Data/
cp ../motif-clustering-Viestra-private/HumanTFs-ccbr-TableS1.csv Data/
cp PWMs/Yin2017/NIHMS1670407-supplement-Supplemental_tables_S1-6.xlsx Data/
families=read_csv("../../../motif-clustering-Viestra-private/HumanTFs-ccbr-TableS1.csv")


code/SELEX-motif-collection/read_excel_tables.R
code/SELEX-motif-collection/extract_motifs_from_excel_Jolma2013.R 
code/SELEX-motif-collection/extract_motifs_from_excel_Jolma2015.R 
code/SELEX-motif-collection/extract_motifs_from_excel_Nitta2015.R  
code/SELEX-motif-collection/extract_motifs_Morgunova_2015.R 
code/SELEX-motif-collection/extract_motifs_from_excel_Yin2017.R 
code/SELEX-motif-collection/extract_motifs_from_excel_fromYimeng_version2.2.R 

code/SELEX-motif-collection/collect_all_motif_metadata.R 


