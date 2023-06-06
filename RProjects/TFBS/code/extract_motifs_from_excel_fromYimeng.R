

#############Read excel tables################

library(tidyverse)
library(universalmotif)

library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(readr)



rm(list=ls())

types=c("pfm_composite_new", "pfm_spacing_new", "20230420")

path="../../PWMs/fromYimeng/pfm_from_newData/"

dir("../../PWMs/fromYimeng/pfm_from_newData/pfm_composite_new/")

datas=tibble()
  
for(t in types){
  data=as_tibble(do.call(rbind,strsplit(dir(paste0(path, t, "/")),"_")))
  data$V8=gsub(".pfm", "",as.character(data$V8))
  data$symbol=paste(data$V1, data$V2,sep="_")
  data=data[,-c(1,2)]
  names(data)=c("ligand", "batch", "seed", "multinomial", "cycle", "short", "symbol")
  data=data %>% relocate(symbol, .before=ligand)
  data$study="fromYimeng"
  data$type=t
  data$organism="Homo_sapiens"
  data$clone=NA
  data$family=NA
  data$experiment="CAP-SELEX"
  data$representative=NA
  data$comment=NA
  data$filename=paste0("PWMs/fromYimeng/pfm_from_newData/", t,"/", dir(paste0(path, t, "/")))
  datas=rbind(datas, data)
}

datas=datas[, c("symbol",	"clone","family",	"organism","study","experiment",
                "ligand",	"batch", "seed",	"multinomial",
                "cycle","representative", "short","type","comment","filename")]

#datas$multinomial=gsub("m", "",as.character(datas$multinomial))
#datas$cycle=gsub("c", "",as.character(datas$cycle))

#write.table(datas, file="../../PWMs/fromYimeng/metadata.csv", row.names = FALSE,sep="\t")
#saveRDS(datas, file="data/fromYimeng.Rds")

dir.create(paste0("../../PWMs_final/fromYimeng/pwms_space/pfm_composite_new/"), recursive=TRUE)
dir.create(paste0("../../PWMs_final/fromYimeng/pwms_space/pfm_spacing_new/"), recursive=TRUE)
dir.create(paste0("../../PWMs/fromYimeng/transfac/pfm_composite_new/"), recursive=TRUE)
dir.create(paste0("../../PWMs/fromYimeng/transfac/pfm_spacing_new/"), recursive=TRUE)

pwms=paste0("../../PWMs_final/fromYimeng/pwms/Homo_sapiens/")
pwms_space=paste0("../../PWMs_final/fromYimeng/pwms_space/Homo_sapiens/")
transfac_path=paste0("../../PWMs_final/fromYimeng/transfac/Homo_sapiens/")
dir.create(pwms, recursive=TRUE)
dir.create(pwms_space, recursive=TRUE)
dir.create(transfac_path, recursive=TRUE)

#write pwms as .cspd
#write pwms as space separated

append=FALSE

datas$filename_final=""
datas$IC=NA
datas$length=NA
datas$IC_universal=NA
datas$consensus=""
#Remove these from metadata, empty TFs
remove=c()


# cspd2meme_failed=["ELK1_HOXA3_CAP-SELEX_TCTGTG40NAGA_AAA_RSCGGTAATKR_1_2_YES", 
#                   "ETV5_HOXA2_CAP-SELEX_TCGGCG40NACG_AY_RSCGGWAATKR_1_3_NO", 
#                   "FOXJ3_TBX21_CAP-SELEX_TCCTCC40NCTC_AAE_WAACAACACMY_1_3_YES",
#                   "HOXB2_ETV1_CAP-SELEX_TGTATA40NGGT_AY_ACCGGAAATGA_2_2_NO", 
#                   "MEIS1_DRGX_CAP-SELEX_TAGGTG40NTGT_AT_TGTCAATTA_1_3_YES", 
#                   "PBX4_HOXA1_CAP-SELEX_TTCTGT40NAAC_AS_NYNATMAATCA_1_2_YES", 
#                   "PBX4_HOXA10_CAP-SELEX_TGACAA40NACT_AS_CCATAAATCA_1_2_NO", 
#                   "PBX4_HOXA10_CAP-SELEX_TGACAA40NACT_AS_NTCGTAAATCA_1_2_YES", 
#                   "TEAD4_EOMES_CAP-SELEX_TGTCGT40NGCG_AY_RGAATGYGTGA_1_3_YES", 
#                   "ETV5_EVX1_CAP-SELEX_TCGTGA40NGTC_AY_RSCGGWAATKR_1_3_NO", 
#                   "ARNT2_Methyl-HT-SELEX_TGTTAG40NACC_KT_AA_1_2_NA", "HES6_Methyl-HT-SELEX_TAACCG40NCAA_KAE_AA_1_4_NA", 
#                   "MYCN_Methyl-HT-SELEX_40NTATAATTAATT_KAL_AA_1_3b0_NA", "ATF4_Methyl-HT-SELEX_TACCGC40NTAC_KW_AA_1_3b0_NA", 
#                   "E2F1_Methyl-HT-SELEX_TGCTTC40NTTA_KW_AA_2_4_NA", "E2F4_Methyl-HT-SELEX_TGGCCG40NGCT_KO_AA_1_4_NA", 
#                   "ZBTB33_HT-SELEX_TAGCTG40NCTA_KR_AA_1_4_NA", "ESX1_RUNX2_TTCGTC40NTGCG_YII_AA_m1_c3b0_short_pfm_spacing_new", 
#                   "ETV5_RUNX2_TCATGA40NCCTG_YII_NCCGCANNCCGGAWGYN_m1_c3b0u_short_pfm_spacing_new", "FLI1_ATF7_TCCTTC40NCGAT_YJI_AA_m1_c3b0_short_pfm_spacing_new",
#                   "HMX3_RUNX2_TACTGG40NGGA_YII_NRACCGCAATTAMN_m1_c3b0u_short_pfm_spacing_new", "HSF4_ONECUT1_TCGCAA40NTGT_YJII_AA_m1_c3b0_short_pfm_spacing_new", 
#                   "PITX1_FOXD2_TCGCGA40NGTT_YWI_AA_m1_c3b0_short_pfm_spacing_new", "POU3F4_ONECUT1_TCCCAA40NGAA_YJII_AA_m1_c3b0_short_pfm_spacing_new",
#                   "POU4F1_ONECUT1_TAACGC40NGTG_YJII_AA_m1_c3b0_short_pfm_spacing_new","POU4F3_ONECUT1_TAAGTT40NGAT_YJII_AA_m1_c3b0_short_pfm_spacing_new",
#                   "POU4F3_TEAD4_TACCAT40NGCG_YJII_AA_m1_c3b0_short_pfm_spacing_new"]

#These are empty
# [1] "PWMs/fromYimeng/pfm_from_newData/pfm_spacing_new/ESX1_RUNX2_TTCGTC40NTGCG_YII_AA_m1_c3b0_short.pfm"
# [1] "PWMs/fromYimeng/pfm_from_newData/pfm_spacing_new/ETV5_RUNX2_TCATGA40NCCTG_YII_NCCGCANNCCGGAWGYN_m1_c3b0u_short.pfm"
# [1] "PWMs/fromYimeng/pfm_from_newData/pfm_spacing_new/FLI1_ATF7_TCCTTC40NCGAT_YJI_AA_m1_c3b0_short.pfm"
# [1] "PWMs/fromYimeng/pfm_from_newData/pfm_spacing_new/HMX3_RUNX2_TACTGG40NGGA_YII_NRACCGCAATTAMN_m1_c3b0u_short.pfm"
# [1] "PWMs/fromYimeng/pfm_from_newData/pfm_spacing_new/HSF4_ONECUT1_TCGCAA40NTGT_YJII_AA_m1_c3b0_short.pfm"
# [1] "PWMs/fromYimeng/pfm_from_newData/pfm_spacing_new/PITX1_FOXD2_TCGCGA40NGTT_YWI_AA_m1_c3b0_short.pfm"
# [1] "PWMs/fromYimeng/pfm_from_newData/pfm_spacing_new/POU3F4_ONECUT1_TCCCAA40NGAA_YJII_AA_m1_c3b0_short.pfm"
# [1] "PWMs/fromYimeng/pfm_from_newData/pfm_spacing_new/POU4F1_ONECUT1_TAACGC40NGTG_YJII_AA_m1_c3b0_short.pfm"
# [1] "PWMs/fromYimeng/pfm_from_newData/pfm_spacing_new/POU4F3_ONECUT1_TAAGTT40NGAT_YJII_AA_m1_c3b0_short.pfm"
# [1] "PWMs/fromYimeng/pfm_from_newData/pfm_spacing_new/POU4F3_TEAD4_TACCAT40NGCG_YJII_AA_m1_c3b0_short.pfm"




for(m in 1:nrow(datas)){
  f=datas$filename[m]
  PWM=read.table(paste0("../../", f))
  #Empty matrix
  if( ncol(PWM) ==1 ) {
    remove=c(remove, m)
    print(f)
    #Zero matrix
  }else if( sum(as.numeric(as.matrix(PWM[,-1]))) ==0 ){
    remove=c(remove, m)    
    print(f)
  }else{
  
    filename=paste0(strsplit(f, "/")[[1]][-c(1,2,3,4)],collapse="/")
    #add extension
    filename_final=paste0("_",datas$type[m] ,".pfm",filename)
    
    
    if(paste0(pwms,filename_final) %in% datas$filename_final){
      print("Warning! Non-unique filenames and IDs")
    }
    
  
    ID=gsub(".pfm", "",filename_final)
    datas$ID[m]=ID
    
    #write if to pwms, tab separated
    write.table(PWM,row.names = FALSE, col.names=FALSE, quote=FALSE,
                file=paste0(pwms,filename_final) ,sep="\t")
    
    
    
    #paste0( strsplit( strsplit(f, "/")[[1]][5], ".pfm")[[1]], "_",datas$type[m])
    write.table(PWM,row.names = FALSE, col.names=FALSE, quote=FALSE,
                file=paste0(pwms_space,filename_final) ,sep=" ")
    
    
    pcm=as.matrix( read.csv(paste0(pwms,filename_final), header=FALSE, sep="\t") )
    colnames(pcm)=NULL
    rownames(pcm)=c("A", "C", "G", "T")
    
    pfm <- TFBSTools::PFMatrix( bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                                profileMatrix=pcm
    )
    
    #pwm <- TFBSTools::toPWM(pfm, type="log2probratio", pseudocounts=0.01)
    icm <- TFBSTools::toICM(pfm, pseudocounts=0.01, schneider=FALSE)
   
    datas$IC[m]=sum(rowSums(icm)) #6.77892 # universalmotif computes this wrong?
    
    #motif length
    datas$length[m]=length(pfm)
    
    motif=universalmotif::read_matrix(file=paste0(pwms,filename_final), sep="\t", header=FALSE)
    
   
    datas[m, "IC_universal"]=motif@icscore
    datas[m, "consensus"]=motif@consensus
    
    
    #Write transfac format
    
    transfac=paste0(transfac_path,filename_final)
    
    write_transfac(motif, file=transfac, overwrite = TRUE, append = FALSE)
    
    #file=paste0("../../PWMs/fromYimeng/pwms_space/",paste0(strsplit(f, "/")[[1]][-c(1,2,3)],collapse="/"))
    
    motif=universalmotif::read_matrix(file=paste0(pwms_space,filename_final), sep=" ", header=FALSE)
    motif@name=paste0( strsplit( strsplit(f, "/")[[1]][5], ".pfm")[[1]], "_",datas$type[m])
  
    transfac_file=paste0(transfac_path,filename_final)
    write_transfac(motif, file=transfac, overwrite = TRUE, append = FALSE)
    
    #Write .scpd format, all into a single file
    PWM=as.matrix(PWM, dimnames=NULL)
    rownames(PWM)=c("A", "C", "G", "T")
    write.table(paste0(">",  ID),   
                append=append, row.names = FALSE, col.names=FALSE, quote=FALSE,
                file=paste0("../../PWMs_final/fromYimeng/Homo_sapiens_all", ".scpd"))
    append=TRUE
    
    write.table(PWM,append=append, row.names = TRUE, col.names=FALSE, quote=FALSE,
                file=paste0("../../PWMs_final/fromYimeng/Homo_sapiens_all", ".scpd"))
    
    
    datas$filename_final[m]=paste0(pwms,filename_final)
    
  
  }
}

datas=datas[-remove,]

#Add family mappings

TF_family_mappings <- read_delim("~/projects/motif-clustering-Viestra-private/HumanTFs-ccbr-TableS1.csv", delim=",")

#commondf <- merge(datas, TF_family_mappings, by = 'symbol', all.x = TRUE)

#There is no singles in this data
#single <- datas$symbol[str_detect(datas$symbol, "_") == FALSE]
#single_index <- which(datas$symbol %in% single)
#Family_names <- rep("Unknown", nrow(datas))
#Family_names[single_index] <- commondf$DBD[single_index]

#pair <- datas$symbol[str_detect(datas$symbol, "_") == TRUE]
#pair_index <- which(datas$symbol %in% pair)
#[pair_index]
pairs <- as.data.frame(do.call(rbind, strsplit(datas$symbol, split="_")))
colnames(pairs)=c("symbol", "symbol1")
pairs$order=seq(nrow(pairs))

commondf_pairs <- merge(pairs, TF_family_mappings, by = 'symbol', all.x = TRUE)[,c(1,2,3,5)]
commondf_pairs <- commondf_pairs[order(commondf_pairs$order), ]
commondf_pairs$order <- NULL

names(commondf_pairs)=c("symbol0", "symbol", "DBD0")
commondf_pairs$order=seq(nrow(commondf_pairs))
commondf_pairs2 = merge(commondf_pairs, TF_family_mappings, by = 'symbol', all.x = TRUE)[,c(1,2,3,4,6)]
commondf_pairs2 <- commondf_pairs2[order(commondf_pairs2$order), ]
commondf_pairs2$order <- NULL
commondf_pairs2 <- commondf_pairs2[, c(2,1,3,4)]
names(commondf_pairs2)=c("symbol0", "symbol1", "DBD0", "DBD1")

pair_families <- paste(commondf_pairs2$DBD0,commondf_pairs2$DBD1, sep = "_")

datas <- add_column(datas, Lambert2018_families = pair_families, .before = 4)
unique_Lambertfamilies <- unique(datas$Lambert2018_families)

datas$filename=NULL
datas <- rename(datas, filename = filename_final)

write.table(datas, file="../../PWMs_final/fromYimeng/metadata.csv", row.names = FALSE,sep="\t")
saveRDS(datas, file="RData/fromYimeng.Rds")


