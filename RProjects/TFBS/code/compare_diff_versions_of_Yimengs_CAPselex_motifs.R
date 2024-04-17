library(readr)

rm(list=ls())



fromYimeng_version1 <- read_delim("~/projects/TFBS/PWMs_final/fromYimeng/metadata.csv", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE) #22

table(fromYimeng_version1$type)
# 20230420 pfm_composite_new   pfm_spacing_new 
# 14               503               192 


fromYimeng_version2 <- read_delim("~/projects/TFBS/PWMs_final_version2/fromYimeng/metadata.csv", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE) #22

table(fromYimeng_version2$type)
#Batch_spacing_final_v2  BatchOverlap_final_v2 
#298                    971 



column_order=c( "ID", "symbol", "clone","family", "Lambert2018_families", "organism", "study","experiment",          
"ligand", "batch","seed", "multinomial","cycle","representative",  "short",               
"type", "comment",  "filename", "IC", "IC_universal", "length","consensus")  

fromYimeng_version1=fromYimeng_version1[,column_order]
fromYimeng_version2=fromYimeng_version2[,column_order]
                  
#Are the old motifs in the new set?

head(fromYimeng_version1$ID)
head(fromYimeng_version2$ID)

#If I remove the composite_new, spacing_new and 20230420 from the the ID are they still unique

#709 -> 707 unique
ID2=fromYimeng_version1 %>% rowwise() %>%
                    mutate(concatenated = paste(symbol, ligand, batch, seed, multinomial, cycle, sep = "_")) %>% select(concatenated) %>%pull(concatenated)

#which appear multiple times, these are exactly the same motifs, can you find these in the new set
double_ind=which(ID2 %in% names(which(table(ID2)>1)))
#"HOXB13_MEIS3_TGCGAT40NGAT_YRIIII_MATAAANNTGTCA_m1_c3b0u" This pair is not in the set
fromYimeng_version2$ID[which(fromYimeng_version2$symbol=="HOXB13_MEIS3")] #Empty
#"HOXD10_TBX6_TCTACT40NACA_YAAII_NYMRTAAAANAGGTGTNA_m1_c3b0" YES this is composite
tmp=fromYimeng_version2[which(fromYimeng_version2$ID %in% names(which(table(ID2)>1))),]

ID2_all=fromYimeng_version1 %>% rowwise() %>%
             mutate(concatenated = paste(symbol, ligand, batch, seed, multinomial, cycle, sep = "_")) %>% select(concatenated) %>%pull(concatenated)

ID2=unique(fromYimeng_version1 %>% rowwise() %>%
  mutate(concatenated = paste(symbol, ligand, batch, seed, multinomial, cycle, sep = "_")) %>% select(concatenated) %>%pull(concatenated))


length(which(fromYimeng_version2$ID %in% ID2)) #558
length(which(ID2 %in% fromYimeng_version2$ID )) #558
length(which(!(ID2 %in% fromYimeng_version2$ID ))) #149

length(which(!(fromYimeng_version2$ID %in% ID2 ))) #711

#Which are not in the final set
length(which(!(ID2 %in% fromYimeng_version2$ID ))) #149

#708-149, There are 559/709 in the new set, 149 motifs missing 


missing_IDs=ID2[which(!(ID2 %in% fromYimeng_version2$ID ))]

missing_symbols=as.data.frame(do.call(rbind, strsplit(missing_IDs, "_"))[,c(1,2)])
names(missing_symbols)=c("TF1", "TF2")

missing_symbols=missing_symbols %>% mutate(symbol=paste(TF1, TF2, sep="_")) %>% select(symbol) %>% pull(symbol)

missing_IDs[which(missing_symbols %in% fromYimeng_version2$symbol)[1]]
#"DLX6_TBX4_TCTCAT40NTGT_YWII_NTAATKAGGTGTNRN_m1_c3b0"

fromYimeng_version2 %>% filter(symbol==missing_symbols[1]) %>% select(ID) %>% pull(ID)
#"DLX6_TBX4_TCTCAT40NTGT_YWII_NAGGTGTTAATTRN_m1_c3b0"

#Are the missing ones such that there are multiple motifs for the same pair
numbers=lapply( as.list(unique(missing_symbols)), function(x) length(which(fromYimeng_version1$symbol %in% x)) )

#There are multiple instances in the first version, can these be found in the second version, not
length(which(unique(missing_symbols)[which(numbers>1)] %in% fromYimeng_version2$symbol))



read.table((fromYimeng_version1 %>% filter(symbol=="DLX6_TBX4") %>% select(filename) %>% pull(filename))[1])
# DLX6_TBX4_TCTCAT40NTGT_YWII_NAGGTGTTAATTRN_m1_c3b0_short_composite_new This is also in the second set
# DLX6_TBX4_TCTCAT40NTGT_YWII_NTAATKAGGTGTNRN_m1_c3b0_short_composite_new


read.table((fromYimeng_version2 %>% filter(symbol=="DLX6_TBX4") %>% select(filename) %>% pull(filename))[1])

# DLX6_TBX4_TCTCAT40NTGT_YWII_NAGGTGTTAATTRN_m1_c3b0 This is also in the first set