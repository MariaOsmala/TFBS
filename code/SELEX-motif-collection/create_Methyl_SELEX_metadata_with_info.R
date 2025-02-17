#Code here describes how to get the Methyl_SELEX_metadata_with_info.csv file

bisulfite_selex <- read_excel("../../PWMs/Yin2017/NIHMS1670407-supplement-Supplemental_tables_S1-6.xlsx", 
                              sheet="S3 bisulfite-SELEX data", col_names=FALSE, skip=21, n_max=543-21,col_types=c(rep("text",7), rep("numeric",10))) #n_max=8996-22, 

#b table, replicate experiments	
bisulfite_selex_replicate <- read_excel("../../PWMs/Yin2017/NIHMS1670407-supplement-Supplemental_tables_S1-6.xlsx", 
                                        sheet="S3 bisulfite-SELEX data", col_names=FALSE, skip=548, n_max=599-548,col_types=c(rep("text",7), rep("numeric",10)))


colnames(bisulfite_selex)=colnames(bisulfite_selex_replicate)=c("symbol",
                                                                "clone",
                                                                "call",
                                                                "methyl-SELEX_call",
                                                                "bisulfite_call", 
                                                                "comment",
                                                                "seed",
                                                                "start_position_of_CG",
                                                                "enrichment_of_mCG",
                                                                "CG_frequency_cycle_3_normal_seq",
                                                                "TG_frequency_cycle_3_normal_seq",
                                                                "CG_frequency_cycle_3_bisulfite_seq",
                                                                "TG_frequency_cycle_3_bisulfite_seq",
                                                                "CG_frequency_cycle_4_normal_seq",
                                                                "TG_frequency_cycle_4_normal_seq",
                                                                "CG_frequency_cycle_4_bisulfite_seq",
                                                                "TG_frequency_cycle_4_bisulfite_seq")

naind=which(is.na(bisulfite_selex$symbol)) 

# Create a helper function to identify breaks in sequences
breaks <- c(1, diff(naind) != 1)

group_ids <- cumsum(breaks) #146

# Split the original vector by the group identifiers
groups <- split(naind, group_ids)

# Print the groups
print(groups)

for(g in groups){
  #g=groups[[4]]
  bisulfite_selex[g,"symbol"]=bisulfite_selex[g[1]-1,"symbol"]
  
}

naind=which(is.na(bisulfite_selex_replicate$symbol)) 

# Create a helper function to identify breaks in sequences
breaks <- c(1, diff(naind) != 1)

group_ids <- cumsum(breaks) #146

# Split the original vector by the group identifiers
groups <- split(naind, group_ids)

# Print the groups
print(groups)

for(g in groups){
  #g=groups[[4]]
  bisulfite_selex_replicate[g,"symbol"]=bisulfite_selex_replicate[g[1]-1,"symbol"]
  
}

#which symbol has asterisk in its name

asterisk_ind=grep("\\*", bisulfite_selex$symbol)
bisulfite_selex$asterisk="NO"
bisulfite_selex$asterisk[asterisk_ind]="YES"
bisulfite_selex$symbol=gsub("\\*","", bisulfite_selex$symbol)

asterisk_ind=grep("\\*", bisulfite_selex_replicate$symbol)
bisulfite_selex_replicate$asterisk="NO"
bisulfite_selex_replicate$asterisk[asterisk_ind]="YES"
bisulfite_selex_replicate$symbol=gsub("\\*","", bisulfite_selex_replicate$symbol)

bisulfite_selex$replicate="NO"
bisulfite_selex_replicate$replicate="YES"

bisulfite_selex_all=rbind(bisulfite_selex, bisulfite_selex_replicate)


#section c: TFs for which bisulfite data was not obtained
#TF_name		call	methyl-SELEX_call		comment	seed for methyl-SELEX

TFs_without_bisulfite_data <- read_excel("../../PWMs/Yin2017/NIHMS1670407-supplement-Supplemental_tables_S1-6.xlsx", 
                                         sheet="S3 bisulfite-SELEX data", col_names=FALSE, skip=605, n_max=703-605,col_types=c(rep("text",7)))

TFs_without_bisulfite_data =TFs_without_bisulfite_data[,-c(2,5)]
colnames(TFs_without_bisulfite_data) <- c("symbol", "call", "methyl-SELEX_call", "comment", "seed") #seed_methyl-SELEX

naind=which(is.na(TFs_without_bisulfite_data$symbol)) 

# Create a helper function to identify breaks in sequences
breaks <- c(1, diff(naind) != 1)

group_ids <- cumsum(breaks) #146

# Split the original vector by the group identifiers
groups <- split(naind, group_ids)

# Print the groups
print(groups)

for(g in groups){
  #g=groups[[4]]
  TFs_without_bisulfite_data[g,"symbol"]=TFs_without_bisulfite_data[g[1]-1,"symbol"]
  
}





Methyl_SELEX_metadata <- PWMs_metadata %>% filter(experiment=="Methyl-HT-SELEX")
#Match seeds in TFs_without_bisulfite_data and Methyl_SELEX_metadata
Methyl_SELEX_metadata <- Methyl_SELEX_metadata %>%
  left_join(select(TFs_without_bisulfite_data, "symbol", "call","methyl-SELEX_call", "comment", "seed" ), by = c('symbol','seed'))


#All possible outcomes
#Table S5 - Contingency table describing the combined effect of the methyl-preference of a TF, 
#methylation of its site, and biological activity of the TF on the methylated state of the bound CpG

#TFs preference: The binding of TFs to DNA is inhibited (MethylMinus), enhanced (MethylPlus) or not affected (Little effect) by CpG methylation.
#state : The state of the cytosine, non-methylated cytosine (C) or methylated cytosine (mC)
#binding: TFs could (YES) or couldn't (No) bind to subsequences with methylated or non-methylated CpG site(s).
#effect on C:  The effect on cytosine after TFs binding to DNA, 
#inducing methylation on normal cytosine (methylation), 
#demethylation of methylated cytosine (demethylation) 
#or keeping the original state (no effect).

possible_outcomes <- read_excel("../../PWMs/Yin2017/NIHMS1670407-supplement-Supplemental_tables_S1-6.xlsx", 
                                sheet="S5 All possible outcomes", col_names=FALSE, skip=8, n_max=32-8,col_types=c(rep("text",5)))
colnames(possible_outcomes)=c("TFs_preference", "state", "binding", "effect_on_C", "outcome")

#Is the methyl plus/minus info unique for each TF
TFs_unique=unique(bisulfite_selex_all$symbol) #289

methyl_info<-list()

#All methylation info for each TF protein
for(TF in TFs_unique){
  #TF=TFs_unique[1]
  methyl_info[[TF]]=unique(bisulfite_selex_all$call[which(bisulfite_selex_all$symbol==TF)])
  
}

max_length <- max(sapply(methyl_info, length))

# Pad the shorter vectors with NA values
lst_padded <- lapply(methyl_info, function(v) {
  length(v) <- max_length
  return(v)
})

# Combine the vectors by row
result <- do.call(rbind, lst_padded)

#which motifs have unique methylation status or na
na_counts <- apply(result, 1, function(row) sum(is.na(row)))

which(na_counts==3)

df_methyl_info=data.frame(call=result[which(na_counts==3),1])
df_methyl_info$symbol=rownames(df_methyl_info)

Methyl_SELEX_metadata <- Methyl_SELEX_metadata %>%
  left_join(select(df_methyl_info, "symbol", "call", ), by = c('symbol'))


which(na_counts==2)
#KLF17    RFX5    BCL6   BCL6B CREB3L1    E2F2    HSF5     YY1 
#64     114     155     156     162     167     200     232 

which(na_counts==1)
#PAX3 PAX7 
#94   96 
which(na_counts==4)
#empty

bisulfite_selex_all[which(bisulfite_selex_all$symbol=="KLF17"),c("symbol", "clone","call", "seed")]
Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="KLF17",c("symbol", "clone", "seed")]
Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="KLF17",c("call.y")]=bisulfite_selex_all[which(bisulfite_selex_all$symbol=="KLF17")[2:1],"call"] 

bisulfite_selex_all[which(bisulfite_selex_all$symbol=="RFX5"),c("symbol", "clone","call", "seed")]
Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="RFX5",c("symbol", "clone", "seed", "call.y")]
Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="RFX5", "call.y"]=bisulfite_selex_all[which(bisulfite_selex_all$symbol=="RFX5")[c(1,3)],"call"]

bisulfite_selex_all[which(bisulfite_selex_all$symbol=="BCL6"),c("symbol", "clone","call", "seed")]
Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="BCL6",c("symbol", "clone", "seed", "call.y")]
Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="BCL6", "call.y"]=bisulfite_selex_all[which(bisulfite_selex_all$symbol=="BCL6")[c(1,2)],"call"]

bisulfite_selex_all[which(bisulfite_selex_all$symbol=="BCL6B"),c("symbol", "clone","call", "seed")]
Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="BCL6B",c("symbol", "clone", "seed", "call.y")]
Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="BCL6B", "call.y"]=bisulfite_selex_all[which(bisulfite_selex_all$symbol=="BCL6B")[c(1,2)],"call"]

bisulfite_selex_all[which(bisulfite_selex_all$symbol=="CREB3L1"),c("symbol", "clone","call", "seed")]
Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="CREB3L1",c("symbol", "clone", "seed", "call.y")]


collect_match_ind=list()

bisulfite= DNAStringSet(bisulfite_selex_all[which(bisulfite_selex_all$symbol=="CREB3L1"),c("seed")]$seed)
methyl_selex= DNAStringSet(Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="CREB3L1",c("seed")]$seed)

for(i in 1:length(methyl_selex)){
  test=pairwiseAlignment(bisulfite, methyl_selex[i])
  collect_match_ind[[i]]=which(test@score==max(test@score))
}
for(i in 1:length(methyl_selex)){
  Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol=="CREB3L1")[i], "call.y"]=unique(bisulfite_selex_all[which(bisulfite_selex_all$symbol=="CREB3L1"),"call"][collect_match_ind[[i]],"call"])
}

# getOption("digits")
# options(digits = 22)
# 
# pairwiseAlignment(bisulfite, methyl_selex[2])
# pairwiseAlignment(bisulfite, methyl_selex[1])
# pairwiseAlignment(bisulfite, methyl_selex[1])
# pairwiseAlignment(bisulfite, methyl_selex[1])
# pairwiseAlignment(bisulfite, methyl_selex[1])
# 
# 
# query = DNAStringSet( as.vector(concordant %>% filter(symbol==concordant_TF) %>% dplyr::select(seed))$seed)
# target = DNAStringSet(as.vector(Methyl_SELEX_metadata$seed) )
# names(target)=paste0(Methyl_SELEX_metadata$symbol,"_", Methyl_SELEX_metadata$clone)
# 
# Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol == concordant_TF),]
# 
# alignment=pairwiseAlignment(target, query[1])
# scores=alignment@score
# names(scores)=Methyl_SELEX_metadata$symbol
# 
# 
# scores_order=order(scores, decreasing=TRUE)
# scores_sorted=sort(scores, decreasing=TRUE)


bisulfite_selex_all[which(bisulfite_selex_all$symbol=="E2F2"),c("symbol", "clone","call", "seed")]
Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="E2F2",c("symbol", "clone", "seed", "call.y")]

Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="E2F2", "call.y"]=bisulfite_selex_all[which(bisulfite_selex_all$symbol=="E2F2")[c(5,1)],"call"]

bisulfite_selex_all[which(bisulfite_selex_all$symbol=="HSF5"),c("symbol", "clone","call", "seed")]
Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="HSF5",c("symbol", "clone", "seed", "call.y")]

Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="HSF5", "call.y"]=bisulfite_selex_all[which(bisulfite_selex_all$symbol=="HSF5")[c(2,1)],"call"]

bisulfite_selex_all[which(bisulfite_selex_all$symbol=="YY1"),c("symbol", "clone","call", "seed")]
Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="YY1",c("symbol", "clone", "seed", "call.y")]

Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="YY1", "call.y"]=bisulfite_selex_all[which(bisulfite_selex_all$symbol=="YY1")[c(2,1)],"call"]

bisulfite_selex_all[which(bisulfite_selex_all$symbol=="PAX3"),c("symbol", "clone","call", "seed")]
Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="PAX3",c("symbol", "clone", "seed", "call.y")]

collect_match_ind=list()

bisulfite= DNAStringSet(bisulfite_selex_all[which(bisulfite_selex_all$symbol=="PAX3"),c("seed")]$seed)
methyl_selex= DNAStringSet(Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="PAX3",c("seed")]$seed)

for(i in 1:length(methyl_selex)){
  test=pairwiseAlignment(bisulfite, methyl_selex[i])
  collect_match_ind[[i]]=which(test@score==max(test@score))
}

Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="PAX3", "call.y"]=c("Multiple effects", "MethylPlus", "MethylMinus")


bisulfite_selex_all[which(bisulfite_selex_all$symbol=="PAX7"),c("symbol", "clone","call", "seed")]
Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="PAX7",c("symbol", "clone", "seed", "call.y")]

collect_match_ind=list()

bisulfite= DNAStringSet(bisulfite_selex_all[which(bisulfite_selex_all$symbol=="PAX7"),c("seed")]$seed)
methyl_selex= DNAStringSet(Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="PAX7",c("seed")]$seed)

for(i in 1:length(methyl_selex)){
  test=pairwiseAlignment(bisulfite, methyl_selex[i])
  collect_match_ind[[i]]=which(test@score==max(test@score))
}

Methyl_SELEX_metadata[Methyl_SELEX_metadata$symbol=="PAX7", "call.y"]=c("MethylMinus", "Multiple Effects",  
                                                                        "MethylPlus",
                                                                        "MethylMinus",  
                                                                        "Multiple Effects",  
                                                                        "MethylPlus")

#Test which motifs do not have CpG

Methyl_SELEX_metadata$no_cpg=""


dinucleotides <- c("CG", "YG", "MG", "SG", "CK", "CR", "CS")
#   "YK", "YR", "YS", "MK", "MR", "MS", "SK", "SR", "SS")

# Function to check for dinucleotides
check_dinucleotides <- function(sequence, dinucleotides) {
  result=FALSE
  for (dinucleotide in dinucleotides) {
    positions <- gregexpr(pattern = dinucleotide, sequence)[[1]]
    if (positions[1] != -1) {
      #print(paste("Dinucleotide", dinucleotide, "found at positions:", paste(positions, collapse=", ")))
      result=TRUE
      
    }
  }
  result
}

cpg_sites=sapply(as.list(Methyl_SELEX_metadata$consensus), check_dinucleotides, dinucleotides=dinucleotides)

Methyl_SELEX_metadata$no_cpg[which(cpg_sites==FALSE)]="No CpG"



write.table(Methyl_SELEX_metadata,file="tables/Methyl_SELEX_metadata_with_info.tsv",sep="\t", row.names = FALSE)



#Below this some trials, not used
stop()


#Are the TFs in the replicate table subset of the main table
TFs_unique_replicate=unique(bisulfite_selex_replicate$symbol) #30
length(which(TFs_unique_replicate %in% TFs_unique)) #28

TFs_unique_replicate[which(!(TFs_unique_replicate %in% TFs_unique))] #all are included, note the asterisk

#Are the TFs in the TFs_without_bisulfite_data missing from the bisulfite selex table
TFs_unique_without_bisulfite_data=unique(TFs_without_bisulfite_data$symbol) #82
length(which(TFs_unique_without_bisulfite_data %in% TFs_unique)) #0

#match metadata and methylation info
as.vector(unique(PWMs_metadata %>% filter(experiment=="Methyl-HT-SELEX") %>% select(symbol)))$symbol[1] #541

Methyl_SELEX_metadata <- PWMs_metadata %>% filter(experiment=="Methyl-HT-SELEX")
write.table(Methyl_SELEX_metadata, file="tables/Methyl_SELEX_metadata.tsv",sep="\t", row.names = FALSE)
write.table(bisulfite_selex,file="tables/bisulfite_selex.tsv",sep="\t", row.names = FALSE)
write.table(bisulfite_selex_replicate,file="tables/bisulfite_selex_replicate.tsv",sep="\t", row.names = FALSE)
write.table(TFs_without_bisulfite_data,file="tables/TFs_without_bisulfite_data.tsv",sep="\t", row.names = FALSE)

methyl_metadata_names=c("methyl_symbol",
                        "methyl_clone",
                        "methyl_call",                              
                        "methyl-SELEX_call",
                        "bisulfite_call",
                        "methyl_comment",                           
                        "methyl_seed",
                        "start_position_of_CG",
                        "enrichment_of_mCG",                 
                        "CG_frequency_cycle_3_normal_seq",
                        "TG_frequency_cycle_3_normal_seq",  
                        "CG_frequency_cycle_3_bisulfite_seq",
                        "TG_frequency_cycle_3_bisulfite_seq",
                        "CG_frequency_cycle_4_normal_seq",    "TG_frequency_cycle_4_normal_seq",   
                        "CG_frequency_cycle_4_bisulfite_seq", "TG_frequency_cycle_4_bisulfite_seq")

methyl_metadata=as.data.frame(matrix("", nrow(Methyl_SELEX_metadata), length(methyl_metadata_names)))
names(methyl_metadata)=methyl_metadata_names

Methyl_SELEX_metadata=cbind(Methyl_SELEX_metadata, methyl_metadata)

concordant=as_tibble(bisulfite_selex[which(bisulfite_selex$comment=="concordant"),])

bisulfite_selex
target = DNAStringSet(as.vector(Methyl_SELEX_metadata$seed) )

collect_match_ind=c()

Methyl_SELEX_metadata$seed_match_score=NA

#TFs_unique=unique(concordant$symbol)
#TFs_unique=TFs_unique[-which(is.na(TFs_unique))]

TFs_unique=unique(bisulfite_selex$symbol)
TFs_unique=TFs_unique[-which(is.na(TFs_unique))]

for(concordant_TF in TFs_unique){
  #concordant_TF=unique(concordant$symbol)[1]
  print(concordant_TF)
  
  #concordant_table=concordant %>% filter(symbol==concordant_TF)
  #query = DNAStringSet( as.vector(concordant_table %>% dplyr::select(seed))$seed)
  
  concordant_table=bisulfite_selex %>% filter(symbol==concordant_TF)
  
  query = DNAStringSet( as.vector(concordant_table %>% dplyr::select(seed))$seed)
  
  
  for(i in 1:length(query)){
    print(i)
    q=query[i]
    
    alignment=pairwiseAlignment(target, q)
    scores=alignment@score
    names(scores)=Methyl_SELEX_metadata$symbol
    scores_order=order(scores, decreasing=TRUE)
    scores_sorted=sort(scores, decreasing=TRUE)
    scores_sorted[which(names(scores_sorted)==concordant_TF)]
    
    
    match_ind=scores_order[which(names(scores_sorted)==concordant_TF)][1]
    if(!is.na(match_ind)){
      if(scores_sorted[which(names(scores_sorted)==concordant_TF)][1]>0){ #Score needs to be higher than 0
        
        #Is the match_ind already filled
        if(!(match_ind %in% collect_match_ind)){
          
          #Do the clone match
          if(Methyl_SELEX_metadata$clone[match_ind]==concordant_table$clone[i]){
            Methyl_SELEX_metadata$seed_match_score[match_ind]=scores_sorted[which(names(scores_sorted)==concordant_TF)][1] #score
            Methyl_SELEX_metadata[match_ind, names(methyl_metadata)]=concordant_table[i,]
            
          }else{
            print("clone does not match")
          }
        }else{
          print("already filled")
        }
        
        collect_match_ind=c(collect_match_ind, match_ind)
        
      }
    }
  }
  
}

write.table(Methyl_SELEX_metadata, file="tables/Methyl_SELEX_metadata.tsv",sep="\t", row.names = FALSE)


concordant_TF=unique(concordant$symbol)[1]
seq1 = DNAString(as.vector(concordant %>% filter(symbol==concordant_TF) %>% dplyr::select(seed))$seed)

BiocManager::install("DECIPHER")

library(DECIPHER)

collect_match_ind=c()

sequence_list=c( as.vector(concordant %>% filter(symbol==concordant_TF) %>% dplyr::select(seed))$seed, as.vector(Methyl_SELEX_metadata$seed))

sequences = DNAStringSet( sequence_list )

query = DNAStringSet( as.vector(concordant %>% filter(symbol==concordant_TF) %>% dplyr::select(seed))$seed)
target = DNAStringSet(as.vector(Methyl_SELEX_metadata$seed) )
names(target)=paste0(Methyl_SELEX_metadata$symbol,"_", Methyl_SELEX_metadata$clone)

Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol == concordant_TF),]

alignment=pairwiseAlignment(target, query[1])
scores=alignment@score
names(scores)=Methyl_SELEX_metadata$symbol


scores_order=order(scores, decreasing=TRUE)
scores_sorted=sort(scores, decreasing=TRUE)

scores_order[which(names(scores_sorted)==concordant_TF)] #168(this is already taken), 169

match_ind=scores_order[which(names(scores_sorted)==concordant_TF)][1]

collect_match_ind=c(collect_match_ind, match_ind)

Methyl_SELEX_metadata[match_ind), 
names(methyl_metadata)]=concordant[]


alignment=pairwiseAlignment(target, query[2])
scores=alignment@score
names(scores)=Methyl_SELEX_metadata$symbol


scores_order=order(scores, decreasing=TRUE)
scores_sorted=sort(scores, decreasing=TRUE)

#Remove those that are already in the collect_match_ind

match_ind=scores_order[which(names(scores_sorted)==concordant_TF)] #168(this is already taken), 169

match_ind=match_ind[-which(match_ind %in% collect_match_ind)][1]

collect_match_ind=c(collect_match_ind, match_ind)




methyl_minus=bisulfite_selex[which(bisulfite_selex$comment=="concordant" & bisulfite_selex$call=="MethylMinus"),]

methyl_minus_TF=unique(methyl_minus$symbol)

which(!(Methyl_SELEX_metadata$symbol %in% methyl_minus_TF))

for( TF in methyl_minus_TF){
  TF=methyl_minus_TF[3]
  print(length(which(Methyl_SELEX_metadata$symbol==TF) )) #This can be 0,1,2,3,4
  #If two, they are from different clones
}

methyl_minus[methyl_minus$symbol==TF,c("seed","start_position_of_CG")]
methyl_minus[methyl_minus$symbol==TF,]

NRTGAYGTYAYN                    6
NRTGAYGTGGYN                    6
NGCCACRTCAYN YKCCACRTCAYC
NRTGAYGTCAYN GRTGAYGTCAYN
NGCCACGTGKMN NGCCACGTGKMN


seq1 = DNAString("NRTGAYGTYAYN")
seq2 = DNAString("NRTGAYGTGGYN")

seq3 = DNAString("NGCCACRTCAYN")
seq4 = DNAString("NRTGAYGTCAYN")
seq5 = DNAString("NGCCACGTGKMN")

# Perform pairwise alignment
pairwiseAlignment(seq1, seq3)
pairwiseAlignment(seq1, seq4) #Highest score
pairwiseAlignment(seq1, seq5)

pairwiseAlignment(seq2, seq3)
pairwiseAlignment(seq2, seq4) #Highest score
pairwiseAlignment(seq2, seq5)


# Print alignment
print(alignment)
print(alignment2)


# Calculate percentage identity
matches = sum(matchPattern(alignment, alignment))
total = width(alignment)
percentage_identity = matches / total * 100

print(paste("Percentage identity: ", percentage_identity))



bisulfite_selex %>% filter(symbol==TF,) #FL eDBD
PWMs_metadata[which(PWMs_metadata$symbol==TF),]
Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==TF),c("seed", "consensus")] #eDBD FL
Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==TF),] #eDBD FL

TF=methyl_minus_TF[2]
Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==TF), names(methyl_metadata)]=methyl_minus[methyl_minus$symbol==TF,]

TF=methyl_minus_TF[3] #Both tables have two hits, eDBD and FL, sort methyl_minus based on the order in metadata
Methyl_SELEX_metadata$clone[which(Methyl_SELEX_metadata$symbol==TF)]

Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==TF), 
                      names(methyl_metadata)]=methyl_minus[which(methyl_minus$symbol==TF)[match(Methyl_SELEX_metadata$clone[which(Methyl_SELEX_metadata$symbol==TF)], 
                                                                                                methyl_minus$clone[methyl_minus$symbol==TF] )],]
#Same info twice, the seeds are equal at the interesting positions, is this even CpG
TF=methyl_minus_TF[4]
Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==TF), names(methyl_metadata)]=rbind(methyl_minus[methyl_minus$symbol==TF,],methyl_minus[methyl_minus$symbol==TF,])

TF=methyl_minus_TF[5] #Same situation as above
Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==TF), names(methyl_metadata)]=rbind(methyl_minus[methyl_minus$symbol==TF,],methyl_minus[methyl_minus$symbol==TF,])

TF=methyl_minus_TF[6]
#Which seeds in metadata are of same length?


same_length=Methyl_SELEX_metadata$length[which(Methyl_SELEX_metadata$symbol==TF)]==unique(nchar(methyl_minus[methyl_minus$symbol==TF,"seed"]$seed))

Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==TF)[which(same_length)],names(methyl_metadata)]=methyl_minus[which(methyl_minus$symbol==TF)[match(Methyl_SELEX_metadata$clone[which(Methyl_SELEX_metadata$symbol==TF)], 
                                                                                                                                                            methyl_minus$clone[methyl_minus$symbol==TF] )],]
TF=methyl_minus_TF[7]


#GTYAYGTGAY

NRTGAYGTYAYN
NRTGAYGTGGYN
PWMs_list[[which(PWMs_metadata$symbol==methyl_minus_TF)[1]]]
PWMs_list[[which(PWMs_metadata$symbol==methyl_minus_TF)[2]]]



PWMs_metadata[which(PWMs_metadata$symbol==methyl_minus_TF),]

Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==methyl_minus_TF),]
#RTCAYGTGMN
#NNNNCGNNNN

NGCCACRTCAYN
NRTGAYGTCAYN
NGCCACGTGKMN


TF=unique((PWMs_metadata %>% filter(experiment=="Methyl-HT-SELEX") %>% select(symbol))$symbol)[5]

#IUPAC nucleotide codes https://www.bioinformatics.org/sms/iupac.html


Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==bisulfite_selex$symbol[1]),"seed"]
"CYAATTAN"
"NTCGTTAN"
"CYAATTAN"
"NTCGTTAN"

bisulfite_selex[1,"seed"]
TTAAYGAW

bisulfite_selex[3,]
NGTYGTAAAWN
bisulfite_selex[3,"seed"]
NGTYGTAAAAN

Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==bisulfite_selex$symbol[2]),]
"NGTCGTAAANN"
"NGYAATAAAAN"


for(TF in unique((PWMs_metadata %>% filter(experiment=="Methyl-HT-SELEX") %>% select(symbol))$symbol)[1:5]){
  print("TF")
  print(TF)
  print("Methyl_SELEX")
  print(Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==TF),])
  print(Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==TF),"seed"])
  print("bisulfite_SELEX")
  print(bisulfite_selex[which(bisulfite_selex$symbol==TF),])
  print(bisulfite_selex[which(bisulfite_selex$symbol==TF),"seed"])
  print("bisulfite_SELEX_replicate")
  print(bisulfite_selex_replicate[which(bisulfite_selex_replicate$symbol==TF),])
  print(bisulfite_selex_replicate[which(bisulfite_selex_replicate$symbol==TF),"seed"])
  print("TFs_without_bisulfite_data")
  print(TFs_without_bisulfite_data[which(TFs_without_bisulfite_data$symbol==TF),])
  #  print(TFs_without_bisulfite_data[which(TFs_without_bisulfite_data$symbol==TF),"seed"])
  
  #PROX1: seed matches with bisulfite selex
  
  tmp=Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==TF),] 
  nrow(tmp)
  seed1=Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==TF),"seed"] 
  seed2=bisulfite_selex[which(bisulfite_selex$symbol==TF),"seed"]$seed 
  if(seed1==seed2) #Seed does not match
    
    #MSX2
    tmp=Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==TF),] 
  nrow(tmp)
  seed1=Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==TF),"seed"] #two
  seed2=bisulfite_selex[which(bisulfite_selex$symbol==TF),"seed"]$seed #only one
  if(seed1==seed2) #Seed does not match
    
    #MSX1
    tmp=Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==TF),] 
  nrow(tmp)
  seed1=Methyl_SELEX_metadata[which(Methyl_SELEX_metadata$symbol==TF),"seed"] #two
  tmp2=bisulfite_selex[which(bisulfite_selex$symbol==TF),]
  seed2=bisulfite_selex[which(bisulfite_selex$symbol==TF),"seed"]$seed #only one
  if(seed1==seed2) #Seed does not match
    
}

PWMs_metadata[which(PWMs_metadata$symbol==as.vector(unique(PWMs_metadata %>% filter(experiment=="Methyl-HT-SELEX") %>% select(symbol)))$symbol[1]), ]
bisulfite_selex[which(bisulfite_selex$symbol==as.vector(unique(PWMs_metadata %>% filter(experiment=="Methyl-HT-SELEX") %>% select(symbol)))$symbol[1]), "seed"]

test=PWMs_metadata %>% filter(experiment=="Methyl-HT-SELEX") %>% left_join(bisulfite_selex, by=c("seed"))

flights2 %>% left_join(planes, by = "tailnum")


which(PWMs_metadata$experiment=="Methyl-HT-SELEX")


