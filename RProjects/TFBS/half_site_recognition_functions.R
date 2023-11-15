find_common_prefix <- function(strings) {
  # If the vector is empty, return an empty string
  if (length(strings) == 0) return("")
  
  # Find the minimum length of strings
  min_length <- min(nchar(strings))
  
  # Initialize common prefix as an empty string
  common_prefix <- ""
  
  # Iterate through each character position
  for (i in 1:min_length) {
    # Get the i-th character from each string
    chars <- substr(strings, i, i)
    
    # If not all characters are the same, break
    if (length(unique(chars)) > 1) {
      break
    }
    
    # Append the common character to the common prefix
    common_prefix <- paste0(common_prefix, chars[1])
  }
  
  return(common_prefix)
}


find_representative_monomer <- function(monomer, monomers, represented_by_representatives, metadata) {
  representative_bool=FALSE
  mono1=monomers[grep(monomer,monomers$ID),]
  
  #There can be representative dimers, remove those
  if(length(which(mono1$type %in% c( "dimeric", "monomeric or dimeric", "dimer of dimers", "putative multimer", "monomer or dimer","trimeric",
                                     "putatively multimeric")))!=0){
    mono1=mono1[-which(mono1$type %in% c( "dimeric", "monomeric or dimeric", "dimer of dimers", "putative multimer", "monomer or dimer","trimeric",
                                          "putatively multimeric")),]
  }
  
  
  #Exclude METHYL-SELEX
  if(length(which(mono1$experiment=="Methyl-HT-SELEX"))!=0){
    mono1=mono1[-which(mono1$experiment=="Methyl-HT-SELEX"),]
  }
  
  if( length(which(mono1$new_representative=="YES"))!=0 ){
    mono1=mono1[which(mono1$new_representative=="YES"),]
    representative_bool=TRUE
  }
  common_prefix=find_common_prefix(mono1$ID[grep(monomer,mono1$ID)])
  if(representative_bool==FALSE){
    indices <- which(sapply(represented_by_representatives, function(x) any(grepl(common_prefix, x))))
    #can be many, Let's take the shortes, or the one with highest information content?
    mono_representative=names(represented_by_representatives)[indices]
    motif1=metadata %>% filter(ID %in% mono_representative) %>%filter(length==min(length)) %>% pull(ID) #Atoh1
    if(length(motif1)>1){
      motif1=metadata %>% filter(ID %in% mono_representative) %>%filter(length==min(length)) %>%filter(IC==max(IC))%>% pull(ID)
    }
    
  }else{
    motif1=mono1 %>%filter(length==min(length)) %>% pull(ID) 
    if(length(motif1)>1){
      motif1=mono1 %>%filter(length==min(length)) %>%filter(IC==max(IC)) %>% pull(ID)
    }
  }
  
  list(motif=motif1, representative_bool=representative_bool)
}


find_monomer <- function(monomer, monomers, metadata) {
  representative_bool=FALSE
  mono1=monomers[grep(monomer,monomers$ID),]
  
  #There can be representative dimers, remove those
  if(length(which(mono1$type %in% c( "dimeric", "monomeric or dimeric", "dimer of dimers", "putative multimer", "monomer or dimer","trimeric",
                                     "putatively multimeric")))!=0){
    mono1=mono1[-which(mono1$type %in% c( "dimeric", "monomeric or dimeric", "dimer of dimers", "putative multimer", "monomer or dimer","trimeric",
                                          "putatively multimeric")),]
  }
  
  
  #Exclude METHYL-SELEX
  if(length(which(mono1$experiment=="Methyl-HT-SELEX"))!=0){
    mono1=mono1[-which(mono1$experiment=="Methyl-HT-SELEX"),]
  }
  
  if( length(which(mono1$new_representative=="YES"))!=0 ){
    mono1=mono1[which(mono1$new_representative=="YES"),]
    representative_bool=TRUE
  }
  
  if(representative_bool==FALSE & nrow(mono1)>1 ){
    
    mono1=mono1 %>%filter(length==min(length))
    motif1=mono1 %>% pull(ID)
    
    if(length(motif1)>1){
      mono1=mono1 %>%filter(length==min(length)) %>%filter(IC==max(IC))
      motif1=mono1 %>% pull(ID)
    }
    
  }else{
    motif1=mono1 %>%filter(length==min(length)) %>% pull(ID) 
    if(length(motif1)>1){
      motif1=mono1 %>%filter(length==min(length)) %>%filter(IC==max(IC)) %>% pull(ID)
    }
  }
  
  list(motif=motif1, representative_bool=representative_bool)
}


