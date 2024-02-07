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
  
  #Include also mouse motifs, do not mind about the upper/lowercase characters in the name
  mono1=monomers[grep(monomer, monomers$symbol, ignore.case=TRUE),]
  
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
  
  #If no motif found, try also dimeric
  if(length(motif1)==0){
    mono1=monomers[grep(monomer, monomers$symbol, ignore.case=TRUE),]
    
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
    
    
  }
  
  list(motif=motif1, representative_bool=representative_bool)
}

monomer_start_end <- function(tomtom_reverse_order, optimal_offset, optimal_orientation, optimal_overlap, dimer_length) {
  
  #tomtom_reverse_order=capselex$First_tomtom_reverse_order[i]
  #optimal_offset=capselex$First_Optimal_offset[i]
  #optimal_orientation=capselex$First_Optimal_Orientation[i]
  #optimal_overlap=capselex$First_Optimal_Overlap
  #dimer_length=capselex$length[i]
  
  #tomtom_reverse_order=capselex$Second_tomtom_reverse_order[i]
  #optimal_offset=capselex$Second_Optimal_offset[i]
  #optimal_orientation=capselex$Second_Optimal_Orientation[i]
  #optimal_overlap=capselex$Second_Optimal_Overlap[i]
  #dimer_length=capselex$length[i]
  
  
  if(tomtom_reverse_order==FALSE){
    offset=optimal_offset
  }else{
    offset=-optimal_offset
  }
  if(optimal_orientation=="-"){
    if(offset<0){
      end=dimer_length
      #If offset -1 no -1 -1+1
      #If offset -2, -1.  -2+1
      start=dimer_length-offset-optimal_overlap+offset+1 
      
    }else{
      end=dimer_length-offset
      start=dimer_length-offset-optimal_overlap+1
    }
    
  }else{
    if(offset<0){
      offset=0
    }
    start=1+offset
    end=offset+optimal_overlap
  }
  #The start needs to be min 1 and end needs to be max the dimer_length of the dimer
  
  start=max(1, start)
  end=min(dimer_length, end)
  
  list(start=start, end=end)
}
