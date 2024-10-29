# q: what does these function do?

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

#Remove potential dimers, then select the shortest motif, and select the one with highest information content
#Do not care about the representatives
find_monomer <- function(monomer, monomers) {
  
  dimeric_bool=FALSE
  
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
  
  
  if(nrow(mono1)>1 ){
    
    mono1=mono1 %>% filter(length==min(length))
    
    if(length(motif1)>1){
      mono1=mono1 %>%filter(IC==max(IC))
     
    }
    
  }
  
  motif1=mono1 %>% pull(ID)
  
  #If no motif found, try also dimeric
  if(length(motif1)==0){
    dimeric_bool=TRUE
    mono1=monomers[grep(monomer, monomers$symbol, ignore.case=TRUE),]
    
    if(length(which(mono1$experiment=="Methyl-HT-SELEX"))!=0){
      mono1=mono1[-which(mono1$experiment=="Methyl-HT-SELEX"),]
    }
    
    if(nrow(mono1)>1 ){
      
      mono1=mono1 %>% filter(length==min(length))
      
      if(nrow(mono1)>1){
        mono1=mono1 %>%filter(IC==max(IC))
        
      }
      
    }
    #If still
    motif1=mono1 %>% pull(ID)
    
  }
  
  list(motif=motif1, dimeric_bool=dimeric_bool)
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


kmer_position_and_orientation <- function(motif,kmer, monomers) {
  #Find the indexes in the motif which match the k-mer 
  
  pcm= as.matrix(read.table(monomers %>% filter(ID==motif$motif) %>% select(filename) %>% pull()))
  rownames(pcm) <- c("A", "C", "G", "T")
  
  pwm <- TFBSTools::toPWM(pcm, type = "log2probratio", pseudocounts = 1)
  
  #Convert the 6-mer to a DNAString object
  six_mer_seq <- Biostrings::DNAString(kmer)
  
  # Score the 6-mer against each position in the PWM
  scores <- numeric()
  
  for (i in 1:(ncol(pwm) - length(six_mer_seq) + 1)) { #1:4
    motif_submatrix <- pwm[, i:(i + length(six_mer_seq) - 1)]
    score <- sum(motif_submatrix[cbind( match( strsplit( as.character(six_mer_seq), "")[[1]], rownames(motif_submatrix)), 1:ncol(motif_submatrix))])
    scores <- c(scores, score)
  }
  
  
  # Reverse complement of the 6-mer
  rev_six_mer_seq <- Biostrings::reverseComplement(six_mer_seq)
  
  # Score the reverse complement
  rev_scores <- numeric()
  
  for (i in 1:(ncol(pwm) - length(rev_six_mer_seq) + 1)) {
    motif_submatrix <- pwm[, i:(i + length(rev_six_mer_seq) - 1)]
    score <- sum(motif_submatrix[cbind(match(strsplit(as.character(rev_six_mer_seq), "")[[1]], rownames(motif_submatrix)), 1:ncol(motif_submatrix))])
    rev_scores <- c(rev_scores, score)
  }
  
  # Combine scores from both orientations
  combined_scores <- c(scores, rev_scores)
  combined_positions <- c(1:length(scores), 1:length(rev_scores))
  
  # Find the best position considering both orientations
  best_combined_position <- which.max(combined_scores)
  best_combined_score <- max(combined_scores)
  
  # Determine if it's forward or reverse
  if (best_combined_position <= length(scores)) {
    orientation <- "Forward"
    position <- best_combined_position
  } else {
    orientation <- "Reverse"
    position <- best_combined_position - length(scores)
  }
  
  result<-list()
  result$position=position
  result$orientation=orientation
  result$score=best_combined_score
  result$best_forward_score=scores[which.max(scores)]
  result$best_reverse_score=rev_scores[which.max(rev_scores)]
  result$best_forward_position=which.max(scores)
  result$best_reverse_position=which.max(rev_scores)
  result
}


find_monomer_match_kmer <- function(monomer, monomers, six_mer) {
  representative_bool=FALSE
  
  alternatives=list()
  alternatives[["BHLHB8"]]="BHLHA15"
  alternatives[["PROX2"]]="PROX1"
  alternatives[["NKX2-4"]]="NKX2-3"
  alternatives[["DBX2"]]="DLX2"
  alternatives[["BHLHB5"]]="BHLHA15"
  alternatives[["ZBTB33"]]="ZBTB33" #?
  alternatives[["PBX2"]]="PBX2" #?
  alternatives[["FOXI2"]]="FOXI1"
  alternatives[["gntR"]]="gntR" #?
  alternatives[["NKX2-6"]]="NKX2-3"
  
  if(monomer %in% names(alternatives)){
    monomer=alternatives[[monomer]]
  }
  
  
  
  #Include also mouse motifs, do not mind about the upper/lowercase characters in the name
  mono1=monomers[grep(monomer, monomers$symbol, ignore.case=TRUE),]
  
  #There can be representative dimers, remove those
  if(length(which(mono1$type %in% c( "dimeric", "monomeric or dimeric", "dimer of dimers", "putative multimer", "monomer or dimer","trimeric",
                                     "putatively multimeric")))!=0){
    mono1=mono1[-which(mono1$type %in% c( "dimeric", "monomeric or dimeric", "dimer of dimers", "putative multimer", "monomer or dimer","trimeric",
                                          "putatively multimeric")),]
  }
  
  #Consider only shortest
  mono1=mono1 %>%filter(length==min(length))
  
  
  #Which motif gives highest score against the k-mer
  pwms=list()
  for(j in 1:nrow(mono1)){
    pcm= as.matrix(read.table(monomers %>% filter(ID==mono1$ID[j]) %>% select(filename) %>% pull()))
    rownames(pcm) <- c("A", "C", "G", "T")
    pwms[[j]] <- TFBSTools::toPWM(pcm, type = "log2probratio", pseudocounts = 1)
  }
  
  
  #Convert the 6-mer to a DNAString object
  six_mer_seq <- Biostrings::DNAString(six_mer)
  # Reverse complement of the 6-mer
  rev_six_mer_seq <- Biostrings::reverseComplement(six_mer_seq)
  
  scores_forward_list=list()
  scores_reverse_list=list()
  
  scores_list=list()
  positions_list=list()
  j=1
  for(pwm in pwms){
    # Score the 6-mer against each position in the PWM
    scores <- numeric()
    # Score the reverse complement
    rev_scores <- numeric()
    
    for (i in 1:(ncol(pwm) - length(six_mer_seq) + 1)) { #1:4
      motif_submatrix <- pwm[, i:(i + length(six_mer_seq) - 1)]
      score <- sum(motif_submatrix[cbind( match( strsplit( as.character(six_mer_seq), "")[[1]], rownames(motif_submatrix)), 1:ncol(motif_submatrix))])
      scores <- c(scores, score)
    }
    
    for (i in 1:(ncol(pwm) - length(rev_six_mer_seq) + 1)) {
      motif_submatrix <- pwm[, i:(i + length(rev_six_mer_seq) - 1)]
      score <- sum(motif_submatrix[cbind(match(strsplit(as.character(rev_six_mer_seq), "")[[1]], rownames(motif_submatrix)), 1:ncol(motif_submatrix))])
      rev_scores <- c(rev_scores, score)
    }
    
    # Combine scores from both orientations
    scores_forward_list[[j]]=scores
    scores_reverse_list[[j]]=rev_scores
    scores_list[[j]] <- c(scores, rev_scores)
    positions_list[[j]] <- c(1:length(scores), 1:length(rev_scores))
    j=j+1
  }
  
  #Find the maximum value in the entire list
  best_combined_score <- max(unlist(scores_list))
  
  # Find the vectors that contain the overall maximum value
  vectors_with_max <- lapply(scores_list, function(x) which(x == best_combined_score))
  motifs_with_max=which(sapply(vectors_with_max, length) > 0)
  
  # this is list
  best_combined_positions <- vectors_with_max[sapply(vectors_with_max, length) > 0]
  
  # Find the best position considering both orientations
  
  mono1$kmer_score=""
  mono1$kmer_start=""
  mono1$kmer_end=""
  mono1$kmer_orientation=""
  
  j=1
  for(motif_ind in motifs_with_max){
  
    #There can be multiple equally good positions, consider only the first one
    if(length(best_combined_positions[[j]])>1){
      best_combined_positions[[j]]=best_combined_positions[[j]][1]
    } 
    
    # Determine if it's forward or reverse
    if (best_combined_positions[[j]] <= length(scores_forward_list[[motif_ind]])) {
      mono1$kmer_score[motif_ind] <- best_combined_score
      mono1$kmer_orientation[motif_ind] <- "Forward"
      mono1$kmer_start[motif_ind] <- best_combined_positions[[j]]
      mono1$kmer_end[motif_ind] <- best_combined_positions[[j]]+length(six_mer_seq)-1
    } else {
      mono1$kmer_score[motif_ind] <- best_combined_score
      mono1$kmer_orientation[motif_ind] <- "Reverse"
      mono1$kmer_start[motif_ind]<- best_combined_positions[[j]] - length(scores_forward_list[[motif_ind]])
      mono1$kmer_end[motif_ind]<- best_combined_positions[[j]] - length(scores_forward_list[[motif_ind]])+length(six_mer_seq)-1
      
    }
    j=j+1
  
  }
  
  mono1=mono1[motif_ind,]
  
  
  #Which motif has seed that best matches with the kmer
  
  # Convert the 6-mer to a DNAString object
  # six_mer_seq <- Biostrings::DNAString(six_mer)
  # 
  # kmers=mono1$seed
  # 
  # # Initialize a list to store the best matches and their scores
  # best_matches <- list()
  # 
  # # Iterate over each k-mer
  # for (kmer in kmers) {
  #   kmer_seq <- Biostrings::DNAString(kmer)
  #   best_score <- -Inf
  #   best_subseq <- ""
  #   
  #   # Slide the 6-mer along the k-mer
  #   for (i in 1:(nchar(kmer) - nchar(six_mer) + 1)) {
  #     subseq <- Biostrings::subseq(kmer_seq, start=i, width=nchar(six_mer))
  #     score <- pwalign::pairwiseAlignment(six_mer_seq, subseq, type="local", scoreOnly=TRUE)
  #     
  #     if (score > best_score) {
  #       best_score <- score
  #       best_subseq <- as.character(subseq)
  #     }
  #   }
  #   
  #   # Store the best subsequence and score for this k-mer
  #   best_matches[[kmer]] <- list(subsequence = best_subseq, score = best_score)
  # }
  # 
  # # Convert the list to a data frame for easier viewing
  # best_matches_df <- do.call(rbind, lapply(names(best_matches), function(x) {
  #   data.frame(kmer = x, best_subsequence = best_matches[[x]]$subsequence, score = best_matches[[x]]$score)
  # }))
  # 
  # # Sort by the best match score
  # sorted_best_matches_df <- best_matches_df[order(-best_matches_df$score),]
  # 
  # #print(sorted_best_matches_df)
  
  #Best k-mers
  
 # mono1=mono1 %>% filter (seed %in% sorted_best_matches_df$kmer[ which(sorted_best_matches_df$score==max(sorted_best_matches_df$score))])
  
  
  
 
  motif1=mono1 %>% pull(ID)
  
  if(length(motif1)>1){
    mono1=mono1 %>%filter(IC==max(IC))
    motif1=mono1 %>% pull(ID)
  }
  
 
  list(motif=motif1, representative_bool=monomers %>% filter(ID==motif1) %>% select(new_representative) %>% pull(new_representative), 
       kmer_start=as.numeric(mono1$kmer_start), kmer_end=as.numeric(mono1$kmer_end), kmer_orientation=mono1$kmer_orientation)
}



kmer_position_in_motif <- function(motif_filename, kmer) {
  #In what order, orientation, and spacing the kmers are in the motif
  
  pcm= as.matrix(read.table(motif_filename))
  rownames(pcm) <- c("A", "C", "G", "T")
  
  #Add NNN to both sides of the motif
  
    
  pcm=cbind(matrix(as.integer(100),nrow=4,ncol=3),pcm,matrix(as.integer(100),nrow=4,ncol=3))
  colnames(pcm)=NULL
  
  pwm <- TFBSTools::toPWM(pcm, type = "log2probratio", pseudocounts = 1)
  
  
  #Convert the 6-mer to a DNAString object
  kmer_seq <- Biostrings::DNAString(kmer)
  
  # Reverse complement of the 6-mer
  rev_kmer_seq <- Biostrings::reverseComplement(kmer_seq)
  
  
  # Score the 6-mer against each position in the PWM
  scores <- numeric()
  # Score the reverse complement
  rev_scores <- numeric()
  
  for (i in 1:(ncol(pwm) - length(kmer_seq) + 1)) { #1:4
    motif_submatrix <- pwm[, i:(i + length(kmer_seq) - 1)]
    score <- sum(motif_submatrix[cbind( match( strsplit( as.character(kmer_seq), "")[[1]],
                                               rownames(motif_submatrix)), 1:ncol(motif_submatrix))])
    scores <- c(scores, score)
  }
  
  for (i in 1:(ncol(pwm) - length(rev_kmer_seq) + 1)) {
    motif_submatrix <- pwm[, i:(i + length(rev_kmer_seq) - 1)]
    score <- sum(motif_submatrix[cbind(match(strsplit(as.character(rev_kmer_seq), "")[[1]], 
                                             rownames(motif_submatrix)), 1:ncol(motif_submatrix))])
    rev_scores <- c(rev_scores, score)
  }
  
  # Combine scores from both orientations
  scores_combined <- c(scores, rev_scores)
  positions_combined <- c(1:length(scores), 1:length(rev_scores))
  
  #Find the maximum value in the entire list
  best_combined_score <- max(scores_combined)
  best_combined_position <- which(scores_combined == best_combined_score)
  
  # Find the best position considering both orientations
  #If two equally good positions, consider the first one (k-mer likely palindromix)
  if(length(best_combined_position)>1){
    best_combined_position=best_combined_position[1]
  } 
  
  # Determine if it's forward or reverse
  if(best_combined_position <= length(scores) ){
    kmer_score<- best_combined_score
    kmer_orientation <- "Forward"
    kmer_start <- best_combined_position-3
    kmer_end <- best_combined_position+length(kmer_seq)-1-3
  } else {
    kmer_score <- best_combined_score
    kmer_orientation <- "Reverse"
    kmer_start <- best_combined_position - length(scores) -3
    kmer_end<- best_combined_position - length(scores)+length(kmer_seq)-1-3
    
  }
  
  
  
  list(kmer_score=kmer_score,
       kmer_start=kmer_start,
       kmer_end=kmer_end,
       kmer_orientation=kmer_orientation)
}


