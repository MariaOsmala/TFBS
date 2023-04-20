#Should I increase the threshold for those which map higher than 350 000

# check_ind=which(match_numbers > 300500)
# 
# for(ci in check_ind){
#   print(ci)
#   #ci=check_ind[1]
#   break_bool=FALSE
#   th=6
#   length_orig=match_numbers[ci]
#   while(break_bool==FALSE){
#     #print(th)
#     if( length(which(all_motif_matches_Teemu[[ci]]$score>th)) < length_orig ){
#       
#       length_new=length(which(all_motif_matches_Teemu[[ci]]$score>th))
#       if(length_new-300000<0){
#         if(300000-length_new > 500){
#           th=th-1
#           break_bool=TRUE
#         }else{
#           break_bool=TRUE
#         }
#       }else if(length_new-300000<500){
#         break_bool=TRUE  
#       } 
#     }else{
#       th=th+1
#     }
#   }
#   
#   length_new=length(which(all_motif_matches_Teemu[[ci]]$score>th))
# 
#   if(length_new==length_orig){
#     #do nothing
#   }else{
#     print(ci)
#     #all_motif_matches_Teemu[[ci]]=
#   }
#   
# }


#hist(all_motif_matches_Teemu[[check_ind[1]]]$score)
#length(which(all_motif_matches_Teemu[[check_ind[1]]]$score>7))
#length(which(all_motif_matches_Teemu[[check_ind[1]]]$score>8))
#length(all_motif_matches_Teemu[[check_ind[1]]]$score>9)

#length(which(match_numbers<300000)) 761
#which(is.na(match(names(which(match_numbers<300000)), names(match_numbers_4) ))) 0

match_numbers_4[match(names(which(match_numbers<300000)), names(match_numbers_4) )][1:5]

ind_match_4=match(names(which(match_numbers<300000)), names(match_numbers_4) )
ind_match=match( names(match_numbers_4), names(which(match_numbers<300000)) )

head(which(match_numbers<300000)[ind_match])
head(match_numbers_4)

str(all_motif_matches_Teemu[[which(match_numbers<300000)[ind_match][1] ]])
names(all_motif_matches_Teemu[ which(match_numbers<300000) ])[ ind_match[1] ] 
"ETV6_HT-SELEX_TGAGTG20NGA_AF_CCGGAASCGGAAGTN_2_3_YES"
str(motif_matches_4_Teemu[[1]])
names(motif_matches_4_Teemu)[[1]]
"ETV6_HT-SELEX_TGAGTG20NGA_AF_CCGGAASCGGAAGTN_2_3_YES"

head(which(match_numbers<300000))
head(which(match_numbers<300000)[ ind_match])

for(name in names( which(match_numbers<300000)[ ind_match] ) ){
  print(name)
  all_motif_matches_Teemu[[name]]=motif_matches_4_Teemu[[name]]
  
}

#Just remove first all < 300000 and add motif_matches_4_Teemu

ind_match_3=match(names(which(match_numbers<300000)), names(match_numbers_3) )
ind_match=match( names(match_numbers_3), names(which(match_numbers<300000)) )

head(which(match_numbers<300000)[ind_match])
head(match_numbers_3)

for(name in names( which(match_numbers<300000)[ ind_match] ) ){
  print(name)
  all_motif_matches_Teemu[[name]]=motif_matches_4_Teemu[[name]]
  
}

match_numbers_4=unlist(lapply(motif_matches_4_Teemu, length))

head(names(which(match_numbers<300000)))

match_numbers_4[ match( names(match_numbers_4), names(which(match_numbers<300000)) )][1:5]

head(match_numbers_4)