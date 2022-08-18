library(tidyverse)
studies=c("Jolma2013", "Nitta2015", "Jolma2015", "Yin2017", "fromYimeng")
PWMs_metadata<-tibble()

for(s in studies){
  print(s)
  tmp=readRDS(file=paste0("data/",s,".Rds"))
  
  PWMs_metadata=rbind(PWMs_metadata,tmp)
}

###### Interesting TFs and pairs

PWMs_metadata %>%
  filter(symbol=="TCF7")

PWMs_metadata %>%
  filter(symbol=="GLI3_RFX3")


PWMs_metadata %>%
  filter(symbol=="GLI3_FLI1")


##################################################3
#How many motifs for each organism

PWMs_metadata %>%
  count(organism)
# 1 Homo_sapiens  3880
# 2 Mus_musculus   134

PWMs_metadata %>%
  filter(organism=="Homo_sapiens") %>%
  count(symbol)
#1434

PWMs_metadata %>%
  filter(organism=="Homo_sapiens" & experiment=="HT-SELEX") %>%
  count(study)
# 1 Jolma2013   709
# 2 Jolma2015    44
# 3 Nitta2015    10
# 4 Yin2017     865


repTFs=nrow(PWMs_metadata %>% #929
              filter(organism=="Homo_sapiens"  & representative=="YES"))

nonrepTFs=nrow(PWMs_metadata %>% #2196
                 filter(organism=="Homo_sapiens" & representative=="NO") )
#no info
NArepTFs=nrow(PWMs_metadata %>% #755
                filter(organism=="Homo_sapiens" & is.na(representative)) )



#How many motifs for each experiment

PWMs_metadata %>%
  filter(experiment=="HT-SELEX") %>%
  count(organism)

# 1 Homo_sapiens  1628
# 2 Mus_musculus   134

PWMs_metadata %>%
  filter(experiment=="CAP-SELEX") %>%
  count(organism)

#1 Homo_sapiens  1323

PWMs_metadata %>%
  filter(experiment=="Methyl-HT-SELEX") %>%
  count(organism)
#1 Homo_sapiens   929

################ Mouse #######################################
nrow(PWMs_metadata %>%
  filter(organism=="Mus_musculus") %>%
  count(symbol))
#82 unique

#representative motifs
PWMs_metadata %>%
       filter(organism=="Mus_musculus" & representative=="YES") %>%
  count(symbol)
#representative motifs for 31 TF, some TFs have two representative motifs
#nonrepresentative motifs for 62 TFs, some TFs have more than one nonrepresentative motif

#are there some TFs which have either representative and nonrepresentative motif but not both
repTFs=as.vector((PWMs_metadata %>% #31
  filter(organism=="Mus_musculus" & representative=="YES") %>%
  count(symbol))[,1]$symbol)

nonrepTFs=as.vector((PWMs_metadata %>% #62
          filter(organism=="Mus_musculus" & representative=="NO") %>%
          count(symbol))[,1]$symbol)

#no info
NArepTFs=nrow(PWMs_metadata %>% #44
                filter(organism=="Mus_musculus" & is.na(representative)) )



length(which(repTFs %in% nonrepTFs)) #11
length(which(nonrepTFs %in% repTFs)) #11

length(which(!repTFs %in% nonrepTFs)) #20
length(which(!nonrepTFs %in% repTFs)) #51

#11 TFs with both representative and nonrepresentative motifs
#51 TFs have only nonrepresentative motif
#20 TFs have only representative motiv

################ Human ##################################################

############### HT-SELEX
nrow(PWMs_metadata %>%
       filter(organism=="Homo_sapiens" & experiment=="HT-SELEX"))
#1628

#count the number from each study
PWMs_metadata %>%
  filter(organism=="Homo_sapiens" & experiment=="HT-SELEX") %>%
  count(study)

# Jolma2013   709
# Jolma2015    44
# Nitta2015    10
# Yin2017     865

#Unique

uniq=PWMs_metadata %>%
  filter(organism=="Homo_sapiens" & experiment=="HT-SELEX") %>%
  count(symbol)

nrow(uniq) #637
uniq$n=as.factor(uniq$n)

ggplot(uniq, aes(x=n)) + 
  geom_histogram(color="black", fill="white",stat="count")+theme_classic()

#are there some TFs which have either representative and nonrepresentative motif but not both
repTFs=as.vector((PWMs_metadata %>% #260
                    filter(organism=="Homo_sapiens" & experiment=="HT-SELEX" & representative=="YES") %>%
                    count(symbol))[,1]$symbol)

nonrepTFs=as.vector((PWMs_metadata %>% #534
                       filter(organism=="Homo_sapiens" & experiment=="HT-SELEX" & representative=="NO") %>%
                       count(symbol))[,1]$symbol)

#no info
NArepTFs=as.vector((PWMs_metadata %>% #34
                       filter(organism=="Homo_sapiens" & experiment=="HT-SELEX" & is.na(representative)) %>%
                       count(symbol))[,1]$symbol)
#studies with TFs with no representative info
PWMs_metadata %>% #34
  +     filter(organism=="Homo_sapiens" & experiment=="HT-SELEX" & is.na(representative)) %>%
  +     count(study)

# 1 Jolma2015    33
# 2 Nitta2015    10
# 3 Yin2017       1


length(which(repTFs %in% nonrepTFs)) #116
length(which(nonrepTFs %in% repTFs)) #116

length(which(!repTFs %in% nonrepTFs)) #94
length(which(!nonrepTFs %in% repTFs)) #368

#116 TFs with both representative and nonrepresentative motifs
#94 TFs have only nonrepresentative motif
#368 TFs have only representative motiv

#NONUNIQUE are there some TFs which have either representative and nonrepresentative motif but not both
repTFs=nrow(PWMs_metadata %>% #370
                    filter(organism=="Homo_sapiens" & experiment=="HT-SELEX" & representative=="YES"))

nonrepTFs=nrow(PWMs_metadata %>% #1214
                       filter(organism=="Homo_sapiens" & experiment=="HT-SELEX" & representative=="NO") )
#no info
NArepTFs=nrow(PWMs_metadata %>% #44
                      filter(organism=="Homo_sapiens" & experiment=="HT-SELEX" & is.na(representative)) )
#studies with TFs with no representative info
PWMs_metadata %>% #34
  +     filter(organism=="Homo_sapiens" & experiment=="HT-SELEX" & is.na(representative)) %>%
  +     count(study)

# 1 Jolma2015    33
# 2 Nitta2015    10
# 3 Yin2017       1


length(which(repTFs %in% nonrepTFs)) #116
length(which(nonrepTFs %in% repTFs)) #116

length(which(!repTFs %in% nonrepTFs)) #94
length(which(!nonrepTFs %in% repTFs)) #368

#116 TFs with both representative and nonrepresentative motifs
#94 TFs have only nonrepresentative motif
#368 TFs have only representative motiv





################ Human ##################################################

############### CAP-SELEX
nrow(PWMs_metadata %>%
       filter(organism=="Homo_sapiens" & experiment=="CAP-SELEX"))
#1323

#count the number from each study
PWMs_metadata %>%
  filter(organism=="Homo_sapiens" & experiment=="CAP-SELEX") %>%
  count(study)

#1 fromYimeng   705
#2 Jolma2015    618

#Unique

PWMs_metadata %>% #315
  filter(organism=="Homo_sapiens" & experiment=="CAP-SELEX" & study=="fromYimeng") %>%
  count(symbol)

uniq=PWMs_metadata %>%
  filter(organism=="Homo_sapiens" & experiment=="CAP-SELEX") %>%
  count(symbol)

nrow(uniq) #806
uniq$n=as.factor(uniq$n)

ggplot(uniq, aes(x=n)) + 
  geom_histogram(color="black", fill="white",stat="count")+theme_classic()

#are there some TFs which have either representative and nonrepresentative motif but not both
repTFs=as.vector((PWMs_metadata %>% #223
                    filter(organism=="Homo_sapiens" & experiment=="CAP-SELEX" & representative=="YES") %>%
                    count(symbol))[,1]$symbol)

nonrepTFs=as.vector((PWMs_metadata %>% #164
                       filter(organism=="Homo_sapiens" & experiment=="CAP-SELEX" & representative=="NO") %>%
                       count(symbol))[,1]$symbol)

#no info
NArepTFs=as.vector((PWMs_metadata %>% #502 , counts for unique
                      filter(organism=="Homo_sapiens" & experiment=="CAP-SELEX" & is.na(representative)) %>%
                      count(symbol))[,1]$symbol)


#studies with TFs with no representative info
PWMs_metadata %>% 
  filter(organism=="Homo_sapiens" & experiment=="CAP-SELEX" & is.na(representative)) %>%
  count(study)

#  1 fromYimeng   705


length(which(repTFs %in% nonrepTFs)) #72
length(which(nonrepTFs %in% repTFs)) #72

length(which(!repTFs %in% nonrepTFs)) #151
length(which(!nonrepTFs %in% repTFs)) #92

#72+151+92=315

#72 TFs with both representative and nonrepresentative motifs
#151 TFs have only nonrepresentative motif
#92 TFs have only representative motiv


#NONUNIQUE are there some TFs which have either representative and nonrepresentative motif but not both
repTFs=nrow(PWMs_metadata %>% #223
                    filter(organism=="Homo_sapiens" & experiment=="CAP-SELEX" & representative=="YES") )
nonrepTFs=nrow(PWMs_metadata %>% #164
                       filter(organism=="Homo_sapiens" & experiment=="CAP-SELEX" & representative=="NO") )

#no info
NArepTFs=nrow(PWMs_metadata %>% #502 , counts for unique
                      filter(organism=="Homo_sapiens" & experiment=="CAP-SELEX" & is.na(representative)) )


################ Human ##################################################

############### Methyl-HT -SELEX
nrow(PWMs_metadata %>%
       filter(organism=="Homo_sapiens" & experiment=="Methyl-HT-SELEX"))
#929

#count the number from each study
PWMs_metadata %>%
  filter(organism=="Homo_sapiens" & experiment=="Methyl-HT-SELEX") %>%
  count(study)

#Yin2017   929


#Unique


uniq=PWMs_metadata %>%
  filter(organism=="Homo_sapiens" & experiment=="Methyl-HT-SELEX") %>%
  count(symbol)

nrow(uniq) #542
uniq$n=as.factor(uniq$n)

ggplot(uniq, aes(x=n)) + 
  geom_histogram(color="black", fill="white",stat="count")+theme_classic()

#are there some TFs which have either representative and nonrepresentative motif but not both
repTFs=as.vector((PWMs_metadata %>% #151
                    filter(organism=="Homo_sapiens" & experiment=="Methyl-HT-SELEX" & representative=="YES") %>%
                    count(symbol))[,1]$symbol)

nonrepTFs=as.vector((PWMs_metadata %>% #440
                       filter(organism=="Homo_sapiens" & experiment=="Methyl-HT-SELEX" & representative=="NO") %>%
                       count(symbol))[,1]$symbol)

#no info
NArepTFs=as.vector((PWMs_metadata %>% #6 , counts for unique
                      filter(organism=="Homo_sapiens" & experiment=="Methyl-HT-SELEX" & is.na(representative)) %>%
                      count(symbol))[,1]$symbol)



length(which(repTFs %in% nonrepTFs)) #54
length(which(nonrepTFs %in% repTFs)) #54

length(which(!repTFs %in% nonrepTFs)) #97
length(which(!nonrepTFs %in% repTFs)) #386


#54 TFs with both representative and nonrepresentative motifs
#97 TFs have only nonrepresentative motif
#386 TFs have only representative motif

#NONUNIQUE are there some TFs which have either representative and nonrepresentative motif but not both
repTFs=nrow(PWMs_metadata %>% #151
                    filter(organism=="Homo_sapiens" & experiment=="Methyl-HT-SELEX" & representative=="YES") )

nonrepTFs=nrow(PWMs_metadata %>% #440
                       filter(organism=="Homo_sapiens" & experiment=="Methyl-HT-SELEX" & representative=="NO") )

#no info
NArepTFs=nrow(PWMs_metadata %>% #6 , counts for unique
                      filter(organism=="Homo_sapiens" & experiment=="Methyl-HT-SELEX" & is.na(representative)) )

