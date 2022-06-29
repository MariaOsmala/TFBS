invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE))


library("tidyverse")

PWM.exps=dir("../../PWMs/")

nro.species=2 #Human and mouse

PWM.summary=data.frame(names=c("Exp", "Species", "Unique.n", "TFs.n", "Family.n"),nrow=length(PWM.exps)*nro.species)

PWM.summary= data.frame(Exp=character(),
                  Species=character(),
                  Unique.n=integer(),
                  TFs.n=integer(),
                  Family.n=integer())


i=2
for (PWM.exp in PWM.exps) {
  print(PWM.exp)
  #PWM.exp=PWM.exps[i]
 
  #Summary from the experiment
  if(length(grep("data",dir(paste0("../../PWMs/", PWM.exp)))) !=0 ){
   PWMs=read_csv(paste0("../../PWMs/", PWM.exp,"/",  dir(paste0("../../PWMs/", PWM.exp))[ grep("data",dir(paste0("../../PWMs/", PWM.exp)))]  )) 
    
    
    species.n=PWMs %>%
      count(Species)
    
    unique.n=nrow(PWMs %>%
                    count(Species, Name)%>%
                    filter(n==1))
    
    uniq.human=nrow(PWMs %>%
      filter(Species=="Homo_sapiens")%>%
      count(Name)%>%
      filter(n==1))
    
    uniq.mouse=nrow(PWMs %>%
                      filter(Species=="Mus_musculus")%>%
                      count(Name)%>%
                      filter(n==1))
    
    
    # TFs.n=PWMs %>%
    #   count(Species, Name)%>%
    #     count(Species)
    
    family.n=PWMs %>%
      count(Species, Family)%>%
      count(Species)
    
    print("Unique match")
    print(unique.n==(uniq.mouse+uniq.human))
    print("Rows match")
    print(unique.n==nrow(PWMs))
    
    
    PWM.summary[i,]=c(PWM.exp, species.n[1,1],  uniq.human, species.n[1,2], family.n[1,2])
    i=i+1
    PWM.summary[i,]=c(PWM.exp, species.n[2,1], uniq.mouse, species.n[2,2], family.n[2,2])
    i=i+1
    
  }else{
    PWM.summary[i,"Exp"]=PWM.exp
    i=i+1
    PWM.summary[i,"Exp"]=PWM.exp
    i=i+1
    
  }
  
  #More info
  TFInfo=read_delim(paste0("../../PWMs/", PWM.exp,"/",  dir(paste0("../../PWMs/", PWM.exp))[ grep("CisBP",dir(paste0("../../PWMs/", PWM.exp)))] , 
                       "/TF_Information.txt" ), delim="\t") 
  
  Motif_IDs<-list()
  for(name in TFInfo$TF_Name){
    Motif_IDs[[name]] <- strsplit(TFInfo[TFInfo$TF_Name==name, "Motif_ID"]$Motif_ID, ",\t|,")[[1]]
  }
  
  
  TFInfo=select(TFInfo, TF_ID, TF_Name, Species_Name, Status, Gene_ID, Family)
  
  TFInfo$Motif_Nro=sapply(Motif_IDs, length)
  
  cisbp <- read_cisbp(system.file("extdata", "cisbp.txt",package = "universalmotif"))
  write_matrix(cisbp, file="../../PWMs/tmp.txt")
  
  #Read motifs, convert to MOODS format
  PWM.list=read_cisbp(paste0("../../PWMs/", PWM.exp,"/",  dir(paste0("../../PWMs/", PWM.exp))[ grep("CisBP",dir(paste0("../../PWMs/", PWM.exp)))] , 
                     "/PWM.txt" ), skip=0)
  #More info
  PWMs=read_delim(paste0("../../PWMs/", PWM.exp,"/",  dir(paste0("../../PWMs/", PWM.exp))[ grep("CisBP",dir(paste0("../../PWMs/", PWM.exp)))] , 
                           "/PWM.txt" ), delim="/t") 
  
}
