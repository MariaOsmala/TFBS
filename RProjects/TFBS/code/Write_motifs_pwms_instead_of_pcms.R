

metadata <- read_delim("~/projects/TFBS/PWMs_final/metadata.csv", 
                       +     delim = "\t", escape_double = FALSE, 
                       +     trim_ws = TRUE)

#Write motifs as pwms

for(i in 1:nrow(metadata)){
  print(i)
  pcm=as.matrix( read.table(paste0( metadata$filename[i]), header=FALSE) )
  dimnames(pcm)=list(c("A", "C", "G", "T"))
  #A 1294  802 6919 6919    0  156 6919 1801
  #C 6919 3055 1044  249   41 1461  493 2373
  #G 3132  905  938   56    0  536 1133 1625
  #T 2804 3864  657    0 6919 6919  462 1120
  
  
  pcm_class <- TFBSTools::PFMatrix(ID=metadata$ID[i], name=metadata$symbol[i], 
                                   #matrixClass="Zipper-Type", 
                                   strand="+", 
                                   bg=c(A=0.25, C=0.25, G=0.25, T=0.25), 
                                   tags=list(family=metadata$Lambert2018.families[i], 
                                             species=metadata$organism[i] 
                                             #tax_group="vertebrates", 
                                             #medline="7592839", 
                                             #type="SELEX", 
                                             #ACC="P53762", 
                                             #pazar_tf_id="TF0000003",
                                             #TFBSshape_ID="11", 
                                             #TFencyclopedia_ID="580"
                                   ),
                                   profileMatrix=pcm
  )
  
  
  
  
  pwm_class <- TFBSTools::toPWM(pcm_class, type="prob", pseudocounts = 0.01, bg=c(A=0.25, C=0.25, G=0.25, T=0.25))
  
  write.table(pwm_class@profileMatrix, file=paste0("../../PFMs_space/",metadata$ID[i], ".pfm"), row.names = FALSE, col.names=FALSE, sep=" ") 
  write.table(pwm_class@profileMatrix, file=paste0("../../PFMs_tab/",metadata$ID[i], ".pfm"), row.names = FALSE, col.names=FALSE, sep="\t") 
  
}


