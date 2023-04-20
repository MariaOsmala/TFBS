library("rtracklayer")
library("TFBSTools")
library("Biostrings")
library("motifmatchr")
library("GenomicRanges")

representatives=read.table("/scratch/project_2006203/motif-clustering-Viestra-private/metadata/new_representatives.tsv", sep="\t", header=TRUE)
test=gsub(".pfm", "", do.call(rbind, strsplit(representatives$filename,"/"))[,5])
rep_motifs=test[which(representatives$new_representative=="YES")]

#remove _pfm_composite_new and _pfm_spacing_new

rep_motifs=gsub("_pfm_composite_new", "", rep_motifs)
rep_motifs=gsub("_pfm_spacing_new", "", rep_motifs)

#Which proteins contain PROX in their names

PROXes=representatives[grep("PROX", representatives$symbol),] #9

#CAP-selex PROXes

CAP_selex_PROX_HOX=PROXes[PROXes$experiment=="CAP-SELEX",]

CAP_selex_PROX_HOX=CAP_selex_PROX_HOX[-grep("MAFK", CAP_selex_PROX_HOX$symbol),]

profileMatrix=as.matrix(read.table(file=paste0("/scratch/project_2006203/TFBS/",CAP_selex_PROX_HOX$filename[1])))
dimnames(profileMatrix)=list(c("A", "C", "G", "T"))
  
i=1


motifs<-PFMatrixList()
CAP_selex_PROX_HOX=PROXes

for(i in 1:nrow(CAP_selex_PROX_HOX)){
  print(i)
  profileMatrix=as.matrix(read.table(file=paste0("/scratch/project_2006203/TFBS/",CAP_selex_PROX_HOX$filename[i])))
  dimnames(profileMatrix)=list(c("A", "C", "G", "T"))
  
  
motifs[[CAP_selex_PROX_HOX$ID[i]]] <- PFMatrix(ID=CAP_selex_PROX_HOX$ID[i], name=CAP_selex_PROX_HOX$symbol[i], 
                #matrixClass="Zipper-Type", 
                strand="+", 
                bg=c(A=0.25, C=0.25, G=0.25, T=0.25), 
                tags=list(family=CAP_selex_PROX_HOX$Lambert2018.families[i], 
                          #species="10090", 
                          #tax_group="vertebrates", 
                          #medline="7592839", 
                          #type="SELEX", 
                          #ACC="P53762", 
                          #pazar_tf_id="TF0000003",
                          T#FBSshape_ID="11", 
                          #TFencyclopedia_ID="580"
                          ),
                profileMatrix=profileMatrix
)
}



#motifmatchr expects input PWMs to use either natural logarithm or log 2. 
#If the input is a PFM, the TFBSTools toPWM is used for making the PWM 
#with the default psueodcounts of 0.8 and base 2 logarithm. For more control of the pseudocounts, 
#simply use the toPWM function to convert your PFMs prior to calling matchMotifs.

#MOODS
#--bg pA pC pG pT      background distribution for computing thresholds from
#                        p-value with --batch (default is 0.25 for all alleles)
#--ps p                total pseudocount added to each matrix column in log-
#                        odds conversion (default = 0.01)

#  --log-base x          logarithm base for log-odds conversion (default
#                       natural logarithm)

#  --lo-bg pA pC pG pT   background distribution for log-odds conversion
#                        (default is 0.25 for all alleles)

library(TFBSTools)
example_pwms <- do.call(PWMatrixList,lapply(motifs, toPWM, 
                                            pseudocounts = 0.01, bg=c(A=0.25, C=0.25, G=0.25, T=0.25)))

# Get motif matches for example motifs in peaks 

test_GRanges=GRanges(seqnames = Rle( "chr1", rep(1, 1) ),
                     ranges = IRanges(start=214150857, end = 214151185 ) )
motif_ix <- matchMotifs(example_pwms, test_GRanges, genome = "hg19", p.cutoff=5e-02, out="positions") #out=scores


#chr1:214150857-214151185
#sort based on the max score

motif_ix_maxscore<-list()

motif_ix_lnscore<-list()

for(i in 5:length(motif_ix)){
  print(i)
  motif_ix[[i]]=motif_ix[[i]][rev(order(motif_ix[[i]]$score))]
  motif_ix_lnscore[[i]]=motif_ix[[i]]
  motif_ix_lnscore[[i]]$score=motif_ix[[i]]$score/log2(exp(1))
  motif_ix_lnscore[[i]]=motif_ix_lnscore[[i]][motif_ix_lnscore[[i]]$score>=2]
  export.bed(motif_ix_lnscore[[i]], paste0("/scratch/project_2006203/TFBS/Results/Kazenwadel2023/",names(motif_ix)[[i]],".bed") )
  motif_ix_maxscore[[names(motif_ix)[[i]] ]]=as.data.frame(motif_ix_lnscore[[i]][1,])
}

best_motif_matches=do.call(rbind, motif_ix_maxscore)
best_motif_matches$pfm=names(motif_ix_maxscore)

#best_motif_matches$ln_score=best_motif_matches$score/log2(exp(1))

#log.p: if TRUE, probabilities p are given as log(p), natural logarithm.

#The logarithm log_b(x) log_10(x) can be computed from he logarithms of x and b with respect to an arbitrary base k:
#log_b(x)=log_e(x)/log_e(b) #log_10(p)=log_e(p)/log_e(10), phyper returns log_e(p)
p <- phyper(k - 1, n, Nm, m, lower.tail=F, log.p=T) / log(10)


#
#mirs <- GRanges( paste('chr', mirInfo[i,1], sep=''),
#                 IRanges(start=mirInfo[i,4], width=1),
#                 strand=whereverYouKeepTheStrand )
#upstream <- flank(mirs, 200)
#upseqs <- Views(Hsapiens, test_GRanges)


bigbed=system("ls /scratch/project_2006203/TFBS/Results/MOODS_Teemu_processed*/MOODS_bigbed/*PROX*top.bed", intern=TRUE)



bigbed=bigbed[which( gsub("_top.bed", "", do.call(rbind, strsplit(bigbed, "/"))[,8]) %in% gsub(".pfm","", do.call(rbind, strsplit( CAP_selex_PROX_HOX$filename, "/"))[,5]) )]

#GRCh37/hg19 chr. 1: 214,150,857-214,151,185

#test region
test_GRanges=GRanges(seqnames = Rle( "chr1", rep(1, 1) ),
                     ranges = IRanges(start=214150857, end = 214151185 ),
                     strand = Rle(strand("*"), rep(1, 1))  )
seqlevelsStyle(TFBS_GRanges_GRCh37)="NCBI"
seqlevels(test_GRanges)=seqlevels(seqinfo(keepSeqlevels(gtf_GRCh37, c(1), pruning.mode="coarse")))
seqinfo(test_GRanges)=seqinfo(keepSeqlevels(gtf_GRCh37, 1, pruning.mode="coarse"))
seqlevelsStyle(test_GRanges) #NCBI
genome(test_GRanges) <- "GRCh37"


#Convert MOODS hits to Hg19
#These are CGh38 coordinates, convert to CGh37
#ch is hg38ToHg19.over.chain
#utils::download.file('ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz',"/projappl/project_2006203/liftOver/hg38ToHg19.over.chain.gz")
#utils::download.file('ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz',"/projappl/project_2006203/liftOver/hg19ToHg38.over.chain.gz")
#system("gunzip /projappl/project_2006203/liftOver/hg19ToHg38.over.chain.gz")
#system("gunzip /projappl/project_2006203/liftOver/hg38ToHg19.over.chain.gz")

ChainFilePath='/projappl/project_2006203/liftOver/hg38ToHg19.over.chain'
ch = rtracklayer::import.chain(ChainFilePath)
seqlevelsStyle(TFBS_GRanges) <- "UCSC" #UCSC

best_motif_matches$score[1]

gtf_GRCh37 <- readRDS("/projappl/project_2006203/chrom_lengths/gtf_GRCh37.Rds")

for(i in 1:length(bigbed) ){
  TFBS_GRanges=import.bed( bigbed[i])  #logarithm base for log-odds conversion (default natural logarithm)
  print(i)
  #
  #GRCh37/hg19 chr. 1: 214,150,857-214,151,185
  
  #min(TFBS_GRanges$score)
  #[1] 6.915028
  #> max(TFBS_GRanges$score)
  #[1] 8.901016
  
  best_motif_matches$score[1]
  
  TFBS_GRanges_GRCh37=unlist(rtracklayer::liftOver(TFBS_GRanges, ch))
  seqlevelsStyle(TFBS_GRanges_GRCh37)<-"NCBI"
  #seqlevelsStyle(gtf_GRCh37) NCBI
  
  #gr = genes(EnsDb.Hsapiens.v75)
  #gr1 = gr[seqnames(gr) %in% c(1:22, "X", "Y")]
  #seqlevels(gr1) = as.character(unique(seqnames(gr1)))
  #keepSeqlevels(gr, c(1:22, "X", "Y"), pruning.mode="coarse")
  
  seqlevels(TFBS_GRanges_GRCh37)=seqlevels(seqinfo(keepSeqlevels(gtf_GRCh37, c(1:22, "X", "Y"), pruning.mode="coarse")))
  seqinfo(TFBS_GRanges_GRCh37)=seqinfo(keepSeqlevels(gtf_GRCh37, c(1:22, "X", "Y"), pruning.mode="coarse"))
  seqlevelsStyle(TFBS_GRanges_GRCh37) #NCBI
  
  genome(TFBS_GRanges_GRCh37) <- "GRCh37"
  
  
  print(findOverlaps( TFBS_GRanges_GRCh37,test_GRanges, ignore.strand=TRUE))
  
}












