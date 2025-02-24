library(readr)
library(tidyverse)
library(ggplot2)




#In reality 3933 motifs, typo in filename
# metadata <- read_delim("~/projects/TFBS/PWMs_final_version2.2/metadata_representatives.tsv", 
#                        delim = "\t", escape_double = FALSE, 
#                        trim_ws = TRUE)

metadata <- read_delim("/projappl/project_2006203/TFBS/PWMs_final_version2.2/metadata_representatives.tsv", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)



#Metadata for artificial half sites, 7 motifs
metadata_halfsites <- read_delim("/scratch/project_2006203/TFBS/PWMs_final_version2.2/fromYimeng/half-site-motifs-20240426/metadata.csv", 
                                 delim = "\t", escape_double = FALSE, 
                                 trim_ws = TRUE)



metadata_halfsites$filename=paste0("../../", metadata_halfsites$filename) 

#How to combine two data.frames that contain different column names, extract only shared columns 

#Extract shared columns

shared_columns <- intersect(names(metadata), names(metadata_halfsites))

#Combine the two data.frames

metadata <- rbind(metadata[shared_columns], metadata_halfsites[, shared_columns])





##"FOXD3" (this should be)  "IRF6"    "MAFF"    "NFATC1"  "ONECUT1" "SOX11"   "TGIF2LX"

# alternatives=list()
# alternatives[["BHLHB8"]]="BHLHA15"
# alternatives[["PROX2"]]="PROX1"
# alternatives[["NKX2-4"]]="NKX2-3"
# alternatives[["DBX2"]]="DLX2"
# alternatives[["BHLHB5"]]="BHLHA15"
# alternatives[["ZBTB33"]]="ZBTB33" #?
# alternatives[["PBX2"]]="PBX2" #?
# alternatives[["FOXI2"]]="FOXI1"
# alternatives[["gntR"]]="gntR" #?
# alternatives[["NKX2-6"]]="NKX2-3"

#print(unique(missing_monomers))
#"PBX4"   "BHLHB8" "DBX2"   "FOXO3A" "gntR"   "NKX2-4" "PBX2"   "PROX2"  "ZBTB33" "BHLHB5" "NKX2-6"




#capselex=readRDS(file="/Users/osmalama/projects/TFBS/RProjects/TFBS/RData/half_site_recognition_in_capselex_motifs_relaxed_tomtom_version2.2_artificial_halfsites.RDS")
capselex=readRDS(file="/scratch/project_2006203/TFBS/RProjects/TFBS/RData/half_site_recognition_in_capselex_motifs_relaxed_tomtom_version2.2_artificial_halfsites.RDS")
composites=capselex %>% filter(study=="fromYimeng" & type=="composite") #1131
table(composites$switched_order, useNA="always")
#FALSE  TRUE  <NA> 
#  552   564    15 
table(composites$first.monomer=="")
#FALSE  TRUE 
#1116    15 
table(composites$second.monomer=="")
#FALSE 
#1131

#Is this because of the missing monomers and failure by tomtom?
#Are there still homodimers aligned against the composite motif?

table(is.na(composites$first_start)) #15
table(is.na(composites$second_start)) #15
table( is.na(composites$first_start) | is.na(composites$second_start)) #15



tmp=composites[which((is.na(composites$First_Optimal_offset) | is.na(composites$Second_Optimal_offset)) & composites$first.monomer==""),]

unique(do.call(rbind, strsplit(tmp$symbol,"_"))[,1])
#"BHLHB8" "DBX2"   "FOXO3A" "gntR"   "NKX2-4" "PBX2"   "PROX2"  "ZBTB33" These are the ones for which the monomer is missing

#For how many composite motifs, the alignment to the monomers can be defined? The alignment can not be of the spacing type.
#Answer: 1023, spacing type for 93 and no type for 15 due to missing monomer
table(composites$`type computation`)

spacing=composites %>% filter(type_computation=="spacing")

#tmp=composites %>% filter(ID=="FOXK2_ELF1_TTGAAC40NTAGC_YADIIII_NNGAAAACCGAAWNN_m1_c3b0") #Difficult composite, monomers do not align

#How many of the new composite motifs are predicted as spacing motifs


#Reorder the first.monomer and second.monomer columns if switched_order==TRUE
#capselex=capselex[, c(1,2,25:50)]

# Function to compute KL divergence
compute_kl_divergence <- function(P, Q) {
  rownames(P)=NULL
  rownames(Q)=NULL
  kl_div <- 0
  # Iterate over each position
  for (i in 1:ncol(P)) {
    # Compute KL divergence for each nucleotide at this position
    for (j in 1:nrow(P)) {
      if (P[j, i] > 0 && Q[j, i] > 0) { # Avoid log(0) and division by zero
        kl_div <- kl_div + P[j, i] * log2(P[j, i] / Q[j, i])
      }else{
        print("log(0) and division by zero")
      }
    }
  }
  return(kl_div)
}

#"ATF3_ONECUT1_TCATTC40NATT_YJII_NTGACGTCANNNNATCGATN_m1_c3b0"

composites = composites %>% filter(type_computation=="composite") #1023



#save_path="~/projects/TFBS/PWMs_final_version2.2/composite_cores_and_overlapping_monomer_flanks/"
save_path="/scratch/project_2006203/TFBS/PWMs_final_version2.2/composite_cores_and_overlapping_monomer_flanks/"
#save_path="/scratch/project_2006203/TFBS/PWMs_final_version2.2/composite_cores_and_overlapping_monomer_flanks_another_way/"

kmer_scores_all=list()
composites_correct_table_all=data.frame()

for(i in 1:nrow(composites)){
  print(i)
  #kmer_scores,composites_correct_table, 
  load(file=paste0(save_path, composites$ID[i], ".RData"))
  kmer_scores_all[[composites$ID[i]]]=kmer_scores
  composites_correct_table_all=rbind(composites_correct_table_all, composites_correct_table)
  
}




#save.image(file="/scratch/project_2006203/TFBS/RProjects/TFBS/RData/kmer_score_computations.RData")

load(file="/scratch/project_2006203/TFBS/RProjects/TFBS/RData/kmer_score_computations.RData")
#Filter out those cases where the overlap does not correspond to flanks

#tmp=composites_correct_table_all %>% filter(overlap_not_between_flanks==TRUE) 

composites_correct_table_all=composites_correct_table_all %>% filter(overlap_not_between_flanks==FALSE)  #977

#Remove those with overlap longer than 12


#Consider first those with overlap length 1-5

table(composites_correct_table_all$overlap)
#1   2   3   4   5   6   7   8   9  10  11  12  14  17 
#58 157 283 197 120  43  35  43  12  16   5   6   1   1 


composites_correct_table_all$JI_core_TF1_flank=NA
composites_correct_table_all$JI_core_TF2_flank=NA



#Select k-mers with 90% or 99% of the maximum score

#Or select top 10 % or 1% of the kmers


for(overlap in 1:12){
  print(overlap)
  kmer_scores=kmer_scores_all[composites_correct_table_all %>% filter(overlap==overlap) %>%pull(composite) ]
  
  for(composite in names(kmer_scores)){
    
    ind=which(composites_correct_table_all$composite==composite)
    #composite=names(kmer_scores)[1]
    core_max=max(kmer_scores[[composite]]$core)
    core_max_kmers=names(which(kmer_scores[[composite]]$core>=0.9*core_max))
    
    TF1_flank_max=max(kmer_scores[[composite]]$TF1_flank)
    TF1_flank_max_kmers=names(which(kmer_scores[[composite]]$TF1_flank>=0.9*TF1_flank_max))
    
    TF2_flank_max=max(kmer_scores[[composite]]$TF2_flank)
    TF2_flank_max_kmers=names(which(kmer_scores[[composite]]$TF2_flank>=0.9*TF2_flank_max))
    
    JI_core_TF1_flank=length(intersect(core_max_kmers,TF1_flank_max_kmers))/length(union(core_max_kmers,TF1_flank_max_kmers))
    JI_core_TF2_flank=length(intersect(core_max_kmers,TF2_flank_max_kmers))/length(union(core_max_kmers,TF2_flank_max_kmers))
    
    composites_correct_table_all$JI_core_TF1_flank[ind]=JI_core_TF1_flank
    composites_correct_table_all$JI_core_TF2_flank[ind]=JI_core_TF2_flank
    
  }

}

save.image(file="/scratch/project_2006203/TFBS/RProjects/TFBS/RData/kmer_score_computations_JI.RData")
load(file="/scratch/project_2006203/TFBS/RProjects/TFBS/RData/kmer_score_computations_JI.RData")

composites_correct_table_all_long=composites_correct_table_all%>% pivot_longer(cols=starts_with("JI_core"),
                                                names_to="JI_type", 
                                                values_to="JI_value")



x_limits <- c(0.5, 12)
x_breaks <- seq(0, 12, by = 1)

x_hist=ggplot(composites_correct_table_all%>% filter(overlap %in% 1:12), aes_string(x = "overlap")) +
  geom_histogram(bins = 30, alpha = 0.4, position = "identity") +
  # geom_density(alpha = 0.4, size = 0.1) +
  guides(fill = FALSE) + theme_minimal()+
  theme(
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.text.x = element_blank(),   # Remove x-axis tick labels
    axis.ticks.x = element_blank()   # Remove x-axis tick marks
  )+
  #theme_minimal(axis.title.x=element_blank(),plot.margin = margin()) +
  #theme(plot.margin = margin())+                     # Adjust max_size to control circle scaling
  scale_x_continuous(limits = x_limits, breaks = x_breaks)   # Set x-axis limits and custom tick marks
  
  #xlim(x_limits)   # Set x-axis limits for alignment


composites_correct_table_all_long$overlap[which(composites_correct_table_all_long$JI_type=="JI_core_TF2_flank")]=composites_correct_table_all_long$overlap[which(composites_correct_table_all_long$JI_type=="JI_core_TF2_flank")]+0.25
composites_correct_table_all_long$overlap[which(composites_correct_table_all_long$JI_type=="JI_core_TF1_flank")]=composites_correct_table_all_long$overlap[which(composites_correct_table_all_long$JI_type=="JI_core_TF1_flank")]-0.25

p=ggplot(composites_correct_table_all_long %>% filter(overlap %in% seq(0.75,12.25,0.25)), 
       aes(x = overlap, y = JI_value, color = JI_type)) + 
  geom_count( position="identity") +                       # Use geom_count to represent point density with size
 
  scale_size_area(max_size = 10) +     # Adjust max_size to control circle scaling
  theme_minimal() +
  labs(size = "Count")+                     # Adjust max_size to control circle scaling
  scale_x_continuous(limits = x_limits, breaks = x_breaks)   # Set x-axis limits and custom tick marks
  

library("cowplot")
aligned_x_hist <- align_plots(x_hist, p, align = "v",axis="tblr")[[1]]


# Arrange plots
plot_grid(
  aligned_x_hist
  , p
  , ncol = 1
  , nrow = 2
  , rel_heights = c(0.2, 1)
  , rel_widths = c(1, 0.2)
)

#For how many composites either both Jaccard indexes are equal or below 0.5

#composites_correct_table_all %>% filter(overlap<=12 & JI_core_TF1_flank<=0.5 & JI_core_TF2_flank<=0.5) %>%nrow() #717 of 975, 73,5 %

#Jaccard index both
composites_correct_table_all %>% filter(overlap<=12 & JI_core_TF1_flank<0.5 & JI_core_TF2_flank<0.5) %>%nrow() #529 of 975, 54 %

#Jaccard index either one
composites_correct_table_all %>% filter(overlap<=12 & (JI_core_TF1_flank<0.5 | JI_core_TF2_flank<0.5)) %>%nrow() #872 of 975, 89 %

#872-529=343

composites_correct_table_all %>% filter(overlap>12) %>%select(JI_core_TF1_flank) #2
composites_correct_table_all %>% filter(overlap>12) %>%select(JI_core_TF2_flank) #2

#x-y plot

font_size=20

pdf(file="../../Figures/Jaccard_index.pdf", width = 7, height = 7)

#dotted as an alternative line type
p=ggplot(composites_correct_table_all %>% filter(overlap<=12), 
       aes(x=JI_core_TF1_flank, y=JI_core_TF2_flank)) + 
  # Use geom_count to represent point density with size
  geom_segment(aes(x=0, xend=0.5, y=0.5, yend=0.5), linetype = "dashed", color = "red", size = 0.5) +  # Add vertical dashed line
  geom_segment(aes(x=0.5, xend=0.5, y=0, yend=0.5), linetype = "dashed", color = "red", size = 0.5) +  # Add horizontal dashed line
  #geom_point()+
  #geom_jitter(width = 0.01, height = 0.01, alpha = 0.7)+
  geom_count( position="identity") +
  labs(title="Jaccard index for high affinity k-mers", 
       x="Between composite core and TF1 flank", 
       y="Between composite core and TF2 flank", 
      size="Count") +
  scale_size_area(max_size = 10) +
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5, size=font_size, face="bold"),  # Center the title
    axis.title.x = element_text(margin = margin(t = 15), size=font_size-8),  # Add space above the x-axis title
    axis.title.y = element_text(margin = margin(r = 15), size=font_size-8),
    axis.text.x = element_text(size = font_size-8),                      # Increase x-axis tick label size
    axis.text.y = element_text(size = font_size-8),                      # Increase y-axis tick label size
    legend.title = element_text(size = font_size-2),                     # Increase legend title size
    legend.text = element_text(size = font_size-8) 
  )+
  coord_fixed()  # Ensure that the aspect ratio is 1:1
print(p)
dev.off()

#Compute Jaccard index JI(core, TF1) and JI(core, TF2)

#JI(core, TF1)=intersection(core-kmers, TF1-flank-kmers)/union(core-kmers, TF1-flank kmers)

composites_correct_table_all%>% filter(JI_core_TF1_flank==1 & JI_core_TF2_flank==1) %>% count(composite) #42
composites_correct_table_all%>% filter(JI_core_TF1_flank==1 & JI_core_TF2_flank==1) %>% count(overlap)
composites_correct_table_all%>% filter(JI_core_TF1_flank==0.5 & JI_core_TF2_flank==0.5 & overlap==3) %>% pull(composite) 
composites_correct_table_all%>% filter(JI_core_TF1_flank==0 & JI_core_TF2_flank==0.5 & overlap==3) %>% pull(composite)
composites_correct_table_all%>% filter(JI_core_TF1_flank==0.25 & JI_core_TF2_flank==0.25) %>%pull(composite)
composites_correct_table_all%>% filter(JI_core_TF1_flank==0 & JI_core_TF2_flank==0 &overlap==3) %>%pull(composite)
composites_correct_table_all%>% filter(sj_divergence_flank_TF1_core>4 & sj_divergence_flank_TF2_core>1) %>% filter(overlap==1) %>%pull(composite)
composites_correct_table_all%>% filter(sj_divergence_flank_TF1_core<0.3 & sj_divergence_flank_TF2_core>4) %>% pull(composite)
composites_correct_table_all%>% filter(sj_divergence_flank_TF1_core<1 & sj_divergence_flank_TF2_core<1) %>% pull(composite)
composites_correct_table_all%>% filter(sj_divergence_flank_TF1_core>3 & sj_divergence_flank_TF2_core>3) %>% pull(composite)

composites_correct_table_all%>% filter(JI_core_TF1_flank==0.5 & JI_core_TF2_flank==0 & overlap==3) %>% pull(composite) 
#"MSX1_GABPA_TTACAA40NTAC_YYIIII_NCCGGAWRYAATTAN_m1_c3b0"

composites_correct_table_all%>% 
  filter(sj_divergence_flank_TF1_core>1 & sj_divergence_flank_TF2_core>1 & JI_core_TF1_flank==0 & JI_core_TF2_flank==0 & overlap==6) %>%
  pull(composite)



#Divide the kl-divergences by the length of the overlap
#Mean of the mean kl-divergences of both alignments, order and see the most similar and most different


#composites_correct_table_all$kl_divergence_mean_TF1= composites_correct_table_all$kl_divergence_mean_TF1/composites_correct_table_all$overlap
#composites_correct_table_all$kl_divergence_mean_TF2= composites_correct_table_all$kl_divergence_mean_TF2/composites_correct_table_all$overlap

#composites_correct_table_all$kl_divergence_mean=rowMeans( composites_correct_table_all[,c("kl_divergence_mean_TF1","kl_divergence_mean_TF2")] )

#sort accordint to the mean kl-divergence

#composites_correct_table_all=composites_correct_table_all[order(composites_correct_table_all$kl_divergence_mean),]

saveRDS(composites_correct_table_all, file="~/projects/TFBS/RProjects/TFBS/RData/composites_correct_table_all.RDS")
#composites_correct_table_all=readRDS(file="~/projects/TFBS/RProjects/TFBS/RData/composites_correct_table_all.RDS")



#Take the minimum of the kl-divergence between the core and the flanks


composites_correct_table_all$kl_divergence_min_TF1=apply(composites_correct_table_all[,c("kl_divergence_core_flank_TF1", "kl_divergence_flank_TF1_core")],1,min)
composites_correct_table_all$kl_divergence_min_TF2=apply(composites_correct_table_all[,c("kl_divergence_core_flank_TF2", "kl_divergence_flank_TF2_core")],1,min)

#composites_correct_table_all$kl_divergence_min_TF1=composites_correct_table_all$kl_divergence_min_TF1/composites_correct_table_all$overlap
#composites_correct_table_all$kl_divergence_min_TF2=composites_correct_table_all$kl_divergence_min_TF2/composites_correct_table_all$overlap


#Divide the sj_divergence by overlap length
#composites_correct_table_all$sj_divergence_flank_TF1_core=composites_correct_table_all$sj_divergence_flank_TF1_core/composites_correct_table_all$overlap
#composites_correct_table_all$sj_divergence_flank_TF2_core=composites_correct_table_all$sj_divergence_flank_TF2_core/composites_correct_table_all$overlap

#Draw a scatter plot KL-divergences

ggplot(composites_correct_table_all, 
       aes(x=kl_divergence_min_TF1, y=kl_divergence_min_TF2)) + geom_point()+
  theme_minimal()

#Draw a scatter plot JS-divergences

#For how many composites, the jsd is higher than 0.5 for either of the flanks

composites_correct_table_all %>% filter(sj_divergence_flank_TF1_core>0.5 | sj_divergence_flank_TF2_core >0.5) %>% nrow() #552 56.5%
composites_correct_table_all %>% filter(sj_divergence_flank_TF1_core>=0.5 | sj_divergence_flank_TF2_core >=0.5) %>% nrow() #552 56.5%

composites_correct_table_all %>% filter(sj_divergence_flank_TF1_core>0.5 & sj_divergence_flank_TF2_core >0.5) %>% nrow() #236 24.2%
composites_correct_table_all %>% filter(sj_divergence_flank_TF1_core>=0.5 & sj_divergence_flank_TF2_core >=0.5) %>% nrow() #236 24.2%

#Are the different core-flank pairs for the same composites when considering JI or sj

group1 = composites_correct_table_all$JI_core_TF1_flank < 0.5 & composites_correct_table_all$JI_core_TF2_flank < 0.5 & composites_correct_table_all$overlap<=12
group2= (composites_correct_table_all$sj_divergence_flank_TF1_core>0.5 | composites_correct_table_all$sj_divergence_flank_TF2_core >0.5)
library(VennDiagram)

venn_list <- list(
  JI = which(group1),  # Indices of rows in Group A
  JSD = which(group2)
)

# Draw the Venn diagram
venn.plot <- venn.diagram(
  x = venn_list,
  category.names = c("JI", "JSD"),
  main="Venn diagram of the composite motifs that differ from \n the flanks based on different measures",
  filename = NULL,
  fill = c( "#009E73", "orange"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  cat.pos = 0
)

# Plot the Venn diagram
pdf(file="../../Figures/venn.pdf", width = 4, height = 4)
grid::grid.draw(venn.plot)
dev.off()


#what is the Jensen-Shannon for these
motifs=c("EGR1_TEAD4_TGCTAA40NAACA_YJII_NCRCCCACRCAMATTCCN_m2_c3b0",
         "BARHL1_ELK3_TCTATT40NCAT_YAFI_NCGGAANCGTTTAN_m1_c3b0",
         "FOXI1_ETV1_TGATCA40NAGG_YRIIII_NGTAAACAGGAWRYN_m1_c3b0",
         #"GSX2_TBX4_TTGTAG40NTCT_YWII_NAGGTGTTAATTRN_m1_c3b0u",
         "HOXC11_TBX21_TTTACA40NGCT_YYI_NAGGTGTNRTAAAWN_m1_c3b0u",
         "MSX1_GABPA_TTACAA40NTAC_YYIIII_NCCGGAWRYAATTAN_m1_c3b0")

highlight_points=composites_correct_table_all %>% filter(composite %in% motifs)
highlight_points$label=toupper(paste0(do.call(rbind, strsplit(highlight_points$first.monomer, "_"))[,1], ":", do.call(rbind, strsplit(highlight_points$second.monomer, "_"))[,1]))

library("ggforce")
library("ggrepel")

pdf(file="../../Figures/Jensen-Shannon.pdf", width = 7, height = 7)
q=ggplot(composites_correct_table_all, 
       aes(x=sj_divergence_flank_TF1_core, y=sj_divergence_flank_TF2_core)) + geom_point( color="grey", size=2)+
  theme_minimal()+
  labs(title="Jensen-Shannon divergence", 
       x="Between composite core and TF1 flank", 
       y="Between composite core and TF2 flank")+
  geom_point(data = highlight_points, color = "black", size=3) +
  #geom_count( position="identity") +                       # Use geom_count to represent point density with size
  geom_segment(aes(x=0.5, xend=max(sj_divergence_flank_TF1_core), y=0.5, yend=0.5), linetype = "dashed", color = "red", size = 0.5) +  # Add vertical dashed line
  geom_segment(aes(x=0.5, xend=0.5, y=0.5, yend=max(sj_divergence_flank_TF2_core)), linetype = "dashed", color = "red", size = 0.5) +  # Add horizontal dashed line
  #geom_circle(data = highlight_points, aes(x0 = sj_divergence_flank_TF1_core, y0 = sj_divergence_flank_TF2_core, r = 0.2), color = "red", size = 1) +  # Draw circles around selected points
  #geom_text(data = highlight_points, aes(label = label), vjust = -1, color = "black") +  # Add labels above the highlighted points
  geom_text_repel(data = highlight_points, aes(label = label), 
                  color = "black", 
                  size = 6, 
                  box.padding = 1,  # Padding around the label to avoid overlap
                  point.padding = 1,  # Padding around the point to keep the label away from it
                  segment.color = "black") +  # Color of the leader line
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5, size=font_size, face="bold"),  # Center the title
    axis.title.x = element_text(margin = margin(t = 15), size=font_size-8),  # Add space above the x-axis title
    axis.title.y = element_text(margin = margin(r = 15), size=font_size-8),
    axis.text.x = element_text(size = font_size-8),                      # Increase x-axis tick label size
    axis.text.y = element_text(size = font_size-8)
   
  )+
  coord_fixed()  # Ensure that the aspect ratio is 1:1
print(q)
dev.off()

composites_correct_table_all%>% filter(sj_divergence_flank_TF1_core>3 & sj_divergence_flank_TF2_core>3) %>% pull(composite)
composites_correct_table_all%>% filter(sj_divergence_flank_TF1_core>4 & sj_divergence_flank_TF2_core>1) %>% pull(composite)
composites_correct_table_all%>% filter(sj_divergence_flank_TF1_core<0.3 & sj_divergence_flank_TF2_core>4) %>% pull(composite)
composites_correct_table_all%>% filter(sj_divergence_flank_TF1_core<1 & sj_divergence_flank_TF2_core<1) %>% pull(composite)
composites_correct_table_all%>% filter(sj_divergence_flank_TF1_core>3 & sj_divergence_flank_TF2_core>3) %>% pull(composite)

# plot(align_plots(p, q, align = "vh",axis="tb"))
library("gridExtra")
pdf(file="../../Figures/Jensen-Shannon-Jaccard-Index.pdf", width = 16, height = 7)
grid.arrange(p,q+theme(
  plot.margin = margin(5.5, 60, 5.5, 5.5)  # Add extra margin on the right side to balance with p's legend
), ncol=2)
dev.off()
# plot_grid(
#   q
#   , p
#   , ncol = 2
#   , nrow = 1
#   , rel_heights = c(1, 1)
#   , rel_widths = c(1, 1)
# )

#save the data

write.table(composites_correct_table_all, file="../../Figures/Jensen-Shannon-Jaccard-Index-data.tsv", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

