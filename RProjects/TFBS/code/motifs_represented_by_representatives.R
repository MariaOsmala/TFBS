

# readLines("../../../SELEX-Dominating-Set-Analysis/data/Selex_all.dat", 
#          n=3)

# graph_data <- read.table("../../../SELEX-Dominating-Set-Analysis/data/Selex_all.dat", header = FALSE, sep = "\t")

library("readr")
library("tidyverse")
library("igraph")
library("readr")

TFBS_path="/Users/osmalama/projects/TFBS/"
representatives <- read_tsv(paste0(TFBS_path, "PWMs_final_version2.2/","metadata_representatives.tsv"), col_names=TRUE)

n=nrow(representatives)
n*(n-1)/2

df_motifsimilarity=read_tsv(paste0(TFBS_path, "PWMs_final_version2.2/", "sstat.tsv")) #7732278 rows

sim_lt=pivot_wider(df_motifsimilarity, id_cols="Query_ID",
                   names_from="Target_ID", values_from="Ssum")
sim = matrix(0,dim(representatives)[1], dim(representatives)[1])
dim(sim)
diag(sim)=1
tril_indices=lower.tri(sim, diag = FALSE)
ltril_indices=lower.tri(matrix(0, nrow(sim_lt),nrow(sim_lt)), diag = TRUE)

sim[tril_indices]=sim_lt[,-1][ltril_indices]


sim_sym = sim + t(sim) - diag(diag(sim), nrow(sim), ncol(sim))

a=colnames(sim_lt)[2]
b=sim_lt[,"Query_ID"]$Query_ID

a=c(a,b)

rownames(sim_sym)=a
colnames(sim_sym)=a

#In this case, if the distance is less than or equal to the threshold, 
#there is an edge between the nodes (represented by 1); otherwise, there is not (represented by 0). 
#The diagonal is set to 0 to remove self-loops?

#Convert this distance matrix into an adjacency matrix using a threshold:

threshold <- 1.5e-5  # replace with your threshold
adj_mat <- ifelse(sim_sym > threshold, 1, 0)
diag(adj_mat) <- 1  # #each query and target will at least be "close" to themselves in this hash. Or should this be 0
head(adj_mat)

g <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", diag = FALSE)

plot(g)


#representatives <- as.list(read.table("../../../SELEX-Dominating-Set-Analysis/solutions/motifsimilarity_05_min_dom_set_list.txt", header = FALSE, sep = " ")$V1)
representatives <- as.list(read.table("../../SELEX-Dominating-Set-Analysis/solutions/SELEX_all_version2.2_min_dom_set_list.txt", header = FALSE, sep = " ")$V1)

represented_by_representatives=lapply(representatives, function(x, g) neighbors(g, x)$name, g )
names(represented_by_representatives)=representatives

#Remove "nonrepresentatives" that are actually representatives

# Initialize a list to keep track of indices where changes occur
indices <- list()

# Function to remove items and capture indices
process_vectors <- function(x) {
  # Find indices of elements to remove
  ind_to_remove <- which(x %in% unlist(representatives))
  
  # Store indices if not empty
  if (length(ind_to_remove) > 0) {
    return(list(cleaned = x[!x %in% unlist(representatives)], removed_indices = ind_to_remove))
  } else {
    return(list(cleaned = x, removed_indices = NULL))
  }
}

# Apply the function to each element of the list
processed_list <- lapply(represented_by_representatives, process_vectors)

# Extract cleaned lists and indices separately
cleaned_list <- lapply(processed_list, `[[`, "cleaned")

table(sapply(represented_by_representatives, length))
table(sapply(cleaned_list, length))

removed_indices <- lapply(processed_list, `[[`, "removed_indices")

name=names(which(lapply(removed_indices, length)!=0)[1])

removed_indices[[name]]

represented_by_representatives[[name]][removed_indices[[name]]] %in% unlist(representatives)

represented_by_representatives=cleaned_list



saveRDS(represented_by_representatives, file=paste0(TFBS_path,"RProjects/TFBS/RData/represented_by_representatives_version2.2.RDS"))

names(represented_by_representatives)=representatives

# Convert the list to a tibble
df <- tibble::enframe(represented_by_representatives, name = "representative", value = "represented")

# Pad the list elements with NAs to make them of equal length
max_length <- max(map_int(df$represented, length)) #122

padded_list <- map(df$represented, ~ c(.x, rep(NA, max_length - length(.x))))
# Convert the list to a dataframe

tmp=do.call(rbind, padded_list) #1031 x 122

tmp2 <- cbind(data.frame(representative=df$representative), data.frame(tmp))

#Does the represented by representatives contain the representative itself? No, because they were removed above.
non_representatives=unique(as.vector(unlist(tmp2[,2:ncol(tmp2)])))
#remove na
non_representatives=non_representatives[!is.na(non_representatives)]
representatives_in_unrepresentatives=unlist(representatives)[which(unlist(representatives) %in% non_representatives)]

#Convert na values to "" 

tmp2[is.na(tmp2)] <- ""

write.table(tmp2, file=paste0(TFBS_path, "PWMs_final_version2.2/","sstat_represented_by_representatives.tsv"), row.names = FALSE, col.names = FALSE, sep="\t")

stop()


TFBS_path="/Users/osmalama/projects/TFBS/"
sstat_clusters=read.table("~/projects/TFBS/PWMs_final_version2.2/sstat_represented_by_representatives.tsv", header=FALSE, sep="\t")


#Order the table according to the number of motifs in the cluster

row_ord=order(apply(sstat_clusters,1,function(x) sum(x!="")), decreasing=TRUE)

sstat_clusters=sstat_clusters[row_ord,]

paste0("Cluster_", 1:nrow(sstat_clusters))

metadata=read.table("~/projects/TFBS/PWMs_final_version2.2/metadata_representatives.tsv", header=TRUE, sep="\t")

#Add cluster as the fourth column of metadata
#The motif can appear in multiple clusters
#what is the maximum number of clusters a motif can appear in? Seems that 5

#Convert data frame rows to list

sstat_list=unlist(as.list(as.data.frame(t(sstat_clusters))))

sstat_list=sstat_list[-which(sstat_list=="")]
table(table(sstat_list))
# 1    2    3    4    5 
# 3218  541  113   60    1

metadata$cluster1=""
metadata$cluster2=""
metadata$cluster3=""
metadata$cluster4=""
metadata$cluster5=""

length(which( table(sstat_list)==1))
length(which( table(sstat_list)==2))
length(which( table(sstat_list)==3))
length(which( table(sstat_list)==4))
length(which( table(sstat_list)==5))

represented_by_representatives_version2.2 <- readRDS("~/projects/TFBS/RProjects/TFBS/RData/represented_by_representatives_version2.2.RDS")

for(motif in names(which( table(sstat_list)==1))){
  #print(motif)
  #motif=names(which( table(sstat_list)==1))[1]
  #Find out the vector in a list that contains the motif
  index <- as.numeric(which(sapply(represented_by_representatives_version2.2, function(x) motif %in% x)))
  if(length(index)==0){
    index=which(names(represented_by_representatives_version2.2)==motif)
  }
  #What is the representative motif
  
  #What is the order of this cluster
  print(paste0("Cluster_",which(sstat_clusters[,1]==names(represented_by_representatives_version2.2)[index])))
  metadata$cluster1[which(metadata$ID==motif)] = paste0("Cluster_",which(sstat_clusters[,1]==names(represented_by_representatives_version2.2)[index]))
  
}

for(motif in names(which( table(sstat_list)==2))){
  #print(motif)
  #motif=names(which( table(sstat_list)==2))[1]
  #Find out the vector in a list that contains the motif
  index <- as.numeric(which(sapply(represented_by_representatives_version2.2, function(x) motif %in% x)))
  index2=which(names(represented_by_representatives_version2.2)==motif)
  
  index=c(index,index2)
  
  #What is the representative motif
  
  #What is the order of this cluster
  print(paste0("Cluster_",which(sstat_clusters[,1] %in% names(represented_by_representatives_version2.2)[index])))
  metadata[which(metadata$ID==motif), c("cluster1", "cluster2")] = paste0("Cluster_",which(sstat_clusters[,1] %in% names(represented_by_representatives_version2.2)[index]))
}

for(motif in names(which( table(sstat_list)==3))){
  #print(motif)
  #motif=names(which( table(sstat_list)==3))[1]
  #Find out the vector in a list that contains the motif
  index <- as.numeric(which(sapply(represented_by_representatives_version2.2, function(x) motif %in% x)))
  index2=which(names(represented_by_representatives_version2.2)==motif)
  
  index=c(index,index2)
  
  #What is the representative motif
  
  #What is the order of this cluster
  print(paste0("Cluster_",which(sstat_clusters[,1] %in% names(represented_by_representatives_version2.2)[index])))
  metadata[which(metadata$ID==motif), c("cluster1", "cluster2", "cluster3")] = paste0("Cluster_",which(sstat_clusters[,1] %in% names(represented_by_representatives_version2.2)[index]))
}

for(motif in names(which( table(sstat_list)==4))){
  #print(motif)
  #motif=names(which( table(sstat_list)==3))[1]
  #Find out the vector in a list that contains the motif
  index <- as.numeric(which(sapply(represented_by_representatives_version2.2, function(x) motif %in% x)))
  index2=which(names(represented_by_representatives_version2.2)==motif)
  
  index=c(index,index2)
  
  #What is the representative motif
  
  #What is the order of this cluster
  print(paste0("Cluster_",which(sstat_clusters[,1] %in% names(represented_by_representatives_version2.2)[index])))
  metadata[which(metadata$ID==motif), c("cluster1", "cluster2", "cluster3", "cluster4")] = paste0("Cluster_",which(sstat_clusters[,1] %in% names(represented_by_representatives_version2.2)[index]))
}

for(motif in names(which( table(sstat_list)==5))){
  #print(motif)
  #motif=names(which( table(sstat_list)==3))[1]
  #Find out the vector in a list that contains the motif
  index <- as.numeric(which(sapply(represented_by_representatives_version2.2, function(x) motif %in% x)))
  index2=which(names(represented_by_representatives_version2.2)==motif)
  
  index=c(index,index2)
  
  #What is the representative motif
  
  #What is the order of this cluster
  print(paste0("Cluster_",which(sstat_clusters[,1] %in% names(represented_by_representatives_version2.2)[index])))
  metadata[which(metadata$ID==motif), c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5")] = paste0("Cluster_",which(sstat_clusters[,1] %in% names(represented_by_representatives_version2.2)[index]))
}


table(apply(metadata[,c(25:29)],1,function(x) sum(x!="")))



write.table(metadata, file=paste0(TFBS_path, "PWMs_final_version2.2/","metadata_representatives_clusters.tsv"), 
            row.names = TRUE, col.names = TRUE, sep="\t")


########## Plot clusters ##############



library(ggplot2)
library(ggseqlogo)
library(grid)
library(dplyr)
library("readr")
library("tidyverse")
library("igraph")
library("readr")

TFBS_path="/Users/osmalama/projects/TFBS/"
#tomtom data

# Read the TSV file into an R data frame
TFBS_path <- "/Users/osmalama/projects/TFBS/"  # Adjust this path accordingly
tomtom <- read.table(paste0(TFBS_path, 'Results/tomtom_final_version2.2/tomtom_relaxed/tomtom.tsv'), sep="\t", header=TRUE)

#save.image("~/projects/TFBS/RProjects/TFBS/RData/motif_webpage.RData")
load("~/projects/TFBS/RProjects/TFBS/RData/motif_webpage.RData")

#which tomtom_IDs are not in metadata

tomtom_IDs=unique(c(tomtom$Query_ID, tomtom$Target_ID))

tomtom_IDs[which(!(tomtom_IDs %in% metadata$ID))]
metadata$ID[which(!(metadata$ID %in% tomtom_IDs))]

#-# missing from the end of the ID
for(wrong_id in tomtom_IDs[which(!(tomtom_IDs %in% metadata$ID))]){
  correct_ID=which(tomtom$Query_ID==wrong_id)
  tomtom$Query_ID[correct_ID]=paste0(wrong_id,"-")
  
  correct_ID=which(tomtom$Target_ID==wrong_id)
  tomtom$Target_ID[correct_ID]=paste0(wrong_id,"-")
  
}

# Define the cluster you are interested in

for(cluster_nro in 1:1232){
cluster <- paste0('Cluster_',cluster_nro)
print(cluster)

# Filter rows in the data frame based on the specific cluster for each group
df1 <- metadata %>% filter(cluster1 == cluster)
df2 <- metadata %>% filter(cluster2 == cluster)
df3 <- metadata %>% filter(cluster3 == cluster)
df4 <- metadata %>% filter(cluster4 == cluster)
df5 <- metadata %>% filter(cluster5 == cluster)

df=rbind(df1, df2, df3, df4, df5)

results_path="/Users/osmalama/projects/TFBS/PWMs_final_version2.2/motif_webpage/"

motifs = df$ID



rows <- (tomtom$Query_ID %in% motifs) & (tomtom$Target_ID %in% motifs)


pairwise_df = tomtom[rows,]        
if(nrow(pairwise_df)>=1){
#unique(c(pairwise_df$Query_ID, pairwise_df$Target_ID))

#Choose seed which mean abs Optimal_offset in all comparisons is the smallest

# Define a custom function to calculate the mean of absolute values
abs_mean <- function(x) {
  mean(abs(x))
}

seed_motif <- pairwise_df %>%
  group_by(Query_ID) %>%
  summarise(Optimal_offset = abs_mean(Optimal_offset)) %>%
  arrange(Optimal_offset) %>%
  slice(1) %>%
  pull(Query_ID)


rows <- (tomtom$Query_ID== seed_motif) & (tomtom$Target_ID %in% motifs)

pairwise_df = tomtom[rows,] 

w <- nchar(pairwise_df$Target_consensus)
left = min(-pairwise_df['Optimal_offset'])

l_offset = -left - pairwise_df['Optimal_offset']
right = max(l_offset + w)
r_offset = right - w - l_offset

alignment_df <- pairwise_df %>%
  select(-Query_ID, -Optimal_offset, -p.value, -E.value, -q.value, -Overlap, -Query_consensus) %>%
  mutate(w = w, l_offset = l_offset, r_offset = r_offset)

# Rename the columns
colnames(alignment_df) <- c('motif', 'consensus', 'strand', 'w', 'l_offset', 'r_offset')

# Reset the row indices
alignment_df <- alignment_df %>% mutate(row_number = row_number()) %>% select(-row_number)

# Merge alignment_df with df on the 'motif' and 'ID' columns
alignment_df <- alignment_df %>%
  left_join(df, by = c("motif" = "ID"))

# Sort by the 'symbol' column
alignment_df <- alignment_df %>%
  arrange(symbol)

# Reset the row indices
alignment_df <- alignment_df %>% mutate(index = row_number())

n = nrow(alignment_df) #number of PWMs
l = min(alignment_df['l_offset'])
r = max(alignment_df['r_offset'] + alignment_df['w'])
w = r - l #length of final PWM
}else{
  print("error")
  
}
representative_row=alignment_df%>% filter(new_representative=="YES")

#Remove representative from the list
alignment_df_orig=alignment_df

alignment_df=alignment_df %>% filter(new_representative!="YES")


# Plotting dimensions
#w <- 10
#n <- nrow(alignment_df)
width_scale <- 1.0

# Placeholder for plots
plots <- list()
label_plots <- list()

#Save motif widths to an array
motif_widths <- c()

#Draw first the representative motif

#representative motif

motif_id = representative_row['motif']
filename=gsub("../../","",representative_row['filename'])
rc = representative_row['strand'] == '-'
left = representative_row[1,'l_offset']$Optimal_offset+1
width = as.numeric(representative_row['w'])
motif_widths=c(motif_widths, width)
motif_pfm = paste0(TFBS_path, filename)
pcm=as.matrix( read.table(paste0( motif_pfm), header=FALSE) )
dimnames(pcm)=list(c("A", "C", "G", "T"))

# Reverse complement if needed
if (rc) {
  # Reverse the matrix = flip the matrix horizontally
  pcm_rev=pcm[,ncol(pcm):1]
  
  #Substitute nucleotides with complements = flip the orders of A&T and C&G
  
  pcm_rev_comp=pcm_rev
  pcm_rev_comp["A", ]=pcm_rev["T",]
  pcm_rev_comp["T", ]=pcm_rev["A",]
  pcm_rev_comp["C", ]=pcm_rev["G",]
  pcm_rev_comp["G", ]=pcm_rev["C",]
  pcm=pcm_rev_comp
}

pcm_final = matrix(1,4, w)*1000
pcm_final[,left:(left+width-1)] = pcm
dimnames(pcm_final)=list(c("A", "C", "G", "T"))

if((left-1)>=1 ){
  pcm_final[,1:(left-1)]=NA
}

if((left+width)<=w){
  pcm_final[,(left+width):w]=NA  
}

library(dplyr)
library(stringr)

ind_symbols=do.call(rbind, strsplit(alignment_df_orig$symbol,"_"))

# Find the maximum length of vectors in the list
max_length <- max(sapply(strsplit(alignment_df_orig$symbol,"_"), length))

# Function to pad shorter vectors with NA
pad_with_na <- function(vec, max_length) {
  length(vec) <- max_length
  return(vec)
}

# Apply the padding function to all list elements
padded_list <- lapply( strsplit(alignment_df_orig$symbol,"_"), pad_with_na, max_length = max_length)

ind_symbols=do.call(rbind, padded_list)

# Convert the list to a data frame
df <- do.call(rbind, lapply(padded_list, as.data.frame))




tmp=t(apply(ind_symbols,1,function(x) str_replace_all(x,'[\\-0-9]+$', '') ))

if(unique(class(tmp)==c("matrix", "array"))& ncol(tmp)!=2){
  tmp=data.frame(V1=t(tmp), V2=rep(NA, length(tmp)))
}

Concatenated <- ifelse(is.na(tmp[,2]), tmp[,1], paste(tmp[,1], tmp[,2], sep = "_"))

# View the resulting data frame

motif_names=names(sort(table(Concatenated)[ which( table(Concatenated) %in% sort(unique(table(Concatenated)), decreasing=TRUE)[1:2])], decreasing=TRUE))
motif_names=motif_names[1:min(4, length(motif_names))]

two_best=paste(motif_names, collapse="/")

cl_title=paste0(cluster, "\n", 
                # Replace the regex pattern in 'symbol' column and get the most common gene families
                two_best,
                # alignment_df_orig %>% #choose the two most common symbol prefixes
                #   mutate(symbol = str_replace_all(symbol, '[\\-0-9]+$', '')) %>%  # Replace trailing numbers and hyphens
                #   count(symbol) %>%  # Count occurrences of each gene family
                #   top_n(2, wt = n) %>%  # Get the top 2 most common
                #   pull(symbol) %>%  # Extract the symbols
                #   paste(collapse = '/'),  # Concatenate them with '/',
                " ",
                alignment_df_orig %>% #Get the most common DBD family and replace spaces with underscores
                  count(Lambert2018_families) %>%  # Count occurrences of each family
                  top_n(1, wt = n) %>%  # Get the most common
                  pull(Lambert2018_families) %>%  # Extract the family name
                  str_replace_all(' ', '_')  # Replace spaces with underscores
                
)

label_plots[[1]]=ggplot() +
  theme_void()+ #empty plot
  annotate("text", x = 1, y = 1, 
           label = paste(representative_row$symbol, "\n",representative_row$experiment, "\n",representative_row$study, "", sep = ""), 
           hjust = 0, size = 5,colour = "red", fontface="bold") +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),  # Remove plot margins
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  coord_cartesian(clip = "off")  # Allow text to extend outside the plot
  
  
  



plots[[1]]=ggplot() +
  geom_logo(pcm_final, method = "prob") + theme_minimal()+ggseqlogo::theme_logo()+guides(scale="none")+
  labs(y=element_blank())+xlim(c(l-1,r+1))+
  #labs(y = paste(representative_row$symbol, "\n",representative_row$experiment, "\n",representative_row$study, "", sep = ""), 
      
       #)+
  theme(
  #   plot.title = element_text(
  #   hjust = 0.5,         # Center-align the title
  #   vjust = 1,           # Adjust vertical position (lower values move title closer to the plot)
  #   size = 14,           # Adjust the size of the title
  #   face = "bold",       # Make the title bold
  #   margin = margin(t = 10, b = 10)  # Add margin above (t) and below (b) the title
  # ),
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 0, colour = "red", size=10, face="bold"), #fontsize='large', fontweight='bold', color='r'
    axis.text.x = element_blank(),  # Remove x-axis tick labels
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),       # Remove x-axis title
    axis.ticks = element_blank(),         # Remove axis ticks
    plot.margin = unit(c(0, 0, 0, 0), "cm") #top, right, bottom, left
    )   # Remove y-axis tick labels






#remove NA columns
pcm_final=pcm_final[,colSums(is.na(pcm_final))<nrow(pcm_final)]


pdf(paste0("~/projects/TFBS/PWMs_final_version2.2/motif_webpage/clusters/logos/pdf/",cluster,".pdf"),
    #width=(ncol(pcm_final))*0.125,
    #height=0.5
    width=ncol(pcm_final),
    height=4
    )

ggp=ggplot() +
  geom_logo(pcm_final, method = "prob") + theme_minimal()+ggseqlogo::theme_logo()+guides(scale="none")+
  labs(y=NULL, x=NULL)+
  # labs(y = paste(representative_row$symbol, "\n",representative_row$experiment, "\n",representative_row$study, "", sep = ""), 
  #      #title=cl_title
  # )+
  theme(
    axis.title.x = element_blank(),       # Remove x-axis title
    axis.title.y = element_blank(),       # Remove y-axis title
    axis.text.x = element_blank(),        # Remove x-axis tick labels
    axis.text.y = element_blank(),        # Remove y-axis tick labels
    axis.ticks = element_blank(),         # Remove axis ticks
    plot.margin = unit(c(0, 0, 0, 0), "cm")  # Remove plot margins
    #fontsize='large', fontweight='bold', color='r'
    )  
plot(ggp)

dev.off()


png(paste0("~/projects/TFBS/PWMs_final_version2.2/motif_webpage/clusters/logos/png/",cluster,".png"),
    width=(ncol(pcm_final))*0.125*120,
    height=0.5*120
    #width=ncol(pcm_final),
    #height=4
)

ggp=ggplot() +
  geom_logo(pcm_final, method = "prob") + theme_minimal()+ggseqlogo::theme_logo()+guides(scale="none")+
  labs(y=NULL, x=NULL)+
  # labs(y = paste(representative_row$symbol, "\n",representative_row$experiment, "\n",representative_row$study, "", sep = ""), 
  #      #title=cl_title
  # )+
  theme(
        axis.title.x = element_blank(),       # Remove x-axis title
        axis.title.y = element_blank(),       # Remove y-axis title
        axis.text.x = element_blank(),        # Remove x-axis tick labels
        axis.text.y = element_blank(),        # Remove y-axis tick labels
        axis.ticks = element_blank(),         # Remove axis ticks
        plot.margin = unit(c(0, 0, 0, 0), "cm")  # Remove plot margins
      )
plot(ggp)
dev.off()


#POU3F2_MEF2C


# Iterate over each row to create sequence logos
for (i in seq_len(nrow(alignment_df))) {
  print(i)
  row <- alignment_df[i, ]
  motif_id = row['motif']
  filename=gsub("../../","",row['filename'])
  width = as.numeric(row['w'])
  motif_widths=c(motif_widths, width)
  rc = row['strand'] == '-'
  left = row[1,'l_offset']$Optimal_offset+1
  motif_pfm = paste0(TFBS_path, filename)
  
  pcm=as.matrix( read.table(paste0( motif_pfm), header=FALSE) )
  dimnames(pcm)=list(c("A", "C", "G", "T"))
  
  # Reverse complement if needed
  if (row$strand == "-") {
    # Reverse the matrix = flip the matrix horizontally
    pcm_rev=pcm[,ncol(pcm):1]
    
    #Substitute nucleotides with complements = flip the orders of A&T and C&G
    
    pcm_rev_comp=pcm_rev
    pcm_rev_comp["A", ]=pcm_rev["T",]
    pcm_rev_comp["T", ]=pcm_rev["A",]
    pcm_rev_comp["C", ]=pcm_rev["G",]
    pcm_rev_comp["G", ]=pcm_rev["C",]
    pcm=pcm_rev_comp
  }
  
  pcm_final = matrix(1,4, w)*1000
  pcm_final[,left:(left+width-1)] = pcm
  dimnames(pcm_final)=list(c("A", "C", "G", "T"))
  
  pcm_class <- TFBSTools::PFMatrix(strand="+",
                                        bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                                        
                                        profileMatrix=pcm_final
  )
  
  
  # pwm_class<- TFBSTools::toPWM(pcm_class, type="prob", pseudocounts = 0.01, bg=c(A=0.25, C=0.25, G=0.25, T=0.25))
  # pwm_class@profileMatrix[,1:(left-1)]=NA
  # pwm_class@profileMatrix[,(left+width):w]=NA
  
  if((left-1)>=1 ){
    pcm_final[,1:(left-1)]=NA
  }
  
  if((left+width)<=w){
    pcm_final[,(left+width):w]=NA  
  }
  
  
  #Probability matrices
  # print(ggplot() + ggseqlogo::geom_logo(  pcm_final, 
  #                                         method="probability", font="roboto_medium", col_scheme="auto"  ) + 
  #         ggseqlogo::theme_logo()  +
  #         labs(title=pwm_class@name)+ theme(title =element_text(size=8, face='bold'))+ guides(scale="none"))
  # 
  # 
  
  
  label_plots[[i+1]]=ggplot() +
    theme_void()+ #empty plot
    annotate("text", x = 1, y = 1, 
             label = paste(row$symbol, "\n",row$experiment, "\n",row$study, "", sep = ""), 
             hjust = 0, size = 5,colour = "black", fontface="bold") +
    theme(
      plot.margin = unit(c(0, 0, 0, 0), "cm"),  # Remove plot margins
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    coord_cartesian(clip = "off")  # Allow text to extend outside the plot
  
  
  
  
  
  # Create the sequence logo plot
  p=ggplot() +
    geom_logo(pcm_final, method = "prob") + theme_minimal()+ggseqlogo::theme_logo()+guides(scale="none")+
    labs(y=element_blank(), x=NULL)+xlim(c(l-1,r+1))+
    #coord_cartesian(xlim = c(l-2,r+2))+
    #labs(y = paste(row$symbol, "\n",row$experiment, "\n",row$study, "", sep = ""), x=NULL)+
    #labs(y = paste("ALX4", "\n","HT-SELEX", "\n","Yin2017", "", sep = ""), x=NULL)+
    theme(
      axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1,  size=10, face="bold"),
      axis.text.x = element_blank(),  # Remove x-axis tick labels
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),       # Remove x-axis title
      axis.ticks = element_blank(),         # Remove axis ticks
      plot.margin = unit(c(0, 0, 0, 0), "cm")  # Remove plot margins)   # Remove y-axis tick labels
    )
  plots[[i+1]] <- p
}








# Combine all plots
library(gridExtra)

library(grid)


# pdf("~/projects/TFBS/PWMs_final_version2.2/motif_webpage/cluster1.pdf", 
#     width <- (w+2)*0.125+2,
#     height <- (n+1)*0.75+1)
# grid.arrange(grobs = plots, ncol = 1,
#              heights = c(2.5, rep(1, length(plots) - 1)) )
# dev.off()

combined_plots <- mapply(function(label, plot) {
  arrangeGrob(nullGrob(),
              label, 
              nullGrob(),
              plot, 
              nullGrob(),
              ncol = 5, 
              widths = c(1,1, 0.5,max(motif_widths)/4,1))  # Adjust widths as needed
}, label_plots, plots, SIMPLIFY = FALSE)

grid_title <- textGrob(cl_title, gp=gpar(fontsize=16, fontface="bold"))
# Arrange all combined plots in a single layout with two columns
final_plot <- arrangeGrob(
  grid_title, 
  arrangeGrob(grobs = combined_plots, ncol = 1),
  ncol=1,
  heights = c(1/n*10, 10))

# Display the final combined plot
pdf(paste0("~/projects/TFBS/PWMs_final_version2.2/motif_webpage/clusters/",cluster,".pdf"), 
    width <- (w+2)*0.125+2+8,
    height <- (n+1)+1)
grid.draw(final_plot)
dev.off()

png(paste0("~/projects/TFBS/PWMs_final_version2.2/motif_webpage/clusters/",cluster,".png"), 
    width <- ((w+2)*0.125+2+8),
    height <- ((n+1)+1), units="in", res=300)
grid.draw(final_plot)
dev.off()

# grid_title <- textGrob(cl_title, gp=gpar(fontsize=16, fontface="bold"))
# 
# pdf(paste0("~/projects/TFBS/PWMs_final_version2.2/motif_webpage/clusters/",cluster,".pdf"), 
#     width <- (w+2)*0.125+2+16,
#     height <- (n+1)+1)
# # Arrange the plots with a title
# grid.arrange(
#   grid_title,  # The title goes here
#   arrangeGrob(grobs = plots, ncol = 1),  # The grid of plots
#   nrow = 2,  # We have two rows now: one for the title, one for the plots
#   heights = c(1/n*10, 10)  # Adjust the heights, more space for plots than title
# )
# dev.off()
# 
# png(paste0("~/projects/TFBS/PWMs_final_version2.2/motif_webpage/clusters/",cluster,".png"), 
#     width <- ((w+2)*0.125+2+8),
#     height <- ((n+1)*0.75+1), units="in", res=300)
# # Arrange the plots with a title
# grid.arrange(
#   textGrob(cl_title, gp=gpar(fontsize=16, fontface="bold")),  # The title goes here
#   arrangeGrob(grobs = plots, ncol = 1),  # The grid of plots
#   nrow = 2,  # We have two rows now: one for the title, one for the plots
#   heights = c(1/n*10, 10)  # Adjust the heights, more space for plots than title
# )
# dev.off()


}


#Write PWMs to uniprobe
#open(TFBS_path + 'Results/SSTAT/consensus_pwms.uniprobe', 'w') 


#Write clusters to a file, sorted by the number of motifs in the cluster
#AC0095	3	Jolma2015 Jolma2015 Jolma2015	TEAD4_SPDEF TEAD4_SPDEF TEAD4_SPDEF

cluster_list=list()

clusters_df=data.frame()

for(cluster_nro in 1:1232){
  cluster <- paste0('Cluster_',cluster_nro)
  print(cluster)
  
  # Filter rows in the data frame based on the specific cluster for each group
  df1 <- metadata %>% filter(cluster1 == cluster)
  df2 <- metadata %>% filter(cluster2 == cluster)
  df3 <- metadata %>% filter(cluster3 == cluster)
  df4 <- metadata %>% filter(cluster4 == cluster)
  df5 <- metadata %>% filter(cluster5 == cluster)
  
  df=rbind(df1, df2, df3, df4, df5)
  cluster_list[[cluster]]=df
  
  clusters_df=rbind(clusters_df, data.frame(cluster=cluster, n=nrow(df), 
                   study=paste0(df$study, collapse=" "), 
                   symbol=paste0(df$symbol, collapse=" ")))
  
}



write.table(clusters_df,"~/projects/TFBS/PWMs_final_version2.2/motif_webpage/clusters.txt",sep="\t", 
            col.names=FALSE, row.names=FALSE, quote=FALSE )

write.table(metadata,"~/projects/TFBS/PWMs_final_version2.2/motif_webpage/metadata.tsv",sep="\t", 
            col.names=TRUE, row.names=FALSE, quote=FALSE )


# Make the archetype motif MEME file

get_ipython().system('uniprobe2meme /home/local/osmalama/projects/TFBS/Results/tomtom/tomtom/consensus_pwms.uniprobe > /home/local/osmalama/projects/TFBS/Results/tomtom/tomtom/consensus_pwms.meme')


# In[28]:


# Output final motif metadata table

motif_annot_df.to_csv(TFBS_path + 'Results/tomtom/tomtom/metadata.tsv', header=True, index=True, sep='\t')


# In[29]:


motif_annot_df.head()


# In[35]:


# Make the HTML files
get_ipython().system('bash /home/local/osmalama/projects/TFBS/Results/tomtom/tomtom/runall-Osmala.make-html')

