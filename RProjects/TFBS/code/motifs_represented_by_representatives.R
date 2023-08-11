

# readLines("../../../SELEX-Dominating-Set-Analysis/data/Selex_all.dat", 
          n=3)

# graph_data <- read.table("../../../SELEX-Dominating-Set-Analysis/data/Selex_all.dat", header = FALSE, sep = "\t")

library("readr")
library("tidyverse")
library("igraph")
library("readr")

TFBS_path="/Users/osmalama/projects/TFBS/"
representatives <- read_tsv(paste0(TFBS_path, "PWMs_final/","metadata_representatives.tsv"), col_names=TRUE)


# similarities=read.table("../../../motif-clustering-Viestra-private/metadata/motifsimilarity_to_DSA.tsv", header = TRUE, sep = "\t")
# 
# readLines(          "../../../motif-clustering-Viestra-private/metadata/motifsimilarity_to_DSA.tsv", n=5)              
# 
# read.tsv("../../../motif-clustering-Viestra-private/metadata/motifsimilarity.tsv", header = FALSE, sep = " ")
# 
# df_motifsimilarity=read_tsv(paste0('../../../motif-clustering-Viestra-private/','/metadata/motifsimilarity.tsv')) #5_423_571 rows


df_motifsimilarity=read_tsv(paste0('../../../motif-clustering-Viestra-private/','/metadata/sstat.tsv')) #5_423_571 rows

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
#there is an edge between the nodes (represented by 1); otherwise, there is not (represented by 0). The diagonal is set to 0 to remove self-loops?

#Convert this distance matrix into an adjacency matrix using a threshold:

threshold <- 1.5e-5  # replace with your threshold
adj_mat <- ifelse(sim_sym > threshold, 1, 0)
diag(adj_mat) <- 1  # #each query and target will at least be "close" to themselves in this hash. Or should this be 0
head(adj_mat)

g <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", diag = FALSE)

plot(g)

#Read in the representatives

node <- "your_node"  # replace with your node

#representatives <- as.list(read.table("../../../SELEX-Dominating-Set-Analysis/solutions/motifsimilarity_05_min_dom_set_list.txt", header = FALSE, sep = " ")$V1)
representatives <- as.list(read.table("../../../SELEX-Dominating-Set-Analysis/solutions/SELEX_all_min_dom_set_list.txt", header = FALSE, sep = " ")$V1)

represented_by_representatives=lapply(representatives, function(x, g) neighbors(g, x)$name, g )

names(represented_by_representatives)=representatives

# Convert the list to a tibble
df <- tibble::enframe(represented_by_representatives, name = "representative", value = "represented")

# Pad the list elements with NAs to make them of equal length
max_length <- max(map_int(df$represented, length)) #122

padded_list <- map(df$represented, ~ c(.x, rep(NA, max_length - length(.x))))
# Convert the list to a dataframe

tmp=do.call(rbind, padded_list) #1031 x 122

tmp2 <- cbind(data.frame(representative=df$representative), data.frame(tmp))

write.table(tmp2, file=paste0(TFBS_path, "PWMs_final/","sstat_represented_by_representatives.tsv"), row.names = FALSE, col.names = FALSE, sep="\t")

