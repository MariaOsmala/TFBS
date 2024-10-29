

# readLines("../../../SELEX-Dominating-Set-Analysis/data/Selex_all.dat", 
#          n=3)

# graph_data <- read.table("../../../SELEX-Dominating-Set-Analysis/data/Selex_all.dat", header = FALSE, sep = "\t")

library("readr")
library("tidyverse")
library("igraph")
library("readr")

TFBS_path="/Users/osmalama/projects/TFBS/"
representatives <- read_tsv(paste0(TFBS_path, "PWMs_final_version2.2/","new_spacings_representatives.tsv"), col_names=TRUE)

n=nrow(representatives)
n*(n-1)/2

df_motifsimilarity=read_tsv(paste0(TFBS_path, "PWMs_final_version2.2/", "sstat_new_spacings.tsv")) #639 015

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
representatives <- as.list(read.table("../../../SELEX-Dominating-Set-Analysis/solutions/SELEX_new_spacing_version2.2_min_dom_set_list.txt", header = FALSE, sep = " ")$V1)

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
table(sapply(cleaned_list, length))

removed_indices <- lapply(processed_list, `[[`, "removed_indices")

name=names(which(lapply(removed_indices, length)!=0)[1])

#removed_indices[name]

represented_by_representatives[[name]][removed_indices[[name]]] %in% unlist(representatives)

represented_by_representatives=cleaned_list



saveRDS(represented_by_representatives, file="RData/represented_by_representatives_version2.2_new_spacings.RDS")



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

write.table(tmp2, file=paste0(TFBS_path, "PWMs_final_version2.2/","sstat_represented_by_representatives_new_spacings.tsv"), 
            row.names = FALSE, col.names = FALSE, sep="\t")

hist(sapply(represented_by_representatives, length)+1)
