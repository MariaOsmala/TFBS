library("readr")
library("tidyverse")
library("igraph")
library("readr")
library("RCy3")

# install.packages(c("readr", "tidyverse", "igraph"))
# if(!"RCy3" %in% installed.packages()){
#   install.packages("BiocManager")
#   BiocManager::install("RCy3")
# }

TFBS_path="/Users/osmalama/projects/TFBS/"

# Motif metadata ----------------------------------------------------------

metadata <- read_tsv(paste0(TFBS_path,"Data/SELEX-motif-collection/metadata_final.tsv"), col_names=TRUE)

n=nrow(metadata)
n*(n-1)/2 #7732278

# Extract SSTAT similarities and create distance matrix ----------------------------------------------

# This has the old motif names

df_motifsimilarity=read_tsv(paste0(TFBS_path, "PWMs_final_version2.2/", "sstat.tsv")) #7732278 rows

sim_lt=pivot_wider(df_motifsimilarity, id_cols="Query_ID",
                   names_from="Target_ID", values_from="Ssum") # 3932 3933
sim = matrix(0,dim(metadata)[1], dim(metadata)[1]) #3933 3933
dim(sim)
diag(sim)=1
tril_indices=lower.tri(sim, diag = FALSE) #3933 3933
ltril_indices=lower.tri(matrix(0, nrow(sim_lt),nrow(sim_lt)), diag = TRUE) #3932 3932

sim[tril_indices]=sim_lt[,-1][ltril_indices]


sim_sym = sim + t(sim) - diag(diag(sim), nrow(sim), ncol(sim))

a=colnames(sim_lt)[2]
b=sim_lt[,"Query_ID"]$Query_ID

a=c(a,b)

rownames(sim_sym)=a
colnames(sim_sym)=a

colnames(sim_sym) %in% metadata$ID %>% table()
rownames(sim_sym) %in% metadata$ID %>% table()

match_ind=match(colnames(sim_sym), metadata$old_ID)
colnames(sim_sym)=metadata$ID[match_ind]

match_ind=match(rownames(sim_sym), metadata$old_ID)
rownames(sim_sym)=metadata$ID[match_ind]

#The SSTAT similarities with itself are missing from the distance matrix, but maybe they are not needed

# Convert this distance matrix into an adjacency matrix using a threshold --------

# If the distance is less than or equal to the threshold,  there is an edge between the nodes (represented by 1); 
# Otherwise, there is not (represented by 0). 
#The diagonal is set to 0 to remove self-loops?

threshold <- 1.5e-5  #
adj_mat <- ifelse(sim_sym > threshold, 1, 0)
diag(adj_mat) <- 0  # #each query and target will at least be "close" to themselves in this hash. Or should this be 0
adj_mat[1:5,1:5]

g <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", diag = FALSE) # 27481 connections

edges=sapply(as_ids(E(g)), strsplit, "\\|")
names(edges)=NULL

edges=do.call(rbind, lapply(edges, function(x) data.frame(V1=x[1],V2=x[2]) ))

unique(c(edges$V1, edges$V2)) %>% length() #3198

# Find isolated nodes
isolated_nodes <- V(g)$name[degree(g) == 0] #735

# 735+3198 = 3933

# Write network in sif format ---------------------------------------------

#nodeA <relationship type> nodeB
#nodeC <relationship type> nodeA
#nodeD <relationship type> nodeE nodeF nodeB
#nodeG
#...
#nodeY <relationship type> nodeZ

names(edges)=c("Source", "Target")

edges$interaction_type="motif-motif"

#order columns 
edges <- edges %>% select(Source, interaction_type, Target)

write.table(edges, file="/Users/osmalama/projects/cytoscape/review/motif-network.sif", 
            quote=FALSE, row.names = FALSE, col.names=FALSE, sep="\t")

#write isolated nodes
write.table(data.frame(isolated_nodes), file="/Users/osmalama/projects/cytoscape/review/motif-network.sif", 
            quote=FALSE,append=TRUE, row.names = FALSE, col.names=FALSE, sep="\t")

# Write network in ncol format --------------------------------------------

write_graph(g, file="/Users/osmalama/projects/cytoscape/review/motif-network.ncol", format = c("NCOL"))

# Now assign edge weights (similarities) from similarity matrix
E(g)$weight <- apply(as_data_frame(g, what = "edges")[, 1:2], 1, function(x) {
  sim_sym[x[1], x[2]]
})

write_graph(g, file="/Users/osmalama/projects/cytoscape/review/motif-network-SSTAT-similarity.ncol", format = c("NCOL"))

write.table(isolated_nodes, file="/Users/osmalama/projects/cytoscape/review/isolated_nodes.ncol", 
            quote=FALSE, row.names = FALSE,
            col.names=FALSE)


# Minimum dominating sets -------------------------------------------------

#Don't use this as these have the old IDS
#representatives <- as.list(
#  read.table("~/projects/SELEX-Dominating-Set-Analysis/solutions/SELEX_all_version2.2_min_dom_set_list.txt", 
#             header = FALSE, sep = " ")$V1)

representatives = as.list(metadata %>% filter(representative=="YES") %>% pull(ID))

# Initialize list to store results
dom_neighborhoods <- list()

for (node in unlist(representatives)) {
  # Get neighbors of the node
  #node=unlist(representatives)
  neighbors_node <- neighbors(g, node, mode = "all")  # or "out"/"in" for directed
  
  # Convert to names (if needed)
  neighbor_names <- V(g)[neighbors_node]$name
  
  # Exclude other dominating set members
  non_dom_neighbors <- setdiff(neighbor_names, representatives)
  
  # Store result
  dom_neighborhoods[[node]] <- non_dom_neighbors
}

dom_sizes=sapply(dom_neighborhoods, length) 
names(dom_sizes)=NULL
table(dom_sizes==0) #TRUE 735


sapply(dom_neighborhoods, length) %>% sum() #3652
# 3652 +1232=4884, this is because non-representative can belong to multiple dom sets
c(names(dom_neighborhoods), unlist(dom_neighborhoods)) %>% unique() %>% length() # 3933

## My way ------------------------------------------------------------------

represented_by_representatives=lapply(representatives, function(x, g) neighbors(g, x)$name, g )
names(represented_by_representatives)=representatives

sapply(represented_by_representatives, length) %>% sum() #3730

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
sapply(represented_by_representatives, length) %>% sum() #3730

table(sapply(cleaned_list, length))
sapply(cleaned_list, length) %>% sum() #3652

represented_by_representatives=cleaned_list
c(names(represented_by_representatives), unlist(represented_by_representatives)) %>% unique() %>% length() 

saveRDS(represented_by_representatives, file=paste0(TFBS_path,"RProjects/TFBS/RData/represented_by_representatives_version2.2.RDS"))

# Separate networks for monomers and CAP-SELEX motifs ---------------------



# What about Methyl-SELEX motifs

monomers <- metadata %>% filter(experiment != "CAP-SELEX") %>% pull(ID) #2035 16883 edges
heterodimers <- metadata %>% filter(experiment == "CAP-SELEX") %>% pull(ID) # 1898

isolated_monomers = isolated_nodes[isolated_nodes %in% monomers]  #173
isolated_heterodimers=isolated_nodes[isolated_nodes %in% heterodimers] #280


## Monomer network ---------------------------------------------------------

### motif-motif network ---------------------------------------------------------

g_monomers <- subgraph(g, V(g)$name %in% monomers) #2023

edges_monomers=sapply(as_ids(E(g_monomers)), strsplit, "\\|")
names(edges_monomers)=NULL
edges_monomers=do.call(rbind, lapply(edges_monomers, function(x) data.frame(V1=x[1],V2=x[2]) ))

unique(c(edges_monomers$V1, edges_monomers$V2)) %>% length() #1857
#Are the missing ones the isolated monomers

#1857 + 173 = 2030

monomers[which(!(monomers %in% c(edges_monomers$V1, edges_monomers$V2, isolated_monomers)))]

#Are these representing heterodimers, seems so

metadata %>% filter(ID %in% monomers[which(!(monomers %in% c(edges_monomers$V1, edges_monomers$V2, isolated_monomers)))]) %>% select(ID,representative, seed, study)

#Add them as isolated monomers

isolated_monomers=c(isolated_monomers, 
                    metadata %>% filter(ID %in% monomers[which(!(monomers %in% c(edges_monomers$V1, edges_monomers$V2, isolated_monomers)))]) %>% pull(ID)
                    )

#Are there additional isolated monomers? No, all are already in the list
add_isolated_monomers <- V(g_monomers)$name[degree(g_monomers) == 0] #735

add_isolated_monomers %in% isolated_monomers %>% table() #All true

### motif-TF network ---------------------------------------------------------

edges_monomer_TF=metadata %>% filter(ID %in% monomers) %>% select(ID,symbol)
names(edges_monomer_TF)=c("Source", "Target")

edges_monomer_TF$interaction_type="motif-TF"

edges_monomer_TF <- edges_monomer_TF %>% select(Source, interaction_type, Target)

### Write network in sif format -----------------------------------------------

names(edges_monomers)=c("Source", "Target")

edges_monomers$interaction_type="motif-motif"

#order columns 
edges_monomers <- edges_monomers %>% select(Source, interaction_type, Target)

write.table(edges_monomers, file="/Users/osmalama/projects/cytoscape/review/monomer-network.sif", 
            quote=FALSE, append=FALSE, row.names = FALSE, col.names=FALSE, sep="\t")

#write isolated nodes
write.table(data.frame(isolated_monomers), file="/Users/osmalama/projects/cytoscape/review/monomer-network.sif", 
            quote=FALSE,append=TRUE, row.names = FALSE, col.names=FALSE, sep="\t")

#write motif-TF edges

write.table(edges_monomer_TF, file="/Users/osmalama/projects/cytoscape/review/monomer-network.sif", 
            quote=FALSE,append=TRUE, row.names = FALSE, col.names=FALSE, sep="\t")

#Write node list


monomer_metadata=metadata %>% filter(ID %in% monomers)

#Add new column, use dplyr 

monomer_metadata <- monomer_metadata %>%
  mutate(node_type = "motif") %>% relocate(node_type, .after = ID)




monomer_TF_metadata <- as.data.frame(matrix(NA, nrow = edges_monomer_TF$Target %>% unique()%>% length(), 
                                            ncol = ncol(monomer_metadata)))

colnames(monomer_TF_metadata) <- colnames(monomer_metadata)

monomer_TF_metadata$ID=edges_monomer_TF$Target %>% unique()
monomer_TF_metadata$node_type="TF"
monomer_metadata <- rbind(monomer_metadata, monomer_TF_metadata)

#monomer_metadata$representative_circles <- ifelse(monomer_metadata$representative == "YES", "YES", "NO")


#monomer_metadata <- mutate(monomer_metadata,
#                           TF_label = ifelse(node_type == "TF", ID, ""))


node_names$organism[which(node_names$organism=="NA")]="Homo_sapiens"
mouse_nodes= node_names %>% filter(organism=="Mus_musculus")

#Need to assign info about the protein family to the TF nodes itself

TF_nodes=node_names %>% filter(node_type=="TF")

TF_nodes$name %>% unique() %>% length() #718, all unique

#Map these to motifs

TF_motifs_list=lapply(TF_nodes$name, function(x) which(node_names$symbol==x))
sapply(TF_motifs_list, length) %>% table()

TF_protein_families_list=lapply(TF_motifs_list, function(x) unique(node_names$Lambert2018_families[x]) )
sapply(TF_protein_families_list, length) %>% table() #unique
TF_protein_families<- unlist(TF_protein_families_list)


write.table(monomer_metadata, file="/Users/osmalama/projects/cytoscape/review/monomer_node_table.tsv", 
            quote=FALSE, row.names = FALSE, col.names=TRUE, sep="\t")

#save.image("~/projects/TFBS/RData/motif_networks.RData")
load("~/projects/TFBS/RData/motif_networks.RData")


# Interact with cytoscape -------------------------------------------------

cytoscapePing()
cytoscapeVersionInfo()

saveSession('/Users/osmalama/projects/cytoscape/sessions/monomer_network') #.cys

#NODES should contain a column of character strings named: id. 
#This name can be overridden by the arg: node.id.list. Additional columns are loaded as node attributes.

#EDGES should contain columns of character strings named: source, target and interaction. 

#These names can be overridden by args: source.id.list, target.id.list, interaction.type.list. 

#Additional columns are loaded as edge attributes. 

#The 'interaction' list can contain a single value to apply to all rows; and if excluded altogether,
#the interaction type wiil be set to "interacts with". 

#NOTE: attribute values of types (num) will be imported as (Double); 
#(int) as (Integer); 
#(chr) as (String); 
#and (logical) as (Boolean). 
#(Lists) will be imported as (Lists) in CyREST v3.9+.

# tmp=rbind(edges_monomers,
#       data.frame(Source=isolated_monomers, interaction_type=NA, Target=NA),
#       edges_monomer_TF)
# 
# 
# createNetworkFromDataFrames(monomer_metadata,
#                             tmp, 
#                             title="Monomer Network", 
#                             collection="Review", 
#                             node.id.list="ID", 
#                             source.id.list="Source", 
#                             target.id.list="Target", 
#                             interaction.type.list="interaction_type")





# Open the network in cytospaces
# Interact programmatically with the network 

#clearVisualPropertyMappings(style.name = getCurrentStyle(), visual.property = "NODE_LABEL")

# a discrete mapping, where a style is defined for each, discrete value. 
# This is useful for categorical data (like type) where there is only 
#a limited set of possible values. 

#This is in contast to the other two other types of mappings: continuous and passthrough. 
# In the case of expression values, for example, we will want to use continuous mapping 
# (e.g., to node color), defining a small set of control points, 
# rather than an explicit color for each possible data value.


## Set the TF node shape to diamond and mouse motif to triangle ----------------------------------------


getNodeShapes() 
column <- 'node_type'
values <- c ('TF',  'motif')
shapes <- c ('DIAMOND', 'ELLIPSE')
setNodeShapeMapping (column, values, shapes)

lockNodeDimensions(TRUE)

node_names <- getTableColumns()



setNodeShapeBypass(node.names= mouse_nodes$name,
                   new.shapes="TRIANGLE")


## Remove node label from motifs -------------------------------------------



motif_nodes <- node_names %>% filter(node_type=="motif")

#This worked
setNodePropertyBypass(node.names= motif_nodes$name,
                      visual.property="NODE_LABEL",
                      new.value="")



## Representative motifs with highlighted shape boundary -------------------


#Get the current node edge color

getNodeColor() %>% unique() "#89D0F5"
setNodeColorDefault("#CCCCCC")

getNodeProperty(node_names$name,visual.property="NODE_BORDER_PAINT") %>% unique() #"#CCCCCC"

getNodeProperty(node_names$name,visual.property="NODE_BORDER_STROKE") %>% unique() #"SOLID"
getNodeProperty(node_names$name,visual.property = "NODE_BORDER_WIDTH" ) %>% unique() #0

getNodeProperty(node_names$name,visual.property = "NODE_BORDER_WIDTH" ) %>% unique() #0

setNodeBorderWidthDefault(3)

#Human
representatives=node_names %>% filter(node_type=="motif" & representative=="YES" 
                                      & organism=="Homo_sapiens")
setNodePropertyBypass(node.names= representatives$name,
                      visual.property="NODE_BORDER_PAINT",
                      new.value="blue")

representatives=node_names %>% filter(node_type=="motif" & representative=="YES" 
                                      & organism=="Mus_musculus")
setNodePropertyBypass(node.names= representatives$name,
                      visual.property="NODE_BORDER_PAINT",
                      new.value="red")


# Color the nodes according to protein families ------------------------






Protein_family_colours2=c()
Protein_family_colours2["Homeodomain"]="#9c4c5b"
Protein_family_colours2["Homeodomain; POU"]= "#9c4c5b"
Protein_family_colours2["Homeodomain; Paired box"]= "#9c4c5b"
Protein_family_colours2["CUT; Homeodomain"]= "#9c4c5b"    
Protein_family_colours2["Ets"]=  "#61c350"
Protein_family_colours2["Ets; AT hook"]=  "#61c350"
Protein_family_colours2["bHLH"]= "#aa53cc" 
Protein_family_colours2["T-box"]=  "#9cb835"                
Protein_family_colours2["bZIP"]="#696add"
Protein_family_colours2["C2H2 ZF"]= "#ceae33"
Protein_family_colours2["C2H2 ZF; AT hook"]= "#ceae33"
Protein_family_colours2["BED ZF" ]="#ceae33"
Protein_family_colours2["CCCH ZF"]= "#ceae33"
Protein_family_colours2["Znf"]= "#ceae33"  
Protein_family_colours2["Forkhead"]="#7752a2"
Protein_family_colours2["Nuclear receptor"]="#409437"
Protein_family_colours2["TEA"]="#d94aa9"
Protein_family_colours2["HMG/Sox"]="#4ec584"
Protein_family_colours2["GCM"]="#e8437e"            
Protein_family_colours2["Paired box"]="#54c9b3"
Protein_family_colours2["IRF"]="#db3750"         
Protein_family_colours2["E2F"]="#3bb7cb"
Protein_family_colours2["RFX"]="#cf402b"
Protein_family_colours2["Rel"]="#72a0da"
Protein_family_colours2["AP-2"]="#dd6629"                 
Protein_family_colours2["GATA"]="#6278c4"
Protein_family_colours2["HSF"]="#d98c2f"
Protein_family_colours2["Myb/SANT"]="#c885df"
Protein_family_colours2["Runt"]="#678728"                 
Protein_family_colours2["MADS box"]="#ab4a92"
Protein_family_colours2["Grainyhead"]="#38793f"
Protein_family_colours2["HMG"]="#b93a68"
Protein_family_colours2["SMAD"]="#5aa77b"                 
Protein_family_colours2["DM"]="#b64444"
Protein_family_colours2["Prospero"]="#2c7c62"
Protein_family_colours2["SAND"]="#e77f60"
Protein_family_colours2["MEIS"]="#94b269"
Protein_family_colours2["p53"]="#965581"
Protein_family_colours2["POU"]="#bcad5e"                 
Protein_family_colours2["TFAP"]="#de88b8"
Protein_family_colours2["CSD"]="#646d2c"                
Protein_family_colours2["CENPB"]="#dd7d7f"
Protein_family_colours2["EBF1"]="#8d751f"
Protein_family_colours2["NRF"]="#a15023"
Protein_family_colours2["RRM"]="#dea66d"                  
Protein_family_colours2["XPA"]="#9c6e3d"       

setNodeColorMapping('Lambert2018_families', 
                    colors=Protein_family_colours2, 
                    mapping.type='d')





#Change the node edge color based on representative column

setNodeColorMapping(
  column = "representative",
  values = c("YES", "NO"),
  colors = c("#FF0000", "#000000"),
  default.color = "#000000"
)


clearNodePropertyBypass(
  #node.names= motif_nodes$name ,
  node.names = as.character(node_names %>% filter(node_type=="motif") %>% select(name) %>% unlist()),
  visual.property="NODE_LABEL"
)

motif_node_names <- node_names %>% filter(node_type=="motif") %>% pull(name)

for(node in motif_node_names) {
  clearNodePropertyBypass(
    node.names = node,
    visual.property = "NODE_LABEL"
  )
}




node_names <- getTableColumns("node", "name")$name
clearNodePropertyBypass(node.names = node_names)

node_data <- getTableColumns("node", c("name", "node_type"))

label_nodes <- node_data
label_nodes$new_name=label_nodes$name
label_nodes$new_name[label_nodes$node_type!="TF"]=NA


setNodeLabelBypass(node.names = label_nodes$name,
                   new.labels = label_nodes$new_name)

# Heterodimer network -----------------------------------------------------

g_heterodimers <- subgraph(g, V(g)$name %in% heterodimers) #1898

edges_heterodimers=sapply(as_ids(E(g_heterodimers)), strsplit, "\\|")
names(edges_heterodimers)=NULL
edges_heterodimers=do.call(rbind, lapply(edges_heterodimers, 
                                     function(x) data.frame(V1=x[1],V2=x[2]) ))

unique(c(edges_heterodimers$V1, edges_heterodimers$V2)) %>% length() #1307
#Are the missing ones the isolated heterodimers, NO

#1307 + 559 = 1866

heterodimers[which(!(heterodimers %in% c(edges_heterodimers$V1, edges_heterodimers$V2, 
                                         isolated_heterodimers)))] #32

#Are these represented by monomers, seems so

metadata %>% filter(ID %in% heterodimers[
  which(!(heterodimers %in% c(edges_heterodimers$V1, edges_heterodimers$V2, isolated_heterodimers)))]
  ) %>% select(ID,representative, seed, study)

#Add them as isolated monomers

isolated_heterodimers=c(isolated_heterodimers, 
                    metadata %>% filter(ID %in% heterodimers[
                      which(!(heterodimers %in% c(edges_heterodimers$V1, edges_heterodimers$V2, isolated_heterodimers)))
                      ]) %>% pull(ID)
)

#Are there additional isolated monomers? No, all are already in the list
add_isolated_heterodimers <- V(g_heterodimers)$name[degree(g_heterodimers) == 0] #735

add_isolated_heterodimers %in% isolated_heterodimers %>% table() #All true


### Write network in sif format ---------------------------------------------

names(edges_heterodimers)=c("Source", "Target")

edges_heterodimers$interaction_type="motif-motif"

#order columns 
edges_heterodimers <- edges_heterodimers%>% select(Source, interaction_type, Target)

write.table(edges_heterodimers, file="/Users/osmalama/projects/cytoscape/review/heterodimer-network.sif", 
            quote=FALSE, row.names = FALSE, col.names=FALSE, sep="\t")

#write isolated nodes
write.table(data.frame(isolated_heterodimers), file="/Users/osmalama/projects/cytoscape/review/heterodimer-network.sif", 
            quote=FALSE,append=TRUE, row.names = FALSE, col.names=FALSE, sep="\t")


# Other stuff -------------------------------------------------------------------



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

