
mypks=mypks[-which(mypks %in% c("stats", "graphics", "grDevices", "utils",  "datasets",  "methods",   "base", "compiler", 
                                "cli", "tools", "rstudioapi","pbdZMQ" )  )]



#BiocManager::install("RCy3")
library(RCy3)
library(igraph)

network=RCy3::importNetworkFromFile("../../ATAC-seq-peaks/CATLAS/edges.tsv.sif")


# SIF (Simple Interaction Format)
# NNF (Nested Network Format)
# XGMML
# GraphML
# PSI-MI Level 1 and 2.5
# Cytoscape.js JSON (can be used with Cytoscape.js, a tool described in future sections)
# CX JSON (for Cytoscape Cyberinfrastructure network exchange)
# The default format is SIF, which is a simple tab-delimited network format that provides node names and edge interactions only

library("jsonlite")
network <- fromJSON("../../ATAC-seq-peaks/CATLAS/edges.tsv.cyjs")
str(network$elements$nodes)

#network$elements$nodes$data 442 rows, id, name, Node_Type
nodes=cbind(network$elements$nodes$data, network$elements$nodes$position) 
edges=network$elements$edges$data
edges$distance=NA

#Compute edge distances

for(i in 1:nrow(edges)){

  x_source=nodes$x[ which(nodes$id==edges$source[i])]
  y_source=nodes$y[ which(nodes$id==edges$source[i])]

  x_target=nodes$x[ which(nodes$id==edges$target[i])]
  y_target=nodes$y[ which(nodes$id==edges$target[i])]

  sqrt((x_source - x_target)^2+(y_source - y_target)^2)

  edges$distance[i]=dist(t(matrix(c(x_source, y_source, x_target, y_target) ,nrow=2, ncol=2)), method="euclidean")
  
}

nodes$name[match(edges$source, nodes$id)]

nodes$name[match(edges$target, nodes$id)]

nodes$name[match(nodes$id, edges$source)]

nodes$name[match(nodes$id, edges$target)]

#relations <- data.frame(from=edges$source, to=edges$target, length=edges$distance)
relations <- data.frame(from=nodes$name[match(edges$source, nodes$id)], to=nodes$name[match(edges$target, nodes$id)], length=edges$distance)

#    to=,
 #   same.dept=,
  #  length=,
#)

nodes=nodes[, -1]

library("igraph")
network=graph_from_data_frame(relations, directed = TRUE, vertices = nodes)
print(network, e = TRUE, v = TRUE)
plot(network)


library("tidyverse")
library("tidygraph")
library("ggraph")
#install.packages(c("tidygraph", "ggraph"))

tbl_network = tbl_graph(nodes=nodes[,-c(7,8)], edges=relations, directed=TRUE)

tbl_network=as_tbl_graph(network)
str(tbl_network)

dend_layout <- tbl_network %>% 
  layout_tbl_graph_dendrogram()

pdf("test.pdf")
ggraph(tbl_network)
dev.off()
gg=ggraph(tbl_network, layout = "dendrogram") + 
  #geom_edge_link(width = 1, colour = "lightgray") +
  #geom_node_point(size = 4, colour = "#00AFBB") +
  #geom_node_text(aes(label = name), repel = TRUE)+
  theme_graph()

plot(gg)


library("navdata")
install.packages("navdata")

install.packages(c("tidygraph", "nycflights13"))
library(tidygraph)
install.packages("nycflights13")
library(nycflights13)

# Create tbl_graph object
g_tbl <- flights %>%
  filter(!is.na(dep_time), !is.na(arr_time)) %>%   # Filter out missing values
  select(origin, dest) %>%                         # Select only the origin and destination columns
  as_tbl_graph(directed=TRUE)                                   # Convert to tbl_graph object

# Print the tbl_graph object

ggraph(g_tbl, layout = "graphopt") + 
  geom_edge_link(width = 1, colour = "lightgray") +
  geom_node_point(size = 4, colour = "#00AFBB") +
  geom_node_text(aes(label = name), repel = TRUE)+
  theme_graph()



data("phone.call2")
phone.net <- tbl_graph(
  nodes = phone.call2$nodes, 
  edges = phone.call2$edges,
  directed = TRUE
)




dist_matrix <- as.dist(get.adjacency(network))
hc <- hclust(dist_matrix, method = "ward.D2")
plot(as.dendrogram(hc))

hc <- as.hclust(distances(network))

temp=as.dist(get.adjacency(network))

dend <- as.dendrogram(hclust(as.dist(get.adjacency(network)), method = "ward.D2"))

network <- read.table("../../ATAC-seq-peaks/CATLAS/edges.tsv.sif", header = FALSE, sep = "\t")

edges <- read.table("../../ATAC-seq-peaks/CATLAS/edges.csv", header = FALSE, sep = "\t")
edges <- read.table("../../ATAC-seq-peaks/CATLAS/edges.csv", header = FALSE, sep = "\t")
edges <- read.table("../../ATAC-seq-peaks/CATLAS/edges.csv", header = FALSE, sep = "\t")
g <- graph_from_data_frame(network, directed = FALSE)



plot(dend)

actors <- data.frame(
  name = c(
    "Alice", "Bob", "Cecil", "David",
    "Esmeralda"
  ),
  age = c(48, 33, 45, 34, 21),
  gender = c("F", "M", "F", "M", "F")
)
relations <- data.frame(
  from = c(
    "Bob", "Cecil", "Cecil", "David",
    "David", "Esmeralda"
  ),
  to = c("Alice", "Bob", "Alice", "Alice", "Bob", "Alice"),
  same.dept = c(FALSE, FALSE, TRUE, FALSE, FALSE, TRUE),
  friendship = c(4, 5, 5, 2, 1, 1), advice = c(4, 5, 5, 4, 2, 3)
)

g <- graph_from_data_frame(relations, directed = TRUE, vertices = actors)
print(g, e = TRUE, v = TRUE)