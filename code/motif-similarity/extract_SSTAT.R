library(dplyr)
library(readr)
library(purrr)
library(ggplot2)
library(tidyr)
library(gridExtra)

# This code extracts the SSTAT similarity scores into a single data frame,
# extracts also similarities with the motif itself
# Processes the SSTAT similarities into a format required by the dominating set analysis

results_path="/scratch/project_2006203/MOSTA-SSTAT/results_final_version2.2_correct"

# Reading the CSV file into df_motif_info data frame
#The number of motifs is 3933
df_motif_info <- read_tsv("/scratch/project_2006203/TFBS/PWMs_final_version2.2/metadata_3993_motifs.csv", col_types = cols()) 

# Initialize an empty data frame for df_sstat
df_sstat <- data.frame()

nrows=c()

# Loop over files and read data
for (i in 0:1000) {
  cat(i, "\n") # Print the iteration number
  
  sstat_path <- paste0(results_path, "/result_", i, ".out")
  sstat <- read_delim(sstat_path, "\t", col_names = c("Query_ID", "Target_ID", "Smax", "Ssum", "imax", "bimaxp", "alphaA", "alphaB"), col_types = cols())
  nrows=c(nrows, nrow(sstat))
  # Append the current data frame to df_sstat
  df_sstat <- bind_rows(df_sstat, sstat)
}

# Writing df_sstat to a TSV file, there is 7 732 278 rows, correct! 
n=3933
n*(n-1)/2 #7 732 278

# Read the distance between the same motif
results_path="/scratch/project_2006203/MOSTA-SSTAT/results_final_version2.2"
sstat <- read_delim(paste0(results_path, "/result_with_itself.out"), "\t", 
                    col_names = c("Query_ID", "Target_ID", "Smax", "Ssum", "imax", "bimaxp", "alphaA", "alphaB"), escape_double = FALSE, trim_ws = TRUE)

# Combine sstat with df_sstat and write to file
sstat_all <- bind_rows(sstat, df_sstat) %>%
  arrange(Query_ID)

write_delim(sstat_all, "/scratch/project_2006203/TFBS/PWMs_final_version2.2/sstat_to_DSA.tsv", "\t", col_names = FALSE) #7 736 211
write_tsv(df_sstat, "/scratch/project_2006203/TFBS/PWMs_final_version2.2/sstat.tsv") #7 732 278


hist(sstat$Ssum)

stop()
# Read the similarity data and draw Figures
df_sstat <- read_delim("/projappl/project_2006203/TFBS/PWMs_final_version2/sstat.tsv", "\t", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE)

# Calculating combinations
n <- nrow(df_motif_info)
k <- 2
P = choose(n, k)
print(P) #n*(n-1)/2=7_926_171

# Handling -Inf in Smax
df_sstat$Smax <- ifelse(is.infinite(df_sstat$Smax), min(df_sstat$Smax[!is.infinite(df_sstat$Smax)]), df_sstat$Smax)

min(df_sstat$Smax)
max(df_sstat$Smax)

min(df_sstat$Ssum)
max(df_sstat$Ssum)


# Histogram for Ssum
ggplot(df_sstat, aes(x = Ssum)) +
  geom_histogram(bins = 100, fill = "blue", color = "black") +
  scale_x_continuous(limits = c(0.00005, 0.0002), breaks = seq(0.00005, 0.0002, length.out = 100)) +
  labs(title = "SSTAT Ssum", x = "Similarity", y = "Frequency") +
  theme_minimal()

#ggsave("../MOSTA-SSTAT/figures/Ssum-hist.png")

# Histogram for Smax
ggplot(df_sstat, aes(x = Smax)) +
  geom_histogram(bins = 50, fill = "blue", color = "black") +
  scale_x_continuous(limits = c(-7, 11), breaks = seq(-7, 11, length.out = 50)) +
  labs(title = "SSTAT Smax", x = "Similarity", y = "Frequency") +
  theme_minimal()

#ggsave("../MOSTA-SSTAT/figures/Smax-hist.png")

