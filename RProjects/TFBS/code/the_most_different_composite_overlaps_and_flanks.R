library(tidyverse)
library(readr)
library(gridExtra)
library(ggplot2)
library(grid)
library(magick)


JSD_JI=read_tsv("../../PWMs_final_version2.2/Jensen-Shannon-Jaccard-Index-data.tsv") #977

#Those with composite overlap of length 3-6
JSD_JI_all=JSD_JI

#JSD_JI=JSD_JI %>% filter(overlap %in% c(3,4,5,6)) #643
JSD_JI=JSD_JI %>% filter(overlap %in% c(2,3)) #440


# Define a threshold for "high" values (e.g., 75th percentile)
high_threshold <- quantile(c(JSD_JI$sj_divergence_flank_TF1_core,
                             JSD_JI$sj_divergence_flank_TF2_core), 0.70)  # Example threshold


# Add a composite score to the data frame to prioritize high values and proximity to diagonal
JSD_JI_filtered <- JSD_JI %>%
  filter(JSD_JI$sj_divergence_flank_TF1_core > high_threshold & JSD_JI$sj_divergence_flank_TF2_core > high_threshold) %>%  # Only keep rows where both are high
  mutate(
    js_mean = (sj_divergence_flank_TF1_core + sj_divergence_flank_TF2_core)/2,       # Composite score for high values
    js_diff = abs(sj_divergence_flank_TF1_core - sj_divergence_flank_TF2_core)   # Measure of how close they are to each other (diagonal closeness)
  ) %>%
  arrange(desc(js_mean), js_diff)  # Sort by high value sum and then by closeness to diagonal



JSD_JI_filtered$label=toupper(paste0( strsplit(JSD_JI_filtered$first.monomer, "_") %>% map_chr(1), ":", strsplit(JSD_JI_filtered$second.monomer, "_") %>% map_chr(1)))

# Plot the scatter plot and highlight the sorted points

#label the points

library(ggplot2)
library(ggrepel)
ggplot(JSD_JI_all, aes(x = sj_divergence_flank_TF1_core, y = sj_divergence_flank_TF2_core)) +
  geom_point(color = "grey") +
  geom_point(data = JSD_JI_filtered, aes(x = sj_divergence_flank_TF1_core, y = sj_divergence_flank_TF2_core), color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +  # Add diagonal line
  theme_minimal() +
  labs(title = "Jensen-Shannon divergence",
       x = "between TF1 flank and composite overlap",
       y = "between TF2 flank and composite overlap") + 
  geom_text_repel(data = JSD_JI_filtered, aes(label = label), 
                          color = "black", 
                             size = 3, 
                           box.padding = 0.5,  # Padding around the label to avoid overlap
                           #point.padding = 1,  # Padding around the point to keep the label away from it
                           segment.color = "black", 
                  max.overlaps = 100)   # Color of the leader line



# Select TFs, overlap

data=JSD_JI_filtered[,c("composite","first.monomer", "second.monomer",             
"TF1_start",
"TF1_end",
"TF1_orientation",
"TF2_start",
"TF2_end",
"TF2_orientation",             
"overlap",
"sj_divergence_flank_TF1_core",
"sj_divergence_flank_TF2_core",
"JI_core_TF1_flank",
"JI_core_TF2_flank",
"js_mean",
"js_diff")]    


data$symbols=paste0( strsplit(data$first.monomer, "_") %>% map_chr(1), ":", strsplit(data$second.monomer, "_") %>% map_chr(1))

alignment_figure_path="/Users/osmalama/spacek/Figures_version2.2_composites_pdf/"

data$Image <-paste0(alignment_figure_path, 
                     paste0(JSD_JI_filtered$composite, ".pdf")
                     ) 

#Motifs for which the spacek alignment is different

data$spacek_image=""
ind=which(data$composite %in% c("ELF2_FOXK1_TTATGG40NCCA_YXI_NAGAAAACCGAAWMN_m2_c3b0u",
                      "HOXD9_TBX20_TACTAA40NTGTT_YWIIII_NAGGTGTTAATKN_m1_c3b0"))
data$spacek_image[ind]=paste0(alignment_figure_path, data$composite[ind], "_spacek.pdf")

data$spacek_image=paste0(alignment_figure_path, data$composite, ".pdf")
#data$spacek_image=paste0(alignment_figure_path, data$composite, "_spacek.pdf")
                                                                                                

data=data[,c("symbols",
       "TF1_start", #
        "TF1_end", #
       "TF1_orientation", #
        "TF2_start", #
        "TF2_end", #
        "TF2_orientation", #
        "sj_divergence_flank_TF1_core",
        "sj_divergence_flank_TF2_core",
        "JI_core_TF1_flank",            
        "JI_core_TF2_flank",
        "js_mean" ,
        "js_diff" ,
        "overlap",
        "Image","spacek_image")]

#remame columns 

colnames(data)=c("Symbols",
      "TF1\nstart", #
      "TF1\nend", #
      "TF1\norient", #
      "TF2\nstart", #
      "TF2\nend", #
      "TF2\norient", #
        "TF1\nJSD",
        "TF2\nJSD",
        "TF1\nJI",
        "TF2\nJI",
        "JSD\nmean",
        "JSD\ndiff",
        "Overlap",
        "Tomtom\nalignment",
        "Spacek\nalignment")

#present numerical values with certain precision


data$`TF1\nJSD`=round(data$`TF1\nJSD`,2)
data$`TF2\nJSD`=round(data$`TF2\nJSD`,2)
data$`TF1\nJI`=round(data$`TF1\nJI`,2)
data$`TF2\nJI`=round(data$`TF2\nJI`,2)
data$`JSD\nmean`=round(data$`JSD\nmean`,2)
data$`JSD\ndiff`=round(data$`JSD\ndiff`,2)




library(flextable)


library(magick)

# Create a blank image with optional text
create_placeholder_image <- function(width = 100, height = 100, text = "No Image", file_path = "placeholder.png") {
  # Create a blank white image
  img <- image_blank(width = width, height = height, color = "white")
  
  # Add optional text
  if (!is.null(text) && text != "") {
    img <- image_annotate(img, text, size = 20, color = "black", gravity = "center")
  }
  
  # Save the image to a file
  image_write(img, path = file_path)
  return(file_path)
}

# Generate a placeholder image
placeholder_path <- create_placeholder_image(width = 100, height = 100, text = "Same", file_path = "placeholder.png")


# Function to handle missing images
colformat_image_safe <- function(flextable_obj, j, width = 0.5, height = 0.5, placeholder_image = "placeholder.png") {
  # Replace missing or invalid file paths with the placeholder
  image_paths <- flextable_obj$x$dataset[[j]]
  valid_images <- sapply(image_paths, function(path) {
    if (is.character(path) && file.exists(path)) {
      TRUE
    } else {
      FALSE
    }
  })
  
  # Replace invalid paths with NA or a placeholder
  flextable_obj$x$dataset[[j]][!valid_images] <- placeholder_image
  
  # Apply colformat_image only for valid image paths
  flextable_obj <- flextable_obj %>%
    colformat_image(j = j, width = width, height = height)
  
  return(flextable_obj)
}

# Example usage with data

#How to change the font size 

data$Symbols=toupper(data$Symbols)
data[data==""]="placeholder.png"

ft=colformat_image_safe(ft)

col_widths=dim_pretty(ft)$widths # 6.879489
# 
col_widths=col_widths/sum(col_widths)
# 
# sum(c(0.13, 0.07, 0.07,0.07, 0.07,0.08, 0.07, 0.10,0.17, 0.17))
# 
# col_widths=c(0.19, 0.06, 0.06,0.06, 0.06,0.07, 0.06, 0.1,0.17, 0.17)*sum(dim_pretty(ft)$widths)
# 
# sum((col_widths[1:7]/5)) #1.235162
# 
# col_widths[2:7]=col_widths[2:7]-0.025
# 
# col_widths[8:9]=col_widths[8:9]+7*0.025/2

# Specify which columns are image paths
image_columns <- c("Tomtom\nalignment", "Spacek\nalignment")

# Calculate column widths
col_widths <- calculate_column_widths(data, image_cols = image_columns)


ft <- flextable(data) %>% fontsize(size=6, part="all") %>% 
  colformat_image(j = "Tomtom\nalignment", width = .85, height = 0.85) %>%
  colformat_image(j = "Spacek\nalignment", width = 0.85, height = 0.85) %>% align(align = "center", part = "all") 


for (col in seq_along(col_widths)) {
  ft <- ft %>% width(j = col, width = col_widths[col])
}

ft <- autofit(ft, part="all", add_w = 0, add_h = 0)



# Save the flextable to a Word document or PDF
library(officer)

# Save to a Word file
save_as_docx(ft, path = "table_with_images_new.docx")

# Optional: Save to a PDF (requires pandoc)
library("pandoc")
pandoc_convert(text=ft,to="pdf", file = "table_with_images_new.pdf")

# Save the flextable to an R Markdown file
save_as_rmarkdown <- function(ft, file = "table_output.Rmd") {
  rmd_content <- c(
    "---",
    "title: \"Table with Images\"",
    "output: pdf_document",
    "---",
    "",
    "```{r, echo=FALSE}",
    "library(flextable)",
    "ft <- ", deparse(substitute(ft)), # Embed the flextable directly
    "ft",
    "```"
  )
  
  writeLines(rmd_content, con = file)
}

rmd_file <- "table_output.Rmd"
save_as_rmarkdown(ft, file = rmd_file)
library(rmarkdown)
# Render the R Markdown file as a PDF
render(rmd_file, output_format = "pdf_document")


data=data[,c("symbols",
        "TF1_start",
        "TF1_end",
        "TF1_orientation",
        "TF2_start",
        "TF2_end",
        "TF2_orientation",
        "overlap",
        "Image", "spacek_image")]                      


# Function to scale and embed image in a cell
create_image_grob <- function(image_path, width = 0.8, height = 0.8) {
  img <- image_read(image_path) # Read the image
  raster <- rasterGrob(as.raster(img), width = unit(width, "npc"), height = unit(height, "npc"))
  return(raster)
}


# Create a list of rows with both text and images
table_rows <- mapply(
  function(text, image_path) {
    list(textGrob(text, x = 0, hjust = 0, gp = gpar(fontsize = 10)),
         create_image_grob(image_path))
  },
  text_data, image_paths,
  SIMPLIFY = FALSE
)

# Generate the table layout with borders
table_with_borders <- function(text_data, image_paths) {
  n_rows <- length(text_data)
  n_cols <- 2 # Fixed for this example (text + image)
  
  # Create an empty grid layout
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(n_rows, n_cols)))
  
  for (row in seq_along(text_data)) {
    # Add text cell
    # Draw borders for each cell
    grid.rect(gp = gpar(lwd = 2), vp = viewport(layout.pos.row = row, layout.pos.col = 1))
    grid.rect(gp = gpar(lwd = 2), vp = viewport(layout.pos.row = row, layout.pos.col = 2))
    grid.text(text_data[row], x = unit(0.5, "npc"), 
              y = unit((n_rows - row + 0.5) / n_rows, "npc"), 
              just = "center",
              vp = viewport(layout.pos.row = row, layout.pos.col = 1))
    
    # Add image cell
    img_grob <- create_image_grob(image_paths[row])
    grid.draw(editGrob(img_grob, vp = viewport(layout.pos.row = row, layout.pos.col = 2)))
    
    
  }
}

# Create the PDF
pdf("/Users/osmalama/table_with_borders.pdf", width = 8, height = 6)
table_with_borders(text_data, image_paths)
dev.off()
