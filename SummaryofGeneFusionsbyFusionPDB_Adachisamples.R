library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

samples <- c(
  "J-PDX1249_T", "J-PDX0919_T", "J-PDX0850_T", "J-PDX0845_T", "J-PDX0826_T",
  "J-PDX0769_T", "J-PDX0598_T", "J-PDX0596_T", "J-PDX0438_T", "J-PDX0381_T",
  "J-PDX0360_T", "J-PDX0224_T", "J-PDX0207_T", "J-PDX1107_T", "J-PDX1047_T",
  "J-PDX1037_T", "J-PDX0987_T", "J-PDX0961_T", "J-PDX0940_T", "J-PDX0927_T",
  "J-PDX0775_T", "J-PDX0747_T", "J-PDX0667_T", "J-PDX0591_T", "J-PDX0590_T",
  "J-PDX0587_T", "J-PDX0560_T", "J-PDX0454_T", "J-PDX0334_T", "J-PDX0326_T",
  "J-PDX0282_T", "J-PDX0277_T", "J-PDX0270_T", "J-PDX0263_T", "J-PDX0220_T",
  "J-PDX0175_T", "J-PDX0095_T", "J-PDX0093_T", "J-PDX0083_T", "J-PDX0079_T",
  "J-PDX0065_T", "J-PDX0013_T", "J-PDX0009_T", "J-PDX1002_T", "J-PDX0958_T",
  "J-PDX0865_T", "J-PDX0804_T", "J-PDX0464_T", "J-PDX0387_T"
)


PDXgf = data.frame(fread("C:/Users/Raul/Documents/Proteomics_project/PDX_samples/Level3_Gene_fusions_in_GCATworkflow.tsv"))

genefusions = unique(data.frame(gene_fusion = apply(cbind(PDXgf$V1,"_",PDXgf$V2), 1, paste, collapse = ""), Sample =PDXgf$V3))

genefusions = genefusions [genefusions$Sample %in% samples,]
mat <- table(genefusions)
# Convert table to data frame

df <- as.data.frame(mat)

# Convert counts to presence/absence
df$Present <- df$Freq > 0



library(ComplexHeatmap)
library(dplyr)

# Create the matrix: genes (rows) × samples (columns)
mat <- table(genefusions)

# Convert to character matrix with mutation presence as "Mutated"
binary_mat <- ifelse(mat > 0, "Fusion", "")

# Define color
col <- c("Fusion" = "black")

# Define the drawing function
alter_fun = list(
  background = function(x, y, w, h) {
    grid::grid.rect(x, y, w, h, gp = grid::gpar(fill = "#f0f0f0", col = NA))
  },
  Fusion = function(x, y, w, h) {
    grid::grid.rect(x, y, w*0.9, h*0.9, gp = grid::gpar(fill = "black", col = NA))
  }
)
pdf("C:/Users/Raul/Documents/Proteomics_project/PDX_samples/GeneFusionsSummary_oncoprint_Level3_Adachisamples.pdf", width = 15, height = 15)

# Create OncoPrint
oncoPrint(binary_mat,
          alter_fun = alter_fun,
          col = col,
          remove_empty_columns = TRUE,
          show_column_names = TRUE,
          show_row_names = TRUE,
          row_names_side = "left",   # Move gene names to the left
          pct_side = "right"   ,      # Move percentages to the right
          pct_gp = grid::gpar(fontsize = 15),  # adjust font if needed
          column_title = "Gene Fusion Landscape in Tumoral samples Level 3 Adachi-sensei samples",
          column_title_gp = grid::gpar(fontsize = 20),  # Title font size
          column_names_gp = grid::gpar(fontsize = 15),  # X axis label font size
          heatmap_legend_param = list(title = "Fusion"))
dev.off()





library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
PDXgf2 = data.frame(fread("C:/Users/Raul/Documents/Proteomics_project/PDX_samples/Level2_Gene_fusions_in_GCATworkflow.tsv"))

PDXgf3 = data.frame(fread("C:/Users/Raul/Documents/Proteomics_project/PDX_samples/Level3_Gene_fusions_in_GCATworkflow.tsv"))


PDXgf = rbind(PDXgf2,PDXgf3)

genefusions = unique(data.frame(gene_fusion = apply(cbind(PDXgf$V1,"_",PDXgf$V2), 1, paste, collapse = ""), Sample =PDXgf$V3))

genefusions = genefusions [genefusions$Sample %in% samples,]



mat <- table(genefusions)
# Convert table to data frame

df <- as.data.frame(mat)

# Convert counts to presence/absence
df$Present <- df$Freq > 0



library(ComplexHeatmap)
library(dplyr)

# Create the matrix: genes (rows) × samples (columns)
mat <- table(genefusions)

# Convert to character matrix with mutation presence as "Mutated"
binary_mat <- ifelse(mat > 0, "Fusion", "")

# Define color
col <- c("Fusion" = "black")

# Define the drawing function
alter_fun = list(
  background = function(x, y, w, h) {
    grid::grid.rect(x, y, w, h, gp = grid::gpar(fill = "#f0f0f0", col = NA))
  },
  Fusion = function(x, y, w, h) {
    grid::grid.rect(x, y, w*0.9, h*0.9, gp = grid::gpar(fill = "black", col = NA))
  }
)
pdf("C:/Users/Raul/Documents/Proteomics_project/PDX_samples/GeneFusionsSummary_oncoprint_Level2_Adachisamples.pdf", width = 15, height = 15)

# Create OncoPrint
oncoPrint(binary_mat,
          alter_fun = alter_fun,
          col = col,
          remove_empty_columns = TRUE,
          show_column_names = TRUE,
          show_row_names = TRUE,
          row_names_side = "left",   # Move gene names to the left
          pct_side = "right"   ,      # Move percentages to the right
          pct_gp = grid::gpar(fontsize = 15),  # adjust font if needed
          column_title = "Gene Fusion Landscape in Tumoral samples Level 2 Adachi-sensei samples",
          column_title_gp = grid::gpar(fontsize = 20),  # Title font size
          column_names_gp = grid::gpar(fontsize = 15),  # X axis label font size
          heatmap_legend_param = list(title = "Fusion"))

dev.off()






PDXgf = data.frame(fread("C:/Users/Raul/Documents/Proteomics_project/PDX_samples/Level1_Gene_fusions_in_GCATworkflow.tsv"))

genefusions = unique(data.frame(gene_fusion = apply(cbind(PDXgf$V1,"_",PDXgf$V2), 1, paste, collapse = ""), Sample =PDXgf$V3))

genefusions = genefusions [genefusions$Sample %in% samples,]
mat <- table(genefusions)
# Convert table to data frame

df <- as.data.frame(mat)

# Convert counts to presence/absence
df$Present <- df$Freq > 0



library(ComplexHeatmap)
library(dplyr)

# Create the matrix: genes (rows) × samples (columns)
mat <- table(genefusions)

# Convert to character matrix with mutation presence as "Mutated"
binary_mat <- ifelse(mat > 0, "Fusion", "")

# Define color
col <- c("Fusion" = "black")

# Define the drawing function
alter_fun = list(
  background = function(x, y, w, h) {
    grid::grid.rect(x, y, w, h, gp = grid::gpar(fill = "#f0f0f0", col = NA))
  },
  Fusion = function(x, y, w, h) {
    grid::grid.rect(x, y, w*0.9, h*0.9, gp = grid::gpar(fill = "black", col = NA))
    }
  )
pdf("C:/Users/Raul/Documents/Proteomics_project/PDX_samples/GeneFusionsSummary_oncoprint_Level1_Adachisamples.pdf", width = 15, height = 15)

# Create OncoPrint
oncoPrint(binary_mat,
          alter_fun = alter_fun,
          col = col,
          remove_empty_columns = TRUE,
          show_column_names = TRUE,
          show_row_names = TRUE,
          row_names_side = "left",
          pct_side = "right",
          pct_gp = grid::gpar(fontsize = 15),
          column_title = "Gene Fusion Landscape in Tumoral samples Level 1 Adachi-sensei samples",
          column_title_gp = grid::gpar(fontsize = 20),  # Title font size
          column_names_gp = grid::gpar(fontsize = 15),  # X axis label font size
          heatmap_legend_param = list(title = "Fusion"))
dev.off()



