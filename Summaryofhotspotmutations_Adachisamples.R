library(data.table)
library(dplyr)
library(tidyr)
PDXmut = data.frame(fread("C:/Users/Raul/Documents/Proteomics_project/PDX_samples/all_samples_mutations_merged_with_sample.tsv"))
head(PDXmut)
hotspot = data.frame(fread("C:/Users/Raul/Documents/Proteomics_project/PDX_samples//hotspots.txt"))
PDXmut = PDXmut[grep("synonymous_variant" ,PDXmut$Consequence, invert = TRUE),]
PDXmut = PDXmut[grep("intron_variant" ,PDXmut$Consequence, invert = TRUE),]
PDXmut = PDXmut[grep("non_coding_transcript_exon_variant" ,PDXmut$Consequence, invert = TRUE),]
PDXmut = PDXmut[grep("non_coding_transcript_variant" ,PDXmut$Consequence, invert = TRUE),]
PDXmut = PDXmut[grep("3_prime_UTR_variant" ,PDXmut$Consequence, invert = TRUE),]
PDXmut = PDXmut[grep("5_prime_UTR_variant" ,PDXmut$Consequence, invert = TRUE),]
  

hotspot_split <-data.frame( hotspot %>%
  separate_rows(Variants, sep = "\\|") %>%
  mutate(Variant_Letter = sub(":.*", "", Variants)))

hotspot_split$ResidueAA = substring(hotspot_split$Residue,1,1)

hotspot_split$mutationslash = apply(cbind(hotspot_split$ResidueAA, "/", hotspot_split$Variant_Letter),1,paste, collapse = "" )
PDXmutGeneswithHotspot = PDXmut[PDXmut$SYMBOL %in% hotspot_split$Gene,]  


PDXmutGeneswithHotspot_atHotspot = data.frame()

for(cont in 1: dim(PDXmutGeneswithHotspot)[1]){
  thiscase = PDXmutGeneswithHotspot[cont,]
  hotspotsthisgene = hotspot_split[hotspot_split$Gene == thiscase$SYMBOL,]
  numeric_part <- gsub("[^0-9]", "", hotspotsthisgene$Residue)

  if(sum(thiscase$Protein_position == numeric_part & thiscase$Amino_acids == hotspotsthisgene$mutationslash) != 0){
    print(cont)
    PDXmutGeneswithHotspot_atHotspot = rbind(PDXmutGeneswithHotspot_atHotspot,thiscase )
  }
  
}
PDXmutGeneswithHotspot_atHotspotonlyTumor = PDXmutGeneswithHotspot_atHotspot[grep("TG3", PDXmutGeneswithHotspot_atHotspot$Sample, invert = T),]

PDXmutGeneswithHotspot_atHotspotonlyTumor = PDXmutGeneswithHotspot_atHotspot

PDXmutGeneswithHotspot_atHotspotonlyTumor$Sample <- sub("TG3$", "", PDXmutGeneswithHotspot_atHotspotonlyTumor$Sample)

samples <- c(
  "J-PDX1249", "J-PDX0919", "J-PDX0850", "J-PDX0845", "J-PDX0826",
  "J-PDX0769", "J-PDX0598", "J-PDX0596", "J-PDX0438", "J-PDX0381",
  "J-PDX0360", "J-PDX0224", "J-PDX0207", "J-PDX1107", "J-PDX1047",
  "J-PDX1037", "J-PDX0987", "J-PDX0961", "J-PDX0940", "J-PDX0927",
  "J-PDX0775", "J-PDX0747", "J-PDX0667", "J-PDX0591", "J-PDX0590",
  "J-PDX0587", "J-PDX0560", "J-PDX0454", "J-PDX0334", "J-PDX0326",
  "J-PDX0282", "J-PDX0277", "J-PDX0270", "J-PDX0263", "J-PDX0220",
  "J-PDX0175", "J-PDX0095", "J-PDX0093", "J-PDX0083", "J-PDX0079",
  "J-PDX0065", "J-PDX0013", "J-PDX0009", "J-PDX1002", "J-PDX0958",
  "J-PDX0865", "J-PDX0804", "J-PDX0464", "J-PDX0387"
)


PDXmutGeneswithHotspot_atHotspotonlyTumor[PDXmutGeneswithHotspot_atHotspotonlyTumor$Sample %in% samples,]



library(ggplot2)

# Convert table to data frame
mat <- table(PDXmutGeneswithHotspot_atHotspotonlyTumor$SYMBOL,
             PDXmutGeneswithHotspot_atHotspotonlyTumor$Sample)
df <- as.data.frame(mat)

# Convert counts to presence/absence
df$Present <- df$Freq > 0

# Plot
ggplot(df, aes(x = Var2, y = Var1, fill = Present)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "black")) +
  labs(x = "Sample", y = "Gene", fill = "Present") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


write.table(PDXmutGeneswithHotspot_atHotspotonlyTumor,"C:/Users/Raul/Dropbox/Papers/DIA-NN/Figures//PDXmutGeneswithHotspot_atHotspotonlyTumor.tsv", row.names = FALSE, sep = "\t", quote = FALSE)


library(ComplexHeatmap)
library(dplyr)

# Create the matrix: genes (rows) × samples (columns)
mat <- table(PDXmutGeneswithHotspot_atHotspotonlyTumor$SYMBOL,
             PDXmutGeneswithHotspot_atHotspotonlyTumor$Sample)

# Convert to character matrix with mutation presence as "Mutated"
binary_mat <- ifelse(mat > 0, "Mutated", "")

# Define color
col <- c("Mutated" = "black")

# Define the drawing function
alter_fun = list(
  background = function(x, y, w, h) {
    grid::grid.rect(x, y, w, h, gp = grid::gpar(fill = "#f0f0f0", col = NA))
  },
  Mutated = function(x, y, w, h) {
    grid::grid.rect(x, y, w*0.9, h*0.9, gp = grid::gpar(fill = "black", col = NA))
  }
)
pdf("C:/Users/Raul/Dropbox/Papers/DIA-NN/Figures/hotspotmutationsSummary_oncoprint.pdf", width = 15, height = 7.5)

# Create OncoPrint
oncoPrint(binary_mat,
          alter_fun = alter_fun,
          col = col,
          remove_empty_columns = TRUE,
          show_column_names = TRUE,
          show_row_names = TRUE,
          row_names_side = "left",   # Move gene names to the left
          pct_side = "right"   ,      # Move percentages to the right
          pct_gp = grid::gpar(fontsize = 10),  # adjust font if needed
          heatmap_legend_param = list(title = "Mutation"),
          column_title = "Hotspot Mutation Landscape in Tumoral samples")

dev.off()


library(ggplot2)

# Convert table to data frame
mat <- table(PDXmutGeneswithHotspot_atHotspotonlyTumor$SYMBOL,
             PDXmutGeneswithHotspot_atHotspotonlyTumor$Sample)
df <- as.data.frame(mat)

# Convert counts to presence/absence
df$Present <- df$Freq > 0

# Plot
ggplot(df, aes(x = Var2, y = Var1, fill = Present)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "black")) +
  labs(x = "Sample", y = "Gene", fill = "Present") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



library(ComplexHeatmap)
library(dplyr)

# Create the matrix: genes (rows) × samples (columns)
mat <- table(PDXmutGeneswithHotspot_atHotspotonlyTumor$SYMBOL,
             PDXmutGeneswithHotspot_atHotspotonlyTumor$Sample)

# Convert to character matrix with mutation presence as "Mutated"
binary_mat <- ifelse(mat > 0, "Mutated", "")

# Define color
col <- c("Mutated" = "black")

# Define the drawing function
alter_fun = list(
  background = function(x, y, w, h) {
    grid::grid.rect(x, y, w, h, gp = grid::gpar(fill = "#f0f0f0", col = NA))
  },
  Mutated = function(x, y, w, h) {
    grid::grid.rect(x, y, w*0.9, h*0.9, gp = grid::gpar(fill = "black", col = NA))
  }
)
pdf("C:/Users/Raul/Dropbox/Papers/DIA-NN/Figures/hotspotmutationsSummary_oncoprint.pdf", width = 15, height = 7.5)

# Create OncoPrint
oncoPrint(binary_mat,
          alter_fun = alter_fun,
          col = col,
          remove_empty_columns = TRUE,
          show_column_names = TRUE,
          show_row_names = TRUE,
          row_names_side = "left",   # Move gene names to the left
          pct_side = "right"   ,      # Move percentages to the right
          pct_gp = grid::gpar(fontsize = 10),  # adjust font if needed
          heatmap_legend_param = list(title = "Mutation"),
          column_title = "Hotspot Mutation Landscape in Tumoral samples")

dev.off()


PDXmutGeneswithHotspot_atHotspotonlyTumor=PDXmutGeneswithHotspot_atHotspotonlyTumor[PDXmutGeneswithHotspot_atHotspotonlyTumor$Sample %in% samples,]

library(ggplot2)

# Convert table to data frame
mat <- table(PDXmutGeneswithHotspot_atHotspotonlyTumor$SYMBOL,
             PDXmutGeneswithHotspot_atHotspotonlyTumor$Sample)
df <- as.data.frame(mat)

# Convert counts to presence/absence
df$Present <- df$Freq > 0

# Plot
ggplot(df, aes(x = Var2, y = Var1, fill = Present)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "black")) +
  labs(x = "Sample", y = "Gene", fill = "Present") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



library(ComplexHeatmap)
library(dplyr)

# Create the matrix: genes (rows) × samples (columns)
mat <- table(PDXmutGeneswithHotspot_atHotspotonlyTumor$SYMBOL,
             PDXmutGeneswithHotspot_atHotspotonlyTumor$Sample)

# Convert to character matrix with mutation presence as "Mutated"
binary_mat <- ifelse(mat > 0, "Mutated", "")

# Define color
col <- c("Mutated" = "black")

# Define the drawing function
alter_fun = list(
  background = function(x, y, w, h) {
    grid::grid.rect(x, y, w, h, gp = grid::gpar(fill = "#f0f0f0", col = NA))
  },
  Mutated = function(x, y, w, h) {
    grid::grid.rect(x, y, w*0.9, h*0.9, gp = grid::gpar(fill = "black", col = NA))
  }
)
pdf("C:/Users/Raul/Dropbox/Papers/DIA-NN/Figures/hotspotmutationsSummary_oncoprint_Adachisamples.pdf", width = 15, height = 7.5)

# Create OncoPrint
oncoPrint(binary_mat,
          alter_fun = alter_fun,
          col = col,
          remove_empty_columns = TRUE,
          show_column_names = TRUE,
          show_row_names = TRUE,
          row_names_side = "left",   # Move gene names to the left
          pct_side = "right"   ,      # Move percentages to the right
          pct_gp = grid::gpar(fontsize = 10),  # adjust font if needed
          heatmap_legend_param = list(title = "Mutation"),
          column_title = "Hotspot Mutation Landscape in Tumoral samples")

dev.off()

