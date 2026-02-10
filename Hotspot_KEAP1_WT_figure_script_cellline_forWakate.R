# Install required packages if not already installed
required_packages <- c(
  "stringr", "ggplot2", "RColorBrewer", "reshape2", "dplyr", "data.table",
  "foreach", "gridExtra", "parallel", "extrafont", "patchwork", "ggrepel",
  "cowplot", "magick", "ggtext", "plyr"
)

# Install missing packages
installed_packages <- rownames(installed.packages())
for (pkg in required_packages) {
  if (!pkg %in% installed_packages) {
    install.packages(pkg)
  }
}

# Load required libraries
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(data.table)


# library(beepr)  # Optional: uncomment if needed






windowsFonts("Helvetica" = windowsFont("Arial"))
my_theme <- function() {
  theme_bw(base_family = "Helvetica") %+replace%
    theme(title = element_text(size=7),
          panel.border = element_blank(),
          # panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "gray20", size = 0,5),
          axis.text = element_text(size =7),
          axis.title = element_text(size = 6),
          legend.key = element_blank(),
          #   legend.key.size = unit (2.25 "cm"),
          legend.margin = margin (0.5,0.5,0.5,0.5),
          legend.text = element_text(size = 6),
          legend.title = element_text(size =7),
          strip.text = element_text(size = 6),
          strip.background = element_blank(),
          complete = TRUE,
          plot.tag = element_text(size = 6))
}
width = 21.59 
heigth = 19.05 * 21.59 /  25.4 
height = width/5*3

letterincheswidth =  19.05/2.54
letterinchesheight =  21.59/2.54 * 3/4






setwd("C:/Users/Raul/Dropbox/Papers/DIA-NN/Figures/")

####################
#This could be its own script, but it also takes just seconds
# it provides the gene name associated to each protein ID from the fasta
#This will help relabel the peptide with the most frequent mutation
headers <- grep("^>", readLines("Data_required/Hotspot_and_Gene_Fusion_Analysis/SHIROKANE_06242025_pdxscript/UP000005640_9606_downloaded03072025_oneline.fasta"), value = TRUE)
headersplit =   str_split(headers, "\\|")
Protein_ID = rep("",length(headersplit) )
Gene_name = rep("",length(headersplit) )
for(cont in 1:length(headersplit)){
  Protein_ID[cont] = headersplit[[cont]][2]
  Gene_name[cont] = sub(".*GN=([^ ]+).*", "\\1", headersplit[[cont]][3])
  
}

Id_genematch = data.frame(Gene = Gene_name, Protein_ID = Protein_ID)

##until here
####################




noncanonical_peptides = data.frame(fread("Data_required/Hotspot_and_Gene_Fusion_Analysis/SHIROKANE_06242025_pdxscript/non_canonical_sequences_justsequences.txt", header = F, sep = "\t"))
noncanonical_peptides = noncanonical_peptides[grep(":", noncanonical_peptides$V1),, drop = FALSE]


# Extraction of the mutated peptides from the output file obtained in the post processsing
#e.g.,  Q13485_D537_V:5_ALQLLVEVLHTMPIADPQPLD_3 --> ALQLLVEVLHTMPIADPQPLD
noncanonical_peptides_sequence = str_split(noncanonical_peptides$V1, "_")
#we want the sequence that is before-last
#in case of different lengths of the rows after splitting by "_", we extract in a loop
lengths = unlist(lapply(noncanonical_peptides_sequence, length))
noncanonical_peptides_sequenceonly = c()
for(cont in 1: length(lengths)){
  noncanonical_peptides_sequenceonly= c(noncanonical_peptides_sequenceonly,  noncanonical_peptides_sequence[[cont]][lengths[cont]-1])
}


#how many mutated peptides were found:
#length(unique(noncanonical_peptides_sequenceonly))


outputDIANN =data.frame(fread("Data_required/Hotspot_and_Gene_Fusion_Analysis/SHIROKANE_06242025_pdxscript/Reports/report_peptidoforms.pg_matrix.tsv"), check.names=FALSE)


#extraction of numeric columns for normalization
numericones  =grep("raw.dia", colnames(outputDIANN))
numericoutputDIANN = outputDIANN[,numericones]
maxnumericoutputDIANN = colSums(numericoutputDIANN,na.rm=T)
normalizednumericoutputDIANN = t(t(numericoutputDIANN)/ maxnumericoutputDIANN)  * 1000000



#normalizednumericoutputDIANN = numericoutputDIANN

#this extarcts the sample names and uses them for naming the 
#numeric columns
#HOWEVER THIS MIGHT ONLY WORK WITH OUR CURRENT CELL LINES
# Extract just the file name (no path, no extension)
colnames(normalizednumericoutputDIANN) <- tools::file_path_sans_ext(basename(colnames(outputDIANN)[numericones]))


normalizednumericoutputDIANN = data.frame(normalizednumericoutputDIANN)
grep("raw.dia", colnames(outputDIANN),invert = TRUE)


#metadata: columns that are not part of the numeric columns
metadata = outputDIANN[,grep("raw.dia", colnames(outputDIANN),invert = TRUE)]
# Then combine everything together
normalizednumericoutputDIANN = cbind(normalizednumericoutputDIANN, metadata)


normalizednumericoutputDIANN = normalizednumericoutputDIANN[grep("KEAP1", normalizednumericoutputDIANN$Protein.Names),]




myColors <- c(
  "#E41A1C",  # Red
  "#377EB8",  # Blue
  "#4DAF4A",  # Green
  "#984EA3",  # Purple
  "#FF7F00",  # Orange
  "#E5C100",  # Muted Yellow
  "#A65628",  # Brown
  "#F781BF",  # Pink
  "#999999",  # Gray
  "#1B9E77",  # Teal
  "#D95F02",  # Dark Orange
  "#7570B3",  # Deep Blue
  "#66C2A5",  # Soft Cyan
  "#0033A0",  # Intense Blue
  "#F4A6D7",  # Pastel Pink
  "#FC8D62",  # Coral
  "#8DD3C7",  # Aqua
  "#FFFFB3",  # Light Yellow
  "#BEBADA",  # Lavender
  "#FB8072",  # Salmon
  "#80B1D3",  # Sky Blue
  "#FDB462",  # Light Orange
  "#B3DE69",  # Light Green
  "#FCCDE5",  # Light Pink
  "#D9D9D9",  # Light Gray
  "#BC80BD"   # Soft Purple
)






#pdf("Fig_KEAP1.pdf", width  = letterincheswidth * 0.5  , height =letterinchesheight * 0.5)
png("wakatefigures/Fig_KEAP1_WT_caseforWakate.png", width  = letterincheswidth * 0.5  , height =letterinchesheight * 0.5,units = "in", res = 1000 )

cont = grep("KEAP1",noncanonLabel )

thecanon = normalizednumericoutputDIANN

thisselectionandcanon = thecanon
thisselectionandcanon$Label = "KEAP1"

# First, define the columns that are numeric (you can also select them by name range if needed)
numeric_cols <- c("HeLa.1", "HeLa.2", "A549.1", "A549.2", "CCRFCEM.1", "CCRFCEM.2",
                  "COLO205.1", "COLO205.2", "H226.1", "H226.2", "H23.1", "H23.2",
                  "RPMI8226.1", "RPMI8226.2", "T47D.1", "T47D.2", "4pool.1", "4pool.2",
                  "NCI7ref.1", "NCI7ref.2")

## Then, group by Label and keep the row with the highest numeric sum
#thisselectionandcanon <- thisselectionandcanon %>%
#  rowwise() %>%
#  mutate(row_sum = sum(c_across(all_of(numeric_cols)), na.rm = TRUE)) %>%
#  ungroup() %>%
#  group_by(Label) %>%
#  slice_max(order_by = row_sum, n = 1, with_ties = FALSE) %>%
#  ungroup() %>%
#  select(-row_sum) 


meltselected_normalizednumericoutputDIANN_thisprot = reshape2::melt(thisselectionandcanon)
meltselected_normalizednumericoutputDIANN_thisprot$variable = factor(meltselected_normalizednumericoutputDIANN_thisprot$variable , levels = unique(meltselected_normalizednumericoutputDIANN_thisprot$variable)[length(unique(meltselected_normalizednumericoutputDIANN_thisprot$variable )):1])


meltselected_normalizednumericoutputDIANN_thisprot$Label <- meltselected_normalizednumericoutputDIANN_thisprot$Label |>
  sub("_", " p.", x = _) |>                             # First "_" becomes ": p."
  sub(":[0-9]+_", ": ", x = _) |>                          # Remove number + "_" after ":"
  sub("_", "", x = _)                                     # Remove second "_"

#meltselected_normalizednumericoutputDIANN_thisprot$Label <- sub("^([^:]+) p\\.", "\\1", meltselected_normalizednumericoutputDIANN_thisprot$Label)

meltselected_normalizednumericoutputDIANN_thisprot$Label[grep(":", meltselected_normalizednumericoutputDIANN_thisprot$Label  , invert = TRUE)] = meltselected_normalizednumericoutputDIANN_thisprot$Label[grep(":", meltselected_normalizednumericoutputDIANN_thisprot$Label  , invert = TRUE)] |>
  sub("p.","",x = _)

#dummy_row <- data.frame(
#  Gene_and_mut = "",
#  Canon = factor("TRUE", levels = c("TRUE", "FALSE")),  # Force "TRUE" into factor levels
#  Label =meltselected_normalizednumericoutputDIANN_thisprot$Label[1],      # Won't appear on the plot
#  Sequence = "",
#  variable = meltselected_normalizednumericoutputDIANN_thisprot$variable[1],   # Won't appear on the plot
#  value = NA    # Won't appear on the plot
#  
#)
# Combine real and dummy data
plot_data <- meltselected_normalizednumericoutputDIANN_thisprot
#plot_data = plot_data[plot_data$Canon == TRUE | plot_data$Gene_and_mut == "KEAP1_G12_A:134_LVVVGAAGVGK",]
plot_data$Status = "Mutated"
plot_data$Status[plot_data$Canon == "TRUE"] = "Not Mutated"
plot_data$Status = factor(plot_data$Status, levels = c("Not Mutated", "Mutated"))
plot_data= plot_data[order(plot_data$Status ),] 
#plot_data$Label[plot_data$Label == "KEAP1 p.G12A: LVVVGAAGVGK"] = "KEAP1 LVVVGAGGVGK"
plot_data$variable = str_split_fixed(str_split_fixed(as.character(plot_data$variable), "_", 5)[,5], ".raw",2)[,1]
plot_data$variable [1:2] = c("HeLa.1", "HeLa.2")
plot_data$Label = factor(plot_data$Label, levels = unique(plot_data$Label))

myColors = c("#303030", myColors)
p = ggplot(plot_data, aes(x = variable, y = value , color = Label)) +
  geom_point(size = 1.5) +
  theme_minimal() +
  coord_flip() +
  xlab("Cell Line") +
  ylab("Normalized Signal") +
  scale_color_manual(values = myColors) +
  ggtitle("KEAP1 WT protein signal") +
  # scale_shape_manual(values = c("Not Mutated" = 16, "Mutated" = 17), drop = FALSE) +
  expand_limits(y = 0) +
  my_theme() + 
  theme(
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6)) 
print(p)

dev.off()



