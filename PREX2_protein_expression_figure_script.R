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
          axis.title = element_text(size = 7),
          legend.key = element_blank(),
          #   legend.key.size = unit (2.25 "cm"),
          legend.margin = margin (0.5,0.5,0.5,0.5),
          legend.text = element_text(size = 7),
          legend.title = element_text(size =7),
          strip.text = element_text(size = 7),
          strip.background = element_blank(),
          complete = TRUE,
          plot.tag = element_text(size = 7))
}
width = 21.59 
heigth = 19.05 * 21.59 /  25.4 
height = width/5*3

letterincheswidth =  21.59/2.54
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


outputDIANN =data.frame(fread("Data_required/Adachi_PDX/Reports/report_peptidoforms.pg_matrix.tsv"), check.names=FALSE)


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


PREX2case  = normalizednumericoutputDIANN[grep("PREX2", normalizednumericoutputDIANN$Genes),,drop= FALSE]
meltselected_normalizednumericoutputDIANN_thisprot = reshape2::melt(PREX2case)
# Apply transformations in order
# Apply transformations in order
meltselected_normalizednumericoutputDIANN_thisprot$variable <- gsub(
  "\\.raw$", "", 
  gsub("_(\\d)", ".\\1", 
       gsub("^J\\.", "J-", meltselected_normalizednumericoutputDIANN_thisprot$variable)
  )
)


meltselected_normalizednumericoutputDIANN_thisprot = meltselected_normalizednumericoutputDIANN_thisprot[length(meltselected_normalizednumericoutputDIANN_thisprot$variable):1,]
meltselected_normalizednumericoutputDIANN_thisprot$variable <- factor(
  meltselected_normalizednumericoutputDIANN_thisprot$variable,
  levels = rev(sort(unique(meltselected_normalizednumericoutputDIANN_thisprot$variable)))
)


png("Fig_PREX2_Protein_expression.png", width  = letterincheswidth  , height =letterinchesheight*1.5,units = "in", res = 1000 )

  p = ggplot(meltselected_normalizednumericoutputDIANN_thisprot, aes(x = variable, y = value )) +
    geom_point(size = 3) +
    theme_minimal() +
    coord_flip() +
    xlab("PDX Samples") +
    ylab("Normalized Signal") +
    scale_color_manual(values = myColors) +
    ggtitle("PREX2 expression across PDX samples") +
    expand_limits(y = 0) +
    my_theme()
  print(p)
  
dev.off()

