# Previously labeled as: PREX2_PDX_Geneexpression_forpaper.R
# PREX2 protein expression across 49 PDX tumor samples

# --- Load Libraries ---
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(data.table)
library(grid)

# --- Plot Theme and Dimensions ---
windowsFonts("Helvetica" = windowsFont("Arial"))
my_theme <- function() {
  theme_bw(base_family = "Helvetica") %+replace%
    theme(title = element_text(size=7),
          panel.border = element_blank(),
          axis.line = element_line(colour = "gray20", size = 0,5),
          axis.text = element_text(size =7),
          axis.title = element_text(size = 7),
          legend.key = element_blank(),
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

# --- Load DIA-NN Protein Group Matrix ---
outputDIANN =data.frame(fread("Data_required/Adachi_PDX/Reports/report_peptidoforms.pg_matrix.tsv"), check.names=FALSE)

# --- Column-Sum Normalization ---
numericones  =grep("raw.dia", colnames(outputDIANN))
numericoutputDIANN = outputDIANN[,numericones]
maxnumericoutputDIANN = colSums(numericoutputDIANN,na.rm=T)
normalizednumericoutputDIANN = t(t(numericoutputDIANN)/ maxnumericoutputDIANN)  * 1000000

# --- Extract Sample Names for Column Labels ---
colnames(normalizednumericoutputDIANN) <- tools::file_path_sans_ext(basename(colnames(outputDIANN)[numericones]))

normalizednumericoutputDIANN = data.frame(normalizednumericoutputDIANN)

# --- Combine Normalized Data with Metadata ---
metadata = outputDIANN[,grep("raw.dia", colnames(outputDIANN),invert = TRUE)]
normalizednumericoutputDIANN = cbind(normalizednumericoutputDIANN, metadata)

# --- Color Palette ---
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

# --- Filter for PREX2 Protein and Reshape ---
PREX2case  = normalizednumericoutputDIANN[grep("PREX2", normalizednumericoutputDIANN$Genes),,drop= FALSE]
meltselected_normalizednumericoutputDIANN_thisprot = reshape2::melt(PREX2case)

# Clean sample names: J. -> J-, _1 -> .1, remove .raw suffix
meltselected_normalizednumericoutputDIANN_thisprot$variable <- gsub(
  "\\.raw$", "",
  gsub("_(\\d)", ".\\1",
       gsub("^J\\.", "J-", meltselected_normalizednumericoutputDIANN_thisprot$variable)
  )
)

# Reverse order for display and set factor levels
meltselected_normalizednumericoutputDIANN_thisprot = meltselected_normalizednumericoutputDIANN_thisprot[length(meltselected_normalizednumericoutputDIANN_thisprot$variable):1,]
meltselected_normalizednumericoutputDIANN_thisprot$variable <- factor(
  meltselected_normalizednumericoutputDIANN_thisprot$variable,
  levels = rev(sort(unique(meltselected_normalizednumericoutputDIANN_thisprot$variable)))
)

# --- Generate PREX2 Expression Figure ---
png("Fig_test_PREX2_PDXtest.png", width  = letterincheswidth  , height =letterinchesheight , units = "in", res = 1000 )

plot_data = data.frame(variable = meltselected_normalizednumericoutputDIANN_thisprot$variable,
                       value = meltselected_normalizednumericoutputDIANN_thisprot$value)
p = ggplot(plot_data, aes(x = variable, y = value)) +
  geom_point(size = 1.5) +
  theme_minimal() +
  xlab("PDX Sample") +
  ylab("Normalized Signal") +
  scale_color_manual(values = myColors) +
  ggtitle("PREX2 Normalized SIgnal in PDX Samples") +
  expand_limits(y = 0) +
  my_theme() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing.y = unit(1, "mm"),
    legend.key.height = unit(2, "mm")
  ) +
  guides(
    color = guide_legend(ncol = 2)
  ) + coord_flip()

print(p)
dev.off()
