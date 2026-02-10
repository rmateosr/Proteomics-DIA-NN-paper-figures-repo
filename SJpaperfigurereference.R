library(stringr)
library(RColorBrewer)
library(foreach)
library(ggplot2)
library(gridExtra)
library(data.table)
library(dplyr)
library(beepr)
library(parallel)
library(extrafont)
library(patchwork)
library(ggrepel)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(magick)
library(ggtext)
library(plyr)
layout <- "
BBBBB
CDD##
"


windowsFonts("Helvetica" = windowsFont("Arial"))
my_theme <- function() {
  theme_bw(base_family = "Helvetica") %+replace%
    theme(title = element_text(size=7),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
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




#!/usr/bin/Rscript
set.seed(12)  
library(stringr)
library(RColorBrewer)
library(foreach)
#library(ggsci)
library(ggplot2)
library(gridExtra)
#library(viridis)
#library(gplots)
library(data.table)
library(dplyr)
library(beepr)
library(parallel)
Sjscorethreshold = 10
scorediffthreshold= 10
c = 0.001
dir.create("C:/Users/Raul/Documents/papers/SJ_paper/Ambiguous_082023/figures/fig3SJincluded/")
metadatafile = "C:/Users/Raul/Documents/papers/SJ_paper/Ambiguous_082023/figures/Final_mutmetadata_ambiguousintegrated.tsv"
mutmetadata  =data.frame(fread(metadatafile))
#if(sum(mutmetadata$relevance == "Ambiguous")>0){


PANCANCERNFE2L2 =readRDS("C:/Users/Raul/Documents/papers/SJ_paper/Ambiguous_082023/figures/mc3.v0.2.8.PUBLIC_NRF2.maf.RDS")
npvalue = 1
resultadosnumericPANCANCER =   data.frame(fread("C:/Users/Raul/Documents/papers/SJ_paper/Ambiguous_082023/figures/Scores_CV.tsv"))

resultadosnumericPANCANCER = unique(resultadosnumericPANCANCER)
resultadosnumericPANCANCER= resultadosnumericPANCANCER[!duplicated(resultadosnumericPANCANCER$TCGA_ID),]
colnames(resultadosnumericPANCANCER)[3] <- "Mutation Labels"
resultadosnumericPANCANCERforplot = resultadosnumericPANCANCER[order(resultadosnumericPANCANCER$score, resultadosnumericPANCANCER$`Mutation Labels`, decreasing = T),]
resultadosnumericPANCANCERforplot$`Mutation Labels` = as.character(resultadosnumericPANCANCERforplot$`Mutation Labels`)
resultadosnumericPANCANCERforplot$`Mutation Labels`[resultadosnumericPANCANCERforplot$`Mutation Labels`=="NFE2L2_Exon2_Mut"] <- "NRF2_Exon2_Mut"
resultadosnumericPANCANCERforplot$`Mutation Labels`[resultadosnumericPANCANCERforplot$`Mutation Labels`=="NFE2L2_Other_Exon_Mut"] <- "NRF2_Other_Exon_Mut"
resultadosnumericPANCANCERforplot$`Mutation Labels`[resultadosnumericPANCANCERforplot$`Mutation Labels`=="Exon2_Skipping"] <- "NRF2_Exon2_Skipping"
old = c("None" ,"NRF2_Exon2_Mut",  "NRF2_Other_Exon_Mut" ,"KEAP1" , "CUL3"  ,"CNA" ,"NRF2_Exon2_Skipping", "Multihit")
old = c("None" ,"NRF2_Exon2_Hotspot_GOF",  "NRF2_Exon2_Mut" , "NRF2_Oth	er_Exon_Mut" ,"KEAP1" ,"KEAP1 OL" ,"KEAP1_Ambiguous",  "CUL3"  ,"CNA" ,"NRF2_Exon2_Skipping", "Multihit")
new = c("None" ,"NRF2 Exon 2 Mut",  "NRF2 Other Exon Mut" ,"KEAP1 Mut" ,"KEAP1_OL" ,"KEAP1_Ambiguous", "CUL3 Mut"  ,"NRF2 CNA" ,"NRF2 Exon 2 Skipping", "Multihit")
new = c("None" ,"NRF2 Exon 2 Hotspot GOF", "NRF2 Exon2 Mut" , "NRF2 Other Exon Mut" ,"KEAP1 Mut" , "KEAP1_OL" ,"KEAP1_Ambiguous", "CUL3 Mut"  ,"NRF2 CNA" ,"NRF2 Exon 2 Skipping", "Multihit")
new = c("None" ,"NRF2 E2 HGOF", "NRF2 E2 Amb" , "NRF2 oE" ,"KEAP1 Mut" , "KEAP1 OLO" ,"KEAP1 Amb", "CUL3"  ,"NRF2 CNA" ,"NRF2 E2 Skip", "Multihit")

new = c("None" ,"NRF2 E2 HGOF", "NRF2 E2 Amb" , "NRF2 oE Amb" ,"KEAP1 Mut" , "KEAP1 OLO" ,"KEAP1 Amb", "CUL3 Amb"  ,"NRF2 CNA" ,"NRF2 E2 Skip", "Multihit")
old = c("None" ,"NRF2_Exon2_Hotspot_GOF",  "NRF2_Exon2_Mut" , "NRF2_Other_Exon_Mut" ,"KEAP1" ,"KEAP1 OL" ,"KEAP1_Ambiguous",  "CUL3"  ,"CNA" ,"NRF2_Exon2_Skipping", "Multihit")
new = c("None" ,"NRF2 E2 H", "NRF2 E2 nH" , "NRF2 oE" ,"KEAP1 Mut" , "KEAP1 OLO" ,"KEAP1 nOLO", "CUL3"  ,"NRF2 CNA" ,"NRF2 E2 Skip", "Multihit")

new = c("None" ,"NRF2 Exon 2 Hotspot", "NRF2 Exon 2 non-Hotspot" , "NRF2 Other Exon Mut" ,"KEAP1 Mut" , "KEAP1 Oncogenic or Likely Oncogenic" ,"KEAP1 non-Oncogenic", "CUL3 mut"  ,"NRF2 CNA" ,"NRF2 Exon 2 Skipping", "Multihit")

new = c("None" ,"NFE2L2 Exon 2 Hotspot", "NFE2L2 Exon 2 non-Hotspot" , "NFE2L2 Other Exon Mut" ,"KEAP1 Mut" , "KEAP1 Oncogenic or Likely Oncogenic" ,"KEAP1 non-Oncogenic", "CUL3 Mut"  ,"NFE2L2 CNA" ,"NFE2L2 Exon 2 Skipping", "Multihit")


######







######  






for(cont in 1:length(old)){
  resultadosnumericPANCANCERforplot$`Mutation Labels`[resultadosnumericPANCANCERforplot$`Mutation Labels` == old[cont]] = new[cont]
  
  
}
#resultadosnumericPANCANCERforplot$`Mutation Labels` =  factor(resultadosnumericPANCANCERforplot$`Mutation Labels`, levels = c("None" ,"NRF2_Exon2_Mut",  "NRF2_Other_Exon_Mut" ,"KEAP1" , "CUL3"  ,"CNA" ,"NRF2_Exon2_Skipping", "Multihit"))

resultadosnumericPANCANCERforplot$`Mutation Labels` =  factor(resultadosnumericPANCANCERforplot$`Mutation Labels`, levels = new[new %in% resultadosnumericPANCANCERforplot$`Mutation Labels`])


resultadosnumericPANCANCERforplot$is_mut_predicted = resultadosnumericPANCANCERforplot$score >scorediffthreshold
resultadosnumericPANCANCERforplot$is_mut_observed = resultadosnumericPANCANCERforplot$`Mutation Labels` != "None"

resultadosnumericPANCANCERforplot$is_mut_predicted = factor(resultadosnumericPANCANCERforplot$is_mut_predicted, levels = c(FALSE, TRUE))


resultadosnumericPANCANCERforplot$is_mut_observed = factor(resultadosnumericPANCANCERforplot$is_mut_observed, levels = c(FALSE, TRUE))


#colnames(resultstest) <- c("Observed Not Mutated", "Observed Mutated")
resultadosnumericPANCANCERforplot$`Mutation labels Ambiguous` = as.character(resultadosnumericPANCANCERforplot$`Mutation Labels`)

#resultadosnumericPANCANCERforplot$`Mutation labels Ambiguous`[rowSums(resultadosnumericPANCANCERforplot[,c(7,9:12)])>1]<- "True Multihit"
#
#resultadosnumericPANCANCERforplot$`Mutation labels Ambiguous`[resultadosnumericPANCANCERforplot$`Mutation labels Ambiguous` == "Multihit"] <- colnames(resultadosnumericPANCANCERforplot[resultadosnumericPANCANCERforplot$`Mutation labels Ambiguous` == "Multihit",c(7,9:12)] )[which(resultadosnumericPANCANCERforplot[resultadosnumericPANCANCERforplot$`Mutation labels Ambiguous` == "Multihit",c(7,9:12)] == TRUE,arr.ind = T)[,2]]
#resultadosnumericPANCANCERforplot$`Mutation labels Ambiguous`[resultadosnumericPANCANCERforplot$`Mutation labels Ambiguous` == "True Multihit"] <- "Multihit"

levelsambiguous = c("NFE2L2 Exon2 Mut","NFE2L2 Exon2 Hotspot GOF ", "NFE2L2 Other Exon Mut" ,"CUL3 Mut", "NFE2L2 CNA" ,"NFE2L2 Exon 2 Skipping", "Multihit")
levelsambiguous = c("NFE2L2 Exon2 non-Hotspot","NFE2L2 Other Exon Mut", "KEAP1 non-Oncogenic" ,"CUL3 Mut", "NFE2L2 CNA" ,"NFE2L2 Exon 2 Skipping", "Multihit")
levelsambiguous = levels(resultadosnumericPANCANCERforplot$`Mutation Labels`) [levels(resultadosnumericPANCANCERforplot$`Mutation Labels`) %in% unique(resultadosnumericPANCANCERforplot$`Mutation Labels`)]

resultadosnumericPANCANCERforplot$`Mutation labels Ambiguous` = factor(resultadosnumericPANCANCERforplot$`Mutation labels Ambiguous`, levels = levelsambiguous)
resultadosnumericPANCANCERforplot$`Mutation Labels` <- resultadosnumericPANCANCERforplot$`Mutation labels Ambiguous`
resultstest = table(resultadosnumericPANCANCERforplot$is_mut_predicted,resultadosnumericPANCANCERforplot$`Mutation Labels`)
rownames(resultstest) <- c("Predicted Not Mutated", "Predicted Mutated")

tt = ttheme_default(base_size = 7)

ggbx =ggplot(resultadosnumericPANCANCERforplot, aes(x = ((1:length(`Mutation Labels`))-1) , y =  score, color = `Mutation Labels`, fill = `Mutation Labels`)) + geom_bar(stat="identity",position="dodge", width= 0.5) +xlab("Sample") + ylab("Score") +ggtitle(paste0("TCGA dataset mutation prediction score SJ classifier"))+
  annotation_custom( tableGrob(resultstest,theme =tt), xmin=length(resultadosnumericPANCANCERforplot$score)/2-1, xmax=length(resultadosnumericPANCANCERforplot$score)[1]/2, ymin=max(resultadosnumericPANCANCERforplot$score)/2-1, ymax=max(resultadosnumericPANCANCERforplot$score/2))
#ggbx = ggbx +scale_fill_nejm()

#resultadosnumericPANCANCERforplot= resultadosnumericPANCANCERforplot[1:100,]
ggbx =ggplot(resultadosnumericPANCANCERforplot, aes(x = ((1:length(`Mutation Labels`))-1) , y =  score, color = `Mutation Labels`, fill = `Mutation Labels`)) + geom_bar(stat="identity",position="dodge", width= 0.18) +xlab("Sample") + ylab("Score") +ggtitle(paste0("TCGA dataset mutation prediction score SJ classifier"))+
  annotation_custom( tableGrob(t(resultstest),theme =tt), xmin=length(resultadosnumericPANCANCERforplot$score)/2-1, xmax=length(resultadosnumericPANCANCERforplot$score)[1]/2, ymin=max(resultadosnumericPANCANCERforplot$score)*0.7 , ymax=max(resultadosnumericPANCANCERforplot$score*0.7))
#ggbx = ggbx +scale_fill_nejm()numericPANCANCERforplot$score)*0.7 , ymax=max(resultadosnumericPANCANCERforplot$score*0.7))
#ggbx = ggbx +scale_fill_n
myColors <- brewer.pal(7,"Set1")
myColors <- brewer.pal(7,"Set1")[c(2,4,5,6,7)]
myColors <- brewer.pal(7,"Set1")[c(1:7)]
myColors <-myColors[-3] 
myColors[1] <-"hotpink"
#myColors = c("Darkgray", myColors)
myColors = c(myColors[1:2], "darkseagreen",myColors[3:length(myColors)])
myColors[6] = "#FFD700"

ggbx = ggbx + scale_fill_manual(labels = c("*NFE2L2* Exon 2 non-Hotspot","*NFE2L2* Other Exon Mut", "*KEAP1* non-Oncogenic" ,"*CUL3* Mut", "*NFE2L2* CNA" ,"*NFE2L2* Exon 2 Skipping", "Multihit"),values = myColors, drop = FALSE) + scale_color_manual(labels =c("*NFE2L2* Exon 2 non-Hotspot","*NFE2L2* Other Exon Mut", "*KEAP1* non-Oncogenic" ,"*CUL3* Mut", "*NFE2L2* CNA" ,"*NFE2L2* Exon 2 Skipping", "Multihit"), values = myColors, drop = FALSE)+ theme_minimal()  + geom_hline(yintercept=10, linetype="dashed", color = "red")+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank()  #remove y axis ticks
  )
p1 = ggbx + my_theme() +  geom_hline(yintercept = 0, color = "black")+ 
  theme(legend.position = c(0.85,0.75) ,axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(), axis.line.y = element_blank())+ scale_x_continuous(expand = c(0, 0), limits = c(-1, NA)) +
  theme(legend.text = element_markdown())
p1



#stat="identity",position="dodge", width=0.02) +xlab("Sample") + ylab("Score") +ggtitle(paste0("TCGA dataset mutation prediction score SJ classifier"))+
#  annotation_custom( tableGrob(resultstest), xPANCANCERNFE2L2min=length(resultadosnumericPANCANCERforplot$score)/2-1, xmax=length(resultadosnumericPANCANCERforplot$score)[1]/2, ymin=max(resultadosnumericPANCANCERforplot$score)/2-1, ymax=max(resultadosnumericPANCANCERforplot$score/2))
#


#!/usr/bin/Rscript
set.seed(12)   
args = commandArgs(trailingOnly=TRUE)
Sjscorethreshold = 10
scorediffthreshold= 10
c = 0.001
library(stringr)
library(RColorBrewer)
library(foreach)
library(ggplot2)
library(gridExtra)
library(data.table)
library(dplyr)
library(beepr)
library(parallel)
metadatafile = "C:/Users/Raul/Documents/papers/SJ_paper/Ambiguous_082023/figures/Final_mutmetadata_ambiguousintegrated.tsv"
mutmetadata  =data.frame(fread(metadatafile))
pthresholds =  "p0001"
PANCANCERNFE2L2 = readRDS("C:/Users/Raul/Documents/papers/SJ_paper/Ambiguous_082023/figures/mc3.v0.2.8.PUBLIC_NRF2.maf.RDS")
resultadosnumericPANCANCER =   data.frame(fread("C:/Users/Raul/Documents/papers/SJ_paper/Ambiguous_082023/figures/Scores_CV.tsv"))


resultadosnumericPANCANCER = unique(resultadosnumericPANCANCER)
resultadosnumericPANCANCER= resultadosnumericPANCANCER[!duplicated(resultadosnumericPANCANCER$TCGA_ID),]
resultadosnumericPANCANCER$mutlabels[resultadosnumericPANCANCER$mutlabels %in% unique(mutmetadata$mutlabels[mutmetadata$relevance == "Not_relevant"])] <- "None"

colnames(resultadosnumericPANCANCER)[3] <- "Mutation Labels"


####Relabeling of mutations#####

resultadosnumericPANCANCERforplot = resultadosnumericPANCANCER[order(resultadosnumericPANCANCER$score, resultadosnumericPANCANCER$`Mutation Labels`, decreasing = T),]
resultadosnumericPANCANCERforplot$`Mutation Labels` = as.character(resultadosnumericPANCANCERforplot$`Mutation Labels`)
resultadosnumericPANCANCERforplot$`Mutation Labels`[resultadosnumericPANCANCERforplot$`Mutation Labels`=="NFE2L2_Exon2_Mut"] <- "NRF2_Exon2_Mut"
resultadosnumericPANCANCERforplot$`Mutation Labels`[resultadosnumericPANCANCERforplot$`Mutation Labels`=="NFE2L2_Other_Exon_Mut"] <- "NRF2_Other_Exon_Mut"
resultadosnumericPANCANCERforplot$`Mutation Labels`[resultadosnumericPANCANCERforplot$`Mutation Labels`=="Exon2_Skipping"] <- "NRF2_Exon2_Skipping"

old = c("None" ,"NRF2_Exon2_Mut",  "NRF2_Other_Exon_Mut" ,"KEAP1" , "CUL3"  ,"CNA" ,"NRF2_Exon2_Skipping", "Multihit")
old = c("None" ,"NRF2_Exon2_Hotspot_GOF",  "NRF2_Exon2_Mut" , "NRF2_Other_Exon_Mut" ,"KEAP1_Mut" ,"KEAP1 OL" ,"KEAP1 Ambiguous",  "CUL3"  ,"CNA" ,"NRF2_Exon2_Skipping", "Multihit")
new = c("None" ,"NRF2 Exon 2 Mut",  "NRF2 Other Exon Mut" ,"KEAP1 Mut" ,"KEAP1_OL" ,"KEAP1_Ambiguous", "CUL3 Mut"  ,"NRF2 CNA" ,"NRF2 Exon 2 Skipping", "Multihit")
new = c("None" ,"NRF2 Exon2 GOF Mut", "NRF2 Exon2 Mut" , "NRF2 Other Exon Mut" ,"KEAP1 Mut" , "KEAP1_OL" ,"KEAP1_Ambiguous", "CUL3 Mut"  ,"NRF2 CNA" ,"NRF2 Exon 2 Skipping", "Multihit")

new = c("None" ,"NFE2L2 Exon 2 Hotspot", "NFE2L2 Exon2 Mut" , "NFE2L2 Other Exon Mut" ,"KEAP1" , "KEAP1_OL" ,"KEAP1_Ambiguous", "CUL3 Mut"  ,"NFE2L2 CNA" ,"NFE2L2 Exon 2 Skip", "Multihit")



#Label multihit only with KEAP1 and NRF2 GOF labels
#currently only KEAP1 cases but might need to include GOF in the future

mutmetadata[mutmetadata$TCGA_ID %in% resultadosnumericPANCANCERforplot$TCGA_ID[resultadosnumericPANCANCERforplot$`Mutation Labels` == "Multihit"],]
resultadosnumericPANCANCERforplot =inner_join(resultadosnumericPANCANCERforplot, mutmetadata)

if(length(unique(resultadosnumericPANCANCERforplot$`Mutation Labels`)) <=3){
  resultadosnumericPANCANCERforplot$`Mutation Labels`[resultadosnumericPANCANCERforplot$`Mutation Labels` == "Multihit" & resultadosnumericPANCANCERforplot$isKEAP1OL == "KEAP1_OL"] <- "KEAP1_OL"
  resultadosnumericPANCANCERforplot$`Mutation Labels`[resultadosnumericPANCANCERforplot$`Mutation Labels` == "Multihit" & resultadosnumericPANCANCERforplot$ismutexon2 ==  "NFE2L2_Exon2_Hotspot_GOF"] <- "NFE2L2_Exon2_Hotspot_GOF"
  resultadosnumericPANCANCERforplot$`Mutation Labels`[resultadosnumericPANCANCERforplot$isKEAP1OL == "KEAP1_OL" & resultadosnumericPANCANCERforplot$ismutexon2 ==  "NFE2L2_Exon2_Hotspot_GOF"] <- "Multihit"
  
}



for(cont in 1:length(old)){
  resultadosnumericPANCANCERforplot$`Mutation Labels`[resultadosnumericPANCANCERforplot$`Mutation Labels` == old[cont]] = new[cont]
}

resultadosnumericPANCANCERforplot$`Mutation Labels` =  factor(resultadosnumericPANCANCERforplot$`Mutation Labels`, levels =  new[new %in% resultadosnumericPANCANCERforplot$`Mutation Labels`])
#######






resultadosnumericPANCANCERforplot$is_mut_predicted = resultadosnumericPANCANCERforplot$score >scorediffthreshold
resultadosnumericPANCANCERforplot$is_mut_observed = resultadosnumericPANCANCERforplot$`Mutation Labels` != "None"

resultadosnumericPANCANCERforplot$is_mut_predicted = factor(resultadosnumericPANCANCERforplot$is_mut_predicted, levels = c(FALSE, TRUE))


resultadosnumericPANCANCERforplot$is_mut_observed = factor(resultadosnumericPANCANCERforplot$is_mut_observed, levels = c(FALSE, TRUE))

resultstest = table(resultadosnumericPANCANCERforplot$is_mut_predicted,resultadosnumericPANCANCERforplot$is_mut_observed)
rownames(resultstest) <- c("Predicted Not Mutated", "Predicted Mutated")
colnames(resultstest) <- c("Observed Not Mutated", "Observed Mutated")
resultadosnumericPANCANCERforplot$Score = resultadosnumericPANCANCERforplot$score
ggbx =ggplot(resultadosnumericPANCANCERforplot, aes(x = ((1:length(`Mutation Labels`))-1) , y =  Score, color = `Mutation Labels`, fill = `Mutation Labels`)) + geom_bar(stat="identity",position="dodge", width=0.02) +xlab("Sample") + ylab("Score") +ggtitle(paste0("TCGA dataset mutation prediction score SJ classifier"))+
  annotation_custom( tableGrob(resultstest), xmin=length(resultadosnumericPANCANCERforplot$score)/2-1, xmax=length(resultadosnumericPANCANCERforplot$score)[1]/2, ymin=max(resultadosnumericPANCANCERforplot$score)/2-1, ymax=max(resultadosnumericPANCANCERforplot$score/2))


myColors <- brewer.pal(7,"Set1")
myColors = c("Darkgray", myColors)
myColors = myColors[c(1,2,4,7)]
myColors[4] = "#FFD700"
ggbx = ggbx + scale_fill_manual(values = myColors, drop = FALSE) + scale_color_manual(values = myColors, drop = FALSE)+ theme_minimal()  + geom_hline(yintercept=10, linetype="dashed", color = "red")+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank()  #remove y axis ticks
  )



resultadosnumericPANCANCERforplot1  =resultadosnumericPANCANCERforplot

resultadosnumericPANCANCERforplot1$`Mutation Labels` = as.character(resultadosnumericPANCANCERforplot1$`Mutation Labels`)
cursiva = c("None", "*NFE2L2*<br>Exon 2<br>Hotspot", "*KEAP1*", "*NFE2L2*<br>Exon 2<br>Skip")


resultadosnumericPANCANCERforplot1$`Mutation Labels`[resultadosnumericPANCANCERforplot1$`Mutation Labels` == "NFE2L2 Exon 2 Hotspot"] <- "NFE2L2\nExon 2\nHotspot"
resultadosnumericPANCANCERforplot1$`Mutation Labels`[resultadosnumericPANCANCERforplot1$`Mutation Labels` == "NFE2L2 Exon 2 Skip"] <- "NFE2L2\nExon 2\nSkip"

resultadosnumericPANCANCERforplot1$`Mutation Labels` = factor(resultadosnumericPANCANCERforplot1$`Mutation Labels`, levels = c("None", "NFE2L2\nExon 2\nHotspot", "KEAP1", "NFE2L2\nExon 2\nSkip" ))



p3 = ggplot(resultadosnumericPANCANCERforplot1, aes(x = `Mutation Labels` , y =  Score, fill = `Mutation Labels`)) + 
  geom_boxplot( color = "black",lwd=0.1,outlier.size = 0.05)+ scale_fill_manual(values = myColors, drop = FALSE) +
  scale_color_manual(values = myColors, drop = FALSE)+ theme_minimal()  + xlab("")+
  geom_hline(yintercept=10, linetype="dashed", color = "red") +  scale_x_discrete(labels = cursiva)+
 # theme(text = element_text(size = 25))   + 
  my_theme() +
  theme(legend.position = "none") + ggtitle("CV Results Classifier")  + theme(axis.text.x = element_markdown())
  



pdf("C:/Users/Raul/Documents/papers/SJ_paper/Ambiguous_082023/figures/fig3SJincluded/Fig3a.pdf", width  = letterincheswidth* 4/12, height =letterinchesheight/2 )
p3 + labs( tag = "a")
dev.off() 


pdf("C:/Users/Raul/Documents/papers/SJ_paper/Ambiguous_082023/figures/fig3SJincluded/Fig3a_notag.pdf", width  = letterincheswidth* 4/12, height =letterinchesheight/2 )
p3 
dev.off() 

png("C:/Users/Raul/Documents/papers/SJ_paper/Ambiguous_082023/figures/fig3SJincluded/Fig3a.png", width  = letterincheswidth* 4/12, height =letterinchesheight/2,units = "in", res = 1000 )
p3 + labs( tag = "a")
dev.off() 


png("C:/Users/Raul/Documents/papers/SJ_paper/Ambiguous_082023/figures/fig3SJincluded/Fig3a_notag.png", width  = letterincheswidth* 4/12, height =letterinchesheight/2,units = "in", res = 1000 )
p3
dev.off() 



#im.convert("C:/Users/Raul/Documents/papers/SJ_paper/Ambiguous_082023/figures/fig3SJincluded/Fig3a.pdf", output = "C:/Users/Raul/Documents/papers/SJ_paper/Ambiguous_082023/figures/fig3SJincluded/Fig3a_convert.png", extra.opts="-density 150")


######
# 
# o = optim(par = c(0.5,0.5), fn = betabinomial, method = "L-BFGS-B", lower=c(0+c,0+c), upper=c(Inf,Inf) , control = list(maxit = 1000000), total = N , altered = K)
# 
# o 
#set.seed(40)
library(stringr)
library(RColorBrewer)
#library(Rmpfr)
#library(doParallel)
library(foreach)
library(ggsci)
library(ggplot2)
library(gridExtra)
library(viridis)
library(gplots)
library(data.table)
library(dplyr)
#set.seed(40)
library(beepr)
library(PRROC)
setwd("C:/Users/Raul/Documents/NRF2/automaticpipeline/automatic_versions/LikelyandOncogenic_KEAP1OL_nomultihit_rankown_ambiguousintegrated/Precision_Recall_Curve/")
#compilation of all cancer types
GSVAsamples = data.frame(fread("C:/Users/Raul/Documents/NRF2/expression_classifier/cpmlog2GSVA_allcancers.tsv"))
GSVAsamples = unique(GSVAsamples$TCGA_ID)
SJsamples =  data.frame( fread("Scores_CV.tsv"))
SJsamples = unique(SJsamples$TCGA_ID)
sharedsamples  =intersect(GSVAsamples, SJsamples)
setwd("C:/Users/Raul/Documents/NRF2/classifier_pipeline_SRA_based//")
GSVAthres = 0.5
gatheredthres = 10
pancancerthres = 10
#datalist = readRDS("TCGAfromSRAinputfortraining.RDS")
mutmetadata  =fread("C:/Users/Raul/Documents/papers/SJ_paper/Ambiguous_082023/figures/Final_mutmetadata_ambiguousintegrated.tsv")
allresultswithpvalue = data.frame()
alltreholdchosenbycancertype = data.frame()
alltreholdchosenbycancertype_original = data.frame()

#SJgatheredsamples =  data.frame( fread(paste0("C:/Users/Raul/Documents/NRF2/classifier_pipeline_SRA_based_individual_cancers//resultsindividualGATHERED_directresultclassifier_",pthresholds[npvalue],".tsv")))





#compilation of all cancer types
setwd("C:/Users/Raul/Documents/NRF2/automaticpipeline/automatic_versions/LikelyandOncogenic_KEAP1OL_nomultihit_rankown_ambiguousintegrated/Precision_Recall_Curve/")


mutmetadata  =fread("C:/Users/Raul/Documents/papers/SJ_paper/Ambiguous_082023/figures/Final_mutmetadata_ambiguousintegrated.tsv")


#allresultswithpvalue = data.frame()
pthresholds = "p0001"
for(npvalue in 1:length(pthresholds)){
  resultsandmeta =  data.frame( fread(paste0("C:/Users/Raul/Documents/papers/SJ_paper/Ambiguous_082023/figures/Scores_CV.tsv")))
  resultsandmeta= resultsandmeta[resultsandmeta$TCGA_ID %in% sharedsamples,]
  resultsandmeta= resultsandmeta[!duplicated(resultsandmeta$TCGA_ID),]
  thresholds = seq( min(resultsandmeta$score),  max(resultsandmeta$score), length.out =100)  
  cancers = unique(resultsandmeta$cancer)
  resultsandmeta$mutornot = resultsandmeta$mutlabels != "None"
  pr <- pr.curve(scores.class0 = resultsandmeta$score, weights.class0 = resultsandmeta$mutornot)
  SJ_PRAUC = (pr$auc.integral)
  
  
  SJ_PRAUC = c()
  pr <- pr.curve(scores.class0 = resultsandmeta$score, weights.class0 = resultsandmeta$mutornot)
  SJ_PRAUC = c(SJ_PRAUC,(pr$auc.integral))
  for(cont in 1:length(cancers)){
    pr <- pr.curve(scores.class0 = resultsandmeta$score[resultsandmeta$cancer == cancers[cont]], weights.class0 = resultsandmeta$mutornot[resultsandmeta$cancer == cancers[cont]])
    SJ_PRAUC = c(SJ_PRAUC,(pr$auc.integral))
  }
  
  TP = c()
  FP = c()
  TN = c()
  FN = c()
  for(nthres in 1:length(thresholds)) {
    resultsandmeta$SJmutornot =   resultsandmeta$score >= thresholds[nthres]
    TP = c(TP,  sum(resultsandmeta$mutornot == TRUE & resultsandmeta$SJmutornot == TRUE))
    TN = c(TN,  sum(resultsandmeta$mutornot == FALSE & resultsandmeta$SJmutornot == FALSE))
    FN = c(FN,  sum(resultsandmeta$mutornot == TRUE & resultsandmeta$SJmutornot == FALSE))
    FP = c(FP,  sum(resultsandmeta$mutornot == FALSE & resultsandmeta$SJmutornot == TRUE))
  }
  results = data.frame( 
    FALSE_not_called = TN,
    FALSE_called = FP,
    TRUE_not_called = FN,
    TRUE_called = TP,
    sensitivity = TP/(TP + FN),
    specificity = TN/(TN + FP),
    Precision = TP/(TP + FP),
    Recall = TP/(TP + FN),
    FPR =   1-TN/(TN + FP),
    threshold = thresholds, 
    cancer  ="PANCANCER", 
    method = paste0( "SJ-based classifier")
  )
  allresults  = results
  cancers =unique(resultsandmeta$cancer)
  for(ncancers in 1:length(cancers)){
    resultsandmetabycancer  =resultsandmeta[resultsandmeta$cancer == cancers[ncancers],]
    
    
    
    thresholds = seq( min(resultsandmeta$score),  max(resultsandmeta$score), length.out =100)   
    TP = c()
    FP = c()
    TN = c()
    FN = c()
    for(nthres in 1:length(thresholds)) {
      resultsandmetabycancer$SJmutornot =   resultsandmetabycancer$score >= thresholds[nthres]
      TP = c(TP,  sum(resultsandmetabycancer$mutornot == TRUE & resultsandmetabycancer$SJmutornot == TRUE))
      TN = c(TN,  sum(resultsandmetabycancer$mutornot == FALSE & resultsandmetabycancer$SJmutornot == FALSE))
      FN = c(FN,  sum(resultsandmetabycancer$mutornot == TRUE & resultsandmetabycancer$SJmutornot == FALSE))
      FP = c(FP,  sum(resultsandmetabycancer$mutornot == FALSE & resultsandmetabycancer$SJmutornot == TRUE))
    }
    results = data.frame(
      FALSE_not_called = TN,
      FALSE_called = FP,
      TRUE_not_called = FN,
      TRUE_called = TP,
      sensitivity = TP/(TP + FN),
      specificity = TN/(TN + FP),
      Precision = TP/(TP + FP),
      Recall = TP/(TP + FN),
      FPR =   1-TN/(TN + FP),
      threshold = thresholds, 
      cancer  =cancers[ncancers], 
      method = paste0( "SJ-based classifier" )
    )
    
    allresults = rbind(allresults, results)
    allresultswithpvalue = rbind(allresultswithpvalue,allresults)
    
  }
  
  
  
  
  
  treholdchosenbycancertype = c()
  for(ncancers in 1){#:length(unique(allresults$cancer))){
    a = allresults[allresults$cancer == unique(allresults$cancer)[ncancers],]
    a$threshold = round(a$threshold,2)
    fmeasure =  (2 * a$Precision * a$Recall) / (a$Precision + a$Recall)
    # treholdchosenbycancertype = c(treholdchosenbycancertype,  allresults$threshold [which.max(fmeasure)])
    # fmeasures = c(fmeasures , sqrt(a$sensitivity * a$specificity)  )
    valuechosen=  unique(allresults$threshold [which.max(fmeasure)])
    # valuechosen = pancancerthres
    #valuechosen = 10
    #  itscoords = a[a$threshold == (round(valuechosen,digits = 2)),1:2]
    itscoords = a[a$threshold == (round(valuechosen,digits = 2)),]
    if(length(valuechosen) == 0  ){
      valuechosen = NA
      itscoords = c(NA,NA)
    } else if (is.na(valuechosen)){
      valuechosen = NA
      itscoords = c(NA,NA)
    }
    treholdchosenbycancertype =   rbind(treholdchosenbycancertype, (c(unique(allresults$cancer)[ncancers],valuechosen,itscoords )))
    
  }
  
  
  
  alltreholdchosenbycancertype = rbind(alltreholdchosenbycancertype, data.frame(treholdchosenbycancertype,Method = paste0("SJ-based classifier") ))
  
}

SJanalyzed  =resultsandmeta$TCGA_ID










totalresults = data.frame(fread("C:/Users/Raul/Documents/NRF2/expression_classifier/cpmlog2GSVA_allcancers.tsv"))
totalresults= totalresults[totalresults$TCGA_ID %in% sharedsamples,]
totalresults = totalresults[!duplicated(totalresults$TCGA_ID),]
#mutmetadata  =fread("C:/Users/Raul/Documents/NRF2/SJclassifierversions/mutationburden2000/TCGA_mut_metadata_mutburden2023_2000_exskipupdated.tsv")
mutmetadata  =fread("C:/Users/Raul/Documents/NRF2/automaticpipeline/automatic_versions/LikelyandOncogenic_KEAP1OL_nomultihit_rankown_ambiguousintegrated/Precision_Recall_Curve/mutmetadata_less_StrictKEAP1OL_rankedbyme_ambiguousintegrated.tsv")

mutmetadata$mutornot = FALSE
mutmetadata$mutornot[mutmetadata$mutlabels!= "None" ] <- TRUE
resultsandmeta = inner_join(totalresults, mutmetadata)















totalresults = data.frame(fread("C:/Users/Raul/Documents/NRF2/expression_classifier/cpmlog2ssGSEA_allcancers.tsv"))
extrainfo= unique(as.data.frame(fread("C:/Users/Raul/Documents/NRF2/GSEA_analysis//extrainforTCGA.tsv"))[,1:3])
totalresults = inner_join(totalresults,extrainfo)
#totalresults$score_ssGSEA = totalresults$score_ssGSEA/totalresults$nozeroes * 100
totalresults = totalresults[,1:3]
totalresults = unique(totalresults)
#totalresults= totalresults[totalresults$TCGA_ID %in% sharedsamples,]
#totalresults = totalresults[!duplicated(totalresults$TCGA_ID),]
mutmetadata  =fread("C:/Users/Raul/Documents/NRF2/SJclassifierversions/mutationburden2000/TCGA_mut_metadata_mutburden2023_2000_exskipupdated.tsv")
mutmetadata  =fread("C:/Users/Raul/Documents/NRF2/automaticpipeline/automatic_versions/LikelyandOncogenic_KEAP1OL_nomultihit_rankown_ambiguousintegrated/Precision_Recall_Curve/mutmetadata_less_StrictKEAP1OL_rankedbyme_ambiguousintegrated.tsv")
mutmetadata  =fread("C:/Users/Raul/Documents/papers/SJ_paper/Ambiguous_082023/figures/Final_mutmetadata_ambiguousintegrated.tsv")

mutmetadata$mutornot = FALSE
mutmetadata$mutornot[mutmetadata$mutlabels!= "None" ] <- TRUE
resultsandmeta = inner_join(totalresults, mutmetadata)
resultsandmeta = resultsandmeta[resultsandmeta$TCGA_ID %in% sharedsamples,]
ssGSEA_PRAUC = c()
pr <- pr.curve(scores.class0 = resultsandmeta$score_ssGSEA, weights.class0 = resultsandmeta$mutornot)
ssGSEA_PRAUC = c(ssGSEA_PRAUC,(pr$auc.integral))

for(cont in 1:length(cancers)){
  pr <- pr.curve(scores.class0 = resultsandmeta$score_ssGSEA[resultsandmeta$cancer_label == cancers[cont]], weights.class0 = resultsandmeta$mutornot[resultsandmeta$cancer_label == cancers[cont]])
  ssGSEA_PRAUC = c(ssGSEA_PRAUC,(pr$auc.integral))
}
thresholds = seq(min(resultsandmeta$score_ssGSEA),max(resultsandmeta$score_ssGSEA), length.out =100) 
TP = c()
FP = c()
TN = c()
FN = c()
for(nthres in 1:length(thresholds)) {
  resultsandmeta$ssGSEAmutornot =   resultsandmeta$score_ssGSEA >= thresholds[nthres]
  TP = c(TP,  sum(resultsandmeta$mutornot == TRUE & resultsandmeta$ssGSEAmutornot == TRUE))
  TN = c(TN,  sum(resultsandmeta$mutornot == FALSE & resultsandmeta$ssGSEAmutornot == FALSE))
  FN = c(FN,  sum(resultsandmeta$mutornot == TRUE & resultsandmeta$ssGSEAmutornot == FALSE))
  FP = c(FP,  sum(resultsandmeta$mutornot == FALSE & resultsandmeta$ssGSEAmutornot == TRUE))
}
results = data.frame( 
  FALSE_not_called = TN,
  FALSE_called = FP,
  TRUE_not_called = FN,
  TRUE_called = TP,
  sensitivity = TP/(TP + FN),
  specificity = TN/(TN + FP),
  Precision = TP/(TP + FP),
  Recall = TP/(TP + FN),
  FPR =   1-TN/(TN + FP),
  threshold = thresholds, 
  cancer  ="PANCANCER", 
  method =  "ssGSEA"
  
)

allresults  = results
ggplot(results, aes(x = FPR,y = sensitivity)) + geom_line()
for(ncancers in 1:length(cancers)){
  resultsandmetabycancer  =resultsandmeta[resultsandmeta$Cancer_Type == cancers[ncancers],]
  
  
  thresholds = seq(min(resultsandmetabycancer$score_ssGSEA),max(resultsandmetabycancer$score_ssGSEA), length.out =100) 
  TP = c()
  FP = c()
  TN = c()
  FN = c()
  for(nthres in 1:length(thresholds)) {
    resultsandmetabycancer$ssGSEAmutornot =   resultsandmetabycancer$score_ssGSEA>= thresholds[nthres]
    TP = c(TP,  sum(resultsandmetabycancer$mutornot == TRUE & resultsandmetabycancer$ssGSEAmutornot == TRUE))
    TN = c(TN,  sum(resultsandmetabycancer$mutornot == FALSE & resultsandmetabycancer$ssGSEAmutornot == FALSE))
    FN = c(FN,  sum(resultsandmetabycancer$mutornot == TRUE & resultsandmetabycancer$ssGSEAmutornot == FALSE))
    FP = c(FP,  sum(resultsandmetabycancer$mutornot == FALSE & resultsandmetabycancer$ssGSEAmutornot == TRUE))
  }
  results = data.frame( 
    FALSE_not_called = TN,
    FALSE_called = FP,
    TRUE_not_called = FN,
    TRUE_called = TP,
    sensitivity = TP/(TP + FN),
    specificity = TN/(TN + FP),
    Precision = TP/(TP + FP),
    Recall = TP/(TP + FN),
    FPR =   1-TN/(TN + FP),
    threshold = thresholds, 
    cancer  =cancers[ncancers], 
    method =  "ssGSEA"
  )
  
  allresults = rbind(allresults, results)
  
  allresultswithpvalue = rbind(allresultswithpvalue,allresults)
  
}




treholdchosenbycancertype = c()
for(ncancers in 1){#:length(unique(allresults$cancer))){
  a = allresults[allresults$cancer == unique(allresults$cancer)[ncancers],]
  a$threshold = round(a$threshold,2)
  fmeasure =  (2 * a$Precision * a$Recall) / (a$Precision + a$Recall)
  # treholdchosenbycancertype = c(treholdchosenbycancertype,  allresults$threshold [which.max(fmeasure)])
  # fmeasures = c(fmeasures , sqrt(a$sensitivity * a$specificity)  )
  valuechosen=  unique(allresults$threshold [which.max(fmeasure)])
  # valuechosen = pancancerthres
  #valuechosen = 10
  #  itscoords = a[a$threshold == (round(valuechosen,digits = 2)),1:2]
  itscoords = a[a$threshold == (round(valuechosen,digits = 2)),]
  if(length(valuechosen) == 0  ){
    valuechosen = NA
    itscoords = c(NA,NA)
  } else if (is.na(valuechosen)){
    valuechosen = NA
    itscoords = c(NA,NA)
  }
  treholdchosenbycancertype =   rbind(treholdchosenbycancertype, (c(unique(allresults$cancer)[ncancers],valuechosen,itscoords )))
  
}



alltreholdchosenbycancertype = rbind(alltreholdchosenbycancertype, data.frame(treholdchosenbycancertype,Method = paste0("ssGSEA") ))





allresultswithpvaluePANCANCER = allresultswithpvalue[allresultswithpvalue$cancer == "PANCANCER",]

allresultswithpvaluePANCANCER$method =factor( allresultswithpvaluePANCANCER$method, levels = unique(allresultswithpvaluePANCANCER$method ))

pdf(paste0("C:/Users/Raul/Documents/NRF2/GSEA_analysis/PrecisionRecallcurve_SJvsGSVA_PANCANCER_TOPS_fmeasurethres_ssGSEAincluded_only_gathered_noGSVA_04052000_bySJfrequency.pdf"), width  = 10)
p =  ggplot(allresultswithpvaluePANCANCER, aes(x = Recall,y = Precision, color = method, alpha = 0.7)) + geom_line(linewidth = 2, alpha = 0.7)+                              
  scale_colour_brewer(palette = "Set1") + xlim(c(0,1)) + ylim(c(0,1)) + theme_minimal()
print(p)

dev.off()


#pdf(paste0("C:/Users/Raul/Documents/NRF2/automaticpipeline/automatic_versions/LikelyandOncogenic_KEAP1OL_nomultihit_rankown_ambiguousintegrated/Precision_Recall_Curve/PR_AUC_bycancer.pdf"))
cancers = unique(allresultswithpvalue$cancer)
cont  =1
#for(cont in 1:length(cancers)){
allresultswithpvaluePANCANCER = allresultswithpvalue[allresultswithpvalue$cancer == cancers[cont],]
allresultswithpvaluePANCANCER$method =factor( allresultswithpvaluePANCANCER$method, levels = unique(allresultswithpvaluePANCANCER$method ))
allresultswithpvaluePANCANCER$Method =factor( allresultswithpvaluePANCANCER$method, levels = unique(allresultswithpvaluePANCANCER$method ))

p2 =  ggplot(allresultswithpvaluePANCANCER, aes(x = Recall,y = Precision, color = Method)) + ggtitle(cancers[cont])+                              
  scale_colour_brewer(palette = "Set1")+ xlim(c(0,1)) + ylim(c(0,1)) + theme_minimal() + ggtitle(paste0("Precision Recall AUC "))  + 
  geom_text(family="Helvetica", aes(x = 0.98, y = 0.95, label = paste0( "SJ classifier AUC = ", round(SJ_PRAUC[cont],3),"\n ssGSEA AUC = ", round(ssGSEA_PRAUC[cont],3))),show.legend = FALSE, color = "black", size  =2.5, hjust = 1) +
  coord_fixed(ratio = 1) + geom_line(linewidth = 1, alpha = 0.7)+ 
my_theme() +   theme(legend.position = "bottom", legend.text = element_text(size = 7)) +   guides(color = guide_legend(nrow = 2)) # Adjust the font size for the text within the legend



p2

pdf("C:/Users/Raul/Documents/papers/SJ_paper/Ambiguous_082023/figures/fig3SJincluded/Fig3c.pdf", width  = letterincheswidth * 3/12, height =letterinchesheight/2 )
p2 + labs( tag = "c")
dev.off() 


pdf("C:/Users/Raul/Documents/papers/SJ_paper/Ambiguous_082023/figures/fig3SJincluded/Fig3c_notag.pdf", width  = letterincheswidth* 3/12, height =letterinchesheight/2 )
p2 
dev.off() 

png("C:/Users/Raul/Documents/papers/SJ_paper/Ambiguous_082023/figures/fig3SJincluded/Fig3c.png", width  = letterincheswidth* 3/12, height =letterinchesheight/2,units = "in", res = 500 )
p2 + labs( tag = "c")
dev.off() 


png("C:/Users/Raul/Documents/papers/SJ_paper/Ambiguous_082023/figures/fig3SJincluded/Fig3c_notag.png", width  = letterincheswidth* 3/12, height =letterinchesheight/2,units = "in", res = 500 )
p2
dev.off() 


  ######
  
  print(p2)
  # ggsave(paste0("C:/Users/Raul/Documents/NRF2/automaticpipeline/automatic_versions/LikelyandOncogenic_KEAP1OL_nomultihit_rankown_ambiguousintegrated/Precision_Recall_Curve/PR_AUC.pdf"), p)
#}
#dev.off()
  p4 = ggplot() + labs( tag = "D") +  theme_bw(base_family = "Helvetica") %+replace%
    theme(title = element_text(size=16),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "white", size = 0,5),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 20),
          legend.key = element_blank(),
          legend.key.size = unit (2.25, "cm"),
          legend.margin = margin (0.5,0.5,0.5,0.5),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 18),
          strip.text = element_text(size = 6),
          strip.background = element_blank(),
          complete = TRUE,
          plot.tag = element_text(size = 50))
  
  
  p4 = ggdraw() +
    draw_image("C:/Users/Raul/Documents/papers/SJ_paper/Ambiguous_082023/figures/igv_snapshot_ptrimmed.png") +
     theme_bw(base_family = "Helvetica") %+replace%
    theme(title = element_text(size=7),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "white", size = 0,5),
          axis.text = element_text(size = 7),
          axis.title = element_text(size = 7),
          legend.key = element_blank(),
          legend.margin = margin (0.5,0.5,0.5,0.5),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          strip.text = element_text(size = 7),
          strip.background = element_blank(),
          complete = TRUE,
          plot.tag = element_text(size = 7)) + ggtitle("TCGA-22-5486-01A: *KEAP1* Region") + theme(plot.title =  element_markdown())
  
  
  
  
  pdf("C:/Users/Raul/Documents/papers/SJ_paper/Ambiguous_082023/figures/fig3SJincluded/Fig3b.pdf", width  = letterincheswidth * 5/12, height =letterinchesheight/2 )
  p4+ labs( tag = "b")
  dev.off() 
  
  
  pdf("C:/Users/Raul/Documents/papers/SJ_paper/Ambiguous_082023/figures/fig3SJincluded/Fig3b_notag.pdf", width  = letterincheswidth * 5/12, height =letterinchesheight/2 )
  p4 
  dev.off() 
  
  png("C:/Users/Raul/Documents/papers/SJ_paper/Ambiguous_082023/figures/fig3SJincluded/Fig3b.png", width  = letterincheswidth * 5/12, height =letterinchesheight/2,units = "in", res = 1000 )
  p4 + labs( tag = "b")
  dev.off() 
  
  
  png("C:/Users/Raul/Documents/papers/SJ_paper/Ambiguous_082023/figures/fig3SJincluded/Fig3b_notag.png", width  = letterincheswidth * 5/12, height =letterinchesheight/2,units = "in", res = 1000 )
  p4
  dev.off() 
  
  
  
  
  
  
  
  
  
  
  
  
  ### tentative
  
  library(stringr)
  library(RColorBrewer)
  library(foreach)
  library(ggplot2)
  library(gridExtra)
  library(data.table)
  library(dplyr)
  library(beepr)
  library(parallel)
  metadatafile = "C:/Users/Raul/Documents/papers/SJ_paper/Ambiguous_082023/figures/mutmetadata_less_StrictKEAP1OL_rankedbyme_ambiguousintegrated.tsv"
  pthresholds =  "p0001"
  PANCANCERNFE2L2 = readRDS("C:/Users/Raul/Documents/papers/SJ_paper/Ambiguous_082023/figures/mc3.v0.2.8.PUBLIC_NRF2.maf.RDS")
  resultadosnumericPANCANCER =   data.frame(fread("C:/Users/Raul/Documents/papers/SJ_paper/Ambiguous_082023/figures/Scores_CV.tsv"))
  
  
  
  FN = resultadosnumericPANCANCER[resultadosnumericPANCANCER$mutlabels != "None" & resultadosnumericPANCANCER$score <10,]
  FN$mutlabels[FN$mutlabels =="NRF2_Exon2_GOF_Mut" ] <- "*NFE2L2* Exon 2 Hotspot"
  FN$mutlabels[FN$mutlabels =="KEAP1_Mut"  ] <- "*KEAP1*"
  FN$mutlabels[FN$mutlabels =="Exon2_Skipping" ] <-"*NFE2L2* Exon 2 Skip"
  
  dfFN = data.frame(table(FN$mutlabels, FN$cancer))
  dfFN$`Alteration Reported` = dfFN$Var1 
  
  dfFN$`Alteration Reported` = factor(dfFN$`Alteration Reported`, levels = levels(dfFN$`Alteration Reported`)[c(2,1,3)])
  p3 = ggplot(resultadosnumericPANCANCERforplot1, aes(x = `Mutation Labels` , y =  Score, fill = `Mutation Labels`)) + 
    geom_boxplot( color = "black",lwd=0.1,outlier.size = 0.05)+ scale_fill_manual(values = myColors, drop = FALSE) +
    scale_color_manual(values = myColors, drop = FALSE)+ theme_minimal()  + xlab("")+
    geom_hline(yintercept=10, linetype="dashed", color = "red") +  scale_x_discrete(labels = cursiva)+
    # theme(text = element_text(size = 25))   + 
    my_theme() +
    theme(legend.position = "none") + ggtitle("CV Results Classifier")  + theme(axis.text.x = element_markdown())
  
  
  tp1 = ggplot(dfFN, aes(x = Var2, y = Freq, fill = `Alteration Reported`)) +
    geom_bar(stat = "identity", position=position_dodge()) + my_theme()  +
    xlab("") + ylab("Mutated Samples not reaching threshold") + theme(legend.position = "bottom") +
    theme(legend.text = element_markdown())+
    scale_fill_manual(values = myColors[-1], drop = FALSE)
  tp1
  
  pdf("C:/Users/Raul/Documents/papers/SJ_paper/Ambiguous_082023/figures/fig3SJincluded/Temptativefig3SJincludedc_1.pdf", width  = letterincheswidth , height =letterinchesheight/2 )
  tp1+ labs( tag = "X")
  dev.off() 
  
  
  
  FN$`Alteration Reported` = factor(FN$mutlabels, levels = (unique(FN$mutlabels)))
 tp2 =  ggplot(FN, aes(x = cancer, y = score, fill = `Alteration Reported`)) + geom_boxplot() +
    my_theme()  +
    xlab("") + ylab("Score") + theme(legend.position = "bottom") +
    theme(legend.text = element_markdown())+
    scale_fill_manual(values = myColors[-1], drop = FALSE)

#  
# 
# pdf("C:/Users/Raul/Documents/papers/SJ_paper/Ambiguous_082023/figures/fig3SJincluded/Temptativefig3SJincludedc_2.pdf", width  = letterincheswidth , height =letterinchesheight/2 )
# tp2+ labs( tag = "X")
# dev.off() 
# 
#  
#   
#  ggplot(FN, aes(x = cancer, y = score, fill = mutlabels)) + geom_boxplot() + my_theme()
#  
#  
#  FP = resultadosnumericPANCANCER[resultadosnumericPANCANCER$mutlabels == "None" & resultadosnumericPANCANCER$score >= 10,]
#  
#  
#  dfFP = data.frame(table(FP$mutlabels, FP$cancer))
#  ggplot(dfFP, aes(x = Var2, y = Freq)) + geom_bar(stat = "identity", position=position_dodge()) + my_theme()
#  
#  tp3 = ggplot(dfFP, aes(x = Var2, y = Freq)) +
#    geom_bar(stat = "identity") + my_theme()  +
#    xlab("") + ylab("Not mutated samples reaching threshold")+
#    theme(legend.text = element_markdown())
#  tp3
#  
#  pdf("C:/Users/Raul/Documents/papers/SJ_paper/Ambiguous_082023/figures/fig3SJincluded/Temptativefig3SJincludedc_3.pdf", width  = letterincheswidth , height =letterinchesheight/2 )
#  tp3+ labs( tag = "X")
#  dev.off() 
#  
#  
#  
#  tp4 =  ggplot(FP, aes(x = cancer, y = score)) + geom_boxplot() +
#    my_theme()  +
#    xlab("") + ylab("Score") +
#    theme(legend.text = element_markdown())+
#    scale_fill_manual(values = myColors[-1], drop = FALSE)
#  
#  tp4
#  pdf("C:/Users/Raul/Documents/papers/SJ_paper/Ambiguous_082023/figures/fig3SJincluded/Temptativefig3SJincludedc_4.pdf", width  = letterincheswidth , height =letterinchesheight/2 )
#  tp4+ labs( tag = "X")
#  dev.off() 
#  
#  
#  
#  FN %>% 
#    group_by(cancer) %>%
#    summarize(Median = median(score), Mean = mean(score))
#  ddply(FN, .(cancer), summarize, mean=mean(score))
#  
#  
#  
#  
#  
#  TP = resultadosnumericPANCANCER[resultadosnumericPANCANCER$mutlabels != "None" & resultadosnumericPANCANCER$score >=10,]
#  
#  df = data.frame(resultadosnumericPANCANCER  %>% 
#                    group_by(cancer) %>% 
#                    summarize(TP =sum(score >= 10 & mutlabels != "None"),
#                              FP =sum(score >= 10 & mutlabels == "None"),
#                              FN =sum(score < 10 & mutlabels != "None"),
#                              TN =sum(score < 10 & mutlabels == "None")))
#  
#  dff = df[rowSums(df[,2:4]) != 0,]
#  dfm  =reshape2::melt(dff)
#  
#  dfm$Outcome = dfm$variable
#  tp5 = ggplot(dfm[dfm$variable != "TN",], aes(x  =cancer, y = value, fill = Outcome)) + geom_bar(stat = "identity", position=position_dodge()) + my_theme()+
#    geom_bar(stat = "identity", position=position_dodge()) + my_theme()  +
#    xlab("") + ylab("Count") + theme(legend.position = "bottom") +
#    theme(legend.text = element_markdown()) +   theme(legend.title = element_blank())
#  tp5
#  
#  pdf("C:/Users/Raul/Documents/papers/SJ_paper/Ambiguous_082023/figures/fig3SJincluded/Temptativefig3SJincludedc_5.pdf", width  = letterincheswidth , height =letterinchesheight/2 )
#  tp5+ labs( tag = "X")
#  dev.off() 
#  
#  
#  
#  #