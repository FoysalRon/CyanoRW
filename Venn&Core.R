getwd()
####Load Library####
library("microbiomeSeq")
library("phyloseq")
library("ggplot2")
library("ggtree")
library("microbiome")
library("hrbrthemes")
library("gcookbook")
library("tidyverse")
library("dplyr")
library("coin")
library("VennDiagram")
library("UpSetR")
library("patchwork")
library("RColorBrewer")
library("VennDiagram")
library("ggpattern")
library('cowplot')
library('ggpubr')
library('vegan')
library('microeco')
library('file2meco')
library("fantaxtic")
library("gridExtra")
library("knitr")
library("magrittr")
library("ggnested")
library("stringr")

####Set colour####
ET_c <- c ("#FFE260", "#ffa600", "#C74B1D", "#7C020F", "#f95d6a", "#d45087", "#a05195", "#665191", "lightskyblue1", "steelblue3", "#2f4b7c",  "#AEC248", "#666666", "#33A02C",  "#2F8006", "#3C5E30", "#113A5C", "#06171E", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")
ET_co <- c("steelblue3", "lightsalmon3", "gold3", "aquamarine4")
ET_col <- c("firebrick", "forestgreen", "royalblue", "chocolate", "gold3")

####Load data (all data)####
OTU<-read.csv("asvtableET.csv", header=TRUE, row.names=1)
TAX<-read.csv("taxonomyET.csv", header=TRUE, row.names=1)
META<-read.csv("metadataET.csv", header=TRUE, row.names=1)
class(OTU)
class(TAX)
OTU_data=otu_table(as.matrix(OTU), taxa_are_rows=TRUE)
TAX_data=tax_table(as.matrix(TAX))
META_data=sample_data(META)
physeq<-merge_phyloseq(phyloseq(OTU_data, TAX_data), META_data)
pseq<-t(physeq)

library(MicrobiotaProcess)
#Rarefy the data
psraw <- prune_samples(sample_sums(physeq)>=sort(rowSums(otu_table(physeq)))[3], physeq)
set.seed(1024)
ps <- rarefy_even_depth(psraw)
ps

####Venn diagram (source)####
library(MicEco)
et_vn <- ps_venn(
  ps,
  group="Source",
  fraction = 0,
  weight = FALSE,
  relative = TRUE,
  plot = TRUE,
  labels = list(cex = 0.75),
  size=16,
  hjust=1,
  fill=c("steelblue3", "salmon3"))
et_vn

####Venn diagram (Location)####
et_vn1 <- ps_venn(
  ps,
  group="Location",
  fraction = 0,
  weight = FALSE,
  relative = TRUE,
  plot = TRUE,
  labels = list(cex = 0.75),
  size=16,
  hjust=1,
  fill=c("steelblue3", "salmon3", "gold3", "aquamarine4"))
et_vn1

####Venn diagram (Season)####
et_vn2 <- ps_venn(
  ps,
  group="Season",
  fraction = 0,
  weight = FALSE,
  relative = TRUE,
  plot = TRUE,
  labels = list(cex = 0.75),
  size=16,
  hjust=1,
  fill=c("steelblue3", "salmon3", "gold3", "aquamarine4"))
et_vn2
####Core taxa####
library('ggvenn')
##Core taxa in location (Venn)
table(meta(ps)$Location, useNA = "always")
pseq.rel <- microbiome::transform(ps, "compositional")
ET <- unique(as.character(meta(ps)$Location))
print(ET)
list_core <- c() # an empty object to store information
for (n in ET){ # for each variable n in location
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(ps, Location == n) # Choose sample from location by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, # 0.001 in atleast 90% samples 
                         prevalence = 0.75)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each location.
  list_core[[n]] <- core_m # add to a list core taxa for each location.
  #print(list_core)
}
print(list_core)
p_venn <- ggvenn(list_core,
                 fill_color = c("steelblue3", "orangered4", "gold3", "aquamarine4", "chocolate"),
                 stroke_size = 1.0, set_name_size = 4, text_size = 5)

##Core taxa in source (Venn)
table(meta(ps)$Source, useNA = "always")
pseq.rel <- microbiome::transform(ps, "compositional")
ET1 <- unique(as.character(meta(ps)$Source))
print(ET1)
list_core <- c() # an empty object to store information
for (n in ET1){ # for each variable n in source
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(ps, Source == n) # Choose sample from source by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, # 0.001 in atleast 90% samples 
                         prevalence = 0.75)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each source.
  list_core[[n]] <- core_m # add to a list core taxa for each source.
  #print(list_core)
}
print(list_core)
p_venn1 <- ggvenn(list_core,
                  fill_color = c("steelblue3", "orangered4", "gold3", "aquamarine4", "chocolate"),
                  stroke_size = 1.0, set_name_size = 4, text_size = 5)


##Core taxa in location (Venn)
table(meta(ps)$Location, useNA = "always")
pseq.rel <- microbiome::transform(ps, "compositional")
ET <- unique(as.character(meta(ps)$Location))
print(ET)
list_core <- c() # an empty object to store information
for (n in ET){ # for each variable n in location
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(ps, Location == n) # Choose sample from location by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, # 0.001 in atleast 90% samples 
                         prevalence = 0.75)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each location.
  list_core[[n]] <- core_m # add to a list core taxa for each location.
  #print(list_core)
}
print(list_core)
p_venn <- ggvenn(list_core,
                 fill_color = c("steelblue3", "orangered4", "gold3", "aquamarine4", "chocolate"),
                 stroke_size = 1.0, set_name_size = 4, text_size = 5)


##Core taxa in season (Venn)
table(meta(ps)$Season, useNA = "always")
pseq.rel <- microbiome::transform(ps, "compositional")
ET2 <- unique(as.character(meta(ps)$Season))
print(ET2)
list_core <- c() # an empty object to store information
for (n in ET2){ # for each variable n in season
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(ps, Season == n) # Choose sample from season by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, # 0.001 in atleast 90% samples 
                         prevalence = 0.75)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each season.
  list_core[[n]] <- core_m # add to a list core taxa for each season.
  #print(list_core)
}
print(list_core)
p_venn2 <- ggvenn(list_core,
                  fill_color = c("steelblue3", "orangered4", "gold3", "aquamarine4", "chocolate"),
                  stroke_size = 1.0, set_name_size = 4, text_size = 5)
p_venn2

##Horizontal alignment##
top_row <- plot_grid(et_vn, et_vn1, et_vn2, labels = "AUTO", rel_widths = c(2, 2, 2), rel_heights = c(2, 2, 2), ncol = 3, align = "h")
top_row
bottom_row <- plot_grid(p_venn1, p_venn, p_venn2, labels = c('D', 'E', 'F'), rel_widths = c(2, 2, 2), rel_heights = c(2, 2, 2), ncol = 3, align = "h")
bottom_row
plot_grid(top_row, bottom_row, nrow = 2, rel_heights = c(1.5, 2))

##Vertical alignment##
top_row <- plot_grid(et_vn, et_vn1, et_vn2, labels = "AUTO", rel_widths = c(2, 2, 2), rel_heights = c(2, 2, 2), ncol = 1, align = "v")
top_row
bottom_row <- plot_grid(p_venn1, p_venn, p_venn2, labels = c('D', 'E', 'F'), rel_widths = c(2, 2, 2), rel_heights = c(2, 2, 2), ncol = 1, align = "v")
bottom_row
plot_grid(top_row, bottom_row, ncol = 2, rel_heights = c(2, 2))

##Save plot
ggsave("venn_location_season_sour.tiff", units = c("in"), width=12, height=12, dpi=300, compression="lzw")

####Core-genera heat map####
#Load data (all data)
OTU<-read.csv("asvtableET.csv", header=TRUE, row.names=1)
TAX<-read.csv("taxonomyET.csv", header=TRUE, row.names=1)
META<-read.csv("metadataET.csv", header=TRUE, row.names=1)
class(OTU)
class(TAX)
OTU_data=otu_table(as.matrix(OTU), taxa_are_rows=TRUE)
TAX_data=tax_table(as.matrix(TAX))
META_data=sample_data(META)
physeq<-merge_phyloseq(phyloseq(OTU_data, TAX_data), META_data)
pseq<-t(physeq)

library(MicrobiotaProcess)
# Rarefy the data
psraw <- prune_samples(sample_sums(physeq)>=sort(rowSums(otu_table(physeq)))[3], physeq)
set.seed(1024)
ps <- rarefy_even_depth(psraw)
ps
# sub-set meta-data
table(meta(ps)$Source)
table(meta(ps)$Source, meta(ps)$Season)
# Filter the data to include only gut samples from M3 subject
ps.m3 <- subset_samples(ps, Source == "Water") 
print(ps.m3)
ps.m3 <- prune_taxa(taxa_sums(ps.m3) > 0, ps.m3)
print(ps.m3)
# Calculate compositional version of the data
# (relative abundances)
ps.m3.rel <- microbiome::transform(ps.m3, "compositional")
taxa_names(ps.m3.rel)[1:2]
ps.m3.rel <- microbiome::add_refseq(ps.m3.rel)
# Check if ref_seq slot is added to phyloseq object
print(ps.m3.rel)
taxa_names(ps.m3.rel)[1:3]
core.taxa.standard <- core_members(ps.m3.rel, detection = 0.0001, prevalence = 50/100)
core.taxa.standard

# first combine genus and species names. 
tax_table(ps.m3.rel)[tax_table(ps.m3.rel) == "Unclassified"] <- NA
tax_table(ps.m3.rel)[tax_table(ps.m3.rel) == "Unclassified"] <- NA
tax_table(ps.m3.rel)[tax_table(ps.m3.rel) == "Unclassified"] <- NA
tax_table(ps.m3.rel)[tax_table(ps.m3.rel) == "Unclassified"] <- NA
tax_table(ps.m3.rel)[tax_table(ps.m3.rel) == "Unclassified"] <- NA
tax_table(ps.m3.rel)[tax_table(ps.m3.rel) == "Unclassified"] <- NA
tax_table(ps.m3.rel)[tax_table(ps.m3.rel) == "Unclassified"] <- NA
tax_table(ps.m3.rel)[, colnames(tax_table(ps.m3.rel))] <- gsub(tax_table(ps.m3.rel)[, colnames(tax_table(ps.m3.rel))],  pattern = "[a-z]Unclassified", replacement = "")

ps.m3.rel.gen <- aggregate_taxa(ps.m3.rel, "Genus")
any(taxa_names(ps.m3.rel.gen) == "Unknown")
# Remove Unknown
ps.m3.rel.gen <- subset_taxa(ps.m3.rel.gen, Genus!="Unknown")

library(RColorBrewer)
prevalences <- seq(.05, 1, .05)
detections <- round(10^seq(log10(1e-5), log10(.2), length = 10), 3)

#Remove big unwanted names
taxa_names(ps.m3.rel.gen) <- gsub("Bacteria_Cyanobacteria_Cyanobacteriia_Chloroplast_Chloroplast_Z_",
                                  "", taxa_names(ps.m3.rel.gen))

p1g <- plot_core(ps.m3.rel.gen, 
                 plot.type = "heatmap", 
                 colours = rev(brewer.pal(5, "RdBu")),
                 prevalences = prevalences, 
                 detections = detections, min.prevalence = .9) +
  xlab("Detection Threshold (Relative Abundance (%))") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90),
        axis.text.y = element_text(size = 12, face = "italic"))
p1g

##Save plot
ggsave("core_genera_all_samples_source_season_source.tiff", units = c("in"), width=9, height=8, dpi=300, compression="lzw")