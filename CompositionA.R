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

#Load data (all data)
OTU<-read.csv("asvtableETCW.csv", header=TRUE, row.names=1)
TAX<-read.csv("taxonomyETC.csv", header=TRUE, row.names=1)
META<-read.csv("metadataETW.csv", header=TRUE, row.names=1)
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

####Composition - phyla####
# Create a data table for ggplot
ps_phylum <- ps %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance (or use ps0.ra)
  psmelt() %>%                                         # Melt to long format for easy ggploting
  filter(Abundance > 0.01)                             # Filter out low abundance taxa
head(ps_phylum)
dat.agr = aggregate(Abundance~Source+Location+Phylum, data=ps_phylum, FUN=mean)
gp1 <- ggplot(dat.agr, aes(x = Source, y = Abundance, fill = Phylum)) +
  facet_grid(cols = vars(Location)) +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing =
          unit(0.1, "lines")) +
  theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=10, face="bold"), axis.text.y =
          element_text(hjust = 1, size=12), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=0.5, size=12, face="bold", margin = margin(t = 6, r = 4, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.4, size=12, face="bold", margin = margin(t = 0, r = 6, b = 0, l = 0)), plot.title
        =element_text(angle = 0, hjust = 0, size=10, face="bold", margin = margin(t = 0, r = 0, b = 8, l = 0)) ) +
  geom_bar(stat = "identity", position="fill", colour = "black") +
  scale_fill_manual(values =colors) +
  theme(legend.text=element_blank()) + theme(axis.text.x=element_markdown()) + scale_y_continuous(label = scales::percent) +
  ylab("Relative Abundance") +
  ggtitle ("Top phyla")
gp1

dat.agr = aggregate(Abundance~Source+Season+Phylum, data=ps_phylum, FUN=mean)
gp2 <- ggplot(dat.agr, aes(x = Source, y = Abundance, fill = Phylum)) +
  facet_grid(cols = vars(Season)) +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing =
          unit(0.1, "lines")) +
  theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=10, face="bold"), axis.text.y =
          element_text(hjust = 1, size=12), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=0.5, size=12, face="bold", margin = margin(t = 6, r = 4, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.4, size=12, face="bold", margin = margin(t = 0, r = 6, b = 0, l = 0)), plot.title
        =element_text(angle = 0, hjust = 0, size=10, face="bold", margin = margin(t = 0, r = 0, b = 8, l = 0)) ) +
  geom_bar(stat = "identity", position="fill", colour = "black") +
  scale_fill_manual(values =colors) +
  theme(legend.text=element_text(size=12, face = "italic"), legend.title = element_text(size=12, face="bold")) + theme(axis.text.x=element_markdown(), legend.text=element_text(face="italic"), legend.key.size = unit(10, "pt")) + scale_y_continuous(label = scales::percent) +
  ylab("Relative Abundance") +
  ggtitle ("Top phyla")
gp2

#Arrange plots
plot_grid(gp1, gp2, nrow = 1, labels = "AUTO", rel_widths = c(1.06, 2.0))
#Save plot
ggsave("phyla_location_season_source_facet.tiff", units = c("in"), width=11, height=7, dpi=300, compression="lzw")


###Top 15 Cyano genera###
#Load data (all data)
OTU<-read.csv("asvtableETCc.csv", header=TRUE, row.names=1)
TAX<-read.csv("taxonomyETCc.csv", header=TRUE, row.names=1)
META<-read.csv("metadataET.csv", header=TRUE, row.names=1)
class(OTU)
class(TAX)
OTU_data=otu_table(as.matrix(OTU), taxa_are_rows=TRUE)
TAX_data=tax_table(as.matrix(TAX))
META_data=sample_data(META)
physeq<-merge_phyloseq(phyloseq(OTU_data, TAX_data), META_data)
pseq<-t(physeq)

#Rarefy data
library(MicrobiotaProcess)
psraw <- prune_samples(sample_sums(physeq)>=sort(rowSums(otu_table(physeq)))[3], physeq)
set.seed(1024)
ps <- rarefy_even_depth(psraw)
ps

#Transform data for meco
ps.rel = transform_sample_counts(ps, function(x) x/sum(x)*100)
meco_dataset <- phyloseq2meco(ps)

ET_cc <- c ("#666666", "#FF7F00", "#CAB2D6", "#6A3D9A", "#B15928", "#FFFF99", "#AEC248", "#33A02C",  "#2F8006", "#3C5E30",  "#FB9A99", "orangered", "gold3", "#113A5C", "#06171E", "#33A02C", "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

#Top 15 genera abundance in pie for the source
t1 <- trans_abund$new(dataset = meco_dataset, taxrank = "Genus", ntaxa = 15, groupmean = "Source")
#smaller font
p1 <- t1$plot_pie(facet_nrow = 2) + scale_fill_manual(values = ET_c) + theme(legend.text=element_text(size=12, face = "italic"), legend.title = element_text(size=12, face="bold")) + theme(legend.key.size = unit(0.5, "cm"))
p1

#Top 15 genera abundance in pie for the location
t2 <- trans_abund$new(dataset = meco_dataset, taxrank = "Genus", ntaxa = 15, groupmean = "Location")
#smaller font
p2 <- t2$plot_pie(facet_nrow = 2) + scale_fill_manual(values = ET_c) + theme(legend.text=element_text(size=12, face = "italic"), legend.title = element_text(size=12, face="bold")) + theme(legend.key.size = unit(0.5, "cm"))
p2

#Top 15 genera abundance in pie for the season
t3 <- trans_abund$new(dataset = meco_dataset, taxrank = "Genus", ntaxa = 15, groupmean = "Season")
#smaller font
p3 <- t3$plot_pie(facet_nrow = 2) + scale_fill_manual(values = ET_c) + theme(legend.text=element_text(size=12, face = "italic"), legend.title = element_text(size=12, face="bold")) + theme(legend.key.size = unit(0.5, "cm"))
p3

##Arrange plots
t_row <- plot_grid(p1, p2, labels = c('A', 'B'), rel_widths = c(1.8, 2), rel_heights = c(2, 2), nrow = 1, align = "h")
t_row
plot_grid(t_row, p3, nrow = 2, labels = c('', 'C'), rel_widths = c(2, 2))
##Save plot
ggsave("Cyano_com_location_season_source.tiff", units = c("in"), width=12, height=12, dpi=300, compression="lzw")

##Cyanobacteria in water
##Alpha-diversity (Cyanobacteria)
#Load data (all data)
OTU<-read.csv("asvtableETCcw.csv", header=TRUE, row.names=1)
TAX<-read.csv("taxonomyETCc.csv", header=TRUE, row.names=1)
META<-read.csv("metadataETW.csv", header=TRUE, row.names=1)
class(OTU)
class(TAX)
OTU_data=otu_table(as.matrix(OTU), taxa_are_rows=TRUE)
TAX_data=tax_table(as.matrix(TAX))
META_data=sample_data(META)
physeq<-merge_phyloseq(phyloseq(OTU_data, TAX_data), META_data)
pseq<-t(physeq)

##Line graph (microEco)
t1 <- trans_abund$new(dataset = meco_dataset, taxrank = "Phylum", ntaxa = 8, group = "Season")
t2w <- t1$plot_line(position = position_dodge(0.3), xtext_type_hor = TRUE)
t2w <- t2w + theme_bw() +
  theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing =
          unit(0.1, "lines")) +
  theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=10, face="bold"), axis.text.y =
          element_text(hjust = 1, size=12), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=0.5, size=12, face="bold"), axis.title.y = element_text(angle = 90, hjust = 0.4, size=12, face="bold"), plot.title = element_text(angle = 0, hjust = 0, size=10, face="bold", margin = margin(t = 0, r = 0, b = 8, l = 0)) ) +
  theme(legend.text=element_text(size=12, face = "bold"), legend.title = element_text(size=12, face="bold", colour = "black")) + theme(axis.text.x=element_markdown(), legend.text=element_text(face="bold"), legend.key.size = unit(10, "pt")) +
  ylab("Relative Abundance (%)") +
  xlab("Season") + ggtitle ("Seasonal changes")

##Save plot
ggsave("phyla_changes_season_source.tiff", units = c("in"), width=10, height=8, dpi=300, compression="lzw")

##Composition (Source)##
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
#Rarefy the data
psraw <- prune_samples(sample_sums(physeq)>=sort(rowSums(otu_table(physeq)))[3], physeq)
set.seed(1024)
ps <- rarefy_even_depth(psraw)
ps

#Transform data for MicroEco (plot to abundant phyla)
ps.rel = transform_sample_counts(ps, function(x) x/sum(x)*100)
meco_dataset <- phyloseq2meco(ps)
#Top 15 phyla abundance in pie for the season
t5 <- trans_abund$new(dataset = meco_dataset, taxrank = "Phylum", ntaxa = 15, groupmean = "Source")
#smaller font
w5 <- t5$plot_pie(facet_nrow = 2) + scale_fill_manual(values = ET_cc) + theme(legend.text=element_text(size=12, face = "italic"), legend.title = element_text(size=12, face="bold")) + theme(legend.key.size = unit(0.5, "cm"))
w5
##Save plot
ggsave("Cyano_sd_wa_season.tiff", units = c("in"), width=6, height=6.5, dpi=300, compression="lzw")

#Top abundant species
# Merges ASVs that have the same taxonomy rank (Genus)
gp = tax_glom(ps, taxrank = "Genus")
# Calculate relative abundance
rb = transform_sample_counts(gp, function(x) {x/sum(x)} )
top50 = names(head(sort(rowSums(otu_table(rb)), decreasing = TRUE), 50))
top50
rb = prune_taxa(top50, rb)
rb
rb %>% psmelt() %>% head(50)
rb2 = psmelt(rb)
rb2
dat.agr = aggregate(Abundance~Source+Genus, data=rb2, FUN=mean)

#Set colour (ET_c for Data_ET and ET_ca for ETR file data analysis)
ET_c <- c ("#FFE260", "#ffa600", "#C74B1D", "#7C020F", "#f95d6a", "#d45087", "#a05195", "#665191", "lightskyblue1", "steelblue3", "#2f4b7c",  "#AEC248", "#33A02C",  "#2F8006", "#3C5E30", "#113A5C", "#06171E", "#33A02C", "#FB9A99", "#FDBF6F", "grey80", "#666666", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")
ET_ca <- c ("#FFE260", "#ffa600", "#C74B1D", "#7C020F", "#f95d6a", "#d45087", "#a05195", "#665191", "lightskyblue1", "steelblue3", "#2f4b7c", "#AEC248", "#2F8006", "#3C5E30", "#113A5C", "#06171E", "#33A02C", "grey85", "#666666", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

g1 <- ggplot(dat.agr, aes(x = Source, y = Abundance, fill = Genus)) +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing =
          unit(0.1, "lines")) +
  theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=10, face="bold"), axis.text.y =
          element_text(hjust = 1, size=12), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=0.5, size=12, face="bold", margin = margin(t = 6, r = 4, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.4, size=12, face="bold", margin = margin(t = 0, r = 6, b = 0, l = 0)), plot.title
        =element_text(angle = 0, hjust = 0, size=10, face="bold", margin = margin(t = 0, r = 0, b = 8, l = 0)) ) +
  geom_bar(stat = "identity", position="fill", colour = "black") +
  scale_fill_manual(values =ET_c) +
  theme(legend.text=element_text(size=12, face = "italic"), legend.title = element_text(size=12, face="bold")) + theme(axis.text.x=element_markdown(), legend.text=element_text(face="italic"), legend.key.size = unit(10, "pt")) + scale_y_continuous(label = scales::percent) +
  ylab("Relative Abundance") +
  xlab("Source") + ggtitle ("Top genera")
g1

#Top genera (location)
dat.agr = aggregate(Abundance~Location+Genus, data=rb2, FUN=mean)
g2 <- ggplot(dat.agr, aes(x = Location, y = Abundance, fill = Genus)) +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing =
          unit(0.1, "lines")) +
  theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=10, face="bold"), axis.text.y =
          element_text(hjust = 1, size=12), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=0.5, size=12, face="bold", margin = margin(t = 6, r = 4, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.4, size=12, face="bold", margin = margin(t = 0, r = 6, b = 0, l = 0)), plot.title
        =element_text(angle = 0, hjust = 0, size=10, face="bold", margin = margin(t = 0, r = 0, b = 8, l = 0)) ) +
  geom_bar(stat = "identity", position="fill", colour = "black") +
  scale_fill_manual(values =ET_c) +
  theme(legend.text=element_text(size=12, face = "italic"), legend.title = element_text(size=12, face="bold")) + theme(axis.text.x=element_markdown(), legend.text=element_text(face="italic"), legend.key.size = unit(10, "pt")) + scale_y_continuous(label = scales::percent) +
  ylab("Relative Abundance") +
  xlab("Location") + ggtitle ("Top genera")
g2

#Top genera (Season)
dat.agr = aggregate(Abundance~Season+Genus, data=rb2, FUN=mean)
g3 <- ggplot(dat.agr, aes(x = Season, y = Abundance, fill = Genus)) +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing =
          unit(0.1, "lines")) +
  theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=10, face="bold"), axis.text.y =
          element_text(hjust = 1, size=12), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=0.5, size=12, face="bold", margin = margin(t = 6, r = 4, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.4, size=12, face="bold", margin = margin(t = 0, r = 6, b = 0, l = 0)), plot.title
        =element_text(angle = 0, hjust = 0, size=10, face="bold", margin = margin(t = 0, r = 0, b = 8, l = 0)) ) +
  geom_bar(stat = "identity", position="fill", colour = "black") +
  scale_fill_manual(values =ET_c) +
  theme(legend.text=element_text(size=12, face = "italic"), legend.title = element_text(size=12, face="bold")) + theme(axis.text.x=element_markdown(), legend.text=element_text(face="italic"), legend.key.size = unit(10, "pt")) + scale_y_continuous(label = scales::percent) +
  ylab("Relative Abundance") +
  xlab("Season") + ggtitle ("Top genera")
g3

#Arrange plots
t_row <- plot_grid(g1, g2, labels = c('A', 'B'), rel_widths = c(1.9, 2), rel_heights = c(2, 2), ncol = 2, align = "h")
m_row <- plot_grid(g3, NULL, labels = c('C', ''), rel_widths = c(1.8, 1.5), rel_heights = c(2, 2), ncol = 2, align = "h")
m_row
plot_grid(t_row, m_row, nrow = 2, rel_widths = c(2, 1.8))

##Save plot
ggsave("composition_location_season_source.tiff", units = c("in"), width=12, height=9, dpi=300, compression="lzw")

##Microbial composition (in water)
#Load data (all data)
OTU<-read.csv("asvtableETW.csv", header=TRUE, row.names=1)
TAX<-read.csv("taxonomyET.csv", header=TRUE, row.names=1)
META<-read.csv("metadataETW.csv", header=TRUE, row.names=1)
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
# Merges ASVs that have the same taxonomy rank (Genus)
gp = tax_glom(ps, taxrank = "Genus")
# Calculate relative abundance
rb = transform_sample_counts(gp, function(x) {x/sum(x)} )
top30 = names(head(sort(rowSums(otu_table(rb)), decreasing = TRUE), 30))
top30
rb = prune_taxa(top30, rb)
rb
rb %>% psmelt() %>% head(30)
rb2 = psmelt(rb)
rb2
dat.agr = aggregate(Abundance~Season+Genus, data=rb2, FUN=mean)
ET_ca <- c ("#FFE260", "#ffa600", "#C74B1D", "#7C020F", "#f95d6a", "#d45087", "#a05195", "#665191", "lightskyblue1", "steelblue3", "#2f4b7c", "#AEC248", "#2F8006", "#3C5E30", "#113A5C", "#06171E", "#33A02C", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "grey85", "#666666", "#6A3D9A", "#FFFF99", "#B15928", "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")
#Top genera (Season)
dat.agr = aggregate(Abundance~Season+Genus, data=rb2, FUN=mean)
g4 <- ggplot(dat.agr, aes(x = Season, y = Abundance, fill = Genus)) +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing =
          unit(0.1, "lines")) +
  theme(strip.text = element_text(colour = 'black'), axis.text.x = element_text(angle = 45, hjust = 1, size=10, face="bold"), axis.text.y =
          element_text(hjust = 1, size=12), axis.title.x = element_text(angle = 0, hjust = 0.5, vjust=0.5, size=12, face="bold", margin = margin(t = 6, r = 4, b = 0, l = 0)), axis.title.y = element_text(angle = 90, hjust = 0.4, size=12, face="bold", margin = margin(t = 0, r = 6, b = 0, l = 0)), plot.title
        =element_text(angle = 0, hjust = 0, size=10, face="bold", margin = margin(t = 0, r = 0, b = 8, l = 0)) ) +
  geom_bar(stat = "identity", position="fill", colour = "black") +
  scale_fill_manual(values =ET_ca) +
  theme(legend.text=element_text(size=12, face = "italic"), legend.title = element_text(size=12, face="bold")) + theme(axis.text.x=element_markdown(), legend.text=element_text(face="italic"), legend.key.size = unit(10, "pt")) + scale_y_continuous(label = scales::percent) +
  ylab("Relative Abundance") +
  xlab("Season") + ggtitle ("Top genera in water")
g4

##Save plot
ggsave("top_genera_wa_season2.tiff", units = c("in"), width=7, height=5.0, dpi=300, compression="lzw")
