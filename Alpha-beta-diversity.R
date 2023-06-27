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

####Alpha-diversity####
#Stat-compare-means (add p-value)
my_comparisons <- list( c("Sediment", "Water"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

#Alpha-diversity analysis (source)
(p = plot_richness(ps, x = "Source", measures=c("Observed", "Shannon")))
p5 <- p + geom_boxplot(data = p$data, aes(x = Source, y = value, fill = Source), alpha = 0.8, lwd =0.8, outlier.shape = NA) + geom_jitter(aes(fill = Source), alpha = 0.6, stroke = 1.0, shape = 21, size = 3.0, position = position_jitterdodge()) + theme_bw(base_size = 12, base_line_size = 1.5) + theme (axis.text.x = element_text(size=12, face = "bold", colour = "black", angle = 45, hjust=1)) + theme (axis.text.y = element_text(size = 12)) + theme(text = element_text(size = 12, face = "bold")) + scale_fill_manual(values = c("steelblue3", "orangered4")) + theme(legend.position = "none") + theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "transparent", color = NA), panel.grid.major = element_blank(), panel.grid.minor = element_line(colour = "white")) + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, label = "p.signif", symnum.args = symnum.args)
p5$layers <- p5$layers[-1]
p5

#Stat-compare-means (add p-value)
my_comparisons <- list( c("Inlet", "Bat"), c("Inlet", "Corner"), c("Bat", "Corner"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

#Alpha-diversity analysis (location)
#Alpha-diversity analysis
(p = plot_richness(ps, x = "Location", measures=c("Observed", "Shannon")))
p6 <- p + geom_boxplot(data = p$data, aes(x = Location, y = value, fill = Location), alpha = 0.8, lwd =0.8, outlier.shape = NA) + geom_jitter(aes(fill = Location), alpha = 0.6, stroke = 1.0, shape = 21, size = 3.0, position = position_jitterdodge()) + theme_bw(base_size = 12, base_line_size = 1.5) + theme (axis.text.x = element_text(size=12, face = "bold", colour = "black", angle = 45, hjust=1)) + theme (axis.text.y = element_text(size = 12)) + theme(text = element_text(size = 12, face = "bold")) + scale_fill_manual(values = c("steelblue3", "orangered4", "gold3")) + theme(legend.position = "none") + theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "transparent", color = NA), panel.grid.major = element_blank(), panel.grid.minor = element_line(colour = "white")) + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, label = "p.signif", symnum.args = symnum.args)
p6$layers <- p6$layers[-1]
p6

#Stat-compare-means (add p-value)
my_comparisons <- list( c("Summer", "Winter"), c("Summer", "Spring"), c("Summer", "Autumn"), c("Winter", "Spring"), c("Autumn", "Spring"), c("Winter", "Autumn"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

#Alpha-diversity analysis (season)
#Alpha-diversity analysis
(p = plot_richness(ps, x = "Season", measures=c("Observed", "Shannon")))
p7 <- p + geom_boxplot(data = p$data, aes(x = Season, y = value, fill = Season), alpha = 0.8, lwd =0.8, outlier.shape = NA) + geom_jitter(aes(fill = Season), alpha = 0.6, stroke = 1.0, shape = 21, size = 3.0, position = position_jitterdodge()) + theme_bw(base_size = 12, base_line_size = 1.5) + theme (axis.text.x = element_text(size=12, face = "bold", colour = "black", angle = 45, hjust=1)) + theme (axis.text.y = element_text(size = 12)) + theme(text = element_text(size = 12, face = "bold")) + scale_fill_manual(values = c("steelblue3", "orangered4", "gold3", "aquamarine4")) + theme(legend.position = "none") + theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "transparent", color = NA), panel.grid.major = element_blank(), panel.grid.minor = element_line(colour = "white")) + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, label = "p.signif", symnum.args = symnum.args)
p7$layers <- p7$layers[-1]
p7

##Alpha-diversity (in water)
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

#Stat-compare-means (add p-value)
my_comparisons <- list( c("Summer", "Winter"), c("Summer", "Z_Spring"), c("Summer", "T_Autumn"), c("Winter", "Z_Spring"), c("T_Autumn", "Z_Spring"), c("Winter", "Autumn"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

#Alpha-diversity analysis (season)
#Alpha-diversity analysis
(p = plot_richness(ps, x = "Season", measures=c("Observed", "Shannon")))
p8 <- p + geom_boxplot(data = p$data, aes(x = Season, y = value, fill = Season), alpha = 0.8, lwd =0.8, outlier.shape = NA) + geom_jitter(aes(fill = Season), alpha = 0.6, stroke = 1.0, shape = 21, size = 3.0, position = position_jitterdodge()) + theme_bw(base_size = 12, base_line_size = 1.5) + theme (axis.text.x = element_text(size=12, face = "bold", colour = "black", angle = 45, hjust=1)) + theme (axis.text.y = element_text(size = 12)) + theme(text = element_text(size = 12, face = "bold")) + scale_fill_manual(values = c("steelblue3", "orangered4", "gold3", "aquamarine4")) + theme(legend.position = "none") + theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "transparent", color = NA), panel.grid.major = element_blank(), panel.grid.minor = element_line(colour = "white")) + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, label = "p.signif", symnum.args = symnum.args)
p8$layers <- p8$layers[-1]
p8

##Arrange plots
top_row <- plot_grid(p5, p6, p7, labels = "AUTO", rel_widths = c(2, 2, 2), rel_heights = c(2, 2, 2), ncol = 1, align = "v")
top_row

##Save plot
ggsave("alpha-div_location_season_source.tiff", units = c("in"), width=8, height=12, dpi=300, compression="lzw")

##Alpha-diversity (Cyanobacteria)
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

#Stat-compare-means (add p-value)
my_comparison <- list(c("Summer", "T_Autumn"), c("Summer", "Winter"), c("Summer", "Z_Spring"), c("T_Autumn", "Winter"), c("T_Autumn", "Z_Spring"), c("Winter", "Z_Spring"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

#Alpha-diversity analysis (season)
#Alpha-diversity analysis
(p = plot_richness(ps, x = "Season", measures=c("Observed", "Shannon")))
p10 <- p + geom_boxplot(data = p$data, aes(x = Season, y = value, fill = Season), alpha = 0.8, lwd =0.8, outlier.shape = NA) + geom_jitter(aes(fill = Season), alpha = 0.6, stroke = 1.0, shape = 21, size = 3.0, position = position_jitterdodge()) + theme_bw(base_size = 12, base_line_size = 1.5) + theme (axis.text.x = element_text(size=12, face = "bold", colour = "black", angle = 45, hjust=1)) + theme (axis.text.y = element_text(size = 12)) + theme(text = element_text(size = 12, face = "bold")) + scale_fill_manual(values = c("steelblue3", "orangered4", "gold3", "aquamarine4")) + theme(legend.position = "none") + theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "transparent", color = NA), panel.grid.major = element_blank(), panel.grid.minor = element_line(colour = "white")) + stat_compare_means(method = "wilcox.test", comparisons = my_comparison, label = "p.signif", symnum.args = symnum.args) + ggtitle ("Cyanobacteria in water")
p10$layers <- p10$layers[-1]
p10

##Save plot
ggsave("alpha-div_cyano_water_season_source.tiff", units = c("in"), width=7, height=6, dpi=300, compression="lzw")

######Beta-diversity######
#Un-weighted
library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
plot(random_tree)
physeq1 = merge_phyloseq(physeq, META, random_tree)
ps0.rar.filtered <- core(physeq, detection = 10, prevalence = 0.05)
summarize_phyloseq(ps0.rar.filtered)

#Source
ordu.unwt.uni <- ordinate(physeq1, "PCoA", "unifrac", weighted=F)
unwt.unifrac <- plot_ordination(physeq1, ordu.unwt.uni, color="Source", shape="Source") 
unwt.unifrac <- unwt.unifrac + ggtitle("Unweighted UniFrac")
unwt.unifrac <- unwt.unifrac + theme_bw()
unw_ord1 <- unwt.unifrac + geom_point(size = 4, alpha = 0.9) + theme_bw(base_size = 20, base_line_size = 1.5) + theme(text = element_text(size = 12, face = "bold"))+ stat_ellipse(aes(group = Source), geom="polygon",level=0.8,alpha=0.2) + labs("Source") + scale_color_manual(values = c("steelblue3", "orangered4", "gold3", "aquamarine4"), name = "Source") + scale_shape_manual(values = c(19, 15), name = "Source") + theme(legend.position = "bottom") + theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "transparent", color = NA), panel.grid.major = element_blank(), panel.grid.minor = element_line(colour = "white"))
unw_ord1

#Location
ordu.unwt.uni <- ordinate(physeq1, "PCoA", "unifrac", weighted=F)
unwt.unifrac <- plot_ordination(physeq1, ordu.unwt.uni, color="Location", shape="Location") 
unwt.unifrac <- unwt.unifrac + ggtitle("Unweighted UniFrac")
unwt.unifrac <- unwt.unifrac + theme_bw()
unw_ord2 <- unwt.unifrac + geom_point(size = 4, alpha = 0.9) + theme_bw(base_size = 20, base_line_size = 1.5) + theme(text = element_text(size = 12, face = "bold"))+ stat_ellipse(aes(group = Location), geom="polygon",level=0.8,alpha=0.2) + labs("Location") + scale_color_manual(values = c("steelblue3", "orangered4", "gold3", "aquamarine4"), name = "Location") + scale_shape_manual(values = c(19, 15, 18, 17), name = "Location") + theme(legend.position = "bottom") + theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "transparent", color = NA), panel.grid.major = element_blank(), panel.grid.minor = element_line(colour = "white"))
unw_ord2

#Season
ordu.unwt.uni <- ordinate(physeq1, "PCoA", "unifrac", weighted=F)
unwt.unifrac <- plot_ordination(physeq1, ordu.unwt.uni, color="Season", shape="Season")
unwt.unifrac <- unwt.unifrac + ggtitle("Unweighted UniFrac")
unwt.unifrac <- unwt.unifrac + theme_bw()
unw_ord3 <- unwt.unifrac + geom_point(size = 4, alpha = 0.9) + theme_bw(base_size = 20, base_line_size = 1.5) + theme(text = element_text(size = 12, face = "bold"))+ stat_ellipse(aes(group = Season), geom="polygon",level=0.8,alpha=0.2) + labs("Season") + scale_color_manual(values = c("steelblue3", "orangered4", "gold3", "aquamarine4"), name = "Season") + scale_shape_manual(values = c(19, 15, 18, 17), name = "Season") + theme(legend.position = "bottom") + theme(panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "transparent", color = NA), panel.grid.major = element_blank(), panel.grid.minor = element_line(colour = "white"))
unw_ord3

#PERMANOVA_ANODIS
library(vegan)
metadf <- data.frame(sample_data(physeq1))
unifrac.dist <- UniFrac(physeq1, weighted = FALSE, normalized = TRUE, parallel = FALSE, fast = TRUE)
permanova <- adonis2(unifrac.dist ~ Source, data = metadf)
permanova

#Beta-dispersion
ps.disper <- betadisper(unifrac.dist, metadf$Location)
permutest(ps.disper, pairwise = TRUE)

#Beta-dispersion
ps.disper <- betadisper(unifrac.dist, metadf$Season)
permutest(ps.disper, pairwise = TRUE)

#PERMANOVA_ANODIS
library(vegan)
metadf <- data.frame(sample_data(physeq1))
unifrac.dist <- UniFrac(physeq1, weighted = TRUE, normalized = TRUE, parallel = FALSE, fast = TRUE)
permanova <- adonis2(unifrac.dist ~ Source, data = metadf)
permanova

#Beta-dispersion
ps.disper <- betadisper(unifrac.dist, metadf$Location)
permutest(ps.disper, pairwise = TRUE)

#Beta-dispersion
ps.disper <- betadisper(unifrac.dist, metadf$Season)
permutest(ps.disper, pairwise = TRUE)

##Arrange plots
m_row <- plot_grid(unw_ord1, unw_ord2, unw_ord3, labels = c('A', 'B', 'C'), rel_widths = c(2, 2, 2), rel_heights = c(2, 2, 2), ncol = 1, align = "v")
m_row
b_row <- plot_grid(w_ord1, w_ord2, w_ord3, labels = c('D', 'E', 'F'), rel_widths = c(2, 2, 2), rel_heights = c(2, 2, 2), ncol = 1, align = "v")
b_row
plot_grid(m_row, b_row, ncol = 2, rel_widths = c(2, 2))

##Save plot
ggsave("beta-div_location_season_source.tiff", units = c("in"), width=10, height=14, dpi=300, compression="lzw")