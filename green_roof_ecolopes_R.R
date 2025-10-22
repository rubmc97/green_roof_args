
#Analyses of the ARGs in Green Roof experiment | ECOLOPES

#Load all the packages
library(vegan)
library(dplyr)
library(tidyr)
library(wesanderson)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(pheatmap)
library(RColorBrewer)
library(grid)
library(tibble)
library(linkET)
library(tidyverse)
library(stringr)

#Import Argo's output and normalize gene counts by total bases:
barcode48 = read.table("SQK-NBD114-96_barcode48_trimmed_filt.sarg.tsv", sep = "\t", header = TRUE, quote = "", fill = TRUE, comment.char = "")
barcode48$total_bases = 597942795
barcode48$CPM = (barcode48$copy / barcode48$total_bases) * 1000000
barcode48$sampleID = "bar48"

barcode49 = read.table("SQK-NBD114-96_barcode49_trimmed_filt.sarg.tsv", sep = "\t", header = TRUE, quote = "", fill = TRUE, comment.char = "")
barcode49$total_bases = 561477907
barcode49$CPM = (barcode49$copy / barcode49$total_bases) * 1000000
barcode49$sampleID = "bar49"

barcode50 = read.table("SQK-NBD114-96_barcode50_trimmed_filt.sarg.tsv", sep = "\t", header = TRUE, quote = "", fill = TRUE, comment.char = "")
barcode50$total_bases = 940634141
barcode50$CPM = (barcode50$copy / barcode50$total_bases) * 1000000
barcode50$sampleID = "bar50"

barcode51 = read.table("SQK-NBD114-96_barcode51_trimmed_filt.sarg.tsv", sep = "\t", header = TRUE, quote = "", fill = TRUE, comment.char = "")
barcode51$total_bases = 1061326664
barcode51$CPM = (barcode51$copy / barcode51$total_bases) * 1000000
barcode51$sampleID = "bar51"

barcode52 <- read.table("SQK-NBD114-96_barcode52_trimmed_filt.sarg.tsv", sep = "\t", header = TRUE, quote = "", fill = TRUE, comment.char = "")
barcode52$total_bases = 1241342219
barcode52$CPM = (barcode52$copy / barcode52$total_bases) * 1000000
barcode52$sampleID = "bar52"

barcode53 <- read.table("SQK-NBD114-96_barcode53_trimmed_filt.sarg.tsv", sep = "\t", header = TRUE, quote = "", fill = TRUE, comment.char = "")
barcode53$total_bases = 1731845091
barcode53$CPM = (barcode53$copy / barcode53$total_bases) * 1000000
barcode53$sampleID = "bar53"

barcode54 <- read.table("SQK-NBD114-96_barcode54_trimmed_filt.sarg.tsv", sep = "\t", header = TRUE, quote = "", fill = TRUE, comment.char = "")
barcode54$total_bases = 2399567608
barcode54$CPM = (barcode54$copy / barcode54$total_bases) * 1000000
barcode54$sampleID = "bar54"

barcode55 <- read.table("SQK-NBD114-96_barcode55_trimmed_filt.sarg.tsv", sep = "\t", header = TRUE, quote = "", fill = TRUE, comment.char = "")
barcode55$total_bases = 2510562365
barcode55$CPM = (barcode55$copy / barcode55$total_bases) * 1000000
barcode55$sampleID = "bar55"

barcode56 <- read.table("SQK-NBD114-96_barcode56_trimmed_filt.sarg.tsv", sep = "\t", header = TRUE, quote = "", fill = TRUE, comment.char = "")
barcode56$total_bases = 2237716099
barcode56$CPM = (barcode56$copy / barcode56$total_bases) * 1000000
barcode56$sampleID = "bar56"

barcode57 <- read.table("SQK-NBD114-96_barcode57_trimmed_filt.sarg.tsv", sep = "\t", header = TRUE, quote = "", fill = TRUE, comment.char = "")
barcode57$total_bases = 1224418604
barcode57$CPM = (barcode57$copy / barcode57$total_bases) * 1000000
barcode57$sampleID = "bar57"

barcode58 <- read.table("SQK-NBD114-96_barcode58_trimmed_filt.sarg.tsv", sep = "\t", header = TRUE, quote = "", fill = TRUE, comment.char = "")
barcode58$total_bases = 2952951589
barcode58$CPM = (barcode58$copy / barcode58$total_bases) * 1000000
barcode58$sampleID = "bar58"

barcode59 <- read.table("SQK-NBD114-96_barcode59_trimmed_filt.sarg.tsv", sep = "\t", header = TRUE, quote = "", fill = TRUE, comment.char = "")
barcode59$total_bases = 2064630817
barcode59$CPM = (barcode59$copy / barcode59$total_bases) * 1000000
barcode59$sampleID = "bar59"

args.detected.gr.cpm = rbind(barcode48, barcode49, barcode50, barcode51, barcode52, barcode53, barcode54, barcode55, barcode56, barcode57, barcode58, barcode59)
write.csv(args.detected.gr.cpm, "args_detected_argo_cpm_green_roofs_0925.csv")

sd.args = read.csv("sd_rooftops.csv") #load metadata
df.args.detected.gr.cpm = merge.data.frame(args.detected.gr.cpm, sd.args, by.y = "sampleID")

df.args.sum = df.args.detected.gr.cpm  %>%
  group_by(sampleID, treatment) %>%
  summarise(total_CPM = sum(CPM, na.rm = TRUE)) #Calculate the total abundance of ARGs per sample

palette.treat = wes_palette("Moonrise2") #Palette for the boxplots

custom.order = c("G1", "G2", "G3", "G4")

# Convert the total_CPM column to numeric
df.args.sum$total_CPM = as.numeric(df.args.sum$total_CPM)

#Plotting ARGs total abundance across treatments
args.total.bp = ggplot(df.args.sum, aes(x = treatment, y = total_CPM, fill = treatment)) +
  geom_boxplot(alpha = 0.75, width = 0.5) +
  geom_jitter(size = 1, shape = 16, alpha = 0.45) +
  ylab(label = "CPM") +
  xlab(label = "Treatment") +
  scale_fill_manual(values = palette.treat) +  # Set custom colors
  theme_pubr() +
  labs(title = "ARGs total") + stat_compare_means(method = "kruskal.test", label.y = max(df.args.sum$total_CPM) * 1.05) +
  theme(plot.title = element_text(size = 11),
        axis.title = element_text(size = 11),  # Axis label sizes
        axis.text = element_text(size = 11)) 
#stat_compare_means(label = "p.adjust", hide.ns = TRUE, method = "t.test",
#comparisons = list(c("G1", "G2"), c("G1", "G3"), c("G1", "G4"),
#c("G2", "G3"), c("G2", "G4"), c("G3", "G4")),
#p.adjust.method = "BH") #Pairwise comparisons if needed
args.total.bp

#Calculate ARGs richness
args.rich = df.args.detected.gr.cpm %>%
  filter(copy > 0) %>%  # Optional: only keep detected ARGs
  group_by(sampleID) %>%
  summarise(richness = n_distinct(subtype))

df.args.rich = merge.data.frame(args.rich , sd.args, by.y = "sampleID")

#Plotting ARGs richness across treatments
args.richness.bp <- ggplot(df.args.rich, aes(x = treatment, y = richness, fill = treatment)) +
  geom_boxplot(alpha = 0.6, width = 0.5) +
  geom_jitter(size = 1, shape = 16, alpha = 0.45) +
  ylab("Richness (unique ARGs detected)") +
  xlab("Treatment") +
  scale_fill_manual(values = palette.treat) +  # Set custom colors
  theme_pubr() +
  stat_compare_means(method = "kruskal.test", label.y = max(df.args.rich$richness) * 1.05) +
  labs(title = "ARGs richness") +
  theme(
    plot.title = element_text(size = 11),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11)
  )
args.richness.bp

fig1AB = ggarrange(args.total.bp, args.richness.bp, ncol = 2, common.legend = T)

#aggregate until the obtention of a data matrix for the PCoA
agg.cpm = args.detected.gr.cpm %>%
  group_by(sampleID, subtype, type) %>%
  summarise(total_CPM = sum(CPM, na.rm = TRUE)) %>%
  ungroup()

mat.agg.cpm = agg.cpm %>%
  pivot_wider(
    names_from = subtype,
    values_from = total_CPM,
    values_fill = 0   # fill missing subtype/sample combinations with 0
  )

mat.agg.cpm = mat.agg.cpm %>%
  tibble::column_to_rownames(var = "sampleID")

#Hellinger transformation
arg.hell = decostand(mat.agg.cpm, method = "hellinger")

#Distance with the Bray-Curtis index
dist.mat = vegdist(arg.hell, method = "bray")
#pcoa.res = cmdscale(dist.mat, k = 2, eig = TRUE)
pcoa.res = ape::pcoa(dist.mat)

row.names(sd.gr) = sd.gr$SampleID
pcoa.df = data.frame(SampleID = rownames(pcoa.res$vectors),
                     Axis1 = pcoa.res$vectors[, 1],
                     Axis2 = pcoa.res$vectors[, 2])
pcoa.arg.df = merge(pcoa.df, sd.gr, by = "SampleID")

sd.gr <- sd.gr[labels(dist.mat), ]
identical(rownames(sd.gr), labels(dist.mat))

#Permanova tests to address the significance of the variables's effect on ARG composition
adonis2(dist.mat ~treatment, data = sd.gr, permutations = 999)
adonis2(dist.mat ~greenwaste, data = sd.gr, permutations = 999) 
adonis2(dist.mat ~plants, data = sd.gr, permutations = 999)
adonis2(dist.mat ~year, data = sd.gr, permutations = 999)
adonis2(dist.mat ~greenwaste:year, data = sd.gr, permutations = 999)

# Create the PERMANOVA table as a tableGrob
adonis.out = paste(
  "Factor            R2      P-value",
  "treatment         0.543   0.006",
  "greenwaste        0.236   0.008",
  "plants            0.481   0.001",
  "year              0.186   0.048",
  "greenwaste:year   0.306   0.024",
  sep = "\n"
)

palette.treat = wes_palette("Moonrise2")

pdf("FigC_pcoa_bray_ARGs.pdf", width = 3, height = 4)

#Plot the PCoA based on ARG composition
pcoa.bray.plot = ggplot(pcoa.arg.df, aes(x = Axis1, y = Axis2, color = treatment)) +  # Adjust "LandUse" as needed
  geom_point(aes(color = factor(treatment), shape = factor(greenwaste)), size = 5, alpha = 0.75) + scale_color_manual(values = palette.treat) +
  labs(x = paste0("PCoA1 (", round(pcoa.res$values$Relative_eig[1] * 100, 2), "%)"),
       y = paste0("PCoA2 (", round(pcoa.res$values$Relative_eig[2] * 100, 2), "%)"),
       title = "PCoA AMR gene composition") + theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  annotate(
    "text",
    x = max(pcoa.arg.df$Axis1) * 0.95,
    y = max(pcoa.arg.df$Axis2) * 0.95,
    label = adonis.out,
    hjust = 1, vjust = 1,
    family = "Arial",
    size = 3.5    # ~10 pt
  )

pcoa.bray.plot

dev.off()

#Aggregate until obtaining a data frame for the heatmap
df.args.detected.gr.cpm = df.args.detected.gr.cpm %>%
  select(type, subtype, CPM, SampleID, treatment, sampleID)

#Square root of the CPM to improve visualization
df.args.detected.gr.cpm$sqrtCPM = sqrt(df.args.detected.gr.cpm$CPM)

treat = read.csv("args_treatment.csv", row.names = 1)
#annotation_row = read.csv("args_merged_heat.csv", row.names = 1)

args.mat = reshape2::dcast(df.args.detected.gr.cpm, SampleID ~ subtype, value.var = "sqrtCPM")

args.mat = tibble::column_to_rownames(.data = args.mat, var = "SampleID")

treat$Treatment = factor(treat$Treatment, levels = c("G1", "G2", "G3", "G4"))

annotation_row <- df.args.detected.gr.cpm %>%
  select(type, subtype) %>%
  distinct(subtype, .keep_all = TRUE) %>%
  tibble::column_to_rownames("subtype")


annotation_row$type = factor(annotation_row$type) 


# Specify the order of 'class' levels
class.order = c("aminocoumarin", "aminoglycoside", "bacitracin", "beta-lactam", "biocide", "bleomycin", "capuramycin", "colistin",  "defensin", "fosfomycin", "glycopeptide", "macrolide-lincosamide-streptogramin", "multidrug", "phenicol", "rifamycin", "tetracycline", "trimethoprim")

# Convert 'class' to a factor with the specified order
annotation_row$type = factor(annotation_row$type, levels = class.order)

palette.treat = wes_palette("Moonrise2", n = length(unique(treat$Treatment)))
col_annotation_colors = setNames(palette.treat, levels(treat$Treatment))

# Create a color palette for the 'class' variable in annotation_row
palette_class = c(
  "firebrick",     # Firebrick red
  "steelblue",     # Steel blue
  "springgreen4",   # Forest green
  "chocolate3",     # Chocolate brown
  "gold2",
  "royalblue3",# Gold
  "sienna",   # Saddle brown
  "olivedrab",     # Olive drab
  "darkgrey",
  "lightpink3",# Dark grey
  "slateblue",     # Slate blue
  "wheat",
  "tomato2",
  "maroon",
  "khaki3",
  "orchid4",
  "tan2" # Tan
)
# Map the palette to the 'class' variable
row_annotation_colors = setNames(palette_class, levels(annotation_row$type))
names(palette_class) = levels(annotation_row$type)

bc_args = vegdist(args.mat, method = "bray")

pdf("FigD_heatmap.pdf", width = 3, height = 4)

#Plotting the heatmap
pheatmap(
  args.mat,
  annotation_col = treat,              
  annotation_row = annotation_row,              
  annotation_colors = list(
    Treatment = col_annotation_colors,      # must match treat column name
    type = row_annotation_colors            # must match annotation_row column name
  ),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  clustering_distance_rows = bc_args,
  scale = "none"
)

dev.off()

#Import physicochemical parameters
abiot = read.csv("abiotic_nuts_green_roofs.csv", header = TRUE, row.names = 1)

top.args = agg.cpm %>%
  group_by(subtype) %>%
  summarise(total_CPM_sum = sum(total_CPM, na.rm = TRUE)) %>%
  arrange(desc(total_CPM_sum))

top.group = agg.cpm %>%
  group_by(type) %>%
  summarise(total_CPM_sum = sum(total_CPM, na.rm = TRUE)) %>%
  arrange(desc(total_CPM_sum))

#Filter the data frame based on the top 10 most abundant ARGs
sub.agg.cpm = agg.cpm %>%
  filter(subtype %in% c("bacA", "sugE", "arr", "bla*", "rox", "aph(3')-II", "tlmA", "rph*", "helR", "lysX"))

sub.agg.cpm_wide = sub.agg.cpm %>%
  group_by(sampleID, subtype) %>%              # sum CPM if duplicates exist
  summarise(total_CPM = sum(total_CPM, na.rm=TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = subtype,
    values_from = total_CPM,
    values_fill = 0
  ) %>%
  column_to_rownames("sampleID")

pdf("FigE_pearson_argo_args_abiot_top10.pdf", width = 3, height = 4)
correlate(sub.agg.cpm_wide, abiot) %>% 
  qcorrplot() +
  geom_square() +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"))

dev.off()

#Plot biomass ~ Cmic
cmic.bp = ggplot(abiot, aes(x = Treatment, y = Cmic, fill = Treatment)) +
  geom_boxplot(alpha = 0.6, width = 0.5) +
  geom_jitter(size = 1, shape = 16, alpha = 0.45) +
  ylab("Cmic (ug/g)") +
  xlab("Treatment") +
  scale_fill_manual(values = palette.treat) +  # Set custom colors
  theme_pubr() +
  stat_compare_means(method = "kruskal.test", p.adjust.methods = "BH") +
  labs(title = "Bacterial biomass") +
  theme(
    plot.title = element_text(size = 11),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11)
  )

####
#Load emu's output of the taxonomic profiling based on the 16S to assess bacterial diversity
emu.g11 = read.table("SQK-NBD114-96_barcode48_16S_rel-abundance.tsv", sep = "\t", header = T)
emu.g11 = emu.g11 %>%
  filter(species != "") %>%
  filter(superkingdom == "Bacteria") #remove everything that isn't bacteria

emu.g12 = read.table("SQK-NBD114-96_barcode49_16S_rel-abundance.tsv", sep = "\t", header = T)
emu.g12 = emu.g12 %>%
  filter(species != "") %>%
  filter(superkingdom == "Bacteria")

emu.g13 = read.table("SQK-NBD114-96_barcode50_16S_rel-abundance.tsv", sep = "\t", header = T)
emu.g13 = emu.g13 %>%
  filter(species != "") %>%
  filter(estimated.counts != "") %>%
  filter(superkingdom == "Bacteria")

emu.g21 = read.table("SQK-NBD114-96_barcode51_16S_rel-abundance.tsv", sep = "\t", header = T)
emu.g21 = emu.g21 %>%
  filter(species != "") %>%
  filter(superkingdom == "Bacteria")

emu.g22 = read.table("SQK-NBD114-96_barcode52_16S_rel-abundance.tsv", sep = "\t", header = T)
emu.g22 = emu.g22 %>%
  filter(species != "") %>%
  filter(superkingdom == "Bacteria")

emu.g23 = read.table("SQK-NBD114-96_barcode53_16S_rel-abundance.tsv", sep = "\t", header = T)
emu.g23 = emu.g22 %>%
  filter(species != "") %>%
  filter(superkingdom == "Bacteria")

emu.g31 = read.table("SQK-NBD114-96_barcode54_16S_rel-abundance.tsv", sep = "\t", header = T)
emu.g31 = emu.g31 %>%
  filter(species != "") %>%
  filter(estimated.counts != "") %>%
  filter(superkingdom == "Bacteria")

emu.g32 = read.table("SQK-NBD114-96_barcode55_16S_rel-abundance.tsv", sep = "\t", header = T)
emu.g32 = emu.g32 %>%
  filter(species != "") %>%
  filter(superkingdom == "Bacteria")

emu.g33 = read.table("SQK-NBD114-96_barcode56_16S_rel-abundance.tsv", sep = "\t", header = T)
emu.g33 = emu.g33 %>%
  filter(species != "") %>%
  filter(superkingdom == "Bacteria")

emu.g41 = read.table("SQK-NBD114-96_barcode57_16S_rel-abundance.tsv", sep = "\t", header = T)
emu.g41 = emu.g41 %>%
  filter(species != "") %>%
  filter(superkingdom == "Bacteria")

emu.g42 = read.table("SQK-NBD114-96_barcode58_16S_rel-abundance.tsv", sep = "\t", header = T)
emu.g42 = emu.g42 %>%
  filter(species != "") %>%
  filter(superkingdom == "Bacteria")

emu.g43 = read.table("SQK-NBD114-96_barcode59_16S_rel-abundance.tsv", sep = "\t", header = T)
emu.g43 = emu.g43 %>%
  filter(species != "") %>%
  filter(superkingdom == "Bacteria")

emu_list <- list(
  G11 = emu.g11,
  G12 = emu.g12,
  G13 = emu.g13,
  G21 = emu.g21,
  G22 = emu.g22,
  G23 = emu.g23,
  G31 = emu.g31,
  G32 = emu.g32,
  G33 = emu.g33,
  G41 = emu.g41,
  G42 = emu.g42,
  G43 = emu.g43
)

# Calculate TSS-normalized richness per sample
richness.gr = map_df(names(emu_list), function(name) {
  df = emu_list[[name]]
  
  # Skip empty samples
  if(sum(df$estimated.counts) == 0) {
    richness_tss <- NA
  } else {
    counts = df$estimated.counts
    names(counts) = df$species
    
    # TSS normalization
    tss = counts / sum(counts)
    
    # Richness: number of species with non-zero normalized abundance
    richness_tss = sum(tss > 0)
  }
  
  tibble(
    sampleID = name,
    treatment = str_extract(name, "^G\\d"),  # Extract treatment (e.g., G1)
    richness = richness_tss
  )
})

#Plot bacterial richness
rich.bp = ggplot(richness.gr, aes(x = treatment, y = richness, fill = treatment)) +
  geom_boxplot(alpha = 0.6, width = 0.5) +
  geom_jitter(size = 1, shape = 16, alpha = 0.45) +
  ylab("Richness (observed species)") +
  xlab("Treatment") +
  scale_fill_manual(values = palette.treat) +  # Set custom colors
  theme_pubr() +
  #stat_compare_means(comparisons = list(c("G1", "G2"), c("G2", "G3"), c("G3", "G4"),
  #                                     c("G1", "G3"), c("G2", "G4"), c("G1", "G4")),
  #                 label = "p.signif", hide.ns = TRUE, method = "wilcox.test",
  #                p.adjust.methods = "BH") +
  stat_compare_means(method = "kruskal.test") +
  labs(title = "Bacterial richness") +
  theme(
    plot.title = element_text(size = 11),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11)
  )

pdf("FigS2_biomass_diversity.pdf", width = 3, height = 4)
ggarrange(cmic.bp, rich.bp, ncol = 2, labels = c("A", "B"), common.legend = T)

####
#Load the output of Centrifuger from the taxonomic classification of the reads containing ARGs
g11.args = read.csv("centrifuger_ARG_with_taxa_barcode48.csv")
g11.args$SampleID = "G1.1"

g12.args = read.csv("centrifuger_ARG_with_taxa_barcode49.csv")
g12.args$SampleID = "G1.2"

g13.args = read.csv("centrifuger_ARG_with_taxa_barcode50.csv")
g13.args$SampleID = "G1.3"

g21.args = read.csv("centrifuger_ARG_with_taxa_barcode51.csv")
g21.args$SampleID = "G2.1"

g22.args = read.csv("centrifuger_ARG_with_taxa_barcode52.csv")
g22.args$SampleID = "G2.2"

g23.args = read.csv("centrifuger_ARG_with_taxa_barcode53.csv")
g23.args$SampleID = "G2.3"

g31.args = read.csv("centrifuger_ARG_with_taxa_barcode54.csv")
g31.args$SampleID = "G3.1"

g32.args = read.csv("centrifuger_ARG_with_taxa_barcode55.csv")
g32.args$SampleID = "G3.2"

g33.args = read.csv("centrifuger_ARG_with_taxa_barcode56.csv")
g33.args$SampleID = "G3.3"

g41.args = read.csv("centrifuger_ARG_with_taxa_barcode57.csv")
g41.args$SampleID = "G4.1"

g42.args = read.csv("centrifuger_ARG_with_taxa_barcode58.csv")
g42.args$SampleID = "G4.2"

g43.args = read.csv("centrifuger_ARG_with_taxa_barcode59.csv")
g43.args$SampleID = "G4.3"

#Create a unified data frame merging the output of Centrifuger
args.all = rbind(g11.args, g12.args, g13.args, g21.args, g22.args, g23.args, g31.args, g32.args, g33.args, g41.args, g42.args, g43.args)

args.all = args.all %>%
  separate(ARG_Hit, into = c("ARG", "ARG2"), sep = ";", fill = "right", extra = "merge", remove = FALSE)

args.all = args.all %>%
  separate(ARG_Hit, into = c("Family_ARG", "Gene", "ARG_code"), sep = "\\|", fill = "right", extra = "merge", remove = FALSE)

args.all = args.all %>%
  select(SampleID, Family_ARG, Gene, clade_name)

#Split the taxonomy
args.all.agg = args.all %>%
  separate(clade_name, 
           into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = "\\|", 
           fill = "right", 
           extra = "drop")

args.agg.genus = args.all.agg %>%
  mutate(Genus = ifelse(is.na(Genus), "Unclassified_Genus", Genus)) %>%
  group_by(SampleID, Family_ARG, Gene, Genus) %>%
  summarise(Gene_count = n(), .groups = "drop")

#Add metadata
df.args.agg.genus = merge.data.frame(args.agg.genus, sd.args, by.y = "SampleID")

#Remove g_ from the names
df.args.agg.genus <- df.args.agg.genus %>%
  mutate(Genus = str_remove(Genus, "^g__"))

#Aggregate
agg.by.genus <- df.args.agg.genus %>%
  group_by(Genus) %>%
  summarise(Total_gene_count = sum(Gene_count), .groups = "drop")

#Top 25 genera
top25 = agg.by.genus %>%
  arrange(desc(Total_gene_count)) %>%  # sort descending by count
  slice_head(n = 25) %>%         # take top 25
  pull(Genus)                     # extract as a vector

sum.df.args.agg = df.args.agg.genus %>%
  group_by(treatment, Family_ARG, Gene) %>%
  summarise(Total_gene_count = sum(Gene_count), .groups = "drop")

g1.args.genus = subset.data.frame(df.args.agg.genus, treatment == "G1")
g1.args.genus = g1.args.genus %>%
  select(Gene_count, Genus, Gene, SampleID)

top25.g1 = g1.args.genus %>% 
  filter(Genus %in% top25)

# Sum all others into "Other"
others_sum = g1.args.genus %>% 
  filter(!Genus %in% top25) %>% 
  group_by(SampleID, Gene) %>%
  summarise(
    Genus = "Other",
    Gene_count = sum(Gene_count),
    .groups = "drop"
  )

g1.args.genus.top25 = bind_rows(top25.g1, others_sum)

#Colour vector for all the genera names across treatments
colors26.dark = c(
  "Unclassified_Genus" = "royalblue1",
  "Streptomyces" = "orange2",
  "Bradyrhizobium" = "tomato1",
  "Mycobacterium" = "steelblue",
  "Allosphingosinicella" = "#59A14F",
  "Lentzea" = "darkgoldenrod1",
  "Blastococcus" = "pink3",
  "Sphingomicrobium" = "plum1",
  "Longimicrobium" = "#2E4057",
  "Nocardioides" = "forestgreen",
  "Penaeicola" = "orchid4",
  "Actinoplanes" = "firebrick",
  "Ramlibacter" = "khaki2",
  "Peristeroidobacter" = "sienna",
  "Phytohabitans" = "gold1",
  "Pseudonocardia" = "cyan1",
  "Microvirga" = "#393B79",
  "Micromonospora" = "#637939",
  "Nonomuraea" = "#8C6D31",
  "Arthrobacter" = "#843C39",
  "Phenylobacterium" = "#3182BD",
  "Actinophytocola" = "#E6550D",
  "Bosea" = "seagreen3",
  "Usitatibacter" = "#756BB1",
  "Virgisporangium" = "maroon",
  "Other" = "#636363"
)

bp.gen.args.g1 = ggplot(g1.args.genus.top25, aes(x = Gene, y = Gene_count, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = colors26.dark,
    guide = guide_legend(
      label.theme = element_text(face = "italic")  # legend text in italics
    )
  ) + ggtitle("G1") +
  labs(x = "Gene", y = "Gene Copies", fill = "Genus") +
  ggpubr::theme_pubr() +
  scale_y_continuous(limits = c(0, 75), breaks = c(0, 25, 50, 75)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")  # x-axis labels in italics
  )

#Same for G2
g2.args.genus = subset.data.frame(df.args.agg.genus, treatment == "G2")
g2.args.genus = g2.args.genus %>%
  select(Gene_count, Genus, Gene, SampleID)

top25.g2 = g2.args.genus %>% 
  filter(Genus %in% top25)

others_sum <- g2.args.genus %>% 
  filter(!Genus %in% top25) %>% 
  group_by(SampleID, Gene) %>%
  summarise(
    Genus = "Other",
    Gene_count = sum(Gene_count),
    .groups = "drop"
  )

# Combine top25 and summed others
g2.args.genus.top25 = bind_rows(top25.g2, others_sum)

bp.gen.args.g2 = ggplot(g2.args.genus.top25, aes(x = Gene, y = Gene_count, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = colors26.dark,
    guide = guide_legend(
      label.theme = element_text(face = "italic")  # legend text in italics
    )
  ) + ggtitle("G2") +
  labs(x = "Gene", y = "Gene Copies", fill = "Genus") +
  ggpubr::theme_pubr() +
  scale_y_continuous(limits = c(0, 75), breaks = c(0, 25, 50, 75)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")  # x-axis labels in italics
  )

#Same for G3
g3.args.genus = subset.data.frame(df.args.agg.genus, treatment == "G3")
g3.args.genus = g3.args.genus %>%
  select(Gene_count, Genus, Gene, SampleID)

top25.g3 = g3.args.genus %>% 
  filter(Genus %in% top25)

others_sum <- g3.args.genus %>% 
  filter(!Genus %in% top25) %>% 
  group_by(SampleID, Gene) %>%
  summarise(
    Genus = "Other",
    Gene_count = sum(Gene_count),
    .groups = "drop"
  )

# Combine top25 and summed others
g3.args.genus.top25 = bind_rows(top25.g3, others_sum)

bp.gen.args.g3 <- ggplot(g3.args.genus.top25, aes(x = Gene, y = Gene_count, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = colors26.dark,
    guide = guide_legend(
      label.theme = element_text(face = "italic")  # legend text in italics
    )
  ) + ggtitle("G3") +
  labs(x = "Gene", y = "Gene Copies", fill = "Genus") +
  ggpubr::theme_pubr() +
  scale_y_continuous(limits = c(0, 175), breaks = c(0, 25, 50, 75, 100, 125, 150, 175)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")  # x-axis labels in italics
  )

#Same for G4
g4.args.genus = subset.data.frame(df.args.agg.genus, treatment == "G4")
g4.args.genus = g4.args.genus %>%
  select(Gene_count, Genus, Gene, SampleID)

top25.g4 = g4.args.genus %>% 
  filter(Genus %in% top25)

others_sum <- g4.args.genus %>% 
  filter(!Genus %in% top25) %>% 
  group_by(SampleID, Gene) %>%
  summarise(
    Genus = "Other",
    Gene_count = sum(Gene_count),
    .groups = "drop"
  )

# Combine top25 and summed others
g4.args.genus.top25 = bind_rows(top25.g4, others_sum)

bp.gen.args.g4 = ggplot(g4.args.genus.top25, aes(x = Gene, y = Gene_count, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = colors26.dark,
    guide = guide_legend(
      label.theme = element_text(face = "italic")  # legend text in italics
    )
  ) + ggtitle("G4") +
  labs(x = "Gene", y = "Gene Copies", fill = "Genus") +
  ggpubr::theme_pubr() +
  scale_y_continuous(limits = c(0, 175), breaks = c(0, 25, 50, 75, 100, 125, 150, 175)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic")  # x-axis labels in italics
  )

pdf("Fig2_ARGs_with_taxa.pdf", width = 12, height = 9)

#Put the plots together
ggarrange(bp.gen.args.g1, bp.gen.args.g2, bp.gen.args.g3, bp.gen.args.g4, ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"), common.legend = T)
