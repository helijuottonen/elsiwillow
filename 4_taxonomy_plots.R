# Heli Juottonen 2019-2022 heli.juottonen@alumni.helsinki.fi
# R scripts for Elsi Hietaranta's willow manuscript
# https://github.com/helijuottonen/elsiwillow

# 4: taxonomy plots

# Disclaimer: I'm a self-taught R user. I try my very best to make sure that the scripts work, but I do not 
# guarantee that they are correct, make sense or follow any kind of best practices. Use at your own risk!

# RUN FIRST: 0_basic_processing.R

#library(vegan)
library(ggplot2) #v.3.3.5
library(cowplot) #v.1.1.1
library(phyloseq) #1.32.0
library(dplyr) #v.1.0.8
library(tidyr) #v.1.2.0

##### BACTERIA

# creating a phyloseq object

# defining taxonomy file
taxfile.bact = "data/pajut_bothruns.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.cons.taxonomy.txt"

#importing the taxonomy table into phyloseq format
mothur_tax.bact <- import_mothur(mothur_constaxonomy_file = taxfile.bact)

# defining OTU table object
otutable.bact.tax <- otu_table(as.matrix(elotub.r), taxa_are_rows=FALSE)
# defining metadata object
metadata.bact <- sample_data(elmetab)

elotub.tax <- merge_phyloseq(otutable.bact.tax, metadata.bact, mothur_tax.bact)

# giving proper names to the columns of the taxonomy tables
colnames(tax_table(elotub.tax)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

# check how the phyloseq object looks (should have 785 taxa and 95 samples)
elotub.tax

# merging OTUs by order

elotub_order <- elotub.tax %>%
  tax_glom(taxrank = "Order") %>%                     # agglomerate at class level
  #transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Order)  

# getting rid of extra metadata columns because they seem to mess up pivot_wider
elotub_order2 <- subset(elotub_order, select = c(2,3,5,6,8,14))

length(unique(elotub_order2[,"Order"]))
# 70 orders - too many for plotting ->

# converting to wide for selecting orders with rel.abund. >1%
elotub_order.wide <- pivot_wider(elotub_order2, names_from = Order, values_from = Abundance)
bcolsums <- colSums(elotub_order.wide[5:74])
brel <- as.data.frame(bcolsums/sum(bcolsums))
names(brel)[1] <- "rel.abundance"
# selecting orders with overall relative abundance over 1%
brel_sub <- subset(brel, rel.abundance > 0.01)
# orders > 1% as a vector
ord_over1 <- row.names(brel_sub)

elotub_order_other <- elotub_order2 %>%
  # renaming orders <1% as 'other'
  mutate(Order_mod = ifelse(Order %in% ord_over1, Order, "other")) %>%
  # converting read counts to relative abundances
  mutate(Rel.Abundance = Abundance/sum(Abundance)) %>%
  # removing orders/rows with no reads
  filter(Abundance !=0)

length(unique(elotub_order_other[,"Order_mod"]))
#count(distinct(elotub_order_other, Order_mod))
# should be 13 = 13 orders left, one of them 'other'

# defining colours for plotting
bpal <- c("#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD", "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", "#8569D5", "#5E738F","#D1A33D", "olivedrab4", "yellowgreen", "turquoise4", "navy", "lightgreen", "honeydew4")

# setting treatment names for plotting
elotub_order_other$treat <- as.factor(elotub_order_other$treat)
levels(elotub_order_other$treat)
levels(elotub_order_other$treat) <- c("open", "bagged")

# setting site names for plotting
elotub_order_other$site <- as.factor(elotub_order_other$site)
levels(elotub_order_other$site)
levels(elotub_order_other$site) <- c("A", "B", "C", "D", "E", "F")

# plotting
elotub_order_plot <- ggplot(elotub_order_other, aes(x = site, y = Rel.Abundance, fill = Order_mod, color=Order_mod)) + 
  geom_bar(stat = "identity", position="fill") +
  #scale_fill_manual(values = bpal, limits=force, breaks = c(ord_over1, "other"), name="Order") +
  scale_fill_manual(values = bpal, limits=force, name="Order") +
  scale_colour_manual(values = bpal, limits=force, guide="none") +
  ylab("Relative abundance") +
  theme_cowplot(12) + 
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  facet_grid(. ~ honeybees + treat, scales="free_x", space="free")
  
elotub_order_plot

#tax.plots <- plot_grid(elotub_order_plot, xxxx, labels=c("a)", "b)"), label_fontface="plain", ncol=1, align="v")

elotub_order_plot
pdf("plots/FigS2a_elotub_taxplot.pdf") 
print(elotub_order_plot)
dev.off()

########### FUNGI

# creating a phyloseq object

# defining taxonomy file
taxfile.fungi = "data/pajubothx.fasta.ITS2.precluster.pick.pick.agc.0.03.cons.taxonomy.txt"

#importing the taxonomy table into phyloseq format
mothur_tax.fungi <- import_mothur(mothur_constaxonomy_file = taxfile.fungi)

# defining OTU table object
otutable.fungi.tax <- otu_table(as.matrix(elotuf.r), taxa_are_rows=FALSE)
# defining metadata object
metadata.fungi <- sample_data(elmetab)

elotuf.tax <- merge_phyloseq(otutable.fungi.tax, metadata.fungi, mothur_tax.fungi)

# giving proper names to the columns of the taxonomy tables
# NOTE 4.11.20: for some reason mothur/Unite databased has given two columns for species in the taxonomy file (not a problem but that's why there is Species and Species2)
colnames(tax_table(elotuf.tax)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Species2")

# check how the phyloseq object looks (should have 930 taxa and 95 samples)
elotuf.tax

# merging OTUs by order

elotuf_order <- elotuf.tax %>%
  tax_glom(taxrank = "Order") %>%                     # agglomerate at class level
  #transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Order)  

elotuf_order2 <- subset(elotuf_order, select = c(2,3,5,6,8,14))

length(unique(elotuf_order2[,"Order"]))
# 93 orders - too many for plotting ->

# converting to wide for selecting orders with rel.abund. >1%
elotuf_order.wide <- pivot_wider(elotuf_order2, names_from = Order, values_from = Abundance)

fcolsums <- colSums(elotuf_order.wide[5:97])
frel <- as.data.frame(fcolsums/sum(fcolsums))
names(frel)[1] <- "rel.abundance"
frel_sub <- subset(frel, rel.abundance > 0.01)
ford_over1 <- row.names(frel_sub)

elotuf_order_other <- elotuf_order2 %>%
  mutate(Order_mod = ifelse(Order %in% ford_over1, Order, "other")) %>%
  mutate(Rel.Abundance = Abundance/sum(Abundance)) %>%
  filter(Abundance !=0)

length(unique(elotuf_order_other[,"Order_mod"]))
#count(distinct(elotuf_order_other, Order_mod))
# 15

#fpal <- c("olivedrab3", "olivedrab4", "orange", "yellow", "lightpink3", "hotpink3", "hotpink4", "indianred", "lightpink4")
  
elotuf_order_other$treat <- as.factor(elotuf_order_other$treat)
levels(elotuf_order_other$treat)
levels(elotuf_order_other$treat) <- c("open", "bagged")

elotuf_order_other$site <- as.factor(elotuf_order_other$site)
levels(elotuf_order_other$site)
levels(elotuf_order_other$site) <- c("A", "B", "C", "D", "E", "F")

# removing o__ etc from the order names
elotuf_order_other$Order_mod <- gsub("o__", "",elotuf_order_other$Order_mod)
elotuf_order_other$Order_mod <- gsub("k__", "",elotuf_order_other$Order_mod)
elotuf_order_other$Order_mod <- gsub("p__", "",elotuf_order_other$Order_mod)
elotuf_order_other$Order_mod <- gsub("_unclassified", " uncl.",elotuf_order_other$Order_mod)

elotuf_order_plot <- ggplot(elotuf_order_other, aes(x = site, y = Rel.Abundance, fill = Order_mod, colour=Order_mod)) + 
  geom_bar(stat = "identity", position="fill") +
  #scale_fill_manual(values = bpal, limits=force, breaks = c(ord_over1, "other"), name="Order") +
  scale_fill_manual(values = bpal, limits=force, name="Order") +
  scale_colour_manual(values = bpal, limits=force, guide="none") +
  ylab("Relative abundance") +
  theme_cowplot(12) + 
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  facet_grid(. ~ honeybees + treat, scales="free_x", space="free")

elotuf_order_plot

elotuf_order_plot
pdf("plots/FigS2b_elotuf_taxplot.pdf") 
print(elotuf_order_plot)
dev.off()

