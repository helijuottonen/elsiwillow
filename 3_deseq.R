# Heli Juottonen 2019-2022 heli.juottonen@alumni.helsinki.fi
# R scripts for Elsi Hietaranta's willow manuscript
# https://github.com/helijuottonen/elsiwillow

# 3: Differentially abundant OTUs between treatments with DESeq2

# Disclaimer: I'm a self-taught R user. I try my very best to make sure that the scripts work, but I do not 
# guarantee that they are correct, make sense or follow any kind of best practices. Use at your own risk!

# RUN FIRST: 0_basic_processing.R

#library(vegan)
library("DESeq2") #v.1.28.1
library(ggplot2) #v.3.3.5
library(phyloseq) #1.32.0
library(cowplot) #v.1.1.1

# DEseq2 manual: https://bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf

# setting theme for plotting
theme_set(theme_bw())

################# BACTERIA

# removing bagged samples from metadata and otu table for honeybee test
elmetab_nb <- subset(elmetab, treat=="ctrl")
elotub_nb = elotub.df2[row.names(elotub.df2) %in% row.names(elmetab_nb),]

### DESeq2 starts

# transpose OTU table (deseq wants OTUs in rows)
elotub.t <- t(elotub.df2)
elotub_nb.t <- t(elotub_nb)

# creating deseq object
# here countData = transposed OTU table
# colData = metadata
# design = the column of metadata that contains the treatment
elotub.treat <- DESeqDataSetFromMatrix(countData= elotub.t, colData=elmetab, design= ~treat)
elotub.hb <- DESeqDataSetFromMatrix(countData= elotub_nb.t, colData=elmetab_nb, design= ~honeybees)

# removing OTUs with less x number of reads
# the purpose here is to remove small sporadic OTUs and keep OTU occurring consistently at least in one sample type
# 6 places x 8 replicates = 48 used as the minimum number of reads
keep.elotub.treat <- rowSums(counts(elotub.treat)) >=48
elotub.treat <- elotub.treat[keep.elotub.treat, ]

keep.elotub.hb <- rowSums(counts(elotub.hb)) >=48
elotub.hb <- elotub.hb[keep.elotub.hb, ]

# add treatment into the deseq object as a factor
elotub.treat$treat <- factor(elotub.treat$treat, levels = c("ctrl", "ps"))
elotub.hb$honeybees <- factor(elotub.hb$honeybees, levels = c("absent", "present"))

# run the deseq analysis
elotub.treat <- DESeq(elotub.treat)
elotub.hb <- DESeq(elotub.hb)

# write deseq results into an object
elotub.treat.res <- results(elotub.treat)
elotub.hb.res <- results(elotub.hb)

# check a summary of the results
summary(elotub.treat.res)
summary(elotub.hb.res)

# write the results into a file

write.csv2(elotub.treat.res, file = "results/bact_deseq_treat.csv")
write.csv2(elotub.hb.res, file = "results/bact_deseq_honeybees.csv")


####plotting results

#reading in results
btreat <- read.csv2("results/bact_deseq_treat.csv")
bhoney <- read.csv2("results/bact_deseq_honeybees.csv")

colnames(btreat)[1] <- "otu"
colnames(bhoney)[1] <- "otu"
rownames(btreat) <- btreat$otu
rownames(bhoney) <- bhoney$otu

# subsetting based on p.adj values
btreat2 <- subset(btreat, btreat$padj<0.05 & abs(log2FoldChange)>1)
bhoney2 <- subset(bhoney, bhoney$padj<0.05 & abs(log2FoldChange)>1)

# adding taxonomy for OTUs:
# exctrating taxonomy from the phyloseq object from the taxonomy script

taxfile.bact = "data/pajut_bothruns.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.cons.taxonomy.txt"

#importing the taxonomy tables into phyloseq format
mothur_tax.bact <- import_mothur(mothur_constaxonomy_file = taxfile.bact)

otutable.bact <- otu_table(as.matrix(elotub.df2), taxa_are_rows=FALSE)
metadata.bact <- sample_data(elmetab)

elotub <- merge_phyloseq(otutable.bact, metadata.bact, mothur_tax.bact)

# giving proper names to the columns of the taxonomy tables
colnames(tax_table(elotub)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

elotub.tax = as(tax_table(elotub), "matrix")
elotub.tax.df = as.data.frame(elotub.tax)

elotub.tax.treat <- elotub.tax.df[row.names(elotub.tax.df) %in% row.names(btreat2),]
elotub.tax.honey <- elotub.tax.df[row.names(elotub.tax.df) %in% row.names(bhoney2),]

btreat2$otutax <- paste(elotub.tax.treat$Phylum, ": ", elotub.tax.treat$Family,  " ", "(", btreat2$otu, ")",  sep="")
bhoney2$otutax <- paste(elotub.tax.honey$Phylum, ": ", elotub.tax.honey$Family,  " ", "(", bhoney2$otu, ")",  sep="")

btreat2$otutax <- gsub("_unclassified", " uncl.", btreat2$otutax)
bhoney2$otutax <- gsub("_unclassified", " uncl.", bhoney2$otutax)
bhoney2$otutax <- gsub("_\\(Subgroup_1\\)", " gr.1", bhoney2$otutax)

# setting the order of OTUs based on log fold change
btreat2$otutax = with(btreat2, reorder(otutax, log2FoldChange, median))
bhoney2$otutax = with(bhoney2, reorder(otutax, log2FoldChange, median))

btreat2$colorlfc <-ifelse(btreat2$log2FoldChange > 1, "indianred", "steelblue")
bhoney2$colorlfc <-ifelse(bhoney2$log2FoldChange > 1, "indianred", "steelblue")

#plots

btreat.plot <- ggplot(btreat2, aes(otutax, log2FoldChange, color=colorlfc)) +
  geom_hline(yintercept = 0) +
  geom_point(size=2) + 
  scale_color_identity() +
  theme_cowplot(11) +
  background_grid(color.major="grey90") +
  theme(axis.title.y = element_blank()) +
  ylim(-6,4) +
  coord_flip()

btreat.plot

bhoney.plot <- ggplot(bhoney2, aes(otutax, log2FoldChange, color=colorlfc)) +
  geom_hline(yintercept = 0) +
  geom_point(size=2) + 
  scale_color_identity() +
  theme_cowplot(11) +
  background_grid(color.major="grey90") +
  theme(axis.title.y = element_blank()) +
  ylim(-7,8) +
  coord_flip()

bhoney.plot

################# FUNGI

# removing bagged samples from metadata and otu table for honeybee test
elmetaf_nb <- subset(elmetaf, treat=="ctrl")
elotuf_nb = elotuf.df2[row.names(elotuf.df2) %in% row.names(elmetaf_nb),]


### deseq starts

# transpose OTU table (deseq wants OTUs in rows)
elotuf.t <- t(elotuf.df2)
elotuf_nb.t <- t(elotuf_nb)

# creating deseq object
# here countData = transposed OTU table
# colData = metadata
# design = the column of metadata that contains the treatment

elotuft.treat <- DESeqDataSetFromMatrix(countData= elotuf.t, colData=elmetaf, design= ~treat)
elotuft.hb <- DESeqDataSetFromMatrix(countData= elotuf_nb.t, colData=elmetaf_nb, design= ~honeybees)

# removing OTUs with less than 10 reads (this number can be changed according to what makes sense)
# the purpose here is to remove small sporadic OTUs and keep OTU occurring consistently at least in one sample type
keep.elotuft.treat <- rowSums(counts(elotuft.treat)) >=48
elotuft.treat <- elotuft.treat[keep.elotuft.treat, ]

keep.elotuft.hb <- rowSums(counts(elotuft.hb)) >=48
elotuft.hb <- elotuft.hb[keep.elotuft.hb, ]

# add treatment into the desed object as a factor
elotuft.treat$treat <- factor(elotuft.treat$treat, levels = c("ctrl", "ps"))
elotuft.hb$honeybees <- factor(elotuft.hb$honeybees, levels = c("absent", "present"))

# run the deseq analysis
elotuft.treat <- DESeq(elotuft.treat)
elotuft.hb <- DESeq(elotuft.hb)

# write deseq results into an object
elotuft.treat.res <- results(elotuft.treat)
elotuft.hb.res <- results(elotuft.hb)

# check a summary of the results
summary(elotuft.treat.res)
summary(elotuft.hb.res)

# write the results into a file
write.csv2(elotuft.treat.res, file = "results/fungi_deseq_treat_new.csv")
write.csv2(elotuft.hb.res, file = "results/fungi_deseq_honeybees_nb_new.csv")

### plotting

#reading in data
ftreat <- read.csv2("results/fungi_deseq_treat_new.csv")
fhoney <- read.csv2("results/fungi_deseq_honeybees_nb_new.csv")

colnames(ftreat)[1] <- "otu"
colnames(fhoney)[1] <- "otu"
rownames(ftreat) <- ftreat$otu
rownames(fhoney) <- fhoney$otu

# subsetting based on p.adj values
ftreat2 <- subset(ftreat, ftreat$padj<0.05 & abs(log2FoldChange)>1)
fhoney2 <- subset(fhoney, fhoney$padj<0.05 & abs(log2FoldChange)>1)

# adding taxonomy for OTUs:
# exctrating taxonomy from the phyloseq object from the taxonomy script

taxfile.fungi = "data/pajubothx.fasta.ITS2.precluster.pick.pick.agc.0.03.cons.taxonomy.txt"

#importing the taxonomy tables into phyloseq format
mothur_tax.fungi <- import_mothur(mothur_constaxonomy_file = taxfile.fungi)

otutable.fungi <- otu_table(as.matrix(elotuf.df2), taxa_are_rows=FALSE)
metadata.fungi <- sample_data(elmetaf)

elotuf <- merge_phyloseq(otutable.fungi, metadata.fungi, mothur_tax.fungi)

# giving proper names to the columns of the taxonomy tables
# NOTE 4.11.20: for some reason mothur/Unite databased has given two columns for species in the taxonomy file (not a problem but that's why there is Species and Species2)
colnames(tax_table(elotuf)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Species2")

elotuf.tax = as(tax_table(elotuf), "matrix")
elotuf.tax.df = as.data.frame(elotuf.tax)

elotuf.tax.treat <- elotuf.tax.df[row.names(elotuf.tax.df) %in% row.names(ftreat2),]
elotuf.tax.honey <- elotuf.tax.df[row.names(elotuf.tax.df) %in% row.names(fhoney2),]

ftreat2$otutax <- paste(elotuf.tax.treat$Class, ": ", elotuf.tax.treat$Genus,  " ", "(", ftreat2$otu, ")",  sep="")
fhoney2$otutax <- paste(elotuf.tax.honey$Class, ": ", elotuf.tax.honey$Genus,  " ", "(", fhoney2$otu, ")",  sep="")

# removing unnecessary things from taxonomic names
ftreat2$otutax <- gsub("p__", "", ftreat2$otutax)
ftreat2$otutax <- gsub("k__", "", ftreat2$otutax)
ftreat2$otutax <- gsub("f__", "", ftreat2$otutax)
ftreat2$otutax <- gsub("c__", "", ftreat2$otutax)
ftreat2$otutax <- gsub("o__", "", ftreat2$otutax)
ftreat2$otutax <- gsub("g__", "", ftreat2$otutax)
ftreat2$otutax <- gsub("_unclassified", " uncl.", ftreat2$otutax)
ftreat2$otutax <- gsub("unclassified_", "uncl. ", ftreat2$otutax)
ftreat2$otutax <- gsub("_fam_Incertae_sedis", " fam.inc.sed.", ftreat2$otutax)
ftreat2$otutax <- gsub("_cls_Incertae_sedis", " cl.inc.sed.", ftreat2$otutax)

fhoney2$otutax <- gsub("p__", "", fhoney2$otutax)
fhoney2$otutax <- gsub("k__", "", fhoney2$otutax)
fhoney2$otutax <- gsub("f__", "", fhoney2$otutax)
fhoney2$otutax <- gsub("c__", "", fhoney2$otutax)
fhoney2$otutax <- gsub("o__", "", fhoney2$otutax)
fhoney2$otutax <- gsub("g__", "", fhoney2$otutax)
fhoney2$otutax <- gsub("_unclassified", " uncl.", fhoney2$otutax)
fhoney2$otutax <- gsub("unclassified_", "uncl. ", fhoney2$otutax)
fhoney2$otutax <- gsub("_fam_Incertae_sedis", " fam.inc.sed.", fhoney2$otutax)
fhoney2$otutax <- gsub("_cls_Incertae_sedis", " cl.inc.sed.", fhoney2$otutax)

# setting the order of OTUs based on log fold change
ftreat2$otutax = with(ftreat2, reorder(otutax, log2FoldChange, median))
fhoney2$otutax = with(fhoney2, reorder(otutax, log2FoldChange, median))

ftreat2$colorlfc <-ifelse(ftreat2$log2FoldChange > 1, "indianred", "steelblue")
fhoney2$colorlfc <-ifelse(fhoney2$log2FoldChange > 1, "indianred", "steelblue")

#plots
ftreat.plot <- ggplot(ftreat2, aes(otutax, log2FoldChange, color=colorlfc)) +
  geom_hline(yintercept = 0) +
  geom_point(size=2) + 
  scale_color_identity() +
  theme_cowplot(11) +
  background_grid(color.major="grey90") +
  theme(axis.title.y = element_blank()) +
  ylim(-6,4) +
  coord_flip()

ftreat.plot

fhoney.plot <- ggplot(fhoney2, aes(otutax, log2FoldChange, color=colorlfc)) +
  geom_hline(yintercept = 0) +
  geom_point(size=2) + 
  scale_color_identity() +
  theme_cowplot(11) +
  background_grid(color.major="grey90") +
  theme(axis.title.y = element_blank()) +
  ylim(-7,8) +
  coord_flip()

fhoney.plot

#### combining plots

# bagging
treat.plots <- plot_grid(btreat.plot, ftreat.plot, labels=c("a)", "b)"), label_fontface="plain", ncol=1, rel_heights=c(0.7,2), align="v")

treat.plots
pdf("plots/Fig2_deseq_treat.pdf")
print(treat.plots)
dev.off()

#honeybees
hb.plots <- plot_grid(bhoney.plot, fhoney.plot, labels=c("a)", "b)"), label_fontface="plain", ncol=1, rel_heights=c(0.8,2), align="v")

hb.plots
pdf("plots/Fig1_deseq_honey.pdf")
print(hb.plots)
dev.off()

# note: text colours of the plots edited in Inkscape afterwards

