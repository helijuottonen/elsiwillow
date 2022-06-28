# Heli Juottonen 2019-2022 heli.juottonen@alumni.helsinki.fi
# R scripts for Elsi Hietaranta's willow manuscript
#  https://github.com/helijuottonen/elsiwillow

# 0: Basic processing: formatting input files for later analyses, subsampling OTU tables
# bacterial and fungal OTU tables have been produced in mothur (see manuscript for details)

# Disclaimer: I'm a self-taught R user. I try my very best to make sure that the scripts work, but I do not 
# guarantee that they are correct, make sense or follow any kind of best practices. Use at your own risk!

# loading necessary packages
library(vegan) #v. 2.5-7

# reading in OTU tables (files as they come from mothur):
# bacteria
elotub <- read.table("data/pajut_bothruns.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.shared", header=T)
# fungi
elotuf <- read.table("data/pajubothx.fasta.ITS2.precluster.pick.pick.agc.shared", header=T)


# reading in the metadata file:
elmeta <- read.csv2("data/Elsi_paju_metadata.csv", header=T, row.names=1)

# converting otu tables into dataframes
elotub.df <- as.data.frame(elotub)
elotuf.df <- as.data.frame(elotuf)

# insert sample names in the OTU tables as row names (needed in later commands):
row.names(elotub.df) <- elotub.df$Group
row.names(elotuf.df) <- elotuf.df$Group

# remove unnecessary columns in the OTU table:
# done because the first three columns (=1:3) of the shared files contain unnecessary values
# subset and -c means we are keeping else except the first 3 columns
elotub.df <- subset(elotub.df, select = -c(1:3))
elotuf.df <- subset(elotuf.df, select = -c(1:3))

# making sure the otu table and metadata contain the same samples
# this is done based on row names so otu table and metadata file must have identical formats of row names

elmetab <- elmeta[row.names(elmeta) %in% row.names(elotub.df),]
elmetaf <- elmeta[row.names(elmeta) %in% row.names(elotuf.df),]

elotub.df2 <- elotub.df[row.names(elotub.df) %in% row.names(elmetab),]
elotuf.df2 <- elotuf.df[row.names(elotuf.df) %in% row.names(elmetaf),]

# ordering otu table and metadata based on row names to make sure the order is the same
elotub.df2 <- elotub.df2[order(row.names(elotub.df2)),]
elotuf.df2 <- elotuf.df2[order(row.names(elotuf.df2)),]
elmetab <- elmetab[order(row.names(elmetab)),]
elmetaf <- elmetaf[order(row.names(elmetaf)),]

# subsampling the data (rarefying):

# defining a function for subsampling based on medians:
# picking the same number of reads from all samples: median of read numbers (=row sums)
# except samples that have less reads than the median: all reads are taken
# note:for comparing alpha diversity among samples it would be better to have the exact same number of reads
# this is a compromise because some samples have so few reads
median_raref <- function(x) {
  x_sums <- rowSums(x)
  x_sums2 <- replace(x_sums, x_sums>(median(x_sums)), median(x_sums))
  set.seed(1732)
  x.r <- rrarefy(x, x_sums2)
  x.rdf <- as.data.frame(x.r)
  return(x.rdf)
}

# sumsampling to median number of reads:
elotub.r.all <- median_raref(elotub.df2)
elotuf.r.all <- median_raref(elotuf.df2)

# removing OTUs with less than 5 reads overall
elotub.r <- elotub.r.all[,colSums(elotub.r.all) > 5]
elotuf.r <- elotuf.r.all[,colSums(elotuf.r.all) > 5]

# converting to relative abundances to even out the effects of median rarefaction, removal of rare OTUs
elotub.st <- decostand(elotub.r, "normalize")
elotuf.st <- decostand(elotuf.r, "normalize")

# OTU tables, metadata ready for further analyses 


