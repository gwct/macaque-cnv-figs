############################################################
# For macaque paper, 08.19
# Counts SV stats and runs statistical tests
# Gregg Thomas
############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(ggplot2)
library(reshape2)
library(cowplot)

source("../lib/read_svs.r")
source("../lib/filter_svs.r")
source("../lib/subset_svs.r")

cat("----------\n")

############################################################
# Functions
spacedOut <- function(str1, pad, str2){
  while(nchar(str1)<pad){
    str1 = paste(str1, " ", sep="")
  }
  return(paste(str1, str2, sep=""))
}
# Nice formatting for string output

svCount <- function(df, pad){
  total = length(df$Type)
  dels = length(df$Type[df$Type=="<DEL>"])
  dels_p = signif(dels / total, digits=4)
  dups = length(df$Type[df$Type=="<DUP>"])
  dups_p = signif(dups / total, digits=4)
  
  cat(spacedOut("Total:", pad, total), "\n")
  cat(spacedOut("Deletions:", pad, paste(dels, " (", dels_p, ")", sep="")), "\n")
  cat(spacedOut("Duplications:", pad, paste(dups, " (", dups_p, ")", sep="")), "\n")
  
  bases = sum(df$Length)
  avg_bases = mean(df$Length)
  del_bases = sum(df$Length[df$Type=="<DEL>"])
  avg_del_bases = mean(df$Length[df$Type=="<DEL>"])
  del_bases_p = signif(del_bases / bases, digits=4)
  dup_bases = sum(df$Length[df$Type=="<DUP>"])
  avg_dup_bases = mean(df$Length[df$Type=="<DUP>"])
  dup_bases_p = signif(dup_bases / bases, digits=4)
  
  cat(spacedOut("Bases:", pad, bases), "\n")
  cat(spacedOut("Deleted bases:", pad, paste(del_bases, " (", del_bases_p, ")", sep="")), "\n")
  cat(spacedOut("Duplicated bases:", pad, paste(dup_bases, " (", dup_bases_p, ")", sep="")), "\n")
  
  cat(spacedOut("Avg. length:", pad, avg_bases), "\n")
  cat(spacedOut("Avg. deletion length:", pad, paste(avg_del_bases, "\n")))
  cat(spacedOut("Avg. duplication length:", pad, paste(avg_dup_bases, "\n")))
  
  overlap = length(df$Num.genes[df$Num.genes!=0])
  genes = sum(df$Num.genes)
  avg_genes = mean(df$Num.genes)
  del_genes = sum(df$Num.genes[df$Type=="<DEL>"])
  avg_del_genes = mean(df$Num.genes[df$Type=="<DEL>"])
  del_genes_p = signif(del_genes / genes, digits=4)
  dup_genes = sum(df$Num.genes[df$Type=="<DUP>"])
  dup_genes_p = signif(dup_genes / genes, digits=4)
  avg_dup_genes = mean(df$Num.genes[df$Type=="<DUP>"])
  
  avg_genes_overlap = mean(df$Num.genes[df$Num.genes!=0])
  avg_genes_overlap_del = mean(df$Num.genes[df$Num.genes!=0 & df$Type=="<DEL>"])
  avg_genes_overlap_dup = mean(df$Num.genes[df$Num.genes!=0 & df$Type=="<DUP>"])
  
  
  cat(spacedOut("SVs overlapping genes:", pad, overlap), "\n")
  cat(spacedOut("Genes overlapped:", pad, genes), "\n")
  #cat(spacedOut("Genes string:", pad, genesstr), "\n")
  cat(spacedOut("Deleted genes:", pad, paste(del_genes, " (", del_genes_p, ")", sep="")), "\n")
  cat(spacedOut("Duplicated genes:", pad, paste(dup_genes, " (", dup_genes_p, ")", sep="")), "\n")
  
  cat(spacedOut("Avg. genes per SV:", pad, avg_genes), "\n")
  cat(spacedOut("Avg. genes per deletion:", pad, paste(avg_del_genes, "\n")))
  cat(spacedOut("Avg. genes per duplication:", pad, paste(avg_dup_genes, "\n")))
  
  cat(spacedOut("Avg. genes per overlap SV:", pad, avg_genes_overlap), "\n")
  cat(spacedOut("Avg. genes per overlap deletion:", pad, paste(avg_genes_overlap_del, "\n")))
  cat(spacedOut("Avg. genes per overlap duplication:", pad, paste(avg_genes_overlap_dup, "\n")))
}
# Counts stuff for a CNV data frame

############################################################
# Input info

readdata = T
readonly = F
savelog = T
filterdata = T
# Run options

minlen = F
maxlen = 100000
# CNV length cutoffs

cat("Input info:\n")
cat("readdata   = ", readdata, "\n")
cat("readonly   = ", readonly, "\n")
cat("filterdata = ", filterdata, "\n")
cat("minlen     = ", minlen, "\n")
cat("maxlen     = ", maxlen, "\n")
cat("savelog    = ", savelog, "\n")
# Some info printed to the screen

logfile = "sv-counts"
if(filterdata){
  logfile = paste(logfile, "-filtered", sep="")
}else{
  logfile = paste(logfile, "-unfiltered", sep="")
}
if(maxlen){
  logfile = paste(logfile, "-", maxlen, ".log", sep="")
}else{
  logfile = paste(logfile, "-all.log", sep="")
}

cat("logfile   -> ", logfile, "\n")
# The logfile

continue = readline("Review input vars above. Continue? (y/n) ")

if(!continue %in% c("Y", "y", "yes", "Yes", "YES")){
  stop("User chose to discontinue. Exiting.")
}
# Input info
######################
# Read the data
cat("----------\n")
cat("Starting counts\n")
pad = 40

if(savelog){
  sink(logfile)
}

if(readdata){
  sv_list = readSVs()
  if(filterdata){
    sv_list = filterSVs(sv_list, minlen, maxlen)
  }
  mq_events = sv_list[[1]]; hu_events = sv_list[[2]];
  # Read and filter data
  
  cat("----------\nSubsetting macaque data...\n")
  mqr = subsetSVs(mq_events)
  mq_events_focal = mqr[[2]];
  mq_alleles = mqr[[4]]; mq_alleles_focal = mqr[[5]];
  
  cat("----------\nSubsetting human data...\n")
  hur = subsetSVs(hu_events)
  hu_events_focal = hur[[2]];
  hu_alleles = hur[[4]]; hu_alleles_focal = hur[[5]];
  # Subset data
  
  mq_denovo = subset(mq_events, Denovo=="Y")
  # Get macaque de novos
}

if(readonly){
  if(savelog){
    sink()
  }
  stop("Read only mode -- exiting.")
}
# Read the data
######################
# The counts for various subsets of data. 
# EVENT = all CNVs in all individuals (ie a deletion at the same position in two individuals is counted as 2 events)
# ALLELE = all CNVs across individuals (ie a deletion at the same position in two individuals is counted as 1 allele)
cat("\nALL MACAQUE EVENTS:\n")
svCount(mq_events, pad)

#cat("\nFOCAL MACAQUE EVENTS:\n")
#svCount(mq_events_focal, pad)

cat("\nALL MACAQUE ALLELES:\n")
svCount(mq_alleles, pad)

#cat("\nFOCAL MACAQUE ALLELES:\n")
#svCount(mq_alleles_focal, pad)

cat("\nDE NOVO MACAQUE ALLELES:\n")
svCount(mq_denovo, pad)

cat("\n----------\n")

cat("\nALL HUMAN EVENTS:\n")
svCount(hu_events, pad)

#cat("\nFOCAL HUMAN EVENTS:\n")
#svCount(hu_events_focal, pad)

cat("\nALL HUMAN ALLELES:\n")
svCount(hu_alleles, pad)

#cat("\nFOCAL HUMAN ALLELES:\n")
#svCount(hu_alleles_focal, pad)

cat("\n----------\n")
# The counts for various subsets of data
######################
# Chi-squared test for genes overlapping CNVs between macaques and humans
cat(" -> Chi-squared test for macaque genes in SVs vs human genes in SVs...\n")
mq_del_genes = sum(mq_alleles$Num.genes[mq_alleles$Type=="<DEL>"])
mq_dup_genes = sum(mq_alleles$Num.genes[mq_alleles$Type=="<DUP>"])
hu_del_genes = sum(hu_alleles$Num.genes[hu_alleles$Type=="<DEL>"])
hu_dup_genes = sum(hu_alleles$Num.genes[hu_alleles$Type=="<DUP>"])

gene_counts = data.frame("Species"=c("Macaque","Human"),
                         "num.del"=c(mq_del_genes,hu_del_genes),
                         "num.dup"=c(mq_dup_genes,hu_dup_genes))
gene_counts_t = t(gene_counts[,2:ncol(gene_counts)])
colnames(gene_counts_t) <- gene_counts[,1]
genes_chi = chisq.test(gene_counts_t)
print(genes_chi)
# Macaque CAFE vs SV genes test

cat("\n----------\n")
# Chi-squared test for genes overlapping CNVs between macaques and humans
######################
# Brandler's de novo counts
cat("Counting Brandler de novos...\n")
brandler_denovo = read.csv("../data/brandler-denovo.csv", header=TRUE)

brandler_denovo = subset(brandler_denovo, SVTYPE %in% c("DEL","DUP"))
names(brandler_denovo)[3] = "Length"
names(brandler_denovo)[11] = "Type"
brandler_denovo$Type = as.character(brandler_denovo$Type)
brandler_denovo$Type[brandler_denovo$Type=="DEL"] = "<DEL>"
brandler_denovo$Type[brandler_denovo$Type=="DUP"] = "<DUP>"

cat("\nDE NOVO BRANDLER ALLELES:\n")
svCount(brandler_denovo, pad)

cat("\n----------\n")
# Brandler's de novo counts
######################
# Chi-squared tests for CAFE data
cat("CAFE DATA\n")
cafe_data = read.csv("../data/cafe-traits.csv", header=TRUE)
#cafe_data = subset(cafe_data, Node.type=="Tip")
cafe_data$genes.changed = cafe_data$Genes.gained + cafe_data$Genes.lost
cafe_data$perc.genes.gained = cafe_data$Genes.gained / cafe_data$genes.changed
cafe_data$perc.genes.lost = cafe_data$Genes.lost / cafe_data$genes.changed
# Read the CAFE data

cat("\nMACAQUE GENES:\n")

mq_changes = cafe_data$genes.changed[cafe_data$Node.ID=="Mmula"]
mq_gains = cafe_data$Genes.gained[cafe_data$Node.ID=="Mmula"]
mq_losses = cafe_data$Genes.lost[cafe_data$Node.ID=="Mmula"]
mq_gains_p = cafe_data$perc.genes.gained[cafe_data$Node.ID=="Mmula"]
mq_losses_p = cafe_data$perc.genes.lost[cafe_data$Node.ID=="Mmula"]

cat(spacedOut("Total changes:", pad, mq_changes), "\n")
cat(spacedOut("Losses:", pad, paste(mq_losses, " (", mq_losses_p, ")", sep="")), "\n")
cat(spacedOut("Gains:", pad, paste(mq_gains, " (", mq_gains_p, ")", sep="")), "\n")
# Macaque CAFE counts

cat("\nHUMAN GENES:\n")

hu_changes = cafe_data$genes.changed[cafe_data$Node.ID=="Hsapi"]
hu_gains = cafe_data$Genes.gained[cafe_data$Node.ID=="Hsapi"]
hu_losses = cafe_data$Genes.lost[cafe_data$Node.ID=="Hsapi"]
hu_gains_p = cafe_data$perc.genes.gained[cafe_data$Node.ID=="Hsapi"]
hu_losses_p = cafe_data$perc.genes.lost[cafe_data$Node.ID=="Hsapi"]

cat(spacedOut("Total changes:", pad, hu_changes), "\n")
cat(spacedOut("Losses:", pad, paste(hu_losses, " (", hu_losses_p, ")", sep="")), "\n")
cat(spacedOut("Gains:", pad, paste(hu_gains, " (", hu_gains_p, ")", sep="")), "\n")
# Human CAFE counts

cat(" -> Chi-squared test for macaque SV genes (ALL) vs human SV genes (ALL)...\n")
mq_hu_all_counts = data.frame("Species"=c("Macaque SV all","Human SV all"),
                              "num.del"=c(mq_del_genes,hu_del_genes),
                              "num.dup"=c(mq_dup_genes,hu_dup_genes))
mq_hu_all_counts_t = t(mq_hu_all_counts[,2:ncol(mq_hu_all_counts)])
colnames(mq_hu_all_counts_t) <- mq_hu_all_counts[,1]
mq_hu_all_chi = chisq.test(mq_hu_all_counts_t)
print(mq_hu_all_chi)
# Macaque SV ALL vs Human SV ALL genes test

cat(" -> Chi-squared test for macaque SV genes (ALL) vs macaque CAFE genes...\n")
mq_sv_all_counts = data.frame("Species"=c("Macaque SV all","Macaque CAFE"),
                         "num.del"=c(mq_del_genes,mq_losses),
                         "num.dup"=c(mq_dup_genes,mq_gains))
mq_sv_all_counts_t = t(mq_sv_all_counts[,2:ncol(mq_sv_all_counts)])
colnames(mq_sv_all_counts_t) <- mq_sv_all_counts[,1]
mq_sv_all_chi = chisq.test(mq_sv_all_counts_t)
print(mq_sv_all_chi)
# Macaque CAFE vs SV ALL genes test

cat(" -> Chi-squared test for macaque SV genes (CAFE) vs macaque CAFE genes...\n")
cat("Reading macaque CAFE genes\n")
mq_cafe = read.csv("../data/macaque-cafe-genes-filtered.csv")
mq_cafe_genes = subset(mq_cafe, SV.key %in% mq_alleles$SV.key)
mq_cafe_genes = subset(mq_cafe_genes, SV.type=="<DEL>" | SV.type=="<DUP>")

mq_cafe_dels = length(mq_cafe_genes$SV.type[mq_cafe_genes$SV.type=="<DEL>"])
mq_cafe_dups = length(mq_cafe_genes$SV.type[mq_cafe_genes$SV.type=="<DUP>"])
mq_cafe_total = mq_cafe_dels + mq_cafe_dups
mq_cafe_del_p = mq_cafe_dels / mq_cafe_total
mq_cafe_dup_p = mq_cafe_dups / mq_cafe_total

mq_sv_cafe_counts = data.frame("Species"=c("Macaque SV CAFE","Macaque CAFE"),
                              "num.del"=c(mq_cafe_dels,mq_losses),
                              "num.dup"=c(mq_cafe_dups,mq_gains))
mq_sv_cafe_counts_t = t(mq_sv_cafe_counts[,2:ncol(mq_sv_cafe_counts)])
colnames(mq_sv_cafe_counts_t) <- mq_sv_cafe_counts[,1]
mq_sv_cafe_chi = chisq.test(mq_sv_cafe_counts_t)
print(mq_sv_cafe_chi)
# Macaque CAFE vs SV CAFE genes test
# Chi-squared tests for CAFE data
######################

if(savelog){
  sink()
}

cat("Done!")