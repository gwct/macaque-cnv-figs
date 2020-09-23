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

##############################
avgInds <- function(df){
  cnv_sums = c()
  del_len_sums = c()
  dup_len_sums = c()
  max_cnvs = 0
  min_cnvs = 99999
  for(i in levels(as.factor(df$Individual))){
    tmp = subset(df, Individual == i)
    num_cnvs = length(tmp[,1])
    tmp_dels = subset(tmp, Type == "<DEL>")
    tmp_dups = subset(tmp, Type == "<DUP>")
    cnv_sums = c(cnv_sums, num_cnvs)
    if(num_cnvs > max_cnvs){
      max_cnvs = num_cnvs
    }
    
    if(num_cnvs < min_cnvs){
      min_cnvs = num_cnvs
    }
    
    del_len_sums = c(del_len_sums, sum(tmp_dels$Length))
    dup_len_sums = c(dup_len_sums, sum(tmp_dups$Length))
  }
  means = data.frame("num.cnvs"=mean(cnv_sums),
                     "max.cnvs"=max_cnvs,
                     "min.cnvs"=min_cnvs,
                     "num.del.bases"=mean(del_len_sums),
                     "num.dup.bases"=mean(dup_len_sums))
  return(means)
}

##############################

svCount <- function(df, pad){
  total = length(df$Type)
  dels = length(df$Type[df$Type=="<DEL>"])
  dels_p = signif(dels / total, digits=4)
  dups = length(df$Type[df$Type=="<DUP>"])
  dups_p = signif(dups / total, digits=4)
  
  del_sens = 0.88
  del_corrected = round(dels / del_sens)
  del_missed = del_corrected - dels
  
  dup_sens = 0.65
  dup_corrected = round(dups / dup_sens)
  dup_missed = dup_corrected - dups
  # Sensitivity estimates from Sudmant et al.
  
  cat(spacedOut("Total:",                   pad, total), "\n")
  cat(spacedOut("Deletions:",               pad, paste(dels, " (", dels_p, ")", sep="")), "\n")
  cat(spacedOut("Deletion sensitivity:",    pad, del_sens), "\n")
  cat(spacedOut("Deletions missed:",        pad, del_missed), "\n")
  cat(spacedOut("Deletions corrected:",     pad, del_corrected), "\n")
  
  cat(spacedOut("Duplications:",            pad, paste(dups, " (", dups_p, ")", sep="")), "\n")
  cat(spacedOut("Duplication sensitivity:", pad, dup_sens), "\n")
  cat(spacedOut("Duplication missed:",      pad, dup_missed), "\n")
  cat(spacedOut("Duplication corrected:",   pad, dup_corrected), "\n")
  
  bases = sum(df$Length)
  avg_bases = mean(df$Length)
  min_bases = min(df$Length)
  max_bases = max(df$Length)
  
  del_bases = sum(df$Length[df$Type=="<DEL>"])
  avg_del_bases = mean(df$Length[df$Type=="<DEL>"])
  del_bases_p = signif(del_bases / bases, digits=4)
  min_del_bases = min(df$Length[df$Type=="<DEL>"])
  max_del_bases = max(df$Length[df$Type=="<DEL>"])
  
  dup_bases = sum(df$Length[df$Type=="<DUP>"])
  avg_dup_bases = mean(df$Length[df$Type=="<DUP>"])
  dup_bases_p = signif(dup_bases / bases, digits=4)
  min_dup_bases = min(df$Length[df$Type=="<DUP>"])
  max_dup_bases = max(df$Length[df$Type=="<DUP>"])
  
  cat(spacedOut("Bases:",                   pad, bases), "\n")
  cat(spacedOut("Deleted bases:",           pad, paste(del_bases, " (", del_bases_p, ")", sep="")), "\n")
  cat(spacedOut("Duplicated bases:",        pad, paste(dup_bases, " (", dup_bases_p, ")", sep="")), "\n")
  
  cat(spacedOut("Avg. length:",             pad, avg_bases), "\n")
  cat(spacedOut("Min. length:",             pad, min_bases), "\n")
  cat(spacedOut("Max. length:",             pad, max_bases), "\n")
  
  cat(spacedOut("Avg. deletion length:",    pad, paste(avg_del_bases, "\n")))
  cat(spacedOut("Min. deletion length:",    pad, paste(min_del_bases, "\n")))
  cat(spacedOut("Max. deletion length:",    pad, paste(max_del_bases, "\n")))
  
  cat(spacedOut("Avg. duplication length:", pad, paste(avg_dup_bases, "\n")))
  cat(spacedOut("Min. duplication length:", pad, paste(min_dup_bases, "\n")))
  cat(spacedOut("Max. duplication length:", pad, paste(max_dup_bases, "\n")))
}
# Counts stuff for a CNV data frame

##############################

geneCount <- function(cnv_df, gene_df){

  cat(spacedOut("\nCNVs that overlap at least 1 gene:",             pad, length(subset(cnv_df, Genes.total > 0)[,1])), "\n")
  cat(spacedOut("Avg. genes per CNV:",                              pad, mean(cnv_df$Genes.total)), "\n")
  cat(spacedOut("Avg. genes per CNV, given it overlaps 1:",         pad, mean(subset(cnv_df, Genes.total > 0)$Genes.total)), "\n")
  
  cat(spacedOut("CNVs that delete at least 1 gene:",                pad, length(subset(cnv_df, Genes.del > 0)[,1])), "\n")
  cat(spacedOut("Avg. genes per deletion:",                         pad, mean(cnv_df$Genes.del)), "\n")
  cat(spacedOut("Avg. genes per deletion, given it overlaps 1:",    pad, mean(subset(cnv_df, Genes.del > 0)$Genes.del)), "\n")
  
  cat(spacedOut("CNVs that duplicate at least 1 gene:",             pad, length(subset(cnv_df, Genes.dup > 0)[,1])), "\n")
  cat(spacedOut("Avg. genes per duplication:",                      pad, mean(cnv_df$Genes.dup)), "\n")
  cat(spacedOut("Avg. genes per duplication, given it overlaps 1:", pad, mean(subset(cnv_df, Genes.dup > 0)$Genes.dup)), "\n")
  
  genes = subset(gene_df, Feature.type == 'gene')
  cat(spacedOut("Total genes overlapped:",                          pad, sum(genes$Total)), "\n")
  cat(spacedOut("Genes overlapped by at least one CNV:",            pad, length(subset(genes, At.least.one == TRUE)[,1])), "\n")
  
  cat(spacedOut("Total genes deleted:",                             pad, sum(genes$Num.del)), "\n")
  cat(spacedOut("Genes deleted by at least one CNV:",               pad, length(subset(genes, At.least.one.del == TRUE)[,1])), "\n")
  
  cat(spacedOut("Total genes duplicated:",                          pad, sum(genes$Num.dup)), "\n")
  cat(spacedOut("Genes duplicated by at least one CNV:",            pad, length(subset(genes, At.least.one.dup == TRUE)[,1])), "\n")
  
  cat(spacedOut("Genes deleted / genes duplicated:",                pad, sum(genes$Num.del) / sum(genes$Num.dup)), "\n")
  
  ############
  
  cat(spacedOut("\nCNVs that overlap at least 1 transcript:",             pad, length(subset(cnv_df, Transcripts.total > 0)[,1])), "\n")
  cat(spacedOut("Avg. transcripts per CNV:",                              pad, mean(cnv_df$Transcripts.total)), "\n")
  cat(spacedOut("Avg. transcripts per CNV, given it overlaps 1:",         pad, mean(subset(cnv_df, Transcripts.total > 0)$Transcripts.total)), "\n")
    
  cat(spacedOut("CNVs that delete at least 1 transcript:",                pad, length(subset(cnv_df, Transcripts.del > 0)[,1])), "\n")
  cat(spacedOut("Avg. transcripts per deletion:",                         pad, mean(cnv_df$Transcripts.del)), "\n")
  cat(spacedOut("Avg. transcripts per deletion, given it overlaps 1:",    pad, mean(subset(cnv_df, Transcripts.del > 0)$Transcripts.del)), "\n")
  
  cat(spacedOut("CNVs that duplicate at least 1 transcript:",             pad, length(subset(cnv_df, Transcripts.dup > 0)[,1])), "\n")
  cat(spacedOut("Avg. transcripts per duplication:",                      pad, mean(cnv_df$Transcripts.dup)), "\n")
  cat(spacedOut("Avg. transcripts per duplication, given it overlaps 1:", pad, mean(subset(cnv_df, Transcripts.dup > 0)$Transcripts.dup)), "\n")
  
  transcripts = subset(gene_df, Feature.type == 'transcript')
  cat(spacedOut("Total transcripts overlapped:",                          pad, sum(transcripts$Total)), "\n")
  cat(spacedOut("Transcripts overlapped by at least one CNV:",            pad, length(subset(transcripts, At.least.one == TRUE)[,1])), "\n")
  
  cat(spacedOut("Total transcripts deleted:",                             pad, sum(transcripts$Num.del)), "\n")
  cat(spacedOut("Transcripts deleted by at least one CNV:",               pad, length(subset(transcripts, At.least.one.del == TRUE)[,1])), "\n")
  
  cat(spacedOut("Total transcripts duplicated:",                          pad, sum(transcripts$Num.dup)), "\n")
  cat(spacedOut("Transcripts duplicated by at least one CNV:",            pad, length(subset(transcripts, At.least.one.dup == TRUE)[,1])), "\n")
  
  cat(spacedOut("Transcripts deleted / transcripts duplicated:",          pad, sum(transcripts$Num.del) / sum(transcripts$Num.dup)), "\n")
  
  ############
  
  cat(spacedOut("\nCNVs that overlap at least 1 exon:",             pad, length(subset(cnv_df, Exons.total > 0)[,1])), "\n")
  cat(spacedOut("Avg. exons per CNV:",                              pad, mean(cnv_df$Exons.total)), "\n")
  cat(spacedOut("Avg. exons per CNV, given it overlaps 1:",         pad, mean(subset(cnv_df, Exons.total > 0)$Exons.total)), "\n")
  
  cat(spacedOut("CNVs that delete at least 1 exon:",                pad, length(subset(cnv_df, Exons.del > 0)[,1])), "\n")
  cat(spacedOut("Avg. exons per deletion:",                         pad, mean(cnv_df$Exons.del)), "\n")
  cat(spacedOut("Avg. exons per deletion, given it overlaps 1:",    pad, mean(subset(cnv_df, Exons.del > 0)$Exons.del)), "\n")
  
  cat(spacedOut("CNVs that duplicate at least 1 exon:",             pad, length(subset(cnv_df, Exons.dup > 0)[,1])), "\n")
  cat(spacedOut("Avg. exons per duplication:",                      pad, mean(cnv_df$Exons.dup)), "\n")
  cat(spacedOut("Avg. exons per duplication, given it overlaps 1:", pad, mean(subset(cnv_df, Exons.dup > 0)$Exons.dup)), "\n")
  
  exons = subset(gene_df, Feature.type == 'exon')
  cat(spacedOut("Total exons overlapped:",                          pad, sum(exons$Total)), "\n")
  cat(spacedOut("Exons overlapped by at least one CNV:",            pad, length(subset(exons, At.least.one == TRUE)[,1])), "\n")
  
  cat(spacedOut("Total exons deleted:",                             pad, sum(exons$Num.del)), "\n")
  cat(spacedOut("Exons deleted by at least one CNV:",               pad, length(subset(exons, At.least.one.del == TRUE)[,1])), "\n")
  
  cat(spacedOut("Total exons duplicated:",                          pad, sum(exons$Num.dup)), "\n")
  cat(spacedOut("Exons duplicated by at least one CNV:",            pad, length(subset(exons, At.least.one.dup == TRUE)[,1])), "\n")
  
  cat(spacedOut("Exons deleted / exons duplicated:",                pad, sum(exons$Num.del) / sum(exons$Num.dup)), "\n")
}
# Reports gene overlaps

############################################################
# Input info

readdata = T
readonly = F
savelog = T
filterdata = T
rm_alus = T
# Run options

minlen = F
maxlen = 100000
# CNV length cutoffs

cat("Input info:\n")
cat("readdata   = ", readdata, "\n")
cat("readonly   = ", readonly, "\n")
cat("filterdata = ", filterdata, "\n")
cat("rm_alus    = ", rm_alus, "\n")
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
if(rm_alus){
  logfile = paste(logfile, "-rmalus", sep="")
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
pad = 55

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
  
  if(rm_alus){
    mq_events = subset(mq_events, Length < 275 | Length > 325)
    hu_events = subset(hu_events, Length < 275 | Length > 325)
  }
  # For Alu stuff
  
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
  
  # if(rm_alus){
  #   mq_events = subset(mq_events, Length < 275 | Length > 325)
  #   mq_events_focal = subset(mq_events_focal, Length < 275 | Length > 325)
  #   mq_alleles = subset(mq_alleles, Length < 275 | Length > 325)
  #   mq_alleles_focal = subset(mq_alleles_focal, Length < 275 | Length > 325)
  #   mq_denovo = subset(mq_denovo, Length < 275 | Length > 325)
  #   
  #   hu_events = subset(hu_events, Length < 275 | Length > 325)
  #   hu_events_focal = subset(hu_events_focal, Length < 275 | Length > 325)
  #   hu_alleles = subset(hu_alleles, Length < 275 | Length > 325)
  #   hu_alleles_focal = subset(hu_alleles_focal, Length < 275 | Length > 325)
  # }
  # # For Alu stuff
  
  genes = readGenes(maxlen=maxlen)
  hu_cnv_genes = genes[[1]]; hu_genes = genes[[2]]; mq_cnv_genes = genes[[3]]; mq_genes = genes[[4]];
  # Read the gene data from filtered CNVs for both species.
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
ind_means = avgInds(mq_events)
cat(spacedOut("Average CNVs per individual:",               pad, ind_means$num.cnvs), "(", ind_means$min.cnvs, "-", ind_means$max.cnvs, ")", "\n")
cat(spacedOut("Average bases deleted per individual:",      pad, ind_means$num.del.bases), "\n")
cat(spacedOut("Average bases duplicated per individual:",   pad, ind_means$num.dup.bases), "\n")
cat(spacedOut("Number of CNVs present in all individuals:", pad, 1755), "\n")
# Report average SVs per individual and number in all individuals.
# 1755 is hard coded because they are filtered out. Change lib/filter_svs.r to get.

svCount(mq_events, pad)

#cat("\nFOCAL MACAQUE EVENTS:\n")
#svCount(mq_events_focal, pad)

cat("\nALL MACAQUE ALLELES:\n")
svCount(mq_alleles, pad)
geneCount(mq_cnv_genes, mq_genes)

cat("\nMACAQUE ALLELE GENE OVERLAPS FOR CNVs UNDER 25kb")
mq_genes_25kb = read.csv("../data/macaque-cnv-gene-overlaps-25kb.csv", header=TRUE)
mq_genes_25kb = parseCNVDF(mq_genes_25kb)
cat(spacedOut("\nCNVs under 25kb that overlap at least 1 exon:",  pad, length(subset(mq_genes_25kb, Exons.total > 0)[,1])), "\n")
cat(spacedOut("Avg. exons per CNV under 25kb:",                   pad, mean(subset(mq_genes_25kb, Exons.total > 0)$Exons.total)), "\n")
# Report gene overlaps

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
geneCount(hu_cnv_genes, hu_genes)

cat("\nHUMAN ALLELE GENE OVERLAPS FOR CNVs UNDER 25kb")
hu_genes_25kb = read.csv("../data/brandler-cnv-gene-overlaps-25kb.csv", header=TRUE)
hu_genes_25kb = parseCNVDF(hu_genes_25kb)
cat(spacedOut("\nCNVs under 25kb that overlap at least 1 exon:",  pad, length(subset(hu_genes_25kb, Exons.total > 0)[,1])), "\n")
cat(spacedOut("Avg. exons per CNV under 25kb:",                   pad, mean(subset(hu_genes_25kb, Exons.total > 0)$Exons.total)), "\n")

#cat("\nFOCAL HUMAN ALLELES:\n")
#svCount(hu_alleles_focal, pad)

cat("\n----------\n")

# The counts for various subsets of data
######################
# Chi-squared test for exons overlapping CNVs between macaques and humans
cat(" -> Chi-squared test for macaque exons in SVs vs human exons in SVs...\n")
mq_exons = subset(mq_genes, Feature.type == "exon")
hu_exons = subset(hu_genes, Feature.type == "exon")

mq_del_exons = sum(mq_exons$Num.del)
mq_dup_exons = sum(mq_exons$Num.dup)
hu_del_exons = sum(hu_exons$Num.del)
hu_dup_exons = sum(hu_exons$Num.dup)

exon_counts = data.frame("Species"=c("Macaque","Human"),
                         "num.del"=c(mq_del_exons,hu_del_exons),
                         "num.dup"=c(mq_dup_exons,hu_dup_exons))
exon_counts_t = t(exon_counts[,2:ncol(exon_counts)])
colnames(exon_counts_t) <- exon_counts[,1]
exons_chi = chisq.test(exon_counts_t)
print(exon_counts_t)
print(exons_chi)
cat("\n----------\n")
# Chi-squared test for exons overlapping CNVs between macaques and humans
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

cat("\nMACAQUE CAFE GENES:")

mq_changes = cafe_data$genes.changed[cafe_data$Node.ID=="Mmula"]
mq_gains = cafe_data$Genes.gained[cafe_data$Node.ID=="Mmula"]
mq_losses = cafe_data$Genes.lost[cafe_data$Node.ID=="Mmula"]
mq_gains_p = cafe_data$perc.genes.gained[cafe_data$Node.ID=="Mmula"]
mq_losses_p = cafe_data$perc.genes.lost[cafe_data$Node.ID=="Mmula"]

cat(spacedOut("Total changes:", pad, mq_changes), "\n")
cat(spacedOut("Losses:",        pad, paste(mq_losses, " (", mq_losses_p, ")", sep="")), "\n")
cat(spacedOut("Gains:",         pad, paste(mq_gains, " (", mq_gains_p, ")", sep="")), "\n")
# Macaque CAFE counts

cat("\n----------\n")

cat("\nHUMAN CAFE GENES:\n")

hu_changes = cafe_data$genes.changed[cafe_data$Node.ID=="Hsapi"]
hu_gains = cafe_data$Genes.gained[cafe_data$Node.ID=="Hsapi"]
hu_losses = cafe_data$Genes.lost[cafe_data$Node.ID=="Hsapi"]
hu_gains_p = cafe_data$perc.genes.gained[cafe_data$Node.ID=="Hsapi"]
hu_losses_p = cafe_data$perc.genes.lost[cafe_data$Node.ID=="Hsapi"]

cat(spacedOut("Total changes:", pad, hu_changes), "\n")
cat(spacedOut("Losses:",        pad, paste(hu_losses, " (", hu_losses_p, ")", sep="")), "\n")
cat(spacedOut("Gains:",         pad, paste(hu_gains, " (", hu_gains_p, ")", sep="")), "\n")
# Human CAFE counts 

cat("\n----------\n")

mq_gene_counts = subset(mq_genes, Feature.type == "gene")
hu_gene_counts = subset(hu_genes, Feature.type == "gene")

mq_del_genes = sum(mq_gene_counts$Num.del)
mq_dup_genes = sum(mq_gene_counts$Num.dup)
hu_del_genes = sum(hu_gene_counts$Num.del)
hu_dup_genes = sum(hu_gene_counts$Num.dup)

cat(" -> Chi-squared test for macaque SV genes (ALL) vs human SV genes (ALL)...\n")
mq_hu_all_counts = data.frame("Species"=c("Macaque SV all","Human SV all"),
                              "num.del"=c(mq_del_genes,hu_del_genes),
                              "num.dup"=c(mq_dup_genes,hu_dup_genes))
mq_hu_all_counts_t = t(mq_hu_all_counts[,2:ncol(mq_hu_all_counts)])
colnames(mq_hu_all_counts_t) <- mq_hu_all_counts[,1]
mq_hu_all_chi = chisq.test(mq_hu_all_counts_t)
print(mq_hu_all_counts_t)
print(mq_hu_all_chi)
# Macaque SV ALL vs Human SV ALL genes 

cat("\n----------\n")

mq_transcript_counts = subset(mq_genes, Feature.type == "transcript")
hu_transcript_counts = subset(hu_genes, Feature.type == "transcript")

mq_del_transcripts = sum(mq_transcript_counts$Num.del)
mq_dup_transcripts = sum(mq_transcript_counts$Num.dup)
hu_del_transcripts = sum(hu_transcript_counts$Num.del)
hu_dup_transcripts = sum(hu_transcript_counts$Num.dup)

cat(" -> Chi-squared test for macaque SV transcripts (ALL) vs macaque CAFE genes...\n")
mq_sv_all_counts = data.frame("Species"=c("Macaque SV all","Macaque CAFE"),
                         "num.del"=c(mq_del_transcripts,mq_losses),
                         "num.dup"=c(mq_dup_transcripts,mq_gains))
mq_sv_all_counts_t = t(mq_sv_all_counts[,2:ncol(mq_sv_all_counts)])
colnames(mq_sv_all_counts_t) <- mq_sv_all_counts[,1]
mq_sv_all_chi = chisq.test(mq_sv_all_counts_t)
print(mq_sv_all_counts_t)
print(mq_sv_all_chi)
# Macaque CAFE vs SV ALL genes test

cat("\n----------\n")

cat(" -> Chi-squared test for macaque SV transcripts (CAFE) vs macaque CAFE genes...\n")
cat("Reading macaque CAFE genes\n")
mq_cafe = read.csv("../data/macaque-cafe-genes.csv")
mq_cafe_dels = sum(mq_cafe$Num.del)
mq_cafe_dups = sum(mq_cafe$Num.dup)


mq_cafe_total = mq_cafe_dels + mq_cafe_dups
mq_cafe_del_p = mq_cafe_dels / mq_cafe_total
mq_cafe_dup_p = mq_cafe_dups / mq_cafe_total

mq_sv_cafe_counts = data.frame("Species"=c("Macaque SV CAFE","Macaque CAFE"),
                              "num.del"=c(mq_cafe_dels,mq_losses),
                              "num.dup"=c(mq_cafe_dups,mq_gains))
mq_sv_cafe_counts_t = t(mq_sv_cafe_counts[,2:ncol(mq_sv_cafe_counts)])
colnames(mq_sv_cafe_counts_t) <- mq_sv_cafe_counts[,1]
mq_sv_cafe_chi = chisq.test(mq_sv_cafe_counts_t)
print(mq_sv_cafe_counts_t)
print(mq_sv_cafe_chi)
# Macaque CAFE vs SV CAFE genes test
# Chi-squared tests for CAFE data
######################

if(savelog){
  sink()
}

cat("Done!")


mq_exon_ratio = matrix(c(441,90,3313,301), nrow=2, dimnames=list(c("total", "exon"), c("dup","del")))
mq_exon_chi = chisq.test(mq_exon_ratio)
print(mq_exon_chi)

hu_exon_ratio = matrix(c(1867,451,15746,2419), nrow=2, dimnames=list(c("total", "exon"), c("dup","del")))
hu_exon_chi = chisq.test(hu_exon_ratio)
print(hu_exon_chi)
