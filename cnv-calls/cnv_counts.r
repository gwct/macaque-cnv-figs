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

logfile = "cnv-counts"
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
# Brandler's de novo counts
cat("Counting Brandler de novos...\n")
brandler_denovo = read.csv("../cnv-calls/brandler-denovo.csv", header=TRUE)

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
# CAFE data
cat("CAFE DATA\n")
cafe_data = read.csv("../cafe-data/cafe-traits.csv", header=TRUE)
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
# CAFE data
# Human CAFE counts 
######################
cat("\n----------\n")

if(savelog){
  sink()
}

cat("Done!")



