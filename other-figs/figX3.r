############################################################
# For macaque paper, 08.19
# Checks some distributions by SV.freq
# Gregg Thomas
############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(ggplot2)
library(reshape2)
library(cowplot)
library(tidyverse)

source("../lib/read_svs.r")
source("../lib/filter_svs.r")
source("../lib/subset_svs.r")
source("../lib/design.r")

cat("----------\n")

############################################################
# Input info

readdata = T
readonly = F
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
# Some info printed to the screen

######################
# Read the data
cat("----------\n")
pad = 55

if(readdata){
  sv_list = readSVs()
  if(filterdata){
    sv_list = filterSVs(sv_list, minlen, maxlen, freq_filter=1)
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
  
  genes = readGenes()
  hu_cnv_genes = genes[[1]]; hu_genes = genes[[2]]; mq_cnv_genes = genes[[3]]; mq_genes = genes[[4]];
  # Read the gene data from filtered CNVs for both species.
}

if(readonly){
  stop("Read only mode -- exiting.")
}
# Read the data
######################
# Length distributions grouped by SV.freq

mq_alleles_5k = subset(mq_alleles, Length <= 5000)
mq_alleles_5k$SV.freq = as.factor(mq_alleles_5k$SV.freq)

len_p = ggplot(mq_alleles_5k, aes(x=Length, fill=SV.freq, color=SV.freq)) + 
  geom_histogram(alpha=0.8, position="identity", bins=50) + facet_wrap(~SV.freq) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="CNV length", y="# CNVs") +
  bartheme() +
  theme(legend.position="none")

print(len_p)
ggsave(filename="figX3A.png", len_p, width=16, height=12, units="in")
# Length distributions grouped by SV.freq
######################
# Length distributions grouped by SV.freq, no SV.freq 1

mq_alleles_5k_no1 = subset(mq_alleles_5k, SV.freq != 1)

len_p = ggplot(mq_alleles_5k_no1, aes(x=Length, fill=SV.freq, color=SV.freq)) + 
  geom_histogram(alpha=0.8, position="identity", bins=50) + facet_wrap(~SV.freq) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="CNV length", y="# CNVs") +
  bartheme() +
  theme(legend.position="none")

print(len_p)
ggsave(filename="figX3B.png", len_p, width=16, height=12, units="in")
# Length distributions grouped by SV.freq, no SV.freq 1
######################
# SV.freq distribution

col=corecol(numcol=1, info=TRUE)

freq_p = ggplot(mq_alleles, aes(x=SV.freq)) + 
  geom_histogram(alpha=0.8, position="identity", bins=32, color="#333333", fill=col) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="CNV frequency", y="# CNVs") +
  bartheme()

print(freq_p)
ggsave(filename="figX3C.png", freq_p, width=6, height=4, units="in")

######################

mq_alleles_freq = subset(mq_alleles, SV.freq < 0.95)

freq_p = ggplot(mq_alleles_freq, aes(x=SV.freq)) + 
  geom_histogram(alpha=0.8, position="identity", bins=30, color="#333333", fill=col) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="CNV frequency", y="# CNVs") +
  bartheme()

print(freq_p)
ggsave(filename="figX3D.png", freq_p, width=6, height=4, units="in")
# SV.freq distribution
######################
# SNP allele freqs by SV freq
cat("----------\nReading CNV SNP count data...\n")
snp_freqs = read.csv("../data/cnv-snp-allele-freqs.csv", header=TRUE)
snp_freqs = subset(snp_freqs, SNP.bin != 0)

snp_freqs_all = subset(snp_freqs, SV.freq=="All")
snp_freqs = subset(snp_freqs, SV.freq != "All")

snp_p = ggplot(snp_freqs, aes(x=SNP.freq, y=SNP.count, fill=SV.freq, color=SV.freq)) + 
  geom_bar(stat="identity", alpha=0.8) + facet_wrap(~SV.freq) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="SNP frequency", y="# SNPs") +
  bartheme() +
  theme(legend.position="none")

print(snp_p)
ggsave(filename="figX3E.png", snp_p, width=16, height=12, units="in")
# SNP allele freqs by SV freq
######################
# SNP allele freqs for all SVs
col=corecol(numcol=1, offset=1, info=TRUE)

snp_p_all = ggplot(snp_freqs_all, aes(x=SNP.freq, y=SNP.count)) + 
  geom_bar(stat="identity", alpha=0.8, fill=col, col="#333333") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="SNP frequency", y="# SNPs") +
  bartheme() +
  theme(legend.position="none")

print(snp_p_all)
ggsave(filename="figX3F.png", snp_p_all, width=6, height=4, units="in")
# SNP allele freqs for all SVs
######################
# SNPs by SV bin
cnv_snps = aggregate(snp_freqs$SNP.count, list(by=snp_freqs$SV.count), FUN=sum)
names(cnv_snps) = c("Num.cnvs", "Num.snps")
cnv_snps_all = subset(cnv_snps, Num.cnvs > 1000)
cnv_snps = subset(cnv_snps, Num.cnvs < 1000)

col=corecol(numcol=1, offset=6, info=TRUE)

cnv_snps_p = ggplot(cnv_snps, aes(x=Num.cnvs, y=Num.snps)) +
  geom_smooth(method="glm", fullrange=T, color="#333333", linetype="dashed", se=F) +
  geom_point(size=3, alpha=0.5) +
  geom_point(data=cnv_snps_all, aes(x=Num.cnvs, y=Num.snps), size=4, color=col) +
  xlab("# CNVs in bin") +
  ylab("# SNPs") +
  bartheme()

print(cnv_snps_p)
ggsave(filename="figX3G.png", cnv_snps_p, width=6, height=4, units="in")
# SNPs by SV bin
######################
# CNVs per individual
mq_events_filtered = subset(mq_events, SV.freq < 0.95)
cnv_inds = aggregate(mq_events_filtered$Individual, list(by=mq_events_filtered$Individual), FUN=length)
names(cnv_inds) = c("Individual", "Num.cnvs")
ind_mean = mean(cnv_inds$Num.cnvs)
cnv_inds = rbind(cnv_inds, c("All", 1755))
cnv_inds$Individual = as.factor(cnv_inds$Individual)
cnv_inds$Num.cnvs = as.numeric(cnv_inds$Num.cnvs)

blah = data.frame()
col=corecol(numcol=1, offset=2, info=TRUE)

cnv_ind_p = ggplot(cnv_inds, aes(x=Individual, y=Num.cnvs)) + 
  #geom_point(alpha=0.8, color=col, size=4) +
  geom_bar(stat="identity", alpha=0.8, fill=col, col="#333333") +
  geom_hline(yintercept=ind_mean, size=1, linetype="dashed", color="#333333") +
  #geom_text(data=blah, aes(6, round(ind_mean), label=paste("Mean CNVs per ind =", ind_mean)), vjust=-2, size=6) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="Individual", y="# CNVs") +
  #coord_flip() +
  bartheme() +
  theme(legend.position="none",
        axis.text.x=element_text(angle=45, hjust=1))#,
        #panel.grid.major = element_line(colour="#ececec", size=0.5))

print(cnv_ind_p)
ggsave(filename="figX3H.png", cnv_ind_p, width=8, height=4, units="in")
# CNVs per individual
######################
