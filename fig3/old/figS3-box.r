############################################################
# For macaque paper, 08.19
# Makes length histograms for human and macaque SVs
# Gregg Thomas
############################################################

binSVs <- function(svs){
# Function to pre-bin SVs by length
  binned_svs = data.frame("bin"=seq(from=0,to=5000,by=100), count=0)
  for(i in 1:nrow(svs)) {
    row <- svs[i,]
    bin = floor(row$Length/100)*100
    binned_svs[binned_svs$bin==bin,]$count = binned_svs[binned_svs$bin==bin,]$count + 1
  }
  return(binned_svs) 
}

############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(ggplot2)
library(reshape2)
library(cowplot)
library(ggridges)
library(OneR)
library(ggforce)
library(ggbeeswarm)

source("../lib/read_svs.r")
source("../lib/filter_svs.r")
source("../lib/subset_svs.r")

cat("----------\n")

############################################################

savefiles = F
color_plots = T
rm_alus = F
# Run options

minlen = F
maxlen = 5000
# CNV length cutoffs

baseout = "fig3"
human_full = F
macaque_filter = T
figs2_opt = F
figs3_opt = F
if(figs2_opt){
  human_full = T
  baseout = paste(baseout, "_S2", sep="")
}else if(figs3_opt){
  macaque_filter = F
  baseout = paste(baseout, "_S3", sep="")
}

if(rm_alus){
  human_full = F
  macaque_filter = T
  baseout = paste(baseout, "_S4", sep="")
}
# Supplemental figure options

############################################################

sv_list = readSVs()
mq_events = sv_list[[1]]; hu_events = sv_list[[2]];
sv_list = filterSVs(sv_list, minlen, maxlen, human_full, macaque_filter)
mq_events = sv_list[[1]]; hu_events = sv_list[[2]];
# Read and filter data

cat("----------\nSubsetting macaque data...\n")
mqr = subsetSVs(mq_events)
mq_svs = mqr[[4]]
mq_svs = subset(mq_svs, Length <= maxlen)

cat("----------\nSubsetting human data...\n")
hur = subsetSVs(hu_events)
hu_svs = hur[[4]]
hu_svs = subset(hu_svs, Length <= maxlen)
# Subset data

#mq_svs$Length = -mq_svs$Length


hu_svs_noalu = subset(hu_svs, Length < 275 | Length > 325)
mq_svs_noalu = subset(mq_svs, Length < 275 | Length > 325)

hu_svs$Label = "Human"
mq_svs$Label = "Macaque"
hu_svs_noalu$Label = "Human no Alu"
mq_svs_noalu$Label = "Macaque no Alu"
# For Alu stuff

#cat("----------\nPre-binning macaque data...\n")
#mq_bins = binSVs(mq_svs)
#mq_bins$Species = "Macaque"
#mq_bins$count = -mq_bins$count

#cat("----------\nPre-binning human data...\n")
#hu_bins = binSVs(hu_svs)
#hu_bins$Species = "Human"

#hu_svs = subset(hu_svs, Length > 200)
sv_alleles = rbind(hu_svs, hu_svs_noalu, mq_svs, mq_svs_noalu)
#sv_bins = rbind(hu_bins, mq_bins)
#test_count = mq_events %>% group_by(Individual) %>% count(Individual)
# Count number of events per individual
# Combine the two datasets
######################

######################
# SV length histogram -- Fig 2A
cat("SV length distribution...\n")


fig3a = ggplot(sv_alleles, aes(x=Label, y=Length, fill=Label)) + 
  geom_quasirandom(size=2, alpha=0.7, width=0.25, color="#d3d3d3") +
  geom_boxplot(outlier.shape=NA, alpha=0.3, width=0.5, color="#000000") +
  labs(y="CNV Length") +
  bartheme() +
  theme(legend.position="none")

if(color_plots){
  fig3a = fig3a + scale_fill_manual(name="", values=c("Macaque"='#490092', "Macaque no Alu"='#490092',"Human"='#920000', "Human no Alu"="#920000"))
}else{
  #fig3a = fig3a + scale_fill_grey(name="", labels=c("Human","Macaque"))
  #fig3a = fig3a + scale_fill_manual(name="", values=c("Macaque"='#d6d6d6',"Human"='#5c5c5c'))
  fig3a = fig3a + scale_fill_manual(name="", values=c("Macaque"='#5c5c5c', "Macaque no Alu"="5c5c5c", "Human"=NA, "Human no Alu"=NA))
}

#print(fig3a)
# ggsave(filename="fig3a-alu-box.png", fig3a, width=6, height=4, units="in")


sv_alleles = rbind(hu_svs_noalu, mq_svs_noalu)
sv_alleles$Label[sv_alleles$Label=="Human no Alu"] = "Human"
sv_alleles$Label[sv_alleles$Label=="Macaque no Alu"] = "Macaque"
fig3a = ggplot(sv_alleles, aes(x=Label, y=Length, fill=Label)) + 
  geom_quasirandom(size=2, alpha=0.7, width=0.25, color="#d3d3d3") +
  geom_boxplot(outlier.shape=NA, alpha=0.3, width=0.5, color="#000000") +
  labs(y="CNV Length", x="") +
  bartheme() +
  theme(legend.position="none")

if(color_plots){
  fig3a = fig3a + scale_fill_manual(name="", values=c("Macaque"='#490092', "Human"='#920000'))
}else{
  #fig3a = fig3a + scale_fill_grey(name="", labels=c("Human","Macaque"))
  #fig3a = fig3a + scale_fill_manual(name="", values=c("Macaque"='#d6d6d6',"Human"='#5c5c5c'))
  fig3a = fig3a + scale_fill_manual(name="", values=c("Macaque"='#5c5c5c', "Human"=NA))
}

print(fig3a)
ggsave(filename="figS3.png", fig3a, width=6, height=4, units="in")

stop()









sv_alleles = subset(sv_alleles, Length > 200)
fig3a = ggplot(sv_alleles, aes(x=Label, y=Length, fill=Label)) + 
  geom_quasirandom(size=2, alpha=0.7, width=0.25, color="#d3d3d3") +
  geom_boxplot(outlier.shape=NA, alpha=0.3, width=0.5, color="#000000") +
  labs(y="CNV Length", x="") +
  bartheme() +
  theme(legend.position="none")

if(color_plots){
  fig3a = fig3a + scale_fill_manual(name="", values=c("Macaque"='#490092', "Human"='#920000'))
}else{
  #fig3a = fig3a + scale_fill_grey(name="", labels=c("Human","Macaque"))
  #fig3a = fig3a + scale_fill_manual(name="", values=c("Macaque"='#d6d6d6',"Human"='#5c5c5c'))
  fig3a = fig3a + scale_fill_manual(name="", values=c("Macaque"='#5c5c5c', "Human"=NA))
}

print(fig3a)
# ggsave(filename="fig3a-noalu-short-box.png", fig3a, width=6, height=4, units="in")

stop()
#if(savefiles){
#  outfile = "fig3a.pdf"
#  cat(" -> ", outfile, "\n")
#  ggsave(filename=outfile, fig3a, width=6, height=4, units="in")
#}
# SV length plot

cat(" -> SV length KS test...\n")
hu_svs = subset(hu_svs, Length > 200)
mq_svs = subset(mq_svs, Length > 200)
length_ks = wilcox.test(hu_svs$Length, mq_svs$Length)
print(length_ks)
# SV length KS test
# SV length histogram
######################
