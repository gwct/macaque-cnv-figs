############################################################
# For macaque paper, 08.19
# Makes length histograms for human and macaque SVs
# Gregg Thomas
############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(ggplot2)
library(reshape2)
library(cowplot)
library(ggridges)

source("../lib/read_svs.r")
source("../lib/filter_svs.r")
source("../lib/subset_svs.r")

cat("----------\n")

############################################################

savefiles = F
color_plots = T
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
# Supplemental figure options

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

#sv_alleles = rbind(hu_svs, mq_svs)
sv_alleles = unique(subset(mq_svs, select=c(SV.key,SV.freq,Species)))
allele_freq_counts = sv_alleles %>% group_by(SV.freq) %>% count(SV.freq)
#allele_freqs$SV.freq = as.character(allele_freqs$SV.freq)

#test_count = mq_events %>% group_by(Individual) %>% count(Individual)
# Count number of events per individual


# Combine the two datasets
######################

######################
# SV length histogram -- Fig 2A
cat("SV allele frequency spectrum...\n")
fig3a = ggplot(sv_alleles, aes(x=SV.freq, group=Species)) + 
  geom_histogram(position="identity", alpha=0.8, bins=30) +
  #geom_density(alpha=0.3) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="CNV length", y="Count") +
  theme_classic() +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=16), 
        axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0),color="black"), 
        axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0),color="black"),
        axis.line=element_line(colour='#595959',size=0.75),
        axis.ticks=element_line(colour="#595959",size = 1),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position="right",
        legend.key.width = unit(0.75,  unit = "cm"),
        legend.spacing.x = unit(0.25, 'cm'),
        legend.title = element_blank(),
        legend.text=element_text(size=12),
        plot.title = element_text(hjust=0.5, size=20),
        plot.margin = unit(c(1,1,1,1), "cm")
  )

if(color_plots){
  fig3a = fig3a + scale_fill_manual(name="", values=c("Macaque"='#490092',"Human"='#920000'))
}else{
  fig3a = fig3a + scale_fill_grey(name="", labels=c("Human","Macaque"))
}

print(fig3a)

