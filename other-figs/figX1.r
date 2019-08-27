############################################################
# For macaque paper, 08.19
# Makes scatterplots of sv length vs gene overlap.
# Gregg Thomas
############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(ggplot2)
library(reshape2)
library(plyr)
library(ggpubr)
library(cowplot)
library("ggtree")

source("../lib/read_svs.r")
source("../lib/filter_svs.r")
source("../lib/subset_svs.r")

cat("----------\n")

############################################################

savefiles = T
# Run options

minlen = F
maxlen = 100000
# CNV length cutoffs

sv_list = readSVs()
sv_list = filterSVs(sv_list, minlen, maxlen)
mq_events = sv_list[[1]]; hu_events = sv_list[[2]];
# Read and filter data

cat("----------\nSubsetting macaque data...\n")
mqr = subsetSVs(mq_events)
mq_svs = mqr[[4]]

cat("----------\nSubsetting human data...\n")
hur = subsetSVs(hu_events)
hu_svs = hur[[4]]
# Subset data

######################
# Macaque
mq_reg = lm(mq_svs$Num.genes ~ mq_svs$Length)

mq_svs$Species = factor(mq_svs$Species, levels=c("Macaque", "Human"))

figX1a = ggplot(mq_svs, aes(Length, Num.genes, color=Species)) +
  geom_point(size=3, alpha=0.2) +
  ggtitle("") +
  labs(x="CNV length (bp)", y="# Genes overlapping CNV") +
  scale_color_manual(name="", values=c("Macaque"='#490092',"Human"='#920000'), drop=FALSE) +
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
        plot.margin = unit(c(1,1,0,1), "cm")
  ) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))
print(figX1a)

######################
# Human
hu_reg = lm(hu_svs$Num.genes ~ hu_svs$Length)

figX1b = ggplot(hu_svs, aes(Length, Num.genes, color=Species)) +
  geom_point(size=3, alpha=0.2) +
  ggtitle("") +
  labs(x="CNV length (bp)", y="") +
  scale_color_manual(name="", values=c("Macaque"='#490092',"Human"='#920000'), drop=FALSE) +
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
        plot.margin = unit(c(1,1,0,0), "cm")
  )
print(figX1b)

######################
# Combine plots for figure

cat("Combining proportion plots...\n")
prow = plot_grid(figX1a + theme(legend.position="none"),
                 figX1b + theme(legend.position="none"),
                 align = 'vh',
                 labels = c("A", "B"),
                 label_size = 24,
                 hjust = -1,
                 nrow = 1)
# Combine panels A and B

legend_b = get_legend(figX1a + theme(legend.direction="horizontal", legend.justification="center", legend.box.just="bottom"))
# Extract the legend from one of the plots

p = plot_grid(prow, legend_b, ncol=1, rel_heights=c(1, 0.1))
# Add the legend underneath the row we made earlier with 10% of the height of the row.
print(p)

if(savefiles){
  outfile = "figX1.pdf"
  cat(" -> ", outfile, "\n")
  ggsave(filename=outfile, p, width=10, height=5, units="in")
}
# Save the figure
######################