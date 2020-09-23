############################################################
# For macaque paper, 08.19
# Compares chromosome lengths and # of CNVs
# Gregg Thomas
############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(ggplot2)
library(reshape2)
library(plyr)
library(grid)
library(ggpubr)
library(cowplot)

source("../lib/read_svs.r")
source("../lib/filter_svs.r")
source("../lib/subset_svs.r")

############################################################

savefiles = F
in_data = read.csv("../data/macaque-cnv-chrome-counts.csv")

######################

######################
# All CNVs and chromosome length
figS1a = ggplot(in_data, aes(Length, Num.CNVs)) +
  geom_point(size=3, color="#666666") +
  geom_smooth(method="lm", fullrange=T, size=0.75, linetype="dashed", alpha=0, color="#333333") +
  labs(x="Chromosome length", y="# CNVs") +
  theme_classic() +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=16), 
        axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0),color="black"), 
        axis.title.x=element_text(margin=margin(t=0,r=0,b=0,l=0),color="black"),
        axis.line=element_line(colour='#595959',size=0.75),
        axis.ticks=element_line(colour="#595959",size = 1),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position="right",
        legend.key.width = unit(0.75,  unit = "cm"),
        legend.spacing.x = unit(0.25, 'cm'),
        legend.title = element_blank(),
        legend.text=element_text(size=12),
        plot.title = element_text(hjust=0.5, size=16),
        plot.margin = unit(c(1,1,0,1), "cm")
  )
print(figS1a)
######################

######################
# Deletions and chromosome length
figS1b = ggplot(in_data, aes(Length, Num.dels)) +
  geom_point(size=3, color="#666666") +
  geom_smooth(method="lm", fullrange=T, size=0.75, linetype="dashed", alpha=0, color="#333333") +
  labs(x="Chromosome length", y="# Deletions") +
  theme_classic() +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=16), 
        axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0),color="black"), 
        axis.title.x=element_text(margin=margin(t=0,r=0,b=0,l=0),color="black"),
        axis.line=element_line(colour='#595959',size=0.75),
        axis.ticks=element_line(colour="#595959",size = 1),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position="right",
        legend.key.width = unit(0.75,  unit = "cm"),
        legend.spacing.x = unit(0.25, 'cm'),
        legend.title = element_blank(),
        legend.text=element_text(size=12),
        plot.title = element_text(hjust=0.5, size=16),
        plot.margin = unit(c(1,1,0,1), "cm")
  )
print(figS1b)
######################

######################
# Duplications and chromosome length
figS1c = ggplot(in_data, aes(Length, Num.dups)) +
  geom_point(size=3, color="#666666") +
  geom_smooth(method="lm", fullrange=T, size=0.75, linetype="dashed", alpha=0, color="#333333") +
  #ggtitle("Gene family changes vs SV duplications") +
  #scale_color_manual(name="", values=c('#490092','#920000'), labels=c("Macaque","Human"), drop=FALSE) +
  #scale_y_continuous(limits=c(0, 5)) +
  labs(x="Chromosome length", y="# Duplications") +
  theme_classic() +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=16), 
        axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0),color="black"), 
        axis.title.x=element_text(margin=margin(t=0,r=0,b=0,l=0),color="black"),
        axis.line=element_line(colour='#595959',size=0.75),
        axis.ticks=element_line(colour="#595959",size = 1),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position="right",
        legend.key.width = unit(0.75,  unit = "cm"),
        legend.spacing.x = unit(0.25, 'cm'),
        legend.title = element_blank(),
        legend.text=element_text(size=12),
        plot.title = element_text(hjust=0.5, size=16),
        plot.margin = unit(c(1,1,0,1), "cm")
  )
print(figS1c)
######################

######################
# Making the figure

p = plot_grid(figS1a, figS1b, figS1c, nrow=1, labels=c("A","B","C"), label_size=24)

print(p)

if(savefiles){
  outfile = "figS1.pdf"
  cat(" -> ", outfile, "\n")
  ggsave(filename=outfile, p, width=14, height=4, units="in")
}

######################