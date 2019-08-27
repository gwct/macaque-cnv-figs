############################################################
# For macaque paper, 08.19
# Compares gene family gains/losses to SVs
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

savefiles = T

in_data = read.csv("macaque-cafe-vs-sv-genes-filtered.csv")
in_data$net.sv = in_data$SV.dups + in_data$SV.dels

######################
# Duplications
figX2a = ggplot(in_data, aes(CAFE.change, SV.dups)) +
  geom_point(size=3, alpha=0.25) +
  geom_smooth(method="lm", fullrange=T, size=0.75, linetype="dashed", alpha=0, color="#333333") +
  labs(x="Gene family change", y="# SV duplications") +
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
print(figX2a)
######################

######################
# Deletions
in_data$SV.dels = abs(in_data$SV.dels)
figX2b = ggplot(in_data, aes(CAFE.change, SV.dels)) +
  geom_point(size=3, alpha=0.25) +
  geom_smooth(method="lm", fullrange=T, size=0.75, linetype="dashed", alpha=0, color="#333333") +
  labs(x="Gene family change", y="# SV deletions") +
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
print(figX2b)
######################

######################
# Net SVs
figX2c = ggplot(in_data, aes(CAFE.change, net.sv)) +
  geom_point(size=3, alpha=0.25) +
  geom_smooth(method="lm", fullrange=T, size=0.75, linetype="dashed", alpha=0, color="#333333") +
  labs(x="Gene family change", y="Net SV change") +
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
print(figX2c)
######################

######################
# Making figure

p = plot_grid(figX2a, figX2b, figX2c, nrow=1, labels=c("A","B","C"), label_size=24)

print(p)

if(savefiles){
  outfile = "figX2.pdf"
  cat(" -> ", outfile, "\n")
  ggsave(filename=outfile, p, width=14, height=4, units="in")
}

######################