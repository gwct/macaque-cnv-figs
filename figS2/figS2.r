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
source("../lib/design.r")

############################################################

savefiles = T
rm_alus = T
in_data = read.csv("../cnv-calls/macaque-cnv-chrome-counts.csv")

if(rm_alus){
  in_data = subset(in_data, Length < 275 | Length > 325)
}
# For Alu stuff

######################

######################
# All CNVs and chromosome length
figS2a = ggplot(in_data, aes(Length, Num.CNVs)) +
  geom_point(size=3, color="#666666") +
  geom_smooth(method="lm", fullrange=T, size=0.75, linetype="dashed", alpha=0, color="#333333") +
  labs(x="Chromosome length", y="# CNVs") +
  bartheme()
print(figS2a)
######################

######################
# Deletions and chromosome length
figS2b = ggplot(in_data, aes(Length, Num.dels)) +
  geom_point(size=3, color="#666666") +
  geom_smooth(method="lm", fullrange=T, size=0.75, linetype="dashed", alpha=0, color="#333333") +
  labs(x="Chromosome length", y="# Deletions") +
  bartheme()
print(figS2b)
######################

######################
# Duplications and chromosome length
figS2c = ggplot(in_data, aes(Length, Num.dups)) +
  geom_point(size=3, color="#666666") +
  geom_smooth(method="lm", fullrange=T, size=0.75, linetype="dashed", alpha=0, color="#333333") +
  #ggtitle("Gene family changes vs SV duplications") +
  #scale_color_manual(name="", values=c('#490092','#920000'), labels=c("Macaque","Human"), drop=FALSE) +
  #scale_y_continuous(limits=c(0, 5)) +
  labs(x="Chromosome length", y="# Duplications") +
  bartheme()
print(figS2c)
######################

######################
# Making the figure

p = plot_grid(figS2a, figS2b, figS2c, nrow=1, labels=c("A","B","C"), label_size=24)

print(p)

if(savefiles){
  outfile = "figS2.pdf"
  cat(" -> ", outfile, "\n")
  ggsave(filename=outfile, p, width=14, height=4, units="in")
}

######################