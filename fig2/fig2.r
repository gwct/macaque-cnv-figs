############################################################
# For macaque paper, 08.19
# Maps CNVs to chromes
# Gregg Thomas
############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(ggplot2)
library(reshape2)
library(plyr)
#library(gridExtra)
library(grid)
library(ggpubr)
library(cowplot)

source("../lib/read_svs.r")
source("../lib/filter_svs.r")
source("../lib/subset_svs.r")

############################################################

savefiles = T
color_plots = T
readcnvs = F
# Run options

minlen = F
maxlen = 100000
# CNV length cutoffs

in_data = read.csv("../data/macaque-cnv-chrome-counts.csv")
in_data$Chromosome = as.character(in_data$Chromosome)

if(readcnvs){
  sv_list = readSVs()
  sv_list = filterSVs(sv_list, minlen, maxlen)
  mq_events = sv_list[[1]]; hu_events = sv_list[[2]];
  # Read and filter data
  
  cat("----------\nSubsetting macaque data...\n")
  mqr = subsetSVs(mq_events)
  mq_svs = mqr[[4]]
}

######################

plots = list()

chr_p = ggplot(data.frame()) +
  scale_y_continuous(limits=c(1,23), breaks=1:22, labels=rev(in_data$label)) +
  scale_x_continuous(limits=c(0,max(in_data$Length-1)), breaks=NULL, expand=c(0,0)) +
  ylab("") +
  xlab("") +
  theme_classic() +
  theme(axis.text=element_text(size=16), 
        axis.title=element_blank(), 
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        legend.position="bottom",
        legend.key.width = unit(0.75,  unit = "cm"),
        legend.spacing.x = unit(0.25, 'cm'),
        legend.title = element_blank(),
        legend.text=element_text(size=12),
        plot.title = element_text(hjust=0.5, size=16),
        plot.margin = unit(c(1,1,0,1), "cm")
  ) +
  guides(color = guide_legend(override.aes = list(size=8)))

if(color_plots){
  chr_p = chr_p + scale_color_manual(name="", labels=c("Deletions","Duplications"), values=c("#006ddb","#db6d00"))
}else{
  #chr_p = chr_p + scale_color_grey(name="", labels=c("Deletions","Duplications"))
  chr_p = chr_p + scale_color_manual(name="", labels=c("Deletions","Duplications"), values=c("#000000","#000000")) +
    theme(legend.position="none")
}

blah = data.frame("blah"=c(1))
# Dummy data frame so the chrome segment doesn't layer
# Initial plot layer

ypt = 22

for(i in 1:nrow(in_data)) {
  row <- in_data[i,]
  cat(ypt, row$Chromosome, "\n")
  cur_cnvs = subset(mq_svs, Chromosome == row$Chromosome)
  cur_cnvs$ypt = ypt
  row$ypt = ypt
  
  cur_cnvs[cur_cnvs$Type=="<DEL>", ]$ypt = cur_cnvs[cur_cnvs$Type=="<DEL>", ]$ypt + 0.13
  cur_cnvs[cur_cnvs$Type=="<DUP>", ]$ypt = cur_cnvs[cur_cnvs$Type=="<DUP>", ]$ypt - 0.13
  # Adjusts the points so deletions appear above the line and duplications below the line
  
  chr_p = chr_p + geom_point(data=cur_cnvs, aes(Pos,ypt,color=Type), shape="|", size=3) +
    geom_segment(data=blah, x=0,xend=row$Length-1,y=ypt,yend=ypt, color="#666666") +
    geom_segment(data=blah, x=0,xend=row$Length-1,y=ypt,yend=ypt, color="#b3b3b3", size=7.2, alpha=0.2) #+
    
  ypt = ypt - 1
}
# Add the chromosome segments one at a time

######################

print(chr_p)
if(savefiles){
  if(color_plots){
    outfile = "fig2.pdf"
  }else{
    outfile = "fig2-grey.pdf"
  }
  cat(" -> ", outfile, "\n")
  ggsave(filename=outfile, chr_p, width=8.5, height=11, units="in")
}
# Save the figure
######################
