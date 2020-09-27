############################################################
# For macaque paper, 08.19
# Makes gene proportion plots with SV and CAFE data
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
library("ggtree")

source("../lib/read_svs.r")
source("../lib/filter_svs.r")
source("../lib/subset_svs.r")
source("../lib/design.r")

cat("----------\n")

############################################################

savefiles = T
color_plots = F
read_data = T
rm_alus = T
# Run options

minlen = F
maxlen = 100000
# CNV length cutoffs

if(read_data){
  sv_list = readSVs()
  sv_list = filterSVs(sv_list, minlen, maxlen)
  mq_events = sv_list[[1]]; hu_events = sv_list[[2]];
  # Read and filter data
  
  if(rm_alus){
    hu_events = subset(hu_events, Length < 275 | Length > 325)
    mq_events = subset(mq_events, Length < 275 | Length > 325)
  }
  # For Alu stuff
  
  cat("----------\nSubsetting macaque data...\n")
  mqr = subsetSVs(mq_events)
  mq_svs = mqr[[4]]
  
  cat("----------\nSubsetting human data...\n")
  hur = subsetSVs(hu_events)
  hu_svs = hur[[4]]
  # Subset data
}

######################
cat("Reading CAFE data...\n")
cafe_data = read.csv("../cafe-data/cafe-traits.csv", header=TRUE)
cafe_data = cafe_data[order(cafe_data$node),]
#cafe_data = subset(cafe_data, Node.type=="Tip")
cafe_data$genes.changed = cafe_data$Genes.gained + cafe_data$Genes.lost
cafe_data$perc.genes.gained = cafe_data$Genes.gained / cafe_data$genes.changed
cafe_data$perc.genes.lost = cafe_data$Genes.lost / cafe_data$genes.changed

cafe_data$branch.str = paste(cafe_data$Genes.gained, "/", cafe_data$Genes.lost)

cat("Plotting tree\n")
tree = read.tree("../cafe-data/cafe-tree.tre")

fig4a = ggtree(tree, size=1, ladderize=F) +
  scale_color_manual(values=c("black","#db6d00")) +
  geom_tiplab(color="#333333", fontface='italic', size=5) +
  ggplot2::xlim(0, 135) +
  theme(plot.title=element_text(size=16, face="bold"))
  #geom_text(aes(x=branch, label=cafe_data$branch.str), size=3.5, 
  #          vjust=-.5, color="#333333")
print(fig4a)

######################
cat("Reading macaque CAFE genes\n")
mq_cafe = read.csv("../cafe-data/macaque-cafe-genes.csv")
mq_cafe_dels = sum(mq_cafe$Num.del)
mq_cafe_dups = sum(mq_cafe$Num.dup)

mq_cafe_total = mq_cafe_dels + mq_cafe_dups
mq_cafe_del_p = mq_cafe_dels / mq_cafe_total
mq_cafe_dup_p = mq_cafe_dups / mq_cafe_total

cat("Reading human CAFE genes\n")
hu_cafe = read.csv("../cafe-data/human-cafe-genes.csv")
hu_cafe_dels = sum(hu_cafe$Num.del)
hu_cafe_dups = sum(hu_cafe$Num.dup)

hu_cafe_total = hu_cafe_dels + hu_cafe_dups
hu_cafe_del_p = hu_cafe_dels / hu_cafe_total
hu_cafe_dup_p = hu_cafe_dups / hu_cafe_total

######################
# CAFE counts
mq_changes = cafe_data$genes.changed[cafe_data$Node.ID=="Mmula"]
mq_gains = cafe_data$Genes.gained[cafe_data$Node.ID=="Mmula"]
mq_losses = cafe_data$Genes.lost[cafe_data$Node.ID=="Mmula"]
mq_gains_p = cafe_data$perc.genes.gained[cafe_data$Node.ID=="Mmula"]
mq_losses_p = cafe_data$perc.genes.lost[cafe_data$Node.ID=="Mmula"]

hu_changes = cafe_data$genes.changed[cafe_data$Node.ID=="Hsapi"]
hu_gains = cafe_data$Genes.gained[cafe_data$Node.ID=="Hsapi"]
hu_losses = cafe_data$Genes.lost[cafe_data$Node.ID=="Hsapi"]
hu_gains_p = cafe_data$perc.genes.gained[cafe_data$Node.ID=="Hsapi"]
hu_losses_p = cafe_data$perc.genes.lost[cafe_data$Node.ID=="Hsapi"]

cafe_genes = data.frame("Label"=c("Macaque gene deletions/duplications","Macaque gene gains/losses",
                                  "Human gene deletions/duplications", "Human gene gains/losses"),
                        "total.svs"=c(mq_cafe_total, mq_changes, hu_cafe_total, hu_changes),
                        "num.del"=c(mq_cafe_dels, mq_losses, hu_cafe_dels, hu_losses),
                        "num.dup"=c(mq_cafe_dups, mq_gains, hu_cafe_dups, hu_gains),
                        "prop.del"=c(mq_cafe_del_p, mq_losses_p, hu_cafe_del_p, hu_losses_p),
                        "prop.dup"=c(mq_cafe_dup_p, mq_gains_p, hu_cafe_dup_p, hu_gains_p))

cafe_prop = subset(cafe_genes, select=c("Label","prop.del","prop.dup"))
cafe_prop_melt = melt(cafe_prop, id.vars="Label")

######################
# Macaque counts
mq_cafe_genes = data.frame("Label"=c("Polymorphic","Fixed",
                                     "Polymorphic","Fixed"),
                           "Type"=c("Deletions","Deletions","Duplications","Duplications"),
                           "prop"=c(mq_cafe_del_p,mq_losses_p,mq_cafe_dup_p,mq_gains_p),
                           "count"=c(mq_cafe_dels,mq_losses,mq_cafe_dups,mq_gains))

######################
# Human counts
hu_cafe_genes = data.frame("Label"=c("Polymorphic","Fixed",
                                     "Polymorphic","Fixed"),
                           "Type"=c("Deletions","Deletions","Duplications","Duplications"),
                           "prop"=c(hu_cafe_del_p,hu_losses_p,hu_cafe_dup_p,hu_gains_p),
                           "count"=c(hu_cafe_dels,hu_losses,hu_cafe_dups,hu_gains))

######################
# Gene proportions plot (Macaque)
fig4b = ggplot(mq_cafe_genes, aes(x=Label, y=prop, fill=Type, label=count)) +
  geom_bar(stat='identity') +
  geom_text(size=6, position=position_stack(vjust=0.5), color="#f2f2f2") +
  scale_x_discrete(limits=c("Polymorphic", "Fixed")) +
  scale_y_continuous(expand=c(0,0)) +
  labs(x="", y="Proportion of changes") +
  coord_flip() + 
  bartheme() +
  theme(legend.position="bottom")
if(color_plots){
  fig4b = fig4b + scale_fill_manual(name="", labels=c("Deletions", "Duplications"), values=c("#009292","#b585e5"))
}else{
  fig4b = fig4b + scale_fill_grey(name="", labels=c("Deletions","Duplications"))
}
print(fig4b)
# Gene proportions plot (Macaque)
######################
# Gene proportions plot (Human)
fig4c = ggplot(hu_cafe_genes, aes(x=Label, y=prop, fill=Type, label=count)) +
  geom_bar(stat='identity') +
  geom_text(size=4, position=position_stack(vjust=0.5), color="#f2f2f2") +
  scale_fill_manual(name="", labels=c("Deletions", "Duplications"), values=c("#009292","#b585e5")) +
  scale_x_discrete(limits=c("Polymorphic", "Fixed")) +
  scale_y_continuous(expand=c(0,0)) +
  ggtitle("Human") +
  labs(x="", y="Proportion of changes") +
  coord_flip() + 
  bartheme() +
  theme(legend.position="bottom")
print(fig4c)
# Gene proportions plot (Human)
######################

cat(" -> Chi-squared test for macaque SVs vs human SVs...\n")
sv_comp_counts = subset(cafe_genes, select=c("Label","num.del","num.dup"))
sv_comp_counts = subset(sv_comp_counts, !Label %in% c("Macaque gene gains/losses", "Human gene gains/losses"))
sv_comp_counts_t = t(sv_comp_counts[,2:ncol(sv_comp_counts)])
colnames(sv_comp_counts_t) <- sv_comp_counts[,1]
sv_comp_chi = chisq.test(sv_comp_counts_t)
print(sv_comp_chi)
# Macaque CAFE vs SV genes test

cat(" -> Chi-squared test for macaque SV vs CAFE genes...\n")
mq_comp_counts = subset(cafe_genes, select=c("Label","num.del","num.dup"))
mq_comp_counts = subset(mq_comp_counts, !Label %in% c("Human gene deletions/duplications", "Human gene gains/losses"))
mq_comp_counts_t = t(mq_comp_counts[,2:ncol(mq_comp_counts)])
colnames(mq_comp_counts_t) <- mq_comp_counts[,1]
mq_comp_chi = chisq.test(mq_comp_counts_t)
print(mq_comp_chi)
# Macaque CAFE vs SV genes test

cat(" -> Chi-squared test for human SV vs CAFE genes...\n")
hu_comp_counts = subset(cafe_genes, select=c("Label","num.del","num.dup"))
hu_comp_counts = subset(hu_comp_counts, !Label %in% c("Macaque gene deletions/duplications", "Macaque gene gains/losses"))
hu_comp_counts_t = t(hu_comp_counts[,2:ncol(hu_comp_counts)])
colnames(hu_comp_counts_t) <- hu_comp_counts[,1]
hu_comp_chi = chisq.test(hu_comp_counts_t)
print(hu_comp_chi)
# Human CAFE vs SV genes test
# Gene proportions plot
######################


######################
# Combine plots for figure
cat("Combining proportion plots...\n")
fig4 = plot_grid(fig4a, fig4b, nrow=2, labels=c("A","B"), label_size=24)
# For macaque proportion plot only

print(fig4)

if(savefiles){
  if(color_plots){
    outfile = "fig4-temp.pdf"
  }else{
    outfile = "fig4-grey-temp.pdf"
  }
  cat(" -> ", outfile, "\n")
  ggsave(filename=outfile, fig4, width=6, height=8, units="in")
}
# Save the figure
######################