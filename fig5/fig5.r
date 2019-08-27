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

cat("----------\n")

############################################################

savefiles = T
color_plots = T
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
cat("Reading CAFE data...\n")
cafe_data = read.csv("../data/cafe-traits.csv", header=TRUE)
#cafe_data = subset(cafe_data, Node.type=="Tip")
cafe_data$genes.changed = cafe_data$Genes.gained + cafe_data$Genes.lost
cafe_data$perc.genes.gained = cafe_data$Genes.gained / cafe_data$genes.changed
cafe_data$perc.genes.lost = cafe_data$Genes.lost / cafe_data$genes.changed

cat("Plotting tree\n")
tree = read.tree("../data/cafe-tree.tre")

fig5a = ggtree(tree, size=1, ladderize=F) +
  scale_color_manual(values=c("black","#db6d00")) +
  geom_tiplab(color="#333333", fontface='italic', size=5) +
  ggplot2::xlim(0, 115) +
  theme(plot.title=element_text(size=16, face="bold"))
print(fig5a)

######################
cat("Reading macaque CAFE genes\n")
mq_cafe = read.csv("../data/macaque-cafe-genes-filtered.csv")
mq_cafe_genes = subset(mq_cafe, SV.key %in% mq_svs$SV.key)
mq_cafe_genes = subset(mq_cafe_genes, SV.type=="<DEL>" | SV.type=="<DUP>")

mq_cafe_dels = length(mq_cafe_genes$SV.type[mq_cafe_genes$SV.type=="<DEL>"])
mq_cafe_dups = length(mq_cafe_genes$SV.type[mq_cafe_genes$SV.type=="<DUP>"])
mq_cafe_total = mq_cafe_dels + mq_cafe_dups
mq_cafe_del_p = mq_cafe_dels / mq_cafe_total
mq_cafe_dup_p = mq_cafe_dups / mq_cafe_total

cat("Reading human CAFE genes\n")
hu_cafe = read.csv("../data/brandler-cafe-genes.csv")
hu_cafe_genes = subset(hu_cafe, SV.key %in% hu_svs$SV.key)
hu_cafe_genes = subset(hu_cafe_genes, SV.type=="<DEL>" | SV.type=="<DUP>")

hu_cafe_dels = length(hu_cafe_genes$SV.type[hu_cafe_genes$SV.type=="<DEL>"])
hu_cafe_dups = length(hu_cafe_genes$SV.type[hu_cafe_genes$SV.type=="<DUP>"])
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
fig5b = ggplot(mq_cafe_genes, aes(x=Label, y=prop, fill=Type, label=count)) +
  geom_bar(stat='identity') +
  geom_text(size=6, position=position_stack(vjust=0.5), color="#f2f2f2") +
  scale_x_discrete(limits=c("Polymorphic", "Fixed")) +
  scale_y_continuous(expand=c(0,0)) +
  labs(x="", y="Proportion of changes") +
  coord_flip() + 
  theme_classic() +
  theme(axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=16),
        axis.title=element_text(size=12), 
        axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0),color="black"), 
        axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0),color="black"),
        axis.line=element_line(colour='#595959',size=0.75),
        axis.ticks=element_line(colour="#595959",size = 1),
        axis.ticks.length=unit(0.2,"cm"),
        axis.ticks.y=element_blank(),
        legend.position="bottom",
        plot.margin = unit(c(1,1,0,0), "cm"),
        plot.title = element_text(size=16, hjust=0.5, color="#333333")
  )

if(color_plots){
  fig5b = fig5b + scale_fill_manual(name="", labels=c("Deletions", "Duplications"), values=c("#009292","#b585e5"))
}else{
  fig5b = fig5b + scale_fill_grey(name="", labels=c("Deletions","Duplications"))
}

print(fig5b)
# Gene proportions plot (Macaque)
######################
# Gene proportions plot (Human)
fig5c = ggplot(hu_cafe_genes, aes(x=Label, y=prop, fill=Type, label=count)) +
  geom_bar(stat='identity') +
  geom_text(size=4, position=position_stack(vjust=0.5), color="#f2f2f2") +
  scale_fill_manual(name="", labels=c("Deletions", "Duplications"), values=c("#009292","#b585e5")) +
  scale_x_discrete(limits=c("Polymorphic", "Fixed")) +
  scale_y_continuous(expand=c(0,0)) +
  ggtitle("Human") +
  labs(x="", y="Proportion of changes") +
  coord_flip() + 
  theme_classic() +
  theme(axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=12), 
        axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0),color="black"), 
        axis.title.x=element_text(margin=margin(t=0,r=0,b=0,l=0),color="black"),
        axis.line=element_line(colour='#595959',size=0.75),
        axis.ticks=element_line(colour="#595959",size = 1),
        axis.ticks.length=unit(0.2,"cm"),
        axis.ticks.y=element_blank(),
        legend.position="bottom",
        plot.margin = unit(c(1,1,0,0), "cm"),
        plot.title = element_text(size=16, hjust=0.5, color="#333333")
  )
print(fig5c)
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

#prow = plot_grid(fig5b + theme(legend.position="none"),
#                 fig5c + theme(legend.position="none"),
#                 align = 'vh',
#                 labels = c("B", "C"),
#                 label_size = 24,
#                 hjust = -1,
#                 nrow = 1)
#legend_b = get_legend(fig5b + theme(legend.direction="horizontal", legend.justification="center", legend.box.just="bottom"))
#ptop = plot_grid(nullGrob(), fig5a, nullGrob(), nrow=1, rel_widths=c(0.1,0.6,0.1))
#
#pcombo = plot_grid(ptop, prow, nrow=2, labels=c("A",""), label_size=24, rel_heights=c(1,0.8))
#fig5 = plot_grid(pcombo, legend_b, ncol=1, rel_heights=c(1, 0.1))
# For macaque and human proportion plots

fig5 = plot_grid(fig5a, fig5b, nrow=2, labels=c("A","B"), label_size=24)
# For macaque proportion plot only

print(fig5)

if(savefiles){
  if(color_plots){
    outfile = "fig5.pdf"
  }else{
    outfile = "fig5-grey.pdf"
  }
  cat(" -> ", outfile, "\n")
  ggsave(filename=outfile, fig5, width=6, height=8, units="in")
}
# Save the figure
######################