############################################################
# For macaque paper, 08.19
# Makes proportion charts for sv types and sv bases
# Gregg Thomas
############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(ggplot2)
library(reshape2)
library(cowplot)

source("../lib/read_svs.r")
source("../lib/filter_svs.r")
source("../lib/subset_svs.r")
source("../lib/design.r")

cat("----------\n")

############################################################

savefiles = F
color_plots = F
rm_alus = T
# Run options

maxlen = 100000
minlen = F
# CNV length cutoffs

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

######################
# SV type proportion plot -- Fig 1B
cat("Geting SV type counts...\n")

mq_dels = length(mq_svs$Type[mq_svs$Type=="<DEL>"])
mq_dups = length(mq_svs$Type[mq_svs$Type=="<DUP>"])
#mq_invs = length(mq_svs$Type[mq_svs$Type=="<INV>"])
mq_total = mq_dels + mq_dups
mq_dels_p = mq_dels / mq_total
mq_dups_p = mq_dups / mq_total

hu_dels = length(hu_svs$Type[hu_svs$Type=="<DEL>"])
hu_dups = length(hu_svs$Type[hu_svs$Type=="<DUP>"])
#hu_invs = length(hu_svs$Type[hu_svs$Type=="<INV>"])
hu_total = hu_dels + hu_dups
hu_dels_p = hu_dels / hu_total
hu_dups_p = hu_dups / hu_total
# Getting alleles for both species less than 100000 bp long and counting SV types. Also removing inversions for humans.

sv_types = data.frame("Species"=c("Macaque","Human"),
                      "total.svs"=c(mq_total, hu_total),
                      "num.del"=c(mq_dels, hu_dels),
                      "num.dup"=c(mq_dups, hu_dups),
                      "prop.del"=c(mq_dels_p, hu_dels_p),
                      "prop.dup"=c(mq_dups_p, hu_dups_p))

sv_types_prop = subset(sv_types, select=c("Species","prop.del","prop.dup"))
sv_types_prop_melt = melt(sv_types_prop, id.vars="Species")
# Organizing type counts

cat("Plotting SV type proportions...\n")
fig1b = ggplot(sv_types_prop_melt, aes(x=Species, y=value, fill=variable)) +
  geom_bar(stat='identity') +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="", y="Proportion of CNV types") +
  coord_flip() + 
  bartheme()
if(color_plots){
  fig1b = fig1b + scale_fill_manual(name="", labels=c("Deletions","Duplications"), values=c("#006ddb","#db6d00"))
}else{
  fig1b = fig1b + scale_fill_grey(name="", labels=c("Deletions","Duplications"))
}
print(fig1b)

#if(savefiles){
#  outfile = "fig1b.pdf"
#  cat(" -> ", outfile, "\n")
#  ggsave(fig1b, filename=outfile, width=6, height=4, units="in")
#}
# SV type proportion plot

cat(" -> Chi-squared test for CNV types...\n")
sv_types_counts = subset(sv_types, select=c("Species","num.del","num.dup"))
sv_types_counts_t = t(sv_types_counts[,2:ncol(sv_types_counts)])
colnames(sv_types_counts_t) <- sv_types_counts[,1]
type_chi = chisq.test(sv_types_counts_t)
print(sv_types_counts_t)
print(type_chi)
# SV type chi-squared test.
# SV type proportion plot -- Fig 1B
######################


######################
# SV bases plot -- Fig 1C
cat("Getting SV base counts...\n")
mq_del_bases = sum(mq_svs$Length[mq_svs$Type=="<DEL>"])
mq_dup_bases = sum(mq_svs$Length[mq_svs$Type=="<DUP>"])
#mq_invs = sum(mq_svs$Type[mq_svs$Type=="<INV>"])
mq_total_bases = mq_del_bases + mq_dup_bases
mq_del_bases_p = mq_del_bases / mq_total_bases
mq_dup_bases_p = mq_dup_bases / mq_total_bases

hu_del_bases = sum(hu_svs$Length[hu_svs$Type=="<DEL>"])
hu_dup_bases = sum(hu_svs$Length[hu_svs$Type=="<DUP>"])
#hu_invs = sum(hu_svs$Type[hu_svs$Type=="<INV>"])
hu_total_bases = hu_del_bases + hu_dup_bases
hu_del_bases_p = hu_del_bases / hu_total_bases
hu_dup_bases_p = hu_dup_bases / hu_total_bases
# Counting base types.

sv_bases = data.frame("Species"=c("Macaque","Human"),
                      "total.svs"=c(mq_total_bases, hu_total_bases),
                      "num.del"=c(mq_del_bases, hu_del_bases),
                      "num.dup"=c(mq_dup_bases, hu_dup_bases),
                      "prop.del"=c(mq_del_bases_p, hu_del_bases_p),
                      "prop.dup"=c(mq_dup_bases_p, hu_dup_bases_p))

sv_bases_prop = subset(sv_bases, select=c("Species","prop.del","prop.dup"))
sv_bases_prop_melt = melt(sv_bases_prop, id.vars="Species")
# Organizing base counts

cat("Plotting SV base proportions...\n")
fig1c = ggplot(sv_bases_prop_melt, aes(x=Species, y=value, fill=variable)) +
  geom_bar(stat='identity') +
  scale_y_continuous(expand = c(0,0)) +
  labs(x="", y="Proportion of bases affected by CNVs") +
  coord_flip() + 
  bartheme()
if(color_plots){
  fig1c = fig1c + scale_fill_manual(name="", labels=c("Deletions","Duplications"), values=c("#006ddb","#db6d00"))
}else{
  fig1c = fig1c + scale_fill_grey(name="", labels=c("Deletions","Duplications"))
}
print(fig1c)

#if(savefiles){
#  outfile = "fig1c.pdf"
#  cat(" -> ", outfile, "\n")
#  ggsave(fig1c, filename=outfile, width=6, height=4, units="in")
#}

cat(" -> Chi-squared test for CNV bases...\n")
sv_bases_counts = subset(sv_bases, select=c("Species","num.del","num.dup"))
sv_bases_counts_t = t(sv_bases_counts[,2:ncol(sv_bases_counts)])
colnames(sv_bases_counts_t) <- sv_bases_counts[,1]
base_chi = chisq.test(sv_bases_counts_t)
print(sv_bases_counts_t)
print(base_chi)
# SV bases chi-squared test.
# SV bases plot
######################


######################
# Combine plots for figure

cat("Combining proportion plots...\n")
prow = plot_grid(fig1b + theme(legend.position="none"),
                 fig1c + theme(legend.position="none"),
                 align = 'vh',
                 labels = c("B", "C"),
                 label_size = 24,
                 hjust = -1,
                 nrow = 1)
# Combine panels B and C


legend_b = get_legend(fig1b + theme(legend.direction="horizontal", legend.justification="center", legend.box.just="bottom"))
# Extract the legend from one of the plots

p = plot_grid(prow, legend_b, ncol=1, rel_heights=c(1, 0.2))
# Add the legend underneath the row we made earlier with 10% of the height of the row.
print(p)

if(savefiles){
  if(color_plots){
    outfile = "fig1bc.pdf"
  }else{
    outfile = "fig1bc-grey.pdf"
  }
  cat(" -> ", outfile, "\n")
  ggsave(filename=outfile, p, width=10, height=4, units="in")
}
# Save the figure.
######################








