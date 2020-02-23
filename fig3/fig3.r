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

sv_alleles = rbind(hu_svs, mq_svs)
#test_count = mq_events %>% group_by(Individual) %>% count(Individual)
# Count number of events per individual
# Combine the two datasets

sv_alleles = subset(sv_alleles, Length < 275 | Length > 325)
# For Alu stuff
######################

######################
# SV length histogram -- Fig 2A
cat("SV length distribution...\n")

fig3a = ggplot(sv_alleles, aes(x=Length, fill=Species, color=Species)) + 
  geom_histogram(alpha=0.8, position="identity", bins=50) +
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
  #fig3a = fig3a + scale_fill_grey(name="", labels=c("Human","Macaque"))
  #fig3a = fig3a + scale_fill_manual(name="", values=c("Macaque"='#d6d6d6',"Human"='#5c5c5c'))
  fig3a = fig3a + scale_fill_manual(name="", values=c("Macaque"='#5c5c5c',"Human"=NA))
  }
fig3a = fig3a + scale_color_manual(name="", values=c("Macaque"=NA,"Human"='#000000'))

print(fig3a)

#if(savefiles){
#  outfile = "fig3a.pdf"
#  cat(" -> ", outfile, "\n")
#  ggsave(filename=outfile, fig3a, width=6, height=4, units="in")
#}
# SV length plot

cat(" -> SV length KS test...\n")
length_ks = ks.test(hu_svs$Length, mq_svs$Length)
print(length_ks)
# SV length KS test
# SV length histogram
######################


######################
# SV deletion length histogram
cat("SV deletion length distribution...\n")
mq_dels = subset(mq_svs, Type=="<DEL>")
hu_dels = subset(hu_svs, Type=="<DEL>")
sv_dels = rbind(hu_dels, mq_dels)

sv_dels = subset(sv_dels, Length < 275 | Length > 325)
# For Alu stuff

fig3b = ggplot(sv_dels, aes(x=Length, fill=Species, color=Species)) + 
  geom_histogram(alpha=0.8, position="identity", bins=50) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="Deletion length", y="Count") +
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
  )

if(color_plots){
  fig3b = fig3b + scale_fill_manual(name="", values=c("Macaque"='#490092',"Human"='#920000'))
}else{
  #fig3b = fig3b + scale_fill_grey(name="", labels=c("Human","Macaque"))
  fig3b = fig3b + scale_fill_manual(name="", values=c("Macaque"='#5c5c5c',"Human"=NA))
}
fig3b = fig3b + scale_color_manual(name="", values=c("Macaque"=NA,"Human"='#000000'))

print(fig3b)

#if(savefiles){
#  outfile = "fig3b.pdf"
#  ggsave(filename=outfile, fig3b, width=6, height=4, units="in")
#}
# SV del length plot

cat(" -> SV deletion length KS test...\n")
del_ks = ks.test(mq_dels$Length, hu_dels$Length)
print(del_ks)
# SV del length KS test
# SV deletion length histogram
######################


######################
# SV duplication length histogram
cat("SV deletion length distribution...\n")
mq_dups = subset(mq_svs, Type=="<DUP>")
hu_dups = subset(hu_svs, Type=="<DUP>")
sv_dups = rbind(hu_dups, mq_dups)

sv_dups = subset(sv_dups, Length < 275 | Length > 325)
# For Alu stuff

fig3c = ggplot(sv_dups, aes(x=Length, fill=Species, color=Species)) + 
  geom_histogram(alpha=0.8, position="identity", bins=50) +
  #scale_fill_manual(name="", labels=c("Macaque","Human"), values=c("#006ddb","#db6d00")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="Duplication length", y="Count") +
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
  )

if(color_plots){
  fig3c = fig3c + scale_fill_manual(name="", values=c("Macaque"='#490092',"Human"='#920000'))
}else{
  #fig3c = fig3c + scale_fill_grey(name="", labels=c("Human","Macaque"))
  fig3c = fig3c + scale_fill_manual(name="", values=c("Macaque"='#5c5c5c',"Human"=NA))
}
fig3c = fig3c + scale_color_manual(name="", values=c("Macaque"=NA,"Human"='#000000'))

print(fig3c)

#if(savefiles){
#  outfile = "fig3c.pdf"
#  ggsave(filename=outfile, fig3c, width=6, height=4, units="in")
#}
# SV del length plot

cat(" -> SV duplication length KS test...\n")
dup_ks = ks.test(mq_dups$Length, hu_dups$Length)
print(dup_ks)
# SV del length KS test
# SV deletion length histogram
######################


######################
# Combine plots for figure

cat("Combining proportion plots...\n")
prow = plot_grid(fig3b + theme(legend.position="none"),
                 fig3c + theme(legend.position="none"),
                 align = 'vh',
                 labels = c("B", "C"),
                 label_size = 24,
                 hjust = -1,
                 nrow = 1)
# Combine panels B and C

pcombo = plot_grid(fig3a + theme(legend.position="none"), prow, nrow=2, labels=c("A","",""), label_size=24, rel_heights=c(1,0.8))
# Add panel A

legend_b = get_legend(fig3b + theme(legend.direction="horizontal", legend.justification="center", legend.box.just="bottom"))
# Extract the legend from one of the plots

p = plot_grid(pcombo, legend_b, ncol=1, rel_heights=c(1, 0.1))
# Add the legend underneath the row we made earlier with 10% of the height of the row.
print(p)

if(savefiles){
  if(color_plots){
    outfile = paste(baseout, ".pdf", sep="")
  }else{
    outfile = paste(baseout, "-grey.png", sep="")
  }
  cat(" -> ", outfile, "\n")
  ggsave(filename=outfile, p, width=10, height=10, units="in")
}
# Save the figure
######################
stop()

















######################
# Some experimenting with other types of plots.
raincloud_theme <- theme(
  text = element_text(size = 10),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  axis.text = element_text(size = 14),
  axis.text.x = element_text(angle = 45, vjust = 0.5),
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 16),
  legend.position = "right",
  plot.title = element_text(lineheight = .8, face = "bold", size = 16),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
  axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"))


test = ggplot(data = sv_alleles, aes(x=Species, y=Length, fill=Species)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .6) +
  geom_point(aes(y = Length, color = Species), 
             position = position_jitter(width = .15), size = .25, alpha = 0.3) +
  stat_boxplot(geom ='errorbar', width=0.1) +
  geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.8, lwd=1) +
  expand_limits(x =c(1.5,2)) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  scale_color_manual(name="", values=c("Macaque"='#490092',"Human"='#920000')) +
  scale_fill_manual(name="", values=c("Macaque"='#490092',"Human"='#920000')) +
  coord_flip() +
  theme_bw()
raincloud_theme


devtools::source_gist("2a1bb0133ff568cbe28d", 
                      filename = "geom_flat_violin.R")
## sourced from github "dgrtwo/geom_flat_violin.R
test = ggplot(data = sv_alleles, 
              mapping = aes(x = Species, 
                            y = Length, 
                            fill = Species)) + 
  geom_flat_violin(scale = "count", 
                   trim = FALSE) + 
  stat_summary(fun.data = mean_sdl, 
               fun.args = list(mult = 1), 
               geom = "pointrange", 
               position = position_nudge(0.05)) + 
  geom_dotplot(binaxis = "y", 
               dotsize = 0.5, 
               stackdir = "down", 
               binwidth = 15, 
               position = position_nudge(-0.025)) + 
  scale_color_manual(name="", values=c("Macaque"='#490092',"Human"='#920000')) +
  scale_fill_manual(name="", values=c("Macaque"='#490092',"Human"='#920000')) +
  theme_bw() + 
  coord_flip() +
  theme(legend.position = "none")



print(test)