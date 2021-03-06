############################################################
# For macaque paper, 08.19
# Makes length histograms for human and macaque SVs
# Gregg Thomas
############################################################

binSVs <- function(svs){
# Function to pre-bin SVs by length
  binned_svs = data.frame("bin"=seq(from=0,to=5000,by=100), count=0)
  for(i in 1:nrow(svs)) {
    row <- svs[i,]
    bin = floor(row$Length/100)*100
    binned_svs[binned_svs$bin==bin,]$count = binned_svs[binned_svs$bin==bin,]$count + 1
  }
  return(binned_svs) 
}

############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(ggplot2)
library(reshape2)
library(cowplot)
library(ggridges)
library(OneR)
library(ggforce)
library(ggbeeswarm)

source("../lib/read_svs.r")
source("../lib/filter_svs.r")
source("../lib/subset_svs.r")
source("../lib/design.r")

cat("----------\n")

############################################################

savefiles = T
color_plots = F
rm_alus = T
# Run options

minlen = F
maxlen = 5000
# CNV length cutoffs

baseout = "fig2"
human_full = F
macaque_filter = T
figs4_opt = F
figs5_opt = F
if(figs4_opt){
  human_full = T
  baseout = paste(baseout, "_S4", sep="")
}else if(figs5_opt){
  macaque_filter = F
  baseout = paste(baseout, "_S5", sep="")
}

if(rm_alus){
  baseout = paste(baseout, "_noalu", sep="")
}
# Supplemental figure options

############################################################

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

if(rm_alus){
  hu_svs = subset(hu_svs, Length < 275 | Length > 325)
  if(!figs5_opt){
    mq_svs = subset(mq_svs, Length < 275 | Length > 325)
  }
}
# For Alu stuff

sv_alleles = rbind(hu_svs, mq_svs)
# Combine the two datasets

#test_count = mq_events %>% group_by(Individual) %>% count(Individual)
# Count number of events per individual
######################

######################
# SV length histogram -- Fig 2A
cat("SV length distribution...\n")

fig2a = ggplot(sv_alleles, aes(x=Length, fill=Species, color=Species)) + 
  geom_histogram(alpha=0.8, position="identity", bins=50, size=0.5) +
  #geom_density(alpha=0.3) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="CNV length", y="# CNVs") +
  bartheme()
if(color_plots){
  fig2a = fig2a + scale_fill_manual(name="", values=c("Macaque"='#490092',"Human"='#920000'))
}else{
  fig2a = fig2a + scale_fill_manual(name="", values=c("Macaque"='#5c5c5c',"Human"=NA))
  }
fig2a = fig2a + scale_color_manual(name="", values=c("Macaque"=NA,"Human"='#000000'))
print(fig2a)
# Full length histogram

#if(savefiles){
#  outfile = "fig2a.pdf"
#  cat(" -> ", outfile, "\n")
#  ggsave(filename=outfile, fig2a, width=6, height=4, units="in")
#}

fig2a_box = ggplot(sv_alleles, aes(x=Species, y=Length, fill=Species)) + 
  geom_quasirandom(size=2, alpha=0.7, width=0.25, color="#d3d3d3") +
  geom_boxplot(outlier.shape=NA, alpha=0.3, width=0.5, color="#000000") +
  labs(y="CNV Length", x="All CNVs") +
  bartheme() +
  theme(legend.position="none")
if(color_plots){
  fig2a_box = fig2a_box + scale_fill_manual(name="", values=c("Macaque"='#490092', "Human"='#920000'))
}else{
  fig2a_box = fig2a_box + scale_fill_manual(name="", values=c("Macaque"='#5c5c5c', "Human"=NA))
}
print(fig2a_box)
# Full lenght boxplot
#ggsave(filename="figS3.png", fig2a, width=6, height=4, units="in")

cat(" -> SV length KS test...\n")
length_ks = ks.test(hu_svs$Length, mq_svs$Length)
print(length_ks)
# Full length KS test
######################

######################
# SV deletion length histogram
cat("SV deletion length distribution...\n")
mq_dels = subset(mq_svs, Type=="<DEL>")
hu_dels = subset(hu_svs, Type=="<DEL>")
sv_dels = rbind(hu_dels, mq_dels)

fig2b = ggplot(sv_dels, aes(x=Length, fill=Species, color=Species)) + 
  geom_histogram(alpha=0.8, position="identity", bins=50, size=0.5) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="Deletion length", y="# Deletions") +
  bartheme()
if(color_plots){
  fig2b = fig2b + scale_fill_manual(name="", values=c("Macaque"='#490092',"Human"='#920000'))
}else{
  #fig2b = fig2b + scale_fill_grey(name="", labels=c("Human","Macaque"))
  fig2b = fig2b + scale_fill_manual(name="", values=c("Macaque"='#5c5c5c',"Human"=NA))
}
fig2b = fig2b + scale_color_manual(name="", values=c("Macaque"=NA,"Human"='#000000'))
print(fig2b)
# Deletion length histogram

#if(savefiles){
#  outfile = "fig2b.pdf"
#  ggsave(filename=outfile, fig2b, width=6, height=4, units="in")
#}

fig2b_box = ggplot(sv_dels, aes(x=Species, y=Length, fill=Species)) + 
  geom_quasirandom(size=2, alpha=0.7, width=0.25, color="#d3d3d3") +
  geom_boxplot(outlier.shape=NA, alpha=0.3, width=0.5, color="#000000") +
  labs(y="CNV Length", x="Deletions") +
  bartheme() +
  theme(legend.position="none")
if(color_plots){
  fig2b_box = fig2b_box + scale_fill_manual(name="", values=c("Macaque"='#490092', "Human"='#920000'))
}else{
  fig2b_box = fig2b_box + scale_fill_manual(name="", values=c("Macaque"='#5c5c5c', "Human"=NA))
}
print(fig2b_box)
# Deletion length boxplot

cat(" -> SV deletion length KS test...\n")
del_ks = ks.test(mq_dels$Length, hu_dels$Length)
print(del_ks)
# Deletion length KS test
######################

######################
# SV duplication length histogram
cat("SV deletion length distribution...\n")
mq_dups = subset(mq_svs, Type=="<DUP>")
hu_dups = subset(hu_svs, Type=="<DUP>")
sv_dups = rbind(hu_dups, mq_dups)

fig2c = ggplot(sv_dups, aes(x=Length, fill=Species, color=Species)) + 
  geom_histogram(alpha=0.8, position="identity", bins=50, size=0.5) +
  #scale_fill_manual(name="", labels=c("Macaque","Human"), values=c("#006ddb","#db6d00")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x="Duplication length", y="# Duplications") +
  bartheme()
if(color_plots){
  fig2c = fig2c + scale_fill_manual(name="", values=c("Macaque"='#490092',"Human"='#920000'))
}else{
  fig2c = fig2c + scale_fill_manual(name="", values=c("Macaque"='#5c5c5c',"Human"=NA))
}
fig2c = fig2c + scale_color_manual(name="", values=c("Macaque"=NA,"Human"='#000000'))
print(fig2c)
# Duplication length histogram

#if(savefiles){
#  outfile = "fig2c.pdf"
#  ggsave(filename=outfile, fig2c, width=6, height=4, units="in")
#}
# SV del length plot

fig2c_box = ggplot(sv_dups, aes(x=Species, y=Length, fill=Species)) + 
  geom_quasirandom(size=2, alpha=0.7, width=0.25, color="#d3d3d3") +
  geom_boxplot(outlier.shape=NA, alpha=0.3, width=0.5, color="#000000") +
  labs(y="CNV Length", x="Duplications") +
  bartheme() +
  theme(legend.position="none")
if(color_plots){
  fig2c_box = fig2c_box + scale_fill_manual(name="", values=c("Macaque"='#490092', "Human"='#920000'))
}else{
  fig2c_box = fig2c_box + scale_fill_manual(name="", values=c("Macaque"='#5c5c5c', "Human"=NA))
}
print(fig2c_box)
# Duplication length boxplot

cat(" -> SV duplication length KS test...\n")
dup_ks = ks.test(mq_dups$Length, hu_dups$Length)
print(dup_ks)
# Duplication length KS test
######################

######################
# Combine plots for figure

cat("Combining histograms...\n")
prow = plot_grid(fig2b + theme(legend.position="none"),
                 fig2c + theme(legend.position="none"),
                 align = 'vh',
                 labels = c("B", "C"),
                 label_size = 24,
                 hjust = -1,
                 nrow = 1)
# Combine panels B and C

pcombo = plot_grid(fig2a + theme(legend.position="none"), prow, nrow=2, labels=c("A","",""), label_size=24, rel_heights=c(1,0.8))
# Add panel A

legend_b = get_legend(fig2b + theme(legend.direction="horizontal", legend.justification="center", legend.box.just="bottom"))
# Extract the legend from one of the plots

p = plot_grid(pcombo, legend_b, ncol=1, rel_heights=c(1, 0.1))
# Add the legend underneath the row we made earlier with 10% of the height of the row.
print(p)

cat("Combining boxplots...\n")
p_box = plot_grid(fig2a_box, fig2b_box, fig2c_box, nrow=1, 
                 align = 'vh',
                 labels = c("A", "B", "C"),
                 label_size = 24,
                 hjust = -1)
# Combine boxplots

if(savefiles){
  if(color_plots){
    outfile = paste(baseout, ".pdf", sep="")
    outfile_box = paste(baseout, "-box.pdf", sep="")
  }else{
    outfile = paste(baseout, "-grey.pdf", sep="")
    outfile_box = paste(baseout, "-box-grey.pdf", sep="")
  }
  cat(" -> ", outfile, "\n")
  ggsave(filename=outfile, p, width=10, height=10, units="in")
  ggsave(filename=outfile_box, p_box, width=10, height=4, units="in")
}
# Save the figure
######################
stop("End of plots.")
















######################
# SV macaque deletion vs. duplication length histogram
fig2a = ggplot(mq_svs, aes(x=Length, fill=Type, color=Type)) + 
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
  fig2a = fig2a + scale_fill_manual(name="", values=c("<DEL>"='#490092',"<DUP>"='#920000'))
}else{
  #fig2a = fig2a + scale_fill_grey(name="", labels=c("Human","Macaque"))
  #fig2a = fig2a + scale_fill_manual(name="", values=c("Macaque"='#d6d6d6',"Human"='#5c5c5c'))
  fig2a = fig2a + scale_fill_manual(name="", values=c("<DEL>"='#5c5c5c',"<DUP>"=NA))
}
fig2a = fig2a + scale_color_manual(name="", values=c("<DEL>"=NA,"<DUP>"='#000000'))

print(fig2a)

cat(" -> SV macaque deletion vs. duplication length KS test...\n")
mq_ks = ks.test(mq_dels$Length, mq_dups$Length)
print(mq_ks)
# SV del length KS test
# SV macaque deletion vs. duplication length histogram
######################


######################
# Some experimenting with other types of plots.


fig2a = ggplot(sv_alleles, aes(x=Species, y=Length)) +
  geom_violin(fill="#999999") + 
  #geom_sina(size=0.5) +
  labs(x="", y="CNV length") +
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

print(fig2a)
## Violin plot

#cat("----------\nPre-binning macaque data...\n")
#mq_bins = binSVs(mq_svs)
#mq_bins$Species = "Macaque"
#mq_bins$count = -mq_bins$count

#cat("----------\nPre-binning human data...\n")
#hu_bins = binSVs(hu_svs)
#hu_bins$Species = "Human"
#sv_bins = rbind(hu_bins, mq_bins)

#fig2a = ggplot(sv_bins, aes(x=bin, y=count, fill=Species)) + 
fig2a = ggplot(sv_bins, aes(x = bin, y = ifelse(Species == "Macaque",-1, 1)*count, fill = Species)) + 
  geom_col() +
  scale_x_continuous(name = "CNV length", limits = c(0, 5000), expand = c(0, 0)) +
  scale_y_continuous(name = "Count", limits=c(-3000,3500), breaks = 1000*(-2:3), labels=c("2000","1000","0","1000","2000","3000")) +
  #scale_fill_manual(values = c("#D55E00", "#0072B2"), guide = "none") +
  #draw_text(x = 70, y = -39, "male", hjust = 0) +
  #draw_text(x = 70, y = 21, "female", hjust = 0) +
  #coord_flip() +
  theme_minimal_grid() +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=16), 
        axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0),color="black"), 
        axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0),color="black"),
        axis.line=element_line(colour='#595959',size=0.75),
        #axis.ticks=element_line(colour="#595959",size = 1),
        #axis.ticks.length=unit(0.2,"cm"),
        legend.position="right",
        legend.key.width = unit(0.75,  unit = "cm"),
        legend.spacing.x = unit(0.25, 'cm'),
        legend.title = element_blank(),
        legend.text=element_text(size=12),
        plot.title = element_text(hjust=0.5, size=20),
        plot.margin = unit(c(1,1,1,1), "cm")
  )

if(color_plots){
  fig2a = fig2a + scale_fill_manual(name="", values=c("Macaque"='#490092',"Human"='#920000'))
}else{
  #fig2a = fig2a + scale_fill_grey(name="", labels=c("Human","Macaque"))
  #fig2a = fig2a + scale_fill_manual(name="", values=c("Macaque"='#d6d6d6',"Human"='#5c5c5c'))
  fig2a = fig2a + scale_fill_manual(name="", values=c("Macaque"='#5c5c5c',"Human"=NA))
}
fig2a = fig2a + scale_color_manual(name="", values=c("Macaque"=NA,"Human"='#000000'))

print(fig2a)
# Top and bottom plot.




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
# Raincould plot