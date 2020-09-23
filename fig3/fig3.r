############################################################
# For macaque paper, 08.19
# Makes de novo SV vs. paternal age plots
# Gregg Thomas
############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(ggplot2)
library(reshape2)
library(plyr)
library(ggpubr)
library(cowplot)

source("../lib/read_svs.r")
source("../lib/filter_svs.r")
source("../lib/subset_svs.r")

cat("----------\n")

############################################################

savefiles = T
color_plots = F
supp = T
rm_alus = T
# Run options

minlen = F
maxlen = 100000
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
mq_events = mqr[[1]];

cat("----------\nSubsetting human data...\n")
hur = subsetSVs(hu_events)
hu_events = hur[[1]];
# Subset data

######################
# Macaque de novos

cat("Counting denovos in macaques...\n")
f1s = c(39243,39317,39313,39242,39226,39230,39234,39316,39227,39231,39318,39314,39229,39319)
# A vector of the focal individuals

mq_denovo = subset(mq_events, Denovo=="Y")
mq_denovo_counts = count(mq_denovo, vars="Individual")
mq_denovo = merge(mq_denovo, mq_denovo_counts, by="Individual")
mq_denovo_counts = ddply(mq_denovo, .(Individual), summarize,
                  Num.dsv=mean(freq),
                  Maternal.age=mean(Maternal.age),
                  Paternal.age=mean(Paternal.age))
# Get the de novos.

cat(" -> Adding F1s with no denovos...\n")
for(ind in f1s){
  if(!ind %in% mq_denovo_counts$Individual){
    cur_mat = mq_events$Maternal.age[mq_events$Individual==ind][1]
    cur_pat = mq_events$Paternal.age[mq_events$Individual==ind][1]
    cur_row = data.frame("Individual"=ind, "Num.dsv"=0, "Maternal.age"=cur_mat, "Paternal.age"=cur_pat)
    mq_denovo_counts = rbind(mq_denovo_counts,cur_row)
  }
}
# Adding the trios with 0 counts back into the frame

cat(" -> Correlating age with denovos...\n")
pat_fit = lm(mq_denovo_counts$Num.dsv ~ mq_denovo_counts$Paternal.age)
mat_fit = lm(mq_denovo_counts$Num.dsv ~ mq_denovo_counts$Maternal.age)

mq_denovo_counts$Species = "Macaque"
mq_denovo_counts$Species = factor(mq_denovo_counts$Species, levels=c("Macaque", "Human"))

cat(" -> Plotting denovos...\n")
fig3a = ggplot(mq_denovo_counts, aes(Paternal.age, Num.dsv, color=Species)) +
  geom_point(size=3, alpha=0.5) +
  #geom_smooth(method="glm", method.args=list(family='poisson'), fullrange=T, size=0.75, linetype="dashed", alpha=0) +
  geom_smooth(method="glm", fullrange=T, size=0.75, linetype="dashed", alpha=0) +
  ggtitle("") +
  scale_y_continuous(limits=c(0, 5)) +
  labs(x="Paternal age (years)", y="# de novo CNVs") +
  bartheme()
  # theme_classic() +
  # theme(axis.text=element_text(size=12), 
  #       axis.title=element_text(size=16), 
  #       axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0),color="black"), 
  #       axis.title.x=element_text(margin=margin(t=0,r=0,b=0,l=0),color="black"),
  #       axis.line=element_line(colour='#595959',size=0.75),
  #       axis.ticks=element_line(colour="#595959",size = 1),
  #       axis.ticks.length=unit(0.2,"cm"),
  #       legend.position="right",
  #       legend.key.width = unit(0.75,  unit = "cm"),
  #       legend.spacing.x = unit(0.25, 'cm'),
  #       legend.title = element_blank(),
  #       legend.text=element_text(size=12),
  #       plot.title = element_text(hjust=0.5, size=20),
  #       plot.margin = unit(c(1,1,0,1), "cm")
  # )

if(color_plots){
  fig3a = fig3a + scale_color_manual(name="", values=c('#490092','#920000'), labels=c("Macaque","Human"), drop=FALSE)
}else{
  fig3a = fig3a + scale_color_grey(name="", labels=c("Macaque","Human"), drop=FALSE) +
    ggtitle("Macaque") +
    theme(plot.margin = unit(c(1,1,1,1), "cm"))
}

print(fig3a)
# Macaque de novo plot

mq_rate = mean(mq_denovo_counts$Num.dsv / 2)
mq_se = ((sd(mq_denovo_counts$Num.dsv) / 2) / sqrt(length(mq_denovo_counts$Individual))) * 1.96
cat(" -> Macaque CNV rate per haploid generation: ", mq_rate, " (95% CI +/-", mq_se, ")\n", sep="")
# Macaque rate per haploid generation

cat(" -> Macaque CNV paternal age correlation: r2 = ", summary(pat_fit)$r.squared, ", df = ", summary(pat_fit)$df[2], ", p = ", summary(pat_fit)$coefficients[2,4], "\n", sep="")
# Paternal age regression

# Macaque de novos
######################

######################
# Human de novos (Brandler's count)
cat("Reading Brandler data...\n")
brandler_svs = read.csv("../data/brandler-denovo.csv", header=TRUE)
brandler_info = read.csv("../data/brandler-samples.csv", header=TRUE)
names(brandler_info)[2] = "ID"
brandler_info = subset(brandler_info, Relationship=="Proband" | Relationship=="Sibling")
brandler_svs$SVLEN = as.numeric(brandler_svs$SVLEN)
# Have to read the Brandler files

if(rm_alus){
  brandler_svs = subset(brandler_svs, SVLEN < 275 | SVLEN > 325)
}
# For Alu stuff

brandler_svs = subset(brandler_svs, SVTYPE %in% c("DEL","DUP"))
if(supp){
  brandler_svs = subset(brandler_svs, VALIDATION == 1)
}
  
brandler_counts = count(brandler_svs, vars="ID")
brandler_svs = merge(brandler_svs, brandler_counts, by="ID")

brandler_counts = merge(brandler_counts, brandler_info[, c("ID", "Mother_chrono_age")], by="ID")
brandler_counts = merge(brandler_counts, brandler_info[, c("ID", "Father_chrono_age")], by="ID")
names(brandler_counts)[2] = "Num.dsv"

for(ind in brandler_info$ID){
  if(!ind %in% brandler_counts$ID){
    cur_mat = brandler_info$Mother_chrono_age[brandler_info$ID==ind][1]
    cur_pat = brandler_info$Father_chrono_age[brandler_info$ID==ind][1]
    cur_row = data.frame("ID"=ind, "Num.dsv"=0, "Mother_chrono_age"=cur_mat, "Father_chrono_age"=cur_pat)
    brandler_counts = rbind(brandler_counts,cur_row)
  }
}
# Adding the trios with 0 counts back into the frame

brandler_counts$Mother_chrono_age = brandler_counts$Mother_chrono_age / 12
brandler_counts$Father_chrono_age = brandler_counts$Father_chrono_age / 12

cat(" -> Correlating age with denovos...\n")
pat_fit_h = lm(brandler_counts$Num.dsv ~ brandler_counts$Father_chrono_age)
mat_fit_h = lm(brandler_counts$Num.dsv ~ brandler_counts$Mother_chrono_age)

cat(" -> Plotting human denovos...\n")
fig3b = ggplot(brandler_counts, aes(Father_chrono_age, Num.dsv, color="Human")) +
  geom_point(size=3, alpha=0.5) +
  #geom_smooth(method="glm", method.args=list(family='poisson'), fullrange=T, size=0.75, linetype="dashed", alpha=0) +
  geom_smooth(method="glm", fullrange=T, size=0.75, linetype="dashed", alpha=0) +
  ggtitle("") +
  scale_y_continuous(limits=c(0, 5)) +
  labs(x="Paternal age (years)", y="") +
  bartheme()
  # theme_classic() +
  # theme(axis.text=element_text(size=12), 
  #       axis.title=element_text(size=16), 
  #       axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0),color="black"), 
  #       axis.title.x=element_text(margin=margin(t=0,r=0,b=0,l=0),color="black"),
  #       axis.line=element_line(colour='#595959',size=0.75),
  #       axis.ticks=element_line(colour="#595959",size = 1),
  #       axis.ticks.length=unit(0.2,"cm"),
  #       legend.position="right",
  #       legend.key.width = unit(0.75,  unit = "cm"),
  #       legend.spacing.x = unit(0.25, 'cm'),
  #       legend.title = element_blank(),
  #       legend.text=element_text(size=12),
  #       plot.title = element_text(hjust=0.5, size=20),
  #       plot.margin = unit(c(1,1,0,1), "cm")
  # )

if(color_plots){
  fig3b = fig3b + scale_color_manual(name="", breaks=c("Macaque","Human"), values=c("Macaque"='#490092',"Human"='#920000'), drop=FALSE)
}else{
  fig3b = fig3b + scale_color_grey(name="", labels=c("Macaque","Human"), drop=FALSE) +
    ggtitle("Human") +
    theme(plot.margin = unit(c(1,1,1,1), "cm"))
}

print(fig3b)

hu_rate = mean(brandler_counts$Num.dsv / 2)
hu_se = ((sd(brandler_counts$Num.dsv) / 2) / sqrt(length(brandler_counts$ID))) * 1.96
cat(" -> Human CNV rate per haploid generation: ", hu_rate, " (95% CI +/-", hu_se, ")\n", sep="")
# Macaque rate per haploid generation

cat(" -> Human CNV paternal age correlation: r2 = ", summary(pat_fit_h)$r.squared, ", df = ", summary(pat_fit_h)$df[2], ", p = ", summary(pat_fit_h)$coefficients[2,4], "\n", sep="")
# Paternal age regression

# Human de novos (Brandler's count)
######################

######################
# Combine plots for figure
cat("Combining plots...\n")
prow = plot_grid(fig3a + theme(legend.position="none"),
                 fig3b + theme(legend.position="none"),
                 align = 'vh',
                 labels = c("A", "B"),
                 label_size = 24,
                 hjust = -1,
                 nrow = 1)


legend_b = get_legend(fig3a + theme(legend.direction="horizontal", legend.justification="center", legend.box.just="bottom"))
# Extract the legend from one of the plots

if(color_plots){
  p = plot_grid(prow, legend_b, ncol=1, rel_heights=c(1, 0.1))
}else{
  p = prow
}
# For color plots: Add the legend underneath the row we made earlier with 10% of the height of the row
# Otherwise, don't add

print(p)

if(savefiles){
  outfile = "fig3"
  if(supp){
    outfile = paste(outfile, "_S7", sep="")
  }
  if(!color_plots){
    outfile = paste(outfile, "-grey", sep="")
  }
  outfile = paste(outfile, ".png", sep="")
  cat(" -> Outfile: ", outfile, "\n")
  ggsave(filename=outfile, p, width=10, height=5, units="in")
}
# Save the figure
######################





