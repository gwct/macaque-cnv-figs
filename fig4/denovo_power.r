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

savefiles = F
color_plots = T
supp = T
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

# Macaque de novos
######################

######################
# Human de novos (Brandler's count)
cat("Reading Brandler data...\n")
brandler_svs = read.csv("../data/brandler-denovo.csv", header=TRUE)
brandler_info = read.csv("../data/brandler-samples.csv", header=TRUE)
names(brandler_info)[2] = "ID"
brandler_info = subset(brandler_info, Relationship=="Proband" | Relationship=="Sibling")
# Have to read the Brandler files

brandler_svs = subset(brandler_svs, SVTYPE %in% c("DEL","DUP"))
if(supp){
  print("HIHI")
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

# Human de novos (Brandler's count)
######################

