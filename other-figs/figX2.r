############################################################
# For macaque paper revisions, 04.20
# Checks distributions of genes overlapping CNVs
# Gregg Thomas
############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(ggplot2)
library(cowplot)
source("../lib/design.r")

cat("----------\n")

############################################################

genHist <- function(df, mode, type) {
  
  if(type=="del"){
    type = "deleted"
    fcol = "#006ddb"
    df$xvar = df$Num.del
    fillstr = "Deletions"
  }else if(type=="dup"){
    type = "duplicated"
    fcol = "#db6d00"
    df$xvar = df$Num.dup
    fillstr = "Duplications"
  }else{
    type = "deleted or duplicated"
    fcol = "#999999"
    df$xvar = df$total
    fillstr = "All"
  }
  
  df = subset(df, xvar > 0)
  xmax = max(df$xvar)
  
  xtitle = paste("# of times ", mode, " ", type, sep="")
  ytitle = paste("# of ", mode, "s", sep="")
  
  p = ggplot(df, aes(x=xvar)) +
    geom_histogram(aes(fill=fillstr), bins=xmax, color="black") +
    xlab(xtitle) +
    ylab(ytitle) +
    scale_x_continuous(breaks=1:xmax) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(name="", values=c("All"='#999999',"Deletions"='#006ddb',"Duplications"='#db6d00')) +
    bartheme() +
    theme(legend.position="none")
  return(p)
}

parseDF <- function(df, mode) {
  df = subset(df, Feature.type == mode)
  df$total = df$Num.del.full + df$Num.del.partial + df$Num.dup.full + df$Num.dup.partial
  df$Num.del = df$Num.del.full + df$Num.del.partial
  df$Num.dup = df$Num.dup.full + df$Num.dup.partial
  return(df)
}

############################################################

mq_gene_counts = read.csv("../data/macaque-cnv-gene-counts.csv", header=T)

mq_genes = parseDF(mq_gene_counts, "gene")
genes_p = genHist(mq_genes, "gene", "all")
print(genes_p)
genes_del_p = genHist(mq_genes, "gene", "del")
print(genes_del_p)
genes_dup_p = genHist(mq_genes, "gene", "dup")
print(genes_dup_p)


cat("Combining gene plots...\n")
prow = plot_grid(genes_del_p, genes_dup_p, align = 'vh', hjust = -1, nrow = 1)
genes_pcombo = plot_grid(genes_p, prow, nrow=2)
print(genes_pcombo)
## Genes


mq_trans = parseDF(mq_gene_counts, "transcript")
trans_p = genHist(mq_trans, "transcript", "all")
print(trans_p)
trans_del_p = genHist(mq_trans, "transcript", "del")
print(trans_del_p)
trans_dup_p = genHist(mq_trans, "transcript", "dup")
print(trans_dup_p)

cat("Combining trancsript plots...\n")
prow = plot_grid(trans_del_p, trans_dup_p, align = 'vh', hjust = -1, nrow = 1)
trans_pcombo = plot_grid(trans_p, prow, nrow=2)
print(trans_pcombo)
## Transcripts


mq_exons = parseDF(mq_gene_counts, "exon")
exons_p = genHist(mq_exons, "exon", "all")
print(exons_p)
exons_del_p = genHist(mq_exons, "exon", "del")
print(exons_del_p)
exons_dup_p = genHist(mq_exons, "exon", "dup")
print(exons_dup_p)


cat("Combining exon plots...\n")
prow = plot_grid(trans_del_p, exons_dup_p, align = 'vh', hjust = -1, nrow = 1)
exons_pcombo = plot_grid(exons_p, prow, nrow=2)
print(exons_pcombo)
## Exons

cat("Combining ALL plots...\n")
full_combo = plot_grid(genes_pcombo, trans_pcombo, exons_pcombo, nrow=3, labels=c("Genes","Transcripts","Exons"), label_size=24, vjust=0.1)
print(full_combo)

cat("Saving...\n")
ggsave(filename="figX2.png", full_combo, width=10, height=20, units="in")


############################################################

# The Chi-square test for ratios of deletions and duplications
cat("\n-------------\n")
cat("Gene overlap Chi-square test:\n")
gene_overlaps = data.frame("Species"=c("Macaque", "Human"), "num.del"=c(952,7804), "num.dup"=c(425,2597))
gene_overlaps_t = t(gene_overlaps[,2:ncol(gene_overlaps)])
colnames(gene_overlaps_t) <- gene_overlaps[,1]
gene_overlaps_chi = chisq.test(gene_overlaps_t)
print(gene_overlaps_chi)
cat("-------------\n")