############################################################
# For macaque paper, 09.20
# Just does some of the Chi-squared tests for gene counts.
# Gregg Thomas
############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

cat("----------\n\n")

############################################################

cat("HUMAN AND MACAQUE GENE RATIOS BY CNV TYPE\n\n")
gene_ratios = matrix(c(1212,5483,226,1132), nrow=2, dimnames=list(c("macaque", "human"), c("del","dup")))
gene_ratio_chi = chisq.test(gene_ratios)
print(gene_ratio_chi)
cat("----------\n\n")
# 
cat("HUMAN AND MACAQUE 10KB UPSTREAM RATIOS BY CNV TYPE\n\n")
up_ratios = matrix(c(469,3132,147,647), nrow=2, dimnames=list(c("macaque", "human"), c("del","dup")))
up_ratio_chi = chisq.test(up_ratios)
print(up_ratio_chi)
cat("----------\n\n")
# 
cat("HUMAN AND MACAQUE 10KB DOWNSTREAM RATIOS BY CNV TYPE\n\n")
down_ratios = matrix(c(537,3174,149,693), nrow=2, dimnames=list(c("macaque", "human"), c("del","dup")))
down_ratio_chi = chisq.test(down_ratios)
print(down_ratio_chi)
cat("----------\n\n")
# 
cat("HUMAN AND MACAQUE TRANSCRIPT RATIOS BY CNV TYPE\n\n")
transcript_ratios = matrix(c(2665,22246,476,4176), nrow=2, dimnames=list(c("macaque", "human"), c("del","dup")))
transcript_ratio_chi = chisq.test(transcript_ratios)
print(transcript_ratio_chi)
cat("----------\n\n")

cat("HUMAN AND MACAQUE EXON RATIOS BY CNV TYPE\n\n")
exon_ratios = matrix(c(936,7597,425,2587), nrow=2, dimnames=list(c("macaque", "human"), c("del","dup")))
exon_ratio_chi = chisq.test(exon_ratios)
print(exon_ratio_chi)
cat("----------\n\n")
# Numbers for these test come from adding full and partial events in TABLE 1
cat("----------\n\n")


cat("MACAQUE TOTAL CNVS BY TYPE TO CNVS OVERLAPPING EXONS BY TYPE\n\n")
mq_exon_ratio = matrix(c(432,90,3214,291), nrow=2, dimnames=list(c("total", "exon"), c("dup","del")))
mq_exon_chi = chisq.test(mq_exon_ratio)
print(mq_exon_chi)
cat("----------\n\n")
# 
cat("HUMAN TOTAL CNVS BY TYPE TO CNVS OVERLAPPING EXONS BY TYPE\n\n")
hu_exon_ratio = matrix(c(1784,443,13745,2267), nrow=2, dimnames=list(c("total", "exon"), c("dup","del")))
hu_exon_chi = chisq.test(hu_exon_ratio)
print(hu_exon_chi)
cat("----------\n\n")
# Numbers for these test come from adding full and partial events in TABLE 2 and from cnv_counts.r