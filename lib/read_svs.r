############################################################
# For macaque paper, 04.19
# Reads SV data files output by sv_parse.py
# Gregg Thomas
############################################################

parseCNVDF <- function(df) {
  df$Genes.total = df$gene.del.full.count + df$gene.del.partial.count + df$gene.dup.full.count + df$gene.dup.partial.count
  df$Genes.del = df$gene.del.full.count + df$gene.del.partial.count
  df$Genes.dup = df$gene.dup.full.count + df$gene.dup.partial.count

  df$Transcripts.total = df$transcript.del.full.count + df$transcript.del.partial.count + df$transcript.dup.full.count + df$transcript.dup.partial.count
  df$Transcripts.del = df$transcript.del.full.count + df$transcript.del.partial.count
  df$Transcripts.dup = df$transcript.dup.full.count + df$transcript.dup.partial.count

  df$Exons.total = df$exon.del.full.count + df$exon.del.partial.count + df$exon.dup.full.count + df$exon.dup.partial.count
  df$Exons.del = df$exon.del.full.count + df$exon.del.partial.count
  df$Exons.dup = df$exon.dup.full.count + df$exon.dup.partial.count
  return(df)
}

readSVs <- function() {
  cat("Reading human data...\n")
  #sv_all_hu = read.csv("../data/brandler-cnvs.csv")
  sv_all_hu = read.csv(unz("../data/brandler-cnvs.zip", "brandler-cnvs.csv"), header = TRUE, sep = ",") 
  sv_all_hu$Species = "Human"
  sv_all_hu$Maternal.age = sv_all_hu$Maternal.age / 12
  sv_all_hu$Paternal.age = sv_all_hu$Paternal.age / 12
  # Converts ages from months to years.

  cat("Reading macaque data...\n")
  sv_all_mq = read.csv("../data/macaque-cnvs.csv")
  sv_all_mq$Species = "Macaque"
  
  return(list(sv_all_mq, sv_all_hu))
}

readGenes <- function() {
  cat("Reading (filtered) human gene data...\n")
  hu_cnv_genes = read.csv("../data/brandler-cnv-gene-overlaps.csv", header = TRUE, sep = ",")
  hu_cnv_genes = parseCNVDF(hu_cnv_genes)
  hu_genes = read.csv("../data/brandler-cnv-gene-counts.csv", header = TRUE, sep = ",")

  cat("Reading (filtered) macaque gene data...\n")
  mq_cnv_genes = read.csv("../data/macaque-cnv-gene-overlaps.csv", header = TRUE, sep = ",")
  mq_cnv_genes = parseCNVDF(mq_cnv_genes)
  mq_genes = read.csv("../data/macaque-cnv-gene-counts.csv", header = TRUE, sep = ",")

  return(list(hu_cnv_genes, hu_genes, mq_cnv_genes, mq_genes))
}