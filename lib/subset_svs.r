############################################################
# For macaque paper, 04.19
# Subset and count SV calls
# Gregg Thomas
############################################################

subsetSVs <- function(svdf){
  svdf_focal = subset(svdf, F1=="Y")
  svdf_denovo = subset(svdf, Denovo=="Y")
  
  #cat("Counts for ALL SV calls:\n")
  #countSV(svdf, svdf_focal, svdf_denovo)
  
  svdf_uniq = subset(svdf, select=c("SV.key", "Type", "Chromosome", "Pos", "Length", "Genes", "Num.genes", "Species", "SV.freq"))
  svdf_uniq = unique(svdf_uniq)
  
  svdf_focal_uniq = subset(svdf_focal, select=c("SV.key", "Type", "Chromosome", "Pos", "Length", "Genes", "Num.genes", "Species", "SV.freq"))
  svdf_focal_uniq = unique(svdf_focal_uniq)
  
  svdf_denovo_uniq = subset(svdf_denovo, select=c("SV.key", "Type", "Chromosome", "Pos", "Length", "Genes", "Num.genes", "Species", "SV.freq"))
  svdf_denovo_uniq = unique(svdf_denovo_uniq)
  
  #cat("Counts for SV alleles:\n")
  #countSV(svdf_uniq, svdf_focal_uniq, svdf_denovo_uniq)
  
  #print(names(svdf_uniq))
  
  results = list(svdf, svdf_focal, svdf_denovo, svdf_uniq, svdf_focal_uniq, svdf_denovo_uniq)
  return(results)
}