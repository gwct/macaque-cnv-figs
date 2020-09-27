############################################################
# For macaque paper, 04.19
# Reads SV data files output by sv_parse.py
# Gregg Thomas
############################################################

readSVs <- function() {
  cat("Reading human data...\n")
  sv_all_hu = read.csv("../cnv-calls/brandler-cnvs.csv", header=T, comment.char="#")
  #sv_all_hu = read.csv(unz("../data/brandler-cnvs.zip", "brandler-cnvs.csv"), header = TRUE, sep = ",") 
  sv_all_hu$Species = "Human"
  sv_all_hu$Maternal.age = sv_all_hu$Maternal.age / 12
  sv_all_hu$Paternal.age = sv_all_hu$Paternal.age / 12
  # Converts ages from months to years.

  cat("Reading macaque data...\n")
  sv_all_mq = read.csv("../cnv-calls/macaque-cnvs.csv", header=T, comment.char="#")
  sv_all_mq$Species = "Macaque"
  
  return(list(sv_all_mq, sv_all_hu))
}
