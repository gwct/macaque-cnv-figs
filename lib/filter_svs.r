############################################################
# For macaque paper, 04.19
# Filters SV calls
# Gregg Thomas
############################################################

filterSVs <- function(sv_list, minlen, maxlen, human_full=F, macaque_filter=T){
  mq_events = sv_list[[1]]; hu_events = sv_list[[2]];
  # Unpack the data list
  
  gtcnv_malesex_filter = 8
  gtcnv_nodup_filter = 12
  gtcnv_dup_filter = 8
  svtyper_filter = 100
  dhffc_filter = 0.7
  dhbfc_filter = 1.3
  qual_filter = 100
  freq_filter = 0.95
  # The filter thresholds
  
  cat("Filtering data:\n")
  
  cat(" -> Removing human SVs present at frequency >", as.character(freq_filter), "\n")
  hu_events = subset(hu_events, SV.freq <= freq_filter)
  
  cat(" -> For human gtCNV male sex chromosome SVs, QUAL >=", as.character(gtcnv_malesex_filter), "\n")
  cat(" -> For human gtCNV non-duplication QUAL >=", as.character(gtcnv_nodup_filter), "\n")
  cat(" -> For human gtCNV duplications QUAL >=", as.character(gtcnv_dup_filter), "\n")
  cat(" -> For human SVtyper QUAL >=", as.character(svtyper_filter), "\n")
  
  sv_gtcnv = subset(hu_events, Genotyper == "gtCNV")
  
  sv_gtcnv_malesex = subset(sv_gtcnv, Sex=="Male" & Chromosome %in% c("Y","X"))
  sv_gtcnv_malesex_sub = subset(sv_gtcnv_malesex, Qual >= gtcnv_malesex_filter)
  
  sv_gtcnv_auto = subset(sv_gtcnv, Sex=="Female" | Sex=="Male" & !Chromosome %in% c("Y","X"))
  
  sv_gtcnv_dup = subset(sv_gtcnv_auto, Type == "<DUP>")
  sv_gtcnv_dup_sub = subset(sv_gtcnv_dup, Qual >= gtcnv_dup_filter)
  
  sv_gtcnv_nodup = subset(sv_gtcnv_auto, Type != "<DUP>")
  sv_gtcnv_nodup_sub = subset(sv_gtcnv_nodup, Qual >= gtcnv_nodup_filter)
  
  sv_svtyper = subset(hu_events, Genotyper == "SVtyper" | Genotyper == "Both")
  sv_svtyper_sub = subset(sv_svtyper, Qual >= svtyper_filter)
  
  if(human_full){
    hu_events = rbind(sv_gtcnv_dup_sub, sv_gtcnv_nodup_sub, sv_gtcnv_malesex_sub, sv_svtyper_sub)
    # All human calls
  }else{
    hu_events = sv_svtyper_sub
    # Only SVtyper human calls
  }
  
  if(maxlen){
    cat(" -> Maximum SV length: ", maxlen, "\n")
    hu_events_over = subset(hu_events, Length > maxlen)
    cat(" --> Human events over max len:", length(hu_events_over$Type), "\n")
    cat(" --> Human alleles over max len:", length(unique(hu_events_over$SV.key)), "\n")
    hu_events = subset(hu_events, Length <= maxlen)
  }
  
  if(minlen){
    cat(" -> Minimum SV length: ", minlen, "\n")
    hu_events_under = subset(hu_events, Length < minlen)
    cat(" --> Human events under min len:", length(hu_events_under$Type), "\n")
    cat(" --> Human alleles under min len:", length(unique(hu_events_under$SV.key)), "\n")
    hu_events = subset(hu_events, Length >= minlen)
  }
  
  cat(" -> Removing human inversions...\n")
  hu_invs = subset(hu_events, Type == "<INV>")
  cat(" --> ", length(hu_invs$Type), " human inversions removed\n")
  hu_events = subset(hu_events, Type != "<INV>")
  # Filter human data
  ##########
  
  if(maxlen){
    cat(" -> Maximum SV length: ", maxlen, "\n")
    events_before = length(mq_events$SV.key)
    alleles_before = length(unique(mq_events$SV.key))
    
    mq_events = subset(mq_events, Length <= maxlen)
    
    events_after = length(mq_events$SV.key)
    alleles_after = length(unique(mq_events$SV.key))
    events_removed = events_before - events_after
    alleles_removed = alleles_before - alleles_after
    cat(" -->", events_removed, "SVs removed at", alleles_removed, "sites.\n")
  }
  
  if(minlen){
    cat(" -> Minimum SV length: ", minlen, "\n")
    events_before = length(mq_events$SV.key)
    alleles_before = length(unique(mq_events$SV.key))
    
    mq_events = subset(mq_events, Length >= minlen)
    
    events_after = length(mq_events$SV.key)
    alleles_after = length(unique(mq_events$SV.key))
    events_removed = events_before - events_after
    alleles_removed = alleles_before - alleles_after
    cat(" -->", events_removed, "SVs removed at", alleles_removed, "sites.\n")
  }
  # Length filter
  
  cat(" -> Removing macaque SVs present at frequency >", as.character(freq_filter), "\n")
  events_before = length(mq_events$SV.key)
  alleles_before = length(unique(mq_events$SV.key))
  
  mq_events = subset(mq_events, SV.freq <= freq_filter)
  
  events_after = length(mq_events$SV.key)
  alleles_after = length(unique(mq_events$SV.key))
  events_removed = events_before - events_after
  alleles_removed = alleles_before - alleles_after
  cat(" -->", events_removed, "SVs removed at", alleles_removed, "sites.\n")
  # Allele frequency filter
  
  if(macaque_filter){
    cat(" -> Minimum SV Qual >=", qual_filter, "\n")
    events_before = length(mq_events$SV.key)
    alleles_before = length(unique(mq_events$SV.key))
    
    mq_events = subset(mq_events, Qual >= qual_filter)
    
    events_after = length(mq_events$SV.key)
    alleles_after = length(unique(mq_events$SV.key))
    events_removed = events_before - events_after
    alleles_removed = alleles_before - alleles_after
    cat(" -->", events_removed, "SVs removed at", alleles_removed, "sites.\n")
    # Quality filter
    
    cat(" -> For macaque DELS, DHFFC <", dhffc_filter, "\n")
    events_before = length(mq_events$SV.key)
    alleles_before = length(unique(mq_events$SV.key))
    
    mq_events = subset(mq_events, Type=="<DUP>" | (DHFFC < dhffc_filter & Type=="<DEL>"))
    
    events_after = length(mq_events$SV.key)
    alleles_after = length(unique(mq_events$SV.key))
    events_removed = events_before - events_after
    alleles_removed = alleles_before - alleles_after
    cat(" -->", events_removed, "SVs removed at", alleles_removed, "sites.\n")
    # DELS
    
    cat(" -> For macaque Dups, DHBFC >", dhbfc_filter, "\n")
    events_before = length(mq_events$SV.key)
    alleles_before = length(unique(mq_events$SV.key))
    
    mq_events = subset(mq_events, Type=="<DEL>" | (DHBFC > dhbfc_filter & Type=="<DUP>"))
    
    events_after = length(mq_events$SV.key)
    alleles_after = length(unique(mq_events$SV.key))
    events_removed = events_before - events_after
    alleles_removed = alleles_before - alleles_after
    cat(" -->", events_removed, "SVs removed at", alleles_removed, "sites.\n")
  }
  # DUPS
  # Read depth filters
  # Filter macaque data
  ##########
  
  return(list(mq_events, hu_events))
}