minsamp <- function(out, group = NA, ignore.min = NULL) {
  
  # Remove all samples N = 2 from dataset before calculating minimum
  if(!is.null(ignore.min)) out >- subset(out, n.samp > ignore.min)
  
  # Break apart by group
  do.call(rbind, lapply(unique(group), function(i) {
    
    # Subset by group
    if(!all(is.na(group))) out <- subset(out, group == i)
    
    # Retrieve minimum number of samples for which error bars still overlap
    if(nrow(out) == 0) 
      
      data.frame(
        group = i,
        min.mean = NA, 
        min.lower.ci = NA, 
        min.upper.ci = NA, 
        min.n = NA
      ) 
    
    else 
      
      data.frame(
        group = i,
        min.mean = min(out$means),
        min.lower.ci = out[which.min(out$means), "lower.ci"],
        min.upper.ci = out[which.min(out$means), "upper.ci"],
        min.n = min(out[out$lower.ci <= out[which.min(out$means), "upper.ci"], "n.samp"])
      )
    
  } ) )
  
}
