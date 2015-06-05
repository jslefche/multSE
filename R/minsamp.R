minsamp = function(output, group = NA, ignore.min = NULL) {
  
  # Remove all samples N = 2 from dataset before calculating minimum
  if(!is.null(ignore.min)) output = subset(output, n.samp > ignore.min)
  
  # Break apart by group
  do.call(rbind, lapply(unique(group), function(i) {
    
    # Subset by group
    if(!all(is.na(group))) output = subset(output, group == i)
  
    # Retrieve minimum number of samples for which error bars still overlap
    if(nrow(output) == 0) 
      
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
        min.mean = min(output$means),
        min.lower.ci = output[which.min(output$means), "lower.ci"],
        min.upper.ci = output[which.min(output$means), "upper.ci"],
        min.n = min(output[output$lower.ci <= output[which.min(output$means), "upper.ci"], "n.samp"])
        )
    
  } ) )
      
}