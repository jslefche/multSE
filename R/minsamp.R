minsamp = function(output, ignore.min = T) {
  
  # Remove all samples N = 2 from dataset before calculating minimum
  if(ignore.min == TRUE) output = subset(output, n.samp > 2)
  
  # Retrieve minimum number of samples for which error bars still overlap
  if(nrow(output) == 0) 
    
    data.frame(
      min.mean = NA, 
      min.lower.ci = NA, 
      min.upper.ci = NA, 
      min.n = NA
      ) else 
  
      data.frame(
        min.mean = min(output$means),
        min.lower.ci = output[which.min(output$means), "lower.ci"],
        min.upper.ci = output[which.min(output$means), "upper.ci"],
        min.n = min(output[output$lower.ci <= output[which.min(output$means), "upper.ci"], "n.samp"])
      )
      
}