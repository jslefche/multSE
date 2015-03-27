# Function to calculate multivariate standard errors of sum-of-squared distances from centroid of multivariate space   
mult.SE = function (d) {
  n = dim(as.matrix(d))[1]
  ss = sum(d^2)/n
  v = ss/(n-1)
  x = sqrt(v/n)
  return(x)
}

# Function to conduct bootstrap and permutation tests
mult.SE.group = function(d, group, nresamp = 10000, ...) {
  
  # Ensure distance matrix is of class==matrix and not class==dist
  D = as.matrix(D)
  
  # Conduct bootstrapping for each group
  do.call(rbind, lapply(unique(group), function(igroup) {
    
    # Subset distance matrix by group
    subset.D = D[group == igroup, group == igroup]
    
    # Conduct boostrapped (replace = T) or permutation-based (replace = F) resampling
    mult.SE.list = lapply(c(T, F), function(replace) {
      
      # Bootstrap subsetted distance matrix by each sample size
      do.call(cbind, lapply(2:ncol(subset.D), function(nsub) {
        # And for each level of number of resamples
        do.call(rbind, lapply(1:nresamp, function(iresamp) {
          
          # Randomly sample distance matrix
          D.samp = subset.D[sample(1:ncol(subset.D), size = nsub, replace), sample(1:ncol(subset.D), size = nsub, replace)]
          # Calculate multivariate SE
          mult.SE(D.samp)
          
        } ) )
      } ) )
    } )
    
    # Calculate means and quantiles
    means = colMeans(mult.SE.list[[2]])
    means.p = colMeans(mult.SE.list[[1]])
    upper.ci = apply(mult.SE.list[[1]], 2, function(x) quantile(x, prob = 0.975, na.rm = T) )
    lower.ci = apply(mult.SE.list[[1]], 2, function(x) quantile(x, prob = 0.025, na.rm = T) )

    # Return data.frame
    data.frame(
      group = igroup,
      n.samp = 1:length(means),
      means = means,
      bias.lower = lower.ci + (means - means.p),
      bias.upper = upper.ci + (means - means.p) )
    
  } ) )

}