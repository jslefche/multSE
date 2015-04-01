# Function to conduct bootstrap and permutation tests
multSE = function(D, nresamp = 10000, group = 1, permanova = F, ... ) {
  
  # Ensure distance matrix is of class==matrix and not class==dist
  D = as.matrix(D)
  
  if(permanova == T) {

    # Conduct permutation (replace = F) and boostrapped (replace = T) resampling
    mult.SE.list = lapply(c(F, T), function(replace) {
      
      # Bootstrap subsetted distance matrix by smallest sample size across all groups
      do.call(cbind, lapply(2:min(table(group)), function(nsub) {
        
        # Define model parameters
        group.resamp = factor(rep(1:length(unique(group)), each = nsub))
        g = length(levels(factor(group.resamp)))
        X = model.matrix( ~ factor(group.resamp))
        H = X %*% solve(t(X) %*% X) %*% t(X)
        
        # And for each level of number of resamples
        do.call(rbind, lapply(1:nresamp, function(iresamp) {
          
          cols = do.call(c, lapply(unique(group), function(igroup) 
            # Randomly sample columns numbers for each group
            sample(which(group == igroup), size = nsub, replace = replace) ) )
          
          # Sample distance matrix from column numbers
          D.samp = D[cols, cols]
          
          # Calculate multivariate SE based on residuals from PERMANOVA
          N = dim(D.samp)[1]
          I = diag(N)
          A = -0.5 * D.samp^2 
          G = A - apply(A, 1, mean) %o% rep(1, N) - rep(1, N) %o% apply(A, 2, mean) + mean(A) 
          MSE = sum(diag((I - H) %*% G))/(N - g)   
          sqrt(MSE/nsub)
          
        } ) )
      } ) )
    } )
    
    #Calculate means and quantiles
    means = colMeans(mult.SE.list[[1]])
    means.p = colMeans(mult.SE.list[[2]])
    lower.ci = apply(mult.SE.list[[2]], 2, function(x) quantile(x, prob = 0.025))
    upper.ci = apply(mult.SE.list[[2]], 2, function(x) quantile(x, prob = 0.975))
    
    # Return data.frame with means and quantiles
    df = data.frame(
      n.samp = 1:length(means),
      means =  means,
      lower.ci = lower.ci + (means - means.p),
      upper.ci = upper.ci + (means - means.p) )
    
  } else {
    
    # Conduct bootstrapping for each group
    df = do.call(rbind, lapply(unique(group), function(igroup) {
      
      # Subset distance matrix by group
      subset.D = D[group == igroup, group == igroup]
      
      if(sum(subset.D) == 0) {
        data.frame(
          group = igroup,
          n.samp = 1,
          means =  0,
          lower.ci = 0,
          upper.ci = 0) 
        
        } else {
        
        # Conduct permutation (replace = F) and boostrapped (replace = T) resampling
        mult.SE.list = lapply(c(F, T), function(replace) {
          
          # Bootstrap subsetted distance matrix by each sample size
          do.call(cbind, lapply(2:ncol(subset.D), function(nsub) {
            # And for each level of number of resamples
            do.call(rbind, lapply(1:nresamp, function(iresamp) {
              
              # Randomly sample distance matrix
              D.samp = subset.D[sample(1:ncol(subset.D), size = nsub, replace), sample(1:ncol(subset.D), size = nsub, replace)]
              
              # Calculate multivariate SE based on SS
              n = dim(as.matrix(D.samp))[1]
              ss = sum(D.samp^2)/n
              v = ss/(n-1)
              sqrt(v/n)
    
            } ) ) 
          } ) )
        } )
        
        #Calculate means and quantiles
        means = colMeans(mult.SE.list[[1]])
        means.p = colMeans(mult.SE.list[[2]])
        lower.ci = apply(mult.SE.list[[2]], 2, function(x) quantile(x, prob = 0.025))
        upper.ci = apply(mult.SE.list[[2]], 2, function(x) quantile(x, prob = 0.975))
        
        # Return data.frame with means and quantiles
        data.frame(
          group = igroup,
          n.samp = 1:length(means),
          means =  means,
          lower.ci = lower.ci + (means - means.p),
          upper.ci = upper.ci + (means - means.p) )
        
        }
        
      } ) )
  }
    
  return(df) 

  }