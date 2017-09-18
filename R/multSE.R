#' multSE: multivariate dissimilarity-based standard errors
#' 
#' @param D a sample-by-species distance matrix
#' @param group an optional vector of groups
#' @param nressamp number of resamples in bootstrapping
#' @param permanova, whether standard errors should be computed across all groups (TRUE) or individually for each group (FALSE)
#' 
#' @return a data.frame with the group (if specified), the means +/- 95% confidence intervals for each level of N (samples)

multSE <- function(D, group, nresamp = 10000, permanova = TRUE) {
  
  # Ensure distance matrix is of class==matrix and not class==dist
  D <- as.matrix(D)
  
  if(permanova == TRUE & !missing(group)) {
    
    # Conduct permutation (replace = F) and boostrapped (replace = T) resampling
    mult.SE.list <- lapply(c(FALSE, TRUE), function(replace) {
      
      # Remove groups with only a single replicate
      if(min(table(group)) == 1) {
        
        remove <- names(table(group)[table(group) == min(table(group))])
        
        D <- D[!group %in% remove, !group %in% remove]
        
        group <- group[!group %in% remove]
        
        message("Groups with 1 replicate have been removed from the analysis!")
        
      }
      
      # Bootstrap subsetted distance matrix by smallest sample size across all groups
      do.call(cbind, lapply(2:min(table(group)), function(nsub) {
        
        # Define model parameters
        group.resamp <- factor(rep(1:length(unique(group)), each = nsub))
        
        g <- length(levels(factor(group.resamp)))
        
        X <- model.matrix( ~ factor(group.resamp))
        
        H <- X %*% solve(t(X) %*% X) %*% t(X)
        
        # And for each level of number of resamples
        do.call(rbind, lapply(1:nresamp, function(iresamp) {
          
          cols <- do.call(c, lapply(unique(group), function(igroup) 
            # Randomly sample columns numbers for each group
            sample(which(group == igroup), size = nsub, replace = replace) ) )
          
          # Sample distance matrix from column numbers
          D.samp <- D[cols, cols]
          
          # Calculate multivariate SE based on residuals from PERMANOVA
          N <- dim(D.samp)[1]
          
          I <- diag(N)
          
          A <- -0.5 * D.samp^2 
          
          G <- A - apply(A, 1, mean) %o% rep(1, N) - rep(1, N) %o% apply(A, 2, mean) + mean(A) 
          
          MSE <- sum(diag((I - H) %*% G))/(N - g)   
          
          sqrt(MSE/nsub)
          
        } ) )
        
      } ) )
      
    } )
    
    #Calculate means and quantiles
    means <- colMeans(mult.SE.list[[1]], na.rm = T)
    
    means.p <- colMeans(mult.SE.list[[2]], na.rm = T)
    
    lower.ci <- apply(mult.SE.list[[2]], 2, function(x) quantile(x, prob = 0.025, na.rm = T))
    
    upper.ci <- apply(mult.SE.list[[2]], 2, function(x) quantile(x, prob = 0.975, na.rm = T))
    
    # Return data.frame with means and quantiles
    df <- data.frame(
      n.samp = 1:length(means),
      means =  means,
      lower.ci = lower.ci + (means - means.p),
      upper.ci = upper.ci + (means - means.p) )
    
  } else {
    
    if(missing(group)) group <- 1
    
    # Conduct bootstrapping for each group
    df <- do.call(rbind, lapply(unique(group), function(igroup) {
      
      # Subset distance matrix by group
      subset.D <- D[group == igroup, group == igroup]
      
      if(sum(subset.D) == 0) {
        data.frame(
          group = as.character(igroup),
          n.samp = 1,
          means =  0,
          lower.ci = NA,
          upper.ci = NA) 
        
      } else {
        
        # Conduct permutation (replace = F) and boostrapped (replace = T) resampling
        mult.SE.list <- lapply(c(F, T), function(replace) {
          
          # Bootstrap subsetted distance matrix by each sample size
          do.call(cbind, lapply(2:ncol(subset.D), function(nsub) {
            # And for each level of number of resamples
            do.call(rbind, lapply(1:nresamp, function(iresamp) {
              
              # Randomly sample distance matrix
              D.samp <- subset.D[sample(1:ncol(subset.D), size = nsub, replace), sample(1:ncol(subset.D), size = nsub, replace)]
              
              # Calculate multivariate SE based on SS
              n <- dim(as.matrix(D.samp))[1]
              
              ss <- sum(D.samp^2)/n
              
              v <- ss/(n-1)
              
              sqrt(v/n)
              
            } ) ) 
          } ) )
        } )
        
        #Calculate means and quantiles
        means <- colMeans(mult.SE.list[[1]], na.rm = T)
        
        means.p <- colMeans(mult.SE.list[[2]], na.rm = T)
        
        lower.ci <- apply(mult.SE.list[[2]], 2, function(x) quantile(x, prob = 0.025, na.rm = T))
        
        upper.ci <- apply(mult.SE.list[[2]], 2, function(x) quantile(x, prob = 0.975, na.rm = T))
        
        # Return data.frame with means and quantiles
        data.frame(
          group = as.character(igroup),
          n.samp = 2:(length(means) + 1),
          means =  means,
          lower.ci = lower.ci + (means - means.p),
          upper.ci = upper.ci + (means - means.p) )
        
      }
      
    } ) )
    
  }
  
  return(df) 
  
}