#######################################################################################
# 
# Function for calculating multivariate standard error separately for each group,
# with double-resampling to obtain error bars, for increasing numbers of replicates.
# (Note also that sample sizes may vary among groups).

# The input will be a distance matrix among all samples, the grouping vector, and the number of
# re-samples (nresamp).
#
# The routine then calculates the means using the permutation approach, while the
# lower and upper quantiles are obtained using the bootstrapping approach
# including an adjustment for the bias in the bootstrap. 
#
# The output consists of three matrices in a list:
# 1. "means" = the means for each group (rows) for each sample size (columns)
# 2. "upper" the upper 0.975 quantile for each group (rows) for each sample size (columns)
# 3. "lower" the lower 0.0.025 quantile for each group (rows) for each sample size (columns)
# 
#######################################################################################
#
MSEgroup.d = function (D, group, nresamp = 1000, ...) {
  
  # Some necessary preliminary functions:            
  mult.SE = function (d) {
    n = dim(as.matrix(d))[1]
    ss = sum(d^2)/n
    v = ss/(n-1)
    x = sqrt(v/n)
    return(x)
  }
  quant.upper = function(x) quantile(x, prob = 0.975, na.rm = TRUE)
  quant.lower = function(x) quantile(x, prob = 0.025, na.rm = TRUE)
  
  # Getting parameters of the problem
  group = factor(group)
  ng = length(levels(group))
  n.i = table(group) 
  
  # Ensure distance matrix is in the form of a matrix (rather than a "distance" object)
  D = as.matrix(D)
  
  # Setting up the matrices for results
  means <- means.b <- lower <- upper <- matrix(rep(0,ng*max(n.i)),ncol = max(n.i), nrow = ng)
  colnames(means) <- colnames(lower) <- colnames(upper) <- 1:max(n.i)
  rownames(means) <- rownames(lower) <- rownames(upper) <- levels(group)
  
  # Loop - do this for each group
  for (igroup in 1:ng) {
    subset.D = D[ group==levels(group)[igroup] , group==levels(group)[igroup] ]
    
    # One matrix is used to store the values under permutation resampling for each sample size
    multSE.store.p = matrix(rep(0,nresamp*n.i[igroup]), ncol = n.i[igroup], nrow = nresamp)
    # One matrix is used to store the values under bootstrap resampling for each sample size
    multSE.store.b = matrix(rep(0,nresamp*n.i[igroup]), ncol = n.i[igroup], nrow = nresamp)
    
    # Bootstrap loop for each sample size.
    for (nsub in 2:n.i[igroup]) {
      for (iresamp in 1:nresamp) {
        ivec.p = sample(1:n.i[igroup], size = nsub, replace = FALSE)
        ivec.b = sample(1:n.i[igroup], size = nsub, replace = TRUE)
        D.perm = subset.D[ivec.p,ivec.p]
        D.boot = subset.D[ivec.b,ivec.b]
        multSE.store.p[iresamp,nsub] = mult.SE(as.dist(D.perm))
        multSE.store.b[iresamp,nsub] = mult.SE(as.dist(D.boot))
      }
    }
    
    # Means and quantiles
    means[igroup,1:n.i[igroup]] = colMeans(multSE.store.p)
    means.b[igroup,1:n.i[igroup]] = colMeans(multSE.store.b)
    upper[igroup,1:n.i[igroup]] = apply(multSE.store.b,MARGIN = 2,quant.upper)
    lower[igroup,1:n.i[igroup]] = apply(multSE.store.b,MARGIN = 2,quant.lower)
    
  } # end of loop for groups
  
  bias =  means - means.b
  lower = lower + bias
  upper = upper + bias
  return(list(means = means, lower = lower, upper = upper))
  
} # end of MSEgroup.d function