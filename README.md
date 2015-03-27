# multSE: Multivariate dissimilarity-based standard error

  Implementation of multivariate dissimilarity-based standard error estimates from:

    Anderson, MJ and J Santanta-Garcon. 2015. "Measures of precision for dissimilarity-based multivariate
    analysis of ecological communities." Ecology Letters: 18(1): 66-73.
    
Version: 0.1 (2015-03-27)

Author: Jon Lefcheck (jslefche@vims.edu)

##Examples

###Load functions
```
# Function from supplements
MSEgroup.d = source("https://github.com/jslefche/multSE/edit/master/R/MSEgroup_d.R")
# New function using vectorization
mult.SE.group = source("https://github.com/jslefche/multSE/edit/master/R/mult_SE_group.R")
```
###Load data
```
pk = read.csv("https://github.com/jslefche/multSE/blob/master/data/PoorKnights.csv")
```
###Run function and plot results
```
# Load vegan library
library(vegan)

# Create species-by-site distance matrix
D = vegdist(pk[,3:49] + 1)

# Run new function
output = mult.SE.group(D, factor(pk$Time), nresamp = 10000)
```
###Plot output
```
# Plot output
library(ggplot2)

ggplot(output, aes(x = n.samp, y = means, group = group)) +
  geom_errorbar(aes(ymax = bias.upper, ymin = bias.lower), width = 0.2)+
  geom_point(aes(shape = group, fill = group), size = 4) + 
  scale_shape_manual(values = c(21, 24:25))+
  scale_fill_manual(values = c("black","grey50","white"))+
  facet_wrap( ~ group) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none", 
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank())
```
###Benchmarks vs. old function
```
system.time(MSEgroup.d(D, factor(pk$Time), nresamp = 10000)) #user: 118.22
system.time(mult.SE.group(D, factor(pk$Time), nresamp = 10000)) #user: 37.3
```
