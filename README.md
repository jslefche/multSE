# multSE: Multivariate dissimilarity-based SE

  Implementation of multivariate dissimilarity-based standard error estimates from:

    Anderson, MJ and J Santanta-Garcon. 2015. "Measures of precision for dissimilarity-based multivariate
    analysis of ecological communities." Ecology Letters: 18(1): 66-73.
    
Version: 0.1 (2015-03-27)

Author: Jon Lefcheck (jslefche@vims.edu)

##Examples

###Load functions
```
library(devtools)
# Function from supplements
MSEgroup.d = source_url("https://raw.githubusercontent.com/jslefche/multSE/master/R/MSEgroup_d.R")[[1]]
# New function using vectorization
mult.SE.group = source_url("https://raw.githubusercontent.com/jslefche/multSE/master/R/mult_SE_group.R")[[1]]
```
###Load data
```
# Poor Knights fish survey data, from supplementary material
pk = read.csv("https://github.com/jslefche/multSE/blob/master/data/PoorKnights.csv")
```
###Calculate multivariate SE and confidence intervals
```
# Load vegan library
library(vegan)

# Create species-by-site distance matrix
D = vegdist(pk[,3:49] + 1)

# Run optimized function to generate multivariate SE for each group
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
![plot](https://github.com/jslefche/jslefche.github.io/blob/master/img/multSE_plot.jpeg?raw=true)
###Benchmarks vs. old function

Can test function provided in supplements to Ecol Letters article, versus new function.

```
library(microbenchmark)

#Run benchmarks
bench = microbenchmark(
  MSEgroup.d(D, factor(pk$Time), nresamp = 10000),
  mult.SE.group(D, factor(pk$Time), nresamp = 10000) )

#And plot
boxplot(bench)
```
