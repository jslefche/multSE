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
tf = tempfile()
download.file("https://raw.githubusercontent.com/jslefche/multSE/master/data/PoorKnights.csv", tf)
pk = read.csv(tf)
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
  labs(x = "Sample size (n)", y = "Multivariate pseudo SE") +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
```
![multSE plot](https://github.com/jslefche/jslefche.github.io/blob/master/img/multSE_plot.jpeg?raw=true)
###Benchmarks vs. old function

Compare `mult.SE.group` to function included in the supplements of ELE paper `MSEgroup.d`.

```
# Calculate system time for old vs. new function 
benchmarks.df = do.call(rbind, lapply(c(10, 100, 1000, 10000), function(i)
  data.frame(
    nresamp = i,
    fn = c("old", "new"),
    time = c(system.time( MSEgroup.d(D, factor(pk$Time), nresamp = i) )[3],
             system.time(mult.SE.group(D, factor(pk$Time), nresamp = i) )[3] ) )
) )

#And plot
ggplot(benchmarks.df, aes(x = nresamp, y = time, group = time, col = time, shape = time)) +
  geom_point(size = 10) +
  scale_color_manual(values = c("red", "black")) + 
  scale_shape_manual(values = c()) +   theme_bw(base_size = 18) +
  labs(x = "Number of resamples", y = "System time (seconds)") +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
```
![multSE benchmark plot](https://github.com/jslefche/jslefche.github.io/blob/master/img/multSE_benchmark.jpeg?raw=true)

As you can see, the new function 'mult.SE.group` is much faster, particularly when the number of resamples is large.
