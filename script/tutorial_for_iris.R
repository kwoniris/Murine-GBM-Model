### Tutorial for Iris ###
## This is an example tutorial from Prof. Fan on how to analyze the dataset using basic R programming skills. 

## load data
load(file="stdeconvolve_resultsall_sansmito_n13_for_iris.RData")

## there are 4 samples (two ORhigh, two ORlow)
table(data$dataset)
## 4052 spots in ORhigh2   
## 3806 spots in ORhigh4    
## 3904 spots in ORlow3    
## 4323 spots in ORlow5

## note that the first two columns correspond to the spot positions
data[1:5,1:5]

## to see what genes are available
ncol(data)
colnames(data[,4:ncol(data)])

## pick a gene 
## note they are mouse genes
g <- 'Ccl5'
g %in% colnames(data)

## plot cells as points
## faceted by the dataset
## color by expression of Slc17a7 on log10 scale
## organize by dataset
library(ggplot2)
ggplot(data, aes(x = x, y = y, col=log10(Ccl5+1))) + 
  geom_point(size = 0.1) + scale_color_gradient(low = 'lightgrey', high='red') + 
  facet_wrap(~dataset, scales="free", ncol = 2) + theme_classic()

## another gene
g <- 'Serpinb2'
g %in% colnames(data)

## pick another gene
ggplot(data, 
       aes(x = x, y = y, col=log10(Serpinb2+1))) + geom_point(size = 0.1) + 
  scale_color_gradient(low = 'lightgrey', high='red') + 
  facet_wrap(~dataset, scales="free", ncol = 2) + theme_classic()

## note in these data, STdeconvolve has also already been performed
## with output as deconGexp and deconProp
## see if you can figure out which deconvolved cell-type corresponds to E cells for example
## and visualize their proportions on the tissue

