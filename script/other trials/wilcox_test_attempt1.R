### Computing Wilcox Test for Igkc gene across ORlow and ORhigh groups

### DATA ###
    ## set wd 
    setwd('/Users/seeunkwon/Desktop/Research/JEFworks Lab/GBM_Mouse_Model')
    ## load relevant packages
    library(tidyverse)
    install.packages("coin")
    library(coin)
    ## load data
    load(file="stdeconvolve_resultsall_sansmito_n13_for_iris.RData")
    ## create new matrix with only spots and gene expression 
    gexp <- data[, 4:ncol(data)]
    ## Notice that in gexp df, first 7858 spots correspond to ORhigh 
    ## and the following 8227 spots correspond to ORlow 


### 2-sample t-test on Igkc Gene Expression for ORlow vs. ORhigh ###
    ## get total Igkc gexp across all 8227 ORlow spots
    igkcORlow <- gexp[7859:nrow(gexp), 'Igkc']
    ## get total Igkc gexp across all 7858 ORhigh spots 
    igkcORhigh <- gexp[1:7858, 'Igkc']
    
    v1 <- gexp[, 'Igkc']
    v2 <- gexp[, 'Ttr']

    ## Perform Wilcox Rank-Sum Test to check whether igkcORlow > igkcORhigh
    res <- wilcox.test(igkcORlow, igkcORhigh, alternative = "g", 
                       na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=T)
    igkc_pval <- res$p.value
    ## use igkc_pval for comparison for future genes 
    
    for (j in 1:1000) {
      vec1 <- gexp[1:7858, j] ## gexp vector across ORhigh spots
      vec2 <- gexp[7859:nrow(gexp), j] ## gexp vector across ORlow spots
      
      ## compute t-test to see whether ORhigh > ORlow for opposite 
      tres <- wilcox.test(vec1, vec2, alternative = "g", 
                          na.rm=TRUE, paired=FALSE, exact=FALSE)
      g_pval <- tres$p.value
      print(g_pval)
    }
    
    vec3 <- gexp[1:7858, 4]
    vec4 <- gexp[7859:nrow(gexp), 4]
    tres <- wilcox.test(vec3, vec4, alternative = "g", 
                        na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=T)


### 2-sample t-test on all genes in data ###
## Goal: to find genes that are higher in expression in ORhigh than ORlow
## create an empty df for p-values
pval_df <- data.frame(matrix(ncol=0, nrow=1))  
## set row name of df as p-value 
rownames(pval_df)[1] <- 'p value'
## initialize variable to be updated 
i <- 1 
## set alpha value for significance 
alpha <- 0.05 # standard alpha value at 95% confididence level
alpha_mod <- 0.05/1000 # corrected alpha value from Bonferroni 
## for loop for going through all genes in data and compute t-test 
for (j in 1:1000) {
  vec1 <- gexp[1:7858, j] ## gexp vector across ORhigh spots
  vec2 <- gexp[7859:nrow(gexp), j] ## gexp vector across ORlow spots
  
  ## compute t-test to see whether ORhigh > ORlow for opposite 
  tres <- wilcox.test(vec1, vec2, alternative = "g", 
                      na.rm=TRUE, paired=FALSE, exact=FALSE)
  ## get p-val of the gene from t-test 
  g_pval <- tres$p.value
  
  ## add p-val only if null is rejected ORhigh > ORlow
  if ((g_pval <= alpha_mod) && (!is.na(g_pval))) {
    pval_df[, i] <- g_pval
    names(pval_df)[i] <- colnames(gexp[j])
    i <- i + 1 ## increment i after adding to list 
  }
}
## note that we have now 11 genes that are differentially expressed in ORhigh vs ORlow 
## with a higher expression in ORhigh than ORlow 

## rearrange pvalues in increasing order (most to least statistically diff)
pval_df <- sort(pval_df, decreasing = FALSE)

## try plotting the 11 genes 
ggplot(data, aes(x = x, y = y, col=log10(Cdk5r2 + 1))) + 
  geom_point(size = 0.1) + scale_color_gradient(low = 'lightgrey', high='red') + 
  facet_wrap(~dataset, scales="free", ncol = 2) + theme_classic() 


## stdeconvolve tutorial 
ggplot(data, aes(x = x, y = y, col=log10(C4b))) + 
  geom_point(size = 0.1) + scale_color_gradient(low = 'lightgrey', high='red') + 
  facet_wrap(~dataset, scales="free", ncol = 2) + theme_classic() 
