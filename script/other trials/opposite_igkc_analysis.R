### Computing t-test for Igkc gene across ORlow and ORhigh groups
  
### DATA ###
  ## set wd 
  setwd('/Users/seeunkwon/Desktop/Research/JEFworks Lab/GBM_Mouse_Model')
  ## load relevant packages
  library(ggplot2)
  library(tidyverse)
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

  ## Perform Welch 2-sample t-test to check whether igkcORlow > igkcORhigh
  res <- t.test(igkcORlow, igkcORhigh, alternative = "g", var.equal=FALSE)
    ## since p-value 2.2e-16, which is less than 0.05, we reject the null 
    ## there is a statistcal difference b/n ORlow and ORhigh for Igkc 
  ## store p-value for Igkc betweeen 2 OR groups as igkc_pval 
  igkc_pval <- res$p.value
  ## use igkc_pval for comparison for future genes 
  
### 2-sample t-test on all genes in data ###
  ## Goal: to find genes that are higher in expression in ORhigh than ORlow
  ## create an empty df for p-values
  pval_df <- data.frame(matrix(ncol=0, nrow=1))  
  ## set row name of df as p-value 
  rownames(pval_df)[1] <- 'p value'
  ## initialize variable to be updated 
  i <- 1 
  ## for loop for going through all genes in data and compute t-test 
  for (j in 1:ncol(gexp)) {
    vec1 <- gexp[1:7858, j] ## gexp vector across ORhigh spots
    vec2 <- gexp[7859:nrow(gexp), j] ## gexp vector across ORlow spots
    
    ## compute t-test to see whether ORhigh > ORlow for opposite 
    tres <- t.test(vec1, vec2, alternative = "g", var.equal=FALSE)  
    ## get p-val of the gene from t-test 
    g_pval <- tres$p.value
    
    ## add p-val only if null is rejected ORhigh > ORlow
    if ((g_pval <= igkc_pval) && (!is.na(g_pval))) {
      pval_df[, i] <- g_pval
      names(pval_df)[i] <- colnames(gexp[j])
      i <- i + 1 ## increment i after adding to list 
    }
  }
  ## note that we have now 11 genes that are differentially expressed in ORhigh vs ORlow 
  ## with a higher expression in ORhigh than ORlow 
  
  ## rearrange pvalues in increasing order (most to least statistically diff)
  pval_df %>% arrange('p value')

  ## try plotting the 11 genes 
  ggplot(data, aes(x = x, y = y, col=log10(Rpl37a+1))) + 
    geom_point(size = 0.1) + scale_color_gradient(low = 'lightgrey', high='red') + 
    facet_wrap(~dataset, scales="free", ncol = 2) + theme_classic() 

  
  ## stdeconvolve tutorial 
  ggplot(data, aes(x = x, y = y, col=log10(C4b))) + 
    geom_point(size = 0.1) + scale_color_gradient(low = 'lightgrey', high='red') + 
    facet_wrap(~dataset, scales="free", ncol = 2) + theme_classic() 
  