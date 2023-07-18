### Continuing Igkc Opposite Gene Expression Analysis ###

### DATA ###
  ## set wd 
  setwd('/Users/seeunkwon/Desktop/Research/JEFworks Lab/GBM_Mouse_Model')
  ## load data
  load(file="stdeconvolve_resultsall_sansmito_n13_for_iris.RData")
  ## note that the first two columns correspond to the spot positions
  data[1:5,1:5]
  ## create new matrix with only spots and gene expression 
  gexp <- data[, 4:ncol(data)]
  
  ## stdeconvolve 
  pos <- data[, c(2,3)]
  
  
### CODE ###
  
### Getting spots around tumor with topic 2 ### 
  
  ## plot the total expression of all genes per spot 
  ggplot(data, aes(x = x, y = y, col=log10(rowSums(gexp)+1))) + 
    geom_point(size = 0.1) + scale_color_gradient(low = 'lightgrey', high='red') + 
    facet_wrap(~dataset, scales="free", ncol = 2) + theme_classic()
  ## filter out non-tumor regions with those spots with gexp >= 3
  testspots <- log10(rowSums(gexp)+1)>3.5
  ## get names of spots that aren't tumor  
  truespots <- names(which(testspots))
  ## get names of spots that are tumor
  tumorspots <- names(which(!testspots))
  ## plot only non-tumor spots 
  ggplot(data[truespots,], aes(x = x, y = y)) + geom_point(size = 0.1, 
    color = "red") + facet_wrap(~dataset, scales="free", ncol = 2) + 
    theme_classic() + ggtitle("Plot of Gene Expression Across Non-Tumor Spots")
  ## plot only tumor spots 
  ggplot(data[tumorspots,], aes(x = x, y = y)) + geom_point(size = 0.1, 
    color = "red") + facet_wrap(~dataset, scales="free", ncol = 2) + 
    theme_classic() + ggtitle("Plot of Gene Expression Across Tumor Region Only")

  ## get spots that are around the tumor, specified as topic 2 
    ## get deconGexp for topic 2 only as a numeric vector
    ## deconGexp2 shows the putative gene exp across selected genes for topic 2
    deconGexp2 <- deconGexp[2,] 
    ## sort deconGexp in decreasing order of gene exp across topic 2
    deconGexp2_desc <- sort(deconGexp2, decreasing=TRUE) 
    
    ## get deconProp for topic 2 
    deconProp2 <- deconProp[, 2]
    ## select for spots with cell-type proportion > 0 
    topic2_spots <- which(deconProp2 > 0) 
    ## notice that 4960 spots (pixels) have >= cell-type proportion for topic 2
    deconProp2 <- deconProp2[topic2_spots]
    ## plot spots in topic 2 
    ggplot(data[topic2_spots,], aes(x = x, y = y)) + geom_point(size = 0.1, 
      color = "red") + facet_wrap(~dataset, scales="free", ncol = 2) + 
      theme_classic() + ggtitle("Plot of Gene Expression for Topic 2")
    
  ## divide topic 2 spots based on ORhigh and ORlow 
    length(topic2_spots) ## 4960 spots in topic 2
    ## get gexp dataframe for topic 2 only 
    gexp_t2 <- gexp[topic2_spots,] ## 4960 spots of 32272 genes total 
    w1 <- data[topic2_spots, 1] == 'ORlow5'
    w2 <- data[topic2_spots, 1] == 'ORlow3'
    ## create new dataframe for all ORlow spots in topic 2 
    gexp_t2_low <- gexp_t2[c(w1, w2), ] ## 2791 spots of 32272 genes total 
    ## create new dataframe for all ORhigh spots in topic 2
    w3 <- data[topic2_spots, 1] == 'ORhigh2'
    w4 <- data[topic2_spots, 1] == 'ORhigh4'
    gexp_t2_hi <- gexp_t2[c(w3, w4), ] ## 2169 spots of 32272 genes total
    ## remove unneeded variables 
    rm(w1, w2, w3, w4)
    
 
### Perform 2-sample t-test on Topic 2 spots only (IGKC) ###
    
    ### 2-sample t-test on Igkc Gene Expression for ORlow vs. ORhigh ###
    ## get total Igkc gexp across all 2791 ORlow spots in topic 2
    igkct2low <- gexp_t2_low[, 'Igkc']
    ## get total Igkc gexp across all 2169 ORhigh spots in topic 2
    igkct2hi <- gexp_t2_hi[, 'Igkc']
    
    ## Perform Welch 2-sample t-test to check whether ORlow > ORhigh
    rest2 <- t.test(igkct2low, igkct2hi, alternative = "g", var.equal=FALSE)
    ## since p-value 6.315e-05, which is less than 0.05, we reject the null 
    ## there is a statistcal difference b/n ORlow and ORhigh for Igkc 
    ## store p-value for Igkc betweeen 2 OR groups as igkc_pval 
    igkcpvalt2 <- rest2$p.value
    ## use igkcpvalt2 for positive comparison for future tests
    
### 2-sample t-test on ALL genes in topic 2 spots ###
    ## Goal: to find genes that are higher in expression in ORhigh than ORlow
    ## in only topic 2 (spots around the tumor)
    ## create an empty df for p-values
    pval_df_t2 <- data.frame(matrix(ncol=0, nrow=1))  
    ## set row name of df as p-value 
    rownames(pval_df_t2)[1] <- 'p value'
    ## initialize variable to be updated 
    i <- 1 
    ## for loop for going through all genes in data and compute t-test 
    for (j in 1:ncol(gexp)) {
      v1 <- gexp_t2_hi[, j] ## gene exp vector across ORhigh spots topic 2
      
      v2 <- gexp_t2_low[, j] ## gene exp vector across ORlow spots topic 2
      
      ## compute t-test to see whether ORhigh > ORlow for opposite 
      test_res <- t.test(v1, v2, alternative = "g", var.equal=FALSE)  
      ## get p-val of the gene from t-test 
      gene_pval <- test_res$p.value
      
      ## add p-val only if null is rejected ORhigh > ORlow
      if ((gene_pval <= igkcpvalt2) && (!is.na(gene_pval))) {
        pval_df_t2[, i] <- gene_pval
        names(pval_df_t2)[i] <- colnames(gexp[j])
        i <- i + 1 ## increment i after adding to list 
      }
    }
    ## we have 395 genes that have ORhigh > ORlow gene exp in topic 2 
    ## rearrange pvalues in increasing order (most to least statistically diff)
    pval_df_t2 %>% arrange('p value')
    
    ## switch rows and cols of dataframe by taking transpose
    pval_df_t2 <- t(pval_df_t2)
    
    ## get list of pvals from dataframe
    pvals_t2 <- pval_df_t2[,1]
    ## get statistical summary of the pvals list from t-test for topic 2
    summary(pvals_t2)
    ## try plotting the pvals on a histogram 
    hist(pvals_t2)
    
    ## check how many overlap with the 11 genes we have previously found from
    ## computing t-test for all spots, not just topic 2
    pval1list <- names(pval_df_t2) ## ch list of genes from t-test for only topic 2
    pval2list <- names(pval_df) ## ch list of genes from t-test for all spots 
    pval3list <- intersect(pval1list, pval2list) ## find overlap genes 
    ## notice 9 genes overlap 
    
    ## try plotting 
    ggplot(data, aes(x = x, y = y, col=log10(Nmi+1))) + 
      geom_point(size = 0.1) + scale_color_gradient(low = 'lightgrey', high='red') + 
      facet_wrap(~dataset, scales="free", ncol = 2) + theme_classic()
    
    
    
    
### Other ###    
    ## plot normalized expression with log10 of 2 highest expressed genes in Topic2
    ggplot(data, aes(x = x, y = y, col=log10(Rims1+1))) + 
      geom_point(size = 0.1) + scale_color_gradient(low = 'lightgrey', high='red') + 
      facet_wrap(~dataset, scales="free", ncol = 2) + theme_classic()

    
    ## other 
    ggplot(data, aes(x = x, y = y, col=log10(Ly6c2+1))) + 
      geom_point(size = 0.1) + scale_color_gradient(low = 'lightgrey', high='red') + 
      facet_wrap(~dataset, scales="free", ncol = 2) + theme_classic()

    
    
    ggplot(data[names(deconProp2),], aes(x = x, y = y, col=deconProp2)) + 
      geom_point(size = 0.1) + scale_color_gradient(low = 'lightgrey', high='red') + 
      facet_wrap(~dataset, scales="free", ncol = 2) + theme_classic()
    
    ggplot(data[names(deconProp2),], aes(x = x, y = y, col=deconProp2>0.1)) + 
      geom_point(size = 0.1) + 
      facet_wrap(~dataset, scales="free", ncol = 2) + theme_classic()
    
    