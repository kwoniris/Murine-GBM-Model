### Computing t-test for Igkc gene across ORlow and ORhigh groups
### Here, we eliminated ORhigh2 from dataset due to experimental limitations.
### @date 07-10-2023

### DATA ###
    ## set wd 
    setwd('/Users/seeunkwon/Desktop/Research/JEFworks Lab/GBM_Mouse_Model')
    ## load relevant packages
    library(tidyverse)

    ## load data
    load(file="stdeconvolve_resultsall_sansmito_n13_for_iris.RData")
  
    ## first 4052 spots correspond to ORhigh2 
    ## remove ORhigh2 from data
    data1 <- data[4053:nrow(data), ]
    
    ## create new gexp matrix with only spots and gene exp 
    gexp1 <- data1[, 4:ncol(data1)]

### 2-sample t-test ### 
    ## IGKC only 
    ## get total Igkc gexp across all 8227 ORlow spots
    igkcORlow <- gexp1[3807:nrow(gexp1), 'Igkc']
    ## get total Igkc gexp across all 3806 ORhigh spots 
    igkcORhigh <- gexp1[1:3806, 'Igkc']
    
    ## perform Welch 2-sample t-test to check whether igkcORlow > igkcORhigh 
    res <- t.test(igkcORlow, igkcORhigh, alternative = "g", var.equal=FALSE, 
                  conf.level=0.95)
    igkc_pval <- res$p.value
    
    ## ALL GENES 
    ## To find genes that are ORlow > ORhigh in expression 
    pvals <- data.frame(matrix(ncol=0, nrow=1))  
    ## set row name of df as p-value 
    rownames(pvals)[1] <- 'p value'
    ## set alpha for 95% confidence level 
    alpha <- 0.05
    ## correct alpha using Bonferroni 
    alpha_mod <- 0.05 / ncol(gexp1)
    ## initialize variable for loop 
    i <- 1
    ## for loop for going through all genes in data and compute t-test 
    for (j in 1:ncol(gexp1)) { 
        vec1 <- gexp1[1:3806, j] ## gexp vector across ORhigh spots
        vec2 <- gexp1[3807:nrow(gexp1), j] ## gexp vector across ORlow spots
        tres <- t.test(vec2, vec1, alternative = "g", var.equal=FALSE, 
                       conf.level=0.95)
        g_pval <- tres$p.value 
        
        ## add p-val only if significant 
        if ((g_pval <= alpha_mod) && (!is.na(g_pval))) {
          pvals[, i] <- g_pval
          names(pvals)[i] <- colnames(gexp1[j])
          i <- i + 1 ## increment i after adding to list 
        }
    }
    ## order pvals from most to least statistically significant
    pvals <- sort(pvals, decreasing = FALSE)
    
### Wilcoxon Rank-Sum Test ### 
    ## IGKC only 
    ## get total Igkc gexp across all 8227 ORlow spots
    igkcORlow <- gexp1[3807:nrow(gexp1), 'Igkc']
    ## get total Igkc gexp across all 3806 ORhigh spots 
    igkcORhigh <- gexp1[1:3806, 'Igkc']
    
    ## perform Welch 2-sample t-test to check whether igkcORlow > igkcORhigh 
    res <- wilcox.test(igkcORlow, igkcORhigh, alternative = "g", 
                       na.rm=TRUE, paired=FALSE, exact=FALSE)
    igkc_pval <- res$p.value

    ## ALL GENES 
    ## To find genes that are ORlow > ORhigh in expression 
    p_wilcox <- data.frame(matrix(ncol=0, nrow=1))  
    ## set row name of df as p-value 
    rownames(p_wilcox)[1] <- 'p value'
    ## initialize variable for loop 
    i <- 1
    ## for loop for going through all genes in data and compute t-test 
    for (j in 1:ncol(gexp1)) { 
      v1 <- gexp1[1:3806, j] ## gexp vector across ORhigh spots
      v2 <- gexp1[3807:nrow(gexp1), j] ## gexp vector across ORlow spots
      tres <- t.test(v2, v1, alternative = "g", var.equal=FALSE, 
                     conf.level=0.95)
      pval_wilcox <- tres$p.value 
      
      ## add p-val only if significant 
      if ((pval_wilcox <= alpha_mod) && (!is.na(pval_wilcox))) {
        p_wilcox[, i] <- pval_wilcox
        names(p_wilcox)[i] <- colnames(gexp1[j])
        i <- i + 1 ## increment i after adding to list 
      }
    }
    ## order pvals from most to least statistically significant
    p_wilcox <- sort(p_wilcox, decreasing = FALSE)
    
### SUBSETTNG DATA FOR TOPIC 2 ONLY ###
    ## get spots that are around the tumor, specified as topic 2 

    ## get deconProp for topic 2 
    deconProp2 <- deconProp[, 2] ## 15973 spots
    ## plot deconProp2 across all samples
    ggplot(data[names(deconProp2), ], aes(x=x, y=y, col=deconProp2)) + 
      geom_point(size = 0.2) + scale_color_gradient(low = 'lightgrey',
                                                    high = 'red') + 
      facet_wrap(~dataset, scales = "free", ncol = 2) + theme_classic()
    
    ## update deconProp2 to those greater than zero 
    deconProp2_nonzero <- deconProp2[deconProp2 > 0] ## 4960 spots total 
    
    ## get spots that are both in deconProp2_nonzero and data1 
    spots_overlap <- intersect(names(deconProp2_nonzero), 
                               rownames(data1)) ## 3846 spots total
    
    ## eliminate ORhigh2 from deconProp2 topic 2 
    deconProp2_curr <- deconProp2_nonzero[spots_overlap]
    
    ## 3846 spots are in topic 2 for ORhigh4, ORlow3, ORlow5
    spots_interest <- names(deconProp2_curr)

    ## subset data for only those spots in topic2 without ORhigh2
    data_t2 <- data1[spots_interest, ]
    gexp1_t2 <- gexp1[spots_interest, ]
    
    ## divide subsetted data for topic 2 based on ORlow and ORhigh
    w1 <- data_t2[, 1] == 'ORlow3'
    w2 <- data_t2[, 1] == 'ORlow5'
    w3 <- data_t2[, 1] == 'ORhigh4'
    ## create gexp for 2822 ORlow spots for 32272 genes in topic 2 
    gexp_t2_ORlow <- gexp1_t2[c(w1,w2), ] 
    ## create gexp for 1024 ORhigh spots for 32272 genes in topic 2
    gexp_t2_ORhigh <- gexp1_t2[w3, ] 
    
### Topic 2: Wilcox and 2-sample t-test ### 
    ## IGKC only 
    ## get total Igkc gexp across all 2822 ORlow spots
    igkc1 <- gexp_t2_ORlow[, 'Igkc']
    ## get total Igkc gexp across all 1024 ORhigh spots 
    igkc2 <- gexp_t2_ORhigh[, 'Igkc']
    
    ## perform Welch 2-sample t-test to check igkcORlow > igkcORhigh 
    res1 <- wilcox.test(igkc1, igkc2, alternative = "g", 
                        na.rm=TRUE, paired=FALSE, exact=FALSE)

    ## ALL GENES 
    ## To find genes that are ORlow > ORhigh in expression 
    pvals <- data.frame(matrix(ncol=0, nrow=1))  
    ## set row name of df as p-value 
    rownames(pvals)[1] <- 'p value'
    ## set alpha for 95% confidence level 
    alpha <- 0.05
    ## correct alpha using Bonferroni 
    alpha_mod <- 0.05 / ncol(gexp1)
    ## initialize variable for loop 
    i <- 1
    ## for loop for going through all genes in data and compute t-test 
    for (j in 1:ncol(gexp1)) { 
      vec1 <- gexp_t2_ORhigh[, j] ## gexp vector across ORhigh spots
      vec2 <- gexp_t2_ORlow[, j]  ## gexp vector across ORlow spots
      res2 <- wilcox.test(vec2, vec1, alternative = "g", 
                          na.rm=TRUE, paired=FALSE, exact=FALSE)
      g_pval <- res2$p.value 
      
      ## add p-val only if significant 
      if ((g_pval <= alpha_mod) && (!is.na(g_pval))) {
        pvals[, i] <- g_pval
        names(pvals)[i] <- colnames(gexp1[j])
        i <- i + 1 ## increment i after adding to list 
      }
    }
    ## order pvals from most to least statistically significant
    pvals <- sort(pvals, decreasing = FALSE)
    
    ## plot gene expression only for the topic 2 regions
    ggplot(data_t2, mapping=aes(x=x, y=y, col = log10(Gm10076+1))) + 
      geom_point(size = 0.1) + 
      scale_color_gradient(low = 'lightgrey', high = 'red') + 
      facet_wrap(~dataset, scales="free", ncol = 3) + theme_classic()
  
### Log Fold Change for Topic 2 ###
    
    ## ORlow gene counts matrix for topic 2
    gexp_t2_ORlow <- gexp1_t2[c(w1,w2), ] 
    ## ORhigh gene counts matrix for topic 2
    gexp_t2_ORhigh <- gexp1_t2[w3, ] 
    
    ## calculate the mean of each gene per OR group
    ## here, we take the average per col to get the avg geneexp across each grp
    low_mean <- apply(gexp_t2_ORlow, 2, mean, na.rm=T)
    hi_mean <- apply(gexp_t2_ORhigh, 2, mean, na.rm=T)
    
    ## calculate fold change for ORhigh / ORlow 
    fc <- hi_mean / low_mean 
    
    ## change inf, na, and NaN values as zero 
    fc[is.na(fc)] <- 0
    fc[is.infinite(fc)] <- 0
    ## order from highest to lowest fold change 
    fc <- sort(fc, decreasing=TRUE)
    ## visualize by plot histogram of fc
    hist(fc)
    
    ggplot(data_t2, mapping=aes(x=x, y=y, col = log10(Chrm3+1))) + 
      geom_point(size = 0.1) + 
      scale_color_gradient(low = 'lightgrey', high = 'red') + 
      facet_wrap(~dataset, scales="free", ncol = 3) + theme_classic()
    
## Try fold change for only those genes in deconGexp to get gen es 
##that are actually expressed and avoid genes that have little to no exp.
    
    ## ORlow gene counts matrix for topic 2, deconGexp genes only 
    foo1 <- gexp_t2_ORlow[, colnames(deconGexp)]
    ## ORhigh gene counts matrix for topic 2, deconGexp genes only 
    foo2 <- gexp_t2_ORhigh[, colnames(deconGexp)]
    
    ## calculate the mean of each gene per OR group
    ## here, we take the average per col to get the avg geneexp across each grp
    avg_low <- apply(foo1, 2, mean, na.rm=T)
    avg_hi <- apply(foo2, 2, mean, na.rm=T)
    
    ## calculate fold change for ORlow / ORhigh
    avg_low / avg_hi 
    ## we can see that Igkc is included as one of the top genes
    
    ## calculate fold change for ORhigh / ORlow
    fc1 <- avg_hi / avg_low
    ## order from highest to lowest fold change 
    fc1 <- sort(fc1, decreasing=TRUE)
    ## plot histogram of fold chagnes
    hist(fc1, main = paste("Histogram of Fold Change for ORhigh / ORlow in Topic 2"))
    summary(fc1)
    
    