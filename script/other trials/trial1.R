### Computing t-test for Igkc gene across ORlow and ORhigh groups
  
## Get 2 Igkc gene expression vectors for ORhigh and ORlow samples 
  ## in gexp df, first 7858 spots correspond to ORhigh 
  ## and the following 8226 spots correspond to ORlow

  ## get Igkc gene expression vector across all 7858 ORhigh spots 
  igkcORhigh <- gexp[1:7858, 'Igkc']
  
  ## get Igkc gene expression vector across all 8226 ORlow spots
  igkcORlow <- gexp[7859:nrow(gexp), 'Igkc']
  
## Create boxplots to visualize distribution of gene exp for each sample
  boxplot(igkcORhigh, igkcORlow, names = c("Igkc Expression in ORhigh", 
                                           "Igkc Expression in ORlow"))
  
## Perform Welch 2-sample t-test to check whether igkcORlow > igkcORhigh
  ## store results in res 
  res <- t.test(igkcORlow, igkcORhigh, alternative = "g", var.equal=FALSE)
  ## since p-value 2.2e-16, which is less than 0.05, we reject the null 
  ## there is a statistcal difference b/n ORlow and ORhigh for Igkc 
  ## store p-value for Igkc betweeen 2 OR groups as igkc_pval variable 
  igkc_pval <- res$p.value
  
  ## use igkc_pval for comparison for future genes 
  
### Compute t-test for all genes in data ##
  ## store p-vals as 
  gexp_pvals <- sapply(1:10, function(j) {
    geneORhigh <- gexp[1:7858, j]
    geneORlow <- gexp[7859:nrow(gexp), j]
    t.test(geneORlow, geneORhigh, alternative = "g", var.equal=FALSE)
    
  })
  ## match each 
colnames(gexp_pvals) <- colnames(gexp[1:10])
gexp_pvals$p.value

  
  
  