## trying out diff correlation metrics to find genes with opposite exp to Igkc

### Spearman's Correlation ###

  ## compute spearman's correlation bw Igkc and gene j 
  ## store as numeric vector of cor values for all 22,607 genes 
  igkc.spear <- sapply(1:ncol(gexp), function(j) {
    cor(gexp[, i], gexp[, j], method = "spearman")
  })
  ## match cos sim vals with respective names of gene j 
  names(igkc.spear) <- colnames(gexp) 
  ## check that cos sim val is 1 when computed with Igkc itself 
  igkc.spear['Igkc']
  ## check how many negative cor vals we get
  sum(igkc.spear < 0) ## 3930 genes have neg spearman cor vals 
  min(igkc.spear)
  which(igkc.spear == min(igkc.spear))
  ## plot histogram of spearman cor vals 
  hist(igkc.spear, xlab = 'Spearman Correlation, Rho', ylab = 'Number of Genes', 
       main = 'Distribution of Spearman with Igkc Gene')
  
  ### Lowest Spearman Correlation ### 
  ## order cosim values in decreasing order 
  spvals <- order(igkc.spear)
  ## get the top 20 lowest spearman corr values 
  low_spvals <- igkc.spear[spvals[1:20]]
  ## Get names of the 20 lowest spearman corr genes 
  names(low_spvals)
  
### Kendall's Correlation ###
  
  ## compute kendall's correlation bw Igkc and gene j 
  ## store as numeric vector of cor values for all 22,607 genes 
  igkc.kend <- sapply(1:ncol(gexp), function(j) {
    cor(gexp[, i], gexp[, j], method = "kendall")
  })
  ## match cos sim vals with respective names of gene j 
  names(igkc.kend) <- colnames(gexp) 
  ## check that cos sim val is 1 when computed with Igkc itself 
  igkc.kend['Igkc']
  ## check how many negative cor vals we get
  sum(igkc.kend < 0) ## 3930 genes have neg spearman cor vals 
  min(igkc.kend)
  which(igkc.kend == min(igkc.kend))
  ## plot histogram of kendall cor vals 
  hist(igkc.kend, xlab = 'Kendall Correlation', ylab = 'Number of Genes', 
       main = 'Distribution of Kendall Correlation with Igkc Gene')
  
  ### Lowest Kendall Correlation ### 
  ## order cosim values in decreasing order 
  kdvals <- order(igkc.kend)
  ## get the top 20 lowest spearman corr values 
  low_kdvals <- igkc.kend[kdvals[1:20]]
  ## Get names of the 20 lowest spearman corr genes 
  names(low_kdvals)
  
  
### Other
  ggplot(data, aes(x = x, y = y, col=log10(Gm11549+1))) + 
    geom_point(size = 0.1) + scale_color_gradient(low = 'lightgrey', high='red') + 
    facet_wrap(~dataset, scales="free", ncol = 2) + theme_classic()
  