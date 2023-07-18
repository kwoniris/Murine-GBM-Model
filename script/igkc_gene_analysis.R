#########################################################################
########################## Igkc Gene Analysis  ##########################
#########################################################################
### Goal: Finding genes with similar vs opposite gene expression to Igkc 
### Data: 4 sample murine GBM model dataset 

### PACKAGES ### 
    ## load the tidyverse package, which includes ggplot2 for visualization
    library(tidyverse)

### DATA ###
    ## set working directory 
    setwd("~/Desktop/Projects/Murine-GBM-Model")
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
    
    ## create new matrix with only spots and gene expression 
    gexp <- data[, 4:ncol(data)]
    
    ## filter out genes and spots with negligible expression 
    genes.to.keep <- which(colSums(gexp) > 0) 
    spots.to.keep <- which(rowSums(gexp) > 0)
    gexp <- gexp[spots.to.keep, genes.to.keep]
    
    ## check dimensions of new matrix to be 16,804 spots by 22,607 genes 
    dim(gexp) 
    
    ## divide data into 4 OR samples
    v1 <- data[, 1] == 'ORlow5'
    gexp.ORlow5 = gexp[v1, ] ## 4323 spots of 22607 genes 
    v2 <- data[, 1] == 'ORlow3'
    gexp.ORlow3 = gexp[v2, ] ## 3904 spots of 22607 genes
    v3 <- data[, 1] == 'ORhigh2'
    gexp.ORhigh2 = gexp[v3, ] ## 4052 spots of 22607 genes 
    v4 <- data[, 1] == 'ORhigh4'
    gexp.ORhigh4 = gexp[v4, ] ## 3806 spots of 22607 genes 

### IGKC Gene ###
    ## getting the gene expression for Igkc gene 
    ## set our gene of interest to Igkc 
    i = 'Igkc'
    ## plot normalized Igkc expression with log10 
    ggplot(data, aes(x = x, y = y, col=log10(Igkc+1))) + 
      geom_point(size = 0.1) + scale_color_gradient(low = 'lightgrey', high='red') + 
      facet_wrap(~dataset, scales="free", ncol = 2) + theme_classic()
    
### Metric #1: Cosine Similarity ###
    ## compute cosine similarity bw Igkc and gene j 
    ## store as numeric vector of cos sim values for all 22,607 genes 
    igkc.cosim <- sapply(1:ncol(gexp), function(j) {
      lsa::cosine(gexp[, i], gexp[, j])
    })
    ## match cos sim vals with respective names of gene j 
    names(igkc.cosim) <- colnames(gexp) 
    ## check that cos sim val is 1 when computed with Igkc itself 
    igkc.cosim['Igkc']
    ## notice that there aren't any negative cos sim vals 
    which(igkc.cosim < 0)

    ## plot histogram of cos sim vals 
    hist(igkc.cosim, xlab = 'Cosine Similarity', ylab = 'Number of Genes', 
         main = 'Distribution of Cosine Similarity with Igkc Gene')
    ## get 95th percentile of cosim values 
    quantile(igkc.cosim, .95) ## too small, 0.0794 
    ## order cosim values in decreasing order 
    cosimvals <- order(-igkc.cosim)
    ## get the top 10 highest cosim values 
    high_cosimvals <- igkc.cosim[cosimvals[2:11]]
    
    ## Top 10 genes with highest cosim vals to Igkc 
    ## Get names of the 10 similar cosim genes 
    names(high_cosimvals)

    ## plot the expression of each gene vs. Igkc on x-y plot
    ## plot the expression of each gene across 4 samples 
    for (x in 1:10) {
      # check which gene 
      name <- names(high_cosimvals[x])
      filename <- paste(name, '_ggplot.png')
      dir <- '/Users/seeunkwon/Desktop/Research/JEFworks Lab/GBM_Mouse_Model/Igkc Similar'
      # create new png file 
      png(file = paste('Igkc_vs_', name, '.png'))
      # plot expression for each gene vs Igkc
      plot(x=gexp[, i], y=gexp[, name], xlab ='Igkc Gene Expression', 
           ylab = paste(name, ' Gene Expression'))
      # store plots as new file 
      dev.off()
      
      # create ggplot for gene expression of that gene across 4 samples
      ggplot(data[geneexp[,name]>0,], aes(x = x, y = y)) + geom_point(size = 0.1) + 
        facet_wrap(~dataset, scales="free", ncol = 2) + theme_classic()
      # save ggplot as new png file for each gene 
      ggsave(filename, plot = last_plot(), device = "png", path = dir)
    }

### Metric #1.1: Cosine Similarity Filtered ### 
    ## compute cosine similarity after filtering out orthogonal cos sim vals 
    
    # create an empty dataframe
    cos.sim.gexp = data.frame(matrix(NA, nrow = nrow(gexp), ncol = 0))
    cos.sim.tot = c()
    k = 0
    
    for(j in 1:ncol(gexp)) { # for-loop over columns of genes 
      x = names(gexp)[j]
      ## create sub-df of gene X and Igkc
      gexp.sub <- gexp[, c(i, x)]
      ## remove rows with any zeros
      gexp.sub <- gexp.sub[apply(gexp.sub!=0, 1, all),]
      ## compute cos similarity bw Igkc and gene X 
      cos.sim.val <- lsa::cosine(gexp.sub[, i], gexp.sub[, x])
      ## filter out NA values and zero (orthogonal) cos sim vals 
      if ((cos.sim.val != 0) & (!is.na(cos.sim.val))) {
        cos.sim.tot <- append(cos.sim.tot, cos.sim.val) ## list of cos sim values
        cos.sim.gexp[, (k+1)] <- gexp[, x]
        names(cos.sim.gexp)[k+1] <- names(gexp)[j]
        names(cos.sim.tot)[k] <- names(gexp)[j]
        k = k+1 ## increment k each time a new val is added to cos.sim.tot
      }
    }
    ## cos.sim.tot represents the list of cosine similarity values 
    ## cos.sim.gexp is the new gene set with nonzero cos sim vs. Igkc 
    
    ## check dimensions of the total cos sim dataframe gene set 
    dim(cos.sim.gexp) ## 20394 genes with cos similarity of nonzero
    
    ## plot modified cosine similarity histogram without orthogonal entries
    hist(cos.sim.tot, xlab = 'Cosine Similarity',
         main = 'Distribution of Cosine Similarity of Igkc Gene',
         cex.main=1.0, cex.lab=1.0, cex.axis=0.75)

### Metric 2: Negative Correlation ###
    ## using negative correlation to find genes with opposite exp to Igkc

    ## create sub gexp matrix with genes with low cos similarity 
    ## low cos sim is set to those less than the mean cos sim 
    cos.sim.low <- which(cos.sim.tot < mean(cos.sim.tot)) 
    ## total of 12817 genes have cos sim lower than the mean value 
    low.gexp <- cos.sim.gexp[, cos.sim.low]
    
    ## create new list for those genes of interest 
    genes.interest = list()
    
    ## create a list of genes with negative correlation 
    for (h in 1:ncol(low.gexp)){
      cor.val = cor(gexp[,'Igkc'], low.gexp[, h], method = 'spearman')
      print(cor.val)
      
      thres.val = -0.035 ## change this threshold value 
      if ((cor.val < thres.val) && (!is.na(cor.val))) {
        genes.interest = append(genes.interest, names(low.gexp)[h])
      }
    }
    ## when the threshold val is -0.035, we have narrowed down to 14 genes 
    ## plot the 14 genes using ggplot2
    ## still, these genes do not show best results that we're looking for. 

### Metric 3: Correlation ###
    ## to find genes with similar correlation vals, just switch to cor() 
    ## instead of lsa::cosine()
    ## here, we used pearson and spearman corr methods 