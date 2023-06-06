# this script describe a function scoring a gene signature. The score is normalised using control genes.
# the code come from a seurat function. View https://www.waltermuskovic.com/2021/04/15/seurat-s-addmodulescore-function/ 
# for more information.

# dds consist in a DESeqDataSet and gene_signature a vector of genes
GetSignatureScore<- function(dds, gene_signature) {
  
  Norm_exprMat <- counts(dds, normalized=T)

  data.avg <- Matrix::rowMeans(x = Norm_exprMat)
  
  # Order genes from lowest to highest mean by row
  data.avg <- data.avg[order(data.avg)]
  
  # the value of nbin is 25 by default
  nbin <- 25
  # create 25 groups of genes. Allows to associateas control genes, genes that are specific
  data.cut <- ggplot2::cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30,
                                  n = nbin,
                                  labels = FALSE,
                                  right = FALSE)
  
  
  # Set the names of the cuts as the gene names
  names(x = data.cut) <- names(x = data.avg)
  data.cut # the genes are order by their average expression
  
  
  # Create an empty list the same length as the number of input gene sets. This will contain the names of the control genes
  cluster.length <- 1 # has to always be one
  ctrl.use <- vector(mode = "list", length = cluster.length) # cluster.length = nbre of gene of interest
  
  
  features <- list(gene_signature)
  ctrl <- 100 # default parameter
  
  # For each of the input gene lists:
  for (i in 1:cluster.length) {
    # Get the gene names from the input gene set as a character vector  
    features.use <- features[[i]]
    
    # Loop through the provided genes (1:num_genes) and for each gene, find ctrl (default=100) genes from the same expression bin (by looking in data.cut):
    for (j in 1:length(x = features.use)) {
      # Within this loop, 'data.cut[features.use[j]]' gives us the expression bin number. We then sample `ctrl` genes from that bin without replacement and add the gene names to ctrl.use.
      ctrl.use[[i]] <- c(ctrl.use[[i]],
                         names(x = sample(x = data.cut[which(x = data.cut == data.cut[features.use[j]])],
                                          size = ctrl, # sélectionne des gènes faisant partie du même bin que celui
                                          replace = FALSE)))
    }
  }
  
  # Test wether the control genes were selected properly
  if(length(ctrl.use[[1]]) != length(gene_signature)*100){
    
    stop('Error : something wrong happended during the selection of control genes.')
    
  }
  
  # all the control genes are place in a matrix
  ctrl.scores <- matrix(data = numeric(length = 1L),
                        nrow = length(x = ctrl.use),
                        ncol = nrow(colData(dds))) # correspond to the number of sample
  
  assay.data <- assay(rlogQ3)
  
  # calculate the average expression of the control genes for each sample
  for (i in 1:length(ctrl.use)) {
    
    # Get control gene names as a vector  
    features.use <- ctrl.use[[i]]
    # For each cell, calculate the mean expression of *all* of the control genes 
    ctrl.scores[i,] <- Matrix::colMeans(x = assay.data[features.use,])
  }
  
  
  
  features.scores <- matrix(data = numeric(length = 1L),
                            nrow = cluster.length,
                            ncol = nrow(colData(dds)))
  
  # Loop through input gene sets and calculate the mean expression of these genes for each cell
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    data.use <- assay.data[features.use, , drop = FALSE]
    features.scores[i, ] <- Matrix::colMeans(x = data.use) # cette formulation est stupide car je calcule la moyenne d'un seul élément (mais ça marche je précise)
  }
  
  features.scores.use <- features.scores - ctrl.scores
  
  return(features.scores.use)
  
}

# the code in the website describe a way to display control genes and their average expression. Since it is a function
# I did not include it.