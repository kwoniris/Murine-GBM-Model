# Function to evaluate similarity between two numeric vectors from a dataset 
## evaluate_similiarity function takes a dataset and determines 
## number of vectors in the dataset 
evaluate_similarity <- function(dataset) {
  num_vectors <- nrow(dataset)
  vector_length <- length(dataset[1, 4:ncol(data)])
  ## from 4 as genes start from 4th column 
  
  ## sim matrix initialized using matrix function 
  similarity_matrix <- matrix(0, nrow = num_vectors, ncol = num_vectors)
  ## will store the similarity scores between each pair of vectors 
  
  ## use nested loops??
  ## compare the corresponding elements using another loop 
  ## iterate over each pair of vectors and calcualte sim score 
  for (i in 1:num_vectors) {
    for (j in 1:num_vectors) {
      ## first, initialize the sim score to zero
      similarity_score <- 0 
      
      for (k in 1:vector_length) {
        similarity_score <- similarity_score + abs(dataset[i, k] - dataset[j, k])
        
      }
      
      similarity_matrix[i, j] <- similarity_score
    }
  }
  
  return(similarity_matrix)
}

## Example Run: 
## if we let vectors as the dataset of all 32,272 geneexp vectors 
## similarity_matrix <- evaluate_similarity(vectors)
## print(similarity_matrix)

## Cons: 
## sim matrix for us is going to be 32,272 by 32,272 
## may not be practical to print or analyze entirely!! 

