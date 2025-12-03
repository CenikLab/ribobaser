# to implement multiplicative imputation
# the method deletes columns /rows , genes /samples which have >80% unobserved or zero values. 
## For Multiplicative to work, samples should be in rows and genes should be in columns
# Ensure genes are rownames instead of a separate column as is the case when importing a csv file 

library(zCompositions)
Imputation <- function(dt1,
                       method = "GBM",
                       output = "p-count",
                       transpose = TRUE
                       ) {
  
  if (!is.numeric(dt1[,1])) {
    gene_names <- dt1[,1]
    dt1 <- dt1[,-1, drop = FALSE]
    rownames(dt1) <- gene_names
  } else {
    gene_names <- rownames(dt1)
  }   
  
  # Optionally transpose input data
  if (transpose) {
    dt1 <- t(dt1)
  }
  
  final_result=  zCompositions::cmultRepl(dt1, method = method, output = output)
  
  # Transpose back if needed
  if (transpose) {
    final_result <- t(final_result)
  }
  return(final_result)
}  

