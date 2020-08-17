# Analysis of expression data

# function to change log2 logits to raw proportions

antilogit2 <- function(x) {
  # x = log2(p/(1-p)) 
  # so p = 2^x / (1 + 2^x)

  return(p = 2^x / (1+ 2^x))
}