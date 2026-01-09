

#' Sample from a 3D simplex using Dirichlet distribution
sample_simplex_dirichlet <- function(n, alpha) {
  stopifnot(length(alpha) == 3)
  x <- matrix(rgamma(n * 3, shape = alpha), ncol = 3, byrow = TRUE)
  x / rowSums(x)
}

read_simplex_input<-functions(file){
  #' Read Simplex Input File
  #'
  #' This function reads a Simplex input file and returns its contents as a character vector.
  #'
  #' @param file A string representing the path to the Simplex input file.
  #'
  #' @return A character vector containing the lines of the Simplex input file.
  #' @export
  #'
  #' @examples
  #' \dontrun{
  #' input_lines <- read_simplex_input("path/to/simplex_input.txt")
  #' }
  
  if(!file.exists(file)){
    stop("The specified file does not exist.")
  }

  lines <- read.csv(file, header = FALSE, stringsAsFactors = FALSE)


  return(lines)
}   