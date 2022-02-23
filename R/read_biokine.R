#' Function to import bio-kine ascii files (.bka)
#'
#' This function allows you to import bio-kine ascii files (.bka). 
#' These are the standard output of BioLogic stopped flow devices.
#' @param path Path to the input .bka file
#' @export
#' @examples
#' read_biokine(path = "my_file.bka")

read_biokine <- function(path = ""){
  comment_lines <- readr::read_lines(path, n_max = 100) %>%
    stringr::str_detect("\"")
  
  comment_lines <- max(which(comment_lines))
  
  data <- readr::read_tsv(file = path, skip = comment_lines, col_names = c("time", "value"))
  
  return(data)
}
