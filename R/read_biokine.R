#' Function to import bio-kine ascii files (.bka)
#'
#' This function allows you to import bio-kine ascii files (.bka). 
#' These are the standard output of BioLogic stopped flow devices.
#' @param path String to count letters in. Defaults to "foo"
#' @export
#' @examples
#' read_biokin()

read_biokine <- function(path = ""){
  comment_lines <- read_lines(path, n_max = 100) %>%
    str_detect("\"")
  
  comment_lines <- max(which(comment_lines))
  
  data <- read_tsv(file = path, skip = comment_lines, col_names = c("time", "value"))
  
  return(data)
}
