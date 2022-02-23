#' Function to import bio-kine ascii files (.bka)
#'
#' This function allows you to import bio-kine ascii files (.bka). 
#' These are the standard output of BioLogic stopped flow devices.
#' @importFrom magrittr %>%
#' @param path Path to the input .bka file
#' @export
#' @examples
#' read_biokine(path = "my_file.bka")

read_biokine <- function(path = ""){
  comment_lines <- readr::read_lines(path, n_max = 100) %>%
    stringr::str_detect("\"")
  
  comment_lines <- max(which(comment_lines))
  
  data <- readr::read_tsv(file = path, skip = comment_lines, col_names = c("time", "value"), col_types = cols(time = col_double(), value = col_double()))
  
  return(data)
}



#' Function to import all bio-kine ascii files (.bka) in a directory
#'
#' This function allows you to import multiple bio-kine ascii files (.bka) in a directory. 
#' These are the standard output of BioLogic stopped flow devices.
#' @importFrom magrittr %>%
#' @param path Path to the input directory containing the .bka files to read
#' @export
#' @examples
#' read_biokine_all(path = "my_dir")

read_biokine_all <- function(path = "."){
  files <- list.files(path, pattern = ".bka$")
  
  data_list <- lapply(files, read_biokine)
  
  names(data_list) <- files
  
  data <- tibble::enframe(data_list) %>%
    tidyr::unnest(value)
  
  return(data)
}
