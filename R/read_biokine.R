#' Function to import bio-kine ascii files (.bka)
#'
#' This function allows you to import bio-kine ascii files (.bka). 
#' These are the standard output of BioLogic stopped flow devices.
#' @param string String to count letters in. Defaults to "foo"
#' @export
#' @examples
#' read_biokin()

read_biokine <- function(string = "foo"){
  return(nchar(string))
}
