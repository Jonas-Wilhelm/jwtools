#' Import bio-kine ascii files (.bka)
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
  
  data <- readr::read_tsv(file = path, skip = comment_lines, col_names = c("time", "value"), col_types = readr::cols(time = readr::col_double(), value = readr::col_double()))
  
  return(data)
}






#' Import all bio-kine ascii files (.bka) in a directory
#'
#' This function allows you to import multiple bio-kine ascii files (.bka) in a directory. 
#' These are the standard output of BioLogic stopped flow devices.
#' @importFrom magrittr %>%
#' @param path Path to the input directory containing the .bka files to read
#' @export
#' @examples
#' read_biokine_all(path = "my_dir")

read_biokine_all <- function(path = "."){
  files <- list.files(path, pattern = ".bka$", full.names = T, recursive = T)
  
  data_list <- lapply(files, read_biokine)
  
  names(data_list) <- list.files(path, pattern = ".bka$", recursive = T)
  
  data <- tibble::enframe(data_list) %>%
    tidyr::unnest(value)
  
  return(data)
}





#' Calculate moving average
#'
#' This function calculates a moving average of a vector for smoothing data 
#' (e.g. time series).
#' @importFrom magrittr %>%
#' @param x Vector for which to calculate moving average
#' @param window Width of the window for the moving average. Number of 
#' neighboring values that will be averaged at each position of x.
#' @export
#' @examples
#' moving_av(x = rnorm(100), window = 5)

moving_av <- function(x, window = 5){
  stats::filter(x, rep(1 / window, window), sides = 2)
}


#' Reduce number of time points by averaging values of multiple time points
#'
#' This function reduces the number of time points in a data frame by averaging 
#' values from multiple time points.
#' @importFrom magrittr %>%
#' @param data A data frame with two columns: time and value
#' @param width Number of time/value pairs that will be averaged into one pair
#' @export
#' @examples
#' timeseries_av(data = tibble(time = 1:100, value = rnorm(100)), window = 5)

timeseries_av <- function(data, window = 5){
  rows <- floor(nrow(data)/window)*window
  data <- data[1:rows,]
  
  var_1 <- colnames(data)[1]
  var_2 <- colnames(data)[2]
  
  data$grouping_column_ts <- rep(1:floor(nrow(data)/window), each = window)
  
  data <- data %>%
    dplyr::group_by(grouping_column_ts) %>%
    dplyr::summarise({{var_1}} := mean(.data[[var_1]]),
                     {{var_2}} := mean(.data[[var_2]])) %>%
    dplyr::select(-grouping_column_ts)
  
  return(data)
}


