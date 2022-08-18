#' Import bio-kine ascii files (.bka)
#'
#' This function allows you to import bio-kine ascii files (.bka). 
#' These are the standard output of BioLogic stopped flow devices.
#' @importFrom magrittr %>%
#' @param path Path to the input .bka file.
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
#' @param path Path to the input directory containing the .bka files to read.
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
#' @param x Vector for which to calculate moving average.
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
#' @param data A data frame with two columns: time and value.
#' @param width Number of time/value pairs that will be averaged into one pair.
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





#' Monte Carlo simulations of linear least squares models
#'
#' Perform Monte Carlo (MC) simulations of non linear least squares models.
#' @importFrom magrittr %>%
#' @param model A nls object to perform MC calculations with.
#' @param runs Number of MC runs to perform (defaults to 1000).
#' @export
#' @examples
#' n <- nls(mpg ~ k * e ^ wt, data = mtcars, start = list(k = 1, e = 2))
#' nls_MC(n)



nls_MC <- function(model, runs = 1000){
  
  sigma        <- broom::glance(model)$sigma
  start        <- summary(model)$param[,1]
  names(start) <- rownames(summary(model)$param)
  params_MC    <- tibble::tibble()
  data         <- broom::augment(model)
  var          <- model$m$formula() %>% as.character() %>% .[2]
  opt          <- model$control
  opt$warnOnly <- FALSE
  
  for(i in 1:runs){
    data_MC <- data
    data_MC[var] <- data$.fitted + rnorm(n = nrow(data), mean = 0, sd = sigma)
    try <- try(
      nls <- nls(model$m$formula(),
                 data    = data_MC, 
                 start   = start, 
                 control = opt),
      silent = T
    )
    if(!"try-error" %in% class(try)){
      params_MC <- dplyr::bind_rows(params_MC, broom::tidy(try) %>% tibble::add_column(run = i))
    }
  }
  if(nrow(params_MC)/length(start) < runs) warning(paste("only", nrow(params_MC)/length(start), "of", runs, "runs converged sucessfully"))
  params_MC <- params_MC %>%
    tibble::add_column(init_estimate = rep(start, nrow(params_MC)/length(start)))
  return(params_MC)
}





#' Monte Carlo simulations of linear least squares models fitted with nlsLM
#' from the minpack.lm package
#'
#' Perform Monte Carlo (MC) simulations of non linear least squares models 
#' fitted with `nlsLM` from the `minpack.lm` package
#' @importFrom magrittr %>%
#' @param model A nlsLM object to perform MC calculations with.
#' @param runs Number of MC runs to perform (defaults to 1000).
#' @export
#' @examples
#' n <- nls(mpg ~ k * e ^ wt, data = mtcars, start = list(k = 1, e = 2))
#' nls_MC(n)



nlsLM_MC <- function(model, runs = 1000){
  
  sigma        <- broom::glance(model)$sigma
  start        <- summary(model)$param[,1]
  names(start) <- rownames(summary(model)$param)
  params_MC    <- tibble::tibble()
  data         <- broom::augment(model)
  var          <- model$m$formula() %>% as.character() %>% .[2]
  opt          <- model$control
  opt$warnOnly <- FALSE
  
  for(i in 1:runs){
    data_MC <- data
    data_MC[var] <- data$.fitted + rnorm(n = nrow(data), mean = 0, sd = sigma)
    try <- try(
      nls <- nls(model$m$formula(),
                 data    = data_MC, 
                 start   = start, 
                 control = opt),
      silent = T
    )
    if(!"try-error" %in% class(try)){
      params_MC <- dplyr::bind_rows(params_MC, broom::tidy(try) %>% tibble::add_column(run = i))
    }
  }
  if(nrow(params_MC)/length(start) < runs) warning(paste("only", nrow(params_MC)/length(start), "of", runs, "runs converged sucessfully"))
  params_MC <- params_MC %>%
    tibble::add_column(init_estimate = rep(start, nrow(params_MC)/length(start)))
  return(params_MC)
}







#' Summarize Monte Carlo simulations of linear least squares models
#'
#' Summarizes the output of the nls_MC() function. For each fitted parameter 
#' mean, original estimate, standard deviation and confidence intervals are computed.
#' @importFrom magrittr %>%
#' @param data A tibble generated by nls_MC().
#' @param conf.level Level for confidence interval calculation (defaults to 0.95).
#' @export
#' @examples
#' n <- nls(mpg ~ k * e ^ wt, data = mtcars, start = list(k = 1, e = 2))
#' s <- nls_MC(n)
#' summarise_nls_MC(s, conf.level = 0.99)



summarise_nls_MC <- function(data, conf.level = 0.95){
  
  lo <- (1-conf.level)/2
  hi <- 1-(1-conf.level)/2
  
  summary <- data %>%
    dplyr::group_by(term) %>%
    dplyr::summarise(
      CI_lo     = quantile(estimate, lo),
      CI_hi     = quantile(estimate, hi),
      std.error = sd(estimate),
      estimate  = mean(estimate),
      init_estimate = mean(init_estimate)) %>%
    dplyr::select(term, estimate, init_estimate, std.error, CI_lo, CI_hi)
  
  return(summary)
}



#' Read TECAN kinetics
#'
#' `read_kinetic_TECAN()` takes a vector of TECAN xlsx kinetic files and vector 
#' of plate layout files (as xlsx), imports data and adds the information 
#' from the plate layouts.
#' @importFrom magrittr %>%
#' @param TECAN_files A vector of TECAN xlsx kinetic files to import.
#' @param layout_files A vector of the respective plate layout files.
#' @param n_cond Number of conditions per well in the plate layout files.
#' @param plate_type Type of plates used. Has to be either "96" or "384".
#' (default is "96")
#' @param pre_trigger Number of pre-trigger time points to discard.
#' @param delay Delay time between starting the kinetic and the first time point 
#' (after pre-trigger) in seconds. Can be either a single value which is used
#' for all plates, or a vector with the same length as `TECAN_files` and 
#' `layout_files` to specify different delay times for each plate.
#' @export
#' @examples
#' read_kinetic_TECAN(TECAN_files = list.files("raw_data", pattern = ".xlsx", full.names = T),
#'                    layout_files = list.files("plate_layouts", pattern = ".xlsx", full.names = T), 
#'                    pre_trigger = 2, delay = c(12,24))



read_kinetic_TECAN <- function(TECAN_files, layout_files, n_cond = 3, 
                               plate_type = "96", pre_trigger = 0, delay = 0){
  
  TECAN_files <- TECAN_files[!grepl("~", TECAN_files)]
  layout_files <- layout_files[!grepl("~", layout_files)]
  
  if(plate_type == "96"){
    cols <- as.character(1:12) %>% stringr::str_pad(2, "left", "0")
    rows <- LETTERS[1:8]
    wells <- c()
    for(i in rows){
      wells <- c(wells, paste0(i,cols))
    }
    range <- paste0("A1:",LETTERS[12],n_cond*8)
  }
  
  if(plate_type == "384"){
    cols <- as.character(1:24) %>% stringr::str_pad(2, "left", "0")
    rows <- LETTERS[1:16]
    wells <- c()
    for(i in rows){
      wells <- c(wells, paste0(i,cols))
    }
    range <- paste0("A1:",LETTERS[24],n_cond*16)
  }
  
  conditions <- tibble::tibble()
  
  for(i in 1:length(layout_files)){
    
    plate_layout <- readxl::read_excel(layout_files[i], 
                                       range = range,
                                       col_names = cols) %>%
      tibble::add_column(cond_n = paste0("cond_", rep(1:n_cond, nrow(.)/n_cond))) %>%
      tibble::add_column(row = rep(rows, each = n_cond)) %>%
      tidyr::pivot_longer(!c("cond_n", "row"), 
                          names_to = "col", 
                          values_to = "cond") %>%
      tibble::add_column(plate = i) %>%
      dplyr::mutate(well = paste0(row, col),
                    p_well = paste(plate, well, sep = "_")) %>%
      tidyr::pivot_wider(names_from = cond_n, 
                         values_from = "cond")
    
    conditions <- dplyr::bind_rows(conditions, plate_layout)
    
  }
  
  conditions <- conditions %>% 
    dplyr::arrange(p_well) %>% 
    stats::na.omit() %>%
    readr::type_convert(col_types = readr::cols())
  
  data <- tibble::tibble()
  
  for(i in 1:length(TECAN_files)) {
    start <- readxl::read_excel(path = TECAN_files[i], range = "A1:A10000", col_names = c("colname"))[[1]] %>% 
      grep(pattern = "Cycle Nr.")
    
    d <- readxl::read_excel(TECAN_files[i], skip = start[1] - 1, n_max = start[2] - start[1] -3, na = "Invalid") %>%
      dplyr::rename_at(dplyr::vars(1:3), ~ c("cycle","time","temp")) %>%
      dplyr::rename_with(.cols = !c("cycle","time","temp"), ~ wells) %>%
      tidyr::pivot_longer(cols = !c("cycle","time","temp"), names_to = "well", values_to = "value") %>%
      tibble::add_column(plate = i) %>%
      dplyr::mutate(p_well = paste(plate, well, sep = "_"))
    
    data <- dplyr::bind_rows(data, d)
    
  }
  
  data <- stats::na.omit(data)
  
  data <- dplyr::left_join(data, conditions, by = c("well", "plate", "p_well"))
  
  delays <- tibble::tibble(plate = 1:length(TECAN_files),
                           delay = delay)
  
  data <- data %>%
    readr::type_convert(col_types = readr::cols()) %>%
    dplyr::filter(cycle > pre_trigger) %>%
    dplyr::group_by(p_well) %>%
    dplyr::mutate(time = time - min(time)) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(delays, by = c("plate")) %>%
    dplyr::mutate(time = time + delay) %>%
    dplyr::select(-delay) %>%
    dplyr::arrange(p_well, time)
  
  
  return(data)
  
}

#' Title
#'
#' Description
#' @importFrom magrittr %>%
#' @export

fancy_scientific <- function(l) {
  l <- format(l, scientific = TRUE)
  l <- gsub("^(.*)e", "'\\1'e", l)
  l <- gsub("e", "%*%10^", l)
  l <- stringr::str_remove(l, "[+]")
  parse(text=l)
}

#' Title
#'
#' Description
#' @importFrom magrittr %>%
#' @export

fancy_scientific2 <- function(l){
  l <- format(l, scientific = TRUE)
  l <- gsub("^(.*)e", "'\\1'e", l)
  l <- gsub("e", "%*%10^", l)
  l <- stringr::str_remove(l, "[+]")
  l <- stringr::str_extract(l, "10\\^-?\\d+")
  #l <- stringr::str_replace(l, "10\\^00", "1")
  parse(text=l)
}


#' Plot plate overview
#'
#' `plot_kinetics_plate()` takes data imported by `read_kinetic_TECAN()` and
#' generates overview plots for each plate in a specified directory
#' @importFrom magrittr %>%
#' @import ggplot2
#' @param data A data frame imported by `read_kinetic_TECAN()`
#' @param plate_type Type of plates used. Has to be either "96" or "384". 
#' (default is "96")
#' @param path Output path for the generated plots
#' @param width Width of the plots in mm, defaults to 400
#' @param height Height of the plots in mm, defaults to 300
#' @export


plot_kinetics_plate <- function(data, plate_type = "96", path = "raw_plots", 
                                width = 450, height = 300){
  
  if(plate_type == "96"){
    cols <- as.character(1:12) %>% stringr::str_pad(2, "left", "0")
    rows <- LETTERS[1:8]
    wells <- c()
    for(i in rows){
      wells <- c(wells, paste0(i,cols))
    }
  }
  
  if(plate_type == "384"){
    cols <- as.character(1:24) %>% stringr::str_pad(2, "left", "0")
    rows <- LETTERS[1:16]
    wells <- c()
    for(i in rows){
      wells <- c(wells, paste0(i,cols))
    }
  }
  
  data$row <- factor(data$row)
  data$row <- forcats::fct_expand(data$row, rows) %>% forcats::fct_relevel(rows)
  data$col <- factor(data$col)
  data$col <- forcats::fct_expand(data$col, cols) %>% forcats::fct_relevel(cols)
  
  for(i in unique(data$plate)){
    q <- data %>%
      dplyr::filter(plate == i) %>%
      ggplot(aes(x = time, y = value))+
      geom_line()+
      facet_grid(row~col, drop = F)+
      labs(title = paste("Plate",i))+
      theme(axis.text.x = element_text(angle = -45, hjust = 0))
    ggsave(paste0(path,"/plate_",i,".pdf"), width = width, height = height, units = "mm")
  }
  
}
