
#' @title Get script name
#' @description Extract the full path of the current script
#' @return Character string of the script path
get_script_name <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else if (!is.null(sys.frames()[[1]]$ofile)) {
    return(normalizePath(sys.frames()[[1]]$ofile))
  } else {
    return("interactive")
  }
}

#' @title Get script directory
#' @description Extract the directory of the current script
#' @return Character string of the script directory
get_script_dir <- function() {
  script_path <- get_script_name()
  if (script_path != "interactive") {
    return(dirname(script_path))
  } else {
    return(getwd())
  }
}

#' @title Get script basename
#' @description Extract the basename of the current script without extension
#' @return Character string of the script basename
get_script_basename <- function() {
  script_path <- get_script_name()
  if (script_path != "interactive") {
    return(tools::file_path_sans_ext(basename(script_path)))
  } else {
    return("interactive")
  }
}

#' @title Get script basename
#' @description Extract the basename of the current script without extension
#' @return Character string of the script basename
#validate_script_input("001_plotAllSampleTracks", directory = "my_data", chromosome = 5)
validate_script_input <- function(script_name, ...) {
  config <- script_configs[[script_name]]
  if (is.null(config)) stop(paste("No configuration found for script:", script_name))
  
  args <- list(...)
  
  for (arg_name in names(config$args)) {
    arg_config <- config$args[[arg_name]]
    value <- args[[arg_name]]
    
    if (is.null(value) && arg_config$required) {
      stop(paste("Required argument", arg_name, "is missing"))
    }
    
    if (is.null(value) && !is.null(arg_config$default)) {
      value <- arg_config$default
    }
    
    if (!is.null(value)) {
      if (typeof(value) != arg_config$type) {
        stop(paste("Argument", arg_name, "should be of type", arg_config$type))
      }
      
      if (!arg_config$validation(value)) {
        stop(arg_config$error_message)
      }
    }
  }
  
  return(TRUE)
}
