#' Input Validation Functions
validate_script_args <- function(...) {
    log_info("Validating script arguments")
    
    # Get script info
    script_info <- get_script_info()
    
    if (is.null(script_info$config)) {
        log_error("No configuration found for script:", script_info$name)
        stop("Invalid script configuration")
    }
    
    # Validate arguments
    validate_args(script_info$config$args, list(...))
}

validate_args <- function(config, args) {
    log_info("Validating arguments")
    
    for (arg_name in names(config)) {
        arg_config <- config[[arg_name]]
        value <- args[[arg_name]]
        
        validate_argument(arg_name, value, arg_config)
    }
    
    TRUE
}

validate_argument <- function(name, value, config) {
    # Check required argument
    if (is.null(value)) {
        if (config$required) {
            log_error("Required argument missing:", name)
            stop(sprintf("Required argument '%s' is missing", name))
        }
        if (!is.null(config$default)) {
            value <- config$default
        }
        return(value)
    }
    
    # Validate type
    if (typeof(value) != config$type) {
        log_error("Invalid type for argument:", name)
        log_error("Expected:", config$type, "Got:", typeof(value))
        stop(sprintf("Invalid type for argument '%s'", name))
    }
    
    # Run custom validation
    if (!config$validation(value)) {
        log_error("Validation failed for argument:", name)
        log_error(config$error_message)
        stop(config$error_message)
    }
    
    value
}
