create_experiment_dir <- function(directory_path, subdirectories){
    if (!dir.exists(directory_path)) {
        dir.create(directory_path, recursive = TRUE)
    } else {
        cat(sprintf("Directory %s already exists.\n", directory_path))
    }
    sapply(file.path(directory_path, subdirectories), function(path_to_create) {
        if (!dir.exists(path_to_create)) {
            dir.create(path_to_create)
        } else {
            cat(sprintf("Directory %s already exists.\n", path_to_create))
        }
    })
}
