# config.R
script_configs <- list(
    "000_setupExperimentDir" = list(
      args = list(
        directory = list(
          type = "string",
          required = TRUE,
          validation = function(value) dir.exists(file.path(Sys.getenv('HOME'), 'data', value)),
          error_message = "Directory does not exist in ~/data/"
        )
      )
    ),
  "003_updateSampleGrid" = list(
    args = list(
      directory = list(
        type = "string",
        required = TRUE,
        validation = function(value) dir.exists(file.path(Sys.getenv('HOME'), 'data', value)),
        error_message = "Directory does not exist in ~/data/"
      )
    )
  ),
  "001_plotAllSampleTracks" = list(
    args = list(
      directory = list(
        type = "string",
        required = TRUE,
        validation = function(value) dir.exists(file.path(Sys.getenv('HOME'), 'data', value)),
        error_message = "Directory does not exist in ~/data/"
      ),
      chromosome = list(
        type = "integer",
        required = FALSE,
        default = 10,
        validation = function(value) value >= 1 && value <= 16,
        error_message = "Chromosome must be between 1 and 16"
      )
    )
  )
)
