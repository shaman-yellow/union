# ==========================================================================
# load union
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

list_unions <- function(path = dirname(find.package("union.utils")), compatible = TRUE)
{
  base <- file.path(path, c("union.publish", "union.project"))
  if (!all(file.exists(base))) {
    stop('!all(file.exists(base)).')
  }
  series <- list.files(path, "union\\.series\\.", full.names = TRUE)
  if (compatible) {
    com <- list.files(path, "utils\\.tool", full.names = TRUE)
  } else {
    com <- NULL
  }
  c(base, series, com)
}

load_unions <- function(path = dirname(find.package("union.utils")), 
  pkgs = list_unions(path))
{
  for (i in pkgs) {
    devtools::load_all(i)
  }
}

new_workflow <- function(name, path = dirname(find.package("union.utils")), try_open = TRUE)
{
  max_series <- .guess_number_from_files(path, "union\\.series\\.")
  seq <- sprintf("%02d", max_series)
  path_max_series <- file.path(path, glue::glue("union.series.{seq}"))
  path_r_max_series <- file.path(path_max_series, "R")
  max_wf <- .guess_number_from_files(path_r_max_series, "workflow_[0-9]+_.*\\.R")
  if (max_wf + 1L >= (max_series + 1) * 10L) {
    message(
      glue::glue(
        "Now series is '{seq}', workflow is '{max_wf}', Need to bulid new 'union.series'"
      )
    )
    seq <- sprintf("%02d", max_series + 1L)
    path_max_series <- file.path(path, glue::glue("union.series.{seq}"))
    path_r_max_series <- file.path(path_max_series, "R")
    new_package(path_max_series, NULL, c("union.project", "data.table"))
  }
  seq_wf <- sprintf("%03d", max_wf + 1L)
  file_new_workflow <- file.path(
    path_r_max_series, glue::glue("workflow_{seq_wf}_{name}.R")
  )
  lines <- readLines(file.path(.expath, "job_templ", "workflow.R"))
  lines <- paste0(lines, collapse = "\n")
  lines <- glue::glue(lines, name = name, .open = ".{{{", .close = "}}}.")
  writeLines(lines, file_new_workflow)
  if (try_open) {
    if (interactive() && requireNamespace("nvimcom", quietly = TRUE)) {
      SendCmdToNvim_lua(glue::glue("vim.cmd([[tabe {file_new_workflow}]])"))
    }
  }
}

.guess_number_from_files <- function(path, pattern) {
  if (!file.exists(path)) {
    stop('!file.exists(path)')
  }
  alls <- list.files(path, pattern)
  num <- as.integer(stringr::str_extract(alls, "[0-9]+"))
  num <- num[!is.na(num)]
  if (length(num)) {
    max(num)
  } else {
    stop('length(num).')
  }
}
