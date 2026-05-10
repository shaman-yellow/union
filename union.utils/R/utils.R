# ==========================================================================
# 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

find_common_cross_groups <- function(df, ...group) {
  dt <- as.data.table(df)
  cols <- setdiff(names(dt), ...group)
  dt[, .key := do.call(paste, c(.SD, sep = "\r")), .SDcols = cols]
  group_sets <- dt[, .(set = list(unique(.key))), by = ...group]
  common <- Reduce(intersect, group_sets$set)
  dt[.key %in% common][, .key := NULL]
}


  if (.Platform$OS.type != "unix") {
    stop('.Platform$OS.type != "unix".')
  }

# ==========================================================
# General Archive Reader with Shell-level ID Filtering
# Supports: .zip, .tar, .tar.gz, .tgz
# Optimized to suppress Broken pipe warnings
# ==========================================================

.shFilter_read_archive_table_by_id <- function(
  path,
  ids,
  id_col,
  pattern = "\\.gz$",
  sep = "\t")
{
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' required.")
  }

  if (!file.exists(path)) {
    stop("Archive path not found.", call. = FALSE)
  }

  path_abs <- normalizePath(path)

  is_zip <- grepl("\\.zip$", path, ignore.case = TRUE)

  is_tar <- grepl(
    "\\.(tar|tar\\.gz|tgz)$",
    path,
    ignore.case = TRUE
  )

  if (is_zip) {

    members <- utils::unzip(
      path_abs,
      list = TRUE
    )$Name

  } else if (is_tar) {

    members <- utils::untar(
      path_abs,
      list = TRUE
    )

  } else {

    stop("Unsupported archive format.")
  }

  target_members <- members[
    grepl(pattern, members)
  ]

  if (length(target_members) == 0L) {
    stop("No members match the pattern.")
  }

  # -----------------------------
  # iterate archive members
  # -----------------------------
  res_list <- lapply(target_members, function(m) {

    message("Processing member: ", m)

    # extraction command
    ext_cmd <- if (is_zip) {

      paste0(
        "unzip -p ",
        shQuote(path_abs),
        " ",
        shQuote(m)
      )

    } else {

      paste0(
        "tar -xOf ",
        shQuote(path_abs),
        " ",
        shQuote(m)
      )
    }

    # detect whether member itself is gz
    if (grepl("\\.gz$", m, ignore.case = TRUE)) {

      stream_cmd <- paste0(
        ext_cmd,
        " | gunzip -c"
      )

    } else {

      stream_cmd <- ext_cmd
    }

    dat <- tryCatch(

      .shFilter_read_table_by_id(
        file = "",
        ids = ids,
        id_col = id_col,
        sep = sep,
        decompress_cmd = stream_cmd
      ),

      error = function(e) {
        warning(
          "Failed processing member: ",
          m,
          "\n",
          conditionMessage(e)
        )
        NULL
      }
    )

    if (!is.null(dat)) {
      dat$file_member <- m
    }

    dat
  })

  res_list <- Filter(
    function(x) !is.null(x) && nrow(x) > 0L,
    res_list
  )

  if (length(res_list) == 0L) {
    return(NULL)
  }

  final_dt <- data.table::rbindlist(
    res_list,
    use.names = TRUE,
    fill = TRUE
  )

  tibble::as_tibble(final_dt)
}

.shFilter_read_table_by_id <- function(
  file,
  ids,
  id_col,
  sep = "\t",
  decompress_cmd = NULL)
{
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' required.")
  }

  if (!file.exists(file) && is.null(decompress_cmd)) {
    stop("File not found.", call. = FALSE)
  }

  # -----------------------------
  # temporary id file
  # -----------------------------
  id_tmp <- tempfile(fileext = ".ids")

  writeLines(as.character(unique(ids)), id_tmp)

  on.exit(unlink(id_tmp), add = TRUE)

  id_tmp_abs <- normalizePath(id_tmp)

  # -----------------------------
  # build stream command
  # -----------------------------
  if (is.null(decompress_cmd)) {

    file_abs <- normalizePath(file)

    if (grepl("\\.gz$", file, ignore.case = TRUE)) {

      stream_cmd <- paste0(
        "gunzip -c ",
        shQuote(file_abs)
      )

    } else {

      stream_cmd <- paste0(
        "cat ",
        shQuote(file_abs)
      )
    }

  } else {

    stream_cmd <- decompress_cmd
  }

  # -----------------------------
  # read header
  # -----------------------------
  header_cmd <- paste0(
    stream_cmd,
    " 2>/dev/null | head -n 1"
  )

  header <- tryCatch(
    system(header_cmd, intern = TRUE),
    error = function(e) character(0)
  )

  if (length(header) == 0L) {
    return(NULL)
  }

  col_names <- strsplit(
    header,
    sep,
    fixed = TRUE
  )[[1L]]

  col_idx <- which(col_names == id_col)

  if (length(col_idx) == 0L) {

    stop(
      "Column '", id_col, "' not found.",
      call. = FALSE
    )
  }

  # -----------------------------
  # awk filter pipeline
  # -----------------------------
  filter_cmd <- paste0(
    stream_cmd,
    " | awk -F ", shQuote(sep),
    " -v target_col=", col_idx,
    " 'NR==FNR {a[$1]; next} ",
    "(FNR==1) || ($target_col in a)' ",
    shQuote(id_tmp_abs),
    " -"
  )

  dat <- tryCatch(
    data.table::fread(
      cmd = filter_cmd,
      sep = sep,
      header = TRUE,
      data.table = FALSE,
      showProgress = FALSE,
      check.names = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(dat) || nrow(dat) == 0L) {
    return(NULL)
  }

  tibble::as_tibble(dat)
}

.check_pkg <- function(pkg)
{
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      paste0(
        "Package '", pkg,
        "' is required but not installed."
      ),
      call. = FALSE
    )
  }
}

