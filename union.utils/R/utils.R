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

.clear_autor_objects <- function(path, verbose = FALSE) {
  objs <- c("autor", "include")
  try(.remove_r_objects(path, objs, objs, backup = FALSE, verbose = verbose))
}

.remove_r_objects <- function(
    path_pkg,
    vec_generic = NULL,
    vec_method = NULL,
    vec_function = NULL,
    recursive = TRUE,
    backup = FALSE,
    dry_run = FALSE,
    verbose = TRUE)
{
  if (!dir.exists(path_pkg)) {
    stop("Package path does not exist.")
  }

  path_r <- file.path(path_pkg, "R")

  if (!dir.exists(path_r)) {
    stop("Cannot find R directory.")
  }

  vec_file <- list.files(
    path_r,
    pattern = "\\.[Rr]$",
    full.names = TRUE,
    recursive = recursive
  )

  vec_generic <- unique(vec_generic)
  vec_method <- unique(vec_method)
  vec_function <- unique(vec_function)

  vec_generic <- vec_generic[!is.na(vec_generic)]
  vec_method <- vec_method[!is.na(vec_method)]
  vec_function <- vec_function[!is.na(vec_function)]

  .get_call_name <- function(expr)
  {
    if (!is.call(expr)) {
      return(NA_character_)
    }

    obj_head <- expr[[1L]]

    if (is.symbol(obj_head)) {
      return(as.character(obj_head))
    }

    if (is.call(obj_head)) {
      return(deparse1(obj_head))
    }

    NA_character_
  }

  .extract_first_arg <- function(expr)
  {
    if (length(expr) < 2L) {
      return(NA_character_)
    }

    obj <- tryCatch(
      eval(expr[[2L]]),
      error = function(e) NA_character_
    )

    as.character(obj)
  }

  .should_remove_expr <- function(expr)
  {
    .walk_call <- function(node)
    {
      if (!is.call(node)) {
        return(FALSE)
      }

      str_call <- .get_call_name(node)

      if (identical(str_call, "setGeneric")) {

        str_name <- .extract_first_arg(node)

        if (str_name %in% vec_generic) {
          return(TRUE)
        }
      }

      if (identical(str_call, "setMethod")) {

        str_name <- .extract_first_arg(node)

        if (str_name %in% vec_method) {
          return(TRUE)
        }
      }

      if (identical(str_call, "<-") || identical(str_call, "=")) {

        obj_lhs <- node[[2L]]
        obj_rhs <- node[[3L]]

        if (
          is.symbol(obj_lhs) &&
            is.call(obj_rhs) &&
            identical(.get_call_name(obj_rhs), "function")
          ) {

          str_name <- as.character(obj_lhs)

          if (str_name %in% vec_function) {
            return(TRUE)
          }
        }
      }

      vec_child <- as.list(node)[-1L]

      any(
        vapply(
          vec_child,
          .walk_call,
          logical(1L)
        )
      )
    }

    .walk_call(expr)
  }

  .process_file <- function(path_file)
  {
    if (verbose) {
      message(
        glue::glue(
          "[READ] {basename(path_file)}"
        )
      )
    }

    vec_line <- readLines(
      path_file,
      warn = FALSE
    )

    obj_expr <- tryCatch(
      parse(
        text = vec_line,
        keep.source = TRUE
      ),
      error = function(e) e
    )

    if (inherits(obj_expr, "error")) {

      message(
        glue::glue(
          "[SKIP] Parse failed: {basename(path_file)}"
        )
      )

      return(NULL)
    }

    lst_srcref <- attr(obj_expr, "srcref")

    if (is.null(lst_srcref)) {

      message(
        glue::glue(
          "[SKIP] Missing srcref: {basename(path_file)}"
        )
      )

      return(NULL)
    }

    vec_remove <- sapply(
      obj_expr,
      .should_remove_expr
    )

    if (!any(vec_remove)) {

      if (verbose) {

        message(
          glue::glue(
            "[KEEP] No target found: {basename(path_file)}"
          )
        )
      }

      return(FALSE)
    }

    vec_delete_line <- rep(
      FALSE,
      length(vec_line)
    )

    vec_idx <- which(vec_remove)

    invisible(
      sapply(
        vec_idx,
        function(idx)
        {
          obj_sr <- lst_srcref[[idx]]

          idx_start <- obj_sr[[1L]]
          idx_end <- obj_sr[[3L]]

          vec_delete_line[
            idx_start:idx_end
          ] <<- TRUE

          if (verbose) {

            message(
              glue::glue(
                "[REMOVE] {basename(path_file)} :: lines {idx_start}-{idx_end}"
              )
            )
          }

          NULL
        }
      )
    )

    vec_new <- vec_line[!vec_delete_line]

    if (dry_run) {

      message(
        glue::glue(
          "[DRY RUN] {basename(path_file)}"
        )
      )

      return(TRUE)
    }

    if (backup) {

      path_backup <- paste0(
        path_file,
        ".bak"
      )

      ok_backup <- file.copy(
        from = path_file,
        to = path_backup,
        overwrite = TRUE
      )

      if (!ok_backup) {

        message(
          glue::glue(
            "[WARN] Backup failed: {basename(path_file)}"
          )
        )
      }
    }

    writeLines(
      vec_new,
      con = path_file
    )

    if (verbose) {

      message(
        glue::glue(
          "[WRITE] {basename(path_file)}"
        )
      )
    }

    TRUE
  }

  lst_res <- lapply(
    vec_file,
    .process_file
  )

  invisible(lst_res)
}

