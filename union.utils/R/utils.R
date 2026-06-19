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

plot_volcano_facet <- function(top_table,
  facet,
  label = "hgnc_symbol",
  use = "adj.P.Val",
  fc = .3,
  cut.p = .05,
  n_top = 10L,
  seed = 2L,
  use.fc = "logFC",
  label.fc = "log2(FC)",
  label.p = paste0("-log10(", use, ")"),
  keep_cols = FALSE,
  pal = NULL,
  mode_fc = 0L,
  f.nudge = .5,
  nudge_y = 0,
  label_by = c("p", "abs_fc"),
  facet_scales = "free_y",
  facet_ncol = NULL,
  show_count = TRUE,
  count_size = 3.2,
  label_size = 3,
  show_legend = TRUE,
  max.overlaps = Inf,
  p_floor = .Machine$double.xmin
)
{
  set.seed(seed)
  label_by <- match.arg(label_by)
  message("Plot facet Volcano.")

  if (!any(label == colnames(top_table))) {
    if (any("rownames" == colnames(top_table))) {
      label <- "rownames"
    }
  }

  vec_need <- unique(c(label, facet, use.fc, use))
  vec_miss <- setdiff(vec_need, colnames(top_table))
  if (length(vec_miss) > 0L) {
    stop(glue::glue(
      "Missing required column(s): {paste(vec_miss, collapse = ', ')}."
    ))
  }

  if (!keep_cols) {
    data <- dplyr::select(
      top_table,
      dplyr::all_of(vec_need)
    )
  } else {
    data <- top_table
  }

  data <- dplyr::filter(
    data,
    !is.na(!!rlang::sym(use)),
    !is.na(!!rlang::sym(use.fc)),
    !is.na(!!rlang::sym(facet))
  )

  data <- dplyr::mutate(
    data,
    .p = !!rlang::sym(use),
    .fc = !!rlang::sym(use.fc),
    .p_plot = pmax(.p, p_floor),
    .neg_log10_p = -log10(.p_plot),
    .abs_fc = abs(.fc),
    change = ifelse(
      .fc > abs(fc) & .p < cut.p,
      "up",
      ifelse(
        .fc < -abs(fc) & .p < cut.p,
        "down",
        "stable"
      )
    ),
    change = factor(change, levels = c("down", "stable", "up"))
  )

  if (is.null(pal)) {
    pal <- c("#053061FF", "#67001FFF")
  }

  if (mode_fc == 0L) {
    xintercept <- c(-abs(fc), abs(fc))
  } else if (mode_fc == 1L) {
    xintercept <- abs(fc)
  } else {
    xintercept <- -abs(fc)
  }

  data_lab <- dplyr::filter(
    data,
    change %in% c("up", "down")
  )

  if (label_by == "p") {
    data_lab <- dplyr::arrange(
      data_lab,
      !!rlang::sym(facet),
      change,
      .p,
      dplyr::desc(.abs_fc)
    )
  } else {
    data_lab <- dplyr::arrange(
      data_lab,
      !!rlang::sym(facet),
      change,
      dplyr::desc(.abs_fc),
      .p
    )
  }

  data_lab <- dplyr::distinct(
    data_lab,
    !!rlang::sym(facet),
    change,
    !!rlang::sym(label),
    .keep_all = TRUE
  )

  data_lab <- dplyr::group_by(
    data_lab,
    !!rlang::sym(facet),
    change
  )

  data_lab <- dplyr::slice_head(
    data_lab,
    n = n_top
  )

  data_lab <- dplyr::ungroup(data_lab)

  if (nrow(data_lab) > 0L) {
    data_lab <- dplyr::group_by(
      data_lab,
      !!rlang::sym(facet)
    )

    data_lab <- dplyr::mutate(
      data_lab,
      .nudge_base = stats::median(.abs_fc, na.rm = TRUE),
      .nudge_base = ifelse(
        is.finite(.nudge_base) & .nudge_base > 0,
        .nudge_base,
        max(.abs_fc, abs(fc), na.rm = TRUE)
      ),
      .nudge_base = ifelse(
        is.finite(.nudge_base) & .nudge_base > 0,
        .nudge_base,
        .5
      ),
      .nudge_x = ifelse(
        change == "up",
        .nudge_base * abs(f.nudge),
        -.nudge_base * abs(f.nudge)
      )
    )

    data_lab <- dplyr::ungroup(data_lab)
  }

  p <- ggplot(
    data,
    aes(
      x = .fc,
      y = .neg_log10_p,
      color = change
    )
  ) +
    geom_point(alpha = .8, stroke = 0, size = 1.5) +
    scale_color_manual(
      values = c(
        down = pal[1L],
        stable = "grey90",
        up = pal[2L]
      ),
      drop = FALSE
    ) +
    geom_hline(
      yintercept = -log10(cut.p),
      linetype = 4L,
      size = .8
    ) +
    geom_vline(
      xintercept = xintercept,
      linetype = 4L,
      size = .8
    ) +
    facet_wrap(
      stats::as.formula(paste0("~", facet)),
      scales = facet_scales,
      ncol = facet_ncol
    ) +
    labs(
      x = label.fc,
      y = label.p
    ) +
    rstyle("theme") +
    geom_blank()

  if (nrow(data_lab) > 0L) {
    p <- p + ggrepel::geom_label_repel(
      data = data_lab,
      nudge_x = data_lab$.nudge_x,
      nudge_y = nudge_y,
      show.legend = FALSE,
      max.overlaps = max.overlaps,
      seed = seed,
      size = label_size,
      aes(
        x = .fc,
        y = .neg_log10_p,
        label = !!rlang::sym(label)
      )
    )
  }

  if (isTRUE(show_count)) {
    message("Use count in volcano plot.")
    data_count <- dplyr::group_by(
      data,
      !!rlang::sym(facet)
    )

    data_count <- dplyr::summarise(
      data_count,
      n_up = sum(change == "up", na.rm = TRUE),
      n_down = sum(change == "down", na.rm = TRUE),
      .groups = "drop"
    )

    data_count <- dplyr::mutate(
      data_count,
      .x = -Inf,
      .y = -Inf,
      .count_label = paste0(
        "Up: ", n_up,
        "\nDown: ", n_down
      )
    )

    p <- p + geom_label(
      data = data_count,
      show.legend = FALSE,
      inherit.aes = FALSE,
      size = count_size,
      label.size = 0,
      alpha = .85,
      hjust = -.05,
      vjust = -.8,
      aes(
        x = .x,
        y = .y,
        label = .count_label
      )
    )
  }

  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }

  p
}

pg_local_recode <- function() {
  conda <- getOption("conda", "~/miniconda3")
  lst <- list(
    fusion = "~/fusion_twas",
    ldscPython = "{conda}/bin/conda run -n ldsc python",
    ldsc = "~/ldsc",
    annovar = "~/disk_sda1/annovar",
    vep = "~/ensembl-vep/vep",
    vep_cache = "~/disk_sda1/.vep",
    python = "{conda}/bin/python3",
    conda = "{conda}/bin/conda",
    conda_env = "{conda}/envs",
    qiime = "{conda}/bin/conda run -n qiime2 qiime",
    musitePython = "{conda}/bin/conda run -n musite python3",
    musitePTM = "~/MusiteDeep_web/MusiteDeep/predict_multi_batch.py",
    musitePTM2S = "~/MusiteDeep_web/PTM2S/ptm2Structure.py",
    hobEnv = "hobpre",
    hobPython = "{conda}/bin/conda run -n hobpre python",
    hobPredict = "~/HOB/HOB_predict.py",
    hobModel = "~/HOB/model",
    hobExtra = "~/HOB/pca_hob.m",
    dl = normalizePath("~/D-GCAN/DGCAN", mustWork = FALSE),
    dl_dataset = normalizePath("~/D-GCAN/dataset", mustWork = FALSE),
    dl_model = normalizePath("~/D-GCAN/DGCAN/model", mustWork = FALSE),
    scfeaPython = "{conda}/bin/conda run -n scFEA python",
    scfea = "~/scFEA/src/scFEA.py",
    scfea_db = "~/scFEA/data",
    musiteModel = normalizePath("~/MusiteDeep_web/MusiteDeep/models", mustWork = FALSE),
    vina = "vina",
    docking_python = "conda run -n docking python",
    mk_prepare_ligand.py = "conda run -n docking mk_prepare_ligand.py",
    mk_prepare_receptor.py = "conda run -n docking mk_prepare_receptor.py",
    prepare_receptor = "~/operation/ADFRsuite_x86_64Linux_1.0/bin/prepare_receptor",
    prepare_gpf.py = "~/operation/ADFRsuite_x86_64Linux_1.0/bin/pythonsh ~/autodock_vina/example/autodock_scripts/prepare_gpf.py",
    autogrid4 = "~/operation/ADFRsuite_x86_64Linux_1.0/bin/autogrid4",
    pymol = "QT_QPA_PLATFORM=xcb pymol",
    scsaEnv = "scsa",
    scsa = "{conda}/bin/conda run -n scsa python3 ~/SCSA/SCSA.py",
    scsa_db = "~/SCSA/whole_v2.db",
    # sirius = .prefix("sirius/bin/sirius", "op"),
    obgen = "obgen"
  )
  envir <- environment()
  lapply(lst, glue::glue, .envir = envir)
}

pg_remote_recode <- function() {
  conda_remote <- getOption("conda_remote", "~/miniconda3")
  lst <- list(
    # vina = "{conda_remote}/bin/conda run -n vina vina",
    vina = "vina",
    qiime = "{conda_remote}/bin/conda run -n qiime2 qiime",
    fastp = "{conda_remote}/bin/conda run -n base fastp",
    bcftools = "{conda_remote}/bin/conda run -n base bcftools",
    elprep = "{conda_remote}/bin/conda run -n base elprep",
    # biobakery_workflows = "{conda_remote}/bin/conda run -n biobakery biobakery_workflows",
    bowtie2 = "{conda_remote}/bin/conda run -n base bowtie2",
    samtools = "{conda_remote}/bin/conda run -n base samtools",
    metaphlan = "{conda_remote}/bin/conda run -n mpa metaphlan",
    Rscript = "{conda_remote}/bin/conda run -n r4-base Rscript",
    merge_metaphlan_tables.py = "{conda_remote}/bin/conda run -n mpa merge_metaphlan_tables.py",
    sirius = "~/operation/sirius/bin/sirius",
    scfeaPython = "{conda_remote}/bin/conda run -n scFEA python",
    scfea = "~/scFEA/src/scFEA.py",
    scfea_db = "~/scFEA/data"
  )
  envir <- environment()
  lapply(lst, glue::glue, .envir = envir)
}

diagnose_object_links <- function(x,
  max_depth = 4L,
  inspect_attributes = TRUE,
  inspect_env = c("record", "children", "none"),
  skip_env_kinds = c("namespace", "package", "global", "base", "empty"),
  skip_attr_names = c("names", "dim", "dimnames", "class", "row.names", "levels"),
  max_env_bindings = 50L,
  max_children_per_node = Inf,
  max_nodes = 5000L,
  calc_serialize = FALSE,
  serialize_depth = 2L)
{
  inspect_env <- match.arg(inspect_env)
  has_lobstr <- base::requireNamespace("lobstr", quietly = TRUE)

  .human_size <- function(n) {
    if (length(n) == 0L || is.na(n)) {
      return(NA_character_)
    }
    units <- c("B", "KB", "MB", "GB", "TB")
    n2 <- as.numeric(n)
    idx <- 1L
    while (n2 >= 1024 && idx < length(units)) {
      n2 <- n2 / 1024
      idx <- idx + 1L
    }
    paste0(format(round(n2, 2L), nsmall = 2L), " ", units[[idx]])
  }

  .safe_object_size <- function(obj) {
    res <- try(as.numeric(utils::object.size(obj)), silent = TRUE)
    if (inherits(res, "try-error")) {
      return(NA_real_)
    }
    res
  }

  .safe_obj_size <- function(obj) {
    if (!has_lobstr) {
      return(NA_real_)
    }
    res <- try(as.numeric(lobstr::obj_size(obj)), silent = TRUE)
    if (inherits(res, "try-error")) {
      return(NA_real_)
    }
    res
  }

  .safe_addr <- function(obj) {
    if (!has_lobstr) {
      return(NA_character_)
    }
    res <- try(as.character(lobstr::obj_addr(obj)), silent = TRUE)
    if (inherits(res, "try-error")) {
      return(NA_character_)
    }
    res
  }

  .short_class <- function(obj) {
    cls <- try(class(obj), silent = TRUE)
    if (inherits(cls, "try-error") || length(cls) == 0L) {
      return(NA_character_)
    }
    paste(cls, collapse = "/")
  }

  .env_kind <- function(env) {
    if (!is.environment(env)) {
      return(NA_character_)
    }
    if (identical(env, emptyenv())) {
      return("empty")
    }
    if (identical(env, globalenv())) {
      return("global")
    }
    if (identical(env, baseenv())) {
      return("base")
    }
    env_name <- environmentName(env)
    if (grepl("^package:", env_name)) {
      return("package")
    }
    if (base::isNamespace(env) || grepl("^namespace:", env_name)) {
      return("namespace")
    }
    "ordinary"
  }

  .env_label <- function(env) {
    if (!is.environment(env)) {
      return(NA_character_)
    }
    env_name <- environmentName(env)
    if (!identical(env_name, "")) {
      return(env_name)
    }
    paste0("<", .env_kind(env), "_env>")
  }

  .safe_serialize_size <- function(obj, depth) {
    if (!isTRUE(calc_serialize) || depth > serialize_depth) {
      return(NA_real_)
    }
    if (is.environment(obj) && .env_kind(obj) %in% skip_env_kinds) {
      return(NA_real_)
    }
    res <- try(length(serialize(obj, NULL, xdr = FALSE)), silent = TRUE)
    if (inherits(res, "try-error")) {
      return(NA_real_)
    }
    as.numeric(res)
  }

  .path_list <- function(parent, nm, idx) {
    if (!is.null(nm) && !is.na(nm) && nzchar(nm)) {
      if (grepl("^[.A-Za-z][.A-Za-z0-9_]*$", nm)) {
        return(paste0(parent, "$", nm))
      }
      return(paste0(parent, "[[\"", nm, "\"]]"))
    }
    paste0(parent, "[[", idx, "]]")
  }

  .path_slot <- function(parent, nm) {
    paste0(parent, "@", nm)
  }

  .path_attr <- function(parent, nm) {
    paste0(parent, " attr(", nm, ")")
  }

  .get_children <- function(obj, path, depth) {
    if (depth >= max_depth) {
      return(list())
    }

    out <- list()

    if (methods::is(obj, "function") && !identical(inspect_env, "none")) {
      env <- environment(obj)
      out[[length(out) + 1L]] <- list(
        path = paste0(path, " <closure_env>"),
        value = env
      )
    }

    if (is.environment(obj) && identical(inspect_env, "children")) {
      kind <- .env_kind(obj)
      if (!(kind %in% skip_env_kinds)) {
        nms <- try(ls(envir = obj, all.names = TRUE), silent = TRUE)
        if (!inherits(nms, "try-error") && length(nms) > 0L) {
          nms <- utils::head(nms, max_env_bindings)
          for (nm in nms) {
            is_active <- try(bindingIsActive(nm, obj), silent = TRUE)
            if (inherits(is_active, "try-error") || isTRUE(is_active)) {
              next
            }
            val <- try(get(nm, envir = obj, inherits = FALSE), silent = TRUE)
            if (inherits(val, "try-error")) {
              next
            }
            out[[length(out) + 1L]] <- list(
              path = paste0(path, "$", nm),
              value = val
            )
          }
        }
      }
    }

    if (isS4(obj)) {
      slots <- slotNames(obj)
      if (length(slots) > 0L) {
        slots <- utils::head(slots, max_children_per_node)
        for (nm in slots) {
          val <- try(methods::slot(obj, nm), silent = TRUE)
          if (inherits(val, "try-error")) {
            next
          }
          out[[length(out) + 1L]] <- list(
            path = .path_slot(path, nm),
            value = val
          )
        }
      }
    } else if (is.list(obj) || is.pairlist(obj)) {
      nms <- names(obj)
      n <- length(obj)
      if (is.infinite(max_children_per_node)) {
        idx <- seq_len(n)
      } else {
        idx <- seq_len(min(n, max_children_per_node))
      }
      for (i in idx) {
        nm <- NA_character_
        if (!is.null(nms) && length(nms) >= i) {
          nm <- nms[[i]]
        }
        out[[length(out) + 1L]] <- list(
          path = .path_list(path, nm, i),
          value = obj[[i]]
        )
      }
    }

    if (isTRUE(inspect_attributes) && !isS4(obj)) {
      attrs <- try(attributes(obj), silent = TRUE)
      if (!inherits(attrs, "try-error") && length(attrs) > 0L) {
        attr_nms <- setdiff(names(attrs), skip_attr_names)
        if (length(attr_nms) > 0L) {
          attr_nms <- utils::head(attr_nms, max_children_per_node)
          for (nm in attr_nms) {
            out[[length(out) + 1L]] <- list(
              path = .path_attr(path, nm),
              value = attrs[[nm]]
            )
          }
        }
      }
    }

    out
  }

  rows <- list()
  visited <- new.env(parent = emptyenv())
  reached_limit <- FALSE

  .walk <- function(obj, path, depth) {
    if (length(rows) >= max_nodes) {
      reached_limit <<- TRUE
      return(invisible(NULL))
    }

    addr <- .safe_addr(obj)
    key <- NA_character_

    is_container <- is.environment(obj) ||
      is.function(obj) ||
      is.list(obj) ||
      isS4(obj) ||
      !is.null(attributes(obj))

    if (isTRUE(is_container) && !is.na(addr)) {
      key <- paste0(typeof(obj), ":", addr)
    }

    duplicated_ref <- FALSE
    if (!is.na(key)) {
      if (exists(key, envir = visited, inherits = FALSE)) {
        duplicated_ref <- TRUE
      } else {
        assign(key, TRUE, envir = visited)
      }
    }

    env_kind <- if (is.environment(obj)) .env_kind(obj) else NA_character_
    env_n_objects <- NA_integer_
    if (is.environment(obj) && !(env_kind %in% skip_env_kinds)) {
      nms <- try(ls(envir = obj, all.names = TRUE), silent = TRUE)
      if (!inherits(nms, "try-error")) {
        env_n_objects <- length(nms)
      }
    }

    object_size_bytes <- .safe_object_size(obj)
    true_size_bytes <- .safe_obj_size(obj)
    serialize_bytes <- .safe_serialize_size(obj, depth)

    rows[[length(rows) + 1L]] <<- list(
      path = path,
      depth = depth,
      type = typeof(obj),
      class = .short_class(obj),
      object_size_bytes = object_size_bytes,
      object_size = .human_size(object_size_bytes),
      true_size_bytes = true_size_bytes,
      true_size = .human_size(true_size_bytes),
      serialize_bytes = serialize_bytes,
      serialize_size = .human_size(serialize_bytes),
      addr = addr,
      duplicated_ref = duplicated_ref,
      is_environment = is.environment(obj),
      is_function = is.function(obj),
      is_s4 = isS4(obj),
      env_kind = env_kind,
      env_label = if (is.environment(obj)) .env_label(obj) else NA_character_,
      env_n_objects = env_n_objects
    )

    if (isTRUE(duplicated_ref)) {
      return(invisible(NULL))
    }

    children <- .get_children(obj, path, depth)
    if (length(children) == 0L) {
      return(invisible(NULL))
    }

    for (child in children) {
      .walk(child$value, child$path, depth + 1L)
    }

    invisible(NULL)
  }

  .walk(x, "x", 0L)

  detail <- do.call(
    rbind,
    lapply(rows, function(z) {
      data.frame(
        path = z$path,
        depth = z$depth,
        type = z$type,
        class = z$class,
        object_size_bytes = z$object_size_bytes,
        object_size = z$object_size,
        true_size_bytes = z$true_size_bytes,
        true_size = z$true_size,
        serialize_bytes = z$serialize_bytes,
        serialize_size = z$serialize_size,
        addr = z$addr,
        duplicated_ref = z$duplicated_ref,
        is_environment = z$is_environment,
        is_function = z$is_function,
        is_s4 = z$is_s4,
        env_kind = z$env_kind,
        env_label = z$env_label,
        env_n_objects = z$env_n_objects,
        stringsAsFactors = FALSE
      )
    })
  )

  size_rank <- detail
  score_size <- size_rank$true_size_bytes
  score_size[is.na(score_size)] <- size_rank$object_size_bytes[is.na(score_size)]
  score_size[is.na(score_size)] <- -Inf
  size_rank <- size_rank[order(-score_size), , drop = FALSE]

  env_refs <- detail[
    detail$is_environment |
      detail$is_function |
      grepl("Environment|env|closure_env|quosure|quo|formula|plot_env",
        detail$path,
        ignore.case = TRUE),
    ,
    drop = FALSE
  ]

  serialize_rank <- detail
  serialize_score <- serialize_rank$serialize_bytes
  serialize_score[is.na(serialize_score)] <- -Inf
  serialize_rank <- serialize_rank[order(-serialize_score), , drop = FALSE]

  summary <- data.frame(
    n_nodes = nrow(detail),
    reached_max_nodes = reached_limit,
    lobstr_available = has_lobstr,
    root_object_size = .human_size(detail$object_size_bytes[[1L]]),
    root_true_size = .human_size(detail$true_size_bytes[[1L]]),
    root_serialize_size = .human_size(detail$serialize_bytes[[1L]]),
    stringsAsFactors = FALSE
  )

  list(
    summary = summary,
    detail = detail,
    top_size = utils::head(size_rank, 30L),
    top_serialize = utils::head(serialize_rank, 30L),
    env_refs = env_refs
  )
}
