# ==========================================================================
# trace S4 method
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -




.autor_s4_trace_runtime <- new.env(parent = emptyenv())
.autor_s4_trace_runtime$current <- NULL
.autor_s4_trace_runtime$n_context <- 0L

.run_s4_trace_enabled <- function()
{
  isTRUE(getOption("autor_s4_trace", FALSE))
}

# .run_make_s4_trace_file <- function(env_trace)
# {
#   path_file <- getOption("autor_s4_trace_file", NULL)
#   if (!is.null(path_file) && length(path_file)) {
#     dir.create(dirname(path_file), FALSE, recursive = TRUE)
#     return(path_file)
#   }
#   dir_trace <- .run_make_s4_trace_dir()
#   file.path(
#     dir_trace,
#     "current.rds"
#   )
# }
#

.run_new_s4_trace_context <- function()
{
  .autor_s4_trace_runtime$n_context <- .autor_s4_trace_runtime$n_context + 1L

  env_trace <- new.env(parent = emptyenv())
  env_trace$id_context <- .autor_s4_trace_runtime$n_context

  env_trace$trace_id <- digest::digest(
    list(
      time = Sys.time(),
      pid = Sys.getpid(),
      host = Sys.info()[["nodename"]],
      n = env_trace$id_context,
      rand = runif(1L)
    ),
    algo = "xxhash64",
    serializeVersion = 3L
  )

  env_trace$records <- list()
  env_trace$stack <- integer(0L)
  env_trace$n_event <- 0L
  env_trace$created_at <- Sys.time()
  env_trace$pid <- Sys.getpid()
  env_trace$host <- Sys.info()[["nodename"]]

  env_trace
}

.run_signature_text <- function(x)
{
  if (is.null(x)) {
    return("")
  }

  if (is(x, "signature")) {
    x <- as.list(x)
  }

  if (is.list(x)) {
    x <- unlist(x, use.names = TRUE)
  }

  if (!length(x)) {
    return("")
  }

  if (is.null(names(x))) {
    return(paste(as.character(x), collapse = "__"))
  }

  paste(
    paste0(names(x), "=", as.character(x)),
    collapse = "__"
  )
}

.run_safe_deparse <- function(x)
{
  paste(
    deparse(x, width.cutoff = 500L),
    collapse = "\n"
  )
}

.run_get_arg_classes <- function(vec_formals, envir)
{
  if (!length(vec_formals)) {
    return(character(0))
  }

  vec_formals <- vec_formals[
    nzchar(vec_formals) &
      vec_formals != "..."
  ]

  if (!length(vec_formals)) {
    return(character(0))
  }

  vapply(
    vec_formals,
    FUN.VALUE = character(1),
    function(name) {
      if (!exists(name, envir = envir, inherits = FALSE)) {
        return("<missing>")
      }

      obj <- try(
        get(name, envir = envir, inherits = FALSE),
        TRUE
      )

      if (inherits(obj, "try-error")) {
        message(glue::glue(
          "Can not get variable from environment, name: {name}"
        ))
        return("<error>")
      }

      paste(class(obj), collapse = "/")
    }
  )
}

.run_start_s4_trace_root <- function(generic)
{
  if (!.run_s4_trace_enabled()) {
    return(FALSE)
  }

  if (is.environment(.autor_s4_trace_runtime$current)) {
    return(FALSE)
  }

  env_trace <- .run_new_s4_trace_context()
  env_trace$root_generic <- generic
  env_trace$root_mode <- "manual"

  .autor_s4_trace_runtime$current <- env_trace

  TRUE
}

.run_end_s4_trace_root <- function(is_root)
{
  if (!isTRUE(is_root)) {
    return(invisible(NULL))
  }

  env_trace <- .autor_s4_trace_runtime$current

  if (is.environment(env_trace)) {
    .run_save_s4_trace_context(env_trace)
  }

  .autor_s4_trace_runtime$current <- NULL

  invisible(NULL)
}

.run_save_s4_trace_context <- function(env_trace)
{
  if (!is.environment(env_trace)) {
    return(invisible(NULL))
  }

  if (!length(env_trace$records)) {
    return(invisible(NULL))
  }

  path_file <- getOption("autor_s4_trace_file", ".s4trace/current.rds")
  dir.create(dirname(path_file), FALSE, recursive = TRUE)

  path_lock <- paste0(path_file, ".lock")

  lock <- try(filelock::lock(path_lock), TRUE)

  if (inherits(lock, "try-error")) {
    message(glue::glue("Can not lock trace file, skip saving trace: {path_file}"))
    return(invisible(NULL))
  }

  on.exit(filelock::unlock(lock), add = TRUE)

  lst_out <- list(
    updated_at = Sys.time(),
    records = env_trace$records
  )

  if (file.exists(path_file)) {
    old <- try(readRDS(path_file), TRUE)

    if (!inherits(old, "try-error") && !is.null(old$records)) {
      lst_out$records <- c(old$records, lst_out$records)
    }
  }

  saveRDS(lst_out, path_file)

  message(glue::glue("S4 trace saved: {path_file}"))

  invisible(path_file)
}

# ==========================================================================
# ==========================================================================
# ==========================================================================

.run_enter_s4_method <- function(generic,
  signature_defined,
  method_hash_raw,
  method_hash_clean,
  formal_names,
  envir)
{
  if (!.run_s4_trace_enabled()) {
    return(NA_integer_)
  }

  env_trace <- .autor_s4_trace_runtime$current
  is_root <- FALSE

  if (!is.environment(env_trace)) {
    env_trace <- .run_new_s4_trace_context()
    env_trace$root_generic <- generic
    env_trace$root_mode <- "auto"

    .autor_s4_trace_runtime$current <- env_trace
    is_root <- TRUE
  }

  env_trace$n_event <- env_trace$n_event + 1L
  id_event <- env_trace$n_event

  id_parent <- if (length(env_trace$stack)) {
    env_trace$stack[[ length(env_trace$stack) ]]
  } else {
    NA_integer_
  }

  str_signature <- .run_signature_text(signature_defined)
  str_method_key <- digest::digest(
    list(
      generic = generic,
      signature = str_signature
    ),
    algo = "xxhash64",
    serializeVersion = 3L
  )

  vec_arg_class <- .run_get_arg_classes(
    formal_names,
    envir = envir
  )

  env_trace$records[[ id_event ]] <- list(
    trace_id = env_trace$trace_id,
    event_id = id_event,
    parent_event_id = id_parent,
    depth = length(env_trace$stack) + 1L,
    is_root = is_root,
    generic = generic,
    signature_defined = signature_defined,
    signature_text = str_signature,
    method_hash_raw = method_hash_raw,
    method_hash_clean = method_hash_clean,
    method_key = str_method_key,
    formal_names = formal_names,
    arg_class = vec_arg_class,
    entered_at = Sys.time(),
    exited_at = NULL,
    status = "entered"
  )

  env_trace$stack <- c(env_trace$stack, id_event)

  id_event
}

.run_exit_s4_method <- function(id_event)
{
  if (!.run_s4_trace_enabled()) {
    return(invisible(NULL))
  }

  env_trace <- .autor_s4_trace_runtime$current

  if (!is.environment(env_trace) || is.na(id_event)) {
    return(invisible(NULL))
  }

  if (length(env_trace$records) >= id_event) {
    env_trace$records[[ id_event ]]$exited_at <- Sys.time()
    env_trace$records[[ id_event ]]$status <- "done"
  }

  if (length(env_trace$stack)) {
    env_trace$stack <- env_trace$stack[ -length(env_trace$stack) ]
  }

  if (!length(env_trace$stack) &&
    identical(env_trace$root_mode, "auto")) {
    .run_save_s4_trace_context(env_trace)
    .autor_s4_trace_runtime$current <- NULL
  }

  invisible(NULL)
}

setReplaceMethod_traceable <- function(f,
  signature = character(),
  definition,
  where = topenv(parent.frame()),
  valueClass = NULL,
  sealed = FALSE)
{
  f_replace <- paste0(f, "<-")

  setMethod_traceable(
    f = f_replace,
    signature = signature,
    definition = definition,
    where = where,
    valueClass = valueClass,
    sealed = sealed
  )
}

setMethod_traceable <- function(f,
  signature = character(),
  definition,
  where = topenv(parent.frame()),
  valueClass = NULL,
  sealed = FALSE)
{
  if (!is.function(definition)) {
    return(methods::setMethod(
      f = f,
      signature = signature,
      definition = definition,
      where = where,
      valueClass = valueClass,
      sealed = sealed
    ))
  }

  vec_formals <- names(formals(definition))

  txt_body_raw <- .run_safe_deparse(body(definition))

  body_clean <- body(definition)
  # body_clean <- .run_strip_semantic_expr(
  #   body(definition),
  #   vec_strip_fun = c(
  #     "methodAdd",
  #     "methodAdd_onExit",
  #     "snapAdd",
  #     "set_lab_legend",
  #     "setLegend",
  #     "step_message"
  #   )
  # )

  txt_body_clean <- .run_safe_deparse(body_clean)

  str_hash_raw <- digest::digest(
    list(
      generic = f,
      signature = signature,
      formals = formals(definition),
      body = txt_body_raw
    ),
    algo = "xxhash64",
    serializeVersion = 3L
  )

  str_hash_clean <- digest::digest(
    list(
      generic = f,
      signature = signature,
      formals = formals(definition),
      body = txt_body_clean
    ),
    algo = "xxhash64",
    serializeVersion = 3L
  )

  expr_enter <- substitute(
    .autor_s4_trace_event_id <- .run_enter_s4_method(
      generic = GENERIC,
      signature_defined = SIGNATURE,
      method_hash_raw = HASH_RAW,
      method_hash_clean = HASH_CLEAN,
      formal_names = FORMALS,
      envir = environment()
    ),
    list(
      GENERIC = f,
      SIGNATURE = signature,
      HASH_RAW = str_hash_raw,
      HASH_CLEAN = str_hash_clean,
      FORMALS = vec_formals
    )
  )

  expr_exit <- quote(
    on.exit(.run_exit_s4_method(.autor_s4_trace_event_id), add = TRUE)
  )

  body_old <- body(definition)

  if (is.call(body_old) && identical(body_old[[ 1L ]], as.name("{"))) {
    body_new <- as.call(c(
      list(as.name("{"), expr_enter, expr_exit),
      as.list(body_old)[ -1L ]
    ))
  } else {
    body_new <- as.call(list(
      as.name("{"),
      expr_enter,
      expr_exit,
      body_old
    ))
  }

  body(definition) <- body_new

  methods::setMethod(
    f = f,
    signature = signature,
    definition = definition,
    where = where,
    valueClass = valueClass,
    sealed = sealed
  )
}

# ==========================================================================
# ==========================================================================
# ==========================================================================

.run_replace_setMethod_text <- function(path_file,
  backup = TRUE)
{
  txt <- readLines(path_file, warn = FALSE)

  txt_new <- txt

  txt_new <- gsub(
    "(^|[^A-Za-z0-9_.:])methods::setReplaceMethod[[:space:]]*\\(",
    "\\1setReplaceMethod_traceable(",
    txt_new,
    perl = TRUE
  )

  txt_new <- gsub(
    "(^|[^A-Za-z0-9_.:])setReplaceMethod[[:space:]]*\\(",
    "\\1setReplaceMethod_traceable(",
    txt_new,
    perl = TRUE
  )

  txt_new <- gsub(
    "(^|[^A-Za-z0-9_.:])methods::setMethod[[:space:]]*\\(",
    "\\1setMethod_traceable(",
    txt_new,
    perl = TRUE
  )

  txt_new <- gsub(
    "(^|[^A-Za-z0-9_.:])setMethod[[:space:]]*\\(",
    "\\1setMethod_traceable(",
    txt_new,
    perl = TRUE
  )

  if (identical(txt, txt_new)) {
    return(invisible(FALSE))
  }

  text_parse <- paste(txt_new, collapse = "\n")
  res <- try(parse(text = text_parse), TRUE)

  if (inherits(res, "try-error")) {
    stop(glue::glue("Parse failed after replacement: {path_file}"))
  }

  if (backup) {
    file.copy(
      path_file,
      paste0(path_file, ".bak"),
      overwrite = TRUE
    )
  }

  writeLines(
    txt_new,
    con = path_file,
    useBytes = TRUE
  )

  message(glue::glue("Replaced setMethod in: {path_file}"))

  invisible(TRUE)
}

replace_setMethod_text_pkgs <- function(path_root,
  backup = TRUE, exclude = "union.utils")
{
  vec_pkg <- list.dirs(
    path_root,
    recursive = FALSE,
    full.names = TRUE
  )

  vec_pkg <- vec_pkg[
    file.exists(file.path(vec_pkg, "DESCRIPTION")) &
      dir.exists(file.path(vec_pkg, "R"))
  ]

  if (!is.null(exclude)) {
    vec_pkg <- vec_pkg[ !basename(vec_pkg) %in% exclude ]
  }

  if (!length(vec_pkg)) {
    stop(glue::glue("No package found under: {path_root}"))
  }

  vec_files <- unlist(lapply(
    vec_pkg,
    function(path_pkg) {
      list.files(
        file.path(path_pkg, "R"),
        pattern = "\\.[Rr]$",
        full.names = TRUE
      )
    }
  ))

  invisible(lapply(
    vec_files,
    .run_replace_setMethod_text,
    backup = backup
  ))
}

.compile_tools <- new.env(parent = emptyenv())

.compile_tools$is_excluded_fun <- function(name)
{
  name %in% c(
    ".run_enter_s4_method",
    ".run_exit_s4_method",
    ".run_s4_trace_enabled",
    ".run_new_s4_trace_context",
    ".run_signature_text",
    ".run_get_arg_classes",
    ".run_save_s4_trace_context",
    ".run_make_s4_trace_file",
    ".run_start_s4_trace_root",
    ".run_end_s4_trace_root",
    "setReplaceMethod_traceable",
    "setMethod_traceable"
  )
}

.compile_tools$get_pkgs <- function() {
  c(basename(list_unions()), "union.utils")
}

.compile_tools$safe_deparse <- function(x)
{
  paste(deparse(x, width.cutoff = 500L), collapse = "\n")
}

.compile_tools$get_call_name <- function(expr)
{
  if (!is.call(expr)) {
    return(NA_character_)
  }

  .compile_tools$safe_deparse(expr[[ 1L ]])
}

.compile_tools$is_ns_call <- function(name)
{
  grepl("::|:::", name)
}

.compile_tools$get_fun_from_pkgs <- function(name, vec_pkg)
{
  if (is.na(name) || !nzchar(name)) {
    return(NULL)
  }

  for (pkg in vec_pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      next
    }

    env_ns <- asNamespace(pkg)

    if (exists(name, envir = env_ns, inherits = FALSE)) {
      obj <- get(name, envir = env_ns, inherits = FALSE)

      if (is.function(obj)) {
        return(obj)
      }
    }
  }

  NULL
}

.compile_tools$is_private_fun <- function(name, vec_pkg)
{
  if (is.na(name) || !nzchar(name)) {
    return(FALSE)
  }

  if (.compile_tools$is_ns_call(name)) {
    return(FALSE)
  }

  !is.null(.compile_tools$get_fun_from_pkgs(name, vec_pkg))
}

.compile_tools$collect_private_calls_expr <- function(expr, vec_pkg)
{
  vec_out <- character(0)

  if (is.call(expr)) {
    name <- .compile_tools$get_call_name(expr)

    if (.compile_tools$is_private_fun(name, vec_pkg) &&
        !.compile_tools$is_excluded_fun(name)) {
      vec_out <- c(vec_out, name)
    }

    vec_child <- unlist(
      lapply(
        as.list(expr),
        .compile_tools$collect_private_calls_expr,
        vec_pkg = vec_pkg
      ),
      use.names = FALSE
    )

    return(unique(c(vec_out, vec_child)))
  }

  if (is.expression(expr)) {
    return(unique(unlist(
      lapply(
        as.list(expr),
        .compile_tools$collect_private_calls_expr,
        vec_pkg = vec_pkg
      ),
      use.names = FALSE
    )))
  }

  character(0)
}

.compile_tools$collect_private_deps_fun <- function(fun, vec_pkg, seen = character(0))
{
  if (!is.function(fun)) {
    return(list(seen = seen, functions = list()))
  }

  vec_calls <- .compile_tools$collect_private_calls_expr(
    body(fun),
    vec_pkg = vec_pkg
  )

  vec_calls <- setdiff(vec_calls, seen)

  if (!length(vec_calls)) {
    return(list(seen = seen, functions = list()))
  }

  lst_nested <- lapply(
    vec_calls,
    function(name) {
      fun_child <- .compile_tools$get_fun_from_pkgs(name, vec_pkg)

      if (is.null(fun_child)) {
        return(list(seen = seen, functions = list()))
      }

      message(glue::glue("Collect private function: {name}"))

      res <- .compile_tools$collect_private_deps_fun(
        fun_child,
        vec_pkg = vec_pkg,
        seen = c(seen, name)
      )

      lst_one <- list()
      lst_one[[ name ]] <- .compile_tools$clean_fun(fun_child)

      list(
        seen = res$seen,
        functions = c(lst_one, res$functions)
      )
    }
  )

  seen_new <- unique(c(
    seen,
    vec_calls,
    unlist(lapply(lst_nested, function(x) x$seen), use.names = FALSE)
  ))

  lst_functions <- do.call(
    c,
    lapply(lst_nested, function(x) x$functions)
  )

  list(
    seen = seen_new,
    functions = lst_functions
  )
}

.compile_tools$get_method_from_record <- function(record)
{
  method <- try(
    methods::selectMethod(
      f = record$generic,
      signature = record$signature_defined
    ),
    TRUE
  )

  if (inherits(method, "try-error")) {
    message(glue::glue(
      "Can not select method: {record$generic} / {record$signature_text}"
    ))
    return(NULL)
  }

  method
}

.compile_tools$collect_private_deps_records <- function(records, vec_pkg)
{
  lst_methods <- list()
  lst_functions <- list()
  seen_fun <- character(0)

  invisible(lapply(
    records,
    function(record) {
      method <- .compile_tools$get_method_from_record(record)

      if (is.null(method)) {
        return(NULL)
      }

      key <- record$method_key

      if (key %in% names(lst_methods)) {
        return(NULL)
      }

      message(glue::glue(
        "Collect S4 method: {record$generic} / {record$signature_text}"
      ))

      lst_methods[[ key ]] <<- .compile_tools$clean_fun(method)

      res <- .compile_tools$collect_private_deps_fun(
        method,
        vec_pkg = vec_pkg,
        seen = seen_fun
      )

      seen_fun <<- unique(c(seen_fun, res$seen))
      lst_functions <<- c(lst_functions, res$functions)

      NULL
    }
  ))

  list(
    methods = lst_methods,
    functions = lst_functions
  )
}

.compile_tools$strip_trace_expr <- function(expr)
{
  if (is.function(expr)) {
    body(expr) <- .compile_tools$strip_trace_expr(body(expr))
    return(expr)
  }

  if (is.call(expr)) {

    name <- .compile_tools$get_call_name(expr)

    if (identical(name, "<-") &&
        length(expr) >= 3L &&
        identical(.compile_tools$safe_deparse(expr[[ 2L ]]), ".autor_s4_trace_event_id")) {
      return(quote(NULL))
    }

    if (identical(name, "on.exit")) {
      txt <- .compile_tools$safe_deparse(expr)

      if (grepl(".run_exit_s4_method", txt, fixed = TRUE)) {
        return(quote(NULL))
      }
    }

    if (name %in% c(
      ".run_enter_s4_method",
      ".run_exit_s4_method"
    )) {
      return(quote(NULL))
    }

    lst_expr <- lapply(
      as.list(expr),
      .compile_tools$strip_trace_expr
    )

    return(as.call(lst_expr))
  }

  if (is.expression(expr)) {
    return(as.expression(lapply(
      as.list(expr),
      .compile_tools$strip_trace_expr
    )))
  }

  expr
}

.compile_tools$as_source_fun <- function(name, fun,
  prefix = ".compiled",
  always_prefix = TRUE)
{
  name <- .compile_tools$make_safe_name(
    name,
    prefix = prefix,
    always_prefix = always_prefix
  )

  fun <- .compile_tools$clean_fun(fun)
  txt <- .compile_tools$safe_deparse(fun)

  paste0(name, " <- ", txt, "\n")
}

.compile_tools$export_source_one <- function(res,
  file = "tmp/compiled_methods.R")
{
  dir.create(
    dirname(file),
    FALSE,
    recursive = TRUE
  )

  txt_helpers <- unlist(lapply(
    names(res$functions),
    function(name) {
      .compile_tools$as_source_fun(
        name,
        res$functions[[ name ]],
        prefix = ".helper",
        always_prefix = TRUE
      )
    }
  ))

  txt_aliases <- .compile_tools$export_helper_aliases(res)

  txt_methods <- unlist(lapply(
    names(res$methods),
    function(name) {
      .compile_tools$as_source_fun(
        name,
        res$methods[[ name ]],
        prefix = ".method",
        always_prefix = TRUE
      )
    }
  ))

  txt <- c(
    "# ============================================================",
    "# Helpers",
    "# ============================================================",
    "",
    txt_helpers,
    "",
    txt_aliases,
    "",
    "# ============================================================",
    "# S4 methods collected from trace",
    "# ============================================================",
    "",
    txt_methods
  )

  writeLines(
    txt,
    con = file,
    useBytes = TRUE
  )

  res_parse <- try(parse(file), TRUE)

  if (inherits(res_parse, "try-error")) {
    message(glue::glue("Parse failed: {file}"))
  } else {
    message(glue::glue("Export compiled source: {file}"))
  }

  invisible(file)
}

.compile_tools$make_safe_name <- function(name,
  prefix = ".compiled",
  always_prefix = TRUE)
{
  name <- as.character(name)
  name <- gsub("[^A-Za-z0-9_.]", "_", name)

  if (always_prefix) {
    name <- paste0(prefix, "_", name)
  }

  if (!grepl("^[A-Za-z.]", name)) {
    name <- paste0(prefix, "_", name)
  }

  if (grepl("^\\.[0-9]", name)) {
    name <- paste0(prefix, "_", sub("^\\.", "", name))
  }

  name
}

.compile_tools$clean_fun <- function(fun)
{
  if (!is.function(fun)) {
    return(fun)
  }

  fun <- .compile_tools$strip_trace_expr(fun)

  fun <- .run_strip_semantic_expr(
    fun,
    vec_strip_fun = c(
      "methodAdd",
      "methodAdd_onExit",
      "snapAdd",
      "set_lab_legend",
      "setLegend",
      "step_message"
    )
  )

  fun
}

.compile_tools$collect_private_deps_trace <- function(path_trace, vec_pkg)
{
  trace <- readRDS(path_trace)

  if (is.null(trace$records)) {
    stop("No records in trace file.")
  }

  res <- .compile_tools$collect_private_deps_records(
    trace$records,
    vec_pkg = vec_pkg
  )

  list(
    records = trace$records,
    methods = res$methods,
    functions = res$functions
  )
}

.compile_tools$get_helper_alias <- function(name)
{
  safe <- .compile_tools$make_safe_name(
    name,
    prefix = ".helper",
    always_prefix = TRUE
  )

  glue::glue("{name} <- {safe}")
}

.compile_tools$export_helper_aliases <- function(res)
{
  if (is.null(res$functions) || !length(res$functions)) {
    return(character(0))
  }

  aliases <- vapply(
    names(res$functions),
    FUN.VALUE = character(1),
    .compile_tools$get_helper_alias
  )

  c(
    "# ---- helper aliases ----",
    aliases,
    ""
  )
}

.compile_tools$method_fun_name <- function(method_key)
{
  .compile_tools$make_safe_name(
    method_key,
    prefix = ".method",
    always_prefix = TRUE
  )
}

.compile_tools$as_records_data <- function(records)
{
  if (!length(records)) {
    return(data.frame())
  }

  data_records <- do.call(
    rbind,
    lapply(
      seq_along(records),
      function(n) {
        x <- records[[ n ]]

        trace_id <- if (!is.null(x$trace_id)) {
          x$trace_id
        } else {
          paste0("legacy_", n)
        }

        event_uid <- paste0(trace_id, "::", x$event_id)

        parent_event_uid <- if (is.na(x$parent_event_id)) {
          NA_character_
        } else {
          paste0(trace_id, "::", x$parent_event_id)
        }

        data.frame(
          record_id = n,
          trace_id = trace_id,
          event_uid = event_uid,
          parent_event_uid = parent_event_uid,
          event_id = x$event_id,
          parent_event_id = x$parent_event_id,
          depth = x$depth,
          generic = x$generic,
          signature_text = x$signature_text,
          method_key = x$method_key,
          method_fun = .compile_tools$method_fun_name(x$method_key),
          method_hash_raw = x$method_hash_raw,
          method_hash_clean = x$method_hash_clean,
          status = x$status,
          stringsAsFactors = FALSE
        )
      }
    )
  )

  data_records
}

.compile_tools$as_method_index <- function(records)
{
  data_records <- .compile_tools$as_records_data(records)

  if (!nrow(data_records)) {
    return(data.frame())
  }

  data_index <- unique(data_records[ , c(
    "generic",
    "signature_text",
    "method_key",
    "method_fun",
    "method_hash_raw",
    "method_hash_clean"
  ) ])

  rownames(data_index) <- NULL

  data_index
}

.compile_tools$as_method_edges <- function(records)
{
  data_records <- .compile_tools$as_records_data(records)

  if (!nrow(data_records)) {
    return(data.frame())
  }

  data_parent <- data_records[ , c(
    "event_uid",
    "method_key",
    "method_fun"
  ) ]

  colnames(data_parent) <- c(
    "parent_event_uid",
    "parent_method_key",
    "parent_method_fun"
  )

  data_edges <- merge(
    data_records,
    data_parent,
    by = "parent_event_uid",
    all.x = FALSE,
    all.y = FALSE
  )

  data_edges <- data_edges[ , c(
    "trace_id",
    "event_uid",
    "parent_event_uid",
    "event_id",
    "parent_event_id",
    "parent_method_key",
    "parent_method_fun",
    "generic",
    "signature_text",
    "method_key",
    "method_fun"
  ) ]

  rownames(data_edges) <- NULL

  data_edges
}

.compile_tools$export_method_index <- function(res,
  file_index = "tmp/method_index.csv",
  file_edges = "tmp/method_edges.csv")
{
  dir.create(
    dirname(file_index),
    FALSE,
    recursive = TRUE
  )

  data_index <- .compile_tools$as_method_index(res$records)
  data_edges <- .compile_tools$as_method_edges(res$records)

  utils::write.csv(
    data_index,
    file = file_index,
    row.names = FALSE
  )

  utils::write.csv(
    data_edges,
    file = file_edges,
    row.names = FALSE
  )

  message(glue::glue("Export method index: {file_index}"))
  message(glue::glue("Export method edges: {file_edges}"))

  invisible(list(
    index = data_index,
    edges = data_edges,
    file_index = file_index,
    file_edges = file_edges
  ))
}

.compile_tools$get_child_methods <- function(idx,
  parent_method_key)
{
  if (is.null(idx$edges) || !nrow(idx$edges)) {
    return(character(0))
  }

  vec_child <- idx$edges$method_fun[
    idx$edges$parent_method_key == parent_method_key
  ]

  unique(vec_child)
}

.compile_tools$get_first_formal_name <- function(fun)
{
  names(formals(fun))[[ 1L ]]
}

.compile_tools$replace_callNextMethod_expr <- function(expr,
  child_methods,
  fallback_name = NULL)
{
  if (is.call(expr)) {

    name <- .compile_tools$get_call_name(expr)

    if (identical(name, "callNextMethod")) {

      if (!length(child_methods)) {
        if (!is.null(fallback_name) && nzchar(fallback_name)) {
          return(as.name(fallback_name))
        }

        return(quote(NULL))
      }

      if (length(child_methods) == 1L) {
        return(as.call(c(
          list(as.name(child_methods[[ 1L ]])),
          as.list(expr)[ -1L ]
        )))
      }

      return(as.call(c(
        list(as.name("list")),
        lapply(
          child_methods,
          function(fun) {
            as.call(c(
              list(as.name(fun)),
              as.list(expr)[ -1L ]
            ))
          }
        )
      )))
    }

    lst_expr <- lapply(
      as.list(expr),
      .compile_tools$replace_callNextMethod_expr,
      child_methods = child_methods,
      fallback_name = fallback_name
    )

    return(as.call(lst_expr))
  }

  if (is.expression(expr)) {
    return(as.expression(lapply(
      as.list(expr),
      .compile_tools$replace_callNextMethod_expr,
      child_methods = child_methods,
      fallback_name = fallback_name
    )))
  }

  expr
}

.compile_tools$resolve_callNextMethod_fun <- function(fun,
  method_key,
  idx)
{
  if (!is.function(fun)) {
    return(fun)
  }

  child_methods <- .compile_tools$get_child_methods(
    idx,
    parent_method_key = method_key
  )

  fallback_name <- .compile_tools$get_first_formal_name(fun)

  if (!length(child_methods) &&
      grepl("callNextMethod", .compile_tools$safe_deparse(body(fun)), fixed = TRUE)) {
    message(glue::glue(
      "Replace orphan callNextMethod with `{fallback_name}` in method: {method_key}"
    ))
  }

  body(fun) <- .compile_tools$replace_callNextMethod_expr(
    body(fun),
    child_methods = child_methods,
    fallback_name = fallback_name
  )

  fun
}

.compile_tools$resolve_callNextMethod_methods <- function(res,
  idx)
{
  if (is.null(res$methods) || !length(res$methods)) {
    return(res)
  }

  res$methods <- setNames(
    lapply(
      names(res$methods),
      function(key) {
        .compile_tools$resolve_callNextMethod_fun(
          res$methods[[ key ]],
          method_key = key,
          idx = idx
        )
      }
    ),
    names(res$methods)
  )

  res
}


