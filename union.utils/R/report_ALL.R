# ==========================================================================
# global config
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

create_remote_project <- function(project = guess_project(), ws = getRemoteWs(), 
  remote = "remote")
{
  cmd <- glue::glue("cd {ws} && mkdir {project}")
  cdRun("ssh ", remote, " '", cmd, "'")
}

getRemoteWs <- function() {
  path <- getOption("remote_working_space")
  if (is.null(path)) {
    stop('is.null(path).')
  }
  path
}

guess_project <- function(path = getwd()) {
  res <- stringr::str_extract(basename(path), "[0-9]+_[a-zA-Z]+[0-9]+")
  if (is.na(res)) {
    stop('is.na(res).')
  }
  res
}

send_job_to_remote <- function(path, to = "rds_jobSave", remote = "remote") {
  if (!is_sshfs_mount(remote)) {
    stop('!is_sshfs_mount(remote).')
  }
  if (!file.exists(path)) {
    stop('!file.exists(path).')
  }
  file.copy(path, file.path(remote, to), TRUE)
}


mark_text <- function(string, color, bold = TRUE, ...) {
  string <- gs(string, "&lt;", "<")
  string <- gs(string, "&gt;", ">")
  ftext <- officer::ftext(
    string, officer::fp_text_lite(color = color, bold = bold, ...)
  )
  paste0("`", officer::to_wml(ftext), "`{=openxml}")
}

mark <- list()

mark$red <- function(string) {
  mark_text(string, color = "#C00000")
}

mark$sig <- mark$red

mark$blue <- function(string) {
  mark_text(string, color = "#2E75B5")
}

mark$th <- mark$blue

mark$green <- function(string) {
  mark_text(string, color = "green")
}


get_file_with_format_name <- function(file, name) {
  filename <- paste0(name, ".", tools::file_ext(file))
  file_new <- file.path(dirname(file), filename)
  file.copy(file, file_new, TRUE)
  url <- glue::glue("file://{normalizePath(file_new)}")
  if (nchar(Sys.which("wl-copy"))) {
    system(glue::glue("echo -n {url} | wl-copy -t text/uri-list"))
  } else {
    stop('nchar(Sys.which("wl-copy")).')
  }
}

gett_file <- function(url) {
  url <- glue::glue("file:/{normalizePath(url)}")
  if (nchar(Sys.which("wl-copy"))) {
    system(glue::glue("echo -n {url} | wl-copy -t text/uri-list"))
  } else {
    stop('nchar(Sys.which("wl-copy")).')
  }
}

gett_files <- function(files) {
  if (!nchar(Sys.which("wl-copy"))) {
    stop("wl-copy not found.", call. = FALSE)
  }

  files <- normalizePath(files, mustWork = TRUE)

  uri <- paste0(
    "file://",
    files,
    collapse = "\n"
  )

  tf <- tempfile(fileext = ".txt")
  writeLines(uri, tf)

  system2(
    "wl-copy",
    c("--type", "text/uri-list"),
    stdin = tf
  )

  invisible(uri)
}

get_contents_refered_from_fields <- function(
  file, ids = "foreword", id_ref = "reference", save_bib = "library.bib",
  to_clipboard = TRUE, lines = NULL)
{
  if (is.null(lines)) {
    lines <- readLines(file)
  }
  fields <- detect_field(lines, c(ids, id_ref))
  if (any(lengths(fields) < 1)) {
    stop('any(lengths(fields) < 1).')
  }
  fun_format <- function(x) paste0(unlist(x), collapse = "\n")
  res <- lapply(ids,
    function(id) {
      res <- parse_references_from_text(
        fun_format(fields[[ id ]]), fun_format(fields[[ id_ref ]])
      )
      indices <- unlist(res$mapping$indices)
      res$reference <- res$reference[ res$reference != "" ]
      if (!all(indices %in% seq_along(res$reference))) {
        # Terror <<- namel(indices, res)
        stop('!all(indices %in% seq_along(res$reference)), not match reference.')
      }
      pmids <- strx(res$reference, "(?<=PMID:[\\s\n\\d]{0,1})[0-9]+")
      if (any(is.na(pmids))) {
        stop('any(is.na(pmids)), some reference do not have pmid, please check manualy')
      }
      bibs <- expect_local_data(
        "tmp", "bib_pmid", get_bibs_by_pmid, list(pmids)
      )
      refs <- .refs(names(bibs))
      list(bibs = bibs, content = glue::glue(res$content))
    })
  contents <- vapply(res, function(x) x$content, character(1))
  if (to_clipboard) {
    gett(contents)
  }
  if (!is.null(save_bib)) {
    bibs <- do.call(c, lapply(res, function(x) x$bibs))
    RefManageR::WriteBib(bibs, save_bib)
  }
}

.run_strip_semantic_layer <- function(path_root,
  vec_strip_fun = c(
    "methodAdd",
    "methodAdd_onExit",
    "snapAdd",
    "set_lab_legend"
  ),
  overwrite = FALSE)
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

  if (!length(vec_pkg)) {
    stop(glue::glue(
      "No package found under: {path_root}"
    ))
  }

  message(glue::glue(
    "Found {length(vec_pkg)} packages."
  ))

  invisible(
    lapply(vec_pkg,
      function(path_pkg) {

        message(glue::glue(
          "Processing package: {basename(path_pkg)}"
        ))

        vec_files <- list.files(
          file.path(path_pkg, "R"),
          pattern = "\\.[Rr]$",
          full.names = TRUE
        )

        invisible(
          lapply(vec_files,
            function(path_file) {

              message(glue::glue(
                "  Processing file: {basename(path_file)}"
              ))

              exprs <- parse(
                path_file,
                keep.source = TRUE
              )

              exprs_new <- lapply(
                as.list(exprs),
                .run_strip_semantic_expr,
                vec_strip_fun = vec_strip_fun
              )

              txt <- vapply(
                exprs_new,
                FUN.VALUE = character(1),
                function(x) {
                  paste(
                    deparse(
                      x,
                      width.cutoff = 500L
                    ),
                    collapse = "\n"
                  )
                }
              )

              if (overwrite) {
                path_out <- path_file
              } else {
                path_out <- paste0(
                  tools::file_path_sans_ext(path_file),
                  ".release.R"
                )
              }

              writeLines(
                txt,
                con = path_out,
                useBytes = TRUE
              )

              message(glue::glue(
                "  Write: {basename(path_out)}"
              ))

              invisible(path_out)
            })
        )
      })
  )
}

.run_strip_semantic_expr <- function(expr,
  vec_strip_fun = c(
    "methodAdd",
    "methodAdd_onExit",
    "snapAdd",
    "set_lab_legend"
  ))
{
  # function closure
  if (is.function(expr)) {

    body(expr) <- .run_strip_semantic_expr(
      body(expr),
      vec_strip_fun = vec_strip_fun
    )

    return(expr)
  }

  # language object
  if (is.call(expr)) {

    fun <- expr[[ 1L ]]

    str_fun <- paste(
      deparse(fun, width.cutoff = 500L),
      collapse = ""
    )

    # strip middleware first
    if (str_fun %in% vec_strip_fun) {

      if (identical(str_fun, "methodAdd") ||
          identical(str_fun, "snapAdd")) {

        if (length(expr) >= 2L) {

          return(
            .run_strip_semantic_expr(
              expr[[ 2L ]],
              vec_strip_fun = vec_strip_fun
            )
          )
        }

        return(quote(NULL))
      }

      if (identical(str_fun, "methodAdd_onExit")) {
        return(quote(invisible(NULL)))
      }

      if (identical(str_fun, "set_lab_legend")) {
        lst_expr <- as.list(expr)

        if (length(lst_expr) >= 4L) {
          lst_expr[ 4L:length(lst_expr) ] <- ""
        }

        return(as.call(lst_expr))
      }
    }

    # strip semantic glue
    if (str_fun %in% c(
      "glue",
      "glue::glue"
    )) {

      txt_expr <- paste(
        deparse(expr, width.cutoff = 500L),
        collapse = "\n"
      )

      if (.is_semantic_text(txt_expr)) {

        message(glue::glue(
          "Strip semantic glue text: {substr(txt_expr, 1L, 80L)}..."
        ))

        return("")
      }
    }

    # recursive rewrite
    lst_expr <- lapply(
      as.list(expr),
      .run_strip_semantic_expr,
      vec_strip_fun = vec_strip_fun
    )

    return(as.call(lst_expr))
  }

  # expression vector
  if (is.expression(expr)) {

    return(as.expression(
      lapply(
        as.list(expr),
        .run_strip_semantic_expr,
        vec_strip_fun = vec_strip_fun
      )
    ))
  }

  # strip semantic character constant
  if (is.character(expr)) {
    txt <- paste0(expr, collapse = "\n")
    if (.is_semantic_text(txt)) {
      message(glue::glue(
          "Strip semantic character: {substr(txt, 1L, 80L)}..."
          ))
      return("")
    }
  }
  expr
}

.is_semantic_text <- function(txt)
{
  if (!nzchar(txt)) {
    return(FALSE)
  }

  # no chinese
  if (!grepl("[\u4e00-\u9fff]", txt, perl = TRUE)) {
    return(FALSE)
  }

  # likely label
  if (nchar(txt) <= 20L &&
      !grepl("[。；：，]", txt)) {
    return(FALSE)
  }

  TRUE
}

.set_gwas_token <- function() {
  token <- "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJzaGFtYW4teWVsbG93QGZveG1haWwuY29tIiwiaWF0IjoxNzc5NTUyMzI0LCJleHAiOjE3ODA3NjE5MjR9.oXetTruXv_kmBjCTc1nkWlHfE6MmdBLShmssblJ-wXzNjdr7hOCI_WDtue6flz3ZySfsXVkqZybod6zcRkap2GFikXRJKZKXjMViVdIaNr0xYoaKghV9dzHSULoiOwy9m-6gIqBozLisVTyjRZSNANOPsBwXOjEp4oS3lHMNnDctAjF66poY0555nSiIa4djFH7FWj7Pml-YuEk2heYqT7_hVr_JNItWLX9IJr8FkrZ_T2snWnyVXvAwCA_KKP8o8TWfmgpDeUhFegSjeVvvidK_QxVjg1j5arOkXfFXHGuEIhGowz2rP04XgJmE6i5IkVxbW48A9qZlquV3hL_xKw"
  expire <- as.Date("2026-06-06")
  if (Sys.Date() >= expire) {
    message(glue::glue("GWAS API token expired, please reset it: <https://api.opengwas.io/profile/>"))
    NULL
  } else {
    token
  }
}

clear_feature <- function(x, name = "key.rds", dir = ".", 
  file = file.path(dir, name))
{
  if (!is(x, "feature")) {
    stop('!is(x, "feature").')
  }
  saveRDS(x, file)
}

load_feature <- function(
  name = "key.rds", dir = ".", 
  file = file.path(dir, name), 
  job_otherwise = file.path("rds_jobSave", "ven.marker.0.rds"), remote = NULL)
{
  if (file.exists(file)) {
    readRDS(file)
  } else {
    if (!is.null(remote)) {
      job_otherwise <- file.path(remote, job_otherwise)
    }
    feature(readRDS(job_otherwise))
  }
}

expect_package <- function(pkg, version, prio_lib = getOption("prio_lib")) {
  if (!requireNamespace(pkg)) {
    stop('!requireNamespace(pkg)')
  }
  if (packageVersion(pkg) >= version) {
    message("Pacakge ", pkg, " as expected.")
    return()
  }
  if (packageVersion(pkg) < version) {
    message(glue::glue("Detach the loaded namespace, search in preferred lib path."))
    unloadNamespace(asNamespace(pkg))
    loadNamespace(pkg, lib.loc = prio_lib)
  }
  if (packageVersion(pkg, lib.loc = prio_lib) < version) {
    stop('packageVersion(pkg, lib.loc = prio_lib) < version.')
  } else {
    message(glue::glue("Successfully loaded the latter R package"))
  }
}

spsv <- function(object, name = NULL, prefix = "tmp") {
  if (is.null(name)) {
    name <- formal_name(rlang::expr_text(substitute(object)))
  }
  fun <- select_savefun(object)
  fun(object, name = name, mkdir = prefix)
}

smart_wrap_expr <- function(plots, size = 3, ..., envir = .GlobalEnv)
{
  calls <- substitute(plots)
  if (as_label(calls[[1]]) != "{") {
    stop('as_label(calls[[1]]) != "{"')
  }
  plots <- lapply(calls[-1],
    function(call) {
      eval(parse(text = as_label(call)), envir = envir)
    })
  smart_wrap(plots, size = size, ...)
}


convert_pdf_in_project <- function(path = "remote", skip = NULL)
{
  dirs <- gs(
    list.files(path, "^r\\.[0-9]+.*\\.r$", full.names = TRUE), 
    "(?<=/)r\\.|\\.r$", "", perl = TRUE
  )
  order <- as.integer(strx(basename(dirs), "[0-9]+"))
  dirs <- dirs[ order(order) ]
  if (!is.null(skip)) {
    message(glue::glue("Skip: \n{bind(dirs[skip], co = '\n')}"))
    dirs <- dirs[-skip]
  }
  targets <- list.files(dirs, "\\.pdf$", full.names = TRUE)
  pbapply::pblapply(targets,
    function(file) {
      newfile <- paste0(tools::file_path_sans_ext(file), ".png")
      res <- try(pdf_convert(file, filenames = newfile, dpi = 300, pages = 1))
      if (inherits(res, "try-error")) {
        sink()
        message(glue::glue("Failed to convert file: {file}"))
      }
    })
}


setup_counting_in_directory <- function(dir, pattern = "^[0-9]+") {
  unlink(
    list.files(dir, pattern, full.names = TRUE), force = TRUE, recursive = TRUE
  )
  options(
    autor_counting_start_dir = dir,
    savedir = list(figs = dir, tabs = dir)
  )
}

output_with_counting_number <- function(plots, envir = .GlobalEnv, 
  fun_wrap = "autor", extra_cmd = NULL)
{
  if (is.null(output <- getOption("autor_counting_start_dir"))) {
    stop('is.null(getOption("autor_counting_start_dir")).')
  }
  if (!dir.exists(output)) {
    stop('!dir.exists(output).')
  }
  calls <- substitute(plots)
  if (!is(calls, "{")) {
    stop('!is(calls, "{").')
  }
  rapp_find_job_name <- function(x) {
    if (is(x, "call") || is(x, "{")) {
      rapp_find_job_name(x[[2]])
    } else if (is(x, "name")) {
      rlang::expr_text(x)
    } else {
      stop("The finally of the 'substitute' is: ", class(x))
    }
  }
  num <- as.integer(guess_number.hb(output, p.pattern = "^[0-9]{2}"))
  fun_num <- function(n) {
    sprintf("%02d", n)
  }
  lapply(calls[-1], 
    function(call) {
      name <- rapp_find_job_name(call)
      job <- try(get(name, envir = .GlobalEnv), TRUE)
      if (inherits(job, "try-error")) {
        .try_loadJob(name, FALSE)
      }
      object <- eval(parse(text = rlang::expr_text(call)))
      outputName <- paste0(fun_num(num), "_", label(object))
      message(glue::glue("Save as: {outputName}"))
      fun_save <- select_savefun(object)
      fun_save(object, name = outputName, mkdir = output)
      expect_file <- file.path(output, outputName)
      if (file.exists(file_pdf <- paste0(expect_file, ".pdf"))) {
        res <- try(
          pdf_convert(file_pdf, filenames = paste0(expect_file, ".png"), dpi = 300, pages = 1)
        )
        if (inherits(res, "try-error")) {
          sink()
          message(glue::glue("Failed to convert file: {file}"))
        }
      }
      num <<- num + 1L
    })
}

methodDefinition_as_setMethod_call <- function(m) {
  stopifnot(methods::is(m, "MethodDefinition"))
  generic <- as.character(m@generic)
  sig <- as.list(m@defined)
  fun <- m@.Data
  call("setMethod", f = generic, signature = unlist(sig), definition = fun)
}

.get_method_defination_in_package <- function(
  lang, env_class, pkgs, envs_search = lapply(pkgs, asNamespace),
  fun_getClass = .guess_class_from_lang, strip = TRUE,
  vec_strip_fun = c("methodAdd", "methodAdd_onExit", "snapAdd", "set_lab_legend")
)
{
  if (!is.language(lang)) {
    stop('!is.language(lang).')
  }
  if (!is(env_class, "environment")) {
    stop('!is(env_class, "environment").')
  }
  if (is(lang, "<-")) {
    objName <- rlang::expr_text(lang[[2]])
    if (is.null(env_class[[objName]])) {
      class <- fun_getClass(lang[[3]], env_class = env_class)
      if (!is.null(class)) {
        env_class[[objName]] <- class
      } else {
        if (grpl(objName, "^metadata")) {
          env_class[[objName]] <- "data.frame"
        } else if (grpl(objName, "^fea")) {
          env_class[[objName]] <- "feature"
        }
      }
    }
    funCall <- lang[[3]]
  } else if (is(lang, "call")) {
    funCall <- lang
  } else if (is(lang, "name")) {
    message(glue::glue("Skip: {rlang::expr_text(lang)}, is a name."))
    return()
  } else if (is(lang, "if")) {
    message(glue::glue("Skip: {rlang::expr_text(lang)}, is a 'if'."))
    return()
  } else if (is(lang, "for")) {
    message(glue::glue("Skip: {rlang::expr_text(lang)}, is a 'for'."))
    return()
  }
  if (!is.call(funCall)) {
    message(glue::glue("Skip: {rlang::expr_text(funCall)}, not a call."))
    return()
  }
  callName <- funCall[[1]]
  if (length(callName) > 1) {
    if (!any(rlang::expr_text(callName[[2]]) == pkgs)) {
      message(
        glue::glue("Skip: {rlang::expr_text(callName)}, other package.")
      )
      return()
    } else {
      fun_name <- rlang::expr_text(callName[[3]])
    }
  } else {
    fun_name <- rlang::expr_text(callName)
  }
  hasThat <- vapply(envs_search, FUN.VALUE = logical(1),
    function(env) {
      exists(fun_name, envir = env, inherits = FALSE)
    })
  if (!any(hasThat)) {
    message(glue::glue("Skip: {fun_name}, not exists."))
    return()
  }
  fun <- get_fun(fun_name, envir = envs_search[ hasThat ][[1]])
  if (isS4(fun)) {
    if (!is(fun, "genericFunction")) {
      message(glue::glue("Skip: {fun_name}, is S4, but not genericFunction."))
      return()
    }
    res <- .expr_resolve_S4(funCall, env_class = env_class)
    f <- try(selectMethod(res$callArgs$fun, signature = res$signature))
    if (inherits(f, "try-error")) {
      stop(glue::glue("`{rlang::expr_text(funCall)}`: Can not found method `{res$fname}` for signature ..."))
    }
    mcall <- methodDefinition_as_setMethod_call(f)
    if (strip) {
      mcall <- .run_strip_semantic_expr(mcall, vec_strip_fun)
    }
    text <- deparse(mcall)
    return(list(text = text, lang = mcall))
  } else {
    if (strip) {
      fun <- .run_strip_semantic_expr(fun, vec_strip_fun)
    }
    text <- deparse(fun)
    text[1] <- paste0(fun_name, " <- ", text[1])
    return(list(text = text, lang = fun))
  }
}

.guess_class_from_lang <- function(lang, pattern = "job_[a-zA-Z0-9_]+", env_class = NULL)
{
  if (is(lang, "name")) {
    from <- rlang::expr_text(lang)
    if (!is.null(env_class[[from]])) {
      return(env_class[[ from ]])
    } else {
      return(NULL)
    }
  }
  code <- rlang::expr_text(lang)
  class <- strx(code, pattern)
  fun_matchClass <- function(string, type = "rds") {
    name <- strx(string, glue::glue("(?<=rds_jobSave/).*(?=.[0-9].{type})"))
    if (!is.na(name) && !is.null(env_class[[ name ]])) {
      env_class[[ name ]]
    } else {
      "job"
    }
  }
  if (isClass(class)) {
    class
  } else {
    if (grpl(code, "copy_job") && identical(lang[[1]], as.name("copy_job"))) {
      fromJob <- rlang::expr_text(lang[[2]])
      if (!is.null(env_class[[fromJob]])) {
        env_class[[ fromJob ]]
      } else {
        NULL
      }
    } else if (grpl(nameCall <- rlang::expr_text(lang[[1]]), "^do_")) {
      class <- glue::glue("job_{s(nameCall, 'do_', '')}")
      if (isClass(class)) {
        class
      } else {
        rlang::abort(glue::glue("Guess class by `do_*` failed: {class}"))
      }
    } else if (grpl(code, "dplyr::")) {
      if (grpl(code, "^dplyr::recode")) {
        "character"
      } else {
        "data.frame"
      }
    } else if (grpl(code, "readRDS.*rds_jobSave")) {
      fun_matchClass(code, "rds")
    } else if (grpl(code, "qs::qread.*rds_jobSave")) {
      fun_matchClass(code, "qs")
    } else {
      toolSub <- c("getsub")
      for (i in toolSub) {
        if (grpl(code, i) && identical(lang[[1]], as.name(i))) {
          fromJob <- rlang::expr_text(lang[[2]])
          if (!is.null(env_class[[fromJob]])) {
            return(env_class[[ fromJob ]])
          }
        }
      }
      NULL
    }
  }
}

.check_bin <- function(x) {
  if (!nchar(Sys.which(x))) {
    stop("Command not found: ", x, call. = FALSE)
  }
}

.wait_wps_doc_ready <- function(file, timeout = 60L) {
  base <- basename(file)

  for (i in seq_len(timeout * 10L)) {
    clients <- suppressWarnings(system2("hyprctl", "clients", stdout = TRUE))
    if (any(grepl(paste0("title: ", base), clients, fixed = TRUE))) {
      system2(
        "hyprctl", c("dispatch", "focuswindow", "class:wps"),
        stdout = FALSE, stderr = FALSE
      )
      Sys.sleep(1)
      return(TRUE)
    }
    Sys.sleep(0.1)
  }

  FALSE
}

wps_pdf <- function(file, timeout = 60L) {
  .check_bin("flatpak")
  .check_bin("hyprctl")
  .check_bin("ydotool")

  file <- normalizePath(file, mustWork = TRUE)

  system2("flatpak", c("run", "com.wps.Office", shQuote(file)), wait = FALSE)

  if (!.wait_wps_doc_ready(file, timeout)) {
    stop("WPS document not ready: ", file, call. = FALSE)
  }

  system2("ydotool", c("key", "56:1", "33:1", "33:0", "56:0"))
  Sys.sleep(0.3)
  system2("ydotool", c("key", "33:1", "33:0"))

  invisible(TRUE)
}

