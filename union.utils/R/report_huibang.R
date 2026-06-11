# ==========================================================================
# huibang
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

create_remote_project.hb <- function(project = guess_project(), ws = getRemoteWs(), 
  remote = "remote")
{
  cmd <- glue::glue("cd {ws} && mkdir {project}")
  cdRun("ssh ", remote, " '", cmd, "'")
}

setup.sshfs <- function(sync = FALSE, project = guess_project(), ws = getRemoteWs(), 
  remote = "remote", path = "remote", mirror = "scripts_mirror")
{
  if (!dir.exists(mirror)) {
    stop('!dir.exists(mirror).')
  }
  if (!dir.exists(path)) {
   stop('!dir.exists("path").')
  }
  if (sync) {
    file_sync <- file.path(.expath, "scripts", "sync.sh")
    cmd <- glue::glue("sh -c 'nohup bash {file_sync} {mirror} {remote}:{ws}/{project} > sync.log 2>&1 &'")
    system(cmd, wait = FALSE)
  }
  # system2("bash", c("-c", cmd), wait = FALSE)
  if (!is_sshfs_mount(glue::glue("../{path}"))) {
    cdRun(glue::glue("nohup sshfs {remote}:{ws} ../{path} >/dev/null 2>&1 &"))
  }
  if (is_sshfs_mount(path)) {
    return(message("The directory has been mount."))
  }
  if (length(list.files(path, include.dirs = TRUE))) {
    stop('length(list.files(path, all.files = TRUE, include.dirs = TRUE)).')
  }
  # umount remote
  cdRun(glue::glue("nohup sshfs {remote}:{ws}/{project} {path} >/dev/null 2>&1 &"))
  repeat {
    Sys.sleep(1)
    if (is_sshfs_mount(path)) {
      return(TRUE)
    }
  }
}

call_nvim_for_remote_setup <- function(file_script, 
  project = guess_project(), mode = "huibang")
{
  if (mode == "huibang") {
    text <- generate_setup_codes.huibang(
      file_script, project, TRUE, path_pkg = "./tmp/f3256d7e/union/union.utils"
    )
  }
  SendCmdToNvim_lua(
    glue::glue(
      "SendCodeToTB([[{text}]])"
    )
  )
}

run_setup_codes.huibang <- function(file_script, project, envir = .GlobalEnv) {
  eval(
    generate_setup_codes.huibang(file_script, project), envir = envir
  )
}

generate_setup_codes.huibang <- function(file_script,
  project, get_text = FALSE, path_pkg = "./union/union.utils", name = NULL
)
{
  file_script <- basename(file_script)
  if (!grpl(file_script, "^r\\.[0-9]+_.*\\.r$")) {
    stop('!grpl(file_script, "^r\\.[0-9]+_.*\\.r$"), not match correct script file name.')
  }
  if (is.null(name)) {
    name <- gs(file_script, "^r\\.|\\.r$", "")
  }
  odir <- paste0("/data/nas1/huanglichuang_OD/project/", project)
  fun_replace <- function(x, envir = parent.frame(1L)) {
    glue::glue(x, .open = ".{{{", .close = "}}}.", envir = envir)
  }
  lang <- quote({
    rm(list = ls()); gc()
    ORIGINAL_DIR <- ".{{{odir}}}."
    output <- file.path(ORIGINAL_DIR, ".{{{name}}}.")
    if (!dir.exists(output)) {
      dir.create(output, recursive = TRUE)
    }
    setwd(ORIGINAL_DIR)

    .libPaths(c('/data/nas2/software/miniconda3/envs/public_R/lib/R/library/', '/data/nas1/huanglichuang_OD/conda/envs/extra_pkgs/lib/R/library/'))

    myPkg <- ".{{{path_pkg}}}."
    if (!dir.exists(myPkg)) {
      stop('Can not found package: ', myPkg)
    }
    devtools::load_all(myPkg)
    load_unions()
    setup.huibang()
  })
  text <- paste0(deparse(lang), collapse = "\n")
  text <- fun_replace(text)
  if (get_text) {
    text
  } else {
    parse(text = text)
  }
}

.upd_pkg_to_remote <- function(...) {
  .send_pkg_to_remote(..., upd = TRUE)
}

table_qc.hb <- function() {
  ftibble(file.path(.expath, "report_qc.csv"))
}

.send_pkg_to_remote <- function(from = "~/union",
  exclude = ".git", to = "remote", upd = FALSE, remoteUntar = TRUE,
  dir_relative = "tmp/f3256d7e", remote = "remote")
{
  if (!is_sshfs_mount(to)) {
    stop('!is_sshfs_mount(to).')
  }
  if (!dir.exists(toDir <- file.path(to, dir_relative))) {
    dir.create(toDir, recursive = TRUE)
  }
  pkg <- basename(from)
  archive_package <- paste0(pkg, ".tar.gz")
  if (!upd && file.exists(file.path(toDir, archive_package))) {
    stop(
      '!upd && file.exists(file.path(to, dir_relative, archive_package))'
    )
  }
  pathPkgFrom <- file.path(from, archive_package)
  # if (file.exists(pathPkgFrom)) {
  #   message(glue::glue('file.exists(pathPkgFrom), remove ...'))
  #   file.remove(pathPkgFrom)
  # }
  cdRun(
    "git ls-files -c -o --exclude-standard -z | xargs -0 tar -czf ", archive_package,
    path = from
  )
  file.copy(
    pathPkgFrom, toDir, TRUE
  )
  exdir <- file.path(toDir, pkg)
  if (dir.exists(exdir)) {
    unlink(exdir, TRUE)
  }
  dir.create(exdir, FALSE)
  if (remoteUntar) {
    ws <- getRemoteWs()
    pr <- guess_project()
    dir_project <- paste0(ws, "/", pr, "/", dir_relative)
    cmd <- glue::glue("cd {dir_project} && tar -xzvf {archive_package} -C {pkg} && command rm {archive_package}")
    cdRun("ssh ", remote, " '", cmd, "'")
  } else {
    untar(
      normalizePath(
        file.path(toDir, archive_package)
      ), exdir = exdir
    )
  }
}

is_sshfs_mount <- function(path = "remote") {
  if (grpl(path, "_mirror$")) {
    return(file.exists(path))
  }
  type <- system(
    glue::glue("findmnt -n -o FSTYPE --target {path}"), intern = TRUE
  )
  type == "fuse.sshfs"
}

push_script.hb <- function(..., .project = guess_project(), 
  .ws = getRemoteWs(), .path = "scripts_mirror", .exlibrary = getOption("remote_R_library", ""))
{
  project <- .project
  ws <- .ws
  path <- .path
  exlibrary <- .exlibrary
  if (!is_sshfs_mount(path)) {
    stop('!is_sshfs_mount(path).')
  }
  dir_project <- file.path(ws, project)
  maxNum <- 1L
  allFiles <- list.files(path)
  vapply(list(...), FUN.VALUE = character(1),
    function(theme) {
      pattern <- glue::glue(
        "^r\\.[0-9]{2}_{{{theme}}}\\.r$", .open = "{{{", .close = "}}}"
      )
      existFiles <- grpf(allFiles, pattern)
      num <- sprintf("%02d", maxNum)
      maxNum <<- maxNum + 1L
      pathScript <- file.path(path, glue::glue("r.{num}_{theme}.r"))
      if (length(existFiles)) {
        if (length(existFiles) > 1) {
          rlang::abort(glue::glue("Theme of {theme} found multiple files: {bind(existFiles)}"))
        } else {
          numReal <- strx(existFiles, "[0-9]+")
          if (numReal != num && sureThat("File exists: {existFiles}, rename to r.{num}_{theme}.r?")) {
            file.rename(file.path(path, existFiles), pathScript)
            return(file.path(pathScript))
          }
          return(file.path(path, existFiles))
        }
      }
      dir_output <- glue::glue("{num}_{theme}")
      script <- readLines(file.path(.expath, "job_templ", "script_setup_huibang.R"))
      script <- glue::glue(
        paste0(script, collapse = "\n"), 
        ORIGINAL_DIR = dir_project, output = dir_output,
        LIBRARY = exlibrary,
        .open = ".{{{", .close = "}}}."
      )
      writeLines(script, pathScript)
      pathScript
    })
}

push_script_runtime.hb <- function(..., .project = guess_project(), 
  .ws = getRemoteWs(), .path = "scripts_mirror", .exlibrary = getOption("remote_R_library", ""))
{
  project <- .project
  ws <- .ws
  path <- .path
  exlibrary <- .exlibrary
  if (!is_sshfs_mount(path)) {
    stop('!is_sshfs_mount(path).')
  }
  dir_project <- file.path(ws, project)
  maxNum <- 1L
  allFiles <- list.files(path)
  vapply(list(...), FUN.VALUE = character(1),
    function(theme) {
      pattern <- glue::glue(
        "^r\\.[0-9]{2}_{{{theme}}}\\.r$", .open = "{{{", .close = "}}}"
      )
      existFiles <- grpf(allFiles, pattern)
      num <- sprintf("%02d", maxNum)
      maxNum <<- maxNum + 1L
      pathScript <- file.path(path, glue::glue("r.{num}_{theme}.r"))
      if (length(existFiles)) {
        if (length(existFiles) > 1) {
          rlang::abort(glue::glue("Theme of {theme} found multiple files: {bind(existFiles)}"))
        } else {
          numReal <- strx(existFiles, "[0-9]+")
          if (numReal != num && sureThat("File exists: {existFiles}, rename to r.{num}_{theme}.r?")) {
            file.rename(file.path(path, existFiles), pathScript)
            return(file.path(pathScript))
          }
          return(file.path(path, existFiles))
        }
      }
      dir_output <- glue::glue("{num}_{theme}")
      file.copy(
        file.path(.expath, "job_templ", "script_runtime_setup_huibang.R"),
        pathScript
      )
      pathScript
    })
}



pdb_packaging.hb <- function(..., dir_save = "material") {
  project <- s(guess_project(), "[0-9]+_", "")
  files <- lapply(list(...),
    function(x) {
      if (!is(x, "job_vina") && !(x@step < 8L)) {
        stop('!is(x, "job_vina") && !(x@step < 8L).')
      }
      x$res_dock_merge$pdb_merge
    })
  files <- unlist(files)
  file_zip <- file.path(dir_save, glue::glue("{project}-{length(files)}对.zip"))
  utils::zip(file_zip, files, flags = "-j")
  gett_file(file_zip)
  return(file_zip)
}

project_packaging.hb <- function(
  file_report,
  overwrite = FALSE,
  overwrite_report = overwrite,
  path = "./remote",
  remote = "remote",
  report_share_to = "~/.var/app/com.tencent.WeChat/xwechat_files",
  wait = TRUE,
  export_pdf = TRUE)
{
  if (!is_sshfs_mount(path)) {
    stop("!is_sshfs_mount(path).", call. = FALSE)
  }
  if (!file.exists(file_report)) {
    stop("!file.exists(file_report).", call. = FALSE)
  }

  ws <- getRemoteWs()
  pr <- guess_project()

  message(glue::glue("Workspace: {ws}\nProject: {pr}"))

  prefix <- strx(pr, "(?<=_)[a-zA-Z]+")
  num_project <- strx(pr, glue::glue("(?<={prefix})[0-9]+"))

  message(glue::glue("Prefix: {prefix}"))
  message(glue::glue("Project number: {prefix}{num_project}"))

  dir_project <- paste0(ws, "/", pr)
  time <- format(Sys.Date(), "%Y%m%d")

  types <- c("scripts", "results", "report")
  names <- setNames(
    as.list(glue::glue("{prefix}_{num_project}_{types}_{time}")),
    types
  )

  all_scripts <- list.files(path, "^r\\.[0-9]+.*\\.r$")
  all_results <- gs(all_scripts, "^r\\.|\\.r$", "")

  fun_bind <- function(x) paste(shQuote(x), collapse = " ")

  cmd_sed <- glue::glue(
    "sed -i \"/^ORIGINAL_DIR\\|^.libPaths/d\" {shQuote(names$scripts)}/*"
  )

  cmd_packaging_scripts <- glue::glue(
    "cd {shQuote(dir_project)} && ",
    "rm -rf {shQuote(names$scripts)} {shQuote(paste0(names$scripts, '.zip'))} && ",
    "mkdir {shQuote(names$scripts)} && ",
    "cp -r {fun_bind(all_scripts)} -t {shQuote(names$scripts)} && ",
    "{cmd_sed} && ",
    "zip -r {shQuote(paste0(names$scripts, '.zip'))} {shQuote(names$scripts)}"
  )

  cmd_packaging_results <- glue::glue(
    "cd {shQuote(dir_project)} && ",
    "rm -rf {shQuote(names$results)} {shQuote(paste0(names$results, '.zip'))} && ",
    "mkdir {shQuote(names$results)} && ",
    "cp -r {fun_bind(all_results)} -t {shQuote(names$results)} && ",
    "zip -r {shQuote(paste0(names$results, '.zip'))} {shQuote(names$results)}"
  )

  file_scripts_zip <- file.path(path, glue::glue("{names$scripts}.zip"))
  file_results_zip <- file.path(path, glue::glue("{names$results}.zip"))

  if (!file.exists(file_scripts_zip) || !file.exists(file_results_zip) || overwrite) {
    message(glue::glue("Packaging scripts: {names$scripts} ..."))
    cdRun("ssh ", remote, " ", shQuote(cmd_packaging_scripts), wait = wait)

    message(glue::glue("Packaging results: {names$results} ..."))
    cdRun("ssh ", remote, " ", shQuote(cmd_packaging_results), wait = wait)
  }

  message("Send report file...")

  toDocx <- file.path(path, glue::glue("{names$report}.docx"))
  toPdf <- file.path(path, glue::glue("{names$report}.pdf"))

  if (!file.exists(toDocx) || overwrite_report) {
    file.copy(file_report, toDocx, overwrite = TRUE)
  }

  if (export_pdf && (!file.exists(toPdf) || overwrite_report)) {
    message("Open report by Flatpak WPS and trigger PDF export...")
    wps_pdf(toDocx)
  }

  text_reply <- glue::glue(
    "已上传分析报告:{dir_project}/{names$report}.docx和对应pdf\n\n",
    "代码压缩包:{dir_project}/{names$scripts}.zip\n\n",
    "结果文件压缩包:{dir_project}/{names$results}.zip"
  )

  # gett(text_reply)
  .gett_report(text_reply, c(toDocx, toPdf))

  if (FALSE && !is.null(report_share_to)) {
    file.copy(toDocx, report_share_to, TRUE)
    file.copy(toPdf, report_share_to, TRUE)
  }

  invisible(list(
    docx = toDocx,
    pdf = toPdf,
    scripts = file_scripts_zip,
    results = file_results_zip,
    reply = text_reply
  ))
}

.gett_report <- function(text, files) {
  gett(text)
  message("Text copied. Paste it first.")
  readline("Press Enter after pasting text...")

  gett_files(files)
  message("Files copied. Paste them now.")

  invisible(TRUE)
}

pkgVersion_remote <- function(pkgs, path = "remote",
  exlibrary = getOption("remote_R_library", ""), remote = "remote")
{
  if (!is_sshfs_mount(path)) {
    stop('!is_sshfs_mount(path).')
  }
  exlibrary <- getOption("remote_R_library", "")
  strs <- glue::glue("'{pkgs}'")
  cmd <- glue::glue("invisible(lapply(c({bind(strs)}), function(x) writeLines(as.character(packageVersion(x)))))")
  lines <- c(exlibrary, cmd)
  dir.create(file.path(path, "tmp"), FALSE)
  file_script <- file.path(path, "tmp", "getPkgInfo.R")
  writeLines(lines, file_script)
  map_file <- file.path("tmp", "getPkgInfo.R")
  ws <- getRemoteWs()
  pr <- guess_project()
  dir_project <- paste0(ws, "/", pr)
  cmd <- glue::glue("cd {dir_project} && Rscript {map_file} ")
  res <- system(paste0("ssh ", remote, " '", cmd, "'"), intern = TRUE)
  if (length(res) != length(pkgs)) {
    stop('length(res) != length(pkgs).')
  }
  res
}

release_remote_package <- function(name = "union", 
  pkg_leader = "union/union.utils", pkg_clear = "union/union.publish",
  from = "tmp/f3256d7e", to = ".", path = "remote", strip = TRUE)
{
  if (!is_sshfs_mount(path)) {
    stop('!is_sshfs_mount(path).')
  }
  codeSetup <- generate_setup_codes.huibang(
    "r.00_release.r", guess_project(),
    TRUE, path_pkg = paste0(from, "/", pkg_leader), name = "."
  )
  fun_replace <- function(x, envir = parent.frame(1L)) {
    glue::glue(x, .open = ".{{{", .close = "}}}.", envir = envir)
  }
  codeCopy <- deparse(substitute({
    dir <- paste0(to, "/", name)
    if (file.exists(dir)) {
      glue::glue("file.exists(dir), overwrite that.")
    }
    file.copy(
      ".{{{paste0(from, '/', name)}}}.", to, TRUE, TRUE
    )
  }))
  codeCopy <- fun_replace(paste0(codeCopy, collapse = "\n"))
  if (strip) {
    codePost <- deparse(substitute({
      .run_strip_semantic_layer(dir, overwrite = TRUE)
      dir <- paste0(to, "/", pkg_clear)
      .clear_autor_objects(dir, TRUE)
    }))
  } else {
    codePost <- ""
  }
  codes <- glue::as_glue(c(codeSetup, "", codeCopy, "", codePost))
  # write script
  tmpdir <- file.path(path, "tmp")
  dir.create(tmpdir, FALSE)
  fileName <- digest::digest("strip", "xxhash32", 3L)
  writeLines(codes, file.path(tmpdir, fileName))
  # run remote
  run_in_project.hb(glue::glue("tmp/{fileName}"))
}

run_remote_output.hb <- function(run = FALSE, skip = NULL,
  files = list.files(path, "^r\\.[0-9]+.*\\.r$", full.names = TRUE),
  order_by_number = TRUE,
  path = "remote", cl = NULL)
{
  if (!is_sshfs_mount(path)) {
    stop('!is_sshfs_mount(path).')
  }
  if (!file.exists(file.path(path, "union"))) {
    stop('!file.exists(file.path(path, "union")), do you forget to release that package?')
  }
  tmpdir <- file.path(path, "tmp")
  dir.create(tmpdir, FALSE)
  if (order_by_number) {
    nums <- as.integer(strx(files, "[0-9]+"))
    files <- files[order(nums)]
  }
  allCodes <- pbapply::pblapply(seq_along(files), cl = cl,
    function(n) {
      file <- files[n]
      message(glue::glue("In script ({n}): {file}"))
      if (n %in% skip) {
        return()
      }
      lines <- readLines(file)
      field_analysis <- grp(lines, "^# FIELD: analysis")
      field_output <- grp(lines, "^# FIELD: output")
      codes <- lines[-(field_analysis:field_output)]
      if (length(field_checkout <- grp(codes, "^# FIELD: checkout"))) {
        if (length(field_checkout) != 1) {
          stop('length(field_checkout) != 1, unknown error.')
        }
        codes <- codes[ -(field_checkout:length(codes)) ]
      }
      fileName <- basename(file)
      writeLines(codes, file.path(tmpdir, fileName))
      if (run) {
        run_in_project.hb(glue::glue("tmp/{fileName}"))
      }
      codes
    })
}

push_overture_as_output.hb <- function(pull = FALSE, push = FALSE,
  ovLoc = getOption("overture_codes_and_location"), override_remote = FALSE,
  dir_check = "./remote_script/push_check", path = "scripts_mirror",
  replace = "take_positions")
{
  if (is.null(ovLoc)) {
    stop('is.null(ovLoc), has not run `project_publish.complex`?')
  }
  if (!is_sshfs_mount(path)) {
    stop('!is_sshfs_mount(path).')
  }
  dir.create(dir_check, FALSE)
  dir_all <- vapply(ovLoc, function(x) x$dir, character(1))
  file_codes <- paste0("r.", dir_all, ".r")
  # setup, be careful, maybe multiple overture into the same file.
  path_codes <- file.path(path, unique(file_codes))
  if (pull) {
    file.copy(path_codes, dir_check, overwrite = TRUE)
  }
  path_codes_local <- file.path(dir_check, unique(file_codes))
  allCodes <- lapply(path_codes_local,
    function(file) {
      codes <- readLines(file)
      posMark <- grp(codes, "^# FIELD: output")
      c(codes[1:(posMark + 1)], "", "setup_counting_in_directory(output)")
    })
  names(allCodes) <- basename(path_codes_local)
  # append the codes
  lapply(seq_along(ovLoc), 
    function(n) {
      file_code <- file_codes[[n]]
      mainCodes <- allCodes[[file_code]]
      code_output <- ovLoc[[n]]$codes
      code_output <- s(code_output, replace, "output_with_counting_number", fixed = TRUE)
      codes <- c(mainCodes, "", code_output)
      allCodes[[file_code]] <<- codes
    })
  pbapply::pblapply(seq_along(path_codes_local), 
    function(n) {
      if (push) {
        path_remote <- file.path(path, basename(path_codes_local[n]))
        writeLines(c(allCodes[[n]], "", ""), path_remote)
      } else {
        writeLines(c(allCodes[[n]], "", ""), path_codes_local[n])
      }
    })
}

push_mirror_to_remote.hb <- function(dir_scripts = "scripts_mirror", path = "remote")
{
  if (!is_sshfs_mount(path)) {
    stop('!is_sshfs_mount(path).')
  }
  files <- list.files(dir_scripts, "^r\\.", full.names = TRUE)
  file.copy(files, path, TRUE)
}

push_mirror_runtime_remote.hb <- function(dir_scripts = "scripts_mirror", 
  project = guess_project(), 
  path = "remote", mode = c("release", "tune")
)
{
  if (!is_sshfs_mount(path)) {
    stop('!is_sshfs_mount(path).')
  }
  mode <- match.arg(mode)
  if (mode == "release") {
    path_pkg <- "./union/union.utils"
  } else {
    path_pkg <- "./tmp/f3256d7e/union/union.utils"
  }
  files <- list.files(dir_scripts, "^r\\.", full.names = TRUE)
  vapply(files, FUN.VALUE = logical(1L),
    function(file) {
      lines <- readLines(file)
      if (length(which <- grp(lines, "# .{{{SETUP}}}.", fixed = TRUE)) == 1L) {
        lines[ which ] <- ""
        codes <- generate_setup_codes.huibang(
          file, project, TRUE, path_pkg = path_pkg
        )
        lines <- append(lines, codes, after = which)
      }
      writeLines(lines, file.path(path, basename(file)))
      TRUE
    })
}

push_checkout_after_output.hb <- function(test = TRUE, 
  path = "scripts_mirror", dir_check = "scripts_checkout")
{
  files <- list.files(path, "^r\\.[0-9]+.*\\.r$", full.names = TRUE)
  fstart <- "# =========================================================================="
  fend <- "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
  env_class <- new.env()
  dir.create(dir_check, FALSE)
  n <- 0L
  lapply(files, 
    function(file) {
      n <<- n + 1L
      message(glue::glue("Processing File: {file}"))
      lines <- readLines(file)
      pattern_field_checkout <- "^# FIELD: checkout"
      text_note <- c("", "", .note_for_reviewer.hb, "", "")
      if (length(pos <- grp(lines, pattern_field_checkout))) {
        message(glue::glue("Field checkout found in code file, remove that."))
        lines <- lines[ seq_len(pos + 2) ]
        lines <- c(lines, text_note)
      } else {
        lines <- c(
          lines, "", "", fstart, "# FIELD: checkout", fend, text_note
        )
      }
      field_analysis <- grp(lines, "^# FIELD: analysis")
      field_output <- grp(lines, "^# FIELD: output")
      mainCodes <- lines[(field_analysis:field_output)]
      mainCodes <- parse(text = mainCodes)
      defines <- lapply(mainCodes,
        function(lang) {
          resolve <- .get_method_defination_in_package(
            lang, env_class, pkgs = basename(list_unions())
          )
          codes <- resolve$text
          if (is.null(codes)) {
            return()
          }
          codes <- codes[ !vapply(codes, is.null, logical(1)) ]
          source <- rlang::expr_text(lang)
          source <- strsplit(source, "\n")[[1]]
          note <- paste0("# ", source)
          c(note, codes, "", "")
        })
      defines <- defines[ !vapply(defines, is.null, logical(1)) ]
      if (length(defines)) {
        lines <- c(
          lines, "if (FALSE) {",
          paste0(stringr::str_pad(" ", 4), unlist(defines)),
          "}", ""
        )
      }
      if (test) {
        writeLines(lines, file.path(dir_check, basename(file)))
      } else {
        writeLines(lines, file)
      }
    })
}

pull_jobs_from_script.hb <- function(files, override = FALSE, 
  test = !override, project = guess_project(), 
  ws = getRemoteWs(), 
  path = "scripts_mirror", dir_save = "remote_script",
  pattern_object = "(?<=\\bclear\\()[a-zA-Z0-9_.]+",
  pattern_level = "(?<=\\bstep)[0-9]+(?=\\()")
{
  if (!is_sshfs_mount(path)) {
    stop('!is_sshfs_mount(path).')
  }
  if (missing(files)) {
    stop('missing(files).')
  }
  if (any(!file.exists(files))) {
    stop('any(!file.exists(files)).')
  }
  dir.create(dir_save, FALSE)
  fileNames <- basename(files)
  localFiles <- vapply(fileNames, FUN.VALUE = character(1),
    function(name) {
      local <- file.path(dir_save, name)
      if (!override && file.exists(local)) {
        stop('file.exists(local): ', local)
      }
      local
    })
  if (!test) {
    file.copy(files, dir_save, overwrite = override)
  }
  scripts <- lapply(localFiles, readLines)
  belongs <- lapply(scripts, 
    function(script) {
      unlist(stringr::str_extract_all(script, pattern_object))
    })
  belongs <- as_df.lst(belongs)
  belongs <- split(belongs$type, belongs$name)
  belongs <- sapply(names(belongs), simplify = FALSE,
    function(oname) {
      file <- unique(belongs[[oname]])
      if (length(file) > 1) {
        message(glue::glue("Detected `{oname}` from multiple file: {bind(file)}"))
        levels <- lapply(file.path(dir_save, file),
          function(file) {
            lines <- grpf(readLines(file), oname, fixed = TRUE)
            max(unlist(stringr::str_extract_all(lines, pattern_level)))
          })
        file <- file[which.max(levels)]
        message(glue::glue("{crayon::yellow(oname)} -> {file}"))
      }
      file
    })
  belongs <- lapply(belongs, 
    function(file) {
      dir <- gs(file, "^r\\.|\\.r$", "")
      list(script = file, dir = dir)
    })
  saveRDS(belongs, ".job_locate_in_script.rds")
  belongs
}

# new_script.hb <- function(theme, num = "guess", project = guess_project(), 
#   ws = getRemoteWs(), 
#   path = "remote", exlibrary = getOption("remote_R_library", ""))
# {
#   if (!is_sshfs_mount(path)) {
#     stop('!is_sshfs_mount(path).')
#   }
#   if (missing(theme)) {
#     stop('missing(theme).')
#   }
#   dir_project <- paste0(ws, "/", project)
#   pattern <- glue::glue(
#     "r\\.[0-9]{2}_{{{theme}}}\\.r", .open = "{{{", .close = "}}}"
#   )
#   existFiles <- list.files(path, pattern)
#   if (length(existFiles)) {
#     stop('length(existFiles).')
#   }
#   if (num == "guess") {
#     num <- guess_number.hb(path)
#   } else if (is.numeric(num)) {
#     num <- sprintf("%02d", as.integer(num))
#   }
#   pathScript <- file.path(path, glue::glue("r.{num}_{theme}.r"))
#   dir_output <- glue::glue("{dir_project}/{num}_{theme}")
#   script <- readLines(file.path(.expath, "job_templ", "script_setup_huibang.R"))
#   script <- glue::glue(
#     paste0(script, collapse = "\n"), 
#     ORIGINAL_DIR = dir_project, output = dir_output,
#     LIBRARY = exlibrary,
#     .open = ".{{{", .close = "}}}."
#   )
#   writeLines(script, pathScript)
#   return(pathScript)
# }

guess_number.hb <- function(path = "remote", p.pattern = "r\\.[0-9]{2}",
  n.pattern = "[0-9]{2}", type = c("files", "dirs"))
{
  type <- match.arg(type)
  if (type == "dirs") {
    alls <- list.dirs(path, recursive = FALSE)
    alls <- alls[ grpl(alls, p.pattern) ]
  } else {
    alls <- list.files(path, p.pattern)
  }
  num <- as.integer(stringr::str_extract(alls, n.pattern))
  num <- num[!is.na(num)]
  if (length(num)) {
    max <- max(num)
  } else {
    max <- 0L
  }
  sprintf("%02d", max + 1)
}

save_small.huibang <- function(name, cutoff = 50, dir = "rdata_smallObject")
{
  dir.create(dir, FALSE)
  file <- file.path(dir, glue::glue("{name}.rdata"))
  message(glue::glue("Save rdata: {file}"))
  save_small(cutoff = cutoff, file = file)
}

setup.huibang <- function() {
  options(
    tibble.print_max = 100,
    pillar.width = 100,
    pillar.max_columns = 15,
    prio_lib = "/data/nas1/huanglichuang_OD/conda/envs/extra_pkgs/lib/R/library/",
    digits = 4,
    warning.length = 5000,
    max.print = 500L,
    path_jobSave = "rds_jobSave",
    future.globals.maxSize = 5e10,
    auto_convert_plots = TRUE,
    wd_prefix = "/data/nas1/huanglichuang_OD/project/",
    db_prefix = "/data/nas1/huanglichuang_OD/project/",
    op_prefix = "/data/nas1/huanglichuang_OD/project/",
    file_batman_compounds_info = "/data/nas2/database/graphban/db/BATMAN_TCM/cids_result.csv",
    path_jobLoadFrom = list(remote = "./rds_jobSave/", local = "./rds_jobSave/lite/"),
    gwas_token = .set_gwas_token(),
    pg_local_recode = list(
      file_mbg = "/data/nas2/database/MR/MBG.allHits.p1e4.txt",
      dir_eqtl = "/data/nas2/database/MR/eQTL_vcf",
      plink_bfile = "/data/nas2/database/MR/g1000_eur/g1000_eur",
      db_drugbank = "/data/nas2/database/graphban/db/Grugbank/id2smiles.csv",
      db_scenic = "/data/nas1/huanglichuang_OD/project/SCENIC",
      # db_scenic = "/data/nas2/database/SCENIC",
      pyscenic = "conda run -n pyscenic pyscenic",
      compass = "conda run -n mebocost compass",
      cellchat_python = "/data/nas1/huanglichuang_OD/conda/envs/extra_pkgs/bin/python",
      rdkit_python = "/data/nas1/huanglichuang_OD/conda/envs/extra_pkgs/bin/python",
      conda = "/data/nas2/software/miniconda3/bin/conda",
      scsaEnv = "scsa",
      mebocostEnv = "mebocost",
      path_mebocost = "/data/nas1/huanglichuang_OD/MEBOCOST",
      scsa = "conda run -n scsa python3 /data/nas1/huanglichuang_OD/SCSA/SCSA.py",
      scsa_db = "/data/nas1/huanglichuang_OD/SCSA/whole_v2.db"
    )
  )
  options("download.file.method" = "wget", "download.file.extra" = "--no-check-certificate")
}

run_in_project_nohup.hb <- function(script, ..., limit_blas = FALSE,
  blas_threads = 1L, login_shell = FALSE, num_file = 1L)
{
  run_in_project.hb(
    script,
    ...,
    wait = FALSE,
    ex1 = "nohup",
    ex2 = glue::glue("> task_nohup_{num_file}.log 2>&1 &"),
    limit_blas = limit_blas,
    blas_threads = blas_threads,
    login_shell = login_shell
  )
}


run_in_project.hb <- function(script = "", remote = "remote",
  fun_map = NULL, wait = TRUE, ex1 = "", ex2 = "",
  limit_blas = FALSE, blas_threads = 1L, login_shell = FALSE,
  ssh_tty = FALSE, ssh_stdin_null = NULL)
{
  if (!is.null(fun_map)) {
    script <- fun_map(script)
  }

  ws <- getRemoteWs()
  pr <- guess_project()
  dir_project <- paste0(ws, "/", pr)

  cmd_env <- ""
  if (isTRUE(limit_blas)) {
    cmd_env <- glue::glue(
      "env ",
      "OPENBLAS_NUM_THREADS={blas_threads} ",
      "OPENBLAS_DEFAULT_NUM_THREADS={blas_threads} ",
      "OMP_NUM_THREADS={blas_threads} ",
      "MKL_NUM_THREADS={blas_threads} ",
      "BLIS_NUM_THREADS={blas_threads} ",
      "VECLIB_MAXIMUM_THREADS={blas_threads} "
    )
  }

  script <- shQuote(script)

  cmd_inner <- glue::glue(
    "cd {shQuote(dir_project)} && ",
    "{ex1} {cmd_env}Rscript {script} {ex2}"
  )

  if (isTRUE(login_shell)) {
    cmd_inner <- glue::glue("bash -l -c {shQuote(cmd_inner)}")
  }

  if (is.null(ssh_stdin_null)) {
    ssh_stdin_null <- !wait
  }

  vec_ssh_opt <- c(
    if (isTRUE(ssh_tty)) "-t" else character(0L),
    if (isTRUE(ssh_stdin_null)) "-n" else character(0L)
  )

  str_ssh_opt <- paste(vec_ssh_opt, collapse = " ")
  str_ssh_opt <- if (nzchar(str_ssh_opt)) paste0(str_ssh_opt, " ") else ""

  cdRun("ssh ", str_ssh_opt, remote, " ", shQuote(cmd_inner), wait = wait)
}


name.hb <- list()

name.hb$check <- function() {
  date <- format(Sys.Date(), "%m%d")
  glue::glue("{s(guess_project(), '[0-9]+_', '')}_关键节点核对_{date}")
}

name.hb$sc <- function() {
  date <- format(Sys.Date(), "%m%d")
  glue::glue("{s(guess_project(), '[0-9]+_', '')}_singleCell_{date}")
}

name.hb$report <- function() {
  date <- format(Sys.Date(), "%m%d")
  path <- glue::glue("/data/nas1/huanglichuang_OD/project/{guess_project()}")
  gett(path)
  message(path)
  glue::glue("{s(guess_project(), '[0-9]+_', '')}_Report_{date}")
}

.note_for_reviewer.hb <- c(
  "# NOTE: 下方代码是以上分析代码中解析出来的，目前只解析一层，没有递归解析",
  "# 递归的话代码会变得非常多，而且会很乱。目前应该够了，method 内部大多都是普通 function",
  "# 查看起来比较方便，可以在加载了我的 R 包后直接输入后查看本体。",
  "# ",
  "# 下方的代码，我在定义的上方写明了在上方哪个分析代码用到了这个本体，",
  "# 希望对您有所帮助"
)

.clHbFo <- list()

.clHbFo$clean_basic <- function(txt_input) {
  # remove {.mark} etc
  # txt_out <- stringr::str_replace_all(txt_input, "\\{\\.[^}]+\\}", "")
  # remove backslash
  txt_out <- stringr::str_replace_all(txt_input, "\\\\", "")
  # normalize spaces
  txt_out <- stringr::str_replace_all(txt_out, "[ \t]+", " ")
  # remove CR
  txt_out <- stringr::str_replace_all(txt_out, "\r", "")
  # replace >
  txt_out <- stringr::str_replace_all(txt_out, "^> ", "- ")
  txt_out <- stringr::str_replace_all(txt_out, "^>$", "")
  txt_out <- stringr::str_trim(txt_out)
  return(txt_out)
}

.clHbFo$normalize_headers <- function(txt_input) {
  vec_lines <- unlist(stringr::str_split(txt_input, "\n"))
  vec_out <- sapply(vec_lines, function(x_line) {
    # detect **header**
    if (stringr::str_detect(x_line, "^\\*\\*.*\\*\\*$")) {
      x_line2 <- stringr::str_replace_all(x_line, "\\*\\*", "")
      x_line2 <- stringr::str_trim(x_line2)
      return(paste0("## ", x_line2))
    }
    # detect markdown header
    if (stringr::str_detect(x_line, "^#+")) {
      str_hash <- stringr::str_extract(x_line, "^#+")
      n_lvl <- nchar(str_hash)
      n_lvl <- min(n_lvl, 3L)
      x_line2 <- stringr::str_replace(x_line, "^#+", "")
      x_line2 <- stringr::str_trim(x_line2)
      return(paste0(strrep("#", n_lvl), " ", x_line2))
    }
    return(x_line)
  }, USE.NAMES = FALSE)
  txt_out <- paste(vec_out, collapse = "\n")
  return(txt_out)
}

.clHbFo$normalize_lists <- function(txt_input) {
  vec_lines <- unlist(stringr::str_split(txt_input, "\n"))
  n_idx_lvl1 <- 0L
  vec_out <- sapply(vec_lines, function(x_line) {
    x_trim <- stringr::str_trim(x_line)

    # detect hierarchical numbering like 1.2.3） or 1.2）
    str_match <- stringr::str_match(
      x_trim,
      "^([0-9]+(\\.[0-9]+)+)([）\\).])"
    )

    if (!is.na(str_match[1])) {
      str_full <- str_match[2]
      n_level <- length(unlist(stringr::str_split(str_full, "\\.")))
      x_new <- stringr::str_replace(
        x_trim,
        "^([0-9]+(\\.[0-9]+)+)([）\\).])",
        ""
      )
      str_indent <- paste(rep(" ", (n_level - 1L) * 4L), collapse = "")
      return(paste0(str_indent, "- ", stringr::str_trim(x_new)))
    }

    # detect simple numbering: （1） or 1） or 1.
    if (stringr::str_detect(x_trim, "^[（(]?[0-9]+[）\\).]")) {
      n_idx_lvl1 <<- n_idx_lvl1 + 1L
      x_new <- stringr::str_replace(
        x_trim,
        "^[（(]?[0-9]+[）\\).]",
        ""
      )
      return(paste0(n_idx_lvl1, ". ", stringr::str_trim(x_new)))
    }

    return(x_line)

  }, USE.NAMES = FALSE)

  txt_out <- paste(vec_out, collapse = "\n")
  return(txt_out)
}

.clHbFo$fix_paragraphs <- function(txt_input) {
  # merge broken lines (not header or list)
  # txt_out <- stringr::str_replace_all(
  #   txt_input,
  #   "([^\\n])\\n([^\\n#\\-])",
  #   "\\1 \\2"
  # )
  # remove excessive blank lines
  txt_out <- stringr::str_replace_all(txt_input, "\n{3,}", "\n\n")
  return(txt_out)
}

.clHbFo$reindex_lists <- function(txt_input) {
  vec_lines <- unlist(stringr::str_split(txt_input, "\n"))
  n_idx <- 1L
  vec_out <- sapply(vec_lines, function(x_line) {
    if (stringr::str_detect(x_line, "^- ")) {
      x_new <- stringr::str_replace(x_line, "^- ", "")
      x_out <- paste0(n_idx, ". ", x_new)
      n_idx <<- n_idx + 1L
      return(x_out)
    }
    n_idx <<- 1L
    return(x_line)
  }, USE.NAMES = FALSE)
  txt_out <- paste(vec_out, collapse = "\n")
  return(txt_out)
}

.clHbFo$split_sections <- function(txt_input) {
  lst_sections <- stringr::str_split(txt_input, "\n(?=#)")[[1]]
  return(lst_sections)
}

.clHbFo$run_clean_md_list <- function(
  txt_input,
  flag_ordered = FALSE
) {
  if (!is.character(txt_input)) {
    message(glue::glue("Input is not character, coercing..."))
    txt_input <- as.character(txt_input)
  }
  txt_step <- .clHbFo$clean_basic(txt_input)
  txt_step <- .clHbFo$normalize_headers(txt_step)
  txt_step <- .clHbFo$split_sections(txt_step)
  txt_step <- unlist(lapply(txt_step, .clHbFo$normalize_lists))
  txt_step <- .clHbFo$fix_paragraphs(txt_step)
  if (flag_ordered) {
    txt_step <- .clHbFo$reindex_lists(txt_step)
    message(glue::glue("List reindexing applied."))
  }
  return(txt_step)
}

.clHbFo$format_file <- function(file, ...) {
  lines <- readLines(file)
  lines <- .clHbFo$run_clean_md_list(lines)
  writeLines(lines, file)
}


