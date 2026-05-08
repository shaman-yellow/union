# ==========================================================================
# workflow of swiss
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_swiss <- setClass("job_swiss", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("http://www.swisstargetprediction.ch/index.php"),
    cite = "[@SwisstargetpredDaina2019]",
    method = "Web tool of `SwissTargetPrediction` used for drug-targets prediction",
    tag = "target:swiss",
    analysis = "SwissTargetPrediction 药物靶点预测"
    ))

job_swiss <- function(smiles, ref = NULL)
{
  if (any(nchar(smiles) > 200)) {
    stop("any(nchar(smiles) > 200)")
  }
  x <- .job_swiss(object = smiles)
  x <- methodAdd(x, "基于 **SwissTargetPrediction** (<http://www.swisstargetprediction.ch/>) 平台对小分子化合物的潜在作用靶点进行预测分析。该平台基于化学结构相似性及已知配体–靶点信息，推测目标化合物可能结合的蛋白靶点，并提供相应的概率评分及靶点类别信息。")
  if (!is.null(ref)) {
    if (is.null(names(smiles))) {
      stop('is.null(names(smiles)), but ref is not NULL')
    }
    x$.feature_compound <- as_feature(names(smiles), ref, nature = "compound")
  }
  return(x)
}

setGeneric("asjob_swiss",
  function(x, ...) standardGeneric("asjob_swiss"))

setMethod("asjob_swiss", signature = c(x = "job_tcmsp"),
  function(x){
    data <- x@tables$step2$ingredients
    .check_columns(data, c("Mol ID", "smiles"))
    job_swiss(nl(data$`Mol ID`, data$smiles, FALSE))
  })

setMethod("asjob_swiss", signature = c(x = "job_pubchemr"),
  function(x){
    if (x@step < 1L) {
      stop("x@step < 1L")
    }
    job_swiss(unlist(x@params$smiles))
  })

setMethod("step0", signature = c(x = "job_swiss"),
  function(x){
    step_message("Prepare your data with function `job_swiss`.")
  })

setMethod("step1", signature = c(x = "job_swiss"),
  function(x, db_file = .prefix("swissTargetPrediction/targets.rds", "db"), tempdir = "download", sleep = 5, port = 4444)
  {
    step_message("Touch the online tools.")
    x$tempdir <- tempdir
    db <- new_db(db_file, ".id")
    db <- not(db, object(x))
    if (length(db@query)) {
      link <- start_drive(download.dir = x$tempdir, port = port)
      link$open()
      lapply(db@query,
        function(query) {
          link$navigate("http://www.swisstargetprediction.ch/index.php")
          ele <- link$findElement("xpath", "//form//div//input[@id='smilesBox']")
          ele$sendKeysToElement(list(query))
          Sys.sleep(3)
          Sys.sleep(sleep)
          ele <- link$findElement("xpath", "//form//div//p//input[@id='submitButton']")
          ele$clickElement()
          Sys.sleep(3)
          Sys.sleep(sleep)
          ele <- FALSE
          n <- 0L
          while ((is.logical(ele) | inherits(ele, "try-error")) & n < 20L) {
            Sys.sleep(1)
            n <- n + 1L
            ele <- try(link$findElement("xpath", "//div//button[@class='dt-button buttons-csv buttons-html5']"), TRUE)
          }
          if (!inherits(ele, "try-error")) {
            ele$clickElement()
          }
          Sys.sleep(1)
          Sys.sleep(sleep)
          Sys.sleep(20)
        })
      link$close()
      end_drive()
      ids <- paste0("query", seq_along(db@query))
      files <- collateFiles(ids, "SwissTargetPrediction.*csv", from = x$tempdir, to = x$tempdir,
        suffix = ".csv")
      data <- ftibble(files)
      names(data) <- db@query
      data <- frbind(data, idcol = TRUE, fill = TRUE)
      db <- upd(db, data)
    }
    targets <- dplyr::filter(as_tibble(db@db), .id %in% object(x))
    targets <- dplyr::rename(targets, smiles = .id, symbols = `Common name`)
    targets <- split_lapply_rbind(targets, seq_len(nrow(targets)), args = list(fill = TRUE),
      function(x) {
        if (grpl(x$symbols, " ")) {
          symbols <- strsplit(x$symbols, " ")[[1]]
          x$symbols <- NULL
          data.frame(x, symbols = symbols, check.names = FALSE)
        } else x
      })
    if (!is.null(names(object(x)))) {
      targets <- dplyr::mutate(
        targets, Name = dplyr::recode(
          smiles, !!!setNames(names(object(x)), as.character(object(x)))
        ), .before = 1
      )
    }
    x <- tablesAdd(x, targets = targets)
    colnames(targets) <- formal_name(colnames(targets))
    x$data_target <- targets
    return(x)
  })

setMethod("step2", signature = c(x = "job_swiss"),
  function(x, cut.p = .1){
    step_message("Filter data")
    if (is.null(names(object(x))) || is.null(x$data_target$Name)) {
      stop('names(object(x)) || is.null(x$data_target$Name).')
    }
    fea <- feature(x, "compound")
    data <- x$data_target
    x$.feature_all <- as_feature(
      lapply(split(data$symbols, data$Name), unique),
      "以 swissTargetPrediction 预测的活性成分的靶点"
    )
    data <- dplyr::filter(data, Probability_ > cut.p)
    x$data_filter <- data
    x$.feature_target <- as_feature(
      split(data$symbols, data$Name), "以 swissTargetPrediction 预测的活性成分的靶点"
    )
    targets <- unique(data$symbols)
    x <- snapAdd(
      x, "本研究中，以 swissTargetPrediction 预测 {snap(fea)} 的作用靶点。设定靶点概率阈值 Probability 为 {cut.p}。共得到 {length(targets)} 个唯一靶点【各成分的靶点统计：{try_snap(data, 'Name', 'symbols')}】。"
    )
    return(x)
  })

setMethod("map", signature = c(x = "job_tcmsp", ref = "job_swiss"),
  function(x, ref){
    refTargets <- ref@tables$step1$targets
    meta <- x@tables$step2$ingredients
    meta <- dplyr::select(meta, smiles, `Mol ID`, `Molecule Name`)
    refTargets <- tbmerge(refTargets, meta, by = "smiles", all.x = TRUE)
    refTargets <- dplyr::select(refTargets, `Mol ID`, `Molecule Name`,
      `Target name` = Target, symbols, probability = `Probability*`)
    refTargets <- dplyr::filter(refTargets, `Mol ID` %in% x@tables$step2$ingredients$`Mol ID`)
    x@tables$step2$compounds_targets <- refTargets
    data <- split_lapply_rbind(refTargets, ~ `Molecule Name`,
      function(x) {
        x <- dplyr::filter(x, probability > 0)
        x <- tibble::add_row(x, symbols = "Others.",
          probability = 0, `Molecule Name` = x$`Molecule Name`[[1]])
        x
      })
    p.targets <- ggplot(data) +
      geom_col(aes(x = reorder(symbols, probability), y = probability, fill = probability)) +
      labs(y = "Targets", x = "Probability") +
      coord_flip() +
      guides(fill = "none") +
      facet_wrap(~ `Molecule Name`, scales = "free_y")
    p.targets <- wrap(p.targets)
    x$p.swissTargets <- .set_lab(p.targets, sig(x), "SwissTargetPrediction-results")
    return(x)
  })
