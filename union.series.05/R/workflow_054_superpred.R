# ==========================================================================
# workflow of superpred
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_superpred <- setClass("job_superpred", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("https://prediction.charite.de/subpages/target_prediction.php"),
    cite = "[@SuperpredUpdaNickel2014]",
    method = "Web tool of `Super-PRED` used for drug-targets relationship prediction",
    tag = "target:sp",
    analysis = "Super-PRED 药物靶点预测"
    ))

setClassUnion("job_PRED", c("job_superpred", "job_swiss"))

job_superpred <- function(smiles, ref = NULL)
{
  x <- .job_superpred(object = smiles)
  x <- methodAdd(x, "基于 SuperPred (<https://prediction.charite.de/>) 平台进一步预测化合物的潜在靶点及药物功能类别。该平台结合机器学习模型、化学指纹特征及已知药物数据库信息，对化合物可能作用的蛋白靶点及 ATC 分类进行综合预测，可用于补充验证候选靶点并评估其潜在药理学属性")
  if (!is.null(ref)) {
    if (is.null(names(smiles))) {
      stop('is.null(names(smiles)), but ref is not NULL')
    }
    x$.feature_compound <- as_feature(names(smiles), ref, nature = "compound")
  }
  return(x)
}

setMethod_traceable("step0", signature = c(x = "job_superpred"),
  function(x){
    step_message("Prepare your data with function `job_superpred`.")
  })

setMethod_traceable("step1", signature = c(x = "job_superpred"),
  function(x, db_file = .prefix("superPred/targets.rds", "db"), tempdir = "download", port = 4444)
  {
    step_message("Touch the online tools.")
    x$tempdir <- tempdir
    db <- new_db(db_file, ".id")
    db <- not(db, object(x))
    if (length(db@query)) {
      link <- start_drive(download.dir = x$tempdir, port = port)
      link$open()
      ##########################
      ##########################
      # The function
      fun_get <- function(query) {
        link$navigate("https://prediction.charite.de/subpages/target_prediction.php")
        ele <- link$findElement("xpath", "//form//div//input[@id='smiles_string']")
        ele$sendKeysToElement(list(query))
        Sys.sleep(3)
        ele <- link$findElement("xpath",
          "//form//div//input[@id='smiles_string']/../div[@class='input-group-append']/button")
        # ele$sendKeysToElement(list("Search", key = "enter"))
        ele$clickElement()
        Sys.sleep(3)
        ele <- link$findElement("xpath", "//table//td/button[@type='submit']")
        # ele$sendKeysToElement(list("Start Calculation", key = "enter"))
        ele$clickElement()
        ele <- FALSE
        n <- 0L
        while ((is.logical(ele) | inherits(ele, "try-error")) & n < 20L) {
          Sys.sleep(1)
          n <- n + 1L
          ele <- try(link$findElement("xpath",
              "//div[@id='targets_wrapper']//div/button[@class='dt-button buttons-csv buttons-html5']"), TRUE)
        }
        if (!inherits(ele, "try-error")) {
          ele$clickElement()
        } else {
          writeLines(c("Target Name,Probability", ","), paste0(x$tempdir, "/", timeName("Targets"), ".csv"))
        }
        Sys.sleep(5)
      }
      ##########################
      ##########################
      ## running body
      groups <- grouping_vec2list(db@query, 5, TRUE)
      pbapply::pblapply(groups,
        function(queries) {
          lapply(queries, fun_get)
          ids <- paste0("query", seq_along(queries))
          files <- collateFiles(ids, "Targets.*csv", from = x$tempdir, to = x$tempdir,
            suffix = ".csv")
          data <- ftibble(files)
          names(data) <- queries
          data <- frbind(data, idcol = TRUE, fill = TRUE)
          db <<- upd(db, data)
          unlink(files, TRUE, TRUE)
        })
      ##########################
      ##########################
      link$close()
      end_drive()
    }
    targets <- dplyr::filter(as_tibble(db@db), .id %in% object(x))
    fun_uniport <- function(ids) {
      e(UniProt.ws::mapUniProt(
          from = "UniProtKB_AC-ID", to = "Gene_Name",
          columns = c("accession", "id"),
          query = list(ids = ids)
          ))
    }
    symbols <- expect_local_data(
      "tmp", "UniProtWs", fun_uniport, list(ids = unique(targets[[ 'UniProt ID' ]]))
    )
    targets <- map(targets, "UniProt ID", symbols, "From", "To", col = "symbols")
    if (!is.null(names(object(x)))) {
      targets <- dplyr::mutate(
        targets, Name = dplyr::recode(
          .id, !!!setNames(names(object(x)), as.character(object(x)))
          ), .before = 1
      )
    }
    targets <- dplyr::mutate(
      targets, Probability_ = as.double(strx(Probability, "[0-9.]+")) / 100
    )
    targets <- set_lab_legend(
      targets,
      glue::glue("{x@sig} targets predicted by Super-Pred"),
      glue::glue("Super-Pred 预测的化合物靶点。")
    )
    x <- tablesAdd(x, targets = targets)
    colnames(targets) <- formal_name(colnames(targets))
    x$data_target <- targets
    return(x)
  })

setMethod_traceable("step2", signature = c(x = "job_superpred"),
  function(x, cut.p = .5)
  {
    step_message("Filter data")
    if (is.null(names(object(x))) || is.null(x$data_target$Name)) {
      stop('names(object(x)) || is.null(x$data_target$Name).')
    }
    fea <- feature(x, "compound")
    data <- x$data_target
    x$.feature_all <- as_feature(
      lapply(split(data$symbols, data$Name), unique),
      "以 Super-PRED 预测的活性成分的靶点"
    )
    data <- dplyr::filter(data, Probability_ > cut.p)
    x$data_filter <- data
    x$.feature_target <- as_feature(
      split(data$symbols, data$Name), "以 Super-PRED 预测的活性成分的靶点"
    )
    targets <- unique(data$symbols)
    x <- snapAdd(
      x, "本研究中，以 Super-PRED 预测 {snap(fea)} 的作用靶点。设定靶点概率阈值 Probability 为 {cut.p * 100}% (&gt; {cut.p * 100}%。共得到 {length(targets)} 个唯一靶点【各成分的靶点统计：{try_snap(data, 'Name', 'symbols')}】。"
    )
    return(x)
  })

setGeneric("asjob_superpred",
  function(x, ...) standardGeneric("asjob_superpred"))

setMethod_traceable("asjob_superpred", signature = c(x = "job_pubchemr"),
  function(x){
    if (x@step < 1L) {
      stop("x@step < 1L")
    }
    job_superpred(unlist(x@params$smiles))
  })

setGeneric("do_herb", 
  function(x, ref, ...) standardGeneric("do_herb"))

setMethod_traceable("do_herb", signature = c(x = "job_pubchemr", ref = "job_superpred"),
  function(x, ref, disease = NULL, disease.score = 5, HLs = NULL, names = NULL, run_step3 = TRUE,
    metadata = NULL)
  {
    if (is.null(names)) {
      names <- object(x)
      if (is.null(names)) {
        stop("is.null(names)")
      }
    } else {
      object(x) <- names
    }
    plots <- list()
    if (TRUE) {
      targets <- dplyr::rename(ref@tables$step1$targets, smiles = .id)
      if (is.null(names(object(x)))) {
        stop("is.null(names(object(x)))")
      }
      dic <- nl(unlist(x$smiles, use.names = FALSE), names(x$smiles))
      dic2 <- nl(unname(object(x)), names(object(x)))
      targets <- dplyr::mutate(targets, pubchem_id = dplyr::recode(smiles, !!!dic),
        name = dplyr::recode(pubchem_id, !!!dic2)
      )
      data <- dplyr::select(targets, name, symbols, probability = Probability, `Model accuracy`)
      data <- dplyr::mutate(data, probability = as.double(strx(probability, "^[0-9]+")),
        `Model accuracy` = as.double(strx(`Model accuracy`, "^[0-9]+")))
      data <- split_lapply_rbind(data, ~ name,
        function(x) {
          x <- head(x, n = 15)
          tibble::add_row(x, symbols = "Others.", probability = 0, name = x$name[[1]],
            `Model accuracy` = 0)
        })
      p.targets <- ggplot(data) +
        geom_col(aes(x = reorder(symbols, probability), y = probability, fill = `Model accuracy`)) +
        labs(y = "Targets", x = "Probability (%)", fill = "Model accuracy (%)") +
        coord_flip() +
        facet_wrap(~ name, scales = "free_y")
      p.targets <- wrap(p.targets)
      p.targets <- .set_lab(p.targets, sig(x), "SuperPred-results")
      plots <- c(plots, namel(p.targets))
    }
    if (is.null(metadata)) {
      data <- tibble::tibble(
        Herb_pinyin_name = "PseudoHerb",
        herb_id = "PseudoHerb_id",
        Herb_cn_name = "PseudoHerb_cn",
        Ingredient.id = unname(object(x)),
        Ingredient.name = names(object(x))
      )
    } else {
      metadata <- dplyr::select(metadata, herb, cid)
      data <- tibble::tibble(cid = unname(object(x)), name = names(object(x)))
      ## the data maybe filtered in step3 (HOB), so not set 'all.x'
      data <- tbmerge(metadata, data, by = "cid", allow.cartesian = TRUE)
      data <- dplyr::mutate(data, herb_id = herb, Herb_cn_name = herb)
      data <- dplyr::rename(data, Herb_pinyin_name = herb,
        Ingredient.id = cid, Ingredient.name = name
      )
    }
    hb <- .job_herb(step = 2L)
    hb@tables$step1$herbs_compounds <- dplyr::select(
      data, herb_id, Ingredient.id, Ingredient.name
    )
    hb@tables$step2$compounds_targets <- dplyr::select(targets,
      Ingredient_id = pubchem_id, Target.name = symbols,
      Target.protein = `Target Name`
    )
    hb@params$herbs_info <- dplyr::rename(
      dplyr::distinct(data, Herb_pinyin_name, Herb_cn_name, herb_id),
      Herb_ = herb_id
    )
    hb@object$herb <- hb@params$herbs_info
    if (run_step3) {
      hb <- suppressMessages(step3(hb, disease = disease, HLs = HLs))
      hb@plots[[ 3 ]] <- c(plots, hb@plots$step3)
      hb@tables[[ 3 ]] <- namel(disease_targets_annotation = hb@tables$step3$disease_targets_annotation)
    }
    return(hb)
  })

setMethod_traceable("asjob_superpred", signature = c(x = "job_herb"),
  function(x, hob_filter = FALSE, ..., tmp = "PubChemR.rds"){
    db <- object(x)$component
    cpds <- x@tables$step1$herbs_compounds
    cpds <- map(cpds, "Ingredient.id", db, "Ingredient_id", "Ingredient_Smile", col = "smile")
    message("Use Smiles to query CIDs.")
    query <- cpds$smile %>% .[!is.na(.) & !grpl(., "Not")]
    if (!is.null(tmp)) {
      db <- new_db(tmp, "Identifier")
      db <- not(db, query)
      if (length(db@query)) {
        info <- try_get_cids.smile(db@query)
        db <- upd(db, info, db@query)
      }
      info <- dplyr::filter(db@db, Identifier %in% !!query)
    } else {
      info <- try_get_cids.smile(query)
    }
    info <- dplyr::distinct(info, Identifier, .keep_all = TRUE)
    if (TRUE) {
      if (!all(query %in% info$Identifier)) {
        print(table(query %in% info$Identifier))
        isThat <- usethis::ui_yeah("Not all smiles got CIDs, continue?")
        if (!isThat) {
          stop("...")
        }
      }
    }
    cpds <- map(cpds, "smile", info, "Identifier", "CID", col = "cid")
    query <- nl(info$CID, info$Identifier, FALSE)
    if (hob_filter) {
      ho <- job_hob(query)
      ho <- suppressMessages(step1(ho, ...))
      ifUse <- as.logical(res(ho)$prediction)
      query <- query[ ifUse ]
      .add_internal_job(.job_hob())
    }
    x <- job_superpred(query)
    if (hob_filter) {
      x$ho <- ho
    }
    x$from_herb <- cpds
    return(x)
  })

setMethod_traceable("map", signature = c(x = "job_herb", ref = "job_superpred"),
  function(x, ref){
    cpds <- ref$from_herb
    data <- dplyr::select(ref@tables$step1$targets,
      smile = .id, Target.name = symbols)
    data <- map(data, "smile", cpds, "smile", "Ingredient.id")
    data <- dplyr::rename(data, Ingredient_id = Ingredient.id)
    x@tables$step2$compounds_targets %<>% dplyr::bind_rows(data)
    x@tables$step2$compounds_targets %<>% dplyr::distinct(Ingredient_id, Target.name, .keep_all = TRUE)
    message("Set `step` to 2L.")
    x@step <- 2L
    return(x)
  })

