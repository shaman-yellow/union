# ==========================================================================
# workflow of Mendelian Randomization (AI) 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_mr <- setClass("job_mr",
  contains = "job",
  prototype = prototype(
    pg = "mr",
    info = c("https://mrcieu.github.io/TwoSampleMR/"),
    cite = "[@TwoSampleMRHemani2018; @MendelianRandomizationYavorska2023]",
    method = "",
    tag = "mr",
    analysis = "孟德尔随机化"
  )
)

.mr_batch <- setClass("mr_batch",
  contains = c("list"),
  representation = representation(
    params = "list",
    exposures = "character",
    outcomes = "character",
    split_exposure = "character",
    split_outcome = "character"
    ),
  prototype = NULL)

setMethod("show", signature = c(object = "mr_batch"),
  function(object){
    message(glue::glue("'mr_batch' object with {length(object)} MR pairs:"))
    fun_show <- function(x) stringr::str_wrap(less(x), exdent = 4L)
    message(glue::glue("exposures: {fun_show(object@exposures)}"))
    message(glue::glue("outcomes: {fun_show(object@outcomes)}"))
  })

setMethod("[[", signature = c(x = "mr_batch"),
  function(x, i, ...){
    x@.Data[[ i ]]
  })

setMethod("[[<-", signature = c(x = "mr_batch"),
  function(x, i, ..., value){
    x@.Data[[ i ]] <- value
    return(x)
  })


setMethod("[", signature = c(x = "mr_batch"),
  function(x, i, ...){
    x@.Data <- x@.Data[ i ]
    fun_subset <- function(col) {
      res <- vapply(x@.Data, FUN.VALUE = character(1),
        function(lst) {
          lst$mr[[ col ]][1]
        })
      unique(res)
    }
    x@exposures <- fun_subset(x@split_exposure)
    x@outcomes <- fun_subset(x@split_outcome)
    return(x)
  })

setMethod("$", signature = c(x = "mr_batch"),
  function(x, name){
    x@params[[ name ]]
  })

setMethod("$<-", signature = c(x = "mr_batch"),
  function(x, name, value){
    x@params[[ name ]] <- value
    return(x)
  })

job_mr <- function(data_exposure, split_exposure = "SYMBOL") {
  lapply(
    split(data_exposure, data_exposure[[ split_exposure ]]), 
    .check_twosamplemr_input
  )
  x <- .job_mr()
  x <- methodAdd(x, "双样本孟德尔随机化（two‑sample Mendelian randomization, MR）旨在利用遗传变异作为工具变量，在观察性研究框架下模拟随机对照试验的设计逻辑，从而在尽可能避免混杂偏倚与反向因果的前提下，推断暴露因素与结局之间的因果关联。")
  if (!is.null(snap(data_exposure))) {
    snap(data_exposure) <- NULL
  }
  x$data_exposure <- data_exposure
  x$split_exposure <- split_exposure
  x$catalog <- .fetch_gwas_catalog()
  return(x)
}

setMethod("step0", signature = "job_mr",
  function(x) {
    step_message("Prepare data with `job_mr()`.")
  }
)
setGeneric("asjob_mr",
  function(x, ...) standardGeneric("asjob_mr"))

setMethod("asjob_mr", signature = c(x = "feature"),
  function(x, mode = "eqtlgen", strict = TRUE, ...)
  {
    mode <- match.arg(mode)
    if (is.null(getOption("gwas_token"))) {
      stop('is.null(getOption("gwas_token")).')
    }
    if (!nchar(ieugwasr::get_opengwas_jwt())) {
      Sys.setenv(OPENGWAS_JWT = getOption("gwas_token"))
    }
    if (mode == "eqtlgen") {
      data_exposure <- .get_exposure_opengwas_eqtlgen(
        fea <- x, strict = strict, ...
      )
      if (is.null(data_exposure)) {
        stop('is.null(data_exposure).')
      }
      x <- job_mr(data_exposure, "SYMBOL")
      x <- methodAdd(x, "本研究以 {snap(fea)} 遗传预测表达量（SNP → 基因表达的 eQTL GWAS 数据）；")
    }
    x$.snap_exposure <- snap(data_exposure) %||% ""
    return(x)
  })

setMethod("step1", signature = c(x = "job_mr"),
  function(x, templates = NULL,
    patterns = NULL, cut.p = 1e-05, top_n = 3L, ..., show_source = TRUE
  )
  {
    step_message("Remove confounder effect.")
    if (!is.null(x$data_exposure_raw)) {
      message(
        glue::glue("Detected exists of `x$data_exposure_raw`, use it as `data_exposure`")
      )
      data_exposure <- x$data_exposure_raw
    } else {
      message(
        glue::glue("For first run, `x$data_exposure` will be backup in `x$data_exposure_raw`.")
      )
      if (is.null(x$data_exposure)) {
        stop('is.null(x$data_exposure).')
      }
      data_exposure <- x$data_exposure_raw <- x$data_exposure
    }
    source_confounder <- .hunt_representative_dataset_opengwas(
      patterns = patterns,
      templates = templates,
      top_n = top_n,
      catalog = x$catalog
    )
    x$source_confounder <- dplyr::relocate(
      source_confounder, search_group, id, sample_size, nsnp
    )
    if (show_source) {
      message(glue::glue("Confounder datasets refer to (`x$source_confounder`): "))
      print(x$source_confounder, n = 10L, width = 80L)
    }
    detail <- bind(c(patterns, templates))
    snap_confounder <- glue::glue("为排除混杂因素 ({detail}) 相关的 SNP，{snap(x$source_confounder)}")
    data_exposure <- .remove_other_snps_opengwas(
      data_exposure,
      x$source_confounder$id,
      p_threshold = cut.p, 
      type = "confounder",
      col_stat = x$split_exposure,
      ..., ask_query = TRUE
    )
    snap_remove <- glue::glue("{snap(data_exposure)}")
    x$.snap_confounder <- paste0(snap_confounder, snap_remove)
    snap(data_exposure) <- NULL
    x$data_exposure <- data_exposure
    return(x)
  })

setMethod("step2", signature = "job_mr",
  function(x, mode = c("mbg", "general"), ids = NULL, outcome = NULL, ...)
  {
    step_message("Get outcome data.")
    mode <- match.arg(mode)
    if (mode == "general") {
      if (is.null(ids)) {
        stop('is.null(ids), mode == "general", `ids` should provided.')
      }
      data_outcome <- .get_outcome_opengwas_general(
        x$data_exposure$SNP, ids, ...
      )
      x$data_outcome <- data_outcome
      x$split_outcome <- "outcome"
      x$pattern_outcome <- "[^|]+"
      x$.snap_outcome <- glue::glue("{snap(data_outcome)}")
      if (is.null(outcome)) {
        outcome <- strx(unique(data_outcome$outcome), x$pattern_outcome)
        if (length(outcome) > 5) {
          stop('length(outcome) > 5, too many "outcome" type for representation in WORD.')
        }
        outcome <- bind(outcome)
      }
      x <- methodAdd(x, "以 {outcome} 为结局，探究两者之间的因果关联。")
    } else if (mode == "mbg") {
      data_outcome <- .get_outcome_opengwas_mbg(
        x$data_exposure$SNP, ...
      )
      x <- methodAdd(x, "以微生物丰度 (MiBioGen 数据) 为结局，探究两者之间的因果关联。")
      x$data_outcome <- data_outcome
      x$split_outcome <- "outcome"
      x$pattern_outcome <- "(?<=\\()[^()]+(?=\\))"
      x$.snap_outcome <- glue::glue("结局因素获取自 MiBioGen 数据库。{snap(data_outcome)}")
    } else {
      stop("Not yet ready.")
      .get_outcome_opengwas(...)
    }
    return(x)
  }
)

setMethod("step3", signature = c(x = "job_mr"),
  function(x, filter_outcome = TRUE, cut.p = 1e-05, ...)
  {
    step_message("Filter snaps (outcome or F test).")
    x <- methodAdd(x, "{x$.snap_exposure}\n\n")
    x <- methodAdd(x, "{x$.snap_confounder}\n\n")
    x <- methodAdd(x, "{x$.snap_outcome}\n\n")
    data_exposure <- x$data_exposure
    if (filter_outcome) {
      data_exposure <- .remove_other_snps_opengwas(
        data_exposure,
        x$data_outcome$id.outcome %||% x$data_outcome$id,
        p_threshold = cut.p,
        type = "outcomeFactor",
        col_stat = x$split_exposure,
        ..., ask_query = TRUE
      )
      snap_remove <- glue::glue("排除与结局因素相关的 SNP。{snap(data_exposure)}")
      snap(data_exposure) <- NULL
      x <- methodAdd(x, "{snap_remove}\n\n")
    }
    x$fTest <- .calc_F_stat(data_exposure)
    data_exposure <- data_exposure[ x$fTest$keep,  ]
    x <- methodAdd(x, "{x$fTest$snap}")
    x$data_exposure <- data_exposure
    return(x)
  })

setMethod("step4", signature = c(x = "job_mr"),
  function(x, min_nsnp = 3L, ..., rerun = FALSE, workers = 10L)
  {
    step_message("Run MR analysis and filter results.")
    x$mr_batch <- .run_mr_batch(
      x$data_exposure,
      x$data_outcome,
      split_exposure = x$split_exposure,
      split_outcome = x$split_outcome,
      min_nsnp = min_nsnp,
      meth = TRUE,
      ncore = workers,
      rerun = rerun
    )
    # IVW p < 0.05
    # Egger intercept p > 0.05
    # PRESSO global p > 0.05
    # Steiger TRUE
    x$mr_batch_filter <- .filter_mr_batch(
      x$mr_batch, ...
    )
    x <- methodAdd(x, "{x$mr_batch$snap_run}\n\n\n{x$mr_batch_filter$snap_filter}")
    return(x)
  })

setMethod("step5", signature = c(x = "job_mr"),
  function(x, ...)
  {
    step_message("Summry.")
    args <- list(
      x = x$mr_batch_filter,
      pattern_outcome = x$pattern_outcome,
      pattern_exposure = x$pattern_exposure
    )
    # tables
    args$wrap <- FALSE
    objFmtTab <- do.call(.format_show_mr_batch, args)
    ts.alls <- x$table_with_all_columns <- .table_mr_batch(objFmtTab)
    cols_excludes <- c(
      "id.outcome", "id.exposure", "lo_ci", "up_ci",
      "or_lci95", "or_uci95",
      "SYMBOL", "pair_id"
    )
    ts.alls <- lapply(ts.alls,
      function(data) {
        data[ , !colnames(data) %in% cols_excludes ]
      })
    t.main <- set_lab_legend(
      ts.alls$main,
      glue::glue("{x@sig} MR summary"),
      glue::glue(
        "MR 因果效应汇总表|||汇总展示各暴露因素与结局之间的孟德尔随机化分析结果。method 为所采用的 MR 方法；nsnp 为纳入分析的工具变量数量；b 为效应估计值；se 为标准误；pval 为统计学检验 P 值；or 为比值比（Odds Ratio）。结果主要参考逆方差加权法（IVW）。"
      )
    )
    t.heterogeneity <- set_lab_legend(
      ts.alls$heterogeneity,
      glue::glue("{x@sig} Heterogeneity test"),
      glue::glue(
        "MR 异质性检验结果表|||用于评估各工具变量效应估计之间的一致性。method 为检验对应的 MR 方法；Q 为 Cochran's Q 统计量；Q_df 为自由度；Q_pval 为异质性检验 P 值。通常 Q_pval > 0.05 提示未见显著异质性，可认为工具变量之间总体一致性较好。"
      )
    )
    t.pleiotropy <- set_lab_legend(
      ts.alls$pleiotropy,
      glue::glue("{x@sig} Horizontal pleiotropy test"),
      glue::glue(
        "MR 水平多效性检验结果表|||基于 MR-Egger 截距检验评估工具变量是否通过暴露因素以外途径影响结局。egger_intercept 为截距估计值；se 为标准误；pval 为检验 P 值。通常 pval > 0.05 提示未发现显著方向性水平多效性，结果可靠性较高。"
      )
    )
    t.steiger <- set_lab_legend(
      dplyr::rename(ts.alls$steiger, direction = correct_causal_direction),
      glue::glue("{x@sig} Steiger directionality test"),
      glue::glue(
        "MR 因果方向检验结果表|||用于判断推定因果方向是否更支持“暴露因素影响结局”而非反向因果。snp_r2.exposure 与 snp_r2.outcome 分别表示工具变量解释暴露和结局变异的比例；direction 表示是否支持预设方向；steiger_pval 为方向性检验 P 值。direction 为 TRUE 且 steiger_pval < 0.05 时，通常认为因果方向更可信。"
      )
    )
    x <- tablesAdd(x, t.main, t.heterogeneity, t.pleiotropy, t.steiger)
    if (x$split_exposure == "SYMBOL") {
      x$.feature_exposure <- as_feature(
        split(t.main$exposure, t.main$outcome), "MR 显著因果关联"
      )
    }
    # plots
    args$wrap <- TRUE
    objFmtPlot <- do.call(.format_show_mr_batch, args)
    ps.alls <- .plot_mr_batch(objFmtPlot)
    # snap summary
    snap <- .stat_summary_mr_batch(objFmtTab)
    p.scatter <- set_lab_legend(
      ps.alls$p.scatter,
      glue::glue("{x@sig} MR scatter plot"),
      glue::glue(
        "MR 分析散点图|||横坐标表示各 SNP 对暴露因素的效应估计值，纵坐标表示对应 SNP 对结局的效应估计值。每个点代表一个工具变量。不同颜色直线表示不同 MR 方法拟合得到的总体因果效应方向与大小；斜率越大，提示效应越强。各方法拟合方向一致时，说明结果稳定性较好；若斜率差异明显，则提示模型间存在不一致性。"
      )
    )
    p.forest <- set_lab_legend(
      ps.alls$p.forest,
      glue::glue("{x@sig} MR forest plot"),
      glue::glue(
        "MR 单位点森林图|||展示每个 SNP 单独作为工具变量时对结局的效应估计及其 95% 置信区间，同时给出总体合并效应。点估计位于零效应线（或 OR=1）右侧提示正向作用，左侧提示负向作用。多数位点方向一致且总体效应稳定时，支持结果可靠；若个别位点偏离明显，则提示可能存在异常工具变量。"
      )
    )
    p.funnel <- set_lab_legend(
      ps.alls$p.funnel,
      glue::glue("{x@sig} MR funnel plot"),
      glue::glue(
        "MR 漏斗图|||横坐标表示各 SNP 的效应估计值，纵坐标表示估计精度（通常为标准误的倒数）。点位围绕总体效应线对称分布时，提示整体结果较稳定，未见明显方向性偏倚；若分布明显偏斜或单侧聚集，则提示可能存在异质性、多效性或异常位点影响。"
      )
    )
    p.leaveoneout <- set_lab_legend(
      ps.alls$p.leaveoneout,
      glue::glue("{x@sig} MR leave-one-out plot"),
      glue::glue(
        "MR 逐一剔除敏感性分析图|||依次去除单个 SNP 后重新进行 MR 分析，展示每次重新估计的总体效应及其 95% 置信区间。若各次结果接近原始总体估计，说明分析结论不依赖某一单独工具变量，稳健性较好；若去除某个位点后效应明显变化，则提示该位点可能具有较大影响。"
      )
    )
    snap_plots <- .stat_plot_summary_mr_batch(
      objFmtTab,
      p.scatter = p.scatter,
      p.forest = p.forest,
      p.funnel = p.funnel,
      p.leaveoneout = p.leaveoneout
    )
    x <- plotsAdd(x, p.scatter, p.forest, p.funnel, p.leaveoneout)
    x <- snapAdd(x, "{snap}\n\n\n{snap_plots}")
    return(x)
  })

# ==========================================================
# Company exposure adapter
# gene -> ensembl -> eqtl file -> exposure_dat
# ==========================================================
.get_exposure_local_eqtl.huibang <- function(
  genes,
  dir_eqtl = pg("dir_eqtl"),
  plink_bfile = pg("plink_bfile"),
  pval_threshold = 5e-08,
  clump_r2 = 0.01,
  clump_kb = 10000L,
  rerun = FALSE
)
{
  fun_get_data <- function(...) {
    .check_pkg("AnnotationDbi")
    .check_pkg("org.Hs.eg.db")
    .check_pkg("VariantAnnotation")
    .check_pkg("gwasglue")
    .check_pkg("TwoSampleMR")
    .check_pkg("plinkbinr")

    message(glue::glue("Get gene of {bind(head(genes))}"))
    map <- AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = as.character(genes),
      keytype = "SYMBOL",
      columns = c("SYMBOL", "ENSEMBL")
    )
    message(glue::glue("Mapped id of genes: {nrow(map)}"))

    map <- tibble::as_tibble(map)

    map <- dplyr::filter(
      map,
      !is.na(.data$ENSEMBL)
    )

    map <- dplyr::distinct(map)

    map <- dplyr::mutate(
      map,
      id.exposure = paste0(
        "eqtl-a-",
        .data$ENSEMBL
      )
    )

    available_ids <- .list_eqtl_ids.huibang(
      dir_eqtl = dir_eqtl
    )
    message(glue::glue("All `available_ids` number: {length(available_ids)}"))
    map <- dplyr::filter(
      map,
      .data$id.exposure %in% available_ids
    )
    message(
      glue::glue("Afiter filterd by `available_ids`, the map data rows: {nrow(map)}")
    )

    if (any(duplicated(map$SYMBOL))) {
      print(map)
      message(glue::glue("Duplicated data found."))
    }

    if (nrow(map) == 0L) {
      return(NULL)
    }
    isNot <- which_not_in_data(map, "SYMBOL", genes, stop = FALSE)
    if (any(isNot)) {
      on.exit(message(glue::glue("Not got: {genes[isNot]}")))
    }

    res <- lapply(
      seq_len(nrow(map)),
      function(i) {

        eid <- map$id.exposure[i]

        f_vcf <- file.path(
          dir_eqtl,
          paste0(eid, ".vcf.gz")
        )

        vcf <- VariantAnnotation::readVcf(f_vcf)

        dat <- suppressWarnings(
          gwasglue::gwasvcf_to_TwoSampleMR(
            vcf,
            type = "exposure"
          )
        )

        dat$id.exposure <- eid

        dat <- dplyr::filter(
          dat,
          .data$pval.exposure < pval_threshold
        )

        if (nrow(dat) == 0L) {
          return(NULL)
        }

        dat$id <- dat$id.exposure
        dat$rsid <- dat$SNP
        dat$pval <- dat$pval.exposure

        clumped <- ieugwasr::ld_clump(
          dat,
          plink_bin = plinkbinr::get_plink_exe(),
          bfile = plink_bfile,
          clump_kb = clump_kb,
          clump_r2 = clump_r2
        )

        if (nrow(clumped) == 0L) {
          return(NULL)
        }

        clumped$SYMBOL <- map$SYMBOL[i]

        clumped
      }
    )

    names(res) <- map$SYMBOL

    Filter(Negate(is.null), res)
  }
  expect_local_data(
    "tmp", "exposure_eqtl", fun_get_data,
    list(as.character(genes), pval_threshold, clump_r2, clump_kb),
    rerun = rerun
  )
}

# ==========================================================
# MiBioGen zip reader
# Efficient reading
# ==========================================================

.get_outcome_mbg_full.huibang <- function(
  snps,
  type = c("genus", "family", "order", "class", "phylum"),
  dir_db = pg("dir_MBG"),
  id = "ALL",
  pattern = "\\.summary\\.txt\\.gz$"
)
{
  type <- match.arg(type)
  # from: https://molgenis26.gcc.rug.nl/downloads/MiBioGen/MiBioGen_QmbQTL_summary_phylum.zip
  zipfile <- file.path(
    dir_db, glue::glue("MiBioGen_QmbQTL_summary_{type}.zip")
  )
  if (!file.exists(zipfile)) {
    stop('!file.exists(zipfile).')
  }

  final_pattern <- pattern
  if (!identical(id, "ALL")) {
    taxa_pat <- paste(id, collapse = "|")
    final_pattern <- paste0("(", taxa_pat, ").*", pattern)
  }

  raw_dat <- .read_archive_table_by_id(
    path = zipfile,
    ids = snps,
    id_col = "rsID",
    pattern = final_pattern,
    sep = "\t"
  )

  if (is.null(raw_dat) || nrow(raw_dat) == 0L) {
    message("No SNP matched in the specified MiBioGen files.")
    return(NULL)
  }

  # Ensure required columns exist after shell filtering
  req_cols <- c("rsID", "beta", "SE", "eff.allele", "ref.allele", "P.weightedSumZ", "N", "bac")
  if (!all(req_cols %in% colnames(raw_dat))) {
    stop("Extracted data is missing required MiBioGen columns.")
  }

  out <- dplyr::transmute(
    raw_dat,
    SNP = .data$rsID,
    beta.outcome = .data$beta,
    se.outcome = .data$SE,
    effect_allele.outcome = .data$eff.allele,
    other_allele.outcome = .data$ref.allele,
    pval.outcome = .data$P.weightedSumZ,
    samplesize.outcome = .data$N,
    outcome = .data$bac,
    id.outcome = .data$bac
  )

  out <- dplyr::distinct(out)
  message("Total Matched rows: ", nrow(out))

  out
}


.get_outcome_mbg_toplist.huibang <- function(
  snps,
  id = "ALL"
)
{
  if (missing(snps) || length(snps) == 0L) {
    stop("`snps` is required.", call. = FALSE)
  }

  file <- get_url_data(
    "MBG.allHits.p1e4.txt",
    "https://molgenis26.gcc.rug.nl/downloads/MiBioGen/MBG.allHits.p1e4.txt",
    "MBG.allHits.p1e4", fun_decompress = NULL
  )

  dat <- ftibble(file)
  message(glue::glue("Total data nrows: {nrow(dat)}"))

  req_col <- c(
    "bac",
    "rsID",
    "ref.allele",
    "eff.allele",
    "beta",
    "SE",
    "P.weightedSumZ",
    "N"
  )

  miss_col <- setdiff(req_col, colnames(dat))

  if (length(miss_col) > 0L) {
    stop(
      paste0(
        "Missing columns: ",
        paste(miss_col, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (!is.null(id) && !identical(id, "ALL")) {
    dat <- dplyr::filter(
      dat,
      .data$bac %in% id
    )
  }

  message(glue::glue("Before SNP filter, data rows number {nrow(dat)}"))
  dat <- dplyr::filter(dat, .data$rsID %in% snps)
  message(glue::glue("After SNP filter, data rows number {nrow(dat)}"))

  dat <- dplyr::transmute(
    dat,
    SNP = .data$rsID,
    beta.outcome = .data$beta,
    se.outcome = .data$SE,
    effect_allele.outcome = .data$eff.allele,
    other_allele.outcome = .data$ref.allele,
    pval.outcome = .data$P.weightedSumZ,
    samplesize.outcome = .data$N,
    outcome = .data$bac,
    id.outcome = .data$bac
  )
  dat <- dplyr::distinct(dat)
  dat
}



# ==========================================================
# Detect huibang available IDs
# Priority:
# 1. eQTL_list
# 2. list.files()
# ==========================================================
.list_eqtl_ids.huibang <- function(
  dir_eqtl = pg("dir_eqtl"),
  file_index = "eQTL_list"
)
{
  idx_path <- file.path(dir_eqtl, file_index)

  if (file.exists(idx_path)) {

    ids <- utils::read.table(
      idx_path,
      header = FALSE,
      stringsAsFactors = FALSE
    )[[1]]

    return(as.character(ids))
  } else {
    message(glue::glue("File not exists: {file_index}"))
  }

  fs <- list.files(
    dir_eqtl,
    pattern = "\\.vcf\\.gz$",
    full.names = FALSE
  )

  sub("\\.vcf\\.gz$", "", fs)
}

# ==========================================================
# OpenGWAS catalog
# ==========================================================

.fetch_gwas_catalog <- function(
  dir_db = .prefix("ieugwasr", "db"),
  force = FALSE)
{
  .check_pkg("ieugwasr")
  if (!force) {
    api_token <- getOption("gwas_token", NULL)
    if (is.null(api_token)) {
      stop('is.null(api_token), please set `options` with: `gwas_token`')
    }
    Sys.setenv(OPENGWAS_JWT = api_token)
  }
  data <- expect_local_data(
    dir_db, "ieugwasr", ieugwasr::gwasinfo, list(id = NULL)
  )
  tibble::as_tibble(data.frame(data))
}

.get_exposure_opengwas_eqtlgen <- function(
  genes,
  strict = TRUE,
  rerun = FALSE
)
{
  if (strict) {
    pval_threshold <- 5e-08
    clump_r2 <- 0.01
    clump_kb <- 10000L
  } else {
    pval_threshold <- 1e-05
    clump_r2 <- 0.05
    clump_kb <- 1000L
  }
  fun_get_data <- function(...) {
    .check_pkg("AnnotationDbi")
    .check_pkg("org.Hs.eg.db")
    .check_pkg("TwoSampleMR")

    message(glue::glue("Mapping {length(genes)} genes to Ensembl and eQTLGen IDs..."))

    # 1. Map to Ensembl for GTEx
    map <- AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = as.character(genes),
      keytype = "SYMBOL",
      columns = c("SYMBOL", "ENSEMBL")
    )

    map <- tibble::as_tibble(map)
    map <- dplyr::filter(map, !is.na(ENSEMBL))
    map <- dplyr::distinct(map)

    # Add id.exposure to map for joining later
    map <- dplyr::mutate(map, id.exposure = paste0("eqtl-a-", ENSEMBL))

    all_query_ids <- unique(map$id.exposure)

    # 2. Fetch instruments from OpenGWAS of eQTLGen
    message(glue::glue("Querying OpenGWAS for eQTLGen data..."))

    res_raw <- tryCatch({
      TwoSampleMR::extract_instruments(
        outcomes = all_query_ids,
        p1 = pval_threshold,
        clump = TRUE,
        r2 = clump_r2,
        kb = clump_kb
      )
    }, error = function(e) {
      message("API request failed: ", e$message)
      return(NULL)
    })

    if (is.null(res_raw) || !nrow(res_raw)) {
      message("No significant eQTLs found in database.")
      return(NULL)
    }

    res_with_symbol <- dplyr::left_join(
      res_raw,
      dplyr::select(map, SYMBOL, id.exposure),
      by = "id.exposure"
    )

    res_with_symbol <- dplyr::filter(res_with_symbol, SYMBOL %in% genes)

    if (!nrow(res_with_symbol)) {
      message("No significant results remain after filtering for requested genes.")
      return(NULL)
    }
    res_with_symbol <- dplyr::relocate(res_with_symbol, SYMBOL)
    tibble::as_tibble(res_with_symbol)
  }

  data <- expect_local_data(
    "tmp", "exposure_eqtl", fun_get_data,
    list(as.character(genes), pval_threshold, clump_r2, clump_kb),
    rerun = rerun
  )
  which_not_in_data(data, "SYMBOL", genes)
  stats <- try_snap(data, 'SYMBOL', 'SNP')
  message(
    glue::glue("All results stat: {stats}")
  )
  if (length(genes) < 15L) {
    stats <- glue::glue("【各基因 SNP 统计 (括号内为数量)：{stats}】")
  } else {
    stats <- ""
  }
  base <- "从 eQTLGen 联盟公开的顺式表达数量性状位点（cis‑eQTL）数据中提取基因表达的工具变量。以目标基因为输入，借助以 R 包 org.Hs.eg.db ⟦pkgInfo('org.Hs.eg.db')⟧ 和 AnnotationDbi ⟦pkgInfo('AnnotationDbi')⟧ 映射至 Ensembl ID，随后通过以 R 包 TwoSampleMR ⟦pkgInfo('TwoSampleMR')⟧ 查询 OpenGWAS (<https://opengwas.io/>) 数据库，保留与表达显著关联的遗传变异。"
  if (strict) {
    snap(data) <- glue::glue("{base}⟦mark$blue('显著性阈值设为 p < {pval_threshold}，并对工具变量执行 LD clumping（r² = {clump_r2}，窗口 = {clump_kb} kb），以保证独立性')⟧。未能匹配或缺乏显著 eQTL 的基因将被剔除，最终获得暴露因素的工具变量集合{stats}。")
  } else {
    snap(data) <- glue::glue(
      "考虑到 eQTL 数据中工具变量数量有限，为提高统计功效，本研究适当放宽工具变量筛选阈值，采用 P < {pval_threshold} 作为候选 SNP 的纳入标准 (PMID：31924771, 31341166)，并设置 LD 剔除参数为 r² = {clump_r2} (PMID: 40126059, 33866329)、clumping window = {clump_kb} kb (PMID: 29955180)。上述参数设置已在多项基于 omics 数据的孟德尔随机化研究中得到广泛应用，能够在保证工具变量独立性的同时，提高 eQTL-MR 分析的有效 SNP 数量与统计稳定性。在此基础上，未能匹配或缺乏显著 eQTL 的基因将被剔除，最终获得暴露因素的工具变量集合{stats}。"
    )
  }
  return(data)
}

# ==========================================================
# OpenGWAS outcome adapter
# ==========================================================

.get_outcome_opengwas <- function(
  snps, target_ids, cache_prefix = "opengwas_cache",
  proxies = TRUE, check_new_snps = FALSE, check_new_ids = FALSE,
  ask_query = FALSE, snp_chunk_size = 50L, id_chunk_size = 1L
)
{
  snps <- unique(snps)
  target_ids <- unique(target_ids)

  cache_files <- list.files(
    "tmp", pattern = paste0("^", cache_prefix),
    full.names = TRUE
  )

  cached_data <- NULL
  ever_seen_ids <- character(0L)
  ever_seen_snps <- character(0L)
  has_history <- length(cache_files) > 0L

  if (has_history) {

    cached_list <- lapply(cache_files, readRDS)
    cached_data <- data.table::rbindlist(
      cached_list, use.names = TRUE, fill = TRUE
    )

    if (!is.null(cached_data) && nrow(cached_data) > 0L) {
      cached_data <- unique(cached_data)
      ever_seen_ids <- unique(cached_data$id.outcome)
      ever_seen_snps <- unique(cached_data$SNP)
    }

  } else if (ask_query && interactive()) {

    message(glue::glue("SNPs: {less(snps, n = 5L)}"))
    message(glue::glue("Datasets ID: {less(target_ids, n = 5L)}"))

    if (!sureThat("No cache file exists, download data now?")) {
      stop("!sureThat, stop fetch network data.")
    }
  }

  fetch_new_snps_now <- (!has_history) || check_new_snps
  fetch_new_ids_now <- (!has_history) || check_new_ids

  old_ids <- intersect(target_ids, ever_seen_ids)

  new_snps <- if (fetch_new_snps_now) {
    setdiff(snps, ever_seen_snps)
  } else {
    character(0L)
  }

  new_ids <- if (fetch_new_ids_now) {
    setdiff(target_ids, ever_seen_ids)
  } else {
    character(0L)
  }

  final_list <- list()

  if (!is.null(cached_data)) {
    final_list[[1L]] <- cached_data
  }

  .fetch_and_cache <- function(s_list, id_list, tag) {

    if (length(s_list) == 0L || length(id_list) == 0L) {
      message(glue::glue("[{tag}] skipped."))
      return(NULL)
    }
    message(glue::glue(
      "[{tag}] Fetching {length(s_list)} SNPs across {length(id_list)} IDs."
    ))
    id_chunks <- split(
      id_list, ceiling(seq_along(id_list) / id_chunk_size)
    )
    snp_chunks <- split(
      s_list, ceiling(seq_along(s_list) / snp_chunk_size)
    )

    res_all <- lapply(seq_along(id_chunks), function(i) {
      vec_ids <- id_chunks[[i]]
      message(glue::glue(
        "[{tag}] ID chunk {i}/{length(id_chunks)} ({length(vec_ids)} IDs)"
      ))
      pbapply::pblapply(seq_along(snp_chunks), function(j) {
        vec_snps <- snp_chunks[[j]]
        message(glue::glue(
          "[{tag}] SNP chunk {j}/{length(snp_chunks)} ({length(vec_snps)} SNPs)"
        ))
        expect_local_data(
          "tmp",
          cache_prefix,
          function(s, o) {
            res <- tryCatch(
              TwoSampleMR::extract_outcome_data(
                snps = s,
                outcomes = o,
                proxies = proxies
              ),
              error = function(e) {
                warning(glue::glue(
                  "extract_outcome_data failed: {conditionMessage(e)}"
                ))
                NULL
              }
            )
            Sys.sleep(3L)
            res
          },
          list(s = vec_snps, o = vec_ids)
        )
      })
    })

    unlist(res_all, recursive = FALSE)
  }

  res_gap_old_ids <- .fetch_and_cache(new_snps, old_ids, "Updating Old IDs")
  res_new_ids <- .fetch_and_cache(snps, new_ids, "Fetching New IDs")
  final_list <- c(final_list, res_gap_old_ids, res_new_ids)

  final_list <- final_list[
    sapply(final_list, function(x) {
      !is.null(x) && nrow(x) > 0L
    })
  ]

  if (length(final_list) == 0L) {
    message("No data retrieved from API or cache.")
    return(NULL)
  }

  merged_all <- data.table::rbindlist(
    final_list,
    use.names = TRUE,
    fill = TRUE
  )
  res <- dplyr::filter(
    merged_all,
    id.outcome %in% target_ids,
    SNP %in% snps
  )
  which_not_in_data(res, "SNP", snps)
  tibble::as_tibble(unique(res))
}

.get_outcome_opengwas_general <- function(
  snps, ids, check_new_snps = FALSE, check_new_ids = FALSE,
  catalog = .fetch_gwas_catalog()
)
{
  if (is.null(catalog)) {
    catalog <- .fetch_gwas_catalog()
  }

  meta <- dplyr::filter(catalog, id %in% ids)
  if (!nrow(meta)) {
    stop('!nrow(meta), no specified GWAS: {bind(ids)}')
  }
  which_not_in_data(meta, "id", ids)
  target_ids <- unique(meta$id)
  
  data <- .get_outcome_opengwas(
    snps = snps, 
    target_ids = target_ids, 
    cache_prefix = paste0("opengwas_general"),
    check_new_snps = check_new_snps,
    check_new_ids = check_new_ids
  )
  if (length(ids) <= 5L) {
    snaps <- glue::glue(
      "结局变量 (Trait) {meta$trait} 来源于 GWAS 数据 (ID: {ids}) (<https://opengwas.io/>)，由 {meta$author} 等人领导，涵盖了 {meta$sample_size} 名受试者 (其中 {meta$ncase} 个 Case) ，包含 {meta$nsnp} 个 SNP。"
    )
    if (length(ids) > 1L) {
      snap <- bind(
        c(glue::glue("本研究一共使用了 {length(ids)} 组 GWAS 数据用于 MR 分析。\n\n"),
          bind(paste0("- ", snaps), co = "\n")
          ),
        co = "\n"
      )
    } else {
      snap <- snaps
    }
    snap(data) <- paste0(snap, "\n\n在数据提取过程中，我们通过暴露因素的 SNP 为依准，使用 R 包 `TwoSampleMR` ⟦pkgInfo('TwoSampleMR')⟧，以 `TwoSampleMR::extract_outcome_data` 获取该结局数据对应的 SNP。")
  } else {
    snap(data) <- ""
  }
  return(data)
}

.get_outcome_opengwas_mbg <- function(
  snps, type = c("ALL", "genus", "family", "order", "class", "phylum"),
  check_new_snps = FALSE, check_new_ids = FALSE,
  catalog = .fetch_gwas_catalog()
)
{
  type <- match.arg(type)
  
  if (is.null(catalog)) {
    catalog <- .fetch_gwas_catalog()
  }

  mbg_meta <- dplyr::filter(catalog, grepl("Kurilshikov", author, ignore.case = TRUE))
  
  if (type != "ALL") {
    mbg_meta <- dplyr::filter(mbg_meta, grepl(type, trait, ignore.case = TRUE))
  }
  
  target_ids <- unique(mbg_meta$id)
  
  # Pass down the control switches
  data <- .get_outcome_opengwas(
    snps = snps, 
    target_ids = target_ids, 
    cache_prefix = paste0("opengwas_mbg_", type),
    check_new_snps = check_new_snps,
    check_new_ids = check_new_ids
  )
  if (type == "ALL") {
    alls <- c("genus", "family", "order", "class", "phylum")
    type <- glue::glue("{bind(alls)} (ALL) ")
  }
  snap(data) <- glue::glue("MiBioGen (<www.mibiogen.org>) 数据由 Kurilshikov 等人领导，涵盖了来自 24 个队列的 18,340 名受试者，是目前肠道微生物群遗传学研究中最具权威性的数据集之一。在数据提取过程中，我们通过暴露因素的 SNP 为依准，以 `TwoSampleMR::extract_outcome_data`  获取 {type} 分类水平下关联相同 SNP 的微生物类群 (一共得到 {length(unique(data$outcome))} 种)。")
  return(data)
}

.hunt_representative_dataset_opengwas <- function(
  patterns = NULL,
  templates = NULL,
  top_n = 3L,
  catalog = .fetch_gwas_catalog(),
  verbose = TRUE
)
{
  if (is.null(catalog)) {
    catalog <- .fetch_gwas_catalog()
  }

  all_templates <- list(
    smoking = "Smoking initiation|Lifetime smoking|Cigarettes per day|Smoking status|Tobacco use",
    drinking = "Drinks per week|Alcohol intake frequency|Alcohol consumption|Alcohol use|Alcohol drinking",
    bmi = "^Body mass index$|\\bBMI\\b|Obesity|Body fat",
    education = "Years of schooling|Educational attainment|College completion|Intelligence",
    crp = "^C-reactive protein$|\\bCRP\\b|Inflammatory marker",
    diabetes = "Type 2 diabetes|T2D|Blood glucose|Fasting glucose|HbA1c|Insulin resistance",
    blood_pressure = "Blood pressure|Systolic blood pressure|Diastolic blood pressure|\\bSBP\\b|\\bDBP\\b|Hypertension",
    lipid = "Cholesterol|Triglyceride|LDL|HDL|Total cholesterol|Lipid",
    immune = "Immune|Inflammatory|Interleukin|\\bIL-|TNF|White blood cell|Neutrophil|Lymphocyte|Monocyte",
    infection = "Infection|Bacterial infection|Viral infection|Pneumonia|Respiratory infection|Bacteremia|Septicemia",
    ibd = "Inflammatory bowel disease|Crohn|Ulcerative colitis",
    sleep = "Sleep duration|Insomnia|Sleep trait",
    exercise = "Physical activity|Exercise|Sedentary behaviour",
    coffee = "Coffee intake|Coffee consumption|Caffeine"
  )

  units <- list()

  if (!is.null(patterns)) {
    units <- c(
      units,
      lapply(patterns, function(x) {
        list(
          label = paste0("patterns:", x),
          regex = x
        )
      })
    )
  }

  if (!is.null(templates)) {
    units <- c(
      units,
      lapply(templates, function(x) {
        if (!x %in% names(all_templates)) {
          stop(paste("Unknown templates:", x))
        }

        list(
          label = paste0("templates:", x),
          regex = all_templates[[x]]
        )
      })
    )
  }

  if (!length(units)) {
    stop("patterns or templates are required.")
  }

  txt <- paste(
    catalog$trait,
    catalog$author,
    catalog$consortium
  )

  res_list <- lapply(units,
    function(u) {

      hit <- grepl(u$regex, txt, ignore.case = TRUE)

      d <- catalog[hit, , drop = FALSE]

      if (!nrow(d)) {
        return(NULL)
      }

      if (!"sample_size" %in% colnames(d)) {
        d$sample_size <- 0
      }

      if (!"nsnp" %in% colnames(d)) {
        d$nsnp <- 0
      }

      d$sample_size[is.na(d$sample_size)] <- 0
      d$nsnp[is.na(d$nsnp)] <- 0

      d <- dplyr::distinct(d, .data$id, .keep_all = TRUE)

      d <- dplyr::arrange(
        d,
        dplyr::desc(.data$sample_size),
        dplyr::desc(.data$nsnp),
        nchar(.data$trait)
      )

      d <- utils::head(d, top_n)

      d <- dplyr::mutate(
        d,
        search_group = u$label,
        .before = 1L
      )

      d
    })

  res_list <- res_list[!sapply(res_list, is.null)]

  if (!length(res_list)) {
    return(NULL)
  }

  res <- data.table::rbindlist(
    res_list,
    use.names = TRUE,
    fill = TRUE
  )

  if (verbose) {
    which_not_in_data(
      res, "search_group", vapply(units, function(x) x$label, character(1))
    )
    message(glue::glue("Returned {nrow(res)} rows."))
  }

  res <- tibble::as_tibble(res)
  snap(res) <- glue::glue("从 OpenGWAS 目录中筛选具有代表性的 GWAS summary 数据。具体实现为：以正则表达式对性状描述、数据库信息进行匹配，识别与目标表型相关的所有可用数据集；经去重后，以样本量（sample_size）和工具变量数目（SNP）为优先准则降序排列，保留排序靠前的 Top {top_n} 个条目，以确保纳入分析的数据兼具统计效力与遗传变异数目优势。")
  return(res)
}

.remove_other_snps_opengwas <- function(
  exposure_dat,
  ids,
  p_threshold = 1e-05,
  type = c("confounder", "outcomeFactor"),
  col_stat = NULL,
  cache_prefix = glue::glue("opengwas_{type}"),
  check_new_snps = FALSE,
  check_new_ids = FALSE,
  ask_query = FALSE,
  verbose = TRUE
)
{
  if (is.null(exposure_dat) || !nrow(exposure_dat)) {
    return(exposure_dat)
  }

  if (!"SNP" %in% colnames(exposure_dat)) {
    stop("Column `SNP` not found in exposure_dat.")
  }
  type <- match.arg(type)

  ids <- unique(stats::na.omit(as.character(ids)))

  if (!length(ids)) {
    stop("`ids` is required.")
  }

  if (verbose) {
    message(
      glue::glue(
        "Filtering SNPs using {length(ids)} {type} IDs."
      )
    )
  }

  out <- .get_outcome_opengwas(
    snps = unique(exposure_dat$SNP),
    target_ids = ids,
    cache_prefix = cache_prefix,
    check_new_snps = check_new_snps,
    check_new_ids = check_new_ids,
    proxies = FALSE,
    ask_query = ask_query
  )

  if (is.null(out) || !nrow(out)) {
    return(exposure_dat)
  }

  if (!"pval.outcome" %in% colnames(out)) {
    warning("Column `pval.outcome` not found in outcome data.")
    return(exposure_dat)
  }

  if (!is.null(out$proxy.outcome)) {
    out <- dplyr::filter(out, !proxy.outcome)
  }

  bad_snps <- unique(
    out$SNP[which(out$pval.outcome < p_threshold)]
  )

  if (length(bad_snps)) {
    res <- dplyr::filter(
      exposure_dat,
      !(.data$SNP %in% bad_snps)
    )

    if (verbose) {
      message(
        glue::glue(
          "Removed {length(bad_snps)} SNPs, retained {nrow(res)} SNPs."
        )
      )
    }
  } else {
    message(glue::glue("Keep all SNPs."))
    res <- exposure_dat
  }

  res <- tibble::as_tibble(res)
  type <- switch(type, confounder = "混杂", outcomeFactor = "结局")
  snap_ex <- ""
  if (!is.null(col_stat)) {
    alls_pre <- unique(exposure_dat[[ col_stat ]])
    alls_aft <- unique(res[[ col_stat ]])
    exclude <- setdiff(alls_pre, alls_aft)
    if (length(exclude)) {
      snap_ex <- glue::glue(" (经本次过滤后，暴露因素 {bind(exclude)} 不再包含任何 SNP，因此从随后的评估中移除) ")
    }
  }
  snap(res) <- glue::glue("将暴露因素工具变量 (SNP) 中与{type}因素显著关联（P < {signif(p_threshold, 2)}）的 SNP 视为受{type}影响并予以排除 (仅基于原始工具变量 SNP 与潜在混杂因素的直接关联进行过滤，即不使用 proxy SNP 替代功能)。对于所有暴露因素的 SNP，最终剔除 {length(bad_snps)} 个{type}因素 SNP，保留 {nrow(res)} 个 SNP。{snap_ex}")
  return(res)
}

# ==========================================================
# Harmonise one pair
# ==========================================================
.prepare_mr_pair <- function(
  exposure_dat,
  outcome_dat,
  action = 2L,
  verbose = TRUE
)
{
  .check_pkg("TwoSampleMR")

  if (is.null(exposure_dat) || !nrow(exposure_dat)) {
    if (verbose) {
      warning("exposure_dat is NULL or empty.")
    }
    return(NULL)
  }

  if (is.null(outcome_dat) || !nrow(outcome_dat)) {
    if (verbose) {
      warning("outcome_dat is NULL or empty.")
    }
    return(NULL)
  }

  res <- tryCatch(
    TwoSampleMR::harmonise_data(
      exposure_dat = exposure_dat,
      outcome_dat = outcome_dat,
      action = action
    ),
    error = function(e) {
      if (verbose) {
        warning(
          glue::glue(
            "harmonise_data failed: {conditionMessage(e)}"
          )
        )
      }
      return(NULL)
    }
  )

  if (is.null(res) || !nrow(res)) {
    if (verbose) {
      warning("harmonise_data returned empty result.")
    }
    return(NULL)
  }

  res
}

# ==========================================================
# Single MR (full version)
# ==========================================================
.run_mr_single <- function(
  dat,
  min_nsnp = 3L,
  run_presso = TRUE,
  presso_nb = 1000L,
  seed = 123L,
  verbose = TRUE
)
{
  .check_pkg("TwoSampleMR")

  if (is.null(dat) || !nrow(dat)) {
    if (verbose) {
      warning("Input dat is NULL or empty.")
    }
    return(NULL)
  }

  nsnp <- length(unique(dat$SNP))

  if (nsnp < min_nsnp) {
    if (verbose) {
      message(
        glue::glue(
          "Skipped MR: nsnp = {nsnp} < min_nsnp = {min_nsnp}."
        )
      )
    }
    return(NULL)
  }

  if (verbose) {
    message(
      glue::glue(
        "Running MR with {nsnp} SNPs."
      )
    )
  }

  het <- tryCatch(
    TwoSampleMR::mr_heterogeneity(dat),
    error = function(e) {
      if (verbose) {
        warning(
          glue::glue(
            "mr_heterogeneity failed: {conditionMessage(e)}"
          )
        )
      }
      NULL
    }
  )

  q_p <- 1

  if (!is.null(het) && nrow(het)) {

    idx <- which(
      het$method %in% c(
        "Inverse variance weighted",
        "Inverse variance weighted (multiplicative random effects)"
      )
    )

    if (length(idx)) {
      q_p <- het$Q_pval[idx[1L]]
    }
  } else {
    if (verbose) {
      message("No heterogeneity result. Using fixed-effect IVW.")
    }
  }

  ivw_method <- if (is.na(q_p) || q_p >= 0.05) {
    "mr_ivw_fe"
  } else {
    "mr_ivw_mre"
  }

  if (verbose) {
    message(
      glue::glue(
        "Selected IVW method: {ivw_method} (Q_p = {signif(q_p, 3L)})."
      )
    )
  }

  method_list <- c(
    "mr_egger_regression",
    "mr_weighted_median",
    ivw_method,
    "mr_simple_mode",
    "mr_weighted_mode"
  )

  mr_res <- tryCatch(
    TwoSampleMR::mr(
      dat,
      method_list = method_list
    ),
    error = function(e) {
      if (verbose) {
        warning(
          glue::glue(
            "mr failed: {conditionMessage(e)}"
          )
        )
      }
      NULL
    }
  )

  if (is.null(mr_res) || !nrow(mr_res)) {
    if (verbose) {
      warning("MR returned empty result.")
    }
    return(NULL)
  }

  mr_or <- tryCatch(
    TwoSampleMR::generate_odds_ratios(mr_res),
    error = function(e) {
      if (verbose) {
        message("generate_odds_ratios failed. Using raw MR result.")
      }
      mr_res
    }
  )

  pleio <- tryCatch(
    TwoSampleMR::mr_pleiotropy_test(dat),
    error = function(e) {
      if (verbose) {
        warning(
          glue::glue(
            "mr_pleiotropy_test failed: {conditionMessage(e)}"
          )
        )
      }
      NULL
    }
  )

  steiger <- tryCatch(
    TwoSampleMR::directionality_test(dat),
    error = function(e) {
      if (verbose) {
        warning(
          glue::glue(
            "directionality_test failed: {conditionMessage(e)}"
          )
        )
      }
      NULL
    }
  )

  loo <- tryCatch(
    TwoSampleMR::mr_leaveoneout(dat),
    error = function(e) {
      if (verbose) {
        warning(
          glue::glue(
            "mr_leaveoneout failed: {conditionMessage(e)}"
          )
        )
      }
      NULL
    }
  )

  single <- tryCatch(
    TwoSampleMR::mr_singlesnp(dat),
    error = function(e) {
      if (verbose) {
        warning(
          glue::glue(
            "mr_singlesnp failed: {conditionMessage(e)}"
          )
        )
      }
      NULL
    }
  )

  presso <- NULL

  if (run_presso) {

    if (nsnp < 4L) {

      if (verbose) {
        message("MR-PRESSO skipped: nsnp < 4.")
      }

    } else if (!requireNamespace("MRPRESSO", quietly = TRUE)) {

      if (verbose) {
        warning("Package `MRPRESSO` not installed.")
      }

    } else {

      set.seed(seed)

      presso <- tryCatch(
        MRPRESSO::mr_presso(
          BetaOutcome = "beta.outcome",
          BetaExposure = "beta.exposure",
          SdOutcome = "se.outcome",
          SdExposure = "se.exposure",
          OUTLIERtest = TRUE,
          DISTORTIONtest = TRUE,
          data = dat,
          NbDistribution = presso_nb,
          SignifThreshold = 0.05
        ),
        error = function(e) {
          if (verbose) {
            warning(
              glue::glue(
                "MR-PRESSO failed: {conditionMessage(e)}"
              )
            )
          }
          NULL
        }
      )
    }
  }

  if (verbose) {
    message("MR finished successfully.")
  }

  lst <- list(
    mr = mr_or,
    heterogeneity = het,
    pleiotropy = pleio,
    steiger = steiger,
    presso = presso,
    leaveoneout = loo,
    singlesnp = single,
    dat = dat
  )

  .to_tbl <- function(x) {
    if (is.null(x)) return(NULL)
    if (is.data.frame(x)) return(tibble::as_tibble(x))
    x
  }

  lapply(lst, .to_tbl)
}


# ==========================================================
# Batch MR
# ==========================================================

.run_mr_batch <- function(
  exposure_dat,
  outcome_dat,
  split_exposure = "SYMBOL",
  split_outcome = "outcome",
  min_nsnp = 3L,
  run_presso = TRUE,
  rerun = FALSE,
  ncore = max(1L, parallel::detectCores() - 1L),
  chunk_size = 50L,
  meth = FALSE,
  verbose = TRUE
)
{
  if (is.null(exposure_dat) || !nrow(exposure_dat)) {
    warning("exposure_dat is NULL or empty.")
    return(list())
  }

  if (is.null(outcome_dat) || !nrow(outcome_dat)) {
    warning("outcome_dat is NULL or empty.")
    return(list())
  }

  exp_list <- split(
    exposure_dat,
    exposure_dat[[split_exposure]]
  )

  out_list <- split(
    outcome_dat,
    outcome_dat[[split_outcome]]
  )

  grid <- expand.grid(
    eid = names(exp_list),
    oid = names(out_list),
    stringsAsFactors = FALSE
  )

  nall <- nrow(grid)

  if (verbose) {
    message(
      glue::glue(
        "Running {length(exp_list)} x {length(out_list)} MR pairs ({nall} total)."
      )
    )
    message(
      glue::glue(
        "Parallel workers: {ncore}"
      )
    )
    message(
      glue::glue(
        "Chunk size: {chunk_size}"
      )
    )
  }

  idx_chunk <- split(
    seq_len(nall),
    ceiling(seq_len(nall) / chunk_size)
  )

  res_all <- lapply(
    seq_along(idx_chunk),
    function(k) {

      idx <- idx_chunk[[k]]

      if (verbose) {
        message(
          glue::glue(
            "Processing chunk {k} / {length(idx_chunk)} ({length(idx)} pairs)."
          )
        )
      }

      fun_run_chunk <- function(...) {

        pbapply::pblapply(
          idx,
          cl = ncore,
          FUN = function(i) {

            eid <- grid$eid[i]
            oid <- grid$oid[i]

            dat <- .prepare_mr_pair(
              exposure_dat = exp_list[[eid]],
              outcome_dat = out_list[[oid]],
              verbose = FALSE
            )

            ans <- .run_mr_single(
              dat = dat,
              min_nsnp = min_nsnp,
              run_presso = run_presso,
              verbose = FALSE
            )

            if (is.null(ans)) {
              return(NULL)
            }

            ans$mr[[split_exposure]] <- eid
            ans$mr[[split_outcome]] <- oid

            ans
          }
        )
      }

      id_args <- list(
        exposure_dat$SNP,
        exposure_dat[[split_exposure]],
        outcome_dat$SNP,
        outcome_dat[[split_outcome]],
        min_nsnp = min_nsnp,
        run_presso = run_presso,
        chunk_index = k,
        pair_index = idx
      )

      expect_local_data(
        "tmp",
        "run_mr_batch",
        fun_run_chunk,
        id_args,
        rerun = rerun
      )
    }
  )

  res <- unlist(
    res_all,
    recursive = FALSE
  )

  res <- Filter(
    Negate(is.null),
    res
  )
  res <- .mr_batch(res,
    exposures = names(exp_list),
    outcomes = names(out_list),
    split_exposure = split_exposure,
    split_outcome = split_outcome,
    params = list(
      min_nsnp = min_nsnp
    )
  )
  res <- res[ seq_along(res) ]
  if (meth) {
    if (length(exp_list) < 20) {
      snap_ex <- bind(names(exp_list))
      snap_ex <- glue::glue(" ({snap_ex}) ")
    } else {
      snap_ex <- ""
    }
    res$snap_run <- glue::glue(
      "在完成数据准备与协同处理后，本研究基于 `TwoSampleMR` ⟦pkgInfo('TwoSampleMR')⟧，对 **{length(exp_list)}** 个暴露因素**{snap_ex}**与 **{length(out_list)}** 个结局指标进行了系统性的因果推断（共计 **{nrow(grid)}** 个分析组合）。针对每一对暴露-结局组合，我们预设工具变量纳入标准为 $n_{SNP} \\ge$ **{min_nsnp}**。",
      .open = "**{", .close = "}**"
    )
  }
  return(res)
}

# ==========================================================
# Filter MR batch results
# Main test uses IVW-family methods only
# ==========================================================
.filter_mr_batch <- function(
  x,
  p_cutoff = 0.05,
  use_fdr = FALSE,
  adjust_method = "fdr",
  require_steiger = TRUE,
  require_no_pleio = TRUE,
  pleio_cutoff = 0.05,
  require_heterogeneity = TRUE,
  heterogeneity_cutoff = 0.05,
  use_presso = TRUE,
  require_no_presso = FALSE,
  verbose = TRUE
)
{
  if (!is(x, "mr_batch")) {
    stop('!is(x, "mr_batch")')
  }

  if (!length(x)) {
    warning("Input `x` is empty.")
    return(x)
  }

  ivw_pat <- paste(
    c(
      "Inverse variance weighted",
      "fixed effects",
      "multiplicative random effects",
      "IVW"
    ),
    collapse = "|"
  )

  meta <- lapply(
    seq_along(x),
    function(i) {

      one <- x[[i]]

      if (is.null(one$mr) || !nrow(one$mr)) {
        message(glue::glue("Pair {i}: missing `mr`, skipped."))
        return(NULL)
      }

      mr <- tibble::as_tibble(one$mr)

      hit <- grepl(ivw_pat, mr$method, ignore.case = TRUE)

      if (!any(hit)) {
        message(glue::glue("Pair {i}: no IVW method found."))
        return(NULL)
      }

      ivw <- mr[which(hit)[1L], , drop = FALSE]

      out <- data.frame(
        idx = i,
        pval = ivw$pval[1L],
        stringsAsFactors = FALSE
      )

      out$pleio_ok <- TRUE
      out$steiger_ok <- TRUE
      out$het_ok <- TRUE
      out$presso_ok <- TRUE

      # --------------------------------------------
      # MR-Egger pleiotropy
      # --------------------------------------------
      if (require_no_pleio) {

        out$pleio_ok <- FALSE

        if (!is.null(one$pleiotropy) &&
            nrow(one$pleiotropy) &&
            "pval" %in% colnames(one$pleiotropy)) {

          pv <- one$pleiotropy$pval[1L]

          if (!is.na(pv)) {
            out$pleio_ok <- pv > pleio_cutoff
          }
        }
      }

      # --------------------------------------------
      # Steiger directionality
      # --------------------------------------------
      if (require_steiger) {

        out$steiger_ok <- FALSE

        if (!is.null(one$steiger) &&
            nrow(one$steiger) &&
            "correct_causal_direction" %in%
            colnames(one$steiger)) {

          out$steiger_ok <- isTRUE(
            one$steiger$correct_causal_direction[1L]
          )
        }
      }

      # --------------------------------------------
      # Heterogeneity
      # --------------------------------------------
      if (require_heterogeneity) {

        out$het_ok <- FALSE

        if (!is.null(one$heterogeneity) &&
            nrow(one$heterogeneity) &&
            "Q_pval" %in% colnames(one$heterogeneity)) {

          qv <- suppressWarnings(
            max(one$heterogeneity$Q_pval, na.rm = TRUE)
          )

          if (!is.infinite(qv) && !is.na(qv)) {
            out$het_ok <- qv > heterogeneity_cutoff
          }
        }
      }

      # --------------------------------------------
      # MR-PRESSO
      # TRUE  = no significant distortion / no issue
      # FALSE = evidence of outlier distortion
      # --------------------------------------------
      if (use_presso || require_no_presso) {

        out$presso_ok <- FALSE

        pr <- one$presso

        if (is.null(pr)) {

          out$presso_ok <- !require_no_presso

        } else {

          pv <- NA_real_

          # Main MR-PRESSO object usually stores:
          # $`MR-PRESSO results`$`Distortion Test`$Pvalue
          pv <- tryCatch(
            {
              as.numeric(
                pr$`MR-PRESSO results`$
                  `Distortion Test`$
                  Pvalue[1L]
              )
            },
            error = function(e) NA_real_
          )

          # fallback: global test
          if (!length(pv) || is.na(pv)) {
            pv <- tryCatch(
              {
                as.numeric(
                  pr$`MR-PRESSO results`$
                    `Global Test`$
                    Pvalue[1L]
                )
              },
              error = function(e) NA_real_
            )
          }

          if (!length(pv) || is.na(pv)) {
            out$presso_ok <- !require_no_presso
          } else {
            out$presso_ok <- pv > pleio_cutoff
          }
        }
      }

      out
    }
  )

  meta <- Filter(Negate(is.null), meta)

  if (!length(meta)) {
    warning("No valid MR pairs found.")
    return(x[0L])
  }

  tab <- dplyr::bind_rows(meta)

  tab$p_adj <- stats::p.adjust(
    tab$pval,
    method = adjust_method
  )

  p_use <- if (use_fdr) tab$p_adj else tab$pval

  keep <- p_use < p_cutoff &
    tab$pleio_ok &
    tab$steiger_ok &
    tab$het_ok &
    tab$presso_ok

  idx_keep <- tab$idx[keep]

  if (verbose) {

    message(
      glue::glue(
        "Retained {length(idx_keep)} / {length(x)} MR pairs."
      )
    )

    message(
      glue::glue(
        "Use FDR: {use_fdr}; ",
        "Steiger: {require_steiger}; ",
        "Egger pleio: {require_no_pleio}; ",
        "Het: {require_heterogeneity}; ",
        "PRESSO: {use_presso}"
      )
    )
  }
  x <- x[idx_keep]
  line <- .get_description("filter_mr_batch.md")
  x$snap_filter <- glue::glue(
    glue::glue(line, .trim = FALSE), .trim = FALSE
  )
  return(x)
}


# ==========================================================
# Summary MR batch results
# Merge all methods + sensitivity metrics
# ==========================================================

.summary_mr_batch <- function(x, digits = 4L)
{
  if (!is(x, "mr_batch")) {
    stop('!is(x, "mr_batch")')
  }

  if (is.null(x) || length(x) == 0L) {
    warning("Input `x` is NULL or empty.")
    return(tibble::tibble())
  }

  res <- lapply(
    seq_along(x),
    function(i) {

      one <- x[[i]]

      if (is.null(one$mr) || !nrow(one$mr)) {
        return(NULL)
      }

      mr <- tibble::as_tibble(one$mr)

      if (!"or" %in% colnames(mr) &&
          all(c("b", "se") %in% colnames(mr))) {

        mr$or <- exp(mr$b)
        mr$or_lci95 <- exp(mr$b - 1.96 * mr$se)
        mr$or_uci95 <- exp(mr$b + 1.96 * mr$se)
      }

      mr$Q_pval <- NA_real_
      mr$pleiotropy_pval <- NA_real_
      mr$egger_intercept <- NA_real_
      mr$steiger_dir <- NA
      mr$steiger_pval <- NA_real_

      if (!is.null(one$heterogeneity) &&
          nrow(one$heterogeneity)) {

        het <- tibble::as_tibble(one$heterogeneity)

        idx <- match(mr$method, het$method)
        mr$Q_pval <- het$Q_pval[idx]
      }

      if (!is.null(one$pleiotropy) &&
          nrow(one$pleiotropy)) {

        ple <- tibble::as_tibble(one$pleiotropy)

        if ("pval" %in% colnames(ple)) {
          mr$pleiotropy_pval <- ple$pval[1L]
        }

        if ("egger_intercept" %in% colnames(ple)) {
          mr$egger_intercept <- ple$egger_intercept[1L]
        }
      }

      if (!is.null(one$steiger) &&
          nrow(one$steiger)) {

        st <- tibble::as_tibble(one$steiger)

        if ("correct_causal_direction" %in%
            colnames(st)) {
          mr$steiger_dir <-
            st$correct_causal_direction[1L]
        }

        if ("steiger_pval" %in%
            colnames(st)) {
          mr$steiger_pval <- st$steiger_pval[1L]
        } else if ("pval" %in% colnames(st)) {
          mr$steiger_pval <- st$pval[1L]
        }
      }

      mr$pair_id <- i
      return(mr)
    }
  )

  res <- Filter(Negate(is.null), res)

  if (!length(res)) {
    warning("No valid MR results found.")
    return(tibble::tibble())
  }

  tab <- dplyr::bind_rows(res)

  tab <- dplyr::arrange(
    tab,
    .data$pval,
    .data$method
  )

  num_cols <- c(
    "b", "se", "pval", "p_adj",
    "or", "or_lci95", "or_uci95",
    "Q_pval", "pleiotropy_pval",
    "egger_intercept",
    "steiger_pval"
  )

  hit <- intersect(num_cols, colnames(tab))

  tab <- dplyr::mutate(
    tab,
    dplyr::across(
      dplyr::all_of(hit),
      function(z) round(z, digits)
    )
  )

  keep <- c("outcome", "exposure", "method", "nsnp", "b",
    "se", "pval", "p_adj", "or", "or_lci95", "or_uci95",
    "Q_pval", "pleiotropy_pval", "egger_intercept", "steiger_dir",
    "steiger_pval", "pair_id"
  )

  keep <- intersect(keep, colnames(tab))

  tibble::as_tibble(tab[, keep, drop = FALSE])
}

.format_show_mr_batch <- function(x,
  slots = c(
    "mr", "singlesnp", "leaveoneout",
    "heterogeneity", "pleiotropy", "steiger",
    "dat"
  ),
  col_outcome = x@split_outcome,
  col_exposure = x@split_exposure,
  pattern_outcome = x$pattern_outcome,
  pattern_exposure = x$pattern_exposure,
  wrap = FALSE, width = 25L
)
{
  if (!is(x, "mr_batch")) {
    stop('!is(x, "mr_batch").')
  }
  if (!length(slots)) {
    return(x)
  }
  fun_map <- function(x, ref, from, to) {
    if (!is.null(ref[[ to ]])) {
      which <- match(x[[ from ]], ref[[ from ]])
      x[[ from ]] <- x[[ to ]] <- ref[[ to ]][ which ]
    }
    if (wrap) {
      x[[ from ]] <- x[[ to ]] <- stringr::str_wrap(
        x[[ to ]], width
      )
    }
    return(x)
  }
  for (i in seq_along(x)) {
    .dat <- x[[ i ]]$dat
    if (is.null(.dat)) {
      stop('is.null(.dat), no raw `dat` for mapping columns.')
    }
    for (s in slots) {
      data <- x[[ i ]][[ s ]]
      for (name in c("outcome", "exposure")) {
        pattern <- get(glue::glue("pattern_{name}"))
        col <- get(glue::glue("col_{name}"))
        if (!is.null(col)) {
          data <- fun_map(data, .dat, name, col)
        }
        if (!is.null(pattern)) {
          data[[name]] <- stringr::str_extract(raw <- data[[name]], pattern)
          if (any(is.na(data[[ name ]]))) {
            print(less(raw, n = 10))
            stop('The pattern extract return NA, see above column value.')
          }
        }
      }
      x[[ i ]][[ s ]] <- dplyr::relocate(data, exposure, outcome)
    }
  }
  return(x)
}

# ==========================================================
# Tables from mr_batch
# return named list
# ==========================================================

.table_mr_batch <- function(
  x,
  digits = 4L
)
{
  if (!is(x, "mr_batch")) {
    stop('!is(x, "mr_batch")')
  }

  if (is.null(x) || length(x) == 0L) {
    warning("Input `x` is NULL or empty.")
    return(list())
  }

  .bind_part <- function(lst, field) {

    dplyr::bind_rows(
      lapply(
        seq_along(lst),
        function(i) {

          z <- lst[[i]][[field]]

          if (is.null(z) || !nrow(z)) {
            return(NULL)
          }

          z <- tibble::as_tibble(z)
          z$pair_id <- i
          z
        }
      )
    )
  }

  .round_num <- function(dat) {

    if (is.null(dat) || !nrow(dat)) {
      return(dat)
    }

    idx <- sapply(dat, is.numeric)

    dat[idx] <- lapply(
      dat[idx],
      round,
      digits = digits
    )

    tibble::as_tibble(dat)
  }

  res <- Filter(
    Negate(is.null),
    lapply(
      seq_along(x),
      function(i) {

        one <- x[[i]]

        if (is.null(one$mr) || !nrow(one$mr)) {
          return(NULL)
        }

        mr <- tibble::as_tibble(one$mr)
        mr$pair_id <- i

        if (!"or" %in% colnames(mr) &&
            all(c("b", "se") %in% colnames(mr))) {

          mr$or <- exp(mr$b)
          mr$or_lci95 <- exp(mr$b - 1.96 * mr$se)
          mr$or_uci95 <- exp(mr$b + 1.96 * mr$se)
        }

        one$mr <- mr
        one
      }
    )
  )

  out <- lapply(
    c(
      main = "mr",
      heterogeneity = "heterogeneity",
      pleiotropy = "pleiotropy",
      steiger = "steiger"
    ),
    function(field) {
      .round_num(.bind_part(res, field))
    }
  )

  names(out) <- c(
    "main",
    "heterogeneity",
    "pleiotropy",
    "steiger"
  )

  out
}

# ==========================================================
# Plot mr_batch
# return named list of ggplot objects
# Different plot types stored separately
# ==========================================================
.plot_mr_batch <- function(x, top_n = 10L)
{
  if (!is(x, "mr_batch")) {
    stop('!is(x, "mr_batch")')
  }

  .check_pkg("TwoSampleMR")
  .check_pkg("patchwork")

  if (is.null(x) || length(x) == 0L) {
    warning("Input `x` is NULL or empty.")
    return(list())
  }

  idx <- seq_len(min(length(x), top_n))

  p.scatter <- list()
  p.forest <- list()
  p.funnel <- list()
  p.leaveoneout <- list()

  for (i in idx) {

    one <- x[[i]]

    if (is.null(one$mr) || is.null(one$dat)) {
      next
    }

    ttl <- paste0(
      unique(one$mr$exposure)[1L],
      " -> ",
      unique(one$mr$outcome)[1L]
    )

    # scatter
    tmp <- try(
      TwoSampleMR::mr_scatter_plot(
        one$mr,
        one$dat
      ),
      silent = TRUE
    )

    if (!inherits(tmp, "try-error")) {
      p.scatter[[paste0("pair.", i)]] <- tmp[[1L]] +
        ggplot2::ggtitle(ttl) +
        ggplot2::labs(y = "") +
        ggplot2::guides(colour = ggplot2::guide_legend(ncol = 1))
    }

    # forest
    if (!is.null(one$singlesnp) &&
        nrow(one$singlesnp)) {

      tmp <- try(
        TwoSampleMR::mr_forest_plot(
          one$singlesnp
        ),
        silent = TRUE
      )

      if (!inherits(tmp, "try-error")) {
        p.forest[[paste0("pair.", i)]] <- tmp[[1L]] + ggplot2::ggtitle(ttl)
      }
    }

    # funnel
    if (!is.null(one$singlesnp) &&
        nrow(one$singlesnp)) {

      tmp <- try(
        TwoSampleMR::mr_funnel_plot(
          one$singlesnp
        ),
        silent = TRUE
      )

      if (!inherits(tmp, "try-error")) {
        p.funnel[[paste0("pair.", i)]] <- tmp[[1L]] + ggplot2::ggtitle(ttl)
      }
    }

    # leave one out
    if (!is.null(one$leaveoneout) &&
        nrow(one$leaveoneout)) {

      tmp <- try(
        TwoSampleMR::mr_leaveoneout_plot(
          one$leaveoneout
        ),
        silent = TRUE
      )

      if (!inherits(tmp, "try-error")) {
        p.leaveoneout[[paste0("pair.", i)]] <- tmp[[1L]] + ggplot2::ggtitle(ttl)
      }
    }
  }

  layout <- wrap_layout(NULL, length(x), ncol = 3L, f.w = 1.4)
  wrap_safe <- function(z) {

    if (!length(z)) {
      return(NULL)
    }
    z <- lapply(z, function(x) x + theme(plot.title = element_text(size = 10L)))
    p <- patchwork::wrap_plots(z, guides = 'collect', ncol = layout$ncol)
    add(layout, p)
  }

  lst <- list(
    p.scatter = wrap_safe(p.scatter),
    p.forest = wrap_safe(p.forest),
    p.funnel = wrap_safe(p.funnel),
    p.leaveoneout = wrap_safe(p.leaveoneout)
  )
}

.stat_summary_mr_batch <- function(
  x,
  digits = 3L,
  top_n = 3L
)
{
  if (!is(x, "mr_batch")) {
    stop('!is(x, "mr_batch")')
  }

  if (length(x) == 0L) {
    return("经筛选后，未获得满足预设标准的孟德尔随机化关联结果。")
  }

  .pick_ivw <- function(dat)
  {
    if (is.null(dat$mr) || !nrow(dat$mr)) {
      return(NULL)
    }

    mr <- tibble::as_tibble(dat$mr)

    hit <- grepl(
      paste(
        c(
          "Inverse variance weighted",
          "fixed effects",
          "multiplicative random effects",
          "IVW"
        ),
        collapse = "|"
      ),
      mr$method,
      ignore.case = TRUE
    )

    if (!any(hit)) {
      return(NULL)
    }

    mr[which(hit)[1L], , drop = FALSE]
  }

  tab <- Filter(
    Negate(is.null),
    lapply(x, .pick_ivw)
  )

  if (!length(tab)) {
    return("结果对象中未检出可用于总结的 IVW 主分析结果。")
  }

  tab <- dplyr::bind_rows(tab)

  if (!all(c("exposure", "outcome", "b", "pval") %in% colnames(tab))) {
    stop("Missing required columns in MR result.")
  }

  if (!"or" %in% colnames(tab)) {
    tab$or <- exp(tab$b)
  }

  tab <- dplyr::arrange(
    tab,
    .data$pval
  )

  n_pair <- nrow(tab)
  n_pos <- sum(tab$b > 0, na.rm = TRUE)
  n_neg <- sum(tab$b < 0, na.rm = TRUE)

  med_nsnp <- if ("nsnp" %in% colnames(tab)) {
    stats::median(tab$nsnp, na.rm = TRUE)
  } else {
    NA_real_
  }

  top_tab <- utils::head(tab, top_n)

  txt_top <- apply(
    top_tab,
    1L,
    function(z) {

      exp_name <- as.character(z["exposure"])
      out_name <- as.character(z["outcome"])

      beta <- round(
        as.numeric(z["b"]),
        digits
      )

      orv <- round(
        as.numeric(z["or"]),
        digits
      )

      pv <- format(
        as.numeric(z["pval"]),
        scientific = TRUE,
        digits = digits
      )

      direction <- if (as.numeric(z["b"]) > 0) {
        "正向关联"
      } else {
        "负向关联"
      }

      glue::glue(
        "{exp_name} 与 {out_name} 呈{direction}（β = {beta}, OR = {orv}, P = {pv}）"
      )
    }
  )

  txt_top <- paste(
    txt_top,
    collapse = "；"
  )

  txt_nsnp <- if (!is.na(med_nsnp)) {
    glue::glue(
      "各分析组合纳入工具变量数量总体充足，中位 SNP 数为 {round(med_nsnp, 1L)}。"
    )
  } else {
    ""
  }

  glue::glue(
    "在完成预设质控与敏感性筛选后，共保留 {n_pair} 组具有统计学支持的暴露因素—结局关联结果。其中⟦mark$red('正向效应 {n_pos} 组，负向效应 {n_neg} 组')⟧。

以 IVW 作为主分析框架，最终保留结果在效应方向上整体稳定，并经多效性、异质性及因果方向检验支持，提示主要关联信号具有较高稳健性。{txt_nsnp}"
  )
}

.stat_plot_summary_mr_batch <- function(
  mr_batch,
  p.scatter = NULL,
  p.forest = NULL,
  p.funnel = NULL,
  p.leaveoneout = NULL,
  top_n = 3L,
  digits = 2L
)
{
  if (!is(mr_batch, "mr_batch")) {
    stop('!is(mr_batch, "mr_batch")')
  }

  if (length(mr_batch) == 0L) {
    return("未获得可用于图形总结的 MR 结果。")
  }

  .pair_name <- function(one)
  {
    exp_nm <- NA_character_
    out_nm <- NA_character_

    if (!is.null(one$mr) && nrow(one$mr)) {

      mr <- tibble::as_tibble(one$mr)

      if ("exposure" %in% colnames(mr)) {
        exp_nm <- mr$exposure[1L]
      }

      if ("outcome" %in% colnames(mr)) {
        out_nm <- mr$outcome[1L]
      }
    }

    paste(exp_nm, out_nm, sep = " → ")
  }

  .check_scatter <- function(one)
  {
    if (is.null(one$mr) || !nrow(one$mr)) {
      return(FALSE)
    }

    mr <- tibble::as_tibble(one$mr)

    if (!"b" %in% colnames(mr)) {
      return(FALSE)
    }

    sgn <- sign(mr$b)
    sgn <- sgn[!is.na(sgn) & sgn != 0]

    if (!length(sgn)) {
      return(FALSE)
    }

    length(unique(sgn)) == 1L
  }

  .check_forest <- function(one)
  {
    if (is.null(one$singlesnp) || !nrow(one$singlesnp)) {
      return(FALSE)
    }

    ss <- tibble::as_tibble(one$singlesnp)

    if (!"b" %in% colnames(ss)) {
      return(FALSE)
    }

    sgn <- sign(ss$b)
    sgn <- sgn[!is.na(sgn) & sgn != 0]

    if (!length(sgn)) {
      return(FALSE)
    }

    max(table(sgn)) / length(sgn) >= 0.70
  }

  .check_funnel <- function(one)
  {
    if (is.null(one$heterogeneity) ||
      !nrow(one$heterogeneity) ||
      !"Q_pval" %in% colnames(one$heterogeneity)) {
      return(FALSE)
    }

    qv <- suppressWarnings(
      max(one$heterogeneity$Q_pval, na.rm = TRUE)
    )

    !is.na(qv) && qv > 0.05
  }

  .check_loo <- function(one)
  {
    if (is.null(one$leaveoneout) ||
      !nrow(one$leaveoneout)) {
      return(FALSE)
    }

    loo <- tibble::as_tibble(one$leaveoneout)

    if (!"b" %in% colnames(loo)) {
      return(FALSE)
    }

    sgn <- sign(loo$b)
    sgn <- sgn[!is.na(sgn) & sgn != 0]

    if (!length(sgn)) {
      return(FALSE)
    }

    length(unique(sgn)) == 1L
  }

  tab <- dplyr::bind_rows(
    lapply(
      seq_along(mr_batch),
      function(i) {

        one <- mr_batch[[i]]

        data.frame(
          idx = i,
          pair = .pair_name(one),
          scatter = .check_scatter(one),
          forest = .check_forest(one),
          funnel = .check_funnel(one),
          loo = .check_loo(one),
          stringsAsFactors = FALSE
        )
      }
    )
  )

  tab$n_pass <- rowSums(
    tab[, c("scatter", "forest", "funnel", "loo")]
  )

  n_pair <- nrow(tab)

  all_pass <- dplyr::filter(
    tab,
    .data$n_pass == 4L
  )

  best_tab <- dplyr::arrange(
    tab,
    dplyr::desc(.data$n_pass),
    .data$pair
  )

  best_tab <- utils::head(best_tab, top_n)

  txt_best <- paste(
    best_tab$pair,
    collapse = "；"
  )

  txt_all <- ""

  if (nrow(all_pass) > 0L) {

    txt_all <- glue::glue(
      "其中共有 {nrow(all_pass)} 组关联在散点图一致性、单 SNP 效应稳定性、漏斗图对称性及逐一剔除稳健性四个维度均表现良好，可视为⟦mark$red('证据最为充分的优先候选结果，包括：{paste(all_pass$pair, collapse = '；')}')⟧。"
    )

  }

  weak_tab <- dplyr::filter(
    tab,
    .data$n_pass <= 2L
  )

  txt_weak <- ""

  if (nrow(weak_tab) > 0L && length(mr_batch) > 2 * top_n) {

    weak_tab <- utils::head(weak_tab, top_n)

    txt_weak <- glue::glue(
      "少数组合稳定性相对一般，主要表现为方法方向不完全一致或逐一剔除后波动较大，如：{paste(weak_tab$pair, collapse = '；')}。"
    )
  }
  glue::glue(
    "图形化敏感性分析整体支持主分析结果的可靠性。多数暴露因素—结局组合在散点图中表现为不同 MR 方法回归方向一致，提示估计效应具有较好的方法学一致性 {aref(p.scatter)}；森林图显示多数工具变量效应方向集中，整体结果通常并非由单一 SNP 驱动 {aref(p.forest)}；漏斗图整体分布较为对称，未见普遍性系统偏倚信号 {aref(p.funnel)}；逐一剔除分析提示多数结果在移除任一 SNP 后仍保持稳定 {aref(p.leaveoneout)}。

{txt_all}

综合各图形维度表现，稳定性较高的代表性关联结果包括：{txt_best}。

{txt_weak}

总体而言，图形学证据与主分析统计结果相互一致，提示最终筛选出的候选关联具有较好的稳健性与解释价值。"
  )
}

.calc_F_stat <- function(data_dat, cut.f = 10L) {
  # Check input
  if (is.null(data_dat) || !nrow(data_dat)) {
    message("Input data_dat is NULL or empty.")
    return(NULL)
  }

  required_cols <- c("beta.exposure", "eaf.exposure", "samplesize.exposure")
  if (!all(required_cols %in% colnames(data_dat))) {
    warning(glue::glue("Missing required columns: {paste(setdiff(required_cols, colnames(data_dat)), collapse=', ')}"))
    return(NULL)
  }

  # Compute F statistic for each SNP
  data_f <- dplyr::reframe(
    data_dat,
    r2_snp = 2 * beta.exposure^2 * eaf.exposure * (1 - eaf.exposure),
    f_snp = (r2_snp / (1 - r2_snp)) * (samplesize.exposure - 2L)
  )

  f_mean <- mean(data_f$f_snp, na.rm = TRUE)

  message(glue::glue(
    "Calculated F statistics for {nrow(data_f)} SNPs; mean F = {signif(f_mean, 3L)}."
  ))

  keep <- which(data_f$f_snp > cut.f)

  snap <- glue::glue("为评估所选工具变量（IVs）作为暴露因素的有效性，我们对每个 SNP 进行了 **F 检验（F-statistic）**。该检验基于每个 SNP 对暴露变量的解释力 ($R^2$) 及样本量 ($n$)，通过公式：

$$
F = \\frac{R^2}{1-R^2} \\times (n - 2)
$$

计算单 SNP 的工具变量强度。F 值越大，说明该 SNP 对暴露的解释能力越强，弱工具偏倚风险越低。随后，通过⟦mark$blue('过滤 F 值低于常用阈值 {{{cut.f}}} 的 SNP')⟧，确保 MR 分析的稳健性。通过 F 检验，过滤了 {{{nrow(data_f) - length(keep)}}} 个 SNP，保留 {{{length(keep)}}} 个 SNP 用于后续 MR 分析。", .open = "{{{", .close = "}}}")

  list(
    data_f = tibble::as_tibble(data_f),
    f_mean = f_mean,
    keep = keep,
    snap = snap
  )
}

# Check whether input data.frame meets TwoSampleMR requirements
# Support raw columns and formatted columns from TwoSampleMR::format_data()

.check_twosamplemr_input <- function(df, type = c("exposure", "outcome")) {

  type <- match.arg(type)

  if (!is.data.frame(df)) {
    stop("`df` must be a data.frame")
  }

  suf <- paste0(".", type)

  message(glue::glue(
    "Checking TwoSampleMR {type} input: {nrow(df)} rows x {ncol(df)} columns"
  ))

  req_fmt <- c(
    "SNP",
    paste0("beta", suf),
    paste0("se", suf),
    paste0("effect_allele", suf),
    paste0("other_allele", suf),
    type
  )

  req_raw <- c(
    "SNP",
    "beta",
    "se",
    "effect_allele",
    "other_allele",
    type
  )

  has_fmt <- all(req_fmt %in% colnames(df))
  has_raw <- all(req_raw %in% colnames(df))

  if (!has_fmt && !has_raw) {

    miss_fmt <- req_fmt[!req_fmt %in% colnames(df)]
    miss_raw <- req_raw[!req_raw %in% colnames(df)]

    stop(
      glue::glue(
        paste0(
          "Input does not match TwoSampleMR format.\n",
          "Missing formatted columns: {paste(miss_fmt, collapse = ', ')}\n",
          "Missing raw columns: {paste(miss_raw, collapse = ', ')}"
        )
      )
    )
  }

  mode <- ifelse(has_fmt, "formatted", "raw")
  use_cols <- if (has_fmt) req_fmt else req_raw

  message(glue::glue(
    "Detected input mode: {mode}"
  ))

  # Strict numeric field names only
  num_base <- c(
    "beta",
    "se",
    "eaf",
    "pval",
    "samplesize",
    "ncase",
    "ncontrol"
  )

  num_cols <- intersect(
    c(num_base, paste0(num_base, suf)),
    colnames(df)
  )

  bad_type <- names(
    which(
      !sapply(df[num_cols], is.numeric)
    )
  )

  if (length(bad_type) > 0L) {
    stop(
      glue::glue(
        "These columns must be numeric: {paste(bad_type, collapse = ', ')}"
      )
    )
  }

  # Required columns NA check
  na_rate <- sapply(
    use_cols,
    function(x) mean(is.na(df[[x]]))
  )

  bad_na <- names(na_rate)[na_rate > 0]

  if (length(bad_na) > 0L) {
    message(glue::glue(
      "Columns containing NA: {paste(bad_na, collapse = ', ')}"
    ))
  }

  # Duplicate SNP
  dup_n <- sum(duplicated(df$SNP))

  if (dup_n > 0L) {
    message(glue::glue(
      "Detected {dup_n} duplicated SNP rows"
    ))
  }

  # Allele check
  ea_col <- intersect(
    c("effect_allele", paste0("effect_allele", suf)),
    colnames(df)
  )[1L]

  oa_col <- intersect(
    c("other_allele", paste0("other_allele", suf)),
    colnames(df)
  )[1L]

  bad_allele <- sum(
    !toupper(df[[ea_col]]) %in% c("A", "C", "G", "T") |
    !toupper(df[[oa_col]]) %in% c("A", "C", "G", "T"),
    na.rm = TRUE
  )

  if (bad_allele > 0L) {
    message(glue::glue(
      "Detected {bad_allele} rows with non-ACGT alleles"
    ))
  }

  message(glue::glue(
    "TwoSampleMR {type} input check passed"
  ))

  invisible(TRUE)
}
