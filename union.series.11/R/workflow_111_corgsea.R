# ==========================================================================
# workflow of corgsea
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_corgsea <- setClass("job_corgsea", 
  contains = c("job"),
  prototype = prototype(
    pg = "corgsea",
    info = c("http://www.gsea-msigdb.org/gsea/downloads.jsp"),
    method = "R package `ClusterProfiler` used for GSEA enrichment",
    cite = "[@ClusterprofilerWuTi2021]",
    tag = "enrich:gsea",
    analysis = "GSEA分析"
    ))

setGeneric("asjob_corgsea",
  function(x, ...) standardGeneric("asjob_corgsea"))

setMethod("asjob_corgsea", signature = c(x = "job_deseq2"),
  function(x, ref, method = "spearman", ...)
  {
    if (x@step < 1L) {
      stop('x@step < 1L.')
    }
    if (any(!ref %in% rownames(object(x)))) {
      stop('any(!ref %in% rownames(object(x))).')
    }
    data <- SummarizedExperiment::assay(x$vst)
    x <- job_corgsea(data, ref, method = method, ...)
    return(x)
  })

setMethod("asjob_corgsea", signature = c(x = "job_limma"),
  function(x, ref, method = "spearman", ...){
    if (x@step < 1L) {
      stop('x@step < 1L.')
    }
    data <- x$normed_data$E
    if (any(!ref %in% rownames(data))) {
      stop('any(!ref %in% rownames(data)).')
    }
    x <- job_corgsea(data, ref, method = method, ...)
    return(x)
  })

job_corgsea <- function(data, ref, method = "spearman",
  rank_strategy = c("auto", "pearson", "spearman", "spearman_pearson_tiebreak"),
  ref_mode = c("individual", "module_score"),
  tie_action = c("switch", "warn", "stop"), min_sd = 0,
  tie_warn_prop = 0.20, tie_stop_prop = 0.50, tiebreak_eps = 1e-8)
{
  method <- match.arg(method, c("pearson", "kendall", "spearman"))
  rank_strategy <- match.arg(rank_strategy)
  ref_mode <- match.arg(ref_mode)
  tie_action <- match.arg(tie_action)

  if (!is(ref, "feature")) {
    stop('!is(ref, "feature").')
  }
  data <- as.matrix(data)
  if (is.null(rownames(data))) {
    stop('is.null(rownames(data)).')
  }
  ref_chr <- as.character(ref)
  missing_ref <- setdiff(ref_chr, rownames(data))
  if (length(missing_ref) > 0L) {
    stop(glue::glue(
      "Reference gene(s) were not found in data rownames: ",
      "{paste(missing_ref, collapse = ', ')}."
    ))
  }
  which.ref <- rownames(data) %in% ref_chr
  if (!any(which.ref)) {
    stop('No reference gene was found in data rownames.')
  }

  vec_sd <- apply(data, 1L, stats::sd, na.rm = TRUE)
  valid_gene <- is.finite(vec_sd) & vec_sd > min_sd
  invalid_ref <- ref_chr[ref_chr %in% rownames(data) & !ref_chr %in% rownames(data)[valid_gene]]
  if (length(invalid_ref) > 0L) {
    stop(glue::glue(
      "Reference gene(s) failed the variance filter: ",
      "{paste(invalid_ref, collapse = ', ')}."
    ))
  }

  data.ref <- data[which.ref & valid_gene, , drop = FALSE]
  data.others <- data[!which.ref & valid_gene, , drop = FALSE]
  if (nrow(data.others) < 10L) {
    stop('nrow(data.others) < 10L.')
  }

  lst_rank <- list()
  lst_diag <- list()

  if (identical(ref_mode, "module_score")) {
    vec_ref <- corgseaFuns$make_module_score(data.ref)
    rank_name <- if (nrow(data.ref) == 1L) rownames(data.ref)[1L] else "Ref_Module"
    res_rank <- corgseaFuns$build_rank(
      vec_ref = vec_ref, mat_other = data.others, rank_name = rank_name,
      method = method, rank_strategy = rank_strategy, tie_action = tie_action,
      tie_warn_prop = tie_warn_prop, tie_stop_prop = tie_stop_prop,
      tiebreak_eps = tiebreak_eps
    )
    lst_rank[[rank_name]] <- res_rank$score
    lst_diag[[rank_name]] <- res_rank$diagnosis
  } else {
    for (i in seq_len(nrow(data.ref))) {
      rank_name <- rownames(data.ref)[i]
      res_rank <- corgseaFuns$build_rank(
        vec_ref = data.ref[i, ], mat_other = data.others, rank_name = rank_name,
        method = method, rank_strategy = rank_strategy, tie_action = tie_action,
        tie_warn_prop = tie_warn_prop, tie_stop_prop = tie_stop_prop,
        tiebreak_eps = tiebreak_eps
      )
      lst_rank[[rank_name]] <- res_rank$score
      lst_diag[[rank_name]] <- res_rank$diagnosis
    }
  }

  data_rank_diagnosis <- do.call(rbind, lst_diag)
  rownames(data_rank_diagnosis) <- NULL
  mat_rank_similarity <- corgseaFuns$rank_similarity(lst_rank)

  x <- .job_corgsea(object = lst_rank)
  x$rank_diagnosis <- data_rank_diagnosis
  x$rank_similarity <- mat_rank_similarity
  x$lst_refine <- list(
    rank_diagnosis = data_rank_diagnosis,
    rank_similarity = mat_rank_similarity
  )
  x@params$cor_method <- method
  x@params$rank_strategy <- rank_strategy
  x@params$ref_mode <- ref_mode
  x@params$tie_action <- tie_action
  x@params$min_sd <- min_sd
  x@params$tie_warn_prop <- tie_warn_prop
  x@params$tie_stop_prop <- tie_stop_prop
  x@params$tiebreak_eps <- tiebreak_eps

  reference_signal_text <- if (identical(ref_mode, "module_score")) {
    "将基因标准化后合成为模块表达信号"
  } else {
    "分别以每个基因的样本间表达变化作为参考信号"
  }
  rank_method_text <- paste(
    unique(corgseaFuns$format_rank_method(data_rank_diagnosis$rank_method_final)),
    collapse = "、"
  )

  x <- methodAdd(
    x,
    glue::glue(
      "以{snap(ref)}相关表达程序为锚点，{reference_signal_text}，计算其与全基因表达谱的相关性，",
      "并据此构建预排序全基因列表。相关性网络分析可用于描述基因表达模式的协同变化 ",
      "(PMID: 19114008)；GSEA 可在全基因连续排序基础上评估预定义基因集是否集中分布于排序列表顶部或底部，",
      "从而识别与目标表达程序协同变化的功能通路 (PMID: 16199517)。",
      "为减少大量并列排序值导致的任意排序影响，本分析对排序列表进行质量评估，",
      "并采用 {rank_method_text} 构建用于富集分析的连续排序统计量。"
    )
  )
  x <- methodAdd(x, corgseaFuns$format_rank_diagnosis_message(data_rank_diagnosis))
  msg_similarity <- corgseaFuns$format_similarity_message(mat_rank_similarity)
  if (!is.null(msg_similarity)) {
    x <- methodAdd(x, msg_similarity)
  }
  return(x)
}

setMethod("step0", signature = c(x = "job_corgsea"),
  function(x){
    step_message("Prepare your data with function `job_corgsea`.")
  })

setMethod("step1", signature = c(x = "job_corgsea"),
  function(x, db, cutoff = .05, pattern = NULL, pvalue = FALSE,
    cutoff.nes = 1, db_anno = NULL, rerun = FALSE,
    min_gs_size = 15L, max_gs_size = 500L, eps = 0,
    collapse_redundant = TRUE, jaccard_cutoff = 0.50,
    mode = c(
      "curated gene sets" = "C2",
      "hallmark gene sets" = "H",
      "positional gene sets" = "C1",
      "regulatory target gene sets" = "C3",
      "computational gene sets" = "C4",
      "ontology gene sets" = "C5",
      "oncogenic signature gene sets" = "C6"
      ), mode_sub = "CP")
  {
    step_message("Custom database for GSEA enrichment.")
    if (missing(db)) {
      mode <- match.arg(mode)
      x <- .set_msig_db(x, mode, mode_sub)
      db <- x$msig_db
    }
    mode_text <- if (length(mode) == 1L) as.character(mode) else "custom gene set"
    if (is.null(db_anno)) {
      db_anno <- x$db_anno
    }
    cli::cli_h1("clusterProfiler::GSEA")
    dir.create("tmp", showWarnings = FALSE, recursive = TRUE)
    use_p <- if (isTRUE(pvalue)) "pvalue" else "p.adjust"

    all.gsea <- pbapply::pbsapply(names(object(x)), simplify = FALSE,
      function(name) {
        glist <- object(x)[[name]]
        args <- list(
          geneList = glist,
          TERM2GENE = db,
          pvalueCutoff = 1,
          pAdjustMethod = "BH",
          minGSSize = min_gs_size,
          maxGSSize = max_gs_size,
          eps = eps,
          seed = TRUE,
          by = "fgsea"
        )
        res.gsea <- expect_local_data(
          "tmp", "gsea", clusterProfiler::GSEA, args, rerun = rerun
        )
        table_gsea_full <- corgseaFuns$format_gsea_table(
          res.gsea, cutoff = cutoff, cutoff.nes = cutoff.nes, use_p = use_p
        )
        if (!is.null(db_anno) && nrow(table_gsea_full) > 0L &&
            all(c("gs_id", "gs_description") %in% colnames(db_anno))) {
          table_gsea_full <- map(
            table_gsea_full, "ID", db_anno, "gs_id", "gs_description",
            col = "Description"
          )
          table_gsea_full <- dplyr::mutate(
            table_gsea_full, Description = stringr::str_wrap(Description, 80)
          )
        }
        if (isTRUE(collapse_redundant)) {
          lst_collapse <- corgseaFuns$collapse_gsea_table(
            table_gsea_full, jaccard_cutoff = jaccard_cutoff
          )
          table_gsea_full <- lst_collapse$full
          table_gsea <- lst_collapse$report
        } else {
          table_gsea <- table_gsea_full
        }
        table_type_text <- if (isTRUE(collapse_redundant)) {
          "代表性通路"
        } else {
          "显著通路"
        }
        table_redundancy_text <- if (isTRUE(collapse_redundant)) {
          "通路代表性依据 leading-edge 基因重叠关系进行归并，以突出相对独立的主要功能方向。"
        } else {
          "结果按统计显著性和标准化富集分数排序展示。"
        }
        table_gsea <- set_lab_legend(
          table_gsea,
          glue::glue("GSEA {table_type_text} of {name} data"),
          glue::glue(
            "基因或基因模块 {name} 的 GSEA {table_type_text}表|||该表展示按 {mode_text} 数据集获得的显著富集通路；",
            "{table_redundancy_text}"
          )
        )
        if (!is.null(db_anno) && nrow(table_gsea) > 0L &&
            all(c("gs_id", "gs_description") %in% colnames(db_anno))) {
          p.gsea <- plot_kegg(table_gsea)
          p.gsea <- .set_lab(
            p.gsea, sig(x), glue::glue("Gene {name} GSEA pathway list of {mode_text}")
          )
          plot_redundancy_text <- if (isTRUE(collapse_redundant)) {
            "该图展示经显著性、效应量和通路相似性筛选后的主要功能方向，"
          } else {
            "该图展示经显著性和效应量筛选后的主要功能方向，"
          }
          p.gsea <- setLegend(
            p.gsea,
            glue::glue(
              "基因或基因模块 {name} 的 GSEA {table_type_text}图。",
              "{plot_redundancy_text}颜色和点大小用于呈现富集显著性与核心基因覆盖情况。"
            )
          )
        } else {
          p.gsea <- NULL
        }
        return(namel(table_gsea, table_gsea_full, p.gsea, res.gsea))
      }
    )

    p.gsea <- lapply(all.gsea, function(x) x$p.gsea)
    table_gsea <- lapply(all.gsea, function(x) x$table_gsea)
    table_gsea_full <- lapply(all.gsea, function(x) x$table_gsea_full)
    if (length(all.gsea) < 6L) {
      snaps <- vapply(names(table_gsea), FUN.VALUE = character(1),
        function(name) {
          corgseaFuns$format_gsea_snap(
            name = name,
            table_report = table_gsea[[name]],
            table_full = table_gsea_full[[name]],
            collapse_redundant = collapse_redundant
          )
        })
      x <- snapAdd(x, paste(snaps, collapse = ""))
    }
    res.gsea <- lapply(all.gsea, function(x) x$res.gsea)
    eps_text <- if (isTRUE(eps == 0)) {
      "极小 P 值估计不设置固定下界，以提高高度显著通路的数值分辨率。"
    } else {
      glue::glue("极小 P 值估计采用 {eps} 作为数值下界，以稳定极端显著结果的统计表示。")
    }
    collapse_text <- if (isTRUE(collapse_redundant)) {
      glue::glue(
        "考虑到 MSigDB 通路集合中存在一定冗余和异质性 (PMID: 26771021)，",
        "显著通路进一步根据 leading-edge 基因的 Jaccard 重叠进行代表性归并，",
        "以突出相对独立的主要功能方向。"
      )
    } else {
      "显著通路未进行相似性归并，结果按统计显著性和标准化富集分数排序展示。"
    }

    x <- methodAdd(
      x,
      glue::glue(
        "使用 {mode_text} 基因集数据库，以 R 包 `clusterProfiler` ⟦pkgInfo('clusterProfiler')⟧ 对预排序全基因列表进行 GSEA ",
        "(PMID: 34557778)。富集结果采用 Benjamini-Hochberg 方法进行多重检验校正，",
        "保留 ⟦mark$blue('{use_p} &lt; {cutoff}，|NES| &gt; {cutoff.nes}')⟧ 的通路；",
        "基因集大小限定为 {min_gs_size}–{max_gs_size} 个基因。",
        "{eps_text}{collapse_text}"
      )
    )
    x@params$res.gsea <- res.gsea
    x@params$db.gsea <- db
    x@params$gsea_filter <- list(
      use_p = use_p,
      cutoff = cutoff,
      cutoff.nes = cutoff.nes,
      min_gs_size = min_gs_size,
      max_gs_size = max_gs_size,
      eps = eps,
      collapse_redundant = collapse_redundant,
      jaccard_cutoff = jaccard_cutoff
    )
    x$db_anno <- db_anno
    if (is.null(x$lst_refine)) {
      x$lst_refine <- list()
    }
    x$lst_refine$gsea_table_full <- table_gsea_full
    x$lst_refine$gsea_table_report <- table_gsea
    x <- tablesAdd(x, table_gsea)
    x <- plotsAdd(x, p.gsea)
    return(x)
  })

setMethod("step2", signature = c(x = "job_corgsea"),
  function(x, top = 10, intersect = TRUE){
    step_message("Select and Visualization")
    x$.feature <- list()
    p.codes <- sapply(names(x$res.gsea), simplify = FALSE, 
      function(name) {
        data <- x@tables$step1$table_gsea[[name]]
        dataTop <- head(data, n = top)
        idTop <- dataTop$ID
        x$.feature[[name]] <<- dataTop$Description
        p.code <- vis(
          x, map = idTop, res.gsea = x$res.gsea[[name]],
          table_gsea = data, .name = name
        )
        p.code
      })
    alls <- names(x$res.gsea)
    snap_ex <- if (x$.args$step1$pvalue) {
      glue::glue(" (按 p.value 排序) ")
    } else {
      glue::glue(" (按 p.adjust 排序) ")
    }
    x <- snapAdd(x, "在 {bind(alls)} 的代表性富集结果中选取 Top {top} 通路进行 GSEA 条码图展示{aref(p.codes)}{snap_ex}。")
    x <- plotsAdd(x, p.codes)
    ins <- lapply(x@tables$step1$table_gsea,
      function(data) {
        head(data$ID, n = top)
      })
    ins <- ins(lst = ins)
    if (length(ins)) {
      insName <- dplyr::filter(x@tables$step1$table_gsea[[1]], ID %in% !!ins)$Description
      x <- snapAdd(x, "不同参考信号的 Top {top} 代表性通路存在 {length(ins)} 个共同功能方向：{bind(insName)}。")
    }
    return(x)
  })

setClassUnion("job_gseaSet", c("job_corgsea", "job_gsea"))

setMethod("vis", signature = c(x = "job_gseaSet"),
  function(x, pattern, map = NULL, res.gsea = NULL, table_gsea = NULL,
    mode = c("kegg", "gsea"), pvalue = FALSE, .name = "", merge = TRUE)
  {
    mode <- match.arg(mode)
    if (is.null(res.gsea)) {
      res.gsea <- x[[ glue::glue("res.{mode}") ]]
    }
    if (is.null(table_gsea)) {
      if (mode == "kegg") {
        table_gsea <- x@tables$step1$table_kegg
      } else if (mode == "gsea") {
        table_gsea <- x@tables$step3$table_gsea
      }
    }
    if (is.null(res.gsea) || is.null(table_gsea)) {
      stop('is.null(res.gsea) || is.null(table_gsea).')
    }
    alls <- table_gsea$ID
    if (is.null(alls)) {
      stop('is.null(alls).')
    }
    if (is.null(map)) {
      whichMapped <- which(grepl(pattern, table_gsea$Description, ignore.case = TRUE))
      map <- alls[ whichMapped ]
    } else {
      whichMapped <- which(table_gsea$ID %in% map)
    }
    if (!length(map)) {
      message(crayon::red("Not match any pathway, skip plot of 'p.code'."))
      p.code <- NULL
    } else {
      if (merge) {
        fun_plot <- function() {
          res.gsea@result <- map(
            res.gsea@result, "ID", table_gsea, "ID", "Description", col = "Description"
          )
          plst <- enrichplot::gseaplot2(res.gsea, map, pvalue_table = pvalue)
          plst[[1]] <- plst[[1]] +
            scale_color_manual(values = color_set(TRUE)) +
            guides(color = guide_legend(ncol = 2)) +
            theme(legend.position = "top")
          print(plst)
        }
        p.code <- grid.grabExpr(fun_plot())
      } else {
        p.code <- sapply(map, simplify = FALSE,
          function(key) {
            title <- dplyr::filter(table_gsea, ID == key)$Description
            grob <- grid.grabExpr(
              print(enrichplot::gseaplot2(res.gsea, key, pvalue_table = pvalue, title = title))
            )
            wrap(grob, 5, 4)
          })
        if (length(map) > 1) {
          layout <- calculate_layout(length(map))
          p.code <- patchwork::wrap_plots(
            lapply(p.code, function(x) x@data), ncol = layout[["cols"]]
          )
          p.code <- wrap_layout(p.code, layout, 3)
        } else {
          p.code <- p.code[[1]]
        }
      }
    }
    ids <- table_gsea$ID[whichMapped]
    des <- table_gsea$Description[whichMapped]
    p.code <- set_lab_legend(
      p.code,
      glue::glue("{sig(x)} GSEA plot {.name}"),
      glue::glue("GSEA 富集条码图 ({.name})|||第一部分为 ES 折线图，峰值位置反映基因集在排序列表中的富集方向与富集强度；正值表示基因集倾向于在列表顶部富集，负值表示基因集倾向于在列表底部富集。第二部分为基因集成员位置图，用竖线标记基因集成员在全基因排序列表中的位置。第三部分为排序后全基因相关统计量的分布，左侧正值表示与参考表达信号正相关，右侧负值表示与参考表达信号负相关。")
    )
    p.code
  })

# ==========================================================================
# helper functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if (!exists("corgseaFuns")) {
  corgseaFuns <- new.env(parent = emptyenv())
}

corgseaFuns$diagnose_rank <- function(vec_score)
{
  vec_score <- vec_score[is.finite(vec_score)]
  n_gene <- length(vec_score)
  if (n_gene == 0L) {
    return(data.frame(
      n_gene = 0L, n_unique = 0L, unique_prop = NA_real_, tie_prop = NA_real_,
      max_tie_block_prop = NA_real_, n_positive = 0L, n_negative = 0L,
      n_zero = 0L, stringsAsFactors = FALSE
    ))
  }
  n_unique <- length(unique(vec_score))
  tab_score <- table(vec_score)
  data.frame(
    n_gene = n_gene,
    n_unique = n_unique,
    unique_prop = n_unique / n_gene,
    tie_prop = 1 - n_unique / n_gene,
    max_tie_block_prop = max(tab_score) / n_gene,
    n_positive = sum(vec_score > 0),
    n_negative = sum(vec_score < 0),
    n_zero = sum(vec_score == 0),
    stringsAsFactors = FALSE
  )
}

corgseaFuns$score_by_cor <- function(vec_ref, mat_other, method)
{
  mat_other <- as.matrix(mat_other)
  vec_score <- stats::cor(
    as.numeric(vec_ref), t(mat_other), method = method,
    use = "pairwise.complete.obs"
  )
  vec_score <- as.numeric(vec_score)
  names(vec_score) <- rownames(mat_other)
  vec_score
}

corgseaFuns$zscore <- function(vec_x)
{
  vec_x <- as.numeric(vec_x)
  sd_x <- stats::sd(vec_x, na.rm = TRUE)
  if (!is.finite(sd_x) || sd_x == 0) {
    return(rep(0, length(vec_x)))
  }
  vec_z <- (vec_x - mean(vec_x, na.rm = TRUE)) / sd_x
  vec_z[!is.finite(vec_z)] <- 0
  vec_z
}

corgseaFuns$build_rank <- function(vec_ref, mat_other, rank_name, method,
  rank_strategy, tie_action, tie_warn_prop, tie_stop_prop, tiebreak_eps)
{
  primary_method <- method
  if (identical(rank_strategy, "pearson")) {
    primary_method <- "pearson"
  } else if (identical(rank_strategy, "spearman") ||
      identical(rank_strategy, "spearman_pearson_tiebreak")) {
    primary_method <- "spearman"
  }

  vec_primary <- corgseaFuns$score_by_cor(vec_ref, mat_other, primary_method)
  vec_score <- vec_primary
  rank_method_final <- primary_method
  rank_reason <- "primary"

  data_diag_initial <- corgseaFuns$diagnose_rank(vec_primary)

  if (identical(rank_strategy, "auto") &&
      is.finite(data_diag_initial$tie_prop) &&
      data_diag_initial$tie_prop > tie_warn_prop &&
      !identical(primary_method, "pearson")) {
    vec_score <- corgseaFuns$score_by_cor(vec_ref, mat_other, "pearson")
    rank_method_final <- "pearson"
    rank_reason <- "auto_switch_from_high_ties"
  }

  if (identical(rank_strategy, "spearman_pearson_tiebreak")) {
    vec_pearson <- corgseaFuns$score_by_cor(vec_ref, mat_other, "pearson")
    genes_common <- intersect(names(vec_score), names(vec_pearson))
    vec_score <- vec_score[genes_common] +
      tiebreak_eps * corgseaFuns$zscore(vec_pearson[genes_common])
    rank_method_final <- "spearman_plus_pearson_tiebreak"
    rank_reason <- "deterministic_tiebreak"
  }

  vec_score <- vec_score[is.finite(vec_score)]
  data_diag_final <- corgseaFuns$diagnose_rank(vec_score)

  if (is.finite(data_diag_final$tie_prop) &&
      data_diag_final$tie_prop > tie_stop_prop &&
      identical(tie_action, "stop")) {
    stop(glue::glue(
      "Rank list of {rank_name} has too many tied values ",
      "({round(data_diag_final$tie_prop * 100, 2)}%)."
    ))
  }

  if (is.finite(data_diag_final$tie_prop) &&
      data_diag_final$tie_prop > tie_warn_prop) {
    warning(glue::glue(
      "Rank list of {rank_name} still has a high tied-rank proportion ",
      "({round(data_diag_final$tie_prop * 100, 2)}%)."
    ))
  }

  data_diag <- data.frame(
    rank_name = rank_name,
    rank_strategy = rank_strategy,
    rank_method_initial = primary_method,
    rank_method_final = rank_method_final,
    rank_reason = rank_reason,
    tie_action = tie_action,
    data_diag_final,
    stringsAsFactors = FALSE
  )

  list(score = sort(vec_score, decreasing = TRUE), diagnosis = data_diag)
}

corgseaFuns$make_module_score <- function(mat_ref)
{
  mat_ref <- as.matrix(mat_ref)
  mat_z <- t(scale(t(mat_ref)))
  mat_z[!is.finite(mat_z)] <- 0
  colMeans(mat_z, na.rm = TRUE)
}

corgseaFuns$rank_similarity <- function(lst_rank)
{
  if (length(lst_rank) < 2L) {
    return(NULL)
  }
  genes_common <- Reduce(intersect, lapply(lst_rank, names))
  if (length(genes_common) < 3L) {
    return(NULL)
  }
  mat_rank <- vapply(lst_rank, FUN.VALUE = numeric(length(genes_common)),
    function(vec_score) {
      vec_score[genes_common]
    })
  stats::cor(mat_rank, method = "spearman", use = "pairwise.complete.obs")
}

corgseaFuns$format_rank_method <- function(vec_method)
{
  vec_method <- as.character(vec_method)
  vec_method[vec_method == "pearson"] <- "Pearson 相关系数"
  vec_method[vec_method == "spearman"] <- "Spearman 秩相关系数"
  vec_method[vec_method == "kendall"] <- "Kendall 秩相关系数"
  vec_method[vec_method == "spearman_plus_pearson_tiebreak"] <-
    "Spearman 秩相关系数结合 Pearson 次级排序"
  vec_method
}

corgseaFuns$format_rank_diagnosis_message <- function(data_diag)
{
  if (is.null(data_diag) || nrow(data_diag) == 0L) {
    return("未生成有效的全基因排序质量评估结果。")
  }
  tie_median <- round(stats::median(data_diag$tie_prop, na.rm = TRUE) * 100, 2)
  tie_max <- round(max(data_diag$tie_prop, na.rm = TRUE) * 100, 2)
  methods <- paste(
    unique(corgseaFuns$format_rank_method(data_diag$rank_method_final)),
    collapse = "、"
  )
  glue::glue(
    "全基因排序质量评估显示，{nrow(data_diag)} 个排序列表的并列排序值比例中位数为 ",
    "{tie_median}%，最高为 {tie_max}%；最终用于富集分析的相关统计量包括 {methods}。"
  )
}

corgseaFuns$format_similarity_message <- function(mat_similarity)
{
  if (is.null(mat_similarity) || ncol(mat_similarity) < 2L) {
    return(NULL)
  }
  vec_upper <- mat_similarity[upper.tri(mat_similarity)]
  vec_upper <- vec_upper[is.finite(vec_upper)]
  if (!length(vec_upper)) {
    return(NULL)
  }
  med_abs <- round(stats::median(abs(vec_upper), na.rm = TRUE), 3)
  max_abs <- round(max(abs(vec_upper), na.rm = TRUE), 3)
  glue::glue(
    "不同参考基因得到的全基因排序列表之间 Spearman 相似度中位数为 {med_abs}，",
    "最高为 {max_abs}；若相似度过高，说明这些参考基因主要反映共同表达程序，",
    "更适合从基因模块层面进行整体功能解释。"
  )
}

corgseaFuns$split_core_genes <- function(vec_core)
{
  vec_core <- as.character(vec_core)
  vec_core[is.na(vec_core)] <- ""
  lapply(strsplit(vec_core, "/", fixed = TRUE), function(x) {
    x <- x[nzchar(x)]
    unique(x)
  })
}

corgseaFuns$format_gsea_table <- function(res.gsea, cutoff, cutoff.nes,
  use_p = "p.adjust")
{
  table_gsea <- dplyr::as_tibble(res.gsea@result)
  if (nrow(table_gsea) == 0L) {
    return(table_gsea)
  }
  if (use_p %in% colnames(table_gsea)) {
    table_gsea <- table_gsea[table_gsea[[use_p]] < cutoff, , drop = FALSE]
  }
  if ("NES" %in% colnames(table_gsea)) {
    table_gsea <- table_gsea[abs(table_gsea$NES) > cutoff.nes, , drop = FALSE]
  }
  table_gsea <- dplyr::as_tibble(table_gsea)
  if (nrow(table_gsea) == 0L) {
    return(table_gsea)
  }
  if ("core_enrichment" %in% colnames(table_gsea)) {
    geneName_list <- corgseaFuns$split_core_genes(table_gsea$core_enrichment)
  } else {
    geneName_list <- vector("list", nrow(table_gsea))
  }
  leading_pct <- rep(NA_real_, nrow(table_gsea))
  if ("leading_edge" %in% colnames(table_gsea)) {
    leading_pct <- suppressWarnings(as.numeric(
      stringr::str_extract(table_gsea$leading_edge, "[0-9]+")
    ))
    leading_pct <- round(leading_pct / 100, 2)
  }
  table_gsea$geneName_list <- geneName_list
  table_gsea$Count <- lengths(geneName_list)
  table_gsea$GeneRatio <- leading_pct
  table_gsea
}

corgseaFuns$collapse_gsea_table <- function(table_gsea, jaccard_cutoff = 0.5)
{
  table_gsea <- dplyr::as_tibble(table_gsea)
  if (nrow(table_gsea) == 0L) {
    table_gsea$is_redundant <- logical(0)
    table_gsea$representative_ID <- character(0)
    table_gsea$representative_Description <- character(0)
    table_gsea$leading_edge_jaccard <- numeric(0)
    return(list(full = table_gsea, report = table_gsea))
  }
  if (!"geneName_list" %in% colnames(table_gsea)) {
    table_gsea$geneName_list <- corgseaFuns$split_core_genes(table_gsea$core_enrichment)
  }
  ord <- order(table_gsea$p.adjust, -abs(table_gsea$NES), table_gsea$ID, na.last = TRUE)
  keep_idx <- integer(0)
  is_redundant <- rep(FALSE, nrow(table_gsea))
  representative_ID <- rep(NA_character_, nrow(table_gsea))
  representative_Description <- rep(NA_character_, nrow(table_gsea))
  leading_edge_jaccard <- rep(NA_real_, nrow(table_gsea))

  for (idx in ord) {
    genes_i <- unique(table_gsea$geneName_list[[idx]])
    direction_i <- sign(table_gsea$NES[idx])
    best_jaccard <- 0
    best_idx <- NA_integer_
    if (length(keep_idx) > 0L && length(genes_i) > 0L) {
      for (idx_keep in keep_idx) {
        direction_keep <- sign(table_gsea$NES[idx_keep])
        if (!identical(direction_i, direction_keep)) {
          next
        }
        genes_keep <- unique(table_gsea$geneName_list[[idx_keep]])
        n_union <- length(union(genes_i, genes_keep))
        if (n_union == 0L) {
          next
        }
        jaccard <- length(intersect(genes_i, genes_keep)) / n_union
        if (jaccard > best_jaccard) {
          best_jaccard <- jaccard
          best_idx <- idx_keep
        }
      }
    }
    if (is.finite(best_jaccard) && best_jaccard >= jaccard_cutoff) {
      is_redundant[idx] <- TRUE
      leading_edge_jaccard[idx] <- best_jaccard
      representative_ID[idx] <- table_gsea$ID[best_idx]
      representative_Description[idx] <- table_gsea$Description[best_idx]
    } else {
      keep_idx <- c(keep_idx, idx)
      leading_edge_jaccard[idx] <- NA_real_
      representative_ID[idx] <- table_gsea$ID[idx]
      representative_Description[idx] <- table_gsea$Description[idx]
    }
  }

  table_gsea$is_redundant <- is_redundant
  table_gsea$representative_ID <- representative_ID
  table_gsea$representative_Description <- representative_Description
  table_gsea$leading_edge_jaccard <- leading_edge_jaccard
  table_report <- table_gsea[!table_gsea$is_redundant, , drop = FALSE]
  list(full = dplyr::as_tibble(table_gsea), report = dplyr::as_tibble(table_report))
}

corgseaFuns$format_gsea_snap <- function(name, table_report, table_full,
  collapse_redundant)
{
  n_full <- nrow(table_full)
  n_report <- nrow(table_report)
  if (n_full == 0L) {
    return(glue::glue("{name} 未获得符合阈值的显著 GSEA 通路。"))
  }
  top_terms <- character(0)
  if ("Description" %in% colnames(table_report) && n_report > 0L) {
    top_terms <- head(as.character(table_report$Description), 3L)
    top_terms <- top_terms[nzchar(top_terms)]
  }
  top_text <- ""
  if (length(top_terms) > 0L) {
    top_text <- glue::glue("代表性通路包括 {paste(top_terms, collapse = '、')}。")
  }
  if (isTRUE(collapse_redundant)) {
    return(glue::glue(
      "{name} 获得 {n_full} 个显著富集通路，经 leading-edge 基因重叠归并后形成 ",
      "{n_report} 个代表性功能方向。{top_text}"
    ))
  }
  glue::glue("{name} 获得 {n_full} 个显著 GSEA 通路。{top_text}")
}
