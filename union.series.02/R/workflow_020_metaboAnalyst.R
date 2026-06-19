# ==========================================================================
# workflow of MetaboAnalyst enrichment analysis
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_metaboAnalyst <- setClass("job_metaboAnalyst",
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "metaboAnalyst",
    info = c("Tutorial: https://www.metaboanalyst.ca"),
    cite = "[@Metaboanalyst4Chong2018]",
    method = "`MetaboAnalystR` was used for metabolite set enrichment analysis.",
    tag = "metaboAnalyst",
    analysis = "MetaboAnalyst 代谢物富集分析"
    ))

job_metaboAnalyst <- function(cpds = NULL,
  cmpd_type = "name",
  libname = "kegg_pathway",
  source_feature = NULL)
{
  x <- .job_metaboAnalyst()
  x$cpds <- metaboAnalystFuns$clean_cpds(cpds)
  x$source_feature <- source_feature
  x$params <- list(
    cmpd_type = cmpd_type,
    libname = libname
  )
  x$lst_refine <- list()
  return(x)
}

setGeneric("asjob_metaboAnalyst",
  function(x, ...) standardGeneric("asjob_metaboAnalyst"))

setMethod("asjob_metaboAnalyst", signature = c(x = "feature"),
  function(x,
    unlist = TRUE,
    cmpd_type = "name",
    libname = "kegg_pathway",
    ...)
  {
    feature_raw <- x
    names_x <- names(x)
    x <- resolve_feature_snapAdd_onExit("x", x)

    if (isTRUE(unlist)) {
      cpds <- unique(as.character(unlist(x, use.names = FALSE)))
      cpds <- cpds[!is.na(cpds) & cpds != ""]
    } else {
      names(x) <- names_x
      cpds <- x
    }

    x <- job_metaboAnalyst(
      cpds = cpds,
      cmpd_type = cmpd_type,
      libname = libname,
      source_feature = feature_raw
    )

    return(x)
  })

setMethod("asjob_metaboAnalyst", signature = c(x = "job_metaboDiff"),
  function(x, ...)
  {
    if (is.null(x$.feature)) {
      stop("`x$.feature` was not found. Please run `step3()` first.")
    }

    asjob_metaboAnalyst(
      x$.feature,
      ...
    )
  })

setMethod("step0", signature = c(x = "job_metaboAnalyst"),
  function(x)
  {
    step_message("Prepare your data with function `job_metaboAnalyst` or `asjob_metaboAnalyst`.")
  })

setMethod("step1", signature = c(x = "job_metaboAnalyst"),
  function(x,
    cpds = NULL,
    cmpd_type = NULL,
    libname = NULL,
    skip = FALSE,
    plot = TRUE,
    plot_top_n = 10L,
    plot_sort_by = c("pvalue", "fdr", "impact"),
    plot_prefix = "metabolites_ORA_dot_",
    cache_dir = NULL,
    cache_path = "tmp",
    default_dpi = 72L,
    ora_p_cutoff = 0.05,
    ora_fdr_cutoff = NULL,
    impact_cutoff = NULL,
    skip_failed_lib = TRUE,
    report_library_status = FALSE,
    run_metaboanalyst_pdf = FALSE,
    verbose = TRUE)
  {
    step_message("MetaboAnalyst metabolite set enrichment analysis.")

    if (isTRUE(skip)) {
      message("Skip MetaboAnalyst metabolite set enrichment analysis.")
      return(x)
    }

    plot_sort_by <- match.arg(plot_sort_by)

    if (is.null(cpds)) {
      cpds <- x$cpds
    }
    cpds <- metaboAnalystFuns$clean_cpds(cpds)
    if (length(cpds) == 0L) {
      stop("`cpds` is empty. Please provide compounds or use `asjob_metaboAnalyst(feature)`.")
    }

    if (is.null(cmpd_type)) {
      cmpd_type <- x$params$cmpd_type
    }
    if (is.null(libname)) {
      libname <- x$params$libname
    }
    libname <- unique(as.character(libname))
    libname <- libname[!is.na(libname) & libname != ""]
    if (length(libname) == 0L) {
      stop("`libname` is empty.")
    }

    if (is.null(cache_dir)) {
      cache_dir <- create_job_cache_dir(
        x,
        name = "metaboAnalyst",
        path = cache_path
      )
    }
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    cache_dir <- normalizePath(cache_dir, winslash = "/", mustWork = TRUE)

    lst_res <- lapply(libname, function(lib_i) {
      metaboAnalystFuns$run_ora_safe(
        cpds = cpds,
        cmpd_type = cmpd_type,
        libname = lib_i,
        plot = isTRUE(run_metaboanalyst_pdf),
        plot_prefix = plot_prefix,
        plot_dir = cache_dir,
        default_dpi = default_dpi,
        verbose = verbose
      )
    })
    names(lst_res) <- libname

    data_status <- metaboAnalystFuns$get_library_status(lst_res)
    lst_success <- lst_res[vapply(lst_res, function(z) isTRUE(z$success), logical(1L))]

    x$dir_cache <- cache_dir
    x$lst_refine$metaboAnalyst_library_status <- data_status

    if (length(lst_success) == 0L) {
      warning("No MetaboAnalyst enrichment result was obtained from the selected metabolite set libraries.")

      if (isTRUE(report_library_status)) {
        t.metaboAnalyst_library_status <- set_lab_legend(
          data_status,
          glue::glue("{x@sig} MetaboAnalyst library status"),
          glue::glue("MetaboAnalyst 数据库运行状态|||该表为内部诊断表，展示不同 MetaboAnalyst 代谢物集合库的运行状态和错误信息。")
        )

        x <- tablesAdd(
          x,
          t.metaboAnalyst_library_status = t.metaboAnalyst_library_status
        )
      }

      if (!isTRUE(skip_failed_lib)) {
        stop("All MetaboAnalyst libraries failed.")
      }

      return(x)
    }

    x$mSet <- lapply(lst_success, function(z) z$result$mSet)
    x$hits <- lapply(lst_success, function(z) z$result$hits)
    x$lst_refine$metaboAnalyst_ora_by_lib <- lapply(lst_success, function(z) z$result)

    data_mapped <- metaboAnalystFuns$combine_result_tables(
      lst_success,
      table_name = "data_mapped"
    )
    data_ora <- metaboAnalystFuns$combine_result_tables(
      lst_success,
      table_name = "data_ora"
    )
    data_enrich <- metaboAnalystFuns$standardize_enrichment_table(
      data_ora,
      p_cutoff = ora_p_cutoff,
      fdr_cutoff = ora_fdr_cutoff,
      impact_cutoff = impact_cutoff
    )

    x$lst_refine$metaboAnalyst_mapping_all <- data_mapped
    x$lst_refine$metaboAnalyst_ora_all <- data_ora
    x$lst_refine$metaboAnalyst_enrichment_all <- data_enrich

    n_input <- length(cpds)
    n_mapped <- if (nrow(data_mapped) == 0L) 0L else length(unique(data_mapped$Query))
    n_ora <- nrow(data_ora)

    enrich_stats <- metaboAnalystFuns$summarize_enrichment_table(
      data_enrich,
      p_cutoff = ora_p_cutoff,
      fdr_cutoff = ora_fdr_cutoff,
      impact_cutoff = impact_cutoff
    )

    vec_success_lib <- names(lst_success)
    text_lib <- paste(vec_success_lib, collapse = "、")
    text_lib_method <- if (length(vec_success_lib) == 1L) {
      glue::glue("选择 {text_lib} 作为代谢物集合库")
    } else {
      glue::glue("分别选择 {text_lib} 作为代谢物集合库")
    }

    text_fdr_method <- if (is.null(ora_fdr_cutoff)) {
      ""
    } else {
      glue::glue("富集结果表同步列示多重检验校正后的 FDR 或 adjusted P value，并以 FDR &lt; {ora_fdr_cutoff} 作为显著性筛选阈值。")
    }

    x <- methodAdd(x, glue::glue(
      "基于差异代谢物列表，采用 R 包 `MetaboAnalystR` ⟦pkgInfo('MetaboAnalystR')⟧ 进行代谢物集合富集分析。输入代谢物首先根据 `{cmpd_type}` 类型进行名称映射和交叉引用，随后{text_lib_method}，并采用超几何检验评估输入代谢物在各代谢物集合中的富集情况。富集显著性主要参考 P value &lt; {ora_p_cutoff}。{text_fdr_method}"
    ))

    if (isTRUE(report_library_status)) {
      t.metaboAnalyst_library_status <- set_lab_legend(
        data_status,
        glue::glue("{x@sig} MetaboAnalyst library status"),
        glue::glue("MetaboAnalyst 数据库运行状态|||该表为内部诊断表，展示不同 MetaboAnalyst 代谢物集合库的运行状态、映射数量、富集结果数量以及错误信息。")
      )

      x <- tablesAdd(
        x,
        t.metaboAnalyst_library_status = t.metaboAnalyst_library_status
      )
    }

    if (nrow(data_mapped) > 0L) {
      t.metaboAnalyst_mapping <- set_lab_legend(
        data_mapped,
        glue::glue("{x@sig} MetaboAnalyst mapped metabolites"),
        glue::glue("MetaboAnalyst 代谢物映射结果|||该表展示输入差异代谢物在 MetaboAnalystR 中的名称映射和数据库交叉引用结果。")
      )

      x <- tablesAdd(
        x,
        t.metaboAnalyst_mapping = t.metaboAnalyst_mapping
      )
    }

    if (nrow(data_enrich) > 0L) {
      t.metaboAnalyst_enrichment <- set_lab_legend(
        data_enrich,
        glue::glue("{x@sig} MetaboAnalyst enrichment result"),
        glue::glue("MetaboAnalyst 富集分析结果|||该表展示基于 {text_lib} 代谢物集合库获得的 ORA 富集结果。ORA 采用超几何检验评估输入代谢物在各代谢物集合中的富集情况，结果表保留原始 P value 和多重检验校正指标。")
      )

      x <- tablesAdd(
        x,
        t.metaboAnalyst_enrichment = t.metaboAnalyst_enrichment
      )
    } else {
      warning("No MetaboAnalyst ORA result was obtained from successful libraries.")
    }

    if (isTRUE(plot) && nrow(data_enrich) > 0L) {
      p.metaboAnalyst_enrichment <- metaboAnalystFuns$plot_enrichment_facet(
        data_enrich,
        top_n = plot_top_n,
        sort_by = plot_sort_by,
        p_cutoff = ora_p_cutoff,
        fdr_cutoff = ora_fdr_cutoff,
        impact_cutoff = impact_cutoff
      )

      text_plot_facet <- if (length(unique(data_enrich$libname)) == 1L) {
        "该图展示代谢物集合富集结果。"
      } else {
        "该图按代谢物集合库分面展示富集结果。"
      }

      p.metaboAnalyst_enrichment <- set_lab_legend(
        wrap(p.metaboAnalyst_enrichment, 6, 5),
        glue::glue("{x@sig} MetaboAnalyst enrichment dot plot"),
        glue::glue("MetaboAnalyst 富集气泡图|||{text_plot_facet}横轴为 -log10(P value)，纵轴为富集代谢物集合名称；点大小表示富集倍数（Enrichment Ratio，按命中数量/期望命中数量计算；若无法计算则以命中代谢物数量显示）；颜色表示 P value，颜色越偏红代表 P value 越小。竖向虚线表示 P value = {ora_p_cutoff} 的显著性参考阈值。")
      )

      x <- plotsAdd(
        x,
        p.metaboAnalyst_enrichment = p.metaboAnalyst_enrichment
      )
    }

    text_p_sig <- if (is.na(enrich_stats$n_p_sig)) {
      "未识别到可用于 P value 统计的列"
    } else {
      glue::glue("{enrich_stats$n_p_sig} 条结果满足 P value &lt; {ora_p_cutoff}")
    }

    text_fdr_sig <- if (is.null(ora_fdr_cutoff) || is.na(enrich_stats$n_fdr_sig)) {
      ""
    } else {
      glue::glue("，{enrich_stats$n_fdr_sig} 条结果满足 FDR &lt; {ora_fdr_cutoff}")
    }

    text_final_sig <- if (is.na(enrich_stats$n_final_sig)) {
      ""
    } else if (!is.null(ora_fdr_cutoff)) {
      glue::glue("；综合设定阈值后共有 {enrich_stats$n_final_sig} 条结果满足筛选条件")
    } else {
      ""
    }

    x <- snapAdd(x, glue::glue(
      "共输入 {n_input} 个差异代谢物用于 MetaboAnalystR 富集分析，在 {text_lib} 代谢物集合库中共有 {n_mapped} 个代谢物完成数据库映射，并获得 {n_ora} 条代谢物集合富集结果；{text_p_sig}{text_fdr_sig}{text_final_sig}。"
    ))

    return(x)
  })

setMethod("step2", signature = c(x = "job_metaboAnalyst"),
  function(x,
    cpds = NULL,
    cmpd_type = NULL,
    organism = "hsa",
    organism_label = "Human",
    node_imp = c("rbc", "dgr"),
    enrichment_method = c("hyperg", "fisher"),
    p_cutoff = 0.05,
    fdr_cutoff = 0.05,
    impact_cutoff = 0.1,
    plot = TRUE,
    plot_top_n = 25L,
    cache_dir = NULL,
    cache_path = "tmp",
    default_dpi = 72L,
    skip_failed = TRUE,
    verbose = TRUE)
  {
    step_message("MetaboAnalyst pathway topology analysis.")

    node_imp <- match.arg(node_imp)
    enrichment_method <- match.arg(enrichment_method)

    if (is.null(cpds)) {
      cpds <- x$cpds
    }
    cpds <- metaboAnalystFuns$clean_cpds(cpds)
    if (length(cpds) == 0L) {
      stop("`cpds` is empty. Please provide compounds or use `asjob_metaboAnalyst(feature)`.")
    }

    if (is.null(cmpd_type)) {
      cmpd_type <- x$params$cmpd_type
    }

    if (is.null(cache_dir)) {
      if (!is.null(x$dir_cache) && dir.exists(x$dir_cache)) {
        cache_dir <- x$dir_cache
      } else {
        cache_dir <- create_job_cache_dir(
          x,
          name = "metaboAnalyst",
          path = cache_path
        )
      }
    }
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    cache_dir <- normalizePath(cache_dir, winslash = "/", mustWork = TRUE)
    x$dir_cache <- cache_dir

    res_pathway <- metaboAnalystFuns$run_pathway_safe(
      cpds = cpds,
      cmpd_type = cmpd_type,
      organism = organism,
      node_imp = node_imp,
      enrichment_method = enrichment_method,
      default_dpi = default_dpi,
      work_dir = cache_dir,
      verbose = verbose
    )

    x$lst_refine$metaboAnalyst_pathway_status <- data.frame(
      organism = organism,
      success = isTRUE(res_pathway$success),
      error = ifelse(isTRUE(res_pathway$success), NA_character_, res_pathway$error),
      stringsAsFactors = FALSE
    )

    if (!isTRUE(res_pathway$success)) {
      warning(glue::glue("MetaboAnalyst pathway topology analysis failed: {res_pathway$error}"))
      if (!isTRUE(skip_failed)) {
        stop(res_pathway$error)
      }
      return(x)
    }

    res_pathway <- res_pathway$result
    data_mapped <- res_pathway$data_mapped
    data_pathway_raw <- res_pathway$data_pathway
    data_pathway <- metaboAnalystFuns$standardize_enrichment_table(
      data_pathway_raw,
      p_cutoff = p_cutoff,
      fdr_cutoff = fdr_cutoff,
      impact_cutoff = impact_cutoff
    )

    data_pathway <- data_pathway[order(data_pathway$pvalue, -data_pathway$impact), , drop = FALSE]

    x$mSet_pathway <- res_pathway$mSet
    x$lst_refine$metaboAnalyst_pathway_analysis <- res_pathway
    x$lst_refine$metaboAnalyst_pathway_mapping <- data_mapped
    x$lst_refine$metaboAnalyst_pathway_result <- data_pathway

    n_input <- length(cpds)
    n_mapped <- if (nrow(data_mapped) == 0L) 0L else length(unique(data_mapped$Query))
    n_result <- nrow(data_pathway)
    n_final <- sum(data_pathway$significant, na.rm = TRUE)

    data_sig <- data_pathway[data_pathway$significant, , drop = FALSE]
    text_sig_path <- if (nrow(data_sig) == 0L) {
      "未获得同时满足 FDR 和 pathway impact 阈值的显著扰动通路"
    } else {
      data_sig <- data_sig[order(data_sig$fdr, -data_sig$impact, data_sig$pvalue), , drop = FALSE]
      paste(head(data_sig$pathway_name, 10L), collapse = "、")
    }

    x <- methodAdd(x, glue::glue(
      "为进一步评估差异代谢物涉及的代谢通路扰动，采用 R 包 `MetaboAnalystR` ⟦pkgInfo('MetaboAnalystR')⟧ 进行 pathway analysis。输入代谢物根据 `{cmpd_type}` 类型进行名称映射和交叉引用，并选择 {organism_label}（{organism}）KEGG pathway library 作为通路背景。通路富集采用 {enrichment_method} 检验进行 over-representation analysis，通路拓扑影响值采用 {node_imp} 方法计算。根据分析方案，将 FDR &lt; {fdr_cutoff} 且 pathway impact &gt; {impact_cutoff} 的通路定义为显著扰动代谢通路。"
    ))

    if (nrow(data_mapped) > 0L) {
      t.metaboAnalyst_pathway_mapping <- set_lab_legend(
        data_mapped,
        glue::glue("{x@sig} MetaboAnalyst pathway mapped metabolites"),
        glue::glue("MetaboAnalyst pathway analysis 代谢物映射结果|||该表展示输入差异代谢物在 MetaboAnalystR pathway analysis 中的名称映射和数据库交叉引用结果。")
      )

      x <- tablesAdd(
        x,
        t.metaboAnalyst_pathway_mapping = t.metaboAnalyst_pathway_mapping
      )
    }

    if (nrow(data_pathway) > 0L) {
      t.metaboAnalyst_pathway <- set_lab_legend(
        data_pathway,
        glue::glue("{x@sig} MetaboAnalyst pathway topology result"),
        glue::glue("MetaboAnalyst pathway analysis 结果|||该表展示基于 {organism_label} KEGG pathway library 的通路富集和拓扑分析结果。P value 来源于 over-representation analysis，FDR 为多重检验校正后的显著性指标，pathway impact 表示匹配代谢物在通路拓扑结构中的相对影响程度。显著扰动代谢通路定义为 FDR &lt; {fdr_cutoff} 且 pathway impact &gt; {impact_cutoff}。")
      )

      x <- tablesAdd(
        x,
        t.metaboAnalyst_pathway = t.metaboAnalyst_pathway
      )
    }

    if (isTRUE(plot) && nrow(data_pathway) > 0L) {
      p.metaboAnalyst_pathway <- metaboAnalystFuns$plot_pathway_impact(
        data_pathway,
        top_n = plot_top_n,
        p_cutoff = p_cutoff,
        fdr_cutoff = fdr_cutoff,
        impact_cutoff = impact_cutoff
      )

      p.metaboAnalyst_pathway <- set_lab_legend(
        p.metaboAnalyst_pathway,
        glue::glue("{x@sig} MetaboAnalyst pathway impact plot"),
        glue::glue("MetaboAnalyst pathway impact 图|||该图展示差异代谢物相关通路的富集显著性和拓扑影响值。横轴为 pathway impact，纵轴为 -log10(P value)，点大小表示命中代谢物数量，颜色表示 FDR。竖向虚线表示 pathway impact = {impact_cutoff}，横向虚线表示 P value = {p_cutoff}。")
      )

      x <- plotsAdd(
        x,
        p.metaboAnalyst_pathway = p.metaboAnalyst_pathway
      )
    }

    x <- snapAdd(x, glue::glue(
      "共输入 {n_input} 个差异代谢物用于 MetaboAnalystR pathway analysis，其中 {n_mapped} 个代谢物完成数据库映射，并获得 {n_result} 条通路分析结果；根据 FDR &lt; {fdr_cutoff} 且 pathway impact &gt; {impact_cutoff} 的标准，共识别到 {n_final} 条显著扰动代谢通路。{ifelse(n_final > 0L, glue::glue('显著通路包括：{text_sig_path}。'), text_sig_path)}"
    ))

    return(x)
  })


# ==========================================================================
# helpers

if (!exists("metaboAnalystFuns", inherits = FALSE)) {
  metaboAnalystFuns <- new.env(parent = emptyenv())
}

metaboAnalystFuns$stop <- function(text)
{
  stop(as.character(glue::glue(text, .envir = parent.frame())), call. = FALSE)
}

metaboAnalystFuns$msg <- function(text, verbose = TRUE)
{
  if (isTRUE(verbose)) {
    message(as.character(glue::glue(text, .envir = parent.frame())))
  }
  invisible(NULL)
}

metaboAnalystFuns$clean_cpds <- function(cpds)
{
  if (is.null(cpds)) {
    return(character())
  }
  if (is.list(cpds)) {
    cpds <- unlist(cpds, use.names = FALSE)
  }
  cpds <- unique(as.character(cpds))
  cpds[!is.na(cpds) & cpds != ""]
}

metaboAnalystFuns$get_map_table <- function(mSet)
{
  if (is.null(mSet$dataSet$map.table)) {
    return(data.frame())
  }

  data_mapped <- tibble::as_tibble(
    data.frame(mSet$dataSet$map.table, check.names = FALSE)
  )

  if ("KEGG" %in% colnames(data_mapped)) {
    data_mapped <- data_mapped[
      !is.na(data_mapped$KEGG) &
        data_mapped$KEGG != "" &
        data_mapped$KEGG != "NA",
      ,
      drop = FALSE
    ]
  }

  data_mapped
}

metaboAnalystFuns$resolve_plot_file <- function(plot_dir,
  plot_base,
  files_before = character())
{
  expected_file <- paste0(plot_base, "dpi72.pdf")
  if (file.exists(expected_file)) {
    return(expected_file)
  }

  files_after <- list.files(
    plot_dir,
    pattern = "\\.pdf$",
    full.names = TRUE
  )
  files_new <- setdiff(files_after, files_before)

  if (length(files_new) == 0L) {
    files_new <- files_after[startsWith(basename(files_after), basename(plot_base))]
  }

  if (length(files_new) == 0L) {
    return(NULL)
  }

  files_new <- files_new[order(file.info(files_new)$mtime, decreasing = TRUE)]
  files_new[1L]
}

metaboAnalystFuns$read_plot_as_binary <- function(file)
{
  if (is.null(file) || !file.exists(file)) {
    return(NULL)
  }

  as_data_binary(.file_fig(file))
}

metaboAnalystFuns$find_col <- function(data_x, candidates)
{
  if (is.null(data_x) || nrow(data_x) == 0L) {
    return(NULL)
  }

  vec_name <- colnames(data_x)
  idx <- match(candidates, vec_name)
  idx <- idx[!is.na(idx)]
  if (length(idx) > 0L) {
    return(vec_name[idx[1L]])
  }

  .clean <- function(x) {
    x <- tolower(as.character(x))
    gsub("[^a-z0-9]+", "", x)
  }

  vec_clean <- .clean(vec_name)
  vec_candidate <- .clean(candidates)
  idx <- match(vec_candidate, vec_clean)
  idx <- idx[!is.na(idx)]
  if (length(idx) > 0L) {
    return(vec_name[idx[1L]])
  }

  return(NULL)
}

metaboAnalystFuns$summarize_ora_table <- function(data_ora,
  p_cutoff = 0.05,
  fdr_cutoff = NULL)
{
  if (is.null(data_ora) || nrow(data_ora) == 0L) {
    return(list(
      p_col = NULL,
      fdr_col = NULL,
      n_p_sig = NA_integer_,
      n_fdr_sig = NA_integer_
    ))
  }

  p_col <- metaboAnalystFuns$find_col(
    data_ora,
    c("Raw p", "Raw.p", "Raw pvalue", "P value", "P.value", "p.value", "pvalue", "P")
  )
  fdr_col <- metaboAnalystFuns$find_col(
    data_ora,
    c("FDR", "False discovery rate", "Adjusted P", "Adjusted.P", "adj.P.Val", "adj.P", "Holm p", "Holm.P")
  )

  if (is.null(p_col)) {
    n_p_sig <- NA_integer_
  } else {
    vec_p <- suppressWarnings(as.numeric(data_ora[[ p_col ]]))
    n_p_sig <- sum(!is.na(vec_p) & vec_p < p_cutoff)
  }

  if (is.null(fdr_col) || is.null(fdr_cutoff)) {
    n_fdr_sig <- NA_integer_
  } else {
    vec_fdr <- suppressWarnings(as.numeric(data_ora[[ fdr_col ]]))
    n_fdr_sig <- sum(!is.na(vec_fdr) & vec_fdr < fdr_cutoff)
  }

  list(
    p_col = p_col,
    fdr_col = fdr_col,
    n_p_sig = n_p_sig,
    n_fdr_sig = n_fdr_sig
  )
}

metaboAnalystFuns$init_data_objects <- function(data_type = "conc",
  anal_type = "msetora", paired = FALSE, default_dpi = 72L)
{
  fun_init <- getExportedValue("MetaboAnalystR", "InitDataObjects")
  vec_formal <- names(formals(fun_init))

  if ("default.dpi" %in% vec_formal) {
    return(e(fun_init(
      data_type,
      anal_type,
      paired,
      default.dpi = default_dpi
    )))
  }

  e(fun_init(data_type, anal_type, paired))
}

metaboAnalystFuns$get_ora_table <- function(mSet)
{
  if (is.null(mSet$analSet$ora.mat)) {
    return(data.frame())
  }

  data_ora <- data.frame(mSet$analSet$ora.mat, check.names = FALSE)
  if (nrow(data_ora) > 0L && !"metabolite_set" %in% colnames(data_ora)) {
    data_ora <- cbind(
      data.frame(metabolite_set = rownames(data_ora), stringsAsFactors = FALSE),
      data_ora,
      stringsAsFactors = FALSE
    )
  }
  rownames(data_ora) <- NULL
  tibble::as_tibble(data_ora)
}


metaboAnalystFuns$is_mset <- function(mSet)
{
  is.list(mSet) && !is.null(mSet$dataSet) && !is.null(mSet$analSet)
}

metaboAnalystFuns$run_ora_safe <- function(cpds,
  cmpd_type = "name",
  libname = "kegg_pathway",
  plot = FALSE,
  plot_prefix = "metabolites_ORA_dot_",
  plot_dir = NULL,
  default_dpi = 72L,
  verbose = TRUE)
{
  res <- tryCatch(
    metaboAnalystFuns$run_ora(
      cpds = cpds,
      cmpd_type = cmpd_type,
      libname = libname,
      plot = plot,
      plot_prefix = plot_prefix,
      plot_dir = plot_dir,
      default_dpi = default_dpi,
      verbose = verbose
    ),
    error = function(e) {
      list(error = conditionMessage(e))
    }
  )

  if (!is.null(res$error)) {
    return(list(
      success = FALSE,
      libname = libname,
      error = res$error,
      result = NULL
    ))
  }

  list(
    success = TRUE,
    libname = libname,
    error = NA_character_,
    result = res
  )
}

metaboAnalystFuns$rbind_fill <- function(lst_data)
{
  lst_data <- lst_data[vapply(lst_data, function(x) {
    is.data.frame(x) && nrow(x) > 0L
  }, logical(1L))]

  if (length(lst_data) == 0L) {
    return(data.frame())
  }

  vec_col <- unique(unlist(lapply(lst_data, colnames), use.names = FALSE))

  lst_data <- lapply(lst_data, function(data_x) {
    vec_missing <- setdiff(vec_col, colnames(data_x))
    if (length(vec_missing) > 0L) {
      for (col_i in vec_missing) {
        data_x[[ col_i ]] <- NA
      }
    }
    data_x[, vec_col, drop = FALSE]
  })

  do.call(rbind, lst_data)
}

metaboAnalystFuns$combine_result_tables <- function(lst_success,
  table_name)
{
  lst_data <- lapply(lst_success, function(z) {
    data_x <- z$result[[ table_name ]]
    if (is.null(data_x) || !is.data.frame(data_x) || nrow(data_x) == 0L) {
      return(data.frame())
    }
    data_x <- as.data.frame(data_x, stringsAsFactors = FALSE, check.names = FALSE)
    data_x$libname <- z$libname
    data_x
  })

  data_out <- metaboAnalystFuns$rbind_fill(lst_data)
  rownames(data_out) <- NULL
  data_out
}

metaboAnalystFuns$get_library_status <- function(lst_res)
{
  data_status <- do.call(rbind, lapply(lst_res, function(z) {
    if (isTRUE(z$success)) {
      data.frame(
        libname = z$libname,
        success = TRUE,
        n_input = length(z$result$cpds),
        n_mapped = nrow(z$result$data_mapped),
        n_result = nrow(z$result$data_ora),
        error = NA_character_,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        libname = z$libname,
        success = FALSE,
        n_input = NA_integer_,
        n_mapped = NA_integer_,
        n_result = NA_integer_,
        error = z$error,
        stringsAsFactors = FALSE
      )
    }
  }))

  rownames(data_status) <- NULL
  data_status
}

metaboAnalystFuns$standardize_enrichment_table <- function(data_ora,
  p_cutoff = 0.05,
  fdr_cutoff = NULL,
  impact_cutoff = NULL)
{
  if (is.null(data_ora) || nrow(data_ora) == 0L) {
    return(data.frame())
  }

  data_out <- as.data.frame(data_ora, stringsAsFactors = FALSE, check.names = FALSE)

  name_col <- metaboAnalystFuns$find_col(
    data_out,
    c("metabolite_set", "Pathway", "Pathway name", "Name", "Mset", "Metabolite Set", "Metabolite.set")
  )
  p_col <- metaboAnalystFuns$find_col(
    data_out,
    c("Raw p", "Raw.p", "Raw pvalue", "P value", "P.value", "p.value", "pvalue", "P")
  )
  fdr_col <- metaboAnalystFuns$find_col(
    data_out,
    c("FDR", "False discovery rate", "Adjusted P", "Adjusted.P", "adj.P.Val", "adj.P", "Holm p", "Holm.P")
  )
  impact_col <- metaboAnalystFuns$find_col(
    data_out,
    c("Impact", "Pathway Impact", "Pathway.Impact", "pathway impact", "impact")
  )
  hit_col <- metaboAnalystFuns$find_col(
    data_out,
    c("Hits", "Hit", "hit.num", "hits", "Observed")
  )

  if (is.null(name_col)) {
    data_out$pathway_name <- rownames(data_out)
  } else {
    data_out$pathway_name <- as.character(data_out[[ name_col ]])
  }

  if (is.null(p_col)) {
    data_out$pvalue <- NA_real_
  } else {
    data_out$pvalue <- suppressWarnings(as.numeric(data_out[[ p_col ]]))
  }

  if (is.null(fdr_col)) {
    data_out$fdr <- NA_real_
  } else {
    data_out$fdr <- suppressWarnings(as.numeric(data_out[[ fdr_col ]]))
  }

  if (is.null(impact_col)) {
    data_out$impact <- NA_real_
  } else {
    data_out$impact <- suppressWarnings(as.numeric(data_out[[ impact_col ]]))
  }

  if (is.null(hit_col)) {
    data_out$hits <- NA_real_
  } else {
    data_out$hits <- suppressWarnings(as.numeric(data_out[[ hit_col ]]))
  }

  data_out$neg_log10_p <- -log10(data_out$pvalue)
  data_out$neg_log10_fdr <- -log10(data_out$fdr)
  data_out$pass_p <- !is.na(data_out$pvalue) & data_out$pvalue < p_cutoff
  data_out$pass_fdr <- if (is.null(fdr_cutoff)) {
    NA
  } else {
    !is.na(data_out$fdr) & data_out$fdr < fdr_cutoff
  }
  data_out$pass_impact <- if (is.null(impact_cutoff)) {
    NA
  } else {
    !is.na(data_out$impact) & data_out$impact > impact_cutoff
  }

  data_out$significant <- data_out$pass_p
  if (!is.null(fdr_cutoff)) {
    data_out$significant <- data_out$significant & data_out$pass_fdr
  }
  if (!is.null(impact_cutoff)) {
    data_out$significant <- data_out$significant & data_out$pass_impact
  }

  data_out
}

metaboAnalystFuns$summarize_enrichment_table <- function(data_enrich,
  p_cutoff = 0.05,
  fdr_cutoff = NULL,
  impact_cutoff = NULL)
{
  if (is.null(data_enrich) || nrow(data_enrich) == 0L) {
    return(list(
      n_p_sig = NA_integer_,
      n_fdr_sig = NA_integer_,
      n_impact_sig = NA_integer_,
      n_final_sig = NA_integer_
    ))
  }

  n_p_sig <- sum(!is.na(data_enrich$pvalue) & data_enrich$pvalue < p_cutoff)
  n_fdr_sig <- if (is.null(fdr_cutoff)) {
    NA_integer_
  } else {
    sum(!is.na(data_enrich$fdr) & data_enrich$fdr < fdr_cutoff)
  }
  n_impact_sig <- if (is.null(impact_cutoff)) {
    NA_integer_
  } else {
    sum(!is.na(data_enrich$impact) & data_enrich$impact > impact_cutoff)
  }
  n_final_sig <- sum(data_enrich$significant, na.rm = TRUE)

  list(
    n_p_sig = n_p_sig,
    n_fdr_sig = n_fdr_sig,
    n_impact_sig = n_impact_sig,
    n_final_sig = n_final_sig
  )
}

metaboAnalystFuns$select_top_enrichment <- function(data_enrich,
  top_n = 10L,
  sort_by = c("pvalue", "fdr", "impact"))
{
  sort_by <- match.arg(sort_by)

  if (is.null(data_enrich) || nrow(data_enrich) == 0L) {
    return(data.frame())
  }

  data_x <- data_enrich[!is.na(data_enrich$pvalue), , drop = FALSE]
  if (nrow(data_x) == 0L) {
    return(data.frame())
  }

  lst_data <- split(data_x, data_x$libname)
  lst_top <- lapply(lst_data, function(data_i) {
    if (sort_by == "fdr" && any(!is.na(data_i$fdr))) {
      data_i <- data_i[order(data_i$fdr, data_i$pvalue), , drop = FALSE]
    } else if (sort_by == "impact" && any(!is.na(data_i$impact))) {
      data_i <- data_i[order(-data_i$impact, data_i$pvalue), , drop = FALSE]
    } else {
      data_i <- data_i[order(data_i$pvalue), , drop = FALSE]
    }
    data_i[seq_len(min(top_n, nrow(data_i))), , drop = FALSE]
  })

  data_top <- metaboAnalystFuns$rbind_fill(lst_top)
  rownames(data_top) <- NULL
  data_top
}

metaboAnalystFuns$plot_enrichment_facet <- function(data_enrich,
  top_n = 10L,
  sort_by = c("pvalue", "fdr", "impact"),
  p_cutoff = 0.05,
  fdr_cutoff = NULL,
  impact_cutoff = NULL)
{
  sort_by <- match.arg(sort_by)

  data_plot <- metaboAnalystFuns$select_top_enrichment(
    data_enrich,
    top_n = top_n,
    sort_by = sort_by
  )

  if (nrow(data_plot) == 0L) {
    metaboAnalystFuns$stop("No enrichment result is available for plotting.")
  }

  data_plot$plot_name <- paste(data_plot$libname, data_plot$pathway_name, sep = "___")
  data_plot <- data_plot[order(data_plot$libname, data_plot$pvalue), , drop = FALSE]
  data_plot$plot_name <- factor(data_plot$plot_name, levels = rev(unique(data_plot$plot_name)))

  if (any(!is.na(data_plot$fdr))) {
    data_plot$color_value <- data_plot$neg_log10_fdr
    color_label <- "-log10(FDR)"
  } else {
    data_plot$color_value <- data_plot$neg_log10_p
    color_label <- "-log10(P value)"
  }

  if (all(is.na(data_plot$hits))) {
    data_plot$hits <- 1
  }

  p <- ggplot2::ggplot(
    data_plot,
    ggplot2::aes(
      x = neg_log10_p,
      y = plot_name,
      size = hits,
      color = color_value
    )
  ) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::facet_wrap(~libname, scales = "free_y") +
    ggplot2::scale_y_discrete(labels = function(z) sub("^.*___", "", z)) +
    ggplot2::geom_vline(
      xintercept = -log10(p_cutoff),
      linetype = 4L,
      linewidth = 0.5
    ) +
    ggplot2::labs(
      x = "-log10(P value)",
      y = NULL,
      size = "Hits",
      color = color_label
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(size = 8)
    )

  p
}

metaboAnalystFuns$run_ora <- function(cpds,
  cmpd_type = "name",
  libname = "kegg_pathway",
  plot = TRUE,
  plot_prefix = "metabolites_ORA_dot_",
  plot_dir = NULL,
  default_dpi = 72L,
  verbose = TRUE)
{
  if (!requireNamespace("MetaboAnalystR", quietly = TRUE)) {
    metaboAnalystFuns$stop(
      "Package `MetaboAnalystR` is required for MetaboAnalyst enrichment analysis."
    )
  }

  if (!is.null(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    plot_dir <- normalizePath(plot_dir, winslash = "/", mustWork = TRUE)
    old_wd <- getwd()
    setwd(plot_dir)
    on.exit(setwd(old_wd), add = TRUE)
  }

  cpds <- metaboAnalystFuns$clean_cpds(cpds)
  if (length(cpds) == 0L) {
    metaboAnalystFuns$stop("No valid compound names were provided.")
  }

  metaboAnalystFuns$msg(
    "Run MetaboAnalystR ORA: {length(cpds)} compound(s), type = {cmpd_type}, lib = {libname}.",
    verbose = verbose
  )

  mSet <- metaboAnalystFuns$init_data_objects(
    data_type = "conc",
    anal_type = "msetora",
    paired = FALSE,
    default_dpi = default_dpi
  )
  if (!metaboAnalystFuns$is_mset(mSet)) {
    metaboAnalystFuns$stop("MetaboAnalystR returned an invalid object after InitDataObjects.")
  }

  mSet <- e(MetaboAnalystR::Setup.MapData(mSet, cpds))
  if (!metaboAnalystFuns$is_mset(mSet)) {
    metaboAnalystFuns$stop("MetaboAnalystR returned an invalid object after Setup.MapData.")
  }

  mSet <- e(MetaboAnalystR::CrossReferencing(mSet, cmpd_type))
  if (!metaboAnalystFuns$is_mset(mSet)) {
    metaboAnalystFuns$stop("MetaboAnalystR returned an invalid object after CrossReferencing.")
  }

  mSet <- e(MetaboAnalystR::CreateMappingResultTable(mSet))
  if (!metaboAnalystFuns$is_mset(mSet)) {
    metaboAnalystFuns$stop("MetaboAnalystR returned an invalid object after CreateMappingResultTable.")
  }

  mSet <- e(MetaboAnalystR::SetMetabolomeFilter(mSet, FALSE))
  if (!metaboAnalystFuns$is_mset(mSet)) {
    metaboAnalystFuns$stop("MetaboAnalystR returned an invalid object after SetMetabolomeFilter.")
  }

  mSet <- e(MetaboAnalystR::SetCurrentMsetLib(mSet, libname, 2L))
  if (!metaboAnalystFuns$is_mset(mSet)) {
    metaboAnalystFuns$stop("MetaboAnalystR returned an invalid object after SetCurrentMsetLib for library `{libname}`.")
  }

  mSet <- e(MetaboAnalystR::CalculateHyperScore(mSet))
  if (!metaboAnalystFuns$is_mset(mSet)) {
    metaboAnalystFuns$stop("MetaboAnalystR returned an invalid object after CalculateHyperScore for library `{libname}`.")
  }

  plot_info <- NULL
  if (isTRUE(plot)) {
    if (is.null(plot_dir)) {
      plot_dir <- file.path(tempdir(), "metaboAnalyst_plot")
    }
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

    plot_base <- file.path(
      plot_dir,
      paste0(plot_prefix, libname, "_")
    )
    files_before <- list.files(
      plot_dir,
      pattern = "\\.pdf$",
      full.names = TRUE
    )

    res_plot <- try(
      e(MetaboAnalystR::PlotEnrichDotPlot(
        mSet,
        "ora",
        plot_base,
        "pdf",
        width = NA
      )),
      silent = TRUE
    )

    if (inherits(res_plot, "try-error")) {
      warning("MetaboAnalystR dot plot generation failed.")
    } else {
      mSet <- res_plot
      plot_file <- metaboAnalystFuns$resolve_plot_file(
        plot_dir = plot_dir,
        plot_base = plot_base,
        files_before = files_before
      )
      plot_fig <- metaboAnalystFuns$read_plot_as_binary(plot_file)
      plot_info <- list(
        dir = plot_dir,
        prefix = plot_base,
        file = plot_file,
        fig = plot_fig
      )
    }
  }

  data_mapped <- metaboAnalystFuns$get_map_table(mSet)
  data_ora <- metaboAnalystFuns$get_ora_table(mSet)

  list(
    mSet = mSet,
    data_mapped = data_mapped,
    data_ora = data_ora,
    hits = mSet$analSet$ora.hits,
    plot_info = plot_info,
    cpds = cpds,
    cmpd_type = cmpd_type,
    libname = libname
  )
}

# --------------------------------------------------------------------------
# Updated helpers for report-ready ORA dot plot and pathway topology analysis

metaboAnalystFuns$parse_numeric <- function(x)
{
  if (is.numeric(x)) {
    return(as.numeric(x))
  }

  vec_x <- trimws(as.character(x))
  vec_x[vec_x %in% c("", "NA", "NaN", "Inf", "-Inf", "NULL", "null")] <- NA_character_
  vec_x <- sub("/.*$", "", vec_x)
  suppressWarnings(as.numeric(vec_x))
}

metaboAnalystFuns$standardize_enrichment_table <- function(data_ora,
  p_cutoff = 0.05,
  fdr_cutoff = NULL,
  impact_cutoff = NULL)
{
  if (is.null(data_ora) || nrow(data_ora) == 0L) {
    return(data.frame())
  }

  data_out <- as.data.frame(data_ora, stringsAsFactors = FALSE, check.names = FALSE)

  name_col <- metaboAnalystFuns$find_col(
    data_out,
    c("metabolite_set", "Pathway", "Pathway name", "Name", "Mset", "Metabolite Set", "Metabolite.set")
  )
  p_col <- metaboAnalystFuns$find_col(
    data_out,
    c("Raw p", "Raw.p", "Raw pvalue", "P value", "P.value", "p.value", "pvalue", "P")
  )
  fdr_col <- metaboAnalystFuns$find_col(
    data_out,
    c("FDR", "False discovery rate", "Adjusted P", "Adjusted.P", "adj.P.Val", "adj.P", "Holm p", "Holm.P")
  )
  impact_col <- metaboAnalystFuns$find_col(
    data_out,
    c("Impact", "Pathway Impact", "Pathway.Impact", "pathway impact", "impact")
  )
  hit_col <- metaboAnalystFuns$find_col(
    data_out,
    c("Hits", "Hit", "hit.num", "hits", "Observed", "Match Status", "Matched")
  )
  expected_col <- metaboAnalystFuns$find_col(
    data_out,
    c("Expected", "expected", "Exp", "exp")
  )
  ratio_col <- metaboAnalystFuns$find_col(
    data_out,
    c("Enrichment Ratio", "Enrichment.Ratio", "EnrichmentRatio", "Rich factor", "Rich.factor", "Ratio", "Fold Enrichment")
  )

  if (is.null(name_col)) {
    data_out$pathway_name <- rownames(data_out)
  } else {
    data_out$pathway_name <- as.character(data_out[[ name_col ]])
  }

  if (is.null(p_col)) {
    data_out$pvalue <- NA_real_
  } else {
    data_out$pvalue <- metaboAnalystFuns$parse_numeric(data_out[[ p_col ]])
  }

  if (is.null(fdr_col)) {
    data_out$fdr <- NA_real_
  } else {
    data_out$fdr <- metaboAnalystFuns$parse_numeric(data_out[[ fdr_col ]])
  }

  if (is.null(impact_col)) {
    data_out$impact <- NA_real_
  } else {
    data_out$impact <- metaboAnalystFuns$parse_numeric(data_out[[ impact_col ]])
  }

  if (is.null(hit_col)) {
    data_out$hits <- NA_real_
  } else {
    data_out$hits <- metaboAnalystFuns$parse_numeric(data_out[[ hit_col ]])
  }

  if (is.null(expected_col)) {
    data_out$expected <- NA_real_
  } else {
    data_out$expected <- metaboAnalystFuns$parse_numeric(data_out[[ expected_col ]])
  }

  if (is.null(ratio_col)) {
    data_out$enrichment_ratio <- NA_real_
  } else {
    data_out$enrichment_ratio <- metaboAnalystFuns$parse_numeric(data_out[[ ratio_col ]])
  }

  idx_ratio_na <- is.na(data_out$enrichment_ratio) |
    !is.finite(data_out$enrichment_ratio)
  idx_can_ratio <- idx_ratio_na &
    !is.na(data_out$hits) &
    !is.na(data_out$expected) &
    data_out$expected > 0
  data_out$enrichment_ratio[idx_can_ratio] <-
    data_out$hits[idx_can_ratio] / data_out$expected[idx_can_ratio]

  data_out$neg_log10_p <- -log10(data_out$pvalue)
  data_out$neg_log10_fdr <- -log10(data_out$fdr)
  data_out$pass_p <- !is.na(data_out$pvalue) & data_out$pvalue < p_cutoff
  data_out$pass_fdr <- if (is.null(fdr_cutoff)) {
    NA
  } else {
    !is.na(data_out$fdr) & data_out$fdr < fdr_cutoff
  }
  data_out$pass_impact <- if (is.null(impact_cutoff)) {
    NA
  } else {
    !is.na(data_out$impact) & data_out$impact > impact_cutoff
  }

  data_out$significant <- data_out$pass_p
  if (!is.null(fdr_cutoff)) {
    data_out$significant <- data_out$significant & data_out$pass_fdr
  }
  if (!is.null(impact_cutoff)) {
    data_out$significant <- data_out$significant & data_out$pass_impact
  }

  data_out
}

metaboAnalystFuns$select_top_enrichment <- function(data_enrich,
  top_n = 10L,
  sort_by = c("pvalue", "fdr", "impact"))
{
  sort_by <- match.arg(sort_by)

  if (is.null(data_enrich) || nrow(data_enrich) == 0L) {
    return(data.frame())
  }

  data_x <- data_enrich[!is.na(data_enrich$pvalue), , drop = FALSE]
  if (nrow(data_x) == 0L) {
    return(data.frame())
  }

  if (!"libname" %in% colnames(data_x)) {
    data_x$libname <- "Metabolite set"
  }

  lst_data <- split(data_x, data_x$libname)
  lst_top <- lapply(lst_data, function(data_i) {
    if (sort_by == "fdr" && any(!is.na(data_i$fdr))) {
      data_i <- data_i[order(data_i$fdr, data_i$pvalue), , drop = FALSE]
    } else if (sort_by == "impact" && any(!is.na(data_i$impact))) {
      data_i <- data_i[order(-data_i$impact, data_i$pvalue), , drop = FALSE]
    } else {
      data_i <- data_i[order(data_i$pvalue), , drop = FALSE]
    }
    data_i[seq_len(min(top_n, nrow(data_i))), , drop = FALSE]
  })

  data_top <- metaboAnalystFuns$rbind_fill(lst_top)
  rownames(data_top) <- NULL
  data_top
}

metaboAnalystFuns$plot_enrichment_facet <- function(data_enrich,
  top_n = 10L,
  sort_by = c("pvalue", "fdr", "impact"),
  p_cutoff = 0.05,
  fdr_cutoff = NULL,
  impact_cutoff = NULL)
{
  sort_by <- match.arg(sort_by)

  data_plot <- metaboAnalystFuns$select_top_enrichment(
    data_enrich,
    top_n = top_n,
    sort_by = sort_by
  )

  if (nrow(data_plot) == 0L) {
    metaboAnalystFuns$stop("No enrichment result is available for plotting.")
  }

  if (!"libname" %in% colnames(data_plot)) {
    data_plot$libname <- "Metabolite set"
  }

  data_plot$plot_name <- paste(data_plot$libname, data_plot$pathway_name, sep = "___")
  data_plot <- data_plot[order(data_plot$libname, data_plot$pvalue), , drop = FALSE]
  data_plot$plot_name <- factor(data_plot$plot_name, levels = rev(unique(data_plot$plot_name)))

  data_plot$size_value <- data_plot$enrichment_ratio
  size_label <- "Enrichment Ratio"
  if (all(is.na(data_plot$size_value))) {
    data_plot$size_value <- data_plot$hits
    size_label <- "Hits"
  }
  if (all(is.na(data_plot$size_value))) {
    data_plot$size_value <- 1
    size_label <- "Value"
  }

  p <- ggplot2::ggplot(
    data_plot,
    ggplot2::aes(
      x = neg_log10_p,
      y = plot_name,
      size = size_value,
      color = pvalue
    )
  ) +
    ggplot2::geom_point(alpha = 0.95) +
    ggplot2::scale_y_discrete(labels = function(z) sub("^.*___", "", z)) +
    ggplot2::geom_vline(
      xintercept = -log10(p_cutoff),
      linetype = 2L,
      linewidth = 0.5,
      color = "grey35"
    ) +
    ggplot2::scale_color_gradient(
      low = "#ff3b1f",
      high = "#f6e8b1",
      name = "P-value"
    ) +
    ggplot2::labs(
      x = "-log10(P value)",
      y = NULL,
      size = size_label,
      title = "Overview of Enriched Metabolite Sets"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(size = 8)
    )

  if (length(unique(data_plot$libname)) > 1L) {
    p <- p + ggplot2::facet_wrap(~libname, scales = "free_y")
  }

  p
}

metaboAnalystFuns$set_kegg_path_lib <- function(mSet,
  organism = "hsa", version = "current")
{
  fun_path <- getExportedValue("MetaboAnalystR", "SetKEGG.PathLib")
  vec_formal <- names(formals(fun_path))

  res <- try(
    {
      if (length(vec_formal) >= 3L || "version" %in% vec_formal) {
        e(fun_path(mSet, organism, version))
      } else {
        e(fun_path(mSet, organism))
      }
    },
    silent = TRUE
  )

  if (inherits(res, "try-error")) {
    metaboAnalystFuns$stop("MetaboAnalystR::SetKEGG.PathLib failed for organism `{organism}`.")
  }

  res
}

metaboAnalystFuns$calculate_ora_score <- function(mSet,
  node_imp = "rbc", enrichment_method = "hyperg")
{
  fun_score <- getExportedValue("MetaboAnalystR", "CalculateOraScore")
  e(fun_score(mSet, node_imp, enrichment_method))
}

metaboAnalystFuns$get_pathway_table <- function(mSet)
{
  data_pathway <- metaboAnalystFuns$get_ora_table(mSet)

  if (nrow(data_pathway) == 0L) {
    return(data_pathway)
  }

  vec_path_names <- try(
    getExportedValue("MetaboAnalystR", "GetORA.pathNames")(mSet),
    silent = TRUE
  )

  if (!inherits(vec_path_names, "try-error") && length(vec_path_names) == nrow(data_pathway)) {
    data_pathway$pathway_name <- as.character(vec_path_names)
  }

  data_pathway
}

metaboAnalystFuns$run_pathway <- function(cpds,
  cmpd_type = "name",
  organism = "hsa",
  node_imp = "rbc",
  enrichment_method = "hyperg",
  default_dpi = 72L,
  work_dir = NULL,
  verbose = TRUE)
{
  if (!is.null(work_dir)) {
    dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
    old_wd <- getwd()
    setwd(work_dir)
    on.exit(setwd(old_wd), add = TRUE)
  }

  if (!requireNamespace("MetaboAnalystR", quietly = TRUE)) {
    metaboAnalystFuns$stop(
      "Package `MetaboAnalystR` is required for MetaboAnalyst pathway analysis."
    )
  }

  cpds <- metaboAnalystFuns$clean_cpds(cpds)
  if (length(cpds) == 0L) {
    metaboAnalystFuns$stop("No valid compound names were provided.")
  }

  metaboAnalystFuns$msg(
    "Run MetaboAnalystR pathway analysis: {length(cpds)} compound(s), type = {cmpd_type}, organism = {organism}.",
    verbose = verbose
  )

  mSet <- metaboAnalystFuns$init_data_objects(
    data_type = "conc",
    anal_type = "pathora",
    paired = FALSE,
    default_dpi = default_dpi
  )
  if (!metaboAnalystFuns$is_mset(mSet)) {
    metaboAnalystFuns$stop("MetaboAnalystR returned an invalid object after InitDataObjects.")
  }

  mSet <- e(MetaboAnalystR::Setup.MapData(mSet, cpds))
  if (!metaboAnalystFuns$is_mset(mSet)) {
    metaboAnalystFuns$stop("MetaboAnalystR returned an invalid object after Setup.MapData.")
  }

  mSet <- e(MetaboAnalystR::CrossReferencing(mSet, cmpd_type))
  if (!metaboAnalystFuns$is_mset(mSet)) {
    metaboAnalystFuns$stop("MetaboAnalystR returned an invalid object after CrossReferencing.")
  }

  mSet <- e(MetaboAnalystR::CreateMappingResultTable(mSet))
  if (!metaboAnalystFuns$is_mset(mSet)) {
    metaboAnalystFuns$stop("MetaboAnalystR returned an invalid object after CreateMappingResultTable.")
  }

  mSet <- metaboAnalystFuns$set_kegg_path_lib(
    mSet,
    organism = organism,
    version = "current"
  )
  if (!metaboAnalystFuns$is_mset(mSet)) {
    metaboAnalystFuns$stop("MetaboAnalystR returned an invalid object after SetKEGG.PathLib.")
  }

  mSet <- e(MetaboAnalystR::SetMetabolomeFilter(mSet, FALSE))
  if (!metaboAnalystFuns$is_mset(mSet)) {
    metaboAnalystFuns$stop("MetaboAnalystR returned an invalid object after SetMetabolomeFilter.")
  }

  mSet <- metaboAnalystFuns$calculate_ora_score(
    mSet,
    node_imp = node_imp,
    enrichment_method = enrichment_method
  )
  if (!metaboAnalystFuns$is_mset(mSet)) {
    metaboAnalystFuns$stop("MetaboAnalystR returned an invalid object after CalculateOraScore.")
  }

  data_mapped <- metaboAnalystFuns$get_map_table(mSet)
  data_pathway <- metaboAnalystFuns$get_pathway_table(mSet)

  list(
    mSet = mSet,
    data_mapped = data_mapped,
    data_pathway = data_pathway,
    hits = mSet$analSet$ora.hits,
    cpds = cpds,
    cmpd_type = cmpd_type,
    organism = organism,
    node_imp = node_imp,
    enrichment_method = enrichment_method
  )
}

metaboAnalystFuns$run_pathway_safe <- function(cpds,
  cmpd_type = "name",
  organism = "hsa",
  node_imp = "rbc",
  enrichment_method = "hyperg",
  default_dpi = 72L,
  work_dir = NULL,
  verbose = TRUE)
{
  res <- tryCatch(
    metaboAnalystFuns$run_pathway(
      cpds = cpds,
      cmpd_type = cmpd_type,
      organism = organism,
      node_imp = node_imp,
      enrichment_method = enrichment_method,
      default_dpi = default_dpi,
      work_dir = work_dir,
      verbose = verbose
    ),
    error = function(e) {
      list(error = conditionMessage(e))
    }
  )

  if (!is.null(res$error)) {
    return(list(
      success = FALSE,
      organism = organism,
      error = res$error,
      result = NULL
    ))
  }

  list(
    success = TRUE,
    organism = organism,
    error = NA_character_,
    result = res
  )
}

metaboAnalystFuns$plot_pathway_impact <- function(data_pathway,
  top_n = 25L,
  p_cutoff = 0.05,
  fdr_cutoff = 0.05,
  impact_cutoff = 0.1)
{
  if (is.null(data_pathway) || nrow(data_pathway) == 0L) {
    metaboAnalystFuns$stop("No pathway analysis result is available for plotting.")
  }

  data_plot <- data_pathway[!is.na(data_pathway$pvalue), , drop = FALSE]
  if (nrow(data_plot) == 0L) {
    metaboAnalystFuns$stop("No valid P value is available for plotting.")
  }

  data_plot <- data_plot[order(data_plot$pvalue, -data_plot$impact), , drop = FALSE]
  data_plot <- data_plot[seq_len(min(top_n, nrow(data_plot))), , drop = FALSE]
  data_plot$plot_fdr <- data_plot$fdr
  if (all(is.na(data_plot$plot_fdr))) {
    data_plot$plot_fdr <- data_plot$pvalue
    color_label <- "P-value"
  } else {
    color_label <- "FDR"
  }
  if (all(is.na(data_plot$hits))) {
    data_plot$hits <- 1
  }
  if (all(is.na(data_plot$impact))) {
    data_plot$impact <- 0
  }

  ggplot2::ggplot(
    data_plot,
    ggplot2::aes(
      x = impact,
      y = neg_log10_p,
      size = hits,
      color = plot_fdr
    )
  ) +
    ggplot2::geom_point(alpha = 0.9) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = pathway_name),
      size = 3,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    ggplot2::geom_vline(
      xintercept = impact_cutoff,
      linetype = 2L,
      linewidth = 0.5,
      color = "grey35"
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(p_cutoff),
      linetype = 2L,
      linewidth = 0.5,
      color = "grey35"
    ) +
    ggplot2::scale_color_gradient(
      low = "#d73027",
      high = "#fee08b",
      name = color_label
    ) +
    ggplot2::labs(
      title = "MetaboAnalyst Pathway Analysis",
      x = "Pathway impact",
      y = "-log10(P value)",
      size = "Hits"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      panel.grid.minor = ggplot2::element_blank()
    )
}
