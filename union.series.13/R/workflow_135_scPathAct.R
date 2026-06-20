# ==========================================================================
# workflow of scPathAct
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# workflow_135_scPathAct_v7.R
# Single-cell pathway activity score workflow.

.job_scPathAct <- setClass("job_scPathAct",
  contains = c("job"),
  prototype = prototype(
    pg = "scPathAct",
    info = c(""),
    cite = "PMID: 23323831",
    method = "",
    tag = "scPathAct",
    analysis = "细胞群通路活性评分"
  )
)

setGeneric("asjob_scPathAct",
  function(x, ...) standardGeneric("asjob_scPathAct"))

setMethod("asjob_scPathAct", signature = c(x = "job_seurat"),
  function(x,
    gene_sets = NULL,
    features = NULL,
    cells = NULL,
    assay = NULL,
    slot = "data",
    layer = NULL,
    sample_col = NULL,
    condition_col = NULL,
    cell_group_col = NULL,
    ...
  )
  {
    seu <- object(x)
    params_source <- params(x)

    assay <- scpaFuns$resolve_assay(seu, assay)
    data_metadata <- scpaFuns$prepare_metadata(seu, params_source)

    vec_cells <- scpaFuns$resolve_cells(
      seu = seu,
      data_metadata = data_metadata,
      cells = cells
    )

    data_metadata <- scpaFuns$align_metadata_to_cells(
      data_metadata = data_metadata,
      cells = vec_cells
    )

    mat_expr <- scpaFuns$get_seurat_expression(
      seu = seu,
      assay = assay,
      slot = slot,
      layer = layer
    )

    vec_cells <- intersect(vec_cells, colnames(mat_expr))
    if (length(vec_cells) == 0L) {
      stop("No selected cells were found in the expression matrix.")
    }

    data_metadata <- scpaFuns$align_metadata_to_cells(
      data_metadata = data_metadata,
      cells = vec_cells
    )

    features_keep <- scpaFuns$resolve_features(
      mat_expr = mat_expr,
      gene_sets = gene_sets,
      features = features
    )

    if (length(features_keep) == 0L) {
      stop("No valid features were retained for scPathAct input.")
    }

    mat_expr <- mat_expr[features_keep, vec_cells, drop = FALSE]

    levels_use <- x$levels
    if (is.null(levels_use)) {
      levels_use <- params_source$levels
    }

    group_by_use <- x$group.by
    if (is.null(group_by_use)) {
      group_by_use <- params_source$group.by
    }

    condition_col <- scpaFuns$resolve_condition_col(
      data_metadata = data_metadata,
      levels = levels_use,
      condition_col = condition_col
    )

    cell_group_col <- scpaFuns$resolve_cell_group_col(
      data_metadata = data_metadata,
      group_by = group_by_use,
      cell_group_col = cell_group_col
    )

    sample_col <- scpaFuns$resolve_sample_col(
      data_metadata = data_metadata,
      sample_col = sample_col
    )

    object_new <- list(
      mat_expr = mat_expr,
      metadata = data_metadata,
      gene_sets = gene_sets,
      assay = assay,
      slot = slot,
      layer = layer,
      cells = vec_cells,
      features = features_keep,
      condition_col = condition_col,
      cell_group_col = cell_group_col,
      sample_col = sample_col,
      levels = levels_use,
      source = list(
        class = class(x),
        pg = x@pg
      )
    )

    y <- .job_scPathAct(object = object_new)
    y@params <- append(y@params, params_source)

    y$metadata <- data_metadata
    y$levels <- levels_use
    y$group.by <- cell_group_col
    y$condition_col <- condition_col
    y$cell_group_col <- cell_group_col
    y$sample_col <- sample_col
    y$assay <- assay
    y$slot <- slot
    y$layer <- layer
    y$gene_sets <- gene_sets
    y$input_n_gene <- nrow(mat_expr)
    y$input_n_cell <- ncol(mat_expr)

    return(y)
  })

setMethod("step0", signature = c(x = "job_scPathAct"),
  function(x)
  {
    step_message("Prepare scRNA-seq pathway activity input with `asjob_scPathAct`.")
    return(x)
  })

setMethod("step1", signature = c(x = "job_scPathAct"),
  function(x,
    db,
    db_anno = NULL,
    mode = c(
      "hallmark gene sets" = "H",
      "curated gene sets" = "C2",
      "positional gene sets" = "C1",
      "regulatory target gene sets" = "C3",
      "computational gene sets" = "C4",
      "ontology gene sets" = "C5",
      "oncogenic signature gene sets" = "C6",
      "all gene sets" = "all"
    ),
    mode_sub = "CP",
    species = "Homo sapiens",
    score_method = c("gsva", "ssgsea"),
    score_unit = c("cell", "sample_cell_group", "sample_condition"),
    aggregation_fun = c("mean"),
    min_gs_size = 10L,
    max_gs_size = 500L,
    kcdf = "auto",
    mx.diff = TRUE,
    abs.ranking = FALSE,
    verbose = TRUE,
    rerun = FALSE,
    ...
  )
  {
    step_message("Pathway activity scoring.")

    obj <- object(x)
    scpaFuns$check_scpathact_object(obj)

    score_method <- match.arg(score_method)
    score_unit <- match.arg(score_unit)
    aggregation_fun <- match.arg(aggregation_fun)
    mode <- match.arg(mode)

    cli::cli_alert_info(glue::glue(
      "Input expression matrix: {nrow(obj$mat_expr)} genes x {ncol(obj$mat_expr)} cells."
    ))
    cli::cli_alert_info(glue::glue(
      "Scoring unit: {scpaFuns$get_score_unit_label(score_unit)}."
    ))

    if (missing(db)) {
      if (!is.null(obj$gene_sets)) {
        db <- scpaFuns$gene_sets_to_db(obj$gene_sets)
        mode_text <- "custom gene set"
      } else {
        mode_sub_use <- mode_sub
        if (!identical(mode, "C2") && identical(mode_sub_use, "CP")) {
          mode_sub_use <- NULL
        }
        cli::cli_alert_info(glue::glue(
          "Preparing MSigDB gene sets: mode = {mode}, sub = {ifelse(is.null(mode_sub_use), 'none', mode_sub_use)}."
        ))
        x <- scpaFuns$set_scpathact_msig_db(
          x = x,
          mode = mode,
          sub = mode_sub_use,
          species = species
        )
        db <- x$msig_db
        mode_text <- scpaFuns$get_mode_text(mode = mode, mode_sub = mode_sub_use)
      }
    } else {
      mode_text <- "custom gene set"
    }

    if (is.null(db_anno)) {
      db_anno <- x$db_anno
    }

    cli::cli_alert_info("Matching gene sets to expression features.")
    lst_gs <- scpaFuns$prepare_gene_sets(
      db = db,
      features = rownames(obj$mat_expr),
      min_gs_size = min_gs_size,
      max_gs_size = max_gs_size
    )

    if (length(lst_gs$gene_sets) == 0L) {
      stop("No gene set passed size filtering after matching expression genes.")
    }

    cli::cli_alert_info(glue::glue(
      "Matched gene sets: {length(lst_gs$gene_sets)} retained from {nrow(lst_gs$gene_set_stat)} total records."
    ))

    data_score_input <- scpaFuns$prepare_pathway_score_input(
      obj = obj,
      score_unit = score_unit,
      aggregation_fun = aggregation_fun
    )

    mat_expr_score <- data_score_input$mat_expr
    data_unit_metadata <- data_score_input$metadata

    cli::cli_alert_info(glue::glue(
      "Running {score_method} on {nrow(mat_expr_score)} genes x {ncol(mat_expr_score)} {scpaFuns$get_score_unit_label(score_unit)} units and {length(lst_gs$gene_sets)} gene sets."
    ))

    args_gsva <- list(
      mat_expr = mat_expr_score,
      gene_sets = lst_gs$gene_sets,
      score_method = score_method,
      kcdf = kcdf,
      mx.diff = mx.diff,
      abs.ranking = abs.ranking,
      verbose = verbose,
      ...
    )

    mat_score <- expect_local_data(
      "tmp", "scPathAct_score", scpaFuns$run_pathway_score,
      args_gsva, rerun = rerun
    )

    data_score_summary <- scpaFuns$summarize_score_matrix(mat_score)
    data_unit_stat <- scpaFuns$summarize_score_units(data_unit_metadata)

    t.gene_set_stat <- set_lab_legend(
      lst_gs$gene_set_stat,
      glue::glue("scPathAct gene set matching summary"),
      glue::glue(
        "通路活性评分基因集匹配表|||该表展示 {mode_text} 基因集中每个通路在单细胞表达矩阵中的基因匹配数量，",
        "并根据 {min_gs_size}–{max_gs_size} 个基因的大小阈值保留用于通路活性评分的基因集。"
      )
    )

    t.score_unit_stat <- set_lab_legend(
      data_unit_stat,
      glue::glue("scPathAct pathway score unit summary"),
      glue::glue(
        "通路活性评分单位汇总表|||该表展示本分析用于计算 pathway-level score 的表达单位数量、",
        "分组字段及每个表达单位包含的细胞数量。"
      )
    )

    t.score_summary <- set_lab_legend(
      data_score_summary,
      glue::glue("scPathAct pathway activity score summary"),
      glue::glue(
        "通路活性评分汇总表|||该表展示每个基因集在 {scpaFuns$get_score_unit_label(score_unit)} 表达单位中的 pathway-level score 分布，",
        "用于概览目标细胞群的通路活性差异。"
      )
    )

    x$msig_db <- db
    x$db_anno <- db_anno
    x$gene_sets_used <- lst_gs$gene_sets
    x$gene_set_stat <- lst_gs$gene_set_stat
    x$pathway_score_input <- mat_expr_score
    x$pathway_score_metadata <- data_unit_metadata
    x$pathway_score <- mat_score
    x$pathway_score_summary <- data_score_summary
    x$score_unit <- score_unit
    x$score_method <- score_method

    if (is.null(x$lst_refine)) {
      x$lst_refine <- list()
    }
    x$lst_refine$msig_db <- db
    x$lst_refine$db_anno <- db_anno
    x$lst_refine$gene_sets_used <- lst_gs$gene_sets
    x$lst_refine$gene_set_stat <- lst_gs$gene_set_stat
    x$lst_refine$pathway_score_input <- mat_expr_score
    x$lst_refine$pathway_score_metadata <- data_unit_metadata
    x$lst_refine$pathway_score <- mat_score
    x$lst_refine$pathway_score_summary <- data_score_summary
    x$lst_refine$pathway_score_unit_stat <- data_unit_stat

    x <- tablesAdd(
      x,
      t.gene_set_stat = t.gene_set_stat,
      t.score_unit_stat = t.score_unit_stat,
      t.score_summary = t.score_summary
    )

    x <- methodAdd(
      x,
      scpaFuns$format_scoring_method_text(
        x = x,
        obj = obj,
        score_method = score_method,
        score_unit = score_unit,
        mode_text = mode_text,
        min_gs_size = min_gs_size,
        max_gs_size = max_gs_size,
        n_gene = nrow(mat_expr_score),
        n_unit = ncol(mat_expr_score),
        n_gene_set = nrow(mat_score),
        kcdf = kcdf
      )
    )

    x <- snapAdd(
      x,
      scpaFuns$format_score_snap(
        mat_score = mat_score,
        score_unit = score_unit,
        data_unit_metadata = data_unit_metadata
      )
    )

    return(x)
  })

setMethod("step2", signature = c(x = "job_scPathAct"),
  function(x,
    compare_by = c("cell_group", "condition"),
    use_p = c("adj.P.Val", "P.Value"),
    cutoff = 0.05,
    cutoff.delta_score = 0,
    top_n = 20L,
    rank_by = c("abs_t", "t", "abs_delta_score", "delta_score"),
    plot_select = c("balanced", "overall"),
    add_plot = TRUE
  )
  {
    step_message("Differential pathway activity analysis.")

    compare_by <- match.arg(compare_by)
    use_p <- match.arg(use_p)
    rank_by <- match.arg(rank_by)
    plot_select <- match.arg(plot_select)

    obj <- object(x)
    mat_score <- x$pathway_score
    if (is.null(mat_score)) {
      stop("Run step1 before differential pathway activity analysis.")
    }

    data_metadata <- x$pathway_score_metadata
    if (is.null(data_metadata)) {
      data_metadata <- obj$metadata
      data_metadata$unit <- data_metadata$cell
    }
    data_metadata <- scpaFuns$align_metadata_to_units(
      data_metadata = data_metadata,
      units = colnames(mat_score)
    )

    compare_col <- if (identical(compare_by, "cell_group")) {
      obj$cell_group_col
    } else {
      obj$condition_col
    }

    if (is.null(compare_col) || !compare_col %in% colnames(data_metadata)) {
      stop("No valid comparison column was found in pathway score metadata.")
    }

    cli::cli_alert_info(glue::glue(
      "Comparing pathway scores by `{compare_col}` across {ncol(mat_score)} score units."
    ))

    data_diff <- scpaFuns$run_limma_diff(
      mat_score = mat_score,
      data_metadata = data_metadata,
      compare_col = compare_col,
      compare_by = compare_by,
      levels = obj$levels
    )

    data_diff <- scpaFuns$annotate_pathway_terms(
      data_diff = data_diff,
      db_anno = x$db_anno
    )
    data_diff <- scpaFuns$add_diff_threshold_status(
      data_diff = data_diff,
      use_p = use_p,
      cutoff = cutoff,
      cutoff.delta_score = cutoff.delta_score
    )

    data_report <- scpaFuns$filter_diff_report(
      data_diff = data_diff,
      use_p = use_p,
      cutoff = cutoff,
      cutoff.delta_score = cutoff.delta_score,
      top_n = top_n,
      rank_by = rank_by,
      select_mode = plot_select
    )

    data_plot <- scpaFuns$select_diff_plot_data(
      data_diff = data_diff,
      use_p = use_p,
      cutoff = cutoff,
      cutoff.delta_score = cutoff.delta_score,
      top_n = top_n,
      rank_by = rank_by,
      select_mode = plot_select
    )

    t.pathway_diff_full <- set_lab_legend(
      data_diff,
      glue::glue("scPathAct differential pathway activity full table"),
      glue::glue(
        "通路活性差异分析完整表|||该表展示不同 `{compare_col}` 分组间 pathway-level score 的差异分析结果；",
        "delta_score 表示目标组相对于参照组的 pathway-level score 差异，",
        "该值对应 limma 输出中的 logFC 列，但在本分析中表示通路评分差异。",
        "结合 {scpaFuns$get_threshold_text(use_p, cutoff, cutoff.delta_score)} 标记差异通路。"
      )
    )

    t.pathway_diff <- set_lab_legend(
      data_report,
      glue::glue("scPathAct differential pathway activity table"),
      glue::glue(
        "通路活性差异分析结果表|||该表展示满足 {scpaFuns$get_threshold_text(use_p, cutoff, cutoff.delta_score)} 的差异通路，",
        "并按 {scpaFuns$get_plot_select_text(top_n, rank_by, plot_select)} 得到代表性结果；",
        "delta_score 表示目标组相对于参照组的 pathway-level score 差异，对应 limma 输出中的 logFC 列。"
      )
    )

    x$pathway_diff_full <- data_diff
    x$pathway_diff <- data_report
    x$pathway_diff_plot <- data_plot
    if (is.null(x$lst_refine)) {
      x$lst_refine <- list()
    }
    x$lst_refine$pathway_diff_full <- data_diff
    x$lst_refine$pathway_diff_report <- data_report
    x$lst_refine$pathway_diff_plot <- data_plot

    x <- tablesAdd(
      x,
      t.pathway_diff = t.pathway_diff,
      t.pathway_diff_full = t.pathway_diff_full
    )

    if (isTRUE(add_plot) && nrow(data_plot) > 0L) {
      p.pathway_t_heatmap <- scpaFuns$plot_pathway_t_heatmap(
        data_plot = data_plot,
        compare_by = compare_by
      )
      p.pathway_t_heatmap <- set_lab_legend(
        p.pathway_t_heatmap,
        glue::glue("scPathAct differential pathway activity t-value plot"),
        glue::glue(
          "GSVA 通路活性差异图|||该图展示每个比较中满足 {scpaFuns$get_threshold_text(use_p, cutoff, cutoff.delta_score)} 后，",
          "按 {scpaFuns$get_plot_select_text(top_n, rank_by, plot_select)} 筛选的代表性通路；横轴为 limma t 值，",
          "红色表示目标组评分升高，绿色表示目标组评分降低，灰色表示未达到设定阈值。"
        )
      )
      x <- plotsAdd(x, p.pathway_t_heatmap = p.pathway_t_heatmap)
    }

    x <- methodAdd(
      x,
      scpaFuns$format_diff_method_text(
        compare_by = compare_by,
        compare_col = compare_col,
        use_p = use_p,
        cutoff = cutoff,
        cutoff.delta_score = cutoff.delta_score,
        top_n = top_n,
        rank_by = rank_by,
        select_mode = plot_select
      )
    )

    if (nrow(data_report) > 0L) {
      x <- snapAdd(
        x,
        scpaFuns$format_diff_snap(
          data_diff = data_diff,
          data_report = data_report,
          compare_col = compare_col,
          use_p = use_p,
          cutoff = cutoff,
          cutoff.delta_score = cutoff.delta_score
        )
      )
    }

    return(x)
  })

if (!exists("scpaFuns")) {
  scpaFuns <- new.env(parent = emptyenv())
}

scpaFuns$set_scpathact_msig_db <- function(x, mode, sub = NULL, species = "Homo sapiens")
{
  fun_data <- function(mode)
  {
    if (packageVersion("msigdbr") < "10.0.0") {
      db_anno <- e(msigdbr::msigdbr(species = species, category = mode))
    } else {
      db_anno <- e(msigdbr::msigdbr(species = species, collection = mode))
    }
    return(db_anno)
  }

  if (length(mode) == 1L && mode != "all") {
    db_anno <- fun_data(mode)
    x <- methodAdd(
      x,
      glue::glue("以 R 包 `msigdbr` ⟦pkgInfo('msigdbr')⟧ 获取 MSigDB 数据库 {mode} 基因集。")
    )
    if (!is.null(sub) && "gs_subcat" %in% colnames(db_anno)) {
      if (identical(sub, "CP")) {
        select <- c("CP:REACTOME", "CP:KEGG", "CP:WIKIPATHWAYS")
      } else {
        select <- as.character(sub)
      }
      db_anno <- dplyr::filter(db_anno, gs_subcat %in% !!select)
      x <- methodAdd(x, glue::glue("该基因集包含多个子集：{try_snap(db_anno, 'gs_subcat', 'gs_name')}。"))
      x <- methodAdd(x, glue::glue("选取 {bind(select)} 子集用于后续分析。"))
    }
  } else {
    if (length(mode) == 1L && mode == "all") {
      mode <- c("H", paste0("C", 1L:8L))
    }
    db_anno <- lapply(mode,
      function(type) {
        fun_data(type)
      })
    db_anno <- dplyr::bind_rows(db_anno)
    x <- methodAdd(
      x,
      glue::glue("以 R 包 `msigdbr` ⟦pkgInfo('msigdbr')⟧ 获取 MSigDB 数据库 {bind(mode)} 基因集。")
    )
  }

  x$mode <- mode
  x$db_anno <- db_anno
  x$msig_db <- dplyr::select(db_anno, gs_id, symbol = gene_symbol)
  return(x)
}

scpaFuns$resolve_assay <- function(seu, assay = NULL)
{
  assay_names <- names(seu@assays)
  if (length(assay_names) == 0L) {
    stop("No assay was found in the Seurat object.")
  }

  if (is.null(assay)) {
    assay <- SeuratObject::DefaultAssay(seu)
  }

  if (!assay %in% assay_names) {
    stop(paste0("Assay `", assay, "` was not found in the Seurat object."))
  }

  return(assay)
}

scpaFuns$prepare_metadata <- function(seu, params = list())
{
  data_metadata <- NULL

  if (!is.null(params$metadata)) {
    data_metadata <- as.data.frame(params$metadata, stringsAsFactors = FALSE)
  }

  if (is.null(data_metadata) || nrow(data_metadata) == 0L) {
    if ("meta.data" %in% slotNames(seu)) {
      data_metadata <- as.data.frame(seu@meta.data, stringsAsFactors = FALSE)
    }
  }

  if (is.null(data_metadata) || nrow(data_metadata) == 0L) {
    data_metadata <- data.frame(
      cell = colnames(seu),
      stringsAsFactors = FALSE
    )
  }

  if (!"cell" %in% colnames(data_metadata)) {
    rn <- rownames(data_metadata)
    if (!is.null(rn) && length(rn) == nrow(data_metadata)) {
      data_metadata$cell <- rn
    }
  }

  if (!"cell" %in% colnames(data_metadata)) {
    stop("Metadata must contain a `cell` column or valid row names.")
  }

  data_metadata$cell <- as.character(data_metadata$cell)
  data_metadata <- data_metadata[!is.na(data_metadata$cell), , drop = FALSE]

  if (anyDuplicated(data_metadata$cell) > 0L) {
    stop("Duplicated cell IDs were found in metadata.")
  }

  rownames(data_metadata) <- data_metadata$cell
  return(data_metadata)
}

scpaFuns$resolve_cells <- function(seu, data_metadata, cells = NULL)
{
  cells_object <- colnames(seu)

  if (!is.null(cells)) {
    cells <- as.character(cells)
    cells <- intersect(cells, cells_object)
    cells <- intersect(cells, data_metadata$cell)
    return(cells)
  }

  cells <- intersect(data_metadata$cell, cells_object)
  if (length(cells) == 0L) {
    cells <- cells_object
  }

  return(cells)
}

scpaFuns$align_metadata_to_cells <- function(data_metadata, cells)
{
  cells <- as.character(cells)
  idx <- match(cells, data_metadata$cell)

  if (any(is.na(idx))) {
    stop("Some selected cells were not found in metadata.")
  }

  data_metadata <- data_metadata[idx, , drop = FALSE]
  rownames(data_metadata) <- data_metadata$cell

  return(data_metadata)
}

scpaFuns$get_seurat_expression <- function(seu,
  assay,
  slot = "data",
  layer = NULL
)
{
  if (!is.null(layer)) {
    mat_expr <- SeuratObject::GetAssayData(
      object = seu,
      assay = assay,
      layer = layer
    )
  } else {
    mat_expr <- tryCatch(
      SeuratObject::GetAssayData(
        object = seu,
        assay = assay,
        slot = slot
      ),
      error = function(e) {
        SeuratObject::GetAssayData(
          object = seu,
          assay = assay,
          layer = slot
        )
      }
    )
  }

  if (is.null(mat_expr) || nrow(mat_expr) == 0L || ncol(mat_expr) == 0L) {
    stop("The extracted expression matrix is empty.")
  }

  return(mat_expr)
}

scpaFuns$resolve_features <- function(mat_expr,
  gene_sets = NULL,
  features = NULL
)
{
  features_available <- rownames(mat_expr)
  if (is.null(features_available)) {
    stop("The expression matrix must have gene names as row names.")
  }

  features_keep <- features_available

  if (!is.null(features)) {
    features <- as.character(features)
    features_keep <- intersect(features_keep, features)
  }

  features_gene_set <- scpaFuns$get_gene_set_features(gene_sets)
  if (length(features_gene_set) > 0L) {
    features_keep <- intersect(features_keep, features_gene_set)
  }

  return(features_keep)
}

scpaFuns$get_gene_set_features <- function(gene_sets)
{
  if (is.null(gene_sets)) {
    return(character(0L))
  }

  if (is.list(gene_sets) && !is.data.frame(gene_sets)) {
    features <- unique(as.character(unlist(gene_sets, use.names = FALSE)))
    features <- features[!is.na(features) & nzchar(features)]
    return(features)
  }

  if (is.data.frame(gene_sets)) {
    col_candidates <- c("gene", "Gene", "genes", "Genes", "symbol", "Symbol", "gene_symbol")
    col_gene <- intersect(col_candidates, colnames(gene_sets))
    if (length(col_gene) > 0L) {
      features <- unique(as.character(gene_sets[[col_gene[1L]]]))
      features <- features[!is.na(features) & nzchar(features)]
      return(features)
    }
  }

  features <- tryCatch(
    unique(as.character(unlist(gene_sets, use.names = FALSE))),
    error = function(e) character(0L)
  )
  features <- features[!is.na(features) & nzchar(features)]

  return(features)
}

scpaFuns$resolve_condition_col <- function(data_metadata,
  levels = NULL,
  condition_col = NULL
)
{
  if (!is.null(condition_col)) {
    if (!condition_col %in% colnames(data_metadata)) {
      stop(paste0("condition_col `", condition_col, "` was not found in metadata."))
    }
    return(condition_col)
  }

  if (is.null(levels) || length(levels) == 0L) {
    return(NULL)
  }

  levels <- as.character(levels)
  col_candidates <- scpaFuns$get_cols_containing_levels(data_metadata, levels)

  if (length(col_candidates) == 0L) {
    return(NULL)
  }

  col_prefer <- c(
    "group", "Group", "condition", "Condition", "disease", "Disease",
    "status", "Status", "phenotype", "Phenotype"
  )
  col_hit <- intersect(col_prefer, col_candidates)

  if (length(col_hit) > 0L) {
    return(col_hit[1L])
  }

  return(col_candidates[1L])
}

scpaFuns$get_cols_containing_levels <- function(data_metadata, levels)
{
  col_candidates <- character(0L)

  for (nm in colnames(data_metadata)) {
    values <- unique(as.character(data_metadata[[nm]]))
    values <- values[!is.na(values)]
    if (all(levels %in% values)) {
      col_candidates <- c(col_candidates, nm)
    }
  }

  return(col_candidates)
}

scpaFuns$resolve_cell_group_col <- function(data_metadata,
  group_by = NULL,
  cell_group_col = NULL
)
{
  if (!is.null(cell_group_col)) {
    if (!cell_group_col %in% colnames(data_metadata)) {
      stop(paste0("cell_group_col `", cell_group_col, "` was not found in metadata."))
    }
    return(cell_group_col)
  }

  if (!is.null(group_by) && length(group_by) > 0L) {
    group_by <- as.character(group_by[1L])
    if (group_by %in% colnames(data_metadata)) {
      return(group_by)
    }
  }

  col_prefer <- c(
    "celltype", "cell_type", "scsa_cell", "annotation", "seurat_clusters",
    "cluster", "clusters"
  )
  col_hit <- intersect(col_prefer, colnames(data_metadata))

  if (length(col_hit) > 0L) {
    return(col_hit[1L])
  }

  return(NULL)
}

scpaFuns$resolve_sample_col <- function(data_metadata, sample_col = NULL)
{
  if (!is.null(sample_col)) {
    if (!sample_col %in% colnames(data_metadata)) {
      stop(paste0("sample_col `", sample_col, "` was not found in metadata."))
    }
    return(sample_col)
  }

  col_prefer <- c(
    "orig.ident", "sample", "Sample", "sample_id", "SampleID",
    "donor", "Donor", "patient", "Patient", "subject", "Subject"
  )
  col_hit <- intersect(col_prefer, colnames(data_metadata))

  if (length(col_hit) == 0L) {
    return(NULL)
  }

  for (nm in col_hit) {
    values <- unique(as.character(data_metadata[[nm]]))
    values <- values[!is.na(values)]
    if (length(values) > 1L) {
      return(nm)
    }
  }

  return(col_hit[1L])
}

scpaFuns$check_scpathact_object <- function(obj)
{
  required_names <- c("mat_expr", "metadata")
  missing_names <- setdiff(required_names, names(obj))

  if (length(missing_names) > 0L) {
    stop(paste0(
      "Invalid scPathAct object. Missing fields: ",
      paste(missing_names, collapse = ", ")
    ))
  }

  if (nrow(obj$mat_expr) == 0L || ncol(obj$mat_expr) == 0L) {
    stop("The scPathAct expression matrix is empty.")
  }

  if (!"cell" %in% colnames(obj$metadata)) {
    stop("The scPathAct metadata must contain a `cell` column.")
  }

  if (!identical(as.character(obj$metadata$cell), colnames(obj$mat_expr))) {
    stop("Metadata cell order must be identical to expression matrix column order.")
  }

  return(TRUE)
}

scpaFuns$gene_sets_to_db <- function(gene_sets)
{
  if (is.data.frame(gene_sets)) {
    nm <- colnames(gene_sets)
    term_col <- intersect(c("gs_id", "term", "pathway", "ID", "id", "gs_name"), nm)
    gene_col <- intersect(c("symbol", "gene", "Gene", "gene_symbol"), nm)
    if (length(term_col) == 0L || length(gene_col) == 0L) {
      stop("Custom gene set data.frame must contain term and gene columns.")
    }
    db <- data.frame(
      gs_id = as.character(gene_sets[[term_col[1L]]]),
      symbol = as.character(gene_sets[[gene_col[1L]]]),
      stringsAsFactors = FALSE
    )
    return(db)
  }

  if (!is.list(gene_sets)) {
    stop("Custom gene_sets must be a named list or a TERM2GENE-like data.frame.")
  }

  if (is.null(names(gene_sets)) || any(!nzchar(names(gene_sets)))) {
    stop("Custom gene set list must be named.")
  }

  db <- do.call(rbind, lapply(names(gene_sets),
    function(nm) {
      data.frame(
        gs_id = nm,
        symbol = as.character(gene_sets[[nm]]),
        stringsAsFactors = FALSE
      )
    }))
  rownames(db) <- NULL

  return(db)
}

scpaFuns$prepare_gene_sets <- function(db,
  features,
  min_gs_size = 10L,
  max_gs_size = 500L
)
{
  db <- as.data.frame(db, stringsAsFactors = FALSE)
  nm <- colnames(db)
  term_col <- intersect(c("gs_id", "term", "pathway", "ID", "id", "gs_name"), nm)
  gene_col <- intersect(c("symbol", "gene", "Gene", "gene_symbol"), nm)

  if (length(term_col) == 0L || length(gene_col) == 0L) {
    stop("Gene set database must contain term and gene columns.")
  }

  term_col <- term_col[1L]
  gene_col <- gene_col[1L]
  db[[term_col]] <- as.character(db[[term_col]])
  db[[gene_col]] <- as.character(db[[gene_col]])
  db <- db[!is.na(db[[term_col]]) & !is.na(db[[gene_col]]), , drop = FALSE]
  db <- db[nzchar(db[[term_col]]) & nzchar(db[[gene_col]]), , drop = FALSE]

  lst_raw <- split(db[[gene_col]], db[[term_col]])
  lst_raw <- lapply(lst_raw, unique)

  gene_sets <- lapply(lst_raw,
    function(x) {
      intersect(x, features)
    })
  n_total <- vapply(lst_raw, length, integer(1L))
  n_match <- vapply(gene_sets, length, integer(1L))
  keep <- n_match >= min_gs_size & n_match <= max_gs_size
  gene_sets_used <- gene_sets[keep]

  data_stat <- data.frame(
    ID = names(gene_sets),
    n_gene_database = as.integer(n_total),
    n_gene_matched = as.integer(n_match),
    keep = as.logical(keep),
    stringsAsFactors = FALSE
  )

  data_stat <- data_stat[order(data_stat$keep, data_stat$n_gene_matched, decreasing = TRUE), , drop = FALSE]
  rownames(data_stat) <- NULL

  return(list(
    gene_sets = gene_sets_used,
    gene_set_stat = data_stat
  ))
}

scpaFuns$get_score_unit_label <- function(score_unit)
{
  if (identical(score_unit, "sample_cell_group")) {
    return("sample-cell-group")
  }
  if (identical(score_unit, "sample_condition")) {
    return("sample-condition")
  }
  return("single-cell")
}

scpaFuns$get_mode_text <- function(mode, mode_sub = NULL)
{
  if (is.null(mode_sub)) {
    return(as.character(mode))
  }
  return(paste0(as.character(mode), "-", as.character(mode_sub)))
}

scpaFuns$prepare_pathway_score_input <- function(obj,
  score_unit = "cell",
  aggregation_fun = "mean"
)
{
  mat_expr <- obj$mat_expr
  data_metadata <- obj$metadata

  if (identical(score_unit, "cell")) {
    data_unit <- data_metadata
    data_unit$unit <- data_unit$cell
    data_unit$n_cell <- 1L
    data_unit$score_unit <- "cell"
    return(list(
      mat_expr = mat_expr,
      metadata = data_unit
    ))
  }

  if (!identical(aggregation_fun, "mean")) {
    stop("Only mean aggregation is currently supported.")
  }

  if (identical(score_unit, "sample_cell_group")) {
    if (is.null(obj$sample_col) || !obj$sample_col %in% colnames(data_metadata)) {
      stop("score_unit = `sample_cell_group` requires a valid sample_col.")
    }
    if (is.null(obj$cell_group_col) || !obj$cell_group_col %in% colnames(data_metadata)) {
      stop("score_unit = `sample_cell_group` requires a valid cell_group_col.")
    }
    group_cols <- unique(c(obj$sample_col, obj$condition_col, obj$cell_group_col))
  } else if (identical(score_unit, "sample_condition")) {
    if (is.null(obj$sample_col) || !obj$sample_col %in% colnames(data_metadata)) {
      stop("score_unit = `sample_condition` requires a valid sample_col.")
    }
    if (is.null(obj$condition_col) || !obj$condition_col %in% colnames(data_metadata)) {
      stop("score_unit = `sample_condition` requires a valid condition_col.")
    }
    group_cols <- unique(c(obj$sample_col, obj$condition_col))
  } else {
    stop("Unsupported score_unit.")
  }

  group_cols <- group_cols[!is.na(group_cols) & nzchar(group_cols)]
  group_cols <- group_cols[group_cols %in% colnames(data_metadata)]

  data_key <- data_metadata[, group_cols, drop = FALSE]
  data_key[] <- lapply(data_key, function(x) as.character(x))
  unit <- do.call(paste, c(data_key, sep = "__"))
  unit[is.na(unit) | !nzchar(unit)] <- "Unknown"

  mat_unit <- scpaFuns$aggregate_expression_by_unit(
    mat_expr = mat_expr,
    unit = unit
  )
  data_unit <- scpaFuns$summarize_unit_metadata(
    data_metadata = data_metadata,
    unit = unit,
    group_cols = group_cols,
    score_unit = score_unit
  )

  data_unit <- scpaFuns$align_metadata_to_units(
    data_metadata = data_unit,
    units = colnames(mat_unit)
  )

  return(list(
    mat_expr = mat_unit,
    metadata = data_unit
  ))
}

scpaFuns$aggregate_expression_by_unit <- function(mat_expr, unit)
{
  unit <- as.character(unit)
  if (length(unit) != ncol(mat_expr)) {
    stop("The length of score-unit labels must match the number of expression columns.")
  }

  unit <- factor(unit, levels = unique(unit))
  lst_idx <- split(seq_along(unit), unit)

  mat_unit <- do.call(cbind, lapply(lst_idx,
    function(idx) {
      scpaFuns$row_mean_safe(mat_expr[, idx, drop = FALSE])
    }))

  rownames(mat_unit) <- rownames(mat_expr)
  colnames(mat_unit) <- names(lst_idx)
  mat_unit <- as.matrix(mat_unit)

  return(mat_unit)
}

scpaFuns$row_mean_safe <- function(mat)
{
  if (inherits(mat, "Matrix")) {
    return(Matrix::rowMeans(mat))
  }
  return(rowMeans(mat))
}

scpaFuns$summarize_unit_metadata <- function(data_metadata,
  unit,
  group_cols,
  score_unit
)
{
  unit <- as.character(unit)
  units <- unique(unit)

  data_unit <- do.call(rbind, lapply(units,
    function(one_unit) {
      idx <- which(unit == one_unit)
      row <- data_metadata[idx[1L], group_cols, drop = FALSE]
      row[] <- lapply(row, as.character)
      row$unit <- one_unit
      row$n_cell <- length(idx)
      row$score_unit <- score_unit
      return(row)
    }))

  rownames(data_unit) <- NULL
  data_unit <- as.data.frame(data_unit, stringsAsFactors = FALSE)
  data_unit$n_cell <- as.integer(data_unit$n_cell)

  return(data_unit)
}

scpaFuns$align_metadata_to_units <- function(data_metadata, units)
{
  units <- as.character(units)
  if (!"unit" %in% colnames(data_metadata)) {
    if ("cell" %in% colnames(data_metadata)) {
      data_metadata$unit <- data_metadata$cell
    } else {
      stop("Pathway score metadata must contain a `unit` column.")
    }
  }

  data_metadata$unit <- as.character(data_metadata$unit)
  idx <- match(units, data_metadata$unit)

  if (any(is.na(idx))) {
    stop("Some pathway score units were not found in metadata.")
  }

  data_metadata <- data_metadata[idx, , drop = FALSE]
  rownames(data_metadata) <- data_metadata$unit

  return(data_metadata)
}

scpaFuns$summarize_score_units <- function(data_unit_metadata)
{
  data_unit_metadata <- as.data.frame(data_unit_metadata, stringsAsFactors = FALSE)
  cols_keep <- intersect(
    c("unit", "score_unit", "n_cell", "orig.ident", "sample", "Sample", "sample_id", "group", "condition", "scsa_cell", "celltype", "cell_type", "seurat_clusters"),
    colnames(data_unit_metadata)
  )

  data_unit <- data_unit_metadata[, cols_keep, drop = FALSE]
  if (!"unit" %in% colnames(data_unit)) {
    data_unit$unit <- data_unit_metadata$unit
  }
  if (!"n_cell" %in% colnames(data_unit)) {
    data_unit$n_cell <- 1L
  }

  rownames(data_unit) <- NULL
  return(data_unit)
}

scpaFuns$format_scoring_method_text <- function(x,
  obj,
  score_method,
  score_unit,
  mode_text,
  min_gs_size,
  max_gs_size,
  n_gene,
  n_unit,
  n_gene_set,
  kcdf
)
{
  assay_text <- obj$assay
  layer_text <- if (!is.null(obj$layer)) {
    paste0("layer = ", obj$layer)
  } else {
    paste0("slot = ", obj$slot)
  }

  unit_text <- if (identical(score_unit, "cell")) {
    "本分析以单个细胞作为通路评分单位，得到每个细胞的 pathway-level score。"
  } else if (identical(score_unit, "sample_cell_group")) {
    glue::glue(
      "本分析先按 `{obj$sample_col}` 与 `{obj$cell_group_col}` 对细胞表达值取均值，",
      "以样本内细胞群作为通路评分单位，得到每个样本-细胞群组合的 pathway-level score。"
    )
  } else {
    glue::glue(
      "本分析先按 `{obj$sample_col}` 与 `{obj$condition_col}` 对细胞表达值取均值，",
      "以样本分组作为通路评分单位，得到每个样本-分组组合的 pathway-level score。"
    )
  }

  text <- glue::glue(
    "基于目标单细胞群的标准化表达矩阵，从 Seurat 对象的 `{assay_text}` assay（{layer_text}）提取表达值。",
    "{unit_text}",
    "采用 R 包 `GSVA` ⟦pkgInfo('GSVA')⟧ 的 {score_method} 方法进行通路活性评分 (PMID: 23323831)，",
    "将 {n_gene} 个基因 × {n_unit} 个评分单位的表达矩阵转换为 {n_gene_set} 个基因集 × {n_unit} 个评分单位的 pathway-level score matrix。",
    "基因集来源为 {mode_text}，并限定匹配后的基因集大小为 {min_gs_size}–{max_gs_size} 个基因；GSVA 核密度参数设置为 kcdf = {kcdf}。"
  )

  return(text)
}

scpaFuns$format_score_snap <- function(mat_score,
  score_unit,
  data_unit_metadata
)
{
  text <- glue::glue(
    "通路活性评分共保留 {nrow(mat_score)} 个基因集，",
    "并在 {ncol(mat_score)} 个 {scpaFuns$get_score_unit_label(score_unit)} 评分单位上获得 pathway-level score。"
  )

  if ("n_cell" %in% colnames(data_unit_metadata) && !identical(score_unit, "cell")) {
    text <- glue::glue(
      "{text} 每个评分单位包含细胞数的中位数为 {stats::median(data_unit_metadata$n_cell, na.rm = TRUE)}。"
    )
  }

  return(as.character(text))
}

scpaFuns$run_pathway_score <- function(mat_expr,
  gene_sets,
  score_method = "gsva",
  kcdf = "auto",
  mx.diff = TRUE,
  abs.ranking = FALSE,
  verbose = TRUE,
  ...
)
{
  version_gsva <- packageVersion("GSVA")
  args_extra <- list(...)
  args_gsva_call <- list(verbose = verbose)

  if ("BPPARAM" %in% names(args_extra)) {
    args_gsva_call$BPPARAM <- args_extra$BPPARAM
    args_extra$BPPARAM <- NULL
  }

  if (isTRUE(verbose)) {
    cli::cli_alert_info(glue::glue(
      "GSVA engine input: {nrow(mat_expr)} genes x {ncol(mat_expr)} units; {length(gene_sets)} gene sets."
    ))
  }

  if (version_gsva >= "1.50.0") {
    if (identical(score_method, "ssgsea")) {
      args_param <- c(
        list(
          exprData = mat_expr,
          geneSets = gene_sets,
          normalize = TRUE,
          verbose = verbose
        ),
        args_extra
      )
      args_param <- scpaFuns$keep_formal_args(GSVA::ssgseaParam, args_param)
      param <- do.call(GSVA::ssgseaParam, args_param)
    } else {
      args_param <- list(
        exprData = mat_expr,
        geneSets = gene_sets,
        kcdf = kcdf,
        absRanking = abs.ranking,
        verbose = verbose
      )
      param_formals <- names(formals(GSVA::gsvaParam))
      if ("maxDiff" %in% param_formals) {
        args_param$maxDiff <- mx.diff
      } else if ("mxDiff" %in% param_formals) {
        args_param$mxDiff <- mx.diff
      }
      if ("sparse" %in% param_formals && inherits(mat_expr, "sparseMatrix") &&
          !("sparse" %in% names(args_extra))) {
        args_param$sparse <- TRUE
      }
      args_param <- c(args_param, args_extra)
      args_param <- scpaFuns$keep_formal_args(GSVA::gsvaParam, args_param)
      param <- do.call(GSVA::gsvaParam, args_param)
    }
    args_gsva_call <- c(list(param), args_gsva_call)
    mat_score <- do.call(GSVA::gsva, args_gsva_call)
  } else {
    args_gsva <- c(
      list(
        expr = mat_expr,
        gset.idx.list = gene_sets,
        method = score_method,
        kcdf = kcdf,
        mx.diff = mx.diff,
        abs.ranking = abs.ranking,
        verbose = verbose
      ),
      args_extra
    )
    mat_score <- do.call(GSVA::gsva, args_gsva)
  }

  mat_score <- as.matrix(mat_score)
  return(mat_score)
}

scpaFuns$keep_formal_args <- function(fun, args)
{
  args <- args[!duplicated(names(args), fromLast = TRUE)]
  formal_names <- names(formals(fun))
  if ("..." %in% formal_names) {
    return(args)
  }
  args <- args[names(args) %in% formal_names]
  return(args)
}

scpaFuns$summarize_score_matrix <- function(mat_score)
{
  data_summary <- data.frame(
    ID = rownames(mat_score),
    mean_score = as.numeric(rowMeans(mat_score, na.rm = TRUE)),
    sd_score = as.numeric(apply(mat_score, 1L, stats::sd, na.rm = TRUE)),
    min_score = as.numeric(apply(mat_score, 1L, min, na.rm = TRUE)),
    max_score = as.numeric(apply(mat_score, 1L, max, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
  data_summary <- data_summary[order(data_summary$sd_score, decreasing = TRUE), , drop = FALSE]
  rownames(data_summary) <- NULL

  return(data_summary)
}

scpaFuns$run_limma_diff <- function(mat_score,
  data_metadata,
  compare_col,
  compare_by = "cell_group",
  levels = NULL
)
{
  if (identical(compare_by, "condition") && !is.null(levels) && length(levels) >= 2L) {
    data_diff <- scpaFuns$run_limma_one_contrast(
      mat_score = mat_score,
      data_group = data_metadata[[compare_col]],
      contrast_name = paste0(levels[1L], "_vs_", levels[2L]),
      positive_level = levels[1L],
      negative_level = levels[2L]
    )
    data_diff$compare_col <- compare_col
    data_diff$compare_by <- compare_by
    return(data_diff)
  }

  groups <- unique(as.character(data_metadata[[compare_col]]))
  groups <- groups[!is.na(groups) & nzchar(groups)]

  lst_diff <- lapply(groups,
    function(group_one) {
      data_diff <- scpaFuns$run_limma_one_vs_rest(
        mat_score = mat_score,
        data_group = data_metadata[[compare_col]],
        group_one = group_one
      )
      data_diff$compare_col <- compare_col
      data_diff$compare_by <- compare_by
      return(data_diff)
    })

  data_diff <- dplyr::bind_rows(lst_diff)
  return(data_diff)
}

scpaFuns$run_limma_one_contrast <- function(mat_score,
  data_group,
  contrast_name,
  positive_level,
  negative_level
)
{
  data_group <- as.character(data_group)
  keep <- data_group %in% c(positive_level, negative_level)
  mat_use <- mat_score[, keep, drop = FALSE]
  data_group <- data_group[keep]
  data_group <- factor(data_group, levels = c(negative_level, positive_level))

  design <- stats::model.matrix(~ data_group)
  fit <- limma::lmFit(mat_use, design)
  fit <- limma::eBayes(fit)
  data_top <- limma::topTable(fit, coef = 2L, number = Inf, sort.by = "t")
  data_top <- scpaFuns$format_limma_table(data_top)
  data_top$contrast <- contrast_name
  data_top$group <- positive_level
  data_top$reference <- negative_level
  return(data_top)
}

scpaFuns$run_limma_one_vs_rest <- function(mat_score,
  data_group,
  group_one
)
{
  data_group <- ifelse(as.character(data_group) == group_one, group_one, "Other")
  data_group <- factor(data_group, levels = c("Other", group_one))

  design <- stats::model.matrix(~ data_group)
  fit <- limma::lmFit(mat_score, design)
  fit <- limma::eBayes(fit)
  data_top <- limma::topTable(fit, coef = 2L, number = Inf, sort.by = "t")
  data_top <- scpaFuns$format_limma_table(data_top)
  data_top$contrast <- paste0(group_one, "_vs_Other")
  data_top$group <- group_one
  data_top$reference <- "Other"
  return(data_top)
}

scpaFuns$format_limma_table <- function(data_top)
{
  data_top <- as.data.frame(data_top, stringsAsFactors = FALSE)
  data_top$ID <- rownames(data_top)
  rownames(data_top) <- NULL

  if ("logFC" %in% colnames(data_top)) {
    data_top$delta_score <- data_top$logFC
  }

  data_top$direction <- ifelse(data_top$delta_score >= 0, "Up", "Down")
  col_order <- intersect(
    c("ID", "contrast", "group", "reference", "delta_score", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "direction"),
    colnames(data_top)
  )
  col_other <- setdiff(colnames(data_top), col_order)
  data_top <- data_top[, c(col_order, col_other), drop = FALSE]

  return(data_top)
}

scpaFuns$annotate_pathway_terms <- function(data_diff, db_anno = NULL)
{
  data_diff <- as.data.frame(data_diff, stringsAsFactors = FALSE)

  data_diff$Pathway_ID <- data_diff$ID
  data_diff$Pathway_Name <- data_diff$ID
  data_diff$Pathway_Description <- data_diff$ID

  if (!is.null(db_anno) && nrow(db_anno) > 0L) {
    db_anno <- as.data.frame(db_anno, stringsAsFactors = FALSE)
    id_col <- intersect(c("gs_id", "ID", "id"), colnames(db_anno))
    name_col <- intersect(c("gs_name", "gs_exact_source", "term", "pathway"), colnames(db_anno))
    desc_col <- intersect(c("gs_description", "description", "Description"), colnames(db_anno))

    if (length(id_col) > 0L) {
      id_col <- id_col[1L]
      data_anno <- data.frame(
        Pathway_ID = as.character(db_anno[[id_col]]),
        stringsAsFactors = FALSE
      )
      if (length(name_col) > 0L) {
        data_anno$Pathway_Name <- as.character(db_anno[[name_col[1L]]])
      }
      if (length(desc_col) > 0L) {
        data_anno$Pathway_Description <- as.character(db_anno[[desc_col[1L]]])
      }
      data_anno <- data_anno[!is.na(data_anno$Pathway_ID) & nzchar(data_anno$Pathway_ID), , drop = FALSE]
      data_anno <- data_anno[!duplicated(data_anno$Pathway_ID), , drop = FALSE]
      idx <- match(data_diff$ID, data_anno$Pathway_ID)

      if ("Pathway_Name" %in% colnames(data_anno)) {
        hit <- !is.na(idx) & !is.na(data_anno$Pathway_Name[idx]) & nzchar(data_anno$Pathway_Name[idx])
        data_diff$Pathway_Name[hit] <- data_anno$Pathway_Name[idx[hit]]
      }
      if ("Pathway_Description" %in% colnames(data_anno)) {
        hit <- !is.na(idx) & !is.na(data_anno$Pathway_Description[idx]) & nzchar(data_anno$Pathway_Description[idx])
        data_diff$Pathway_Description[hit] <- data_anno$Pathway_Description[idx[hit]]
      }
    }
  }

  data_diff$Pathway <- scpaFuns$clean_pathway_label(data_diff$Pathway_Name)
  desc_ok <- !is.na(data_diff$Pathway_Description) & nzchar(data_diff$Pathway_Description) &
    data_diff$Pathway_Description != data_diff$Pathway_ID
  data_diff$Pathway_Label <- data_diff$Pathway
  data_diff$Pathway_Label[desc_ok] <- data_diff$Pathway_Description[desc_ok]
  data_diff$Pathway_Label <- scpaFuns$clean_pathway_label(data_diff$Pathway_Label)

  col_first <- intersect(
    c("ID", "Pathway", "Pathway_Label", "Pathway_Name", "Pathway_Description", "contrast", "group", "reference"),
    colnames(data_diff)
  )
  col_other <- setdiff(colnames(data_diff), col_first)
  data_diff <- data_diff[, c(col_first, col_other), drop = FALSE]

  return(data_diff)
}

scpaFuns$clean_pathway_label <- function(x)
{
  x <- as.character(x)
  x <- gsub("^HALLMARK_", "HALLMARK ", x)
  x <- gsub("^REACTOME_", "REACTOME ", x)
  x <- gsub("^KEGG_", "KEGG ", x)
  x <- gsub("^WP_", "WIKIPATHWAYS ", x)
  x <- gsub("_", " ", x)
  x <- stringr::str_squish(x)
  return(x)
}

scpaFuns$get_threshold_text <- function(use_p,
  cutoff,
  cutoff.delta_score = 0
)
{
  if (!is.null(cutoff.delta_score) && is.finite(cutoff.delta_score) &&
      cutoff.delta_score > 0) {
    return(glue::glue("{use_p} &lt; {cutoff} 且 |delta_score| ≥ {cutoff.delta_score}（pathway-level score 差异阈值）"))
  }
  return(glue::glue("{use_p} &lt; {cutoff}"))
}

scpaFuns$add_diff_threshold_status <- function(data_diff,
  use_p = "adj.P.Val",
  cutoff = 0.05,
  cutoff.delta_score = 0
)
{
  data_diff <- as.data.frame(data_diff, stringsAsFactors = FALSE)

  if (use_p %in% colnames(data_diff)) {
    data_diff$pass_p <- !is.na(data_diff[[use_p]]) & data_diff[[use_p]] < cutoff
  } else {
    data_diff$pass_p <- TRUE
  }

  if ("delta_score" %in% colnames(data_diff) && !is.null(cutoff.delta_score) &&
      is.finite(cutoff.delta_score) && cutoff.delta_score > 0) {
    data_diff$pass_delta_score <- !is.na(data_diff$delta_score) &
      abs(data_diff$delta_score) >= cutoff.delta_score
  } else {
    data_diff$pass_delta_score <- TRUE
  }

  data_diff$pass_filter <- data_diff$pass_p & data_diff$pass_delta_score
  return(data_diff)
}

scpaFuns$filter_diff_by_thresholds <- function(data_diff,
  use_p = "adj.P.Val",
  cutoff = 0.05,
  cutoff.delta_score = 0
)
{
  if (!"pass_filter" %in% colnames(data_diff)) {
    data_diff <- scpaFuns$add_diff_threshold_status(
      data_diff = data_diff,
      use_p = use_p,
      cutoff = cutoff,
      cutoff.delta_score = cutoff.delta_score
    )
  }
  data_diff <- data_diff[data_diff$pass_filter, , drop = FALSE]
  return(data_diff)
}

scpaFuns$filter_diff_report <- function(data_diff,
  use_p = "adj.P.Val",
  cutoff = 0.05,
  cutoff.delta_score = 0,
  top_n = 20L,
  rank_by = "abs_t",
  select_mode = c("balanced", "overall")
)
{
  select_mode <- match.arg(select_mode)
  if (nrow(data_diff) == 0L) {
    return(data_diff)
  }

  data_sig <- scpaFuns$filter_diff_by_thresholds(
    data_diff = data_diff,
    use_p = use_p,
    cutoff = cutoff,
    cutoff.delta_score = cutoff.delta_score
  )

  if (nrow(data_sig) == 0L) {
    return(data_sig)
  }

  if (identical(select_mode, "balanced")) {
    data_sig$rank_score <- scpaFuns$get_direction_rank_score(data_sig, rank_by)
    top_each <- as.integer(ceiling(top_n / 2L))
    lst_contrast <- split(data_sig, data_sig$contrast)
    data_sig <- dplyr::bind_rows(lapply(lst_contrast,
      function(data_one) {
        data_up <- data_one[data_one$t > 0, , drop = FALSE]
        data_down <- data_one[data_one$t < 0, , drop = FALSE]
        if (nrow(data_up) > 0L) {
          data_up <- dplyr::arrange(data_up, dplyr::desc(rank_score))
          data_up <- utils::head(data_up, top_each)
        }
        if (nrow(data_down) > 0L) {
          data_down <- dplyr::arrange(data_down, dplyr::desc(rank_score))
          data_down <- utils::head(data_down, top_each)
        }
        dplyr::bind_rows(data_up, data_down)
      }))
  } else {
    data_sig$rank_score <- scpaFuns$get_diff_rank_score(data_sig, rank_by)
    data_sig <- dplyr::arrange(data_sig, contrast, dplyr::desc(rank_score))
    data_sig <- dplyr::group_by(data_sig, contrast)
    data_sig <- dplyr::slice_head(data_sig, n = top_n)
    data_sig <- dplyr::ungroup(data_sig)
  }

  data_sig <- dplyr::arrange(data_sig, contrast, dplyr::desc(t))
  data_sig <- as.data.frame(data_sig, stringsAsFactors = FALSE)

  return(data_sig)
}

scpaFuns$get_direction_rank_score <- function(data_diff, rank_by)
{
  if (identical(rank_by, "delta_score") || identical(rank_by, "abs_delta_score")) {
    return(abs(data_diff$delta_score))
  }
  return(abs(data_diff$t))
}

scpaFuns$get_diff_rank_score <- function(data_diff, rank_by)
{
  if (identical(rank_by, "t")) {
    return(data_diff$t)
  }
  if (identical(rank_by, "delta_score")) {
    return(data_diff$delta_score)
  }
  if (identical(rank_by, "abs_delta_score")) {
    return(abs(data_diff$delta_score))
  }
  return(abs(data_diff$t))
}

scpaFuns$get_rank_label <- function(rank_by)
{
  if (identical(rank_by, "t")) {
    return("limma t 值")
  }
  if (identical(rank_by, "delta_score")) {
    return("通路评分差异")
  }
  if (identical(rank_by, "abs_delta_score")) {
    return("通路评分差异绝对值")
  }
  return("limma t 值绝对值")
}

scpaFuns$get_plot_select_text <- function(top_n, rank_by, select_mode)
{
  if (identical(select_mode, "balanced")) {
    top_each <- as.integer(ceiling(top_n / 2L))
    return(glue::glue(
      "{scpaFuns$get_rank_label(rank_by)} 在上调和下调方向分别筛选不超过 {top_each} 个通路"
    ))
  }
  return(glue::glue("{scpaFuns$get_rank_label(rank_by)} 筛选前 {top_n} 个通路"))
}

scpaFuns$select_diff_plot_data <- function(data_diff,
  use_p = "adj.P.Val",
  cutoff = 0.05,
  cutoff.delta_score = 0,
  top_n = 20L,
  rank_by = "abs_t",
  select_mode = c("balanced", "overall")
)
{
  if (nrow(data_diff) == 0L) {
    return(data_diff)
  }

  select_mode <- match.arg(select_mode)

  data_plot <- scpaFuns$filter_diff_report(
    data_diff = data_diff,
    use_p = use_p,
    cutoff = cutoff,
    cutoff.delta_score = cutoff.delta_score,
    top_n = top_n,
    rank_by = rank_by,
    select_mode = select_mode
  )

  if (!"pass_filter" %in% colnames(data_plot)) {
    data_plot <- scpaFuns$add_diff_threshold_status(
      data_diff = data_plot,
      use_p = use_p,
      cutoff = cutoff,
      cutoff.delta_score = cutoff.delta_score
    )
  }
  data_plot$Group <- ifelse(
    data_plot$pass_filter & data_plot$t > 0, "Up",
    ifelse(data_plot$pass_filter & data_plot$t < 0, "Down", "Not")
  )
  data_plot$Group <- factor(data_plot$Group, levels = c("Down", "Not", "Up"))

  data_plot <- dplyr::arrange(data_plot, contrast, dplyr::desc(t))
  data_plot <- dplyr::group_by(data_plot, contrast)
  data_plot$plot_label <- make.unique(stringr::str_wrap(data_plot$Pathway_Label, 60L))
  data_plot$plot_label <- factor(data_plot$plot_label, levels = rev(unique(data_plot$plot_label)))
  data_plot <- dplyr::ungroup(data_plot)
  data_plot <- as.data.frame(data_plot, stringsAsFactors = FALSE)

  return(data_plot)
}

scpaFuns$plot_pathway_t_heatmap <- function(data_plot, compare_by)
{
  data_plot <- as.data.frame(data_plot, stringsAsFactors = FALSE)
  data_plot$Group <- factor(data_plot$Group, levels = c("Down", "Not", "Up"))

  p <- ggplot2::ggplot(data_plot, ggplot2::aes(x = t, y = plot_label, fill = Group)) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
    ggplot2::scale_fill_manual(
      values = c(Down = "#00B050", Not = "#BFBFBF", Up = "#FF0000"),
      drop = FALSE
    ) +
    ggplot2::labs(
      x = "t value of GSVA score",
      y = NULL,
      fill = "Group"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 7),
      legend.position = "right"
    )

  if (length(unique(data_plot$contrast)) > 1L) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(contrast), scales = "free_y")
  }

  return(p)
}

scpaFuns$format_diff_method_text <- function(compare_by,
  compare_col,
  use_p,
  cutoff,
  cutoff.delta_score = 0,
  top_n,
  rank_by,
  select_mode
)
{
  compare_text <- if (identical(compare_by, "condition")) {
    glue::glue("以 `{compare_col}` 定义的疾病/处理分组为比较对象")
  } else {
    glue::glue("以 `{compare_col}` 定义的每个细胞群分别对比其他剩余细胞群")
  }

  select_text <- if (identical(select_mode, "balanced")) {
    top_each <- as.integer(ceiling(top_n / 2L))
    glue::glue(
      "在每个比较中，分别从评分升高和评分降低方向各选取不超过 {top_each} 个代表性通路用于结果展示"
    )
  } else {
    glue::glue(
      "在每个比较中，根据 {scpaFuns$get_rank_label(rank_by)} 选取前 {top_n} 个代表性通路用于结果展示"
    )
  }

  threshold_text <- scpaFuns$get_threshold_text(
    use_p = use_p,
    cutoff = cutoff,
    cutoff.delta_score = cutoff.delta_score
  )

  text <- glue::glue(
    "基于通路活性评分矩阵，{compare_text}，采用 R 包 `limma` ⟦pkgInfo('limma')⟧ 对 pathway-level score 进行差异分析。",
    "差异结果采用 Benjamini-Hochberg 方法进行多重检验校正，并以 {threshold_text} 筛选差异通路。",
    "{select_text}。delta_score 表示目标组相对于参照组的 pathway-level score 差异，对应 limma 输出中的 logFC 列；当设置 delta_score 阈值时，该阈值用于限制通路评分差异幅度。"
  )

  return(text)
}

scpaFuns$format_diff_snap <- function(data_diff,
  data_report,
  compare_col,
  use_p,
  cutoff,
  cutoff.delta_score = 0
)
{
  if (nrow(data_report) == 0L) {
    return("")
  }

  contrasts <- unique(as.character(data_report$contrast))
  contrasts <- contrasts[!is.na(contrasts) & nzchar(contrasts)]

  texts <- vapply(head(contrasts, 3L), FUN.VALUE = character(1L),
    function(one_contrast) {
      data_all <- data_diff[data_diff$contrast == one_contrast, , drop = FALSE]
      data_sig <- scpaFuns$filter_diff_by_thresholds(
        data_diff = data_all,
        use_p = use_p,
        cutoff = cutoff,
        cutoff.delta_score = cutoff.delta_score
      )
      n_up <- sum(data_sig$delta_score > 0, na.rm = TRUE)
      n_down <- sum(data_sig$delta_score < 0, na.rm = TRUE)
      n_total <- nrow(data_sig)

      data_one <- data_report[data_report$contrast == one_contrast, , drop = FALSE]
      data_up <- data_one[data_one$delta_score > 0, , drop = FALSE]
      data_down <- data_one[data_one$delta_score < 0, , drop = FALSE]

      up_text <- ""
      if (nrow(data_up) > 0L) {
        data_up <- data_up[order(data_up$t, decreasing = TRUE), , drop = FALSE]
        up_text <- glue::glue(
          "评分升高最明显的代表性通路为 {data_up$Pathway_Label[1L]}（t = {round(data_up$t[1L], 2L)}，delta_score = {round(data_up$delta_score[1L], 3L)}）"
        )
      }

      down_text <- ""
      if (nrow(data_down) > 0L) {
        data_down <- data_down[order(data_down$t, decreasing = FALSE), , drop = FALSE]
        down_text <- glue::glue(
          "评分降低最明显的代表性通路为 {data_down$Pathway_Label[1L]}（t = {round(data_down$t[1L], 2L)}，delta_score = {round(data_down$delta_score[1L], 3L)}）"
        )
      }

      lead_text <- paste(c(up_text, down_text), collapse = "；")
      lead_text <- gsub("^；|；$", "", lead_text)
      if (!nzchar(lead_text)) {
        one <- data_one[order(abs(data_one$t), decreasing = TRUE), , drop = FALSE][1L, , drop = FALSE]
        lead_text <- glue::glue(
          "t 值绝对值最大的代表性通路为 {one$Pathway_Label}（t = {round(one$t, 2L)}，delta_score = {round(one$delta_score, 3L)}）"
        )
      }

      threshold_text <- scpaFuns$get_threshold_text(
        use_p = use_p,
        cutoff = cutoff,
        cutoff.delta_score = cutoff.delta_score
      )
      sig_text <- glue::glue(
        "{one_contrast} 共检出 {n_total} 条差异通路（Up = {n_up}，Down = {n_down}，{threshold_text}）"
      )
      if (n_total == 0L) {
        sig_text <- glue::glue(
          "{one_contrast} 在 {threshold_text} 阈值下未检出显著差异通路"
        )
      }

      return(as.character(glue::glue("{sig_text}；{lead_text}。")))
    })

  text <- paste(texts, collapse = "")
  return(text)
}
