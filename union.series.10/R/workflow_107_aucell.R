# ==========================================================================
# workflow of aucell
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_aucell <- setClass("job_aucell", 
  contains = c("job"),
  prototype = prototype(
    pg = "aucell",
    info = c("https://www.bioconductor.org/packages/release/bioc/vignettes/AUCell/inst/doc/AUCell.html"),
    cite = "[@SCENIC_single_Aibar_2017]",
    method = "",
    tag = "aucell",
    analysis = "AUCell 识别细胞的基因集活性"
    ))

setGeneric("asjob_aucell",
   function(x, ...) standardGeneric("asjob_aucell"))

aucellFuns <- new.env(parent = emptyenv())

setMethod("asjob_aucell", signature = c(x = "job_seurat"),
  function(x, sets, name = names(sets)[1], join = TRUE,
    sets_feature = NULL, assay = SeuratObject::DefaultAssay(object(x)),
    id_kegg = NULL, ...)
  {
    mtx <- Seurat::GetAssayData(object(x), assay = assay, layer = "data")
    if (is.null(mtx) || is.null(rownames(mtx))) {
      stop('is.null(mtx) || is.null(rownames(mtx)).')
    }
    if (!is.null(id_kegg)) {
      data_kegg <- geneFuns$get_kegg_pathway_gene_table(
        pathway_ids = id_kegg,
        species = "hsa"
      )
      sets <- split(data_kegg$symbol, data_kegg$pathway_name)
      sets <- as_feature(sets, name)
    }
    if (!is(sets, "feature")) {
      stop('!is(sets, "feature").')
    }
    if (is(sets, "feature_char")) {
      sets <- join(sets, name)
    }
    if (missing(name)) {
      message(glue::glue("Missing `name`, will use: {name}"))
    }
    snap <- stat_features(sets, name, join, "sets")
    methodAdd_onExit("x", "{snap}以该基因集作为 AUCell 输入。")
    pr <- params(x)
    if (is.null(pr$metadata)) {
      stop('is.null(pr$metadata).')
    }
    x <- job_aucell(mtx, sets, ...)
    x@params <- append(x@params, pr)
    x$.feature_genesSets <- sets
    x$name <- name
    return(x)
  })

job_aucell <- function(mtx, sets)
{
  if (!is(mtx, "dgCMatrix")) {
    stop('!is(mtx, "dgCMatrix").')
  }
  if (!is(sets, "GeneSetCollection")) {
    sets <- as_collection(sets)
  }
  gids <- unique(unlist(GSEABase::geneIds(sets)))
  isIns <- gids %in% rownames(mtx)
  message(glue::glue("Has genes: {try_snap(isIns)}"))
  if (all(!isIns)) {
    stop('all(!isIns).')
  }
  x <- .job_aucell(object = mtx)
  x <- methodAdd(x, "AUCell 是一种基于单细胞转录组数据评估基因集活性的分析方法，其主要目的是在单细胞分辨率下量化预定义基因集（如信号通路、转录因子靶基因集或细胞状态特征基因）的活跃程度。该方法通过对每个细胞内基因表达进行排序，并计算目标基因集在高表达基因中的富集面积（AUC score），从而评估该基因集在不同细胞中的相对活性。")
  x <- methodAdd(x, "以 R 包 `AUCell` ⟦pkgInfo('AUCell')⟧ {cite_show('SCENIC_single_Aibar_2017')} 识别单细胞数据集的基因集调控活性。")
  x$sets <- sets
  x$gids <- gids
  return(x)
}

setMethod("step0", signature = c(x = "job_aucell"),
  function(x){
    step_message("Prepare your data with function `job_aucell`.")
  })

setMethod("step1", signature = c(x = "job_aucell"),
  function(x, workers = NULL, group.by = x$group.by,
    rerun = FALSE, fun_name = function(x) s(x, "^HALLMARK_", ""),
    score_methods = "AUCell", score_cutoff = 0.2, top_cell_prop = 0.1,
    group_test = TRUE, group_col = "group", box_score_methods = "AUCell",
    box_max_facet_plot = 50L, box_stat_only = NULL)
  {
    step_message("Running gene set scoring...")
    box_stat_only_expr <- if (is.null(box_stat_only)) {
      "NULL"
    } else {
      as.character(isTRUE(box_stat_only))
    }
    if (is.remote(x)) {
      x <- run_job_remote(x, wait = 3L,
        {
          x <- step1(
            x,
            workers = "{workers}"
          )
        }
      )
      return(x)
    }

    score_methods <- aucellFuns$resolve_score_methods(score_methods)
    fun_show <- function(string) stringr::str_wrap(gs(string, "_", " "), 20)

    fun_aucell <- function(...) {
      if (!is.null(workers)) {
        workers <- e(BiocParallel::MulticoreParam(workers))
      }
      e(AUCell::AUCell_run(object(x), x$sets, BPPARAM = workers))
    }

    res_aucell <- expect_local_data(
      "tmp", "AUcell", fun_aucell,
      list(colnames(object(x)), sig(x), names(x$sets), x$gids),
      rerun = rerun
    )

    x$res_aucell <- e(AUCell::getAUC(res_aucell))
    if (!is.null(fun_name)) {
      rownames(x$res_aucell) <- fun_name(rownames(x$res_aucell))
    }

    lst_score <- list(AUCell = x$res_aucell)
    other_methods <- setdiff(score_methods, c("AUCell", "Scoring"))

    if (length(other_methods) > 0L) {
      fun_score <- function(...) {
        aucellFuns$run_score_methods(
          mtx = object(x),
          sets = x$sets,
          methods = other_methods,
          score_cutoff = score_cutoff
        )
      }

      lst_other_score <- expect_local_data(
        "tmp", "gene_set_multi_scores", fun_score,
        list(
          colnames(object(x)),
          rownames(object(x)),
          sig(x),
          names(x$sets),
          x$gids,
          other_methods,
          score_cutoff
        ),
        rerun = rerun
      )

      if (length(lst_other_score) > 0L) {
        if (!is.null(fun_name)) {
          lst_other_score <- lapply(
            lst_other_score,
            aucellFuns$rename_score_matrix,
            fun_name = fun_name
          )
        }
        lst_score <- c(lst_score, lst_other_score)
      }
    }

    if ("Scoring" %in% score_methods) {
      lst_score$Scoring <- aucellFuns$get_composite_scoring_score(
        lst_score,
        component_methods = setdiff(score_methods, "Scoring")
      )
    }

    score_methods_done <- names(lst_score)
    data_score_long <- aucellFuns$as_score_long(
      lst_score,
      score_cutoff = score_cutoff,
      top_cell_prop = top_cell_prop
    )
    data_key_cells <- dplyr::filter(data_score_long, is_top_cell)

    if (is.null(x$lst_all_others)) {
      x$lst_all_others <- list()
    }
    x$lst_all_others$gene_set_scoring <- list(
      scores = lst_score,
      score_methods = score_methods_done,
      score_long = data_score_long,
      key_cells = data_key_cells,
      score_cutoff = score_cutoff,
      top_cell_prop = top_cell_prop
    )


    other_methods_done <- setdiff(score_methods_done, "AUCell")

    if (length(other_methods_done) == 0L) {
      x <- methodAdd(
        x,
        "采用 AUCell 对预定义候选基因集在单细胞水平进行活性评估。AUCell 通过对每个细胞内基因表达排序，计算候选基因集在高表达基因中的富集面积，从而获得单细胞层面的基因集活性评分。"
      )
    } else {
      other_methods_report <- setdiff(other_methods_done, "Scoring")
      if ("Scoring" %in% score_methods_done) {
        x <- methodAdd(
          x,
          "采用多算法单细胞基因集评分策略对预定义候选基因集在单细胞水平进行活性评估。本分析实际计算的评分方法包括：{bind(score_methods_done)}。其中，AUCell 用于计算基于基因表达排序的基因集活性评分；{bind(other_methods_report)} 用于从不同算法角度补充评估候选基因集活性；Scoring 为基于上述评分结果构建的综合评分，即对各方法的同一基因集评分进行 0-1 归一化后取平均，用于评估不同算法所得候选基因集活性的整体一致性。"
        )
      } else {
        x <- methodAdd(
          x,
          "采用多算法单细胞基因集评分策略对预定义候选基因集在单细胞水平进行活性评估。本分析实际计算的评分方法包括：{bind(score_methods_done)}。其中，AUCell 用于计算基于基因表达排序的基因集活性评分；{bind(other_methods_done)} 用于从不同算法角度补充评估候选基因集活性。"
        )
      }
    }

    if ("Scoring" %in% score_methods_done) {
      x <- methodAdd(
        x,
        "基于单细胞评分结果，记录每个候选基因集评分位于前 {top_cell_prop * 100}% 的细胞作为高评分细胞；对于综合 Scoring 结果，同时保留 score &gt;= {score_cutoff} 的阈值标记，用于辅助展示候选基因集在细胞群体中的分布。"
      )
    } else {
      x <- methodAdd(
        x,
        "基于单细胞评分结果，记录每个候选基因集评分位于前 {top_cell_prop * 100}% 的细胞作为高评分细胞，用于辅助展示候选基因集在细胞群体中的分布。"
      )
    }



    if (isTRUE(group_test)) {
      data_box_score <- aucellFuns$prepare_group_score_data(
        data_score_long = data_score_long,
        metadata = x$metadata,
        group_col = group_col,
        celltype_col = group.by,
        group_levels = x$levels,
        score_methods = box_score_methods,
        score_methods_done = score_methods_done
      )

      if (is.null(data_box_score) || nrow(data_box_score$data) == 0L) {
        warning(
          "Skip cell type-specific group comparison of gene set scores because valid score data was not available, cell type annotation was missing, or more than one gene set was provided.",
          call. = FALSE
        )
      } else if (!exists(".map_boxplot2", mode = "function")) {
        warning(
          "Skip group comparison of gene set scores because `.map_boxplot2()` was not found.",
          call. = FALSE
        )
      } else {
        p.box_score <- .map_boxplot2(
          data_box_score$data,
          pvalue = TRUE,
          x = "group",
          y = "value",
          ids = "var",
          xlab = "Group",
          ylab = "Gene set score",
          test = "wilcox.test",
          max_facet_plot = box_max_facet_plot,
          stat_only = box_stat_only
        )

        pvalue_box_score <- attr(p.box_score, "pvalue")
        compare_box_score <- attr(p.box_score, "compare")

        layout <- wrap_layout(NULL, length(pvalue_box_score))
        p.box_score <- set_lab_legend(
          add(layout, wrap(p.box_score)),
          glue::glue("{x@sig} boxplot of gene set enrichment score"),
          glue::glue(
            "细胞类型内基因集富集评分箱形图|||该图展示各细胞类型内候选基因集评分在 {bind(data_box_score$group_levels)} 组间的细胞水平分布差异。",
            "每个分面对应一个细胞类型；箱线图基于单细胞评分绘制，组间比较采用 Wilcoxon 秩和检验。"
          )
        )

        x <- plotsAdd(x, p.box_score)

        x$lst_all_others$gene_set_scoring$group_test <- list(
          data = data_box_score$data,
          plot = p.box_score,
          pvalue = pvalue_box_score,
          compare = compare_box_score,
          group_col = group_col,
          celltype_col = data_box_score$celltype_col,
          gene_set = data_box_score$gene_set,
          group_levels = data_box_score$group_levels,
          score_methods = data_box_score$score_methods
        )

        if (length(data_box_score$score_methods) == 1L && data_box_score$score_methods[1L] == "AUCell") {
          x <- methodAdd(
            x,
            "基于 metadata 中的 {group_col} 分组信息和 {data_box_score$celltype_col} 细胞类型注释，进一步在各细胞类型内分别比较候选基因集（{data_box_score$gene_set}）AUCell 评分在 {data_box_score$group_levels[1]} 组和 {data_box_score$group_levels[2]} 组之间的细胞水平分布差异。组间比较采用 Wilcoxon 秩和检验，P &lt; 0.05 认为差异具有统计学意义。"
          )
        } else {
          x <- methodAdd(
            x,
            "基于 metadata 中的 {group_col} 分组信息和 {data_box_score$celltype_col} 细胞类型注释，进一步在各细胞类型内分别比较候选基因集（{data_box_score$gene_set}）评分在 {data_box_score$group_levels[1]} 组和 {data_box_score$group_levels[2]} 组之间的细胞水平分布差异。本分析纳入的评分方法包括 {bind(data_box_score$score_methods)}。组间比较采用 Wilcoxon 秩和检验，P &lt; 0.05 认为差异具有统计学意义。"
          )
        }

        if (exists(".stat_compare_by_pvalue", mode = "function")) {
          x <- snapAdd(
            x,
            .stat_compare_by_pvalue(
              p.box_score,
              data_box_score$group_levels,
              "基因集",
              mode = "enrichment"
            )
          )
        } else {
          x <- snapAdd(
            x,
            "基于细胞类型内的单细胞评分分布，对候选基因集评分进行组间 Wilcoxon 秩和检验，并以箱线图展示其组间差异{aref(p.box_score)}。"
          )
        }
      }
    }

    if (!is.null(group.by)) {
      metadata <- x$metadata
      fun_mean <- function(...) {
        auc <- .get_auc_from_job_aucell(x)
        data <- cbind(
          metadata[, group.by, drop = FALSE],
          as.data.frame(auc)
        )
        data <- as.data.table(data)
        data <- data.table:::`[.data.table`(
          data, , lapply(.SD, mean), by = group.by,
          .SDcols = setdiff(names(data), group.by)
        )
        tibble::as_tibble(data)
      }

      data <- expect_local_data(
        "tmp", "aucell_mean", fun_mean,
        list(rownames(x$res_aucell), metadata$cell, group.by, x$gids),
        rerun = rerun
      )

      x$res_aucell_mean <- data
      data_score_mean <- aucellFuns$summarize_score_by_group(
        lst_score = lst_score,
        metadata = metadata,
        group.by = group.by
      )
      data_score_mean_long <- aucellFuns$as_score_mean_long(
        lst_score_mean = data_score_mean,
        group.by = group.by
      )

      x$lst_all_others$gene_set_scoring$score_mean <- data_score_mean
      x$lst_all_others$gene_set_scoring$score_mean_long <- data_score_mean_long


      x <- methodAdd(
        x,
        "基于细胞注释信息（{group.by}）对同一细胞群内所有细胞的基因集评分取平均，以评估不同细胞类型的整体候选基因集活性。"
      )

      if (nrow(data_score_mean_long) > 0L) {
        p.score_heatmap <- aucellFuns$plot_score_heatmap(
          data_score_mean_long,
          group.by = group.by
        )
        p.score_heatmap <- set_lab_legend(
          p.score_heatmap,
          glue::glue("{x@sig} Gene set score dot heatmap"),
          glue::glue(
            "基因集评分点图|||该图展示已计算基因集评分在不同细胞群中的平均活性分布。",
            "横坐标为评分方法，纵坐标为细胞群（{group.by}）；",
            "点颜色表示同一评分方法内、同一候选基因集在不同细胞群之间标准化后的相对平均评分，点大小表示平均评分。",
            "当包含多个候选基因集时，不同候选基因集以分面形式展示。"
          )
        )
        x <- plotsAdd(x, p.score_heatmap)
        x$lst_all_others$gene_set_scoring$p_score_heatmap <- p.score_heatmap
        x <- snapAdd(
          x,
          "基于细胞注释信息计算候选基因集在不同细胞群中的平均评分，并以点图展示其相对分布{aref(p.score_heatmap)}。"
        )
      }

      if (ncol(data) < 21L) {
        layout <- wrap_layout(NULL, ncol(data) - 1L, 3)
        data <- tidyr::pivot_longer(
          data, -!!rlang::sym(group.by),
          names_to = "Function", values_to = "Activity"
        )
        if (!is.null(fun_show)) {
          data <- dplyr::mutate(data, Function = fun_show(Function))
        }
        p.aucell_mean <- ggplot(data, aes(x = reorder(!!rlang::sym(group.by), Activity), y = Activity)) +
          geom_col() +
          facet_wrap(~ Function, ncol = layout$ncol) +
          labs(x = "Cell types", y = "Activity") +
          coord_flip() +
          theme_minimal()
        p.aucell_mean <- set_lab_legend(
          add(layout, p.aucell_mean),
          glue::glue("{x@sig} Mean AUCell Activity"),
          glue::glue("各细胞类型 AUCell 功能活性|||每个分面代表一个独立的功能通路或生物学过程（Function），横坐标表示不同细胞群体（按平均活性值排序），纵坐标表示该群体的平均 AUCell 活性评分（Activity）。")
        )
        x <- snapAdd(
          x, "对于每个功能基因集，计算对应细胞群体的平均 AUCell 活性分数，并以分面柱状图形式展示{aref(p.aucell_mean)}。"
        )
        x <- plotsAdd(x, p.aucell_mean)
      }
    }
    return(x)
  })

aucellFuns$resolve_score_methods <- function(score_methods)
{
  base_methods <- c("AUCell", "UCell", "singscore", "ssGSEA", "AddModuleScore")
  all_methods <- c(base_methods, "Scoring")
  score_methods <- unique(as.character(score_methods))
  score_methods <- trimws(score_methods)

  if (length(score_methods) == 0L || any(is.na(score_methods))) {
    stop("`score_methods` must contain at least one method name.")
  }

  if (any(tolower(score_methods) %in% c("all", "multi", "multiple"))) {
    return(all_methods)
  }

  idx <- match(tolower(score_methods), tolower(all_methods))
  if (any(is.na(idx))) {
    stop("Unsupported score method: ", bind(score_methods[is.na(idx)]))
  }

  score_methods <- all_methods[idx]

  if ("Scoring" %in% score_methods) {
    score_methods <- unique(c(base_methods, "Scoring"))
  }

  score_methods
}


aucellFuns$get_gene_list <- function(sets, features)
{
  lst_gene <- GSEABase::geneIds(sets)
  if (is.null(names(lst_gene))) {
    names(lst_gene) <- names(sets)
  }
  lst_gene <- lapply(lst_gene, function(x) unique(as.character(x)))
  lst_gene <- lapply(lst_gene, function(x) intersect(x, features))
  is_valid <- vapply(lst_gene, length, integer(1L)) > 0L
  if (all(!is_valid)) {
    stop("No gene set has genes matched to the expression matrix.")
  }
  if (any(!is_valid)) {
    warning("Gene sets without matched genes were removed: ", bind(names(lst_gene)[!is_valid]))
  }
  lst_gene[is_valid]
}

aucellFuns$rename_score_matrix <- function(mat, fun_name = NULL)
{
  if (!is.null(fun_name)) {
    rownames(mat) <- fun_name(rownames(mat))
  }
  mat
}

aucellFuns$normalize_score_matrix <- function(score, set_names, cell_names,
  method_name)
{
  if (is.null(score)) {
    return(NULL)
  }
  if (is.data.frame(score)) {
    score <- as.matrix(score)
  }
  if (!is.matrix(score)) {
    score <- as.matrix(score)
  }
  rownames(score) <- sub("_UCell$", "", rownames(score))
  colnames(score) <- sub("_UCell$", "", colnames(score))
  rownames(score) <- sub("^AddModuleScore_", "", rownames(score))
  colnames(score) <- sub("^AddModuleScore_", "", colnames(score))
  if (all(cell_names %in% rownames(score))) {
    score <- t(score[cell_names, , drop = FALSE])
  } else if (all(cell_names %in% colnames(score))) {
    score <- score[, cell_names, drop = FALSE]
  } else {
    stop("Cannot align score matrix cells for method: ", method_name)
  }
  if (!is.null(set_names)) {
    rownames(score) <- sub("_UCell$", "", rownames(score))
    keep <- rownames(score) %in% set_names
    if (any(keep)) {
      score <- score[keep, , drop = FALSE]
    }
  }
  score
}

aucellFuns$run_score_methods <- function(mtx, sets, methods,
  score_cutoff = 0.2)
{
  lst_gene <- aucellFuns$get_gene_list(sets, rownames(mtx))
  methods <- aucellFuns$resolve_score_methods(methods)
  methods <- setdiff(methods, "AUCell")
  lst_score <- list()
  for (method in methods) {
    cli::cli_alert_info("Calculate gene set score by {.val {method}}.")
    score <- switch(
      method,
      UCell = aucellFuns$run_ucell_score(mtx, lst_gene),
      singscore = aucellFuns$run_singscore_score(mtx, lst_gene),
      ssGSEA = aucellFuns$run_ssgsea_score(mtx, lst_gene),
      AddModuleScore = aucellFuns$run_addmodulescore_score(mtx, lst_gene),
      stop("Unsupported score method: ", method)
    )
    score <- aucellFuns$normalize_score_matrix(
      score,
      set_names = names(lst_gene),
      cell_names = colnames(mtx),
      method_name = method
    )
    lst_score[[ method ]] <- score
  }
  lst_score
}

aucellFuns$run_ucell_score <- function(mtx, lst_gene)
{
  if (!requireNamespace("UCell", quietly = TRUE)) {
    stop("Package `UCell` is required for UCell scoring.")
  }
  fun <- get("ScoreSignatures_UCell", envir = asNamespace("UCell"))
  res <- try(fun(mtx, features = lst_gene), silent = TRUE)
  if (inherits(res, "try-error")) {
    res <- try(fun(expr.mat = mtx, features = lst_gene), silent = TRUE)
  }
  if (inherits(res, "try-error")) {
    stop("Failed to calculate UCell score.")
  }
  res
}

aucellFuns$run_singscore_score <- function(mtx, lst_gene)
{
  if (!requireNamespace("singscore", quietly = TRUE)) {
    stop("Package `singscore` is required for singscore scoring.")
  }
  data_expr <- as.matrix(mtx)
  data_rank <- singscore::rankGenes(data_expr)
  lst_score <- lapply(lst_gene, function(genes) {
    res <- singscore::simpleScore(data_rank, upSet = genes)
    if ("TotalScore" %in% colnames(res)) {
      return(res$TotalScore)
    }
    res[[ ncol(res) ]]
  })
  score <- do.call(rbind, lst_score)
  rownames(score) <- names(lst_gene)
  colnames(score) <- colnames(mtx)
  score
}

aucellFuns$run_ssgsea_score <- function(mtx, lst_gene)
{
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    stop("Package `GSVA` is required for ssGSEA scoring.")
  }
  data_expr <- as.matrix(mtx)
  res <- try(
    GSVA::gsva(
      expr = data_expr,
      gset.idx.list = lst_gene,
      method = "ssgsea",
      kcdf = "Gaussian",
      verbose = FALSE
    ),
    silent = TRUE
  )
  if (inherits(res, "try-error")) {
    param <- GSVA::ssgseaParam(data_expr, lst_gene)
    res <- GSVA::gsva(param, verbose = FALSE)
  }
  res
}

aucellFuns$run_addmodulescore_score <- function(mtx, lst_gene)
{
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package `Seurat` is required for AddModuleScore scoring.")
  }
  obj <- Seurat::CreateSeuratObject(counts = mtx)
  obj <- try(
    SeuratObject::SetAssayData(obj, assay = SeuratObject::DefaultAssay(obj), layer = "data", new.data = mtx),
    silent = TRUE
  )
  if (inherits(obj, "try-error")) {
    obj <- Seurat::CreateSeuratObject(counts = mtx)
    obj <- SeuratObject::SetAssayData(obj, assay = SeuratObject::DefaultAssay(obj), slot = "data", new.data = mtx)
  }
  prefix <- "AddModuleScore_"
  obj <- Seurat::AddModuleScore(
    object = obj,
    features = lst_gene,
    name = prefix,
    assay = SeuratObject::DefaultAssay(obj),
    slot = "data"
  )
  score_cols <- grep(paste0("^", prefix), colnames(obj@meta.data), value = TRUE)
  score_cols <- score_cols[seq_along(lst_gene)]
  score <- as.matrix(obj@meta.data[, score_cols, drop = FALSE])
  colnames(score) <- names(lst_gene)
  score
}

aucellFuns$scale_score_matrix01 <- function(mat)
{
  mat <- as.matrix(mat)
  out <- matrix(
    0,
    nrow = nrow(mat),
    ncol = ncol(mat),
    dimnames = dimnames(mat)
  )

  for (i in seq_len(nrow(mat))) {
    score <- as.numeric(mat[i, ])
    score_min <- min(score, na.rm = TRUE)
    score_max <- max(score, na.rm = TRUE)

    if (!is.finite(score_min) || !is.finite(score_max) || score_max == score_min) {
      out[i, ] <- 0
    } else {
      out[i, ] <- (score - score_min) / (score_max - score_min)
    }
  }

  out
}

aucellFuns$get_composite_scoring_score <- function(lst_score,
  component_methods = c("AUCell", "UCell", "singscore", "ssGSEA", "AddModuleScore"))
{
  component_methods <- intersect(component_methods, names(lst_score))
  component_methods <- setdiff(component_methods, "Scoring")

  if (length(component_methods) == 0L) {
    stop("No component score is available for composite Scoring.")
  }

  common_sets <- Reduce(
    intersect,
    lapply(component_methods, function(method) rownames(lst_score[[ method ]]))
  )
  common_cells <- Reduce(
    intersect,
    lapply(component_methods, function(method) colnames(lst_score[[ method ]]))
  )

  if (length(common_sets) == 0L || length(common_cells) == 0L) {
    stop("Cannot align component scores for composite Scoring.")
  }

  lst_scaled <- lapply(component_methods, function(method) {
    mat <- lst_score[[ method ]]
    mat <- mat[common_sets, common_cells, drop = FALSE]
    aucellFuns$scale_score_matrix01(mat)
  })

  score <- Reduce(`+`, lst_scaled) / length(lst_scaled)
  rownames(score) <- common_sets
  colnames(score) <- common_cells

  score
}


aucellFuns$as_score_long <- function(lst_score, score_cutoff = 0.2,
  top_cell_prop = 0.1)
{
  lst_data <- lapply(names(lst_score), function(method) {
    mat <- lst_score[[ method ]]
    data <- as.data.frame(t(mat))
    data$cell <- rownames(data)
    data <- tidyr::pivot_longer(
      data,
      -cell,
      names_to = "Function",
      values_to = "Score"
    )
    data$Method <- method
    data
  })
  data <- dplyr::bind_rows(lst_data)
  data <- dplyr::group_by(data, Method, Function)
  data <- dplyr::mutate(
    data,
    top_cutoff = stats::quantile(Score, probs = 1 - top_cell_prop, na.rm = TRUE),
    is_top_cell = Score >= top_cutoff,
    score_cutoff = score_cutoff,
    pass_score_cutoff = ifelse(Method == "Scoring", Score >= score_cutoff, NA)
  )
  dplyr::ungroup(data)
}


aucellFuns$prepare_group_score_data <- function(data_score_long, metadata,
  group_col = "group", celltype_col = NULL, group_levels = NULL,
  score_methods = "AUCell", score_methods_done = NULL)
{
  if (is.null(metadata) || !is.data.frame(metadata)) {
    return(NULL)
  }
  if (is.null(celltype_col) || length(celltype_col) != 1L || is.na(celltype_col)) {
    return(NULL)
  }
  if (!all(c("cell", group_col, celltype_col) %in% colnames(metadata))) {
    return(NULL)
  }
  if (is.null(data_score_long) || nrow(data_score_long) == 0L) {
    return(NULL)
  }
  if (!all(c("cell", "Method", "Function", "Score") %in% colnames(data_score_long))) {
    return(NULL)
  }

  gene_sets <- unique(as.character(data_score_long$Function))
  gene_sets <- gene_sets[!is.na(gene_sets)]
  if (length(gene_sets) != 1L) {
    warning(
      "Skip cell type-specific Wilcoxon test because the current boxplot helper is restricted to one gene set.",
      call. = FALSE
    )
    return(NULL)
  }

  score_methods_done <- unique(score_methods_done)
  score_methods_done <- score_methods_done[!is.na(score_methods_done)]
  if (length(score_methods_done) == 0L) {
    score_methods_done <- unique(as.character(data_score_long$Method))
  }

  if (length(score_methods) == 1L && identical(score_methods, "all")) {
    score_methods_use <- score_methods_done
  } else {
    score_methods_use <- intersect(score_methods, score_methods_done)
  }
  if (length(score_methods_use) == 0L) {
    return(NULL)
  }

  data_meta <- metadata[, c("cell", group_col, celltype_col), drop = FALSE]
  colnames(data_meta) <- c("cell", "group", "celltype")
  data_meta$cell <- as.character(data_meta$cell)
  data_meta$group <- as.character(data_meta$group)
  data_meta$celltype <- as.character(data_meta$celltype)
  data_meta <- data_meta[
    !is.na(data_meta$cell) &
      !is.na(data_meta$group) &
      !is.na(data_meta$celltype),
    , drop = FALSE
  ]

  group_values <- unique(data_meta$group)
  group_values <- group_values[!is.na(group_values)]
  if (!is.null(group_levels) && length(group_levels) == 2L && all(group_levels %in% group_values)) {
    group_levels_use <- as.character(group_levels)
  } else if (length(group_values) == 2L) {
    group_levels_use <- as.character(group_values)
  } else {
    return(NULL)
  }

  data <- dplyr::filter(data_score_long, Method %in% score_methods_use)
  data <- dplyr::left_join(data, data_meta, by = "cell")
  data <- dplyr::filter(
    data,
    !is.na(group),
    !is.na(celltype),
    group %in% group_levels_use
  )
  if (nrow(data) == 0L) {
    return(NULL)
  }

  data$group <- factor(data$group, levels = group_levels_use)

  data_valid <- lapply(split(data, data$celltype, drop = TRUE), function(dat) {
    n_by_group <- table(dat$group)
    if (!all(group_levels_use %in% names(n_by_group))) {
      return(NULL)
    }
    if (any(n_by_group[group_levels_use] < 1L)) {
      return(NULL)
    }
    dat
  })
  data <- dplyr::bind_rows(data_valid)
  if (is.null(data) || nrow(data) == 0L) {
    return(NULL)
  }

  data$value <- data$Score
  if (length(score_methods_use) == 1L) {
    data$var <- as.character(data$celltype)
  } else {
    data$var <- paste(data$Method, data$celltype, sep = " | ")
  }
  data$var <- factor(data$var, levels = unique(data$var))

  list(
    data = data,
    score_methods = score_methods_use,
    group_levels = group_levels_use,
    celltype_col = celltype_col,
    gene_set = gene_sets[1L]
  )
}

aucellFuns$summarize_score_by_group <- function(lst_score, metadata, group.by)
{
  lapply(lst_score, function(mat) {
    data <- cbind(
      metadata[, group.by, drop = FALSE],
      as.data.frame(t(mat))
    )
    data <- as.data.table(data)
    data <- data.table:::`[.data.table`(
      data, , lapply(.SD, mean), by = group.by,
      .SDcols = setdiff(names(data), group.by)
    )
    tibble::as_tibble(data)
  })
}

aucellFuns$as_score_mean_long <- function(lst_score_mean, group.by)
{
  lst_data <- lapply(names(lst_score_mean), function(method) {
    data <- lst_score_mean[[ method ]]
    if (is.null(data) || nrow(data) == 0L) {
      return(NULL)
    }
    data <- tidyr::pivot_longer(
      data,
      -!!rlang::sym(group.by),
      names_to = "Function",
      values_to = "Score"
    )
    data$Method <- method
    data
  })
  lst_data <- lst_data[!vapply(lst_data, is.null, logical(1L))]
  if (length(lst_data) == 0L) {
    return(tibble::tibble())
  }
  data <- dplyr::bind_rows(lst_data)
  data <- dplyr::group_by(data, Method, Function)
  data <- dplyr::mutate(
    data,
    Activity = aucellFuns$scale_score(Score)
  )
  dplyr::ungroup(data)
}

aucellFuns$scale_score <- function(x)
{
  x <- as.numeric(x)
  if (length(x) == 0L || all(is.na(x))) {
    return(rep(0, length(x)))
  }
  sd_x <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(sd_x) || sd_x == 0) {
    return(rep(0, length(x)))
  }
  as.numeric(scale(x))
}

aucellFuns$scale_score01 <- function(x)
{
  x <- as.numeric(x)
  if (length(x) == 0L || all(is.na(x))) {
    return(rep(0.5, length(x)))
  }
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)
  if (!is.finite(x_min) || !is.finite(x_max) || x_max == x_min) {
    return(rep(0.5, length(x)))
  }
  (x - x_min) / (x_max - x_min)
}

aucellFuns$plot_score_heatmap <- function(data, group.by)
{
  if (is.null(data) || nrow(data) == 0L) {
    stop("`data` is empty.")
  }

  data$CellGroup <- data[[ group.by ]]
  data$Method <- factor(data$Method, levels = unique(data$Method))
  data$Function <- factor(data$Function, levels = unique(data$Function))
  data$CellGroup <- factor(data$CellGroup, levels = rev(unique(data$CellGroup)))

  data <- dplyr::group_by(data, Method, Function)
  data <- dplyr::mutate(
    data,
    MeanScore01 = aucellFuns$scale_score01(Score)
  )
  data <- dplyr::ungroup(data)

  p <- ggplot(
      data,
      aes(x = Method, y = CellGroup)
    ) +
    geom_point(aes(size = MeanScore01, color = Activity), alpha = 0.9) +
    labs(
      x = "Scoring methods",
      y = "Cell groups",
      color = "Scaled mean score",
      size = "Mean score"
    ) +
    theme_bw() +
    scale_color_gradient2(low = "blue", mid = "lightyellow", high = "red") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_line(linewidth = 0.2),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )

  if (length(unique(data$Function)) > 1L) {
    p <- p + facet_wrap(~ Function)
  }

  p
}


setMethod("quantile", signature = c(x = "job_aucell"),
  function(x, cols = NULL, cut = .75,
    gather = c("merge", "intersect", "respective"), 
    name = x$name, group.by = x$group.by, ...)
  {
    gather <- match.arg(gather)
    if (is.null(cols)) {
      auc <- .get_auc_from_job_aucell(x)
      cols <- colnames(auc)
    }
    data <- x$res_aucell_mean
    if (is.null(data)) {
      stop('is.null(data), no `x$res_aucell_mean` data.')
    }
    celltypes <- quantile(
      data.frame(data), cols, get = group.by, cut = cut, ...
    )
    if (gather == "intersect") {
      if (length(cols) > 1) {
        celltypes <- ins(lst = celltypes)
      } else {
        celltypes <- unlist(celltypes)
      }
      snap <- glue::glue("{bind(cols)} 基因集的平均 AUCell 活性为 Top {(1 - cut) * 100}% 的细胞类型")
      as_feature(as.character(celltypes), snap, nature = "cell")
    } else if (gather == "respective") {
      snap <- glue::glue("基因集的平均 AUCell 活性为 Top {(1 - cut) * 100}% 的细胞类型")
      as_feature(
        lapply(celltypes, as.character), snap, nature = "cell"
      )
    } else if (gather == "merge") {
      snap <- glue::glue("{name} 基因集的平均 AUCell 活性为 Top {(1 - cut) * 100}% 合并后的细胞类型")
      celltypes <- unique(as.character(unlist(celltypes)))
      as_feature(setNames(list(celltypes), name), snap, nature = "cell")
    }
  })

setMethod("step2", signature = c(x = "job_aucell"),
  function(x, group.by = "seurat_clusters")
  {
    step_message("Annotate for clusters.")
    auc <- .get_auc_from_job_aucell(x)
    dt <- as.data.table(auc)
    dt$cluster <- x$metadata[[ group.by ]]
    if (is.null(dt$cluster)) {
      stop('is.null(dt$cluster), no value in of column `group.by` in metadata.')
    }
    cluster_mean <- dt[ ,
      lapply(.SD, mean),
      by = cluster,
      .SDcols = setdiff(names(dt), "cluster")
      ]
    mean_scaled <- t(scale(t(cluster_mean[, -1])))
    x <- methodAdd(x, "依据细胞分群结果 ({group.by}) 对同一亚群内细胞的 AUC 分数取平均，以获得各细胞亚群的整体功能活性特征。为消除不同基因集之间评分尺度差异，对每个亚群的平均 AUC 矩阵按 Cluster 进行 Z-score 标准化处理。")
    mean_scaled <- tibble::as_tibble(
      mean_scaled
    )
    mean_scaled <- dplyr::mutate(
      mean_scaled, cluster = cluster_mean$cluster, 
      .before = 1
    )
    mean_scaled <- tidyr::pivot_longer(
      mean_scaled, -cluster,
      names_to = "Function", values_to = "Activity"
    )
    annotation <- dplyr::group_by(mean_scaled, cluster)
    annotation <- dplyr::summarise(
      annotation, Function = Function[ which.max(Activity) ]
    )
    x <- methodAdd(x, "对于每个细胞亚群，进一步筛选其标准化活性值最高的功能通路，并将该通路作为该亚群的主要功能注释（Annotation）。")
    data <- dplyr::filter(
      mean_scaled, Function %in% !!annotation$Function
    )
    data <- map(
      data, "cluster", annotation, "cluster", "Function", col = "Annotation"
    )
    args <- list(
      .data = data, .row = quote(Function), .column = quote(cluster),
      .value = quote(Activity), group_by = quote(Annotation),
      cluster_columns = TRUE, column_names_rot = 45,
      cluster_rows = TRUE,
      row_names_max_width = grobWidth(textGrob(data$Function, gpar(fontsize = 10, fontface = 1)))
    )
    rm(dt, auc)
    p.hp <- wrap_scale_heatmap(
      funPlot(heatmap_with_group, args),
      data$cluster, data$Function, pre_width = 6
    )
    p.hp <- set_lab_legend(
      p.hp,
      glue::glue("{x@sig} Cluster functional enrichment score heatmap"),
      glue::glue("功能富集得分热图|||热图展示不同细胞亚群与代表性功能状态之间的对应关系，其中横坐标细胞亚群，纵坐标表示表示功能基因集（Function），颜色梯度表示标准化后的相对活性强弱（Activity），暖色代表该亚群中该功能相对激活，冷色代表相对低活性。热图的 Cluster 对应有注释类型。热图行列均基于功能活性模式进行层次聚类。")
    )
    types <- unique(annotation$Function)
    x <- snapAdd(
      x, "如图{aref(p.hp)} (热图仅展示有 cluster 注释的功能的活性)，⟦mark$red('AUCell 亚群功能富集一共注释了 {length(types)} 种类型的功能，分别为：{bind(types)}')⟧。"
    )
    x$metadata <- map(
      x$metadata, group.by, annotation, "cluster", "Function",
      col = "AUCell_Function"
    )
    x <- plotsAdd(x, p.hp)
    return(x)
  })

setMethod("step3", signature = c(x = "job_aucell"),
  function(x, use.trait, data_trait = x$metadata, rerun = FALSE)
  {
    step_message("Correlation with trait data.")
    auc <- .get_auc_from_job_aucell(x)
    types <- unique(x$metadata[[ "AUCell_Function" ]])
    data_trait <- data_trait[, use.trait, drop = FALSE]
    data_aucell <- as.data.frame(auc[, colnames(auc) %in% types])
    fun_cor <- function(...) {
      cli::cli_alert_info("safe_fortify_cor")
      safe_fortify_cor(data_trait, data_aucell)
    }
    x$cor_trait_aucell <- expect_local_data(
      "tmp", "aucell_trait_activity_cor", fun_cor, rerun = rerun,
      list(
        x$metadata$cell, colnames(data_aucell), 
        colnames(data_trait)
      )
    )
    snap_cor <- .stat_ggcor_table_list(
      x$cor_trait_aucell, "Trait", "Function"
    )
    p.cor_trait_aucell <- .ggcor_add_general_style(ggcor::quickcor(x$cor_trait_aucell))
    p.cor_trait_aucell <- set_lab_legend(
      wrap_scale_heatmap(p.cor_trait_aucell, length(use.trait), length(types), raw = FALSE),
      glue::glue("{x@sig} trait correlation with AUCell"),
      glue::glue("{bind(use.trait)} 与 AUCell 活性关联分析热图|||热图中颜色表示相关系数的大小，颜色越深表示相关系数越高。P 值以 * 标注 ({.md_p_significant})。")
    )
    x <- snapAdd(x, "对表型 ({bind(use.trait)}) 与 AUCell Function 活性之间关联分析，如图{aref(p.cor_trait_aucell)}，{snap_cor}")
    x <- plotsAdd(x, p.cor_trait_aucell)
    return(x)
  })

.get_auc_from_job_aucell <- function(x) {
  if (!identical(colnames(x$res_aucell), x$metadata$cell)) {
    stop('!identical(colnames(x$res_aucell), x$metadata$cell).')
  }
  t(x$res_aucell)
}

setMethod("map", signature = c(x = "job_aucell", ref = "job_seurat"),
  function(x, ref, use.trait = NULL, use.function = NULL,
    group.by = ref$group.by, pal = NULL, .name = "seurat")
  {
    fun_show <- function(string) stringr::str_wrap(gs(string, "_", " "), 20)
    fun_rename_title <- function(lst) {
      lst$title <- fun_show(lst$title)
      lst
    }
    if (x@step == 1L) {
      if (is.null(use.function)) {
        auc <- .get_auc_from_job_aucell(x)
        use.function <- colnames(auc)
      }
      meta <- dplyr::select(x$metadata, cell)
      meta <- cbind(meta, auc)
      meta <- data.frame(meta[, -1L, drop = FALSE], row.names = meta$cell)
      object(ref) <- SeuratObject::AddMetaData(object(ref), meta)
      layout <- wrap_layout(NULL, length(use.function))
      ps.map <- e(Seurat::FeaturePlot(object(ref), 
          features = use.function, combine = FALSE
          ))
      if (!is(ps.map, "list")) {
        ps.map <- list(ps.map)
      }
      ps.map <- lapply(ps.map, 
        function(x) {
          x + theme(plot.title = element_text(face = "plain", size = 10))
        })
      if (!is.null(fun_show)) {
        ps.map <- lapply(
          ps.map, .set_ggplot_content,
          fun = fun_rename_title,
          slot = "labels"
        )
      }
      p.map <- add(layout, ps.map, TRUE)
      p.map <- set_lab_legend(p.map,
        glue::glue("{x@sig} AUCell Activity UMAP mapping"),
        glue::glue("AUCell 功能活性 UMAP 图||| {bind(use.function)} 的 AUCell 功能活性 UMAP 图。")
      )
      x[[ glue::glue("map_{.name}") ]] <- namel(p.map, metadata = meta)
    } else if (x@step > 1L) {
      meta <- x$metadata[, !duplicated(colnames(x$metadata))]
      meta <- dplyr::select(meta, cell, AUCell_Function)
      col_r_trait <- NULL
      col_trait <- NULL
      if (!is.null(use.trait)) {
        if (is.null(x$cor_trait_aucell)) {
          stop('is.null(x$cor_trait_aucell), but !is.null(use.trait)')
        }
        data_cor <- dplyr::filter(
          x$cor_trait_aucell, .row.names %in% use.trait
        )
        lst_cor <- split(data_cor, ~ .row.names)
        col_r_trait <- glue::glue("r_{use.trait}")
        for (i in seq_along(lst_cor)) {
          meta <- map(
            meta, "AUCell_Function", lst_cor[[i]], ".col.names", 
            "r", col = col_r_trait[ i ]
          )
        }
        allAvai <- c(colnames(meta(ref)), colnames(meta))
        col_trait <- use.trait[ use.trait %in% allAvai ]
        if (any(isNot <- !use.trait %in% col_trait)) {
          warning("Can not got trait from `ref`: ", bind(use.trait[ isNot ]))
        }
      }
      meta <- data.frame(meta[, -1L, drop = FALSE], row.names = meta$cell)
      object(ref) <- SeuratObject::AddMetaData(object(ref), meta)
      group <- c("AUCell_Function", group.by, col_trait)
      ps.map <- e(Seurat::DimPlot(
          object(ref), pt.size = if (dim(object(ref))[2] > 30000) .3 else .5,
          group.by = group,
          cols = color_set(), combine = FALSE
          ))
      if (!is.null(pal) && !is.null(use.trait)) {
        whichTrait <- which(group %in% use.trait)
        for (i in whichTrait) {
          ps.map[[i]] <- ps.map[[i]] + scale_color_manual(values = pal)
        }
      }
      if (!is.null(col_r_trait)) {
        ps2.map <- e(Seurat::FeaturePlot(object(ref), 
            features = col_r_trait, combine = FALSE
            ))
        if (!is(ps2.map, "list")) {
          ps2.map <- list(ps2.map)
        }
        for (i in seq_along(col_r_trait)) {
          ps2.map[[i]] <- ps2.map[[i]] + .scale_for_cor_palette("color")
        }
        group <- c(group, col_r_trait)
        ps.map <- c(ps.map, ps2.map)
      }
      legend_ex <- ""
      if (!is.null(use.trait)) {
        legend_ex <- glue::glue("其中，{bind(col_r_trait)} 对应为 AUCell_Function 与 {bind(use.trait)} 的关联分析的相关系数。")
      }
      layout <- z7(wrap_layout(NULL, length(group), ncol = 2), 1.7, 1)
      p.map <- add(layout, ps.map, TRUE)
      p.map <- set_lab_legend(p.map,
        glue::glue("{x@sig} Cell Function UMAP mapping"),
        glue::glue("AUCell 功能活性 UMAP 图|||依次对应为 {bind(group)} 的 UMAP 图。{legend_ex}")
      )
      x[[ glue::glue("map_{.name}") ]] <- namel(p.map, metadata = meta)
    }
    return(x)
  })

setMethod("clear", signature = c(x = "job_aucell"),
  function(x, save = TRUE, lite = TRUE, suffix = NULL, name = substitute(x, parent.frame(1)))
  {
    eval(name)
    if (save) {
      callNextMethod(
        x, save = save, lite = FALSE, suffix = suffix, name = name
      )
    }
    object(x) <- NULL
    x$res_aucell <- NULL
    if (lite) {
      callNextMethod(
        x, save = FALSE, lite = TRUE, suffix = suffix, name = name
      )
    }
    return(x)
  })

setMethod("map", signature = c(x = "job_seurat", ref = "job_aucell"),
  function(x, ref, type = "AUC", scale = FALSE, run_focus = TRUE)
  {
    if (ref@step < 1L) {
      stop('ref@step < 1L.')
    }
    if (type == "AUC") {
      res <- t(ref$res_aucell)
      colnames(res) <- paste0("AUC_", colnames(res))
    }
    res <- res[match(rownames(object(x)@meta.data), rownames(res)), , drop = FALSE]
    if (scale) {
      res <- scale(res)
    }
    object(x)@meta.data <- object(x)@meta.data[, !colnames(object(x)@meta.data) %in% colnames(res)]
    object(x)@meta.data <- cbind(object(x)@meta.data, res)
    if (ncol(res) <= 20 && run_focus) {
      x <- focus(x, colnames(res), name = "AUCell", cols = c("skyblue", "blue", "black"))
    }
    return(x)
  })

setMethod("set_remote", signature = c(x = "job_aucell"),
  function(x, wd = glue::glue("~/aucell_{x@sig}")){
    x$wd <- wd
    rem_dir.create(wd, wd = ".")
    return(x)
  })

# ==========================================================================

geneFuns <- new.env(parent = baseenv())

geneFuns$check_pkg <- function(pkg)
{
  if (!base::requireNamespace(pkg, quietly = TRUE)) {
    base::stop(glue::glue("Package `{pkg}` is required but not installed."))
  }

  base::invisible(TRUE)
}

geneFuns$set_kegg_pathway_id <- function(pathway_id, species = "hsa")
{
  vec_id <- base::trimws(base::as.character(pathway_id))
  vec_id <- base::sub("^path:", "", vec_id)

  vec_id <- base::ifelse(
    base::grepl("^\\d{5}$", vec_id),
    base::paste0(species, vec_id),
    vec_id
  )

  vec_id <- base::ifelse(
    base::grepl("^map\\d{5}$", vec_id),
    base::paste0(species, base::sub("^map", "", vec_id)),
    vec_id
  )

  vec_id
}

geneFuns$get_cache_file <- function(cache_dir, key)
{
  if (base::is.null(cache_dir)) {
    return(NULL)
  }

  if (!base::dir.exists(cache_dir)) {
    base::dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }

  key <- base::gsub("[^A-Za-z0-9_.-]+", "_", key)
  base::file.path(cache_dir, base::paste0(key, ".rds"))
}

geneFuns$get_kegg_pathway_table <- function(species = "hsa",
  cache_dir = "cache/kegg", refresh = FALSE)
{
  geneFuns$check_pkg("KEGGREST")

  cache_file <- geneFuns$get_cache_file(
    cache_dir,
    base::paste0("kegg_pathway_table_", species)
  )

  if (!refresh && !base::is.null(cache_file) && base::file.exists(cache_file)) {
    return(base::readRDS(cache_file))
  }

  base::message(glue::glue("Download KEGG pathway table: {species}."))
  vec_pathway <- KEGGREST::keggList("pathway", species)

  data_pathway <- data.frame(
    pathway_id = base::sub("^path:", "", base::names(vec_pathway)),
    pathway_name = base::as.character(vec_pathway),
    stringsAsFactors = FALSE
  )

  data_pathway$pathway_name <- base::sub(" - .*$", "", data_pathway$pathway_name)

  if (!base::is.null(cache_file)) {
    base::saveRDS(data_pathway, cache_file)
  }

  data_pathway
}

geneFuns$search_kegg_pathway <- function(query, species = "hsa",
  cache_dir = "cache/kegg", refresh = FALSE, ignore.case = TRUE)
{
  data_pathway <- geneFuns$get_kegg_pathway_table(
    species = species,
    cache_dir = cache_dir,
    refresh = refresh
  )

  vec_keep <- base::grepl(
    query,
    data_pathway$pathway_name,
    ignore.case = ignore.case
  )

  data_pathway[vec_keep, , drop = FALSE]
}

geneFuns$get_kegg_pathway_entry <- function(pathway_id, species = "hsa",
  cache_dir = "cache/kegg", refresh = FALSE)
{
  geneFuns$check_pkg("KEGGREST")

  pathway_id <- geneFuns$set_kegg_pathway_id(pathway_id, species = species)

  cache_file <- geneFuns$get_cache_file(
    cache_dir,
    base::paste0("kegg_pathway_entry_", pathway_id)
  )

  if (!refresh && !base::is.null(cache_file) && base::file.exists(cache_file)) {
    return(base::readRDS(cache_file))
  }

  base::message(glue::glue("Download KEGG pathway entry: {pathway_id}."))
  lst_entry <- KEGGREST::keggGet(pathway_id)

  if (base::length(lst_entry) == 0L) {
    base::stop(glue::glue("No KEGG entry was returned: {pathway_id}."))
  }

  entry <- lst_entry[[ 1L ]]

  if (!base::is.null(cache_file)) {
    base::saveRDS(entry, cache_file)
  }

  entry
}

geneFuns$get_gene_table_from_kegg_entry <- function(entry, pathway_id)
{
  if (base::is.null(entry$GENE) || base::length(entry$GENE) == 0L) {
    return(data.frame(
      pathway_id = character(),
      pathway_name = character(),
      entrez_id = character(),
      symbol = character(),
      description = character(),
      stringsAsFactors = FALSE
    ))
  }

  vec_gene <- entry$GENE
  vec_desc <- base::as.character(vec_gene)
  vec_entrez <- base::names(vec_gene)

  if (base::is.null(vec_entrez) || base::length(vec_entrez) == 0L ||
      base::all(base::is.na(vec_entrez)) || base::all(vec_entrez == "")) {
    vec_entrez <- vec_desc[base::seq(1L, base::length(vec_desc), by = 2L)]
    vec_desc <- vec_desc[base::seq(2L, base::length(vec_desc), by = 2L)]
  }

  vec_entrez <- base::sub("^.+:", "", base::as.character(vec_entrez))
  vec_symbol <- base::trimws(base::sub(",.*$", "", base::sub(";.*$", "", vec_desc)))

  pathway_name <- NA_character_

  if (!base::is.null(entry$NAME) && base::length(entry$NAME) > 0L) {
    pathway_name <- base::paste(base::as.character(entry$NAME), collapse = " ")
    pathway_name <- base::sub(" - .*$", "", pathway_name)

    if (base::length(pathway_name) == 0L || base::is.na(pathway_name) ||
        pathway_name == "") {
      pathway_name <- NA_character_
    }
  }

  data_gene <- data.frame(
    pathway_id = base::rep(pathway_id, base::length(vec_entrez)),
    pathway_name = base::rep(pathway_name, base::length(vec_entrez)),
    entrez_id = vec_entrez,
    symbol = vec_symbol,
    description = vec_desc,
    stringsAsFactors = FALSE
  )

  data_gene <- data_gene[!base::is.na(data_gene$entrez_id), , drop = FALSE]
  data_gene <- data_gene[data_gene$entrez_id != "", , drop = FALSE]
  data_gene <- data_gene[!base::duplicated(data_gene$entrez_id), , drop = FALSE]
  row.names(data_gene) <- NULL

  data_gene
}


geneFuns$get_kegg_pathway_gene_table <- function(pathway_ids,
  species = "hsa", cache_dir = "cache/kegg", refresh = FALSE)
{
  vec_pathway_id <- geneFuns$set_kegg_pathway_id(
    pathway_ids,
    species = species
  )

  lst_gene <- base::lapply(vec_pathway_id, function(pathway_id) {
    entry <- geneFuns$get_kegg_pathway_entry(
      pathway_id = pathway_id,
      species = species,
      cache_dir = cache_dir,
      refresh = refresh
    )

    geneFuns$get_gene_table_from_kegg_entry(
      entry = entry,
      pathway_id = pathway_id
    )
  })

  data_gene <- do.call(rbind, lst_gene)

  if (base::is.null(data_gene) || base::nrow(data_gene) == 0L) {
    base::warning("No genes were obtained from KEGG pathways.")
    return(data.frame(
      pathway_id = character(),
      pathway_name = character(),
      entrez_id = character(),
      symbol = character(),
      description = character(),
      stringsAsFactors = FALSE
    ))
  }

  data_gene <- data_gene[!base::duplicated(
    base::paste(data_gene$pathway_id, data_gene$entrez_id, sep = "|||")
  ), , drop = FALSE]

  row.names(data_gene) <- NULL
  data_gene
}

geneFuns$get_kegg_pathway_gene_list <- function(pathway_ids,
  species = "hsa", gene_id_type = c("symbol", "entrez"),
  cache_dir = "cache/kegg", refresh = FALSE, min_genes = 5L)
{
  gene_id_type <- base::match.arg(gene_id_type)

  data_gene <- geneFuns$get_kegg_pathway_gene_table(
    pathway_ids = pathway_ids,
    species = species,
    cache_dir = cache_dir,
    refresh = refresh
  )

  if (base::nrow(data_gene) == 0L) {
    return(list())
  }

  id_col <- switch(
    gene_id_type,
    symbol = "symbol",
    entrez = "entrez_id"
  )

  data_gene$gene_id <- base::as.character(data_gene[[ id_col ]])
  data_gene <- data_gene[!base::is.na(data_gene$gene_id), , drop = FALSE]
  data_gene <- data_gene[data_gene$gene_id != "", , drop = FALSE]

  data_gene$gene_set_name <- base::ifelse(
    !base::is.na(data_gene$pathway_name) & data_gene$pathway_name != "",
    base::paste0(data_gene$pathway_id, " | ", data_gene$pathway_name),
    data_gene$pathway_id
  )

  lst_gene <- base::split(data_gene$gene_id, data_gene$gene_set_name)
  lst_gene <- base::lapply(lst_gene, unique)

  vec_keep <- base::vapply(
    lst_gene,
    function(vec_gene) base::length(vec_gene) >= min_genes,
    logical(1L)
  )

  lst_gene[vec_keep]
}

