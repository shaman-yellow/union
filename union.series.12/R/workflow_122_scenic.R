# ==========================================================================
# workflow of scenic
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_scenic <- setClass("job_scenic", 
  contains = c("job"),
  prototype = prototype(
    pg = "scenic",
    info = c(
      "https://pyscenic.readthedocs.io/en/latest/",
      "https://arboreto.readthedocs.io/en/latest/",
      "https://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Setup.html",
      "https://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Running.html",
      "https://github.com/aertslab/SCENICprotocol/blob/master/notebooks/PBMC10k_SCENIC-protocol-CLI.ipynb"
      ),
    cite = "",
    method = "",
    tag = "scenic",
    analysis = "SCENIC 转录因子分析"
    ))

setGeneric("asjob_scenic",
  function(x, ...) standardGeneric("asjob_scenic"))

setMethod("step0", signature = c(x = "job_scenic"),
  function(x){
    step_message("Prepare your data with function `job_scenic`.")
  })

setMethod("asjob_scenic", signature = c(x = "job_seurat"),
  function(x, sketch = TRUE, ncells = 10000, filter = TRUE,
    org = c("human"), dir_db = pg("db_scenic"), overwrite = FALSE)
  {
    dir_cache <- create_job_cache_dir(x, "scene")
    args <- list(
      sketch = sketch, ncells = ncells, filter = filter, org = org
    )
    hash <- digest::digest(args, "xxhash64", serializeVersion = 3)
    file_expr <- file.path(dir_cache, glue::glue("expr_{hash}.loom"))
    methodAdd_onExit(
      "x", "为在复杂细胞群体中识别关键转录因子及其调控网络，为后续功能解析提供依据，本研究采用 pySCENIC (0.12.1, PMID: 32561888) 进行转录调控分析。"
    )
    methodAdd_onExit("x", "该方法以基因表达矩阵为输入，首先基于共表达关系（如 GRNBoost2）推断转录因子与靶基因之间的调控网络，其次结合 motif 注释对网络进行剪枝以保留高置信调控关系，最后通过 AUCell 在单细胞水平量化各调控模块（regulon）的活性，从而识别关键转录因子及其调控特征。")
    if (sketch) {
      if (ncol(object(x)) > ncells * 1.5) {
        ncells_pre <- ncol(object(x))
        message(glue::glue("Use `Seurat::SketchData` to sampling cells."))
        Seurat::DefaultAssay(object(x)) <- "RNA"
        if (!length(Seurat::VariableFeatures(object(x)))) {
          object(x) <- e(Seurat::NormalizeData(object(x)))
          object(x) <- e(Seurat::FindVariableFeatures(object(x)))
        }
        object(x) <- e(Seurat::SketchData(
            object(x), ncells = ncells, seed = x$seed
            ))
        counts <- SeuratObject::GetAssayData(
          object(x), assay = "sketch", layer = "counts"
        )
        methodAdd_onExit("x", "在为 GRNBoost2 提供构建调控网络的表达数据之前，为保证计算效率的同时尽可能保留数据的整体结构信息，本研究采用了基于 “sketch” 的代表性抽样策略 (`Seurat::SketchData`，使用的算法为 LeverageScore)。⟦mark$blue('在抽样前，数据集包含 {ncells_pre} 个细胞，抽样后，子集包含 { ncells } 个细胞')⟧。")
      }
    } else {
      counts <- SeuratObject::GetAssayData(
        object(x), assay = "RNA", layer = "counts"
      )
    }
    if (filter) {
      ngenes_pre <- nrow(counts)
      keep_genes <- Matrix::rowSums(counts > 0) > (0.01 * ncol(counts))
      counts <- counts[keep_genes, ]
      methodAdd_onExit("x", "对基因进行预过滤以降低数据稀疏性并减少计算负担，具体而言，⟦mark$blue('仅保留在至少 1% 细胞中被检测到（表达值大于 0）的基因。在过滤前，数据集包含 {ngenes_pre} 个基因，过滤后，数据集包含 {sum(keep_genes)} 个基因')⟧。")
    }
    fun_save <- function(data, file) {
      if (file.exists(file) && overwrite) {
        file.remove(file)
      }
      if (!file.exists(file)) {
        loom <- SCopeLoomR::build_loom(file, data)
        loom$close_all()
      }
    }
    fun_save(counts, file_expr)
    used_genes <- rownames(counts)
    file_expr_all <- file_expr
    if (sketch) {
      file_expr_all <- file.path(dir_cache, glue::glue("expr_all_{hash}.loom"))
      counts_alls <- SeuratObject::GetAssayData(
        object(x), assay = "RNA", layer = "counts"
      )
      counts_alls <- counts_alls[ rownames(counts_alls) %in% used_genes, ]
      message(glue::glue("Set sketch. Save all cells for AUCells running."))
      fun_save(counts_alls, file_expr_all)
    }
    pr <- params(x)
    x <- .job_scenic()
    x$used_genes <- used_genes
    x$dir_cache <- dir_cache
    x$hash <- hash
    x$file_expr <- file_expr
    x$file_expr_all <- file_expr_all
    x$sketch <- sketch
    x@params <- append(x@params, pr)
    if (org == "human") {
      pattern <- "hg38"
      pattern_mutate <- "hgnc"
    }
    x$org <- org
    if (!dir.exists(dir_db)) {
      dir.create(dir_db, FALSE)
    }
    fun_get <- function(pattern, fun_download) {
      fun_list <- function() list.files(dir_db, pattern, recursive = TRUE, full.names = TRUE)
      files <- fun_list()
      if (org == "human" && !length(files) && sureThat("No any {pattern} files, download it?")) {
        fun_download()
        files <- fun_list()
      }
      files
    }
    x$file_feathers <- fun_get(
      glue::glue("{pattern}.*\\.rankings.feather$"), .download_feather_files
    )
    x$file_motifs <- fun_get(
      glue::glue("{pattern_mutate}.*\\.tbl$"), .download_motif2tf_files
    )
    x$file_tf <- fun_get(
      glue::glue("{pattern}\\.txt$"), .download_tf_files
    )
    lapply(c("file_tf", "file_feathers", "file_motifs"),
      function(name) {
        if (!length(x[[name]])) {
          stop(glue::glue("Not found any file for {name}?"))
        }
      })
    return(x)
  })

setMethod("step1", signature = c(x = "job_scenic"),
  function(x, workers = 16)
  {
    step_message("Run pyscenic grn")
    # cmd <- "pyscenic grn {f_loom_path_scenic} {f_tfs} -o adj.csv --num_workers {workers}"
    x$file_log <- "task.log"
    x$file_adj <- file.path(x$dir_cache, glue::glue("adj_{x$hash}.tsv"))
    cmd <- glue::glue(
      "{pg('pyscenic')} grn {x$file_expr} {x$file_tf} -o {x$file_adj} --num_workers {workers} --seed {x$seed}"
    )
    .expect_file_with_cmd_run(
      cmd, NULL, x$file_adj, x$file_log, check = FALSE
    )
    return(x)
  })

setMethod("step2", signature = c(x = "job_scenic"),
  function(x, workers = 3L)
  {
    step_message("Run pyscenic ctx")
    x$file_reg <- file.path(x$dir_cache, glue::glue("reg_{x$hash}.tsv"))
    cmd <- .glue_cmd(
      "{pg('pyscenic')} ctx {x$file_adj} {bind(x$file_feathers, co = ' ')} ",
      "--annotations_fname {x$file_motifs} --expression_mtx_fname {x$file_expr} ",
      "--output {x$file_reg} --mask_dropouts --num_workers {workers}"
    )
    .expect_file_with_cmd_run(cmd, x$file_adj, x$file_reg, x$file_log)
    return(x)
  })

setMethod("step3", signature = c(x = "job_scenic"),
  function(x, workers = 3L){
    step_message("Run pyscenic aucell")
    x$file_aucell <- file.path(x$dir_cache, glue::glue("aucell_{x$hash}.loom"))
    if (!is.null(x$sketch) && x$sketch) {
      x <- methodAdd(x, "在 GRNBoost2，cisTarget 的算法部分，使用了上述子集细胞进行运算。而在 AUCell 部分，我们将使用全部细胞来计算 AUCell 活性评分，从而得到所有细胞的 AUCell 活性。")
    }
    cmd <- .glue_cmd(
      "{pg('pyscenic')} aucell {x$file_expr_all} {x$file_reg} ",
      "--output {x$file_aucell} --num_workers {workers} --seed {x$seed} "
    )
    .expect_file_with_cmd_run(cmd, x$file_reg, x$file_aucell, x$file_log)
    return(x)
  })

setMethod("map", signature = c(x = "job_scenic", ref = "job_seurat"),
  function(x, ref, name = NULL, name.by = NULL, ntop = 3,
    group.by = ref$group.by, .name = "seurat")
  {
    step_message("Mapping results into seurat (heatmap)")
    init(snap(x)) <- TRUE
    init(meth(x)) <- TRUE
    name_save <- glue::glue("map_{.name}")
    lst_results <- list()
    if (!is.null(name)) {
      if (is.null(name.by)) {
        stop('is.null(name.by), no default.')
      }
      message(
        glue::glue("Get sub-object of cells {bind(name)}, within: {name.by}")
      )
      ref <- getsub(ref, !!rlang::sym(name.by) %in% name)
      x <- snapAdd(x, "聚焦于 {bind(name)} 的 regulon 活性。", step = .name)
    }
    object <- .get_seurat_with_regulon(x, ref)
    SeuratObject::Idents(object) <- group.by
    object <- e(Seurat::ScaleData(object))
    p.hp <- e(Seurat::DoHeatmap(
      object, rownames(object), group.by = x$group.by
    ))
    data <- tibble::as_tibble(.get_ggplot_content(p.hp[[1]]))
    # remove void columns generated by `DoHeatmap`
    data <- dplyr::filter(data, !is.na(Expression))
    lst_results$data <- data
    data <- dplyr::mutate(
      data, Cell = paste0(Identity, "_", Cell)
    )
    args <- list(
      .data = data, .row = quote(Feature), .column = quote(Cell),
      .value = quote(Expression), group_by = quote(Identity),
      split_by = quote(Identity),
      cluster_columns = TRUE, show_column_names = FALSE,
      cluster_rows = TRUE, show_row_names = FALSE,
      palette_value = c("#053061FF", "White", "#67001FFF")
    )
    p.hp <- wrap_scale_heatmap(
      funPlot(heatmap_with_group, args),
      data$Cell, data$Feature, pre_width = 3, w.size = .002, 
      h.size = .005, min_width = 12, min_height = 8
    )
    p.hp <- set_lab_legend(
      p.hp,
      glue::glue("{x@sig} AUCell all regulon heatmap"),
      glue::glue("不同细胞亚群的所有 regulon 活性谱热图|||横坐标表示单细胞，按照细胞亚群进行排序；顶部颜色注释条表示对应的细胞分组（Identity）。纵坐标表示 regulon。颜色梯度表示 regulon 活性水平（Expression），其中暖色表示相对高活性，冷色表示相对低活性。")
    )
    lst_results$p.hp <- p.hp
    if (!is.null(ntop)) {
      features <- .get_representative_regulon(data, ntop)
      args$.data <- dplyr::filter(data, Feature %in% !!features)
      args$show_row_names <- TRUE
      p.hp_top <- wrap(
        funPlot(heatmap_with_group, args), 12, 8
      )
      p.hp_top <- set_lab_legend(
        p.hp_top,
        glue::glue("{x@sig} AUCell top regulon heatmap"),
        glue::glue("不同细胞亚群的 Top regulon 活性谱热图|||横坐标表示单细胞，按照细胞亚群进行排序；顶部颜色注释条表示对应的细胞分组（Identity）。纵坐标表示筛选获得的代表性 regulon。颜色梯度表示 regulon 活性水平（Expression），其中暖色表示相对高活性，冷色表示相对低活性。")
      )
      lst_results$p.hp_top <- p.hp_top
      x <- snapAdd(
        x, "为解析不同细胞亚群间潜在转录调控网络的异质性，首先在各细胞群（Identity）内计算每个 regulon 的平均活性水平（Expression），以表征对应转录因子调控程序在该亚群中的整体活化状态。随后，分别筛选各细胞群中活性最高的前 {ntop} 个 regulon 与活性最低的前 {ntop} 个 regulon，前者代表该亚群特异性增强或持续激活的关键调控网络，后者则提示相对沉默或受抑制的调控程序。最终，将所有细胞亚群筛选得到的 regulon 合并并去重，形成具有生物学代表性的候选调控集合，用于热图展示及亚群功能状态比较{aref(p.hp_top)}，从而识别驱动细胞命运转变、功能分化及微环境适应的核心转录调控模式。", step = .name
      )
    }
    x[[ name_save ]] <- lst_results
    return(x)
  })

.get_representative_regulon <- function(df, n = 5L) {

  # Summarise mean expression of each feature in each identity
  stat <- dplyr::group_by(
    df,
    Identity,
    Feature
  )

  stat <- dplyr::summarise(
    stat,
    Expression = mean(Expression, na.rm = TRUE),
    .groups = "drop"
  )

  # Select top n and bottom n features per identity
  stat2 <- dplyr::group_by(
    stat,
    Identity
  )

  stat2 <- dplyr::group_modify(
    stat2,
    function(.x, .y) {

      top_df <- dplyr::slice_max(
        .x,
        order_by = Expression,
        n = n,
        with_ties = FALSE
      )

      bottom_df <- dplyr::slice_min(
        .x,
        order_by = Expression,
        n = n,
        with_ties = FALSE
      )

      dplyr::bind_rows(top_df, bottom_df)
    }
  )

  # Return feature list
  unique(stat2$Feature)
}


setMethod("diff", signature = c(x = "job_scenic"),
  function(x, ref, name = NULL, split.by = ref$group.by, compare.by = "group",
    ident.1 = x$levels[1], ident.2 = x$levels[2],
    use.p = c("p_val_adj", "p_val"), 
    cut.fc = .5, maxShow = 10, .name = "seurat")
  {
    step_message("Differential analysis")
    init(snap(x)) <- TRUE
    init(meth(x)) <- TRUE
    name_save <- glue::glue("diff_{.name}")
    lst_results <- list()
    if (!is(ref, "job_seurat")) {
      stop('!is(ref, "job_seurat").')
    }
    if (!is.null(name)) {
      split.by <- NULL
    }
    use.p <- match.arg(use.p)
    x <- copy_job(x)
    object <- .get_seurat_with_regulon(x, ref)
    # find markers
    if (!is.null(split.by)) {
      stop("The code for this module has not been written yet.")
      SeuratObject::Idents(object) <- split.by
      cellTypes <- levels(
        SeuratObject::Idents(object)
      )
      cli::cli_alert_info("Seurat::FindMarkers")
      diff_regulon <- lapply(cellTypes, 
        function(type) {
          set <- SeuratObject:::subset.Seurat(object, !!rlang::sym(split.by) == type)
          SeuratObject::Idents(set) <- compare.by
          diff_regulon <- Seurat::FindMarkers(
            set, ident.1 = ident.1, ident.2 = ident.2, 
            min.pct = .1, min.cells.group = 1
          )
          tibble::as_tibble(markers)
        })
    } else {
      SeuratObject::Idents(object) <- compare.by
      diff_regulon <- Seurat::FindMarkers(
        object, ident.1 = ident.1, ident.2 = ident.2, min.pct = .1
      )
      diff_regulon <- as_tibble(diff_regulon, idcol = "regulon")
      x <- methodAdd(
        x, "将 Regulon 的 AUCell 活性数据作为 Assay 对象导入 Seurat，对 {name} (n = {ncol(object)}) 以 Seurat::FindMarkers 差异分析 ({ident.1} vs {ident.2}) (⟦mark$blue('{detail(use.p)} &lt; 0.05，|Log2FC| > {cut.fc}')⟧)。", step = .name
      )
    }
    diff_regulon <- dplyr::filter(
      diff_regulon, !!rlang::sym(use.p) < .05, abs(avg_log2FC) > !!cut.fc
    )
    diff_regulon <- set_lab_legend(
      diff_regulon,
      glue::glue("{x@sig} {name} diff Regulon data"),
      glue::glue("{name}中的差异 Regulon")
    )
    # x <- tablesAdd(x, diff_regulon)
    lst_results$diff_regulon <- diff_regulon
    if (nrow(diff_regulon) < 50) {
      object@meta.data <- dplyr::mutate(
        object@meta.data, .group_id = paste0(
          !!rlang::sym(compare.by), "__", !!rlang::sym(ref$group.by)
        )
      )
      dotData <- e(Seurat::DotPlot(
          object, features = diff_regulon$regulon,
          group.by = ".group_id"))
      dotData <- .get_ggplot_content(dotData)
      # features.plot, id
      dotData <- tidyr::separate(
        tibble::as_tibble(dotData), 
        id, into = c("Group", "Cell"), sep = "__", remove = FALSE
      )
      dotData <- dplyr::rename(
        dotData, Regulon = features.plot, Activity_Scaled = avg.exp.scaled,
        Cell_type = id
      )
      p.hp <- funPlot(heatmap_with_group,
        list(
          .data = dotData, .row = quote(Regulon), .column = quote(Cell_type),
          .value = quote(Activity_Scaled), cluster_columns = FALSE,
          cluster_rows = TRUE, group_by = quote(Group)
          ))
      p.hp <- set_lab_legend(
        p.hp,
        glue::glue("{x@sig} Diff regulon in {name}"),
        glue::glue("{name}中的差异 Regulon|||纵坐标表示 Regulon，横坐标表示细胞类型，热图颜色代表标准化后的 Regulon 活性。")
      )
      # x <- plotsAdd(x, p.hp)
      lst_results$p.hp <- p.hp
    } else {
      stop("Too many significant regulon for draw dotplot.")
    }
    levels <- c(ident.1, ident.2)
    data <- head(diff_regulon, n = maxShow)
    cp <- .setup_compare_pvalue_with_table(
      data, "regulon", use.p, "avg_log2FC",
      levels = levels
    )
    snap <- .stat_compare_by_pvalue(cp, levels, "", mode = "active")
    x <- snapAdd(
      x, "对 Regulon 的 AUCell 活性差异分析一共筛选到 {nrow(diff_regulon)} 个显著差异的 Regulon。对显著性结果按 Log2FC 降序排序，排名前 {nrow(data)} 的 Regulon 中，{snap}", step = .name
    )
    x[[ name_save ]] <- lst_results
    x <- stepPostModify(x, formal = FALSE)
    return(x)
  })

.get_seurat_with_regulon <- function(x, ref) {
  # get data
  loom <- SCopeLoomR::open_loom(x$file_aucell)
  meta_regulons <- SCopeLoomR::get_regulons(loom, "Regulons")
  regulons_AUC <- SCopeLoomR::get_regulons_AUC(loom, "RegulonsAUC")
  SCopeLoomR::close_loom(loom)
  # format
  meta_regulons <- apply(
    meta_regulons, 1, function(x) names(x)[x], simplify = FALSE
  )
  matrix_auc <- e(AUCell::getAUC(regulons_AUC))
  message(glue::glue("Matrix auc, dim: {bind(dim(matrix_auc))}"))
  # set in seurat
  object <- object(ref)
  message(
    glue::glue("Object Seurat, dim: {bind(dim(object))}")
  )
  matrix_auc <- matrix_auc[, SeuratObject::Cells(object)]
  message(
    glue::glue("Get from matrix auc, dim: {bind(dim(matrix_auc))}")
  )
  object[["regulon"]] <- e(SeuratObject::CreateAssayObject(data = matrix_auc))
  SeuratObject::DefaultAssay(object) <- "regulon"
  return(object)
}

# regulons_thresholds <- SCopeLoomR::get_regulon_thresholds(loom)

.expect_file_with_cmd_run <- function(cmd, file_pre, file_aft, file_log = "task.log",
  check = TRUE, fun_check = NULL)
{
  if (is.null(fun_check)) {
    fun_check <- function(file) {
      length(readLines(file, 10)) == 10
    }
  }
  if (check) {
    if (!file.exists(file_pre) || !fun_check(file_pre)) {
      stop(glue::glue("{file_pre} is invalid, the previous program has not yet output results?"))
    }
  }
  if (!file.exists(file_aft)) {
    cmd_run <- glue::glue("sh -c 'date && nohup {cmd} > {file_log} 2>&1 &'")
    message(crayon::yellow(glue::glue("Running in the background, please wait for the result output to complete before proceeding to the next step.")))
    system(cmd_run, wait = FALSE)
  } else {
    message(glue::glue("file.exists: {file_aft}"))
  }
}

.download_feather_files <- function(path) {
  url <- "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/"
  system(glue::glue("wget -r -np -nH --cut-dirs=7 -R 'index.html*' -P {path} {url} --no-check-certificate"))
}

.download_motif2tf_files <- function(path) {
  url <- "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
  system(glue::glue("wget -P {path} {url} --no-check-certificate"))
}

.download_tf_files <- function(path) {
  url <- "https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt"
  system(glue::glue("wget -P {path} {url} --no-check-certificate"))
}

.glue_cmd <- function(..., envir = parent.frame(1)) {
  strings <- vapply(list(...), FUN.VALUE = character(1),
    function(string) {
      glue::glue(string, .envir = envir)
    })
  paste0(strings, collapse = " ")
}

setMethod("status", signature = c(x = "job_scenic"),
  function(x){
    system(glue::glue("tail -f {x$file_log}"))
  })

