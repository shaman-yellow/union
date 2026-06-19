# ==========================================================================
# workflow of infercnv
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_infercnv <- setClass("job_infercnv",
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("https://bioconductor.org/packages/release/bioc/html/infercnv.html"),
    method = "Package inferCNV used for CNV anlysis and cancer cell prediction",
    tag = "scrna:cancer",
    analysis = "InferCNV 变异拷贝数分析"
    ))

setGeneric("asjob_infercnv",
  function(x, ...) standardGeneric("asjob_infercnv"))

setMethod("asjob_infercnv", signature = c(x = "job_seurat"),
  function(x, ref, groups = NULL, ..., recluster = TRUE, subset = "seurat_clusters",
    group.by = x$group.by, outdir = create_job_cache_dir(x, "infercnv"),
    max_cells_per_group = NULL, seed = 123L
  )
  {
    if (missing(ref)) {
      stop('missing(ref).')
    }
    if (!is.null(seed)) {
      message(glue::glue("Set seed as: {seed}"))
      set.seed(seed)
    }
    if (!is.null(max_cells_per_group)) {
      if (!is.numeric(max_cells_per_group) || length(max_cells_per_group) != 1L ||
        is.na(max_cells_per_group) || max_cells_per_group <= 0) {
        stop('`max_cells_per_group` should be a positive numeric scalar or NULL.')
      }
      max_cells_per_group <- as.integer(max_cells_per_group)
    }
    metadata <- dplyr::select(
      as_tibble(object(x)@meta.data), rownames, !!rlang::sym(group.by), !!rlang::sym(subset)
    )
    if (any(!ref %in% metadata[[ group.by ]])) {
      stop('any(!ref %in% metadata[[ group.by ]]).')
    }
    if (is.null(groups)) {
      groups <- unique(as.character(metadata[[group.by]]))
      groups <- groups[ !groups %in% ref ]
    }
    allGroups <- unique(c(groups, ref))
    metadata <- dplyr::filter(metadata, !!rlang::sym(group.by) %in% allGroups)
    dir.create(outdir, FALSE)
    if (recluster && subset == "seurat_clusters" && !is.null(groups)) {
      hash <- digest::digest(
        list(x@sig, rownames(object(x)), groups)
      )
      file_submeta <- add_filename_suffix(
        file.path(outdir, "submetadata.rds"), hash
      )
      if (!file.exists(file_submeta)) {
        sub <- asjob_seurat_sub(x, !!rlang::sym(group.by) %in% groups)
        sub <- step1(sub)
        sub <- step2(sub)
        sub <- step3(sub, ...)
        metadata_sub <- as_tibble(object(sub)@meta.data)
        snap <- snap(metadata_sub) <- snap(sub)
        saveRDS(metadata_sub, file_submeta)
      } else {
        metadata_sub <- readRDS(file_submeta)
        snap <- snap(metadata_sub)
      }
      metadata <- map(
        metadata, "rownames", metadata_sub, "rownames", subset, col = subset
      )
    }
    if (!is.null(subset)) {
      metadata <- dplyr::mutate(
        metadata, group_subset = paste0(
          !!rlang::sym(group.by), "_", !!rlang::sym(subset)
          ),
        group_subset = ifelse(
          !!rlang::sym(group.by) %in% !!ref, 
          as.character(!!rlang::sym(group.by)), group_subset
        )
      )
      message(glue::glue("\n{showStrings(metadata$group_subset, trunc = FALSE)}"))
      metadata <- dplyr::select(metadata, rownames, group_subset)
      group.by <- "group_subset"
    }
    counts <- e(SeuratObject::LayerData(object(x), "count"))
    rownames(counts) <- gname(rownames(counts))
    counts <- counts[ !duplicated(rownames(counts)), ]
    message(glue::glue("Before Cells filter: {bind(dim(counts))}"))
    counts <- counts[, colnames(counts) %in% metadata$rownames]
    message(glue::glue("After Cells filter: {bind(dim(counts))}"))
    ranges <- get_gene_ranges(rownames(counts))
    counts <- counts[rownames(counts) %in% ranges$symbols, ]
    ranges <- ranges[match(rownames(counts), ranges$symbols), ]
    genes <- dplyr::select(
      ranges, symbols, seqnames, start, end
    )
    tmp.metadata <- file.path(outdir, "metadata.tsv")
    write_tsv(metadata, tmp.metadata, col.names = FALSE)
    tmp.genes <- file.path(outdir, "genes.tsv")
    write_tsv(genes, tmp.genes, col.names = FALSE)
    cell_groups <- unique(as.character(metadata[[ group.by ]]))
    if (!requireNamespace("infercnv")) {
      setup_jags_from_libpaths()
    }
    fun_create <- function(...) {
      e(infercnv::CreateInfercnvObject(
          raw_counts_matrix = counts,
          gene_order_file = tmp.genes,
          annotations_file = tmp.metadata,
          ref_group_names = ref,
          max_cells_per_group = max_cells_per_group
          ))
    }
    params <- params(x)
    obj.cnv <- expect_local_data(
      "tmp", "infercnv_create", fun_create,
      list(
        sig(x), params$metadata$cell, metadata$rownames,
        metadata[[group.by]], ref, groups, subset, max_cells_per_group, seed
      )
    )
    x <- .job_infercnv(object = obj.cnv)
    x@params <- append(x@params, params)
    x$outdir <- outdir
    x$tmp.metadata <- tmp.metadata
    x$tmp.genes <- tmp.genes
    x$meta_used <- metadata
    x$infercnv_group_col <- group.by
    x$infercnv_subset <- subset
    x$infercnv_ref <- ref
    x$infercnv_groups <- groups
    x$max_cells_per_group <- max_cells_per_group
    x$infercnv_seed <- seed
    x <- methodAdd(x, "采用 `inferCNV` ⟦pkgInfo('infercnv')⟧ 基于单细胞转录组表达谱推断染色体尺度拷贝数变异，并以预设正常参考细胞作为基线，评估候选肿瘤细胞群的 CNV 信号。")
    if (!is.null(max_cells_per_group)) {
      x <- methodAdd(
        x,
        glue::glue(
          "考虑到候选肿瘤细胞数量较大，为降低大规模 inferCNV 计算负担并避免单一高丰度细胞组过度主导热图展示，",
          "分析在构建 inferCNV 对象时对每个注释组最多随机纳入 {max_cells_per_group} 个细胞；",
          "细胞数未超过该阈值的注释组全部保留。抽样过程固定随机种子，以保证结果可复现。"
        )
      )
    }
    if (exists("snap") && !is.function(snap)) {
      message(glue::glue("Add seurat recluster snap."))
      x <- snapAdd(x, snap)
    }
    x <- snapAdd(x, "使用 `inferCNV` 识别肿瘤细胞的染色体拷贝数变异，选择 {bind(ref)} 为参考细胞 (正常细胞)，识别 {bind(groups)} 中的拷贝数变异。分析中，{bind(groups)}以 {subset} 标识次级聚类。")
    return(x)
  })

setMethod("step0", signature = c(x = "job_infercnv"),
  function(x){
    step_message("Prepare your data with function `asjob_infercnv`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_infercnv"),
  function(x, workers = 4, cutoff = .1, hmm = FALSE,
    analysis_mode = c("samples", "subclusters"), force_rerun = FALSE, ...)
  {
    step_message("Run inferCNV.")
    message(crayon::yellow(glue::glue("In linux, too many cells may result in blank heatmap.")))
    analysis_mode <- match.arg(analysis_mode)
    # https://github.com/broadinstitute/infercnv/issues/362
    if (is.remote(x)) {
      x <- run_job_remote(x, wait = 3L, ...,
        {
          x <- step1(
            x, workers = "{workers}",
            cutoff = "{cutoff}", hmm = "{hmm}",
            analysis_mode = "{analysis_mode}",
            force_rerun = "{force_rerun}"
          )
        }
      )
      return(x)
    }
    if (force_rerun && dir.exists(x$outdir)) {
      unlink(x$outdir, TRUE, TRUE)
    }
    options(scipen = 100)
    if (!requireNamespace("infercnv")) {
      setup_jags_from_libpaths()
    }
    fun_run <- function(...) {
      e(infercnv::run(
          object(x), num_threads = workers, cluster_by_groups = TRUE,
          # cutoff = 1 works well for Smart-seq2
          # and cutoff = 0.1 works well for 10x Genomics
          cutoff = cutoff, denoise = TRUE, HMM = hmm, 
          analysis_mode = analysis_mode,
          out_dir = x$outdir, save_rds = FALSE, save_final_rds = FALSE
          ))
    }
    object(x) <- expect_local_data(
      "tmp", "infercnv_run", fun_run,
      list(
        sig(x), x$metadata$cell, x$meta_used$rownames,
        x$meta_used[[x$infercnv_group_col]], cutoff, hmm, analysis_mode,
        x$max_cells_per_group
      )
    )
    x <- methodAdd(
      x,
      glue::glue(
        "inferCNV 推断中以 {bind(x$infercnv_ref)} 作为正常参考细胞，",
        "以 {bind(x$infercnv_groups)} 作为待评估细胞群；",
        "分析采用 {analysis_mode} 模式进行 CNV 信号估计，并在去噪后输出染色体尺度表达偏移热图。"
      )
    )
    return(x)
  })

setMethod("step2", signature = c(x = "job_infercnv"),
  function(x){
    step_message("Got results.")
    if (is.remote(x)) {
      dir <- file.path(x$map_local, x$outdir)
      if (!dir.exists(dir)) {
        get_file_from_remote(
          x$outdir, x$wd, x$map_local, recursive = TRUE
        )
      }
    } else {
      dir <- x$outdir
    }
    if (!requireNamespace("infercnv")) {
      setup_jags_from_libpaths()
    }
    p.infer <- .file_fig(.cut_png_blank_space(file.path(dir, "infercnv.png")))
    p.infer <- set_lab_legend(
      as_data_binary(p.infer),
      glue::glue("{x@sig} infercnv heatmap"),
      glue::glue("CNV 层次聚类热图|||{.infercnv_heatmap_method()}")
    )
    x <- plotsAdd(x, p.infer = p.infer)
    return(x)
  })

setMethod("step3", signature = c(x = "job_infercnv"),
  function(x, k = 10, clear = FALSE, p_cutoff = 0.05,
    p_adjust_method = "BH", min_cnv_z = 0.5, ref_high_quantile = 0.75,
    obs_high_quantile = 0.5, min_score_delta = 0,
    min_annotation_high_prop = 0.3, annotation_high_quantile = 0.6,
    min_annotation_n = 100L, require_pvalue = FALSE,
    use_obs_rank = TRUE, use_annotation_score = TRUE)
  {
    step_message("Kmean and rank-percentile CNV-burden scoring...")
    if (is.null(object(x))) {
      stop('is.null(object(x)).')
    }
    if (!requireNamespace("infercnv")) {
      setup_jags_from_libpaths()
    }
    if (!is.numeric(ref_high_quantile) || length(ref_high_quantile) != 1L ||
      is.na(ref_high_quantile) || ref_high_quantile <= 0 || ref_high_quantile >= 1) {
      stop('`ref_high_quantile` should be a numeric scalar between 0 and 1.')
    }
    if (!is.numeric(obs_high_quantile) || length(obs_high_quantile) != 1L ||
      is.na(obs_high_quantile) || obs_high_quantile <= 0 || obs_high_quantile >= 1) {
      stop('`obs_high_quantile` should be a numeric scalar between 0 and 1.')
    }
    if (!is.logical(use_obs_rank) || length(use_obs_rank) != 1L ||
      is.na(use_obs_rank)) {
      stop('`use_obs_rank` should be TRUE or FALSE.')
    }
    if (!is.numeric(annotation_high_quantile) ||
      length(annotation_high_quantile) != 1L ||
      is.na(annotation_high_quantile) || annotation_high_quantile <= 0 ||
      annotation_high_quantile >= 1) {
      stop('`annotation_high_quantile` should be a numeric scalar between 0 and 1.')
    }
    if (!is.numeric(min_annotation_high_prop) ||
      length(min_annotation_high_prop) != 1L ||
      is.na(min_annotation_high_prop) || min_annotation_high_prop < 0 ||
      min_annotation_high_prop > 1) {
      stop('`min_annotation_high_prop` should be a numeric scalar between 0 and 1.')
    }
    if (!is.numeric(min_annotation_n) || length(min_annotation_n) != 1L ||
      is.na(min_annotation_n) || min_annotation_n < 1) {
      stop('`min_annotation_n` should be a positive numeric scalar.')
    }
    min_annotation_n <- as.integer(min_annotation_n)
    if (!is.logical(use_annotation_score) ||
      length(use_annotation_score) != 1L || is.na(use_annotation_score)) {
      stop('`use_annotation_score` should be TRUE or FALSE.')
    }
    if (!is.logical(require_pvalue) || length(require_pvalue) != 1L ||
      is.na(require_pvalue)) {
      stop('`require_pvalue` should be TRUE or FALSE.')
    }

    .show_values <- function(x)
    {
      x <- unique(as.character(x))
      x <- x[!is.na(x) & nzchar(x)]
      if (length(x) == 0L) {
        return("无")
      }
      bind(x)
    }

    expr <- object(x)@expr.data
    obs <- unlist(
      object(x)@observation_grouped_cell_indices, use.names = FALSE
    )
    refs <- unlist(
      object(x)@reference_grouped_cell_indices, use.names = FALSE
    )

    clusters <- x$clusters <- kmeanMiniBatch(t(expr[, obs]), k)
    x$kmean_group_indices <- split(obs, paste0("C", clusters))
    x$data_group <- tibble::tibble(
      cells = colnames(expr)[obs],
      clusters = paste0("C", unname(clusters))
    )
    x$data_group <- dplyr::left_join(
      x$data_group, x$meta_used, by = c("cells" = "rownames")
    )

    # CNV burden is a continuous score summarizing genome-wide deviation
    # from the neutral inferCNV baseline. It is used as effect-size evidence.
    data_cell_score <- colMeans((expr - 1) ^ 2)
    expr_obs <- data_cell_score[obs]
    expr_refs <- data_cell_score[refs]

    data <- tibble::tibble(
      group = c(paste0("C", clusters), rep("Ref", length(expr_refs))),
      expr = c(expr_obs, expr_refs)
    )
    data <- dplyr::mutate(
      data,
      type = ifelse(group == "Ref", "Reference", "Observation"),
      type = factor(type, c("Reference", "Observation"))
    )

    ref_center <- stats::median(expr_refs, na.rm = TRUE)
    ref_scale <- stats::mad(expr_refs, constant = 1, na.rm = TRUE)
    if (is.na(ref_scale) || ref_scale <= 0) {
      ref_scale <- stats::IQR(expr_refs, na.rm = TRUE) / 1.349
    }
    if (is.na(ref_scale) || ref_scale <= 0) {
      ref_scale <- stats::sd(expr_refs, na.rm = TRUE)
    }
    if (is.na(ref_scale) || ref_scale <= 0) {
      ref_scale <- 1e-08
    }

    dataPvalue <- test_unbalance_groups(
      data, "Ref", "group", "expr", alternative = "greater"
    )
    dataPvalue$qvalue <- stats::p.adjust(dataPvalue$pvalue, method = p_adjust_method)

    data_group_score <- dplyr::summarise(
      dplyr::filter(data, type == "Observation"),
      cnv_score_median = stats::median(expr, na.rm = TRUE),
      cnv_score_mean = mean(expr, na.rm = TRUE),
      cnv_score_max = max(expr, na.rm = TRUE),
      n_cell = dplyr::n(),
      .by = group
    )
    obs_high_cutoff <- stats::quantile(
      data_group_score$cnv_score_median,
      probs = obs_high_quantile,
      na.rm = TRUE,
      names = FALSE
    )
    n_group_score <- nrow(data_group_score)
    data_group_score$cnv_score_rank <- rank(
      data_group_score$cnv_score_median,
      ties.method = "average"
    )
    if (n_group_score <= 1L) {
      data_group_score$cnv_score_rank_prop <- 1
    } else {
      data_group_score$cnv_score_rank_prop <-
        (data_group_score$cnv_score_rank - 1) / (n_group_score - 1)
    }
    ref_tail_cutoff <- stats::quantile(
      expr_refs,
      probs = ref_high_quantile,
      na.rm = TRUE,
      names = FALSE
    )
    ref_z_cutoff <- ref_center + min_cnv_z * ref_scale

    dataPvalue <- dplyr::left_join(dataPvalue, data_group_score, by = "group")
    dataPvalue <- dplyr::mutate(
      dataPvalue,
      ref_score_median = ref_center,
      ref_score_scale = ref_scale,
      ref_z_cutoff = ref_z_cutoff,
      ref_tail_cutoff = ref_tail_cutoff,
      obs_high_cutoff = obs_high_cutoff,
      cnv_score_delta = cnv_score_median - ref_center,
      cnv_z = cnv_score_delta / ref_scale,
      pass_pvalue = qvalue < p_cutoff,
      pass_ref_z = cnv_z >= min_cnv_z & cnv_score_delta >= min_score_delta,
      pass_ref_tail = cnv_score_median >= ref_tail_cutoff &
        cnv_score_delta >= min_score_delta,
      pass_ref_effect = pass_ref_z | pass_ref_tail,
      pass_obs_rank = if (use_obs_rank) {
        cnv_score_rank_prop >= obs_high_quantile
      } else {
        TRUE
      },
      pass_stat = if (require_pvalue) {
        pass_pvalue
      } else {
        TRUE
      },
      is_cnv_high = pass_stat & pass_ref_effect & pass_obs_rank,
      cnv_call = ifelse(is_cnv_high, "CNV-high", "CNV-low"),
      sig = cnv_call,
      type = factor("Observation", levels(data$type))
    )

    x$dataPvalue <- set_lab_legend(
      dataPvalue,
      glue::glue("{x@sig} kmean cluster CNV-burden result"),
      glue::glue(
        "CNV pattern 分组判定表|||该表展示 `kmean` CNV pattern 分组的 CNV burden、",
        "相对 Reference 的稳健偏移程度、Observation 内部 CNV burden 分位排名、",
        "校正后显著性以及 CNV-high/CNV-low 判定结果。"
      )
    )

    x$data_group <- dplyr::left_join(
      x$data_group,
      dplyr::select(dataPvalue, group, cnv_call, is_cnv_high),
      by = c("clusters" = "group")
    )
    x$data_group <- dplyr::mutate(
      x$data_group,
      isCancer = ifelse(is_cnv_high, "Malignant cell", "CNV-low cell")
    )

    vec_group_col <- setdiff(
      colnames(x$data_group),
      c("cells", "clusters", "cnv_call", "is_cnv_high", "isCancer")
    )
    annotation_col <- vec_group_col[1L]
    data_annotation_call <- dplyr::summarise(
      x$data_group,
      cnv_high_prop = mean(isCancer == "Malignant cell", na.rm = TRUE),
      n_infercnv_cell = dplyr::n(),
      .by = !!rlang::sym(annotation_col)
    )
    data_annotation_cell_score <- tibble::tibble(
      cells = colnames(expr)[obs],
      cnv_score = as.numeric(expr_obs)
    )
    data_annotation_cell_score <- dplyr::left_join(
      data_annotation_cell_score,
      x$meta_used,
      by = c("cells" = "rownames")
    )
    data_annotation_score <- dplyr::summarise(
      data_annotation_cell_score,
      cnv_score_median = stats::median(cnv_score, na.rm = TRUE),
      cnv_score_mean = mean(cnv_score, na.rm = TRUE),
      .by = !!rlang::sym(annotation_col)
    )
    data_annotation_call <- dplyr::left_join(
      data_annotation_call,
      data_annotation_score,
      by = annotation_col
    )
    n_annotation_score <- nrow(data_annotation_call)
    data_annotation_call$annotation_score_rank <- rank(
      data_annotation_call$cnv_score_median,
      ties.method = "average"
    )
    if (n_annotation_score <= 1L) {
      data_annotation_call$annotation_score_rank_prop <- 1
    } else {
      data_annotation_call$annotation_score_rank_prop <-
        (data_annotation_call$annotation_score_rank - 1) /
        (n_annotation_score - 1)
    }
    annotation_high_cutoff <- stats::quantile(
      data_annotation_call$cnv_score_median,
      probs = annotation_high_quantile,
      na.rm = TRUE,
      names = FALSE
    )
    data_annotation_call <- dplyr::mutate(
      data_annotation_call,
      annotation_high_cutoff = annotation_high_cutoff,
      pass_annotation_prop = cnv_high_prop >= min_annotation_high_prop,
      pass_annotation_score = if (use_annotation_score) {
        annotation_score_rank_prop >= annotation_high_quantile
      } else {
        TRUE
      },
      pass_annotation_n = n_infercnv_cell >= min_annotation_n,
      infercnv_call = ifelse(
        pass_annotation_prop & pass_annotation_score & pass_annotation_n,
        "CNV-high", "CNV-low"
      ),
      isCancer = ifelse(
        infercnv_call == "CNV-high",
        "Malignant cell", "CNV-low cell"
      )
    )

    x$data_annotation_call <- set_lab_legend(
      data_annotation_call,
      glue::glue("{x@sig} annotation-group CNV-burden call"),
      glue::glue(
        "注释组 CNV 判定表|||该表按 inferCNV 注释组汇总 CNV-high 细胞比例，",
        "并基于实际纳入 inferCNV 的该组细胞计算组级 CNV burden 中位数及其在各注释组中的相对分层。",
        "仅当组内 CNV-high 细胞比例达到 {min_annotation_high_prop}、",
        "组级 CNV burden 位于注释组高分层（分位阈值 {annotation_high_quantile}），",
        "且纳入细胞数不少于 {min_annotation_n} 时，才将该注释组判定为潜在肿瘤细胞群；",
        "未达到标准的注释组保留为 CNV-low 细胞群。"
      )
    )

    data_kmean_annotation <- dplyr::summarise(
      x$data_group,
      n_cell = dplyr::n(),
      .by = c(!!rlang::sym(annotation_col), clusters, cnv_call)
    )
    data_kmean_annotation <- dplyr::mutate(
      data_kmean_annotation,
      annotation_group = as.character(!!rlang::sym(annotation_col)),
      kmean_group = as.character(clusters)
    )
    data_annotation_total <- dplyr::summarise(
      x$data_group,
      n_annotation_cell = dplyr::n(),
      .by = !!rlang::sym(annotation_col)
    )
    data_annotation_total <- dplyr::mutate(
      data_annotation_total,
      annotation_group = as.character(!!rlang::sym(annotation_col))
    )
    data_kmean_grid <- base::expand.grid(
      annotation_group = unique(data_kmean_annotation$annotation_group),
      kmean_group = unique(as.character(dataPvalue$group)),
      stringsAsFactors = FALSE
    )
    data_kmean_grid <- dplyr::left_join(
      data_kmean_grid,
      dplyr::select(dataPvalue, kmean_group = group, cnv_call),
      by = "kmean_group"
    )
    data_kmean_annotation <- dplyr::left_join(
      data_kmean_grid,
      dplyr::select(
        data_kmean_annotation,
        annotation_group, kmean_group, cnv_call, n_cell
      ),
      by = c("annotation_group", "kmean_group", "cnv_call")
    )
    data_kmean_annotation$n_cell[is.na(data_kmean_annotation$n_cell)] <- 0L
    data_kmean_annotation <- dplyr::left_join(
      data_kmean_annotation,
      dplyr::select(data_annotation_total, annotation_group, n_annotation_cell),
      by = "annotation_group"
    )
    data_kmean_annotation <- dplyr::mutate(
      data_kmean_annotation,
      prop_in_annotation = n_cell / n_annotation_cell,
      clusters = kmean_group
    )
    data_kmean_annotation[[annotation_col]] <- data_kmean_annotation$annotation_group
    x$data_kmean_annotation_link <- set_lab_legend(
      data_kmean_annotation,
      glue::glue("{x@sig} kmeans-annotation CNV link table"),
      glue::glue(
        "k-means 与注释组对应表|||该表展示每个 inferCNV 注释组中不同 k-means CNV pattern 分组的细胞数及比例，",
        "未观察到对应细胞的组合以 0 保留，用于说明抽象的 CNV pattern 分组如何对应到 Seurat cluster / annotation group。"
      )
    )

    data_annotation_plot <- dplyr::left_join(
      data_annotation_cell_score,
      dplyr::select(
        data_annotation_call,
        !!rlang::sym(annotation_col), infercnv_call, isCancer,
        annotation_high_cutoff
      ),
      by = annotation_col
    )
    data_annotation_plot <- dplyr::mutate(
      data_annotation_plot,
      annotation_group = as.character(!!rlang::sym(annotation_col))
    )
    data_annotation_call_plot <- dplyr::mutate(
      data_annotation_call,
      annotation_group = as.character(!!rlang::sym(annotation_col))
    )
    vec_annotation_level <- data_annotation_call[[annotation_col]][
      order(data_annotation_call$cnv_score_median)
    ]
    data_annotation_plot$annotation_group_plot <- factor(
      data_annotation_plot$annotation_group,
      levels = vec_annotation_level
    )
    data_annotation_call_plot$annotation_group_plot <- factor(
      data_annotation_call_plot$annotation_group,
      levels = vec_annotation_level
    )
    data_kmean_annotation$annotation_group_plot <- factor(
      data_kmean_annotation$annotation_group,
      levels = vec_annotation_level
    )

    y_pos <- max(data$expr, na.rm = TRUE)
    p.violin <- ggplot(data) +
      geom_violin(aes(x = reorder(group, expr), y = expr, fill = group)) +
      geom_text(data = dataPvalue, aes(x = group, label = cnv_call, y = y_pos)) +
      ggforce::facet_row(~ type, space = "free", scales = "free_x") +
      scale_fill_manual(values = color_set()) +
      theme_minimal()
    p.violin <- wrap(p.violin, 8, 3.5)
    p.violin <- set_lab_legend(
      p.violin,
      glue::glue("{x@sig} kmean cluster CNV burden violin plot"),
      glue::glue(
        "CNV burden 分层小提琴图|||基于 inferCNV 去噪后的表达偏移矩阵，",
        "按细胞计算 CNV burden（相对中性基线 1 的均方偏移），并对 k-means CNV pattern 分组与 Reference 进行比较。",
        "CNV-high 判定以 CNV burden 为主要效应量证据，结合 Reference 稳健基线、Reference 高分位阈值以及 Observation 内部 CNV burden 分位排名进行筛选。"
      )
    )

    p.annotation_violin <- ggplot(data_annotation_plot) +
      geom_violin(aes(
        x = annotation_group_plot,
        y = cnv_score,
        fill = infercnv_call
      )) +
      geom_hline(
        yintercept = ref_center,
        linetype = "dashed",
        linewidth = .3
      ) +
      geom_hline(
        yintercept = ref_tail_cutoff,
        linetype = "dotted",
        linewidth = .3
      ) +
      labs(x = NULL, y = "CNV burden") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    p.annotation_violin <- wrap(p.annotation_violin, 8, 3.8)
    p.annotation_violin <- set_lab_legend(
      p.annotation_violin,
      glue::glue("{x@sig} annotation-group CNV burden violin plot"),
      glue::glue(
        "注释组 CNV burden 分布图|||该图直接展示每个 inferCNV 注释组中实际纳入细胞的 CNV burden 分布，",
        "虚线表示 Reference 细胞群的 CNV burden 中位数，点线表示 Reference 高分位阈值。",
        "该图用于直观评估 Seurat cluster / annotation group 层面的 CNV 信号强弱。"
      )
    )

    p.kmean_annotation <- ggplot(data_kmean_annotation) +
      geom_tile(aes(x = kmean_group, y = annotation_group_plot, fill = prop_in_annotation)) +
      geom_text(aes(
        x = kmean_group,
        y = annotation_group_plot,
        label = sprintf("%.2f", prop_in_annotation)
      ), size = 2.6, color = "white") +
      labs(x = "k-means CNV pattern group", y = NULL, fill = "Proportion") +
      theme_minimal()
    p.kmean_annotation <- wrap(p.kmean_annotation, 8, 4)
    p.kmean_annotation <- set_lab_legend(
      p.kmean_annotation,
      glue::glue("{x@sig} kmeans-annotation CNV pattern heatmap"),
      glue::glue(
        "k-means 与注释组对应热图|||该图展示每个 Seurat cluster / annotation group 中细胞分配到各 k-means CNV pattern 分组的比例。",
        "没有细胞落入的组合显示为 0，用于追踪 CNV-high k-means 分组主要来源于哪些注释组，从而使肿瘤细胞判定过程更加透明。"
      )
    )

    p.annotation_scatter <- ggplot(data_annotation_call_plot) +
      geom_vline(
        xintercept = min_annotation_high_prop,
        linetype = "dashed",
        linewidth = .3
      ) +
      geom_hline(
        yintercept = annotation_high_cutoff,
        linetype = "dashed",
        linewidth = .3
      ) +
      geom_point(aes(
        x = cnv_high_prop,
        y = cnv_score_median,
        size = n_infercnv_cell,
        color = infercnv_call
      ), alpha = .85) +
      geom_text(aes(
        x = cnv_high_prop,
        y = cnv_score_median,
        label = annotation_group
      ), vjust = -0.8, size = 3) +
      labs(
        x = "Proportion of CNV-high k-means cells",
        y = "Annotation-group median CNV burden",
        size = "Cells"
      ) +
      theme_minimal()
    p.annotation_scatter <- wrap(p.annotation_scatter, 7, 4.2)
    p.annotation_scatter <- set_lab_legend(
      p.annotation_scatter,
      glue::glue("{x@sig} annotation-group CNV-high decision plot"),
      glue::glue(
        "注释组 CNV 判定散点图|||该图以每个 Seurat cluster / annotation group 为单位，",
        "同时展示 CNV-high k-means 细胞比例和组级 CNV burden 中位数。",
        "垂直虚线表示 CNV-high 细胞比例阈值，水平虚线表示注释组 CNV burden 高分层阈值；",
        "位于双阈值高区且满足最低细胞数要求的注释组被判定为 CNV-high 潜在肿瘤细胞群。"
      )
    )

    x <- tablesAdd(
      x,
      t.cnv_kmeans = x$dataPvalue,
      t.cnv_annotation_call = x$data_annotation_call,
      t.kmean_annotation_link = x$data_kmean_annotation_link
    )
    x <- plotsAdd(
      x,
      p.violin = p.violin,
      p.annotation_violin = p.annotation_violin,
      p.kmean_annotation = p.kmean_annotation,
      p.annotation_scatter = p.annotation_scatter
    )

    vec_high_k <- dataPvalue$group[dataPvalue$is_cnv_high]
    vec_low_k <- dataPvalue$group[!dataPvalue$is_cnv_high]
    vec_high_anno <- data_annotation_call[[annotation_col]][
      data_annotation_call$infercnv_call == "CNV-high"
    ]
    vec_low_anno <- data_annotation_call[[annotation_col]][
      data_annotation_call$infercnv_call == "CNV-low"
    ]

    .show_annotation_detail <- function(data, call)
    {
      data <- dplyr::filter(data, infercnv_call == call)
      if (nrow(data) == 0L) {
        return("无")
      }
      data <- dplyr::arrange(data, dplyr::desc(cnv_score_median))
      vec_label <- paste0(
        as.character(data[[annotation_col]]),
        " (CNV-high=", sprintf("%.1f%%", 100 * data$cnv_high_prop),
        ", burden=", signif(data$cnv_score_median, 3),
        ", n=", data$n_infercnv_cell, ")"
      )
      bind(vec_label)
    }

    txt_high_anno_detail <- .show_annotation_detail(
      data_annotation_call, "CNV-high"
    )
    txt_low_anno_detail <- .show_annotation_detail(
      data_annotation_call, "CNV-low"
    )

    x <- methodAdd(
      x,
      glue::glue(
        "为避免仅由大样本量导致的统计显著性造成过度判定，本研究未直接依据 p 值将所有候选细胞定义为肿瘤细胞。",
        "参考 scRNA-seq CNV 推断识别恶性细胞的常用思路，分析基于 inferCNV 表达偏移矩阵计算每个细胞的 CNV burden，",
        "并以参考细胞群建立稳健基线；随后结合相对 Reference 的 CNV burden 增幅、Reference 高分位阈值以及 Observation 细胞群内部的 CNV burden 分位排名，",
        "筛选 CNV-high 的候选肿瘤细胞群。在注释组层面，进一步同时考虑组内 CNV-high 细胞比例、",
        "该组实际纳入 inferCNV 细胞的组级 CNV burden 分层以及最低细胞数要求，",
        "以减少由少量高 CNV 细胞或抽样波动造成的过度判定。未达到 CNV-high 标准的细胞群保留为 CNV-low 细胞群，不额外定义为肿瘤细胞。",
        "该策略与基于 scRNA-seq 推断 CNV/非整倍体信号区分恶性细胞和非恶性细胞的分析思路一致（PMID: 33462507; 39200223）。"
      )
    )
    x <- snapAdd(
      x,
      glue::glue(
        "inferCNV CNV burden 分层显示，CNV-high k-means 分组为：{.show_values(vec_high_k)}；",
        "CNV-low k-means 分组为：{.show_values(vec_low_k)}。",
        "按注释组汇总后，判定为潜在肿瘤细胞的 cluster 为：{.show_values(vec_high_anno)}；",
        "保留为 CNV-low 细胞群的 cluster 为：{.show_values(vec_low_anno)}。",
        "其中，潜在肿瘤细胞 cluster 的判定依据为：{txt_high_anno_detail}；",
        "CNV-low cluster 的对应指标为：{txt_low_anno_detail}。",
        "本步骤同时输出 k-means CNV pattern 与 Seurat cluster / annotation group 的对应热图{aref(p.kmean_annotation)}、",
        "注释组 CNV burden 分布图和注释组判定散点图{aref(p.annotation_scatter)}，用于展示判定过程。"
      )
    )

    if (clear) {
      object(x) <- NULL
    }
    return(x)
  })

setMethod("map", signature = c(x = "job_seurat", ref = "job_infercnv"),
  function(x, ref, from = x$group.by, to = "infercnv_cell",
    map_by = c("cell", "annotation_group"), group_col = NULL,
    subset = NULL, min_prop = 0.5)
  {
    if (ref@step < 3) {
      stop('ref@step < 3.')
    }
    if (!requireNamespace("infercnv")) {
      setup_jags_from_libpaths()
    }
    map_by <- match.arg(map_by)
    res <- ref$data_group

    if ("isCancer" %in% colnames(res)) {
      res$isCancer <- ifelse(
        res$isCancer == "Malignant cell",
        "Malignant cell", "CNV-low cell"
      )
    } else {
      group.sig <- dplyr::filter(ref$dataPvalue, is_cnv_high)$group
      res <- dplyr::mutate(
        res,
        isCancer = ifelse(
          clusters %in% group.sig,
          "Malignant cell", "CNV-low cell"
        )
      )
    }
    message(glue::glue("Is Malignant: {try_snap(res$isCancer == 'Malignant cell')}"))

    if (map_by == "cell") {
      isCancer <- res$isCancer[ match(rownames(object(x)@meta.data), res$cells) ]
      message(glue::glue("Got cell-level annotation: {try_snap(!is.na(isCancer))}"))
    } else {
      if (!is.null(ref$data_annotation_call)) {
        data_group_call <- ref$data_annotation_call
        if (is.null(group_col)) {
          vec_group_col <- setdiff(
            colnames(data_group_call),
            c(
              "cnv_high_prop", "n_infercnv_cell", "cnv_score_median",
              "cnv_score_mean", "annotation_score_rank",
              "annotation_score_rank_prop", "annotation_high_cutoff",
              "pass_annotation_prop", "pass_annotation_score",
              "pass_annotation_n", "infercnv_call", "isCancer"
            )
          )
          group_col <- vec_group_col[1L]
        }
        data_group_call$isCancer <- ifelse(
          data_group_call$isCancer == "Malignant cell",
          "Malignant cell", "CNV-low cell"
        )
      } else {
        if (is.null(group_col)) {
          vec_group_col <- setdiff(
            colnames(res), c("cells", "clusters", "cnv_call", "is_cnv_high", "isCancer")
          )
          if (length(vec_group_col) == 0L) {
            stop('No inferCNV annotation group column was found in `ref$data_group`.')
          }
          group_col <- vec_group_col[1L]
        }
        if (!group_col %in% colnames(res)) {
          stop('`group_col` was not found in `ref$data_group`.')
        }
        data_group_call <- dplyr::summarise(
          res,
          malignant_prop = mean(isCancer == "Malignant cell"),
          n_infercnv_cell = dplyr::n(),
          .by = !!rlang::sym(group_col)
        )
        data_group_call <- dplyr::mutate(
          data_group_call,
          isCancer = ifelse(
            malignant_prop >= min_prop,
            "Malignant cell", "CNV-low cell"
          )
        )
      }

      if (!group_col %in% colnames(data_group_call)) {
        stop('`group_col` was not found in annotation-group inferCNV result.')
      }
      metadata <- as_tibble(object(x)@meta.data)
      if (!is.null(subset)) {
        subset_use <- subset
      } else {
        subset_use <- ref$infercnv_subset
      }
      if (group_col %in% colnames(metadata)) {
        vec_metadata_group <- as.character(metadata[[group_col]])
      } else if (group_col == "group_subset" && !is.null(subset_use) &&
        subset_use %in% colnames(metadata) && from %in% colnames(metadata)) {
        vec_metadata_group <- paste0(metadata[[from]], "_", metadata[[subset_use]])
        vec_metadata_group <- ifelse(
          metadata[[from]] %in% ref$infercnv_ref,
          as.character(metadata[[from]]), vec_metadata_group
        )
      } else if (ref$infercnv_group_col %in% colnames(metadata)) {
        vec_metadata_group <- as.character(metadata[[ref$infercnv_group_col]])
      } else {
        stop('Unable to reconstruct inferCNV annotation groups in Seurat metadata.')
      }
      isCancer <- data_group_call$isCancer[
        match(vec_metadata_group, data_group_call[[group_col]])
      ]
      message(glue::glue("Got annotation-group-level annotation: {try_snap(!is.na(isCancer))}"))
      x$infercnv_group_call <- data_group_call
    }

    object(x)@meta.data[[ to ]] <- ifelse(
      is.na(isCancer) | isCancer != "Malignant cell",
      as.character(object(x)@meta.data[[ from ]]),
      "Malignant cell"
    )
    palette <- .set_palette_in_ending(
      object(x)@meta.data[[ from ]], "Malignant cell"
    )
    object(x)@meta.data[[ to ]] <- factor(
      object(x)@meta.data[[ to ]], names(palette)
    )
    p.map_cancer <- e(Seurat::DimPlot(
        object(x), reduction = "umap", label = FALSE,
        group.by = to, cols = palette
        ))
    p.map_cancer <- wrap(as_grob(p.map_cancer), 7, 4)
    p.map_cancer <- set_lab_legend(
      p.map_cancer,
      glue::glue("{x@sig} cell type annotation with inferCNV idenfication"),
      glue::glue("细胞注释结果的 UMAP 图|||`inferCNV` 的细胞级 CNV-high 判定结果映射至 Seurat 细胞注释。不同颜色代表不同细胞簇类型，横纵坐标为UMAP的两个维度。")
    )
    x$p.map_cancer <- p.map_cancer
    if (map_by == "cell") {
      x <- snapAdd(
        x,
        "将 `inferCNV` 的细胞级 CNV-high 判定结果映射至 Seurat 细胞注释中；未达到 CNV-high 标准的细胞保留原有细胞注释。",
        step = class(ref), add = FALSE
      )
    } else {
      x <- snapAdd(
        x,
        glue::glue(
          "将 `inferCNV` 的注释组级 CNV-high 判定结果映射至 Seurat 细胞注释中；",
          "达到 CNV-high 标准的注释组标记为 Malignant cell，未达到标准的注释组保留原有细胞注释。"
        ),
        step = class(ref), add = FALSE
      )
    }
    x$.map_heading <- glue::glue("Seurat-InferCNV 癌细胞注释")
    return(x)
  })

test_unbalance_groups <- function(data, ref, col.group = "group",
  col.value = "value", fun = function(...) wilcox.test(...)$p.value, ...)
{
  fun <- match.fun(fun)
  data <- split(data, data[[col.group]] == ref)
  if (any(vapply(data, nrow, integer(1)) <= 1)) {
    stop('any(vapply(data, nrow, integer(1)) <= 1).')
  }
  fun_test <- function(x, ...) {
    fun(x, data[["TRUE"]][[col.value]], ...)
  }
  res <- dplyr::reframe(
    data[["FALSE"]], pvalue = fun_test(
      !!rlang::sym(col.value), ...
      ),
    max = max(!!rlang::sym(col.value), na.rm = TRUE),
    min = min(!!rlang::sym(col.value), na.rm = TRUE),
    .by = !!rlang::sym(col.group)
  )
  res <- dplyr::mutate(res,
    sig = ifelse(
      pvalue < .001, "***", ifelse(
        pvalue < .01, "**", ifelse(pvalue < .05, "*", "")
      )
    )
    # sig = paste0(sig, " (", signif(pvalue, 4), ")")
  )
  res
}

kmeanMiniBatch <- function(mtx, k = 10, batch = 100, force = FALSE, ...) {
  hash <- digest::digest(list(mtx, k, batch, force))
  file_cache <- add_filename_suffix("kmean.rds", hash)
  if (file.exists(file_cache)) {
    message(glue::glue("Read from cache: {file_cache}"))
    clusters <- readRDS(file_cache)
  } else {
    if (nrow(mtx) < 1e6 && !force) {
      message("Use `kmeans` for small size data.")
      clusters <- kmeans(mtx, k)
    } else {
      if (!requireNamespace("ClusterR", quietly = TRUE)) {
        install.packages("ClusterR")
      }
      object <- e(ClusterR::MiniBatchKmeans(
        mtx, k, batch_size = batch, num_init = 5
      ))
      clusters <- e(ClusterR::predict_MBatchKMeans(mtx, object$centroids))
    }
    saveRDS(clusters, file_cache)
  }
  if (is(clusters, "kmeans")) {
    clusters <- clusters$cluster
  }
  clusters
}

.cut_png_blank_space <- function(file_png, threshold = .99, 
  max = .01)
{
  img <- png::readPNG(file_png)
  dims <- dim(img)
  is_rgb <- length(dims) >= 3 && (dims[3] %in% c(3, 4))
  if (is_rgb) {
    rgb <- img[,,1:3]
    is_white <- (rgb[,,1] >= threshold) & (rgb[,,2] >= threshold) & (rgb[,,3] >= threshold)
  } else {
    gray <- img[,,1]
    is_white <- (gray >= threshold)
  }
  white_rows <- apply(is_white, 1, all)
  n <- 1L
  status <- white_rows[1]
  group <- rep(0L, length(white_rows))
  for (i in seq_along(white_rows)) {
    if (status != white_rows[i]) {
      n <- n + 1L
      status <- !status
    }
    group[i] <- n
  }
  groups <- split(seq_along(white_rows), group)
  ratio <- lengths(groups) / length(white_rows)
  allStatus <- as.logical(seq_along(ratio) %% 2)
  if (!white_rows[1]) {
    allStatus <- !allStatus
  }
  areThat <- ratio > max & allStatus
  maxNum <- floor(length(white_rows) * max)
  groups[areThat] <- lapply(groups[areThat], 
    function(x) {
      head(x, maxNum)
    })
  keep <- unlist(groups)
  if (is_rgb) {
    cropped_img <- img[keep, , , drop = FALSE]
  } else {
    cropped_img <- img[keep, , drop = FALSE]
  }
  newfile <- add_filename_suffix(file_png, "crop")
  png::writePNG(cropped_img, newfile)
  return(newfile)
}

.infercnv_heatmap_method <- function() {
  "正常细胞的表达值绘制在顶部热图中，肿瘤细胞的表达值绘制在底部热图中，基因在染色体上从左到右排列。正常细胞表达数据实际上从肿瘤细胞表达数据中减去，从而得出差异，其中染色体区域扩增显示为红色块，染色体区域缺失显示为蓝色块。参考<https://github.com/broadinstitute/inferCNV/wiki/Interpreting-the-figure>"
}

setMethod("set_remote", signature = c(x = "job_infercnv"),
  function(x, wd = glue::glue("~/infercnv_{x@sig}")){
    x$wd <- wd
    rem_dir.create(wd, wd = ".")
    return(x)
  })


get_gene_ranges <- function(symbols, version = c("hg38", 
  "hg19"), gname = TRUE)
{
  raw <- symbols
  if (gname) {
    symbols <- gname(symbols)
  }
  version <- match.arg(version)
  name_fun <- glue::glue("TxDb.Hsapiens.UCSC.{version}.knownGene")
  if (!requireNamespace(name_fun, quietly = TRUE)) {
    BiocManager::install(name_fun)
  }
  db <- get_fun(name_fun, asNamespace(name_fun), "S4")
  ranges <- e(GenomicFeatures::genes(db))
  # if (!requireNamespace("EnsDb.Hsapiens.v86")) {
  #   BiocManager::install("EnsDb.Hsapiens.v86")
  #   # EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
  # }
  entrez <- e(AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keytype = "ALIAS",
    keys = symbols,
    column = c("ENTREZID")
  ))
  hasThats <- !is.na(entrez) & (entrez %in% ranges$gene_id)
  entrez <- entrez[ hasThats ]
  ranges <- ranges[ match(entrez, ranges$gene_id) ]
  ranges <- data.frame(ranges)
  ranges$symbols <- raw[ hasThats ]
  return(ranges)
}

setup_jags_from_libpaths <- function(
  extra_libpaths = NULL,
  prefer_env = NULL,
  set_libpaths = TRUE,
  set_path = FALSE,
  verbose = TRUE)
{
  .clean_path <- function(path)
  {
    gsub("/+$", "", normalizePath(path, mustWork = FALSE))
  }

  .prefix_from_r_library <- function(path)
  {
    path <- .clean_path(path)

    if (basename(path) == "library" &&
        basename(dirname(path)) == "R" &&
        basename(dirname(dirname(path))) == "lib") {
      return(dirname(dirname(dirname(path))))
    }

    return(NA_character_)
  }

  .prefix_from_package_path <- function(path)
  {
    path <- .clean_path(path)

    if (basename(dirname(path)) == "library" &&
        basename(dirname(dirname(path))) == "R" &&
        basename(dirname(dirname(dirname(path)))) == "lib") {
      return(dirname(dirname(dirname(dirname(path)))))
    }

    return(NA_character_)
  }

  .is_valid_jags_moddir <- function(path)
  {
    file.exists(file.path(path, "basemod.so")) &&
      file.exists(file.path(path, "bugs.so"))
  }

  if (!is.null(extra_libpaths) && set_libpaths) {
    .libPaths(c(extra_libpaths, .libPaths()))
  }

  vec_lib <- unique(.clean_path(.libPaths()))

  vec_prefix <- stats::na.omit(vapply(
    vec_lib,
    .prefix_from_r_library,
    character(1L)
  ))

  path_rjags <- unlist(lapply(
    vec_lib,
    function(path_lib) {
      path_pkg <- file.path(path_lib, "rjags")
      if (dir.exists(path_pkg)) {
        return(path_pkg)
      }
      character()
    }
  ))

  if (length(path_rjags) > 0L) {
    vec_prefix_pkg <- stats::na.omit(vapply(
      path_rjags,
      .prefix_from_package_path,
      character(1L)
    ))
    vec_prefix <- c(vec_prefix_pkg, vec_prefix)
  }

  path_jags <- Sys.which("jags")
  if (nchar(path_jags)) {
    vec_prefix <- c(dirname(dirname(.clean_path(path_jags))), vec_prefix)
  }

  vec_prefix <- unique(.clean_path(vec_prefix))

  if (!is.null(prefer_env)) {
    idx_prefer <- grepl(prefer_env, vec_prefix)
    vec_prefix <- c(vec_prefix[idx_prefer], vec_prefix[!idx_prefer])
    vec_prefix <- unique(vec_prefix)
  }

  data_candidate <- do.call(
    rbind,
    lapply(
      vec_prefix,
      function(prefix) {
        data.frame(
          prefix = prefix,
          moddir = c(
            file.path(prefix, "lib/JAGS/modules-4"),
            file.path(prefix, "lib64/JAGS/modules-4"),
            file.path(prefix, "JAGS/modules-4")
          ),
          stringsAsFactors = FALSE
        )
      }
    )
  )

  data_candidate$valid <- vapply(
    data_candidate$moddir,
    .is_valid_jags_moddir,
    logical(1L)
  )

  if (!any(data_candidate$valid)) {
    msg <- paste(
      "No valid JAGS module directory was found.",
      "Checked candidates:",
      paste(data_candidate$moddir, collapse = "\n  "),
      sep = "\n  "
    )
    stop(msg, call. = FALSE)
  }

  data_valid <- data_candidate[data_candidate$valid, , drop = FALSE]
  path_jags_mod <- data_valid$moddir[1L]
  env_prefix <- data_valid$prefix[1L]
  path_r_lib <- file.path(env_prefix, "lib/R/library")

  options(jags.moddir = path_jags_mod)

  Sys.setenv(
    JAGS_MODULE_PATH = path_jags_mod,
    JAGS_HOME = env_prefix,
    CONDA_PREFIX = env_prefix
  )

  if (set_path) {
    path_bin <- file.path(env_prefix, "bin")
    path_lib <- file.path(env_prefix, "lib")

    Sys.setenv(
      PATH = paste(path_bin, Sys.getenv("PATH"), sep = .Platform$path.sep),
      LD_LIBRARY_PATH = paste(path_lib, Sys.getenv("LD_LIBRARY_PATH"), sep = .Platform$path.sep)
    )
  }

  if (set_libpaths && dir.exists(path_r_lib)) {
    .libPaths(c(path_r_lib, .libPaths()))
  }

  if (verbose) {
    message("JAGS prefix: ", env_prefix)
    message("JAGS module dir: ", path_jags_mod)
    message("R library: ", path_r_lib)
    message("options('jags.moddir'): ", getOption("jags.moddir"))
  }

  invisible(list(
    prefix = env_prefix,
    moddir = path_jags_mod,
    r_library = path_r_lib,
    candidates = data_candidate
  ))
}

