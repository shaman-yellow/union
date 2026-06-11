# ==========================================================================
# workflow of seurat5n
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_seurat5n <- setClass("job_seurat5n", 
  contains = c("job_seurat"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    pg = "seurat5n",
    info = c("https://satijalab.org/seurat/articles/integration_introduction"),
    cite = "[@DictionaryLearHaoY2024]",
    method = "",
    tag = "seurat5n",
    analysis = "Seurat 集成单细胞数据分析"
    ))

job_seurat5n <- function(dirs, names = NULL, mode = c("sc", "st"), 
  st.filename = "filtered_feature_bc_matrix.h5", ...)
{
  # https://satijalab.org/seurat/articles/parsebio_sketch_integration
  mode <- match.arg(mode)
  n <- 0L
  object <- pbapply::pblapply(dirs,
    function(dir) {
      n <<- n + 1L
      project <- names[n]
      if (mode == "sc") {
        suppressMessages(job_seurat(dir, project = project, ...))@object
      } else {
        suppressMessages(job_seuratSp(dir, filename = st.filename, ...))@object
      }
    })
  x <- .job_seurat5n(object = object)
  # if (!is.null(names)) {
  #   x <- snapAdd(x, "读取 {bind(names)} 样本的数据集。")
  # }
  object(x) <- e(SeuratObject:::merge.Seurat(object(x)[[1]], object(x)[-1]))
  object(x)[[ "percent.mt" ]] <- e(Seurat::PercentageFeatureSet(object(x), pattern = "^MT-"))
  p.qc_pre <- plot_qc.seurat(x)
  x$p.qc_pre <- p.qc_pre
  x <- methodAdd(x, "以 R 包 `Seurat` ⟦pkgInfo('Seurat')⟧ 进行单细胞数据质量控制 (QC) 和下游分析。依据 <{x@info}> 为指导对单细胞数据预处理。")
  return(x)
}

setGeneric("asjob_seurat5n",
  function(x, ...) standardGeneric("asjob_seurat5n"))

setMethod("asjob_seurat5n", signature = c(x = "job_seurat"),
  function(x, split = "orig.ident", assay = "RNA")
  {
    object(x)[[assay]] <- split(object(x)[[assay]], f = object(x)@meta.data[[split]])
    x <- .job_seurat5n(object = object(x))
    object(x)[[ "percent.mt" ]] <- e(Seurat::PercentageFeatureSet(object(x), pattern = "^MT-"))
    SeuratObject::Idents(object(x)) <- split
    p.qc_pre <- plot_qc.seurat(x)
    x$p.qc_pre <- p.qc_pre
    x <- methodAdd(x, "以 R 包 `Seurat` ⟦pkgInfo('Seurat')⟧ 进行单细胞数据质量控制 (QC) 和下游分析。依据 <{x@info}> 为指导对单细胞数据预处理。")
    return(x)
  })

setMethod("step0", signature = c(x = "job_seurat5n"),
  function(x){
    step_message("Prepare your data with function `job_seurat5n`.")
  })

setMethod("map", signature = c(x = "job_seurat", ref = "df"),
  function(x, ref, by.x = "orig.ident", by.ref = "sample", 
    get = "group", col = get)
  {
    object(x)@meta.data[[col]] <- dplyr::recode(
      object(x)@meta.data[[by.x]], !!!setNames(ref[[get]], ref[[by.ref]])
    )
    return(x)
  })

setMethod("step1", signature = c(x = "job_seurat5n"),
  function(x, min.features, max.features, max.count, max.percent.mt = 5)
  {
    step_message("Quality control (QC).")
    if (!is.null(min.features)) {
      ncell <- ncol(object(x))
      ngene <- nrow(object(x))
      object(x) <- e(SeuratObject:::subset.Seurat(
          object(x), subset = nFeature_RNA > min.features &
            nFeature_RNA < max.features & percent.mt < max.percent.mt &
            nCount_RNA < max.count
          ))
      p.qc_aft <- plot_qc.seurat(x)
      x$p.qc_aft <- p.qc_aft <- set_lab_legend(
        p.qc_aft,
        glue::glue("{x@sig} After Quality control"),
        glue::glue("数据过滤后的 QC 图|||{.seurat_qc_note}") #__REVISE__ set_lab_legend 2026-03-23_21:59:08
      )
      p.qc_pre <- set_lab_legend(
        x$p.qc_pre,
        glue::glue("{x@sig} before Quality control"),
        glue::glue("质量控制 (QC) 图 (数据过滤前) |||{.seurat_qc_note}") #__REVISE__ set_lab_legend 2026-03-23_21:51:34
      )
      x <- plotsAdd(x, p.qc_pre = p.qc_pre, p.qc_aft = p.qc_aft)
      x <- methodAdd(
        x, "前期质量控制{aref(p.qc_pre)}，一个细胞至少应有 {min.features} 个基因，并且基因数量小于 {max.features}。线粒体基因的比例小于 {max.percent.mt}%。保留总基因表达量小于 {max.count} 细胞。过滤前，所有样本共包含 {ncell} 个细胞，{ngene} 个基因。过滤后{aref(p.qc_aft)}，⟦mark$red('所有样本共包含{ncol(object(x))}个细胞，{nrow(object(x))} 个基因用于后续分析。')⟧" #__REVISE__ methodAdd 2026-03-23_22:06:48
      )
      # x <- methodAdd(x, "一个细胞至少应有 {min.features} 个基因，并且基因数量小于 {max.features}。线粒体基因的比例小于 {max.percent.mt}%。根据上述条件，获得用于下游分析的高质量细胞。")
    }
    return(x)
  })

setMethod("step2", signature = c(x = "job_seurat5n"),
  function(x, ndims = 20, sct = FALSE, jk = FALSE, workers = 5L)
  {
    step_message("Run standard anlaysis workflow or `SCTransform`.")
    if (is.remote(x)) {
      if (is.null(workers)) {
        stop('is.null(workers).')
      }
      x <- run_job_remote(x, wait = 1,
        {
          x <- step2(x, ndims = "{ndims}", sct = "{sct}", workers = "{workers}")
        }
      )
      return(x)
    }
    if (sct) {
      if (!is.null(workers)) {
        old_plan <- future::plan()
        on.exit(future::plan(old_plan), add = TRUE)
        if (parallelly::supportsMulticore()) {
          future::plan(future::multicore, workers = workers)
        } else {
          future::plan(future::sequential)
          message(
            "Multicore is not supported in this R session. Run by Rscript for reliable SCTransform parallelization."
          )
        }
        message(glue::glue("future workers: {future::nbrOfWorkers()}"))
      }
      object(x) <- e(Seurat::SCTransform(
          object(x),
          method = "glmGamPoi",
          vars.to.regress = "percent.mt",
          verbose = TRUE,
          assay = SeuratObject::DefaultAssay(object(x))
          ))
      message(glue::glue("Shift assays to {SeuratObject::DefaultAssay(object(x))}"))
      x <- methodAdd(
        x, "使用 `Seurat::SCTransform` (默认参数) 对数据集归一化 (<https://satijalab.org/seurat/articles/sctransform_vignette>) 。" #__REVISE__ methodAdd 2026-03-23_22:41:41
      )
    } else {
      object(x) <- e(Seurat::NormalizeData(object(x)))
      object(x) <- e(Seurat::FindVariableFeatures(object(x)))
      object(x) <- e(Seurat::ScaleData(object(x)))
      x <- methodAdd(
        x, "执行标准 Seurat 分析工作流 (`NormalizeData`, `FindVariableFeatures`, `ScaleData`)。"
      )
      p.varfeature <- e(Seurat::VariableFeaturePlot(object(x)))
      p.varfeature <- set_lab_legend(
        wrap(p.varfeature),
        glue::glue("{x@sig} Variable Feature Plot"), #__REVISE__ set_lab_legend 2026-03-23_22:13:35
        glue::glue("高变基因图|||红色代表高变基因，横坐标为基因在所有细胞中的表达水平（log10对数值），纵坐标为基因在所有细胞中的表达水平的标准差，数值越大，表示该基因在细胞中的表达水平越不稳定。生物学差异（如细胞类型、状态等差异）通常会导致某些基因在不同细胞之间表现出较大变异，因此更有可能提供关于生物学现象的信息。") #__REVISE__ set_lab_legend 2026-03-23_22:15:27
      )
      x <- plotsAdd(x, p.varfeature)
    }
    object(x) <- e(Seurat::RunPCA(object(x)))
    x <- methodAdd(x, "随后 PCA 聚类 (`RunPCA`)。")
    if (jk && !sct) {
      object(x) <- e(Seurat::JackStraw(object(x), dims = ndims))
      object(x) <- e(Seurat::ScoreJackStraw(object(x)))
      p.jackPlot <- Seurat::JackStrawPlot(object(x))
      p.jackPlot <- set_lab_legend(
        wrap(p.jackPlot),
        glue::glue("{x@sig} Jack Straw plot"),
        glue::glue("Jackstraw 置换检验 |||通过对原始数据进行多次置换，构建一个零假设分布，然后将实际观测到的主成分得分与该零假设分布进行比较。每个点表示基因在某个主成分上的投影得分与随机背景的比较，大于或等于实际观测主成分得分的比例就是 p 值。p &lt; 0.05 通常认为在该显著性水平下，实际观测到的主成分得分显著高于随机情况下的得分，说明该主成分具有统计学意义，不是由随机因素导致的。通过量化主成分的显著性强度，与均匀分布（虚线）比较，判断哪些主成分更具有统计学意义，富含低p值基因较多的主成分更有统计学意义。")
      )
      x <- plotsAdd(x, p.jackPlot)
      x <- methodAdd(x, "通过 Jackstraw 函数置换检验重新聚类以检验 PC 的选择结果{aref(p.jackPlot)}（P &lt; 0.05）。")
    }
    p.pca_rank <- e(Seurat::ElbowPlot(object(x), ndims))
    # add Seurat::PCAPlot
    p.pca_rank <- set_lab_legend(
      wrap(pretty_elbowplot(p.pca_rank), 4, 4),
      glue::glue("{x@sig} Standard deviations of PCs"),
      glue::glue("主成分 (PC) 的标准化方差 (Standard deviations)|||横坐标为主成分数目，纵坐标代表基于每个主成分对方差解释率的排名（每个主成分的解释方差是其特征值（eigenvalue），表示它解释了总变异的比例），图中每个点表示一个主成分的方差解释比例。") #__REVISE__ set_lab_legend 2026-03-23_22:01:31
    )
    x <- methodAdd(
      x, "使用 ElbowPlot 函数绘制肘图{aref(p.pca_rank)}，帮助确定用于下游分析的主成分以进行后续分析。" #__REVISE__ methodAdd 2026-03-23_22:27:28
    )
    x <- plotsAdd(x, p.pca_rank)
    # x <- snapAdd(x, "数据归一化，PCA 聚类 (Seurat 标准工作流，见方法章节) 后。")
    return(x)
  })

setMethod("step3", signature = c(x = "job_seurat5n"),
  function(x, dims = 1:15, resolution = .2,
    use = c("HarmonyIntegration", "CCAIntegration", "RPCAIntegration"), ...)
  {
    step_message("Identify clusters of cells")
    use <- match.arg(use)
    if (!is.null(x$JoinLayers) && x$JoinLayers) {
      message("Job is 'job_seurat5n', but 'JoinLayers' has been performed.")
      object(x) <- e(Seurat::FindNeighbors(object(x), dims = dims, reduction = use))
      object(x) <- e(
        Seurat::FindClusters(object(x), resolution = resolution, ...)
      )
      object(x) <- e(
        Seurat::RunUMAP(
          object(x), dims = dims, reduction = use,
          n.neighbors = 50L,
          min.dist = 0.45,
          ...
        )
      )
    } else {
      if (is.null(x$.before_IntegrateLayers)) {
        object(x) <- e(Seurat::FindNeighbors(object(x), dims = dims, reduction = "pca"))
        object(x) <- e(Seurat::FindClusters(object(x), resolution = resolution,
            cluster.name = "unintegrated_clusters"))
        object(x) <- e(Seurat::RunUMAP(object(x), dims = dims,
            reduction = "pca", reduction.name = "umap_unintegrated",
            n.neighbors = 50L,
            min.dist = 0.45
            ))
        x$.before_IntegrateLayers <- TRUE
      }
      p.umapUint <-  e(Seurat::DimPlot(object(x), reduction = "umap_unintegrated",
          group.by = c("orig.ident", "unintegrated_clusters"), cols = color_set(TRUE)))
      p.umapUint <- set_lab_legend(
        wrap(p.umapUint, 10, 5),
        glue::glue("{x@sig} UMAP Unintegrated"),
        glue::glue("去除批次效应之前的 UMAP 聚类图|||不同颜色代表不同cluster。横纵坐标是 UMAP 降维的两个维度。UMAP能够将高维空间中的数据映射到低维空间中，并保留数据集的局部特性。") #__REVISE__ set_lab_legend 2026-03-23_22:44:03
      )
      x <- plotsAdd(x, p.umapUint)
      ## integrated
      methods <- list(CCAIntegration = Seurat::CCAIntegration,
        HarmonyIntegration = Seurat::HarmonyIntegration,
        RPCAIntegration = Seurat::RPCAIntegration
      )
      use <- match.arg(use, names(methods))
      object <- object(x)
      res <- try(e(Seurat::IntegrateLayers(object = object,
            method = methods[[ use ]], orig.reduction = "pca",
            new.reduction = use, verbose = FALSE,
            normalization.method = if (object@active.assay == "SCT") "SCT" else "LogNormalize")))
      if (!inherits(res, "try-error")) {
        object(x) <- res
      } else {
        warning("Got error while perform `Seurat::IntegrateLayers`, return the job.")
        return(x)
      }
      object(x)[["RNA"]] <- e(SeuratObject::JoinLayers(object(x)[["RNA"]]))
      ## SeuratObject::DefaultDimReduc, search in case of UMAP
      object(x)@reductions$umap_unintegrated <- NULL
      ## method of job_seurat
      x <- callNextMethod(
        x, dims, resolution, reduction = use, ...
      )
      x <- methodAdd(
        x, "结果显示{aref(x@plots$step2$p.pca_rank)}，前 {max(dims)} 个 PCs 以后方差增量减缓逐渐趋于稳定，选择前 {max(dims)} 个 PCs 进行后续聚类分析。", add = FALSE
      )
      x <- methodAdd(
        x, "以 `Seurat::IntegrateLayers` 集成数据，去除批次效应 (使用 {use} 方法)。"
      )
    }
    p.umapInt <-  e(Seurat::DimPlot(object(x),
        group.by = c("orig.ident", "seurat_clusters"), cols = color_set(TRUE)))
    p.umapInt <- set_lab_legend(
      wrap(p.umapInt, 10, 5),
      glue::glue("{x@sig} UMAP Integrated"), #__REVISE__ set_lab_legend 2026-03-23_22:45:48
      glue::glue("去除批次效应之后的 UMAP 聚类图|||不同颜色代表不同cluster。横纵坐标是 UMAP 降维的两个维度。UMAP能够将高维空间中的数据映射到低维空间中，并保留数据集的局部特性。")
    )
    p.umapLabel <-  e(Seurat::DimPlot(object(x),
        group.by = c("seurat_clusters"), 
        cols = color_set(TRUE), label = TRUE))
    p.umapLabel <- set_lab_legend(
      p.umapLabel,
      glue::glue("{x@sig} UMAP with label"),
      glue::glue("UMAP 聚类图|||UMAP 图中带有数字注释了细胞簇属于哪个聚类，有利于分辨。不同颜色代表不同cluster。横纵坐标是 UMAP 降维的两个维度。UMAP能够将高维空间中的数据映射到低维空间中，并保留数据集的局部特性。")
    )
    x$checks$ps.checks <- .plot_check_clustering(object(x))
    x <- plotsAdd(x, p.umapInt, p.umapLabel)
    x <- methodAdd(x, "在 1-{max(dims)} PC 维度下，以 `Seurat::FindNeighbors` 构建 Nearest-neighbor Graph。随后在 {resolution} 分辨率下，以 `Seurat::FindClusters` 函数识别细胞群并以 `Seurat::RunUMAP` 进行 UMAP 聚类。")
    nBefore <- length(levels(object(x)@meta.data$unintegrated_clusters))
    nAfter <- length(levels(object(x)@meta.data$seurat_clusters))
    x <- methodAdd(x, "在去除批次效应前，UMAP 图{aref(p.umapUint)}中各样本保持离散。`Seurat::FindClusters` 共找到 {nBefore} 个细胞簇。")
    x <- methodAdd(x, "去除批次效应后{aref(p.umapInt)}，`Seurat::FindClusters` 找到 {nAfter} 个细胞簇，且各样本相互均匀混合，即批次效应已被良好地处理。选择去除批次效应后的数据集进行后续分析。")
    x$JoinLayers <- TRUE
    return(x)
  })

.plot_check_clustering <- function(seu) {
  p.checks_feature <- Seurat::FeaturePlot(
    object = seu,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    reduction = "umap",
    order = TRUE
  )
  p.checks_dim <- Seurat::DimPlot(
    object = seu,
    group.by = "orig.ident",
    reduction = "umap"
  )
  list(p.checks_feature = wrap(p.checks_feature), p.checks_dim = wrap(p.checks_dim))
}

setMethod("asjob_limma", signature = c(x = "job_seurat"),
  function(x, features, cells, cell_groups, group.by = x$group.by,
    slot = "data", fun_norm = function(x) log2(x + 1), gname = TRUE)
  {
    metadata <- dplyr::mutate(
      as_tibble(object(x)@meta.data), 
      sample = rownames, group = !!rlang::sym(group.by), .before = 1
    )
    if (missing(cells) && !missing(cell_groups)) {
      cells <- metadata[[ group.by ]] %in% cell_groups
    }
    metadata <- metadata[ cells, ]
    assay <- object(x)@assays[[ object(x)@active.assay ]]
    data <- SeuratObject::LayerData(object = assay, layer = slot)
    genes <- data.frame(gene = rownames(data))
    if (gname) {
      genes <- dplyr::mutate(genes, gene = gname(gene))
    }
    if (is.character(features)) {
      features <- genes$gene %in% features
    }
    genes <- genes[ features, , drop = FALSE]
    data <- data[ features, cells ]
    data <- data.frame(data, check.names = FALSE)
    data <- fun_norm(data)
    object <- new_from_package(
      "EList", "limma", list(E = data, targets = metadata, genes = genes)
    )
    validObject(object)
    x <- .job_limma()
    x$normed_data <- object
    x$from_seurat <- TRUE
    return(x)
  })


scFuns <- new.env(parent = emptyenv())

scFuns$is_10x_matrix_dir <- function(dir_path)
{
  if (!dir.exists(dir_path)) {
    return(FALSE)
  }

  has_barcode <- any(file.exists(file.path(
    dir_path,
    c("barcodes.tsv", "barcodes.tsv.gz")
  )))

  has_matrix <- any(file.exists(file.path(
    dir_path,
    c("matrix.mtx", "matrix.mtx.gz")
  )))

  has_feature <- any(file.exists(file.path(
    dir_path,
    c("genes.tsv", "genes.tsv.gz", "features.tsv", "features.tsv.gz")
  )))

  has_barcode && has_matrix && has_feature
}


scFuns$find_10x_matrix_dirs <- function(dir_root)
{
  vec_dir <- list.dirs(dir_root, recursive = TRUE, full.names = TRUE)
  vec_keep <- vapply(vec_dir, scFuns$is_10x_matrix_dir, logical(1L))

  vec_dir[vec_keep]
}

scFuns$check_tar_md5 <- function(file_tar)
{
  file_md5 <- paste0(file_tar, ".md5")

  if (!file.exists(file_md5)) {
    return(NA)
  }

  vec_line <- readLines(file_md5, warn = FALSE)

  if (length(vec_line) == 0L) {
    return(NA)
  }

  chr_expected <- regmatches(
    vec_line[[1L]],
    regexpr("[0-9a-fA-F]{32}", vec_line[[1L]])
  )

  if (length(chr_expected) == 0L || is.na(chr_expected) || !nzchar(chr_expected)) {
    return(NA)
  }

  chr_observed <- unname(tools::md5sum(file_tar))

  identical(tolower(chr_observed), tolower(chr_expected))
}

scFuns$get_10x_matrix_dirs <- function(
  dir_root,
  extract_tar = TRUE,
  check_md5 = FALSE,
  verbose = TRUE
)
{
  if (!dir.exists(dir_root)) {
    stop(glue::glue("Directory does not exist: {dir_root}"), call. = FALSE)
  }

  dir_root <- normalizePath(dir_root, winslash = "/", mustWork = TRUE)

  vec_sample_dir <- list.dirs(dir_root, recursive = FALSE, full.names = TRUE)
  vec_empty_dir <- vec_sample_dir[vapply(
    vec_sample_dir,
    function(x) length(list.files(x, all.files = TRUE, no.. = TRUE)) == 0L,
    logical(1L)
  )]

  vec_tar <- list.files(
    dir_root,
    pattern = "\\.tar$",
    recursive = TRUE,
    full.names = TRUE
  )

  data_matrix_before <- scFuns$make_10x_matrix_table(dir_root)

  n_already <- 0L
  n_extracted <- 0L
  n_failed <- 0L
  n_md5_failed <- 0L

  if (verbose) {
    message(glue::glue("Root directory: {dir_root}"))
    message(glue::glue("Sample directories: {length(vec_sample_dir)}"))
    message(glue::glue("Empty sample directories: {length(vec_empty_dir)}"))
    message(glue::glue("Tar files found: {length(vec_tar)}"))
    message(glue::glue(
      "Valid 10x matrix directories before extraction: {sum(data_matrix_before$is_valid)}"
    ))
  }

  if (verbose && length(vec_empty_dir) > 0L) {
    message(glue::glue(
      "Empty directories: {paste(basename(vec_empty_dir), collapse = ', ')}"
    ))
  }

  for (file_tar in vec_tar) {
    vec_candidate <- scFuns$get_tar_matrix_candidates(file_tar)
    vec_existing_matrix <- scFuns$get_existing_matrix_from_candidates(vec_candidate)

    if (length(vec_existing_matrix) > 0L) {
      n_already <- n_already + 1L

      if (verbose) {
        message(glue::glue("Already extracted: {basename(file_tar)}"))
      }

      next
    }

    if (!extract_tar) {
      next
    }

    if (check_md5) {
      val_md5 <- scFuns$check_tar_md5(file_tar)

      if (identical(val_md5, FALSE)) {
        n_md5_failed <- n_md5_failed + 1L
        n_failed <- n_failed + 1L

        if (verbose) {
          message(glue::glue("Skip tar due to MD5 mismatch: {basename(file_tar)}"))
        }

        next
      }

      if (is.na(val_md5) && verbose) {
        message(glue::glue("MD5 file unavailable or unreadable: {basename(file_tar)}"))
      }
    }

    ok_extract <- tryCatch(
      {
        withCallingHandlers(
          utils::untar(file_tar, exdir = dirname(file_tar)),
          warning = function(w) {
            if (verbose) {
              message(glue::glue(
                "Warning while extracting {basename(file_tar)}: {conditionMessage(w)}"
              ))
            }

            invokeRestart("muffleWarning")
          }
        )

        TRUE
      },
      error = function(e) {
        if (verbose) {
          message(glue::glue(
            "Failed to extract {basename(file_tar)}: {conditionMessage(e)}"
          ))
        }

        FALSE
      }
    )

    vec_existing_matrix <- scFuns$get_existing_matrix_from_candidates(vec_candidate)

    if (ok_extract && length(vec_existing_matrix) > 0L) {
      n_extracted <- n_extracted + 1L

      if (verbose) {
        message(glue::glue("Extracted: {basename(file_tar)}"))
      }
    } else {
      n_failed <- n_failed + 1L

      if (verbose) {
        message(glue::glue(
          "No valid 10x matrix directory found after extraction: {basename(file_tar)}"
        ))
      }
    }
  }

  data_matrix <- scFuns$make_10x_matrix_table(dir_root)

  n_valid_sample <- length(unique(data_matrix$sample_id[data_matrix$is_valid]))
  n_missing_sample <- length(unique(data_matrix$sample_id[!data_matrix$is_valid]))

  if (verbose) {
    message(glue::glue("Already extracted tar files: {n_already}"))
    message(glue::glue("Newly extracted tar files: {n_extracted}"))
    message(glue::glue("Failed tar files: {n_failed}"))
    message(glue::glue("MD5 failed tar files: {n_md5_failed}"))
    message(glue::glue("Valid 10x matrix directories after extraction: {sum(data_matrix$is_valid)}"))
    message(glue::glue("Samples with valid matrix: {n_valid_sample}"))
    message(glue::glue("Samples without valid matrix: {n_missing_sample}"))
  }

  if (verbose && n_missing_sample > 0L) {
    vec_missing_sample <- unique(data_matrix$sample_id[!data_matrix$is_valid])

    message(glue::glue(
      "Samples without valid matrix: {paste(vec_missing_sample, collapse = ', ')}"
    ))
  }

  attr(data_matrix, "summary") <- list(
    n_sample_dir = length(vec_sample_dir),
    n_empty_sample_dir = length(vec_empty_dir),
    n_tar = length(vec_tar),
    n_matrix_before = sum(data_matrix_before$is_valid),
    n_already = n_already,
    n_extracted = n_extracted,
    n_failed = n_failed,
    n_md5_failed = n_md5_failed,
    n_matrix_after = sum(data_matrix$is_valid),
    n_valid_sample = n_valid_sample,
    n_missing_sample = n_missing_sample,
    empty_sample_dir = basename(vec_empty_dir),
    missing_sample = unique(data_matrix$sample_id[!data_matrix$is_valid])
  )

  data_matrix
}

scFuns$get_tar_matrix_candidates <- function(file_tar)
{
  dir_tar <- dirname(file_tar)
  dir_expected <- sub("\\.tar$", "", file_tar)

  vec_entry <- tryCatch(
    utils::untar(file_tar, list = TRUE),
    error = function(e) {
      message(glue::glue(
        "Failed to list tar entries: {basename(file_tar)}; {conditionMessage(e)}"
      ))

      character(0L)
    }
  )

  if (length(vec_entry) == 0L) {
    return(unique(c(dir_expected, dir_tar)))
  }

  vec_entry <- gsub("^\\./", "", vec_entry)
  vec_entry <- gsub("^/+", "", vec_entry)
  vec_entry <- vec_entry[!is.na(vec_entry)]
  vec_entry <- vec_entry[nzchar(vec_entry)]

  if (length(vec_entry) == 0L) {
    return(unique(c(dir_expected, dir_tar)))
  }

  vec_split <- strsplit(vec_entry, "/", fixed = TRUE)

  vec_top <- vapply(
    vec_split,
    function(x) {
      if (length(x) == 0L || is.na(x[[1L]]) || !nzchar(x[[1L]])) {
        return(NA_character_)
      }

      x[[1L]]
    },
    character(1L)
  )

  vec_top <- unique(vec_top[!is.na(vec_top) & nzchar(vec_top)])
  vec_candidate <- file.path(dir_tar, vec_top)

  unique(c(dir_expected, vec_candidate, dir_tar))
}


scFuns$get_existing_matrix_from_candidates <- function(vec_candidate)
{
  vec_candidate <- vec_candidate[dir.exists(vec_candidate)]

  if (length(vec_candidate) == 0L) {
    return(character(0L))
  }

  vec_matrix <- unlist(
    lapply(
      vec_candidate,
      function(dir_path) {
        vec_found <- scFuns$find_10x_matrix_dirs(dir_path)

        if (scFuns$is_10x_matrix_dir(dir_path)) {
          vec_found <- unique(c(dir_path, vec_found))
        }

        vec_found
      }
    ),
    use.names = FALSE
  )

  unique(vec_matrix)
}


scFuns$make_10x_matrix_table <- function(dir_root)
{
  dir_root <- normalizePath(dir_root, winslash = "/", mustWork = TRUE)

  vec_sample_dir <- list.dirs(dir_root, recursive = FALSE, full.names = TRUE)

  if (length(vec_sample_dir) == 0L) {
    return(data.frame(
      sample_id = character(0L),
      sample_dir = character(0L),
      matrix_dir = character(0L),
      matrix_name = character(0L),
      matrix_rel_path = character(0L),
      n_matrix = integer(0L),
      n_tar = integer(0L),
      tar_file = character(0L),
      is_valid = logical(0L),
      is_empty_sample = logical(0L),
      status = character(0L),
      stringsAsFactors = FALSE
    ))
  }

  lst_table <- lapply(
    vec_sample_dir,
    function(dir_sample) {
      sample_id <- basename(dir_sample)

      vec_matrix_dir <- scFuns$find_10x_matrix_dirs(dir_sample)
      vec_matrix_dir <- normalizePath(
        vec_matrix_dir,
        winslash = "/",
        mustWork = TRUE
      )

      vec_tar <- list.files(
        dir_sample,
        pattern = "\\.tar$",
        recursive = TRUE,
        full.names = TRUE
      )

      is_empty_sample <- length(list.files(
        dir_sample,
        all.files = TRUE,
        no.. = TRUE
      )) == 0L

      tar_file <- paste(basename(vec_tar), collapse = "; ")

      if (length(vec_matrix_dir) == 0L) {
        status <- if (is_empty_sample) {
          "empty_sample_dir"
        } else if (length(vec_tar) > 0L) {
          "tar_without_valid_matrix"
        } else {
          "no_tar_no_matrix"
        }

        return(data.frame(
          sample_id = sample_id,
          sample_dir = normalizePath(dir_sample, winslash = "/", mustWork = TRUE),
          matrix_dir = NA_character_,
          matrix_name = NA_character_,
          matrix_rel_path = NA_character_,
          n_matrix = 0L,
          n_tar = length(vec_tar),
          tar_file = tar_file,
          is_valid = FALSE,
          is_empty_sample = is_empty_sample,
          status = status,
          stringsAsFactors = FALSE
        ))
      }

      vec_matrix_rel_path <- substring(
        vec_matrix_dir,
        nchar(normalizePath(dir_sample, winslash = "/", mustWork = TRUE)) + 2L
      )

      data.frame(
        sample_id = sample_id,
        sample_dir = normalizePath(dir_sample, winslash = "/", mustWork = TRUE),
        matrix_dir = vec_matrix_dir,
        matrix_name = basename(vec_matrix_dir),
        matrix_rel_path = vec_matrix_rel_path,
        n_matrix = length(vec_matrix_dir),
        n_tar = length(vec_tar),
        tar_file = tar_file,
        is_valid = TRUE,
        is_empty_sample = is_empty_sample,
        status = "valid_matrix",
        stringsAsFactors = FALSE
      )
    }
  )

  data_table <- do.call(rbind, lst_table)
  rownames(data_table) <- NULL

  data_table
}

# ==========================================================================

seuFuns <- new.env(parent = parent.frame())

seuFuns$.resolve_integrated_reduction_seurat <- function(object, use = NULL)
{
  vec_reduction <- Seurat::Reductions(object)

  if (!is.null(use)) {
    vec_use <- as.character(use)
    vec_hit <- vec_use[vec_use %in% vec_reduction]

    if (length(vec_hit) > 0L) {
      return(vec_hit[[1L]])
    }

    stop(glue::glue(
      "None of the requested integrated reductions exists: {paste(vec_use, collapse = ', ')}."
    ))
  }

  vec_candidate <- c(
    "HarmonyIntegration",
    "CCAIntegration",
    "RPCAIntegration",
    "integrated.dr",
    "harmony",
    "integrated"
  )
  vec_hit <- vec_candidate[vec_candidate %in% vec_reduction]

  if (length(vec_hit) == 0L) {
    stop(glue::glue(
      "Cannot find integrated reduction. Available reductions: {paste(vec_reduction, collapse = ', ')}."
    ))
  }

  return(vec_hit[[1L]])
}

seuFuns$.check_umap_dims_seurat <- function(object, reduction, dims)
{
  n_dim <- ncol(Seurat::Embeddings(object, reduction = reduction))

  if (max(dims) > n_dim) {
    stop(glue::glue(
      "Requested dims exceed available dimensions in '{reduction}': max(dims) = {max(dims)}, available = {n_dim}."
    ))
  }

  return(invisible(TRUE))
}

seuFuns$.backup_reduction_seurat <- function(object, reduction, backup_prefix = "backup")
{
  if (!reduction %in% Seurat::Reductions(object)) {
    return(object)
  }

  time_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")
  backup_name <- paste0(reduction, "_", backup_prefix, "_", time_tag)

  object@reductions[[backup_name]] <- object@reductions[[reduction]]
  message(glue::glue("Backup reduction '{reduction}' to '{backup_name}'."))

  return(object)
}

seuFuns$.drop_reduction_seurat <- function(object, reduction)
{
  if (reduction %in% Seurat::Reductions(object)) {
    object@reductions[[reduction]] <- NULL
  }

  return(object)
}

redraw_umap <- function(x, dims = 1L:15L, resolution = .2,
  use = "HarmonyIntegration", n_neighbors = 20L, min_dist = 0.08, spread = 1.5,
  seed_use = 2026L, replace = FALSE, replace_reduction = replace,
  backup = FALSE, suffix = "redraw", ...)
{
  if (!methods::is(x, "job_seurat5n")) {
    stop('!is(x, "job_seurat5n").')
  }

  dims <- as.integer(dims)
  object_now <- object(x)

  if (!"pca" %in% Seurat::Reductions(object_now)) {
    stop("Cannot redraw unintegrated UMAP because reduction 'pca' is missing.")
  }

  integrated_reduction <- seuFuns$.resolve_integrated_reduction_seurat(
    object = object_now,
    use = use
  )

  seuFuns$.check_umap_dims_seurat(object_now, "pca", dims)
  seuFuns$.check_umap_dims_seurat(object_now, integrated_reduction, dims)

  if (replace_reduction) {
    name_umap_unintegrated <- "umap_unintegrated"
    name_umap_integrated <- "umap"
  } else {
    name_umap_unintegrated <- paste0("umap_unintegrated_", suffix)
    name_umap_integrated <- paste0("umap_integrated_", suffix)
  }

  if (replace_reduction && backup) {
    object_now <- seuFuns$.backup_reduction_seurat(object_now, name_umap_unintegrated)
    object_now <- seuFuns$.backup_reduction_seurat(object_now, name_umap_integrated)
  }

  object_now <- seuFuns$.drop_reduction_seurat(object_now, name_umap_unintegrated)
  object_now <- seuFuns$.drop_reduction_seurat(object_now, name_umap_integrated)

  message(glue::glue(
    "Redraw unintegrated UMAP: pca, dims = {min(dims)}-{max(dims)}, n.neighbors = {n_neighbors}, min.dist = {min_dist}, spread = {spread}."
  ))

  object_now <- Seurat::RunUMAP(
    object_now,
    reduction = "pca",
    dims = dims,
    reduction.name = name_umap_unintegrated,
    reduction.key = "UMAPU_",
    n.neighbors = n_neighbors,
    min.dist = min_dist,
    spread = spread,
    seed.use = seed_use,
    verbose = TRUE,
    ...
  )

  message(glue::glue(
    "Redraw integrated UMAP: {integrated_reduction}, dims = {min(dims)}-{max(dims)}, n.neighbors = {n_neighbors}, min.dist = {min_dist}, spread = {spread}."
  ))

  object_now <- Seurat::RunUMAP(
    object_now,
    reduction = integrated_reduction,
    dims = dims,
    reduction.name = name_umap_integrated,
    reduction.key = "UMAP_",
    n.neighbors = n_neighbors,
    min.dist = min_dist,
    spread = spread,
    seed.use = seed_use,
    verbose = TRUE,
    ...
  )

  object(x) <- object_now

  vec_meta <- colnames(object(x)@meta.data)
  cluster_unintegrated <- if ("unintegrated_clusters" %in% vec_meta) {
    "unintegrated_clusters"
  } else {
    "seurat_clusters"
  }

  p.umapUint <- Seurat::DimPlot(
    object(x),
    reduction = name_umap_unintegrated,
    group.by = c("orig.ident", cluster_unintegrated),
    cols = color_set(TRUE),
    raster = TRUE
  )

  p.umapUint <- set_lab_legend(
    wrap(p.umapUint, 10L, 5L),
    glue::glue("{x@sig} UMAP Unintegrated"),
    glue::glue("去除批次效应之前的 UMAP 聚类图|||不同颜色代表不同样本或整合前 cluster。横纵坐标是 UMAP 降维的两个维度。UMAP 能够将高维空间中的数据映射到低维空间中，并保留数据集的局部特性。")
  )

  p.umapInt <- Seurat::DimPlot(
    object(x),
    reduction = name_umap_integrated,
    group.by = c("orig.ident", "seurat_clusters"),
    cols = color_set(TRUE),
    raster = TRUE
  )

  p.umapInt <- set_lab_legend(
    wrap(p.umapInt, 10L, 5L),
    glue::glue("{x@sig} UMAP Integrated"),
    glue::glue("去除批次效应之后的 UMAP 聚类图|||不同颜色代表不同样本或整合后 cluster。横纵坐标是 UMAP 降维的两个维度。UMAP 能够将高维空间中的数据映射到低维空间中，并保留数据集的局部特性。")
  )

  p.umapLabel <- Seurat::DimPlot(
    object(x),
    reduction = name_umap_integrated,
    group.by = "seurat_clusters",
    cols = color_set(TRUE),
    label = TRUE,
    repel = TRUE,
    raster = TRUE
  )

  p.umapLabel <- set_lab_legend(
    p.umapLabel,
    glue::glue("{x@sig} UMAP with label"),
    glue::glue("UMAP 聚类图|||UMAP 图中带有数字注释了细胞簇属于哪个聚类，有利于分辨。不同颜色代表不同 cluster。横纵坐标是 UMAP 降维的两个维度。UMAP 能够将高维空间中的数据映射到低维空间中，并保留数据集的局部特性。")
  )

  if (replace) {
    if (is.null(x$plots$step3)) {
      x$plots$step3 <- list()
    }

    x$plots$step3$p.umapUint <- p.umapUint
    x$plots$step3$p.umapInt <- p.umapInt
    x$plots$step3$p.umapLabel <- p.umapLabel
  } else {
    x <- plotsAdd(
      x, p.umapUint_new, p.umapInt_new, p.umapLabel_new, step = 3L
    )
  }

  if (exists(".plot_check_clustering", mode = "function", inherits = TRUE)) {
    x$checks$ps.checks <- .plot_check_clustering(object(x))
  }

  x$checks$redraw_umap <- list(
    dims = dims,
    resolution_not_used = resolution,
    integrated_reduction = integrated_reduction,
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    spread = spread,
    seed_use = seed_use,
    replace = replace,
    replace_reduction = replace_reduction,
    umap_unintegrated = name_umap_unintegrated,
    umap_integrated = name_umap_integrated
  )

  return(x)
}


# ==========================================================================


.seurat_qc_note <- "上方小提琴图：每个样本对应一个‘小提琴’，小提琴的宽度代表相应数据的密度，宽度越大表示在该区域内的数据点越密集，更多数据点集中于此区域；宽度越小则表示密度越小，即数据相对较少；过滤标准：nCount 和 nFeature 过高可能是双细胞，过低可能是细胞碎片；percent.mt（线粒体基因表达比例，是细胞内线粒体基因表达量占所有基因表达量的比例）表明细胞状态，值过高可能是细胞正在经历压力或死亡。下方点图：每一个点代表一个细胞，不同颜色代表不同样本；左图横坐标为总基因表达数，纵坐标为线粒体基因比例；右图横坐标为 nCount（总基因表达数），纵坐标为 Feature（总基因数）；正常情况下，nCount 越多那么 nFeature 就越高，呈现出正相关关系，因此检测到的基因表达数应与检测到的基因数目在细胞间高度相关，而线粒体基因比例则不相关（若呈正相关，横坐标越大纵坐标也越大；若呈负相关，横坐标越大纵坐标应越小）。"

