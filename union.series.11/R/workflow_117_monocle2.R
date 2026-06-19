# ==========================================================================
# workflow of monocle2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_monocle2 <- setClass("job_monocle2", 
  contains = c("job"),
  prototype = prototype(
    pg = "monocle2",
    info = c("https://cole-trapnell-lab.github.io/monocle-release/docs/"),
    cite = "",
    method = "",
    tag = "monocle",
    analysis = "拟时序轨迹分析"
    ))

setGeneric("asjob_monocle2",
  function(x, ...) standardGeneric("asjob_monocle2"))

setMethod("asjob_monocle2", signature = c(x = "job_seurat"),
  function(x, compare = .guess_levels_from_job_seurat(x),
    compare.by = "group", group.by = x$group.by, nfeatures = 1000, min.pct = .1)
  {
    if (!requireNamespace("monocle", quietly = TRUE)) {
      stop('!requireNamespace("monocle").')
    }
    metadata <- object(x)@meta.data
    if (SeuratObject::DefaultAssay(object(x)) == "SCT") {
      message(
        glue::glue(
          "Default Assay is SCT, to found variables, use data in assay 'RNA'"
        )
      )
      object <- object(x)
      SeuratObject::DefaultAssay(object) <- "RNA"
      object <- e(Seurat::NormalizeData(object))
      object <- e(Seurat::FindVariableFeatures(object))
      VariableFeatures <- e(Seurat::VariableFeatures(object))
    } else {
      object(x) <- e(Seurat::FindVariableFeatures(
        object(x), nfeatures = nfeatures
      ))
      VariableFeatures <- e(Seurat::VariableFeatures(object(x)))
    }
    if (!length(VariableFeatures)) {
      stop('!length(VariableFeatures).')
    }
    diff_genes <- NULL
    if (!is.null(compare)) {
      if (length(compare) != 2) {
        stop('length(compare) != 2.')
      }
      Seurat::Idents(object(x)) <- compare.by
      if (SeuratObject::DefaultAssay(object(x)) == "SCT" && !is.null(x$seurat_subset) && x$seurat_subset) {
        SeuratObject::DefaultAssay(object(x)) <- "RNA"
        object(x) <- e(Seurat::NormalizeData(object(x)))
        object(x) <- e(Seurat::ScaleData(object(x)))
      }
      diff_genes <- e(
        Seurat::FindMarkers(object(x), compare[1], compare[2], min.pct = min.pct)
      )
      snap(diff_genes) <- glue::glue("使用 Seurat::FindMarkers 默认参数进行组间比较 {bind(compare, co = ' vs ')}")
    }
    # metadata$Cluster <- object(x)@active.ident
    cells <- unique(metadata[[group.by]])
    snapAdd_onExit("x", "将 {bind(cells)} 进行拟时间轴轨迹分析。")
    counts <- object(x)@assays$RNA$counts
    phenoData <- new('AnnotatedDataFrame', data = metadata)
    phenoData$Size_Factor <- rep(NA_real_, ncol(counts))
    featureData = new(
      'AnnotatedDataFrame',
      data = data.frame(gene_short_name = row.names(counts), row.names = row.names(counts))
    )
    cli::cli_alert_info("new('CellDataSet', ...)")
    # don't use `newCellDataSet`
    object <- new(
      "CellDataSet", assayData = Biobase::assayDataNew("environment", exprs = counts),
      phenoData = phenoData, featureData = featureData, 
      lowerDetectionLimit = .01,
      expressionFamily = VGAM::negbinomial.size(),
      dispFitInfo = new.env(hash = TRUE)
    )
    validObject(object)
    x <- .job_monocle2(object = object)
    x$VariableFeatures <- VariableFeatures
    x$group.by <- group.by
    x$diff_genes <- diff_genes
    x <- methodAdd(x, "基于 **Monocle2** 的单细胞转录组拟时序分析，本研究利用其核心的反转图嵌入（Reversed Graph Embedding）算法，在低维空间中重构所选细胞的连续转录状态轨迹。Monocle2 可用于描述复杂生物过程中的细胞状态转换和分支结构（PMID: 28825705）；拟时序分析也可用于研究不同疾病状态或实验条件相关的连续转录变化，而不必局限于严格的发育分化过程（PMID: 37949861）。因此，本流程将根据 ordering genes 的来源对轨迹进行解释：若使用高变基因，则主要反映所选细胞的无监督转录异质性；若使用分组差异基因，则主要反映分组或疾病状态相关的连续转录状态谱。")
    return(x)
  })

setMethod("step0", signature = c(x = "job_monocle2"),
  function(x){
    step_message("Prepare your data with function `job_monocle2`.")
  })

setMethod("step1", signature = c(x = "job_monocle2"),
  function(x, not_run = FALSE)
  {
    step_message("Detect features.")
    if (!not_run) {
      object(x) <- e(BiocGenerics::estimateSizeFactors(object(x)))
      object(x) <- e(BiocGenerics::estimateDispersions(object(x)))
      object(x) <- e(monocle::detectGenes(object(x), min_expr = 1))
    }
    x <- methodAdd(x, "以 R 包 `monocle` ⟦pkgInfo('monocle')⟧ 对所选细胞进行细胞拟时序轨迹分析。")
    return(x)
  })

setMethod("step2", signature = c(x = "job_monocle2"),
  function(x, mode = c("diff", "var"), top = 300,
    try_sig = TRUE, cut.fc = .5, use.p = c("p_val", "p_val_adj"), 
    group = "group", not_run = FALSE)
  {
    step_message("DDRTree.")
    require(DDRTree)
    use.p <- match.arg(use.p)
    if (!missing(mode) && length(mode) > 1) {
      order.by <- mode
      x$order_mode <- "custom_genes"
      x$order_by <- order.by
      x <- methodAdd(x, "使用外部指定的 {length(order.by)} 个 ordering genes，并通过 `monocle::setOrderingFilter` 对细胞轨迹排序。该模式下拟时序反映指定基因集所代表的生物学过程或细胞状态连续变化，需结合输入基因集的来源和功能进行解释。")
      # message(glue::glue("Input genes will be `head` by `top` number."))
      # order.by <- head(order.by, n = top)
    } else {
      mode <- match.arg(mode)
      if (mode == "var") {
        order.by <- x$VariableFeatures
        x$order_mode <- "variable_features"
        x$order_by <- order.by
        x <- methodAdd(x, "以 `Seurat::FindVariableFeatures` 选取 {length(order.by)} 个高变基因，并使用 `monocle::setOrderingFilter` 对细胞轨迹排序。该模式属于无监督 ordering gene 策略，主要用于展示所选细胞内部由高变表达程序驱动的整体转录异质性和潜在状态连续谱，适合探索细胞亚型或未知状态变化。")
      } else {
        # use variable features, or use DEGs with control vs model?
        # [@An_atlas_of_epi_Han_G_2024] 38418883
        message(glue::glue("Use Top {top} DEGs as ordering principle."))
        if (!is.null(x$diff_genes)) {
          data <- x$diff_genes
          snap_ex <- ""
          if (try_sig) {
            dataSig <- dplyr::filter(
              data, abs(avg_log2FC) > cut.fc, !!rlang::sym(use.p) < .05
            )
            if (nrow(dataSig) < top) {
              top <- nrow(dataSig)
              message(glue::glue("Too less significant genes ({top})..."))
            } else {
              message(glue::glue("Ordering significant genes by `avg_log2FC` (n_top: {top})"))
              dataSig <- dplyr::arrange(
                dataSig, dplyr::desc(abs(avg_log2FC)), !!rlang::sym(use.p)
              )
            }
            snap_ex <- glue::glue("(⟦mark$blue('{detail(use.p)} &lt; 0.05, |avg_log2FC| &gt; {cut.fc}')⟧)")
            data <- dataSig
          }
          order.by <- head(rownames(data), n = top)
          x$order_mode <- "differential_genes"
          x$order_by <- order.by
          if (any(!order.by %in% rownames(object(x)))) {
            stop('any(!order.by %in% rownames(object(x))).')
          }
          x <- methodAdd(x, "参考已发表单细胞轨迹分析中使用差异表达基因作为 ordering genes 的策略{cite_show('An_atlas_of_epi_Han_G_2024')}(PMID: 38418883)，本研究根据分组差异表达基因排序选取 Top {top} ({snap(x$diff_genes)}){snap_ex}，并使用 `monocle::setOrderingFilter` 对细胞轨迹排序。该模式属于分组/疾病状态导向的 ordering gene 策略，更适合解析该组间比较相关转录状态在所选细胞中的连续分布；因此所得拟时序主要解释为分组相关细胞状态连续谱，而不直接等同于严格的发育分化方向。")
        } else {
          stop('!is.null(x$diff_genes).')
        }
      }
    }
    if (!not_run) {
      object(x) <- e(monocle::setOrderingFilter(object(x), ordering_genes = order.by))
      object(x) <- e(monocle::reduceDimension(object(x), reduction_method = "DDRTree"))
      object(x) <- e(monocle::orderCells(object(x)))
    }
    return(x)
  })

setMethod("step3", signature = c(x = "job_monocle2"),
  function(x, use = c("Pseudotime", "State", x$group.by), 
    extra = "group", root = NULL, not_run = FALSE)
  {
    step_message("Plot cell Trajectory")
    use <- c(use, extra)
    if (!is.character(use)) {
      stop('!is.character(use).')
    }
    if (!is.null(root) && !not_run) {
      object(x) <- e(monocle::orderCells(object(x), root_state = root))
    }
    cli::cli_alert_info("monocle::plot_cell_trajectory")
    lst <- pbapply::pbsapply(use, simplify = FALSE,
      function(type) {
        p <- monocle::plot_cell_trajectory(object(x), color_by = type)
        if (type == x$group.by) {
          p + guides(color = guide_legend(nrow = 2))
        } else p
      })
    message(glue::glue("Finished plot cell trajectory."))
    p.traj <- smart_wrap(lst, 5, max_ratio = 2)
    snaps <- c(
      Pseudotime = "细胞拟时间连续状态轨迹图，不同颜色代表细胞在轨迹上的相对进程位置；",
      State = "Monocle2 推断的轨迹状态图，不同颜色代表细胞处于不同轨迹状态或分支区域；",
      cell = "不同细胞类型或细胞亚群的轨迹分布图，不同颜色代表不同细胞注释；",
      group = "不同样本分组的轨迹分布图，不同颜色代表细胞所属分组。"
    )
    snaps <- snaps[ seq_along(use) ]
    snaps <- setNames(
      snaps, c("Pseudotime", "State", x$group.by, extra)
    )
    snaps <- bind(snaps[match(use, names(snaps))], co = "")
    p.traj <- set_lab_legend(
      p.traj,
      glue::glue("{x@sig} cell trajectories"),
      glue::glue(
        "细胞拟时轨迹图。|||横纵坐标分别为拟时序降维后的两个维度，图中每个圆点代表一个细胞，黑色圆圈内的数字代表 Monocle2 推断的轨迹状态节点。Pseudotime 表示细胞在所构建轨迹上的相对进程位置；其生物学方向取决于 root_state 设置及所选 ordering genes，需要结合细胞注释、分组信息和关键基因表达共同解释。从左到右、从上到下，各子图分别为：{snaps}"
      )
    )
    x$use <- use
    x <- plotsAdd(x, p.traj)
    x <- methodAdd(x, "使用 `monocle::plot_cell_trajectory` 函数绘制细胞拟时轨迹图，并分别展示 Pseudotime、State、细胞注释及分组信息在轨迹空间中的分布。需要说明的是，Monocle2 计算的 Pseudotime 是沿重构轨迹的相对距离，若未明确指定或验证 root_state，其方向不应被直接解释为真实时间或严格分化顺序；在差异基因排序模式下，更适合解释为分组或疾病状态相关的连续转录状态变化。")
    return(x)
  })

setMethod("step4", signature = c(x = "job_monocle2"),
  function(x, ref, use = x$use, recode = NULL, ...)
  {
    step_message("Plot genes in pseudotime.")
    set.seed(x$seed)
    if (!is(ref, "feature")) {
      stop('!is(ref, "feature").')
    }
    regenes <- genes <- unique(resolve_feature(ref))
    if (!is.null(recode)) {
      regenes <- dplyr::recode(
        genes, !!!setNames(names(recode), unname(recode))
      )
      fun_recode <- function(data) {
        dplyr::mutate(data,
          # f_id = dplyr::recode(f_id, !!!recode),
          # gene_short_name = dplyr::recode(gene_short_name, !!!recode),
          feature_label = dplyr::recode(feature_label, !!!recode)
        )
      }
      fun_recode_layer <- function(layers) {
        for (i in seq_along(layers)) {
          data <- layers[[i]]$data
          if (!is.null(data) && !is.null(data$feature_label)) {
            layers[[i]]$data <- dplyr::mutate(
              layers[[i]]$data, feature_label = dplyr::recode(feature_label, !!!recode)
            )
          }
        }
        return(layers)
      }
    }
    if (length(regenes) > 10) {
      stop('length(regenes) > 10, too many input.')
    }
    if (any(notGot <- !regenes %in% rownames(object(x)))) {
      stop(glue::glue("Not got: {bind(regenes[notGot])}"))
    }
    object <- object(x)[regenes, ]
    cli::cli_alert_info("monocle::plot_genes_in_pseudotime")
    plot_pseudotime <- function(object, color_by, ...) {
      suppressMessages(require(monocle))
      monocle::plot_genes_in_pseudotime(object, color_by = color_by, ...)
    }
    p.geneInPseudo <- pbapply::pbsapply(use,
      function(type) {
        args <- list(object = object, color_by = type, ...)
        p <- callr::r(plot_pseudotime, args, libpath = .libPaths(), show = TRUE)
        if (!is.null(recode)) {
          p <- .set_ggplot_content(p, fun_recode)
          p <- .set_ggplot_content(p, fun_recode_layer, "layers")
        }
        wrap(p, 5, 1.5 * length(regenes))
      }, simplify = FALSE)
    p.geneInPseudo <- set_lab_legend(
      p.geneInPseudo,
      glue::glue("{x@sig} genes in trajectorie of {use}"),
      glue::glue("基因在细胞轨迹图中的表达量变化|||横坐标为细胞的伪时间排序，纵轴表示基因的表达量，每一个点代表一个细胞，颜色代表图例所示的类型 ({use}) 。")
    )
    x <- plotsAdd(x, p.geneInPseudo)
    x <- methodAdd(x, "使用 `monocle::plot_genes_in_pseudotime` 绘制 {snap(ref)} 在所选细胞拟时序轨迹中的表达水平变化，用于观察基因是否沿连续转录状态呈现动态变化。")
    return(x)
  })

setMethod("step5", signature = c(x = "job_monocle2"),
  function(x, point = NULL, list_branches = NULL,
    maxShow = 50, workers = 1, 
    rerun = FALSE, features = NULL)
  {
    step_message("Pseudotime heatmap.")
    if (!is.null(features)) {
      if (!is(features, "feature")) {
        stop('!is(features, "feature").')
      }
      genes <- resolve_feature(features)
      x <- snapAdd(x, "以 `monocle::plot_genes_branched_heatmap` 函数，绘制拟时序相关基因表达热图，展示关键分支点上 {snap(features)} 表达动态变化。")
    } else {
      run_diff <- function(cds, cores) {
        require(monocle)
        require(VGAM)
        monocle::differentialGeneTest(cds = cds, cores = cores)
      }
      fun_cache <- function(...) {
        genes <- x$VariableFeatures
        callr::r(
          run_diff, list(cds = object(x)[genes, ], cores = workers),
          libpath = .libPaths(), show = FALSE
        )
      }
      args <- x$.args[ names(x$.args) %in% paste0("step", 1:4) ]
      diff_test_pseudotime <- expect_local_data(
        "tmp", "monocle_diff", fun_cache, list(args), rerun = rerun
      )
      x <- methodAdd(x, "以 `monocle::differentialGeneTest` 根据 Pseudotime 鉴定高变基因中 (n = {length(x$VariableFeatures)}) 与拟时序连续状态相关的动态变化基因。")
      x <- snapAdd(x, "根据 Pseudotime 一共鉴定到 {nrow(data)} (⟦mark$blue('qval &lt; 0.05')⟧) 个与拟时序连续状态相关的动态变化基因，并以热图展示{aref(p.hp)} (Top {length(genes)})。")
      diff_test_pseudotime <- tibble::as_tibble(diff_test_pseudotime)
      x$diff_test_pseudotime <- diff_test_pseudotime <- dplyr::arrange(diff_test_pseudotime, qval)
      data <- dplyr::filter(diff_test_pseudotime, qval < .05)
      genes <- head(data$gene_short_name, n = maxShow)
    }
    if (length(genes) < 4L) {
      n_cluster <- 1L
    } else {
      n_cluster <- 4L
    }
    fun_heatmap <- function(cds, point, n_cluster, list_branches) {
      require(monocle)
      if (!is.null(point)) {
        if (is.null(list_branches)) {
          monocle::plot_genes_branched_heatmap(
            cds, branch_point = point, num_clusters = n_cluster,
            show_rownames = TRUE, return_heatmap = TRUE
          )
        } else {
          lapply(list_branches,
            function(branch_states) {
              if (length(branch_states) != 2) {
                stop('length(branch_states) != 2.')
              }
              branch_labels <- paste0("Cell state ", branch_states)
              message(
                glue::glue("Compare cell state: {paste(branch_states, collapse = ' vs ')}")
              )
              monocle::plot_genes_branched_heatmap(
                cds, branch_point = point,
                branch_states = branch_states,
                branch_labels = branch_labels,
                num_clusters = n_cluster,
                show_rownames = TRUE, return_heatmap = TRUE
              )
            })
        }
      } else {
        monocle::plot_pseudotime_heatmap(
          cds, num_clusters = n_cluster,
          show_rownames = TRUE, return_heatmap = TRUE
        )
      }
    }
    p.hp <- callr::r(
      fun_heatmap, list(
        cds = object(x)[genes, ], point = point,
        n_cluster = n_cluster, list_branches
        ),
      libpath = .libPaths(), show = TRUE
    )
    if (!is.null(point)) {
      x$hp_raw <- p.hp
      if (is.null(list_branches)) {
        grob <- p.hp$ph_res$gtable
        p.hp <- wrap_scale_heatmap(
          grob, 12, length(genes), 
          pre_height = 3L, raw = FALSE
        )
      } else {
        layout <- wrap_layout(NULL, length(list_branches), f.w = 1.5)
        p.hp <- patchwork::wrap_plots(
          lapply(p.hp, function(x) x$ph_res$gtable), ncol = layout$ncol
        )
        p.hp <- add(layout, p.hp)
      }
      p.hp <- set_lab_legend(
        p.hp,
        glue::glue("{x@sig} pseudotime genes in heatmap with branch"),
        glue::glue("Monocle 拟时分析热图|||拟时间分支节点上的基因表达变化。图中每一行代表一个基因，每一列代表一个拟时序点，颜色表示基因表达水平（例如从蓝色低表达到红色高表达），通过聚类将具有相似表达变化模式的基因归入同一个 Cluster。该热图用于展示基因沿连续细胞状态或分支区域的表达动态，而不应单独解释为严格发育时间。")
      )
    } else {
      p.hp <- set_lab_legend(
        wrap(p.hp$gtable, 5, 8),
        glue::glue("{x@sig} pseudotime genes in heatmap"),
        glue::glue("Monocle 拟时分析热图|||图中每一行代表一个基因，每一列代表一个拟时序点，颜色表示基因表达水平（例如从蓝色低表达到红色高表达），通过聚类将具有相似表达变化模式的基因归入同一个 Cluster。该热图用于展示基因沿连续细胞状态的表达动态，而不应单独解释为严格发育时间。")
      )
    }
    x <- plotsAdd(x, p.hp)
    return(x)
  })

setMethod("mutate", signature = c(x = "job_monocle2"),
  function(x, ...){
    object(x)@phenoData@data <- dplyr::mutate(object(x)@phenoData@data, ...)
    return(x)
  })

# diff_genes <- differentialGeneTest(cds[expressed_genes,], 
#                                   fullModelFormulaStr = "~Status")
# sig_genes <- subset(diff_genes, qval < 0.01)
# ordering_genes <- rownames(sig_genes)[
#   order(abs(sig_genes$fold_change), decreasing = TRUE)[1:300]
# ]

setMethod("set_remote", signature = c(x = "job_monocle2"),
  function(x, wd)
  {
    x$wd <- wd
    return(x)
  })
