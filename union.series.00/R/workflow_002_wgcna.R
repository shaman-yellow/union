# ==========================================================================
# workflow of WGCNA
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_wgcna <- setClass("job_wgcna", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html"),
    cite = "[@WgcnaAnRPacLangfe2008]",
    method = "R package `WGCNA` used for gene co-expression analysis",
    tag = "wgcna",
    analysis = "WGCNA 分析"
    ))

setGeneric("asjob_wgcna",
  function(x, ...) standardGeneric("asjob_wgcna"))

setMethod("asjob_wgcna", signature = c(x = "job_seurat"),
  function(x, features = NULL, cells = NULL)
  {
    step_message("Use SeuratObject:::subset.Seurat to subset the data.")
    hasIt <- features %in% rownames(object(x))
    message("Features found:")
    print(prop.table(table(hasIt)))
    sub <- e(suppressWarnings(SeuratObject:::subset.Seurat(object(x),
        features = features[ hasIt ], cells = cells
        )))
    log_counts <- as_tibble(sub[[ SeuratObject::DefaultAssay(sub) ]]@scale.data)
    metadata <- as_tibble(sub@meta.data)
    gene_annotation <- tibble::tibble(gene = rownames(sub))
    job_wgcna(metadata, log_counts, gene_annotation, "gene")
  })

job_wgcna <- function(metadata, log_counts,
  gene_annotation)
{
  elist <- new_elist(metadata, log_counts, gene_annotation)
  datExpr0 <- as_wgcData(elist)
  .job_wgcna(object = elist, params = list(datExpr0 = datExpr0))
}

setMethod("step0", signature = c(x = "job_wgcna"),
  function(x){
    step_message("Prepare your data with function `job_wgcna`. ",
      "Note that the ", crayon::red("first column"),
      " in each data (metadata, counts, genes annotation)",
      " were used as ID column of ",
      "corresponding content. \n",
      crayon::red("metadata:"), " with 'sample' and 'group' (optional).\n",
      crayon::red("counts:"), " with id column and then expression columns.\n",
      "genes: ", crayon::red("Traits data"), " could be placed herein; ",
      "a step would serve all numeric data in `genes` as trait data then ",
      "perform statistic."
    )
  })

setMethod("step1", signature = c(x = "job_wgcna"),
  function(x, mutate_name = TRUE){
    step_message("Cluster sample tree.",
    "This do:",
    "generate `x@params$raw_sample_tree`; `x@plots[[ 1 ]]`"
    )
    if (mutate_name) {
      dat <- params(x)$datExpr0
      rownames(dat) <- paste0(object(x)$targets$group, "_", seq_along(object(x)$targets$group))
      raw_sample_tree <- draw_sampletree(dat)
      grob_raw_sample_tree <- as_grob(
        expression(draw_sampletree(dat)), environment()
      )
    } else {
      raw_sample_tree <- draw_sampletree(x$datExpr0)
      grob_raw_sample_tree <- as_grob(
        expression(draw_sampletree(x$datExpr0)), environment()
      )
    }
    p.raw_sample_tree <- wrap(
      grob_raw_sample_tree, 
      min(nrow(x$datExpr0) * .4, 20), 
      min(nrow(x$datExpr0), 10)
    )
    p.raw_sample_tree <- set_lab_legend(
      p.raw_sample_tree,
      glue::glue("{x@sig} sample clustering"),
      glue::glue("样本聚类树|||基于样本整体表达谱计算样本间距离，并采用层次聚类展示样本间的整体相似性，用于评估样本分布情况及排查潜在离群样本。纵轴 Height 表示聚类距离，分支连接高度越低，说明样本表达模式越相近；连接高度越高，则提示样本间整体差异相对较大。")
    )
    x <- plotsAdd(x, p.raw_sample_tree)
    x$raw_sample_tree <- raw_sample_tree
    x <- methodAdd(
      x, "**WGCNA** 是一种基于基因表达数据构建加权共表达网络的系统生物学分析方法，其主要目的是识别具有高度协同表达特征的基因模块，并探索其与临床性状或生物学表型之间的关联。该方法通过计算基因间表达相关性构建共表达网络，并根据拓扑重叠矩阵（TOM）对基因进行模块划分，从而筛选与目标表型显著相关的关键模块及枢纽基因（hub genes）。进一步结合功能富集分析，可揭示相关模块在疾病发生发展或生物学过程中的潜在功能与调控机制，为关键基因筛选及机制研究提供依据。"
    )
    x <- methodAdd(x, "以 R 包 `WGCNA` ⟦pkgInfo('WGCNA')⟧ 对数据作共表达分析{cite_show('WgcnaAnRPacLangfe2008')}。分析方法参考 <{x@info}>。")
    return(x)
  })

setMethod("step2", signature = c(x = "job_wgcna"),
  function(x, height = NULL, size = 10L)
  {
    step_message("Cut sample tree with `height` and `size`. ",
      "This do: ",
      "clip `x@object`; generate `x@params$datExpr`; ",
      "generate `x@params$allTraits`. "
    )
    if (!is.null(size) && !is.null(height)) {
      iskeep <- cut_tree(x$raw_sample_tree, height, size)
      message(glue::glue("Keep: {try_snap(iskeep)}\nDrop:\n{showStrings(which(!iskeep))}"))
      datExpr <- exclude(params(x)$datExpr0, iskeep)
      x$datExpr <- datExpr
      object(x) <- clip_data(object(x), datExpr)
      x <- snapAdd(
        x, "以 `WGCNA::cutreeStatic` (cutHeight = {height}, minSize = {size}) 剪切聚类树，滤掉样本 {showStrings(rownames(x$datExpr0)[!iskeep])}。"
      )
    } else {
      x$datExpr <- x$datExpr0
    }
    x$allTraits <- as_wgcTrait(object(x))
    return(x)
  })

setMethod("step3", signature = c(x = "job_wgcna"),
  function(x, cores = 4, powers = 1:50, ...)
  {
    step_message("Analysis of network topology for soft-thresholding powers. ",
      "This do: ",
      "Generate x@params$sft; plots in `x@plots[[ 3 ]]`. "
    )
    if (is.remote(x)) {
      object <- object(x)
      object(x) <- NULL
      x <- run_job_remote(x, wait = 3L, ...,
        {
          x <- step3(x, powers = seq_len("{max(powers)}"), cores = "{cores}")
        }
      )
      object(x) <- object
    } else {
      e(WGCNA::enableWGCNAThreads(cores))
      sft <- cal_sft(params(x)$datExpr, powers = powers)
      x$sft <- sft
      p.sft <- wrap(plot_sft(sft), 10, 5)
      p.sft <- set_lab_legend(
        p.sft,
        glue::glue("{x@sig} soft thresholding powers"),
        glue::glue("WGCNA 软阈值筛选曲线|||通过 WGCNA 的 pickSoftThreshold 方法评估不同软阈值 power 下网络的无尺度拓扑拟合程度和平均连接度。左图展示软阈值与无尺度拓扑模型拟合指数 signed R² 的关系，红色横线表示参考筛选标准；右图展示不同软阈值下的平均连接度变化。综合选择能够使网络接近无尺度拓扑特征，同时保留适当基因连接度的 power 值用于后续共表达网络构建。")
      )
      x <- plotsAdd(x, p.sft)
      x <- methodAdd(x, "以 `WGCNA::pickSoftThreshold` 预测最佳 soft thresholding powers。")
    }
    return(x)
  })

setMethod("step4", signature = c(x = "job_wgcna"),
  function(x, cores = 4, power = x$sft$powerEstimate, 
    mergeCutHeight = .15, minModuleSize = 200L,
    inherit = TRUE, force = FALSE, ...)
  {
    step_message("One-step network construction and module detection.
      Extra parameters would passed to `cal_module`.
      This do: Generate `x@params$MEs`; plots (net) in `x@plots[[ 4 ]]`.
      By default, red{{x@params$sft$powerEstimate}} is used
      as `power` for WGCNA calculation.
      "
    )
    if (is.remote(x)) {
      message(glue::glue("Note: `...` can not passed to remote."))
      object <- object(x)
      object(x) <- NULL
      x <- run_job_remote(
        x, wait = 3, inherit_last_result = inherit, ...,
        {
          x <- step4(x, power = "{power}", cores = "{cores}")
        }
      )
      object(x) <- object
    } else {
      e(WGCNA::enableWGCNAThreads(cores))
      fun_module <- function(mergeCutHeight, minModuleSize, power)
      {
        net <- cal_module(
          x$datExpr, power,
          minModuleSize = minModuleSize, mergeCutHeight = mergeCutHeight, ...
        )
      }
      net <- expect_local_data(
        "tmp", "module", fun_module, list(
          mergeCutHeight, minModuleSize, power
        ), rerun = force
      )
      if (!is(net, "wgcNet")) {
        net <- .wgcNet(net)
      }
      x$MEs <- get_eigens(net)
      ME_genes <- net$colors
      ME_genes <- split(names(ME_genes), unname(ME_genes))
      names(ME_genes) <- paste0("ME", WGCNA::labels2colors(names(ME_genes)))
      x$ME_genes <- ME_genes
      x <- snapAdd(x, "使用 `WGCNA::blockwiseModules` 函数，设定最小模块基因数为 {minModuleSize} (minModuleSize) 以过滤过小模块，并通过合并切割聚类树高度为 {mergeCutHeight} (mergeCutHeight) 的分支模块，以 power {power} (soft thresholding powers) 创建基因共表达模块 【各模块基因数 (括号中为数目)：{try_snap(ME_genes)}】。")
      p.net <- set_lab_legend(
        wrap(net, 7, 6),
        glue::glue("{x@sig} co-expression module"),
        glue::glue("WGCNA 基因共表达模块聚类图|||基于基因间表达相关性构建加权共表达网络，并根据拓扑重叠矩阵对基因进行层次聚类与模块划分。上方树状图表示基因之间的共表达相似性，分支越接近说明基因表达模式越相似；下方不同颜色代表识别得到的不同共表达模块，同一颜色中的基因具有相似的表达变化趋势，可用于后续模块-性状关联分析及关键模块筛选。")
      )
      x <- plotsAdd(x, p.net)
      x <- methodAdd(x, "选择 power 为 {power}, 以 `WGCNA::blockwiseModules` 创建共表达网络，检测基因模块。")
    }
    return(x)
  })

setMethod("step5", signature = c(x = "job_wgcna"),
  function(x, traits = NULL, group_levels = NULL, cut.p = .05, 
    cut.cor = .3, native = TRUE)
  {
    step_message("Correlation test for modules with trait data. ",
      "This do:",
      "Generate plots in `x@plots[[ 5 ]]`; ",
      "tables in `x@tables[[ 5 ]]`"
    )
    if (!is.null(x$traits)) {
      message("Use `x$traits` for correlation.")
      traits <- x$traits
    }
    if (is.null(traits) && !is.null(group_levels) && !is.null(object(x)$targets[[ "group" ]])) {
      message(glue::glue("Use 'group' in `object(x)$targets`: {bind(group_levels)}."))
      traits <- object(x)$targets
      traits$group <- as.integer(factor(traits$group, levels = group_levels))
      x <- snapAdd(x, "将 'group' 设置为数值变量 ({bind(group_levels)} 依次为 {bind(seq_along(group_levels))}) 与基因共表达模块关联分析。")
    }
    if (!is.null(traits)) {
      .check_columns(traits, c("sample"))
      message("Match rownames in expression data.")
      traits <- traits[match(rownames(x@params$datExpr), traits$sample), ]
      rownames <- traits$sample
      message("The numeric columns will calculate correlation with expression data.")
      traits <- dplyr::select_if(traits, is.numeric)
      traits <- data.frame(traits)
      rownames(traits) <- rownames
      x$allTraits <- .wgcTrait(traits)
    }
    if (is.null(params(x)$allTraits)) {
      stop("is.null(params(x)$allTraits) == TRUE")
    }
    if (ncol(params(x)$allTraits) == 0) {
      stop("ncol(params(x)$allTraits) == 0, no data in `allTraits`.")
    }
    useMEs <- x$MEs[, colnames(x$MEs) != "MEgrey" ]
    if (ncol(params(x)$allTraits) == 1L) {
      traitName <- colnames(x$allTraits)
      cor <- e(WGCNA::cor(useMEs, x$allTraits, use = "p"))
      pvalue <- e(WGCNA::corPvalueStudent(cor, nrow(useMEs)))
      if (!identical(rownames(cor), rownames(pvalue))) {
        stop('!identical(rownames(cor), rownames(pvalue))')
      }
      x$corp_group <- dplyr::bind_cols(cor, pvalue)
      colnames(x$corp_group) <- c("cor", "pvalue")
      x$corp_group <- dplyr::mutate(x$corp_group, MEs = rownames(!!cor), .before = 1)
      x$corp_group <- dplyr::arrange(x$corp_group, dplyr::desc(abs(cor)))
      x$corp_group <- set_lab_legend(
        x$corp_group,
        glue::glue("{x@sig} correlation of module with {traitName}"),
        glue::glue("共表达模块与 {traitName} 的关联性")
      )
      if (native) {
        stop("...")
      } else {
        data <- dplyr::mutate(x$corp_group, group = !!traitName)
        fun_palette <- fun_color(
          values = data$cor, category = "div", rev = TRUE
        )
        p.corhp <- e(
          tidyHeatmap::heatmap(data, MEs, group, cor, palette_value = fun_palette)
        )
        p.corhp <- tidyHeatmap::layer_text(
          p.corhp, .value = signif(pvalue, 4)
        )
      }
      p.corhp <- set_lab_legend(
        wrap(p.corhp, 5),
        glue::glue("{x@sig} correlation heatmap"),
        glue::glue("WGCNA 模块-表型相关性热图|||展示 WGCNA 共表达模块与疾病表型 ({traitName}) 之间的相关性分析结果。每一行代表一个共表达模块，热图颜色表示模块特征基因与表型之间的相关方向和强度，红色表示正相关，蓝色表示负相关，颜色越深说明相关性越强。方格中的数值为相关系数，括号中为对应的 P 值，可用于筛选与疾病表型显著相关的关键模块。")
      )
      x <- plotsAdd(x, p.corhp)
      sigModules <- dplyr::filter(
        x$corp_group, abs(cor) > cut.cor, pvalue < cut.p
        )$MEs
      x <- snapAdd(
        x, "筛选显著关联的共表达模块的基因 (pvalue &lt; {cut.p}, cor &gt; {cut.cor})。"
      )
      x$.feature <- as_feature(
        x$ME_genes[ names(x$ME_genes) %in% sigModules ], x,
        analysis = glue::glue("WGCNA 与 {traitName} 显著关联的共表达模块的基因")
      )
    } else {
      hps_corp <- new_heatdata(useMEs, x$allTraits)
      hps_corp <- callheatmap(hps_corp)
      x <- plotsAdd(x, hps_corp = hps_corp)
      x <- tablesAdd(x, corp = hps_corp@data_long)
    }
    return(x)
  })

setMethod("step6", signature = c(x = "job_wgcna"),
  function(x, use.trait = NULL, use = c("adj.pvalue", "pvalue")){
    step_message("Calculate gene significance (GS) and module membership (MM).",
      "This do:",
      "Generate `x@params$mm`, `x@params$gs`; ",
      "tables (filter by pvalue < 0.05) `x@tables[[ 6 ]]`"
    )
    use <- match.arg(use)
    useMEs <- x$MEs[, colnames(x$MEs) != "MEgrey" ]
    mm <- cal_corp(params(x)$datExpr, useMEs, "gene", "module")
    mm.s <- mutate(as_tibble(mm), adj.pvalue = p.adjust(pvalue, "BH"))
    mm.s <- dplyr::filter(mm.s, !!rlang::sym(use) < .05)
    gs <- cal_corp(params(x)$datExpr, params(x)$allTraits, "gene", "trait")
    gs.s <- dplyr::mutate(tibble::as_tibble(gs), adj.pvalue = p.adjust(pvalue, "BH"))
    gs.s <- dplyr::filter(gs.s, !!rlang::sym(use) < .05)
    if (!is.null(use.trait)) {
      gs.s <- dplyr::filter(gs.s, trait %in% dplyr::all_of(use.trait))
    }
    if (FALSE) {
      p.mm_gs <- new_upset(gs = gs.s$gene, mm = mm.s$gene)
      show(p.mm_gs)
      p.mm_gs <- wrap(recordPlot(), 3, 3)
      dev.off()
      x <- plotsAdd(x, p.mm_gs)
    }
    x$mm <- mm
    x$gs <- gs
    x$ins.mm_gs <- intersect(gs.s$gene, mm.s$gene)
    x <- tablesAdd(x, mm = mm.s, gs = gs.s)
    return(x)
  })

cut_tree <- function(tree, height, size) {
  clust <- e(WGCNA::cutreeStatic(tree, height, size))
  clust > 0
}

cal_sft <- function(data, powers = c(c(1:10), seq(12, 20, by = 2))) 
{
  if (!is(data, "wgcData")) {
    stop("is(data, \"wgcData\") == FALSE")
  }
  sft <- e(WGCNA::pickSoftThreshold(data, powerVector = powers, verbose = 5))
  sft
}

plot_sft <- function(sft) 
{
  p1 <- ggplot(sft$fitIndices, aes(x = Power, y = -sign(slope) * SFT.R.sq)) +
    geom_line(color = "darkred", size = 2, lineend = "round") +
    labs(x = "Soft Threshold (power)",
      y = "Scale Free Topology Model Fit, signed R^2") +
    theme_classic()
  p2 <- ggplot(sft$fitIndices, aes(x = Power, y = mean.k.)) +
    geom_line(color = "darkgreen", size = 2, lineend = "round") +
    labs(x = "Soft Threshold (power)",
      y = "Mean Connectivity") +
    theme_classic()
  require(patchwork)
  p1 + p2
}

cal_module <- function(data, power, save_tom = "tom", minModuleSize = 200L, mergeCutHeight = .15, ...)
{
  if (!is(data, "wgcData")) {
    stop("is(data, \"wgcData\") == FALSE")
  }
  require(WGCNA)
  net <- e(WGCNA::blockwiseModules(
      data, power = power,
      TOMType = "unsigned", reassignThreshold = 0, 
      minModuleSize = minModuleSize, mergeCutHeight = mergeCutHeight,
      numericLabels = TRUE, pamRespectsDendro = FALSE, loadTOM = TRUE,
      saveTOMs = TRUE, saveTOMFileBase = save_tom, verbose = 3, ...
      ))
  .wgcNet(net)
}

setMethod("set_remote", signature = c(x = "job_wgcna"),
  function(x, wd = glue::glue("~/wgcna_{x@sig}")){
    x$wd <- wd
    rem_dir.create(wd, wd = ".")
    return(x)
  })

