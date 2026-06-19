# ==========================================================================
# workflow of mebocost
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_mebocost <- setClass("job_mebocost", 
  contains = c("job"),
  prototype = prototype(
    pg = "mebocost",
    info = c("https://github.com/kaifuchenlab/MEBOCOST"),
    cite = "",
    method = "",
    tag = "mebocost",
    analysis = "MEBOCOST 细胞代谢通讯分析"
    ))

setGeneric("asjob_mebocost",
  function(x, ...) standardGeneric("asjob_mebocost"))

meboFuns <- new.env(parent = emptyenv())

setMethod("asjob_mebocost", signature = c(x = "job_seurat"),
  function(x, dir_cache = create_job_cache_dir(x, "mebocost"), 
    conda_env = pg("mebocostEnv"), group.by = x$group.by)
  {
    message("Convert data.")
    cli::cli_alert_info("Convert to Assay")
    # object <- scCustomize::Convert_Assay(
    #   object(x), "RNA", convert_to = "V3"
    # )
    requireNamespace("Seurat")
    assay <- SeuratObject::DefaultAssay(object(x))
    if (is(object(x)[[ assay ]], "Assay5")) {
      object(x)[[ assay ]] <- as(object(x)[[ assay ]], "Assay")
    }
    hash <- digest::digest(
      list(cells = colnames(object(x)), genes = rownames(object(x))), 
      "xxhash64", serializeVersion = 3
    )
    file_anndata <- file.path(dir_cache, glue::glue("anndata_{hash}.h5ad"))
    if (!file.exists(file_anndata)) {
      activate_env(conda_env, pg("conda"))
      cli::cli_alert_info("sceasy::convertFormat")
      sceasy::convertFormat(
        object(x), from = "seurat", to = "anndata", assay = assay,
        outFile = file_anndata
      )
    } else {
      message(glue::glue('file.exists(file_anndata): {file_anndata}'))
    }
    metadata <- as_tibble(object(x)@meta.data, idcol = "cell")
    levels <- .guess_levels_from_job_seurat(x)
    group.by <- eval(group.by)
    x <- .job_mebocost()
    x$dir_cache <- dir_cache
    x$file_anndata <- file_anndata
    x$metadata <- metadata
    x$levels <- levels
    x$group.by <- group.by
    x$conda_env <- conda_env
    x$hash <- hash
    return(x)
  })

setMethod("step0", signature = c(x = "job_mebocost"),
  function(x){
    step_message("Prepare your data with function `job_mebocost`.")
  })

setMethod("activate", signature = c(x = "job_mebocost"),
  function(x){
    if (is.null(x$conda_env)) {
      stop('is.null(x$conda_env).')
    }
    activate_env(x$conda_env, pg("conda"))
    x$scanpy <- e(reticulate::import("scanpy"))
    x$mebocost <- e(reticulate::import("mebocost")$mebocost)
    return(x)
  })

setMethod("step1", signature = c(x = "job_mebocost"),
  function(x, workers = 10L, species = "human",
    path_mebocost = pg("path_mebocost"), try_compass = FALSE,
    cutoff_exp = 0, cutoff_met = 0, cutoff_prop = 0.01,
    sensor_type = "All")
  {
    step_message("Create mebocost object.")

    x <- activate(x)
    anndata <- x$scanpy$read_h5ad(x$file_anndata)
    x$file_config <- file.path(path_mebocost, "mebocost.conf")

    if (!is.null(x$metadata$group)) {
      x$compare.by <- "group"
    } else {
      x$compare.by <- NULL
    }

    x$mebocost_param <- list(
      species = species,
      cutoff_exp = cutoff_exp,
      cutoff_met = cutoff_met,
      cutoff_prop = cutoff_prop,
      sensor_type = sensor_type,
      workers = workers
    )

    cli::cli_alert_info("Run: mebocost.create_obj")

    object(x) <- x$mebocost$create_obj(
      adata = anndata,
      group_col = x$group.by,
      condition_col = x$compare.by,
      met_est = "mebocost",
      config_path = x$file_config,
      exp_mat = NULL,
      cell_ann = NULL,
      species = species,
      met_pred = NULL,
      met_enzyme = NULL,
      met_sensor = NULL,
      met_ann = NULL,
      scFEA_ann = NULL,
      compass_met_ann = NULL,
      compass_rxn_ann = NULL,
      cutoff_exp = cutoff_exp,
      cutoff_met = cutoff_met,
      cutoff_prop = cutoff_prop,
      sensor_type = sensor_type,
      thread = as.integer(workers)
    )

    if (try_compass) {
      x$file_avgExp <- file.path(x$dir_cache, glue::glue("avgExp_{x$hash}.tsv"))

      if (file.exists(x$file_avgExp)) {
        message(glue::glue("file.exists: {x$file_avgExp}"))
      } else {
        pd <- reticulate::import("pandas")
        np <- reticulate::import("numpy")

        avg_exp <- x$scanpy$get$aggregate(
          anndata, by = list("group", x$group.by), func = "mean"
        )
        cols <- avg_exp$var_names$to_list()
        rows <- paste0(avg_exp$obs[["group"]], " ~ ", avg_exp$obs[[x$group.by]])
        avg_exp <- as.data.frame(avg_exp$layers[["mean"]])
        rownames(avg_exp) <- rows
        colnames(avg_exp) <- cols
        avg_exp <- t(avg_exp)
        avg_exp <- expm1(avg_exp)

        data.table::fwrite(
          avg_exp, x$file_avgExp, sep = "\t", row.names = TRUE
        )
      }
    }

    .describe_cutoff_value <- function(value, label) {
      if (identical(value, "auto")) {
        return(glue::glue(
          "`{label} = auto`，即由 MEBOCOST 根据所有细胞中非零表达/丰度值的分布自动确定阈值。"
        ))
      }

      if (is.numeric(value) && identical(as.numeric(value), 0)) {
        return(glue::glue(
          "`{label} = 0`，即以非零检出作为有效表达/丰度判断标准，保留所有非零信号。"
        ))
      }

      glue::glue(
        "`{label} = {value}`，即使用固定数值 {value} 作为有效表达/丰度判断阈值。"
      )
    }

    .describe_prop_value <- function(value) {
      value_pct <- round(as.numeric(value) * 100, 4L)

      glue::glue(
        "`cutoff_prop = {value}`，即要求代谢物或传感器在对应细胞群体中至少 {value_pct}% 的细胞达到有效检出阈值。"
      )
    }

    .describe_sensor_type <- function(sensor_type) {
      if (length(sensor_type) == 1L && identical(sensor_type, "All")) {
        return("`sensor_type = All`，即纳入 MEBOCOST 数据库中的全部代谢物传感器类型。")
      }

      glue::glue(
        "`sensor_type = {paste(sensor_type, collapse = ', ')}`，即仅纳入指定类型的代谢物传感器。"
      )
    }

    text_cutoff_exp <- .describe_cutoff_value(cutoff_exp, "cutoff_exp")
    text_cutoff_met <- .describe_cutoff_value(cutoff_met, "cutoff_met")
    text_cutoff_prop <- .describe_prop_value(cutoff_prop)
    text_sensor_type <- .describe_sensor_type(sensor_type)

    if (identical(cutoff_exp, 0) && identical(cutoff_met, 0) &&
        isTRUE(abs(cutoff_prop - 0.01) < .Machine$double.eps^0.5)) {
      text_threshold_reason <- paste0(
        "考虑到代谢物传感器、转运体及相关受体基因在单细胞转录组中通常存在较低检出率，",
        "本研究参考既往 MEBOCOST 应用研究的参数设置 (PMID: 40474982)，将 `cutoff_exp` 和 `cutoff_met` 设为 0，",
        "以保留所有非零代谢物或传感器信号，并将 `cutoff_prop` 设为 0.01，",
        "要求至少 1% 细胞达到有效检出阈值。"
      )
    } else if (identical(cutoff_exp, "auto") || identical(cutoff_met, "auto")) {
      text_threshold_reason <- paste0(
        "本研究使用 MEBOCOST 的自动阈值模式或混合阈值模式进行通讯事件筛选；",
        "其中 `auto` 表示根据数据中非零表达/丰度值的分布自动确定有效检出阈值。"
      )
    } else {
      text_threshold_reason <- paste0(
        "本研究根据数据中代谢物和传感器的检出特征设置固定阈值，",
        "以避免在单细胞转录组中因低检出率导致潜在代谢通讯事件被过度过滤。"
      )
    }

    x <- methodAdd(x, "**MEBOCOST** 是一种基于单细胞转录组数据推断代谢物介导细胞间通讯的计算方法，其主要目的是系统性解析不同细胞类型之间通过代谢物–传感器轴所形成的潜在相互作用网络。该方法结合代谢物相关酶表达信息与传感器基因表达谱，推断细胞产生或释放特定代谢物的潜力及其被其他细胞感知的可能性，从而构建代谢通讯关系，并在此基础上识别具有统计显著性的通讯事件及关键代谢物–传感器组合。\n\n")

    x <- methodAdd(x, "将单细胞数据集导入 Python 工具 MEBOCOST 进行细胞间代谢通讯分析。分析物种参数设置为 `{species}`，细胞群体依据 `{x$group.by}` 进行定义；若样本 metadata 中存在分组变量 `group`，则进一步按照该分组分别推断不同条件下的代谢通讯。")

    x <- methodAdd(x, "MEBOCOST 基于内置的代谢物–酶和代谢物–传感器先验知识库，首先聚合代谢物相关酶的表达信息以估计不同细胞群体的代谢物产生潜力，并计算接收细胞中对应传感器基因的表达水平。随后，对于每一组 Sender–Receiver 细胞类型及代谢物–传感器组合，计算代谢物介导的细胞间通讯得分。")

    x <- methodAdd(x, "本研究的 MEBOCOST 初始化参数为：{text_cutoff_exp} {text_cutoff_met} {text_cutoff_prop} {text_sensor_type}")

    x <- methodAdd(x, text_threshold_reason)
    return(x)
  })

setMethod("step2", signature = c(x = "job_mebocost"),
  function(x, species = "homo_sapiens", workers = 1L)
  {
    step_message("Run compass")
    if (x$.args$step1$try_compass && !is.null(x$file_avgExp)) {
      x$dir_compass <- file.path(x$dir_cache, "compass_output")
      x$dir_compass_tmp <- file.path(x$dir_cache, "compass_tmp")
      input <- glue::glue("--data {x$file_avgExp} ")
      output <- glue::glue("--output-dir {x$dir_compass} --temp-dir {x$dir_compass_tmp} ")
      setting <- glue::glue("--num-thread {workers} --species {species} --calc-metabolites --lambda 0")
      system(glue::glue("{pg('compass')} {input} {setting} {output}"))
    }
    return(x)
  })

setMethod("step3", signature = c(x = "job_mebocost"),
  function(x, min_cell_number = 10L, cut.p = .05, use.p = "permutation_test_fdr", rerun = FALSE)
  {
    step_message("Infer communication.")
    x <- activate(x)
    x$use.p <- use.p
    fun_infer <- function(min_cell_number, use.p, cut.p, 
      metadata, ...)
    {
      object(x)$infer_commu(
        n_shuffle = 1000L,
        seed = as.integer(x$seed),
        Return = FALSE,
        thread = as.integer(x$.args$step1$workers),
        save_permuation = TRUE,
        min_cell_number = as.integer(min_cell_number),
        pval_method = use.p,
        pval_cutoff = cut.p
      )
      object(x)
    }
    object(x) <- expect_local_data(
      x$dir_cache, "mebocost", fun_infer, list(
        min_cell_number, use.p, cut.p, x$metadata,
        x$.args$step1$cutoff_exp,
        x$.args$step1$cutoff_met,
        x$.args$step1$cutoff_prop
      ),
      fun_read = x$mebocost$load_obj, fun_save = x$mebocost$save_obj, 
      ext = "pk", rerun = rerun
    )
    x$data_original <- object(x)$original_result
    x$data_cutoff_check <- dplyr::bind_rows(lapply(
        c(0, 0.001, 0.005, 0.01, 0.03, 0.05, 0.10, 0.15),
        function(cutoff_prop) {
          run_check_cutoff_prop(
            data_original = x$data_original,
            cutoff_prop = cutoff_prop,
            pval_method = use.p,
            pval_cutoff = 0.05
          )
        }))
    message(glue::glue("Check cutoff data:\n"))
    print(x$data_cutoff_check)
    t.commu_res <- as_tibble(object(x)$commu_res)
    t.commu_res <- dplyr::filter(t.commu_res, !!rlang::sym(x$use.p) < !!cut.p)
    t.commu_res <- .mutate_get_chain_in_mebocost_table(t.commu_res)
    t.commu_res <- set_lab_legend(
      t.commu_res,
      glue::glue("{x@sig} cell metabolic communication results"),
      glue::glue("酶和传感器共表达检测的显著代谢物介导的细胞间通讯")
    )
    ps.heatmaps <- sapply(x$levels, simplify = FALSE,
      function(group) {
        p <- vis(x, "commu_dotmap", group = group, cut.p = cut.p)
        p <- set_lab_legend(p,
          glue::glue("{x@sig} group {group} significant communication heatmap"),
          glue::glue("Group: {group} 细胞间代谢通讯气泡图|||展示显著的“代谢物–感受器”对在不同细胞对之间的通讯关系，其中纵轴为代谢物及其对应感受器的组合，横轴为具体的发送细胞与接收细胞配对，每个气泡代表一条显著的代谢通讯事件；气泡大小表示通讯强度，数值越大代表该通讯关系越强，气泡颜色表示统计显著性水平。")
        )
      })
    snap_commu <- .stat_table_by_pvalue(
      t.commu_res, n = 5, split = "Condition", use.p = x$use.p, 
      colName = "Chain", target = "细胞代谢通讯", by = "组中检测到"
    )
    x <- snapAdd(x, "通过 MEBOCOST `infer_commu` 在单细胞数据集中一共检测到{aref(ps.heatmaps)} {snap_commu}")
    p.eventnum_bar <- vis(x, "eventnum_bar")
    p.eventnum_bar <- set_lab_legend(
      p.eventnum_bar,
      glue::glue("{x@sig} bar plot of communication events"),
      glue::glue("通讯事件柱状图|||柱状图展示发送方与接收方的通讯数量，横轴为各细胞类型，纵轴为通讯事件数。")
    )
    x$cut.p <- cut.p
    if (use.p == "permutation_test_fdr") {
      x <- methodAdd(x, "对每一对 Sender-Receiver 细胞类型及每一组代谢物-传感器对，计算酶-传感器共表达得分（Sender细胞中代谢物聚合酶表达均值 × Receiver 细胞中传感器表达均值）作为原始通讯强度。通过 1000 次细胞标签置换检验构建零分布，计算经验 p 值，并经 Benjamini‑Hochberg 法进行 FDR 校正，⟦mark$blue('以 FDR &lt; {cut.p} 为阈值筛选出显著的代谢物-传感器结合概率')⟧。")
    } else if (use.p == "permutation_test_pval") {
      x <- methodAdd(x, "对每一对 Sender-Receiver 细胞类型及每一组代谢物-传感器对，计算酶-传感器共表达得分（Sender细胞中代谢物聚合酶表达均值 × Receiver 细胞中传感器表达均值）作为原始通讯强度。通过 1000 次细胞标签置换检验构建零分布，计算经验 p 值，⟦mark$blue('以 p &lt; {cut.p} 为阈值筛选出显著的代谢物-传感器结合概率')⟧。")
    }
    x <- tablesAdd(x, t.commu_res)
    x <- plotsAdd(x, ps.heatmaps, p.eventnum_bar)
    return(x)
  })

setMethod("step4", signature = c(x = "job_mebocost"),
  function(x, flux_pass = TRUE, sig_mccc_only = TRUE, 
    cut.p = x$cut.p, cut.fc = .5, rerun = FALSE)
  {
    step_message("Differential analysis.")
    x <- activate(x)
    compare <- paste0(x$levels[1], "_vs_", x$levels[2])
    cli::cli_alert_info("CommDiff")
    if (flux_pass && is.null(x$is_compass_run)) {
      flux_pass <- FALSE
    }
    fun_diff <- function(...) {
      object(x)$CommDiff(
        comps = as.list(compare),
        sig_mccc_only = sig_mccc_only,
        flux_pass = flux_pass,
        thread = as.integer(x$.args$step1$workers)
      )
      object(x)
    }
    object(x) <- expect_local_data(
      x$dir_cache, "diff", fun_diff, list(flux_pass, sig_mccc_only),
      fun_read = x$mebocost$load_obj, fun_save = x$mebocost$save_obj, 
      ext = "pk", rerun = rerun
    )
    x <- methodAdd(x, "对 {bind(x$levels)} 的细胞代谢通讯进行组间比较。")
    ts.diff_commu <- sapply(compare, simplify = FALSE,
      function(com) {
        data <- tibble::as_tibble(object(x)$diffcomm_res[[com]])
        data <- dplyr::filter(data, abs(Log2FC) > !!cut.fc, !!rlang::sym(x$use.p) < !!cut.p)
        data <- .mutate_get_chain_in_mebocost_table(data)
        data <- set_lab_legend(
          data,
          glue::glue("{x@sig} {com} differential communication"),
          glue::glue("{com} 组间细胞代谢通讯差异分析表格。")
        )
      })
    x <- tablesAdd(x, ts.diff_commu)
    p.diff_flow <- vis(
      x, "diff_flow", compare = compare, cut.p = cut.p, cut.fc = cut.fc
    )
    p.diff_flow <- set_lab_legend(
      p.diff_flow,
      glue::glue("{x@sig} Cellular differential metabolic communication network"),
      glue::glue("细胞差异代谢通讯网络|||{.mebocost_network_note}")
    )
    maxShow <- 5L
    snap_diff <- vapply(names(ts.diff_commu), FUN.VALUE = character(1),
      function(name) {
        data <- ts.diff_commu[[name]]
        data <- dplyr::arrange(data, dplyr::desc(abs(Log2FC)))
        data <- head(data, n = maxShow)
        cp <- .setup_compare_pvalue_with_table(
          data, "Chain", x$use.p, "Log2FC", levels = x$levels
        )
        snap <- .stat_compare_by_pvalue(cp, x$levels, "", mode = "communication")
        glue::glue("对显著性结果按 Log2FC 降序排序，如图{aref(p.diff_flow)}，排名前{nrow(data)}的细胞代谢通讯中，{snap}")
      })
    snap_diff <- bind(snap_diff, co = "\n\n")
    x <- snapAdd(x, "⟦mark$red('对 {bind(x$levels)} 通讯差异分析，一共检测到 {nrow(ts.diff_commu[[1]])} 个显著差异的细胞代谢通讯')⟧。{snap_diff}")
    x <- plotsAdd(x, p.diff_flow)
    return(x)
  })


setMethod("step5", signature = c(x = "job_mebocost"),
  function(x, use.score = c("scale", "raw"),
    axis_level = c("sender_metabolite_receiver", "metabolite_receiver"),
    group_by = NULL,
    axis = c("Sender", "Metabolite_Name", "Sensor", "Receiver"),
    cut.p = NULL,
    significance_cap = 10,
    plot_top_n = 30L,
    plot_sort_by = c("score", "table"))
  {
    step_message("Prioritize candidate metabolic communication axes.")

    use.score <- match.arg(use.score)
    axis_level <- match.arg(axis_level)
    plot_sort_by <- match.arg(plot_sort_by)

    if (is.null(group_by)) {
      group_by <- meboFuns$resolve_axis_group_by(axis_level)
    }

    score_info <- meboFuns$resolve_score_column(x, use.score)
    cut.p <- meboFuns$resolve_p_cutoff(x, cut.p)

    raw <- x@tables$step4$ts.diff_commu[[1L]]

    meboFuns$check_axis_score_input(
      raw,
      axis = axis,
      group_by = group_by,
      score_col = score_info$col,
      p_col = x$use.p
    )

    data <- meboFuns$prepare_axis_event_data(
      raw = raw,
      axis = axis,
      score_col = score_info$col,
      p_col = x$use.p,
      cut.p = cut.p,
      significance_cap = significance_cap
    )

    if (nrow(data) == 0L) {
      warning("No significant positive communication event remained for axis prioritization.")
      return(x)
    }

    data_axis <- meboFuns$summarize_axis_consensus(
      data,
      group_by = group_by,
      score_col = score_info$col
    )

    data_annotation <- meboFuns$summarize_axis_annotation(
      data,
      group_by = group_by
    )
    data_axis <- dplyr::left_join(data_axis, data_annotation, by = group_by)

    data_metabolite <- meboFuns$summarize_metabolite_context(data)
    data_axis <- dplyr::left_join(data_axis, data_metabolite, by = "Metabolite_Name")

    data_axis <- meboFuns$append_overall_score(data_axis)

    data_axis <- dplyr::arrange(data_axis, dplyr::desc(overall_score))

    fshow <- function(value) {
      strx(value, "[^_]+")
    }

    name_axis <- bind(fshow(group_by), co = "——")

    x <- methodAdd(
      x,
      "基于 MEBOCOST 显著差异代谢通讯事件，进一步对 {name_axis} 候选通讯轴进行统一优先级排序。MEBOCOST 是用于推断 metabolite-mediated cell-cell communication (mCCC) 的单细胞分析方法，可整合单细胞转录组、代谢酶、代谢物感受器及代谢通量信息，识别发送细胞产生或释放代谢物并被接收细胞感知的潜在通讯关系 (PMID: 40568942)。\n\n⟦mark$blue('考虑到谷氨酰胺等基础营养代谢物在免疫细胞增殖、细胞因子产生、吞噬及杀菌功能中具有广泛作用，相关代谢物在免疫单细胞通讯推断中可能呈现较高背景活跃度 (PMID: 11533304; PMID: 30360490)。因此，本研究不直接采用通讯事件数量、累积通讯权重或单一通讯强度作为最终排序依据，而采用多证据秩整合策略对候选轴进行优先级评估')⟧。该策略同时整合通讯强度、差异幅度、统计显著性、方向一致性、感受器支持度和细胞通讯对特异性；其中细胞通讯对特异性用于降低广泛覆盖多个 Sender–Receiver 组合的基础代谢物对排序的过度影响。该设计借鉴细胞通讯分析中 magnitude 与 specificity 相区分的思想，以及 TF-IDF/IDF 对广泛出现特征进行降权的原则；多证据秩整合则参考 Rank Product 与 Robust Rank Aggregation 等生物信息学排序整合方法 (PMID: 33024107; PMID: 15327980; PMID: 22247279)。"
    )

    snap_score <- meboFuns$get_overall_score_note(
      model = x$levels[1L],
      group_by = group_by,
      axis_level = axis_level,
      score_label = score_info$label,
      significance_cap = significance_cap
    )

    x <- methodAdd(
      x,
      "候选代谢通讯轴统一优先级评分定义如下：\n\n{snap_score}\n\n"
    )

    if (length(group_by) == 2L) {
      data_plot <- utils::head(data_axis, plot_top_n)
      if (plot_sort_by == "score") {
        data_plot <- dplyr::arrange(data_plot, dplyr::desc(overall_score))
      }

      p.score <- ggplot(
          data_plot,
          aes(x = !!rlang::sym(group_by[1L]), y = !!rlang::sym(group_by[2L]))
        ) +
        geom_point(
          aes(size = overall_score, fill = axis_specificity_rank),
          shape = 21, alpha = 0.85
        ) +
        labs(
          x = fshow(group_by[1L]),
          y = fshow(group_by[2L]),
          size = "Overall score",
          fill = "Cell-pair specificity rank"
        ) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

      p.score <- set_lab_legend(
        wrap_scale_heatmap(
          p.score,
          data_plot[[group_by[1L]]],
          data_plot[[group_by[2L]]]
        ),
        glue::glue("{x@sig} Dotplot for prioritized metabolic communication axes"),
        glue::glue(
          "候选代谢通讯轴综合优先级气泡图|||该图展示基于 MEBOCOST 差异通讯事件获得的候选代谢通讯轴 overall score。",
          "横纵坐标对应于 {bind(fshow(group_by))}；",
          "气泡大小表示 overall score，填充颜色表示细胞通讯对特异性百分位秩。",
          "overall score 基于通讯强度、组间差异幅度、统计显著性、方向一致性、感受器支持度及细胞通讯对特异性进行多证据秩整合，",
          "用于在保留 MEBOCOST 原始差异通讯证据的同时，降低广泛基础代谢背景对候选轴排序的影响。",
          "为保证图形可读性，图中展示 overall score 排名靠前的候选轴。"
        )
      )

      x <- plotsAdd(x, p.score)
    } else if (length(group_by) == 3L &&
        all(c("Sender", "Metabolite_Name", "Receiver") %in% group_by)) {
      data_plot <- utils::head(data_axis, plot_top_n)
      if (plot_sort_by == "score") {
        data_plot <- dplyr::arrange(data_plot, dplyr::desc(overall_score))
      }

      data_plot$axis_label <- apply(
        as.data.frame(data_plot[, group_by, drop = FALSE]),
        1L,
        function(value) {
          paste(as.character(value), collapse = " -> ")
        }
      )
      data_plot$axis_label <- factor(
        data_plot$axis_label,
        levels = rev(unique(data_plot$axis_label))
      )

      p.score <- ggplot(
          data_plot,
          aes(x = overall_score, y = axis_label)
        ) +
        geom_segment(
          aes(x = 0, xend = overall_score, y = axis_label, yend = axis_label),
          size = 0.35,
          color = "grey80"
        ) +
        geom_point(
          aes(size = n_sensor, fill = axis_specificity_rank),
          shape = 21,
          alpha = 0.9
        ) +
        labs(
          x = "Overall score",
          y = "Sender -> Metabolite -> Receiver",
          size = "Sensor support",
          fill = "Cell-pair specificity rank"
        ) +
        theme_bw() +
        theme(
          axis.text.y = element_text(size = 8),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank()
        )

      p.score <- set_lab_legend(
        p.score,
        glue::glue("{x@sig} Ranked plot for prioritized sender-metabolite-receiver axes"),
        glue::glue(
          "候选三元代谢通讯轴综合优先级排序图|||该图展示基于 MEBOCOST 差异通讯事件获得的 Sender–Metabolite–Receiver 候选通讯轴 overall score。",
          "纵轴直接展示发送细胞、代谢物和接收细胞构成的三元通讯轴，横轴表示 overall score；",
          "气泡大小表示支持该通讯轴的感受器数量，填充颜色表示细胞通讯对特异性百分位秩。",
          "overall score 基于通讯强度、组间差异幅度、统计显著性、方向一致性、感受器支持度及细胞通讯对特异性进行多证据秩整合，",
          "因此可在同一坐标尺度下比较不同三元通讯轴的综合优先级，同时避免将基础代谢物的广泛连接度直接等同于关键候选轴优先级。",
          "为保证图形可读性，图中展示 overall score 排名靠前的候选轴。"
        )
      )

      x <- plotsAdd(x, p.score)
    } else {
      p.score <- NULL
    }

    t.overallScore <- set_lab_legend(
      data_axis,
      glue::glue("{x@sig} Prioritized metabolic communication axis table"),
      glue::glue(
        "MEBOCOST 候选代谢通讯轴综合优先级表|||该表展示基于显著差异通讯事件计算得到的候选代谢通讯轴 overall score。",
        "overall_score 为本步骤统一采用的候选轴综合优先级评分，由通讯强度、差异幅度、统计显著性、方向一致性、感受器支持度及细胞通讯对特异性六个百分位秩取平均得到；",
        "axis_specificity 表示代谢物在 Sender–Receiver 细胞通讯对中的分布特异性，数值越高说明该代谢物越集中于较少细胞通讯对；",
        "main_sender 和 main_sensor 分别展示该候选轴中贡献最高的发送细胞和代谢物感受器，用于辅助解释通讯来源和接收端机制。"
      )
    )

    x <- tablesAdd(x, t.overallScore)

    x$data_overall_score <- data_axis

    x$.feature_all_metabolites <- as_feature(
      unique(data_axis$Metabolite_Name), 
      "MEBOCOST 主要代谢物", nature = "compounds"
    )

    snap <- paste(
      as.character(unlist(data_axis[1L, group_by, drop = FALSE], use.names = FALSE)),
      collapse = " -> "
    )

    for (i in group_by) {
      x[[ glue::glue(".feature_{tolower(fshow(i))}") ]] <- as_feature(
        data_axis[[i]][1L],
        "关键通讯轴",
        nature = fshow(i)
      )
    }

    if (!is.null(p.score)) {
      x <- snapAdd(
        x,
        "基于候选代谢通讯轴 overall score，如图{aref(p.score)}，优先级最高的 {name_axis} 通讯轴为 {snap}。"
      )
    } else {
      x <- snapAdd(
        x,
        "基于候选代谢通讯轴 overall score，优先级最高的 {name_axis} 通讯轴为 {snap}。"
      )
    }

    return(x)
  })

meboFuns$resolve_score_column <- function(x, use.score)
{
  if (identical(use.score, "scale")) {
    score_col <- glue::glue("Scaled_Commu_Score_{x$levels[1L]}")
    score_label <- glue::glue("`Scaled_Commu_Score_{x$levels[1L]}`")
  } else {
    score_col <- glue::glue("Commu_Score_{x$levels[1L]}")
    score_label <- glue::glue("`Commu_Score_{x$levels[1L]}`")
  }

  list(
    col = as.character(score_col),
    label = as.character(score_label),
    mode = use.score
  )
}

meboFuns$resolve_p_cutoff <- function(x, cut.p = NULL)
{
  if (!is.null(cut.p)) {
    return(cut.p)
  }

  if (!is.null(x$.args) &&
      !is.null(x$.args$step4) &&
      !is.null(x$.args$step4$cut.p)) {
    return(x$.args$step4$cut.p)
  }

  if (!is.null(x$cut.p)) {
    return(x$cut.p)
  }

  0.05
}

meboFuns$resolve_axis_group_by <- function(axis_level)
{
  if (identical(axis_level, "sender_metabolite_receiver")) {
    return(c("Sender", "Metabolite_Name", "Receiver"))
  }

  if (identical(axis_level, "metabolite_receiver")) {
    return(c("Metabolite_Name", "Receiver"))
  }

  stop("Unsupported `axis_level`: ", axis_level)
}

meboFuns$check_axis_score_input <- function(data, axis, group_by,
  score_col, p_col)
{
  if (length(axis) < 2L) {
    stop("`axis` must contain at least two columns.")
  }

  if (anyDuplicated(axis) > 0L) {
    stop("`axis` cannot contain duplicated column names.")
  }

  if (anyDuplicated(group_by) > 0L) {
    stop("`group_by` cannot contain duplicated column names.")
  }

  if (!all(group_by %in% axis)) {
    stop("All `group_by` columns must be included in `axis`.")
  }

  if (!"Metabolite_Name" %in% group_by) {
    stop("`group_by` must include `Metabolite_Name` for metabolite-axis annotation.")
  }

  needed <- unique(c(axis, group_by, "Log2FC", score_col, p_col))
  missing <- setdiff(needed, colnames(data))

  if (length(missing) > 0L) {
    stop(
      "Missing columns in axis score input: ",
      paste(missing, collapse = ", ")
    )
  }

  invisible(TRUE)
}

meboFuns$prepare_axis_event_data <- function(raw, axis, score_col, p_col,
  cut.p, significance_cap)
{
  data <- dplyr::filter(raw, !!rlang::sym(p_col) < !!cut.p)

  data <- dplyr::select(
    data,
    dplyr::all_of(axis),
    Log2FC,
    !!rlang::sym(score_col),
    !!rlang::sym(p_col)
  )

  data <- data[stats::complete.cases(data[, axis, drop = FALSE]), , drop = FALSE]

  data <- dplyr::mutate(
    data,
    abs_log2fc = abs(Log2FC),
    direction = sign(Log2FC),
    event_weight = !!rlang::sym(score_col) * abs_log2fc,
    neg_log10_p = pmin(
      -log10(pmax(!!rlang::sym(p_col), .Machine$double.xmin)),
      significance_cap
    )
  )

  dplyr::filter(
    data,
    !is.na(!!rlang::sym(score_col)),
    is.finite(!!rlang::sym(score_col)),
    !!rlang::sym(score_col) > 0,
    !is.na(abs_log2fc),
    is.finite(abs_log2fc)
  )
}

meboFuns$summarize_axis_consensus <- function(data, group_by, score_col)
{
  sensor_by <- unique(c(group_by, "Sensor"))

  data_sensor <- dplyr::group_by(data, !!!rlang::syms(sensor_by))
  data_sensor <- dplyr::summarize(
    data_sensor,
    sensor_comm_score = mean(!!rlang::sym(score_col), na.rm = TRUE),
    sensor_abs_log2fc = mean(abs_log2fc, na.rm = TRUE),
    sensor_event_weight = mean(event_weight, na.rm = TRUE),
    sensor_neg_log10_p = mean(neg_log10_p, na.rm = TRUE),
    sensor_n_chain = dplyr::n(),
    sensor_n_sender = dplyr::n_distinct(Sender),
    sensor_direction_consistency = abs(sum(direction * abs_log2fc, na.rm = TRUE)) /
      sum(abs_log2fc, na.rm = TRUE),
    .groups = "drop"
  )

  data_axis <- dplyr::group_by(data_sensor, !!!rlang::syms(group_by))
  data_axis <- dplyr::summarize(
    data_axis,
    axis_strength = mean(sensor_comm_score, na.rm = TRUE),
    axis_abs_log2fc = mean(sensor_abs_log2fc, na.rm = TRUE),
    axis_weight_mean = mean(sensor_event_weight, na.rm = TRUE),
    axis_significance = mean(sensor_neg_log10_p, na.rm = TRUE),
    n_sensor = dplyr::n_distinct(Sensor),
    sensor_support_mean = mean(sensor_n_sender, na.rm = TRUE),
    sensor_support_max = max(sensor_n_sender, na.rm = TRUE),
    sensor_direction_consistency = mean(sensor_direction_consistency, na.rm = TRUE),
    .groups = "drop"
  )

  data_raw <- dplyr::group_by(data, !!!rlang::syms(group_by))
  data_raw <- dplyr::summarize(
    data_raw,
    raw_weight_sum = sum(event_weight, na.rm = TRUE),
    raw_weight_mean = mean(event_weight, na.rm = TRUE),
    n_chain = dplyr::n(),
    n_sender = dplyr::n_distinct(Sender),
    raw_direction_consistency = abs(sum(direction * abs_log2fc, na.rm = TRUE)) /
      sum(abs_log2fc, na.rm = TRUE),
    .groups = "drop"
  )

  data_axis <- dplyr::left_join(data_axis, data_raw, by = group_by)

  dplyr::mutate(
    data_axis,
    direction_consistency = raw_direction_consistency
  )
}

meboFuns$summarize_axis_annotation <- function(data, group_by)
{
  group_sender <- unique(c(group_by, "Sender"))
  data_sender <- dplyr::group_by(data, !!!rlang::syms(group_sender))
  data_sender <- dplyr::summarize(
    data_sender,
    sender_weight = mean(event_weight, na.rm = TRUE),
    .groups = "drop"
  )
  data_sender <- dplyr::arrange(
    data_sender,
    !!!rlang::syms(group_by),
    dplyr::desc(sender_weight)
  )
  data_sender <- dplyr::group_by(data_sender, !!!rlang::syms(group_by))
  data_sender <- dplyr::slice(data_sender, 1L)
  data_sender <- dplyr::ungroup(data_sender)

  if ("Sender" %in% group_by) {
    data_sender <- dplyr::select(
      data_sender,
      dplyr::all_of(group_by),
      main_sender_weight = sender_weight
    )
    data_sender <- dplyr::mutate(data_sender, main_sender = Sender)
  } else {
    data_sender <- dplyr::select(
      data_sender,
      dplyr::all_of(group_by),
      main_sender = Sender,
      main_sender_weight = sender_weight
    )
  }

  group_sensor <- unique(c(group_by, "Sensor"))
  data_sensor <- dplyr::group_by(data, !!!rlang::syms(group_sensor))
  data_sensor <- dplyr::summarize(
    data_sensor,
    sensor_weight = mean(event_weight, na.rm = TRUE),
    .groups = "drop"
  )
  data_sensor <- dplyr::arrange(
    data_sensor,
    !!!rlang::syms(group_by),
    dplyr::desc(sensor_weight)
  )
  data_sensor <- dplyr::group_by(data_sensor, !!!rlang::syms(group_by))
  data_sensor <- dplyr::slice(data_sensor, 1L)
  data_sensor <- dplyr::ungroup(data_sensor)

  if ("Sensor" %in% group_by) {
    data_sensor <- dplyr::select(
      data_sensor,
      dplyr::all_of(group_by),
      main_sensor_weight = sensor_weight
    )
    data_sensor <- dplyr::mutate(data_sensor, main_sensor = Sensor)
  } else {
    data_sensor <- dplyr::select(
      data_sensor,
      dplyr::all_of(group_by),
      main_sensor = Sensor,
      main_sensor_weight = sensor_weight
    )
  }

  dplyr::left_join(data_sender, data_sensor, by = group_by)
}

meboFuns$summarize_metabolite_context <- function(data)
{
  data_pair <- dplyr::distinct(data, Sender, Receiver)
  total_cell_pair <- nrow(data_pair)
  if (!is.finite(total_cell_pair) || total_cell_pair < 1L) {
    total_cell_pair <- 1L
  }

  data <- dplyr::group_by(data, Metabolite_Name)

  data_context <- dplyr::summarize(
    data,
    metabolite_n_receiver = dplyr::n_distinct(Receiver),
    metabolite_n_sender = dplyr::n_distinct(Sender),
    metabolite_n_sensor = dplyr::n_distinct(Sensor),
    metabolite_n_cell_pair = dplyr::n_distinct(paste(Sender, Receiver, sep = "||")),
    metabolite_n_chain = dplyr::n(),
    metabolite_weight_sum = sum(event_weight, na.rm = TRUE),
    metabolite_weight_mean = mean(event_weight, na.rm = TRUE),
    .groups = "drop"
  )

  dplyr::mutate(
    data_context,
    total_cell_pair = total_cell_pair,
    metabolite_cell_pair_fraction = metabolite_n_cell_pair / total_cell_pair,
    axis_specificity = log((1 + total_cell_pair) / (1 + metabolite_n_cell_pair))
  )
}

meboFuns$rank01 <- function(x)
{
  out <- rep(NA_real_, length(x))
  ok <- !is.na(x) & is.finite(x)

  if (sum(ok) == 0L) {
    out[] <- 0
    return(out)
  }

  if (sum(ok) == 1L || length(unique(x[ok])) == 1L) {
    out[ok] <- 1
    out[!ok] <- 0
    return(out)
  }

  out[ok] <- (rank(x[ok], ties.method = "average") - 1) / (sum(ok) - 1)
  out[!ok] <- 0

  out
}

meboFuns$append_overall_score <- function(data)
{
  data <- dplyr::mutate(
    data,
    sensor_support_raw = log1p(n_sensor),
    communication_rank = meboFuns$rank01(axis_strength),
    effect_rank = meboFuns$rank01(axis_abs_log2fc),
    significance_rank = meboFuns$rank01(axis_significance),
    direction_rank = meboFuns$rank01(direction_consistency),
    sensor_support_rank = meboFuns$rank01(sensor_support_raw),
    axis_specificity_rank = meboFuns$rank01(axis_specificity),
    overall_score = rowMeans(
      cbind(
        communication_rank,
        effect_rank,
        significance_rank,
        direction_rank,
        sensor_support_rank,
        axis_specificity_rank
      ),
      na.rm = TRUE
    ),
    evidence_pattern = dplyr::case_when(
      n_sensor >= 2L & axis_specificity_rank >= 0.5 ~ "multi_sensor_specific_axis",
      n_sensor >= 2L ~ "multi_sensor_broad_axis",
      n_sensor == 1L & axis_specificity_rank >= 0.5 ~ "single_sensor_specific_axis",
      TRUE ~ "single_sensor_broad_axis"
    )
  )

  data
}

meboFuns$diagnose_axis_score_bias <- function(x, use.score = c("scale", "raw"),
  axis_level = c("sender_metabolite_receiver", "metabolite_receiver"),
  group_by = NULL,
  axis = c("Sender", "Metabolite_Name", "Sensor", "Receiver"),
  cut.p = NULL,
  significance_cap = 10)
{
  use.score <- match.arg(use.score)
  axis_level <- match.arg(axis_level)

  if (is.null(group_by)) {
    group_by <- meboFuns$resolve_axis_group_by(axis_level)
  }

  score_info <- meboFuns$resolve_score_column(x, use.score)
  cut.p <- meboFuns$resolve_p_cutoff(x, cut.p)

  raw <- x@tables$step4$ts.diff_commu[[1L]]

  meboFuns$check_axis_score_input(
    raw,
    axis = axis,
    group_by = group_by,
    score_col = score_info$col,
    p_col = x$use.p
  )

  data <- meboFuns$prepare_axis_event_data(
    raw = raw,
    axis = axis,
    score_col = score_info$col,
    p_col = x$use.p,
    cut.p = cut.p,
    significance_cap = significance_cap
  )

  data_axis <- meboFuns$summarize_axis_consensus(
    data,
    group_by = group_by,
    score_col = score_info$col
  )

  data_annotation <- meboFuns$summarize_axis_annotation(
    data,
    group_by = group_by
  )
  data_axis <- dplyr::left_join(data_axis, data_annotation, by = group_by)

  data_metabolite <- meboFuns$summarize_metabolite_context(data)
  data_axis <- dplyr::left_join(data_axis, data_metabolite, by = "Metabolite_Name")

  data_axis <- meboFuns$append_overall_score(data_axis)

  dplyr::arrange(data_axis, dplyr::desc(overall_score))
}

meboFuns$get_overall_score_note <- function(model, group_by, axis_level,
  score_label, significance_cap)
{
  axis_text <- if (identical(axis_level, "sender_metabolite_receiver") ||
      all(c("Sender", "Metabolite_Name", "Receiver") %in% group_by)) {
    paste0(
      "本分析将候选通讯轴定义为 Sender–Metabolite–Receiver 三元组，",
      "即保留发送细胞、代谢物和接收细胞三类方向性信息；",
      "Sensor 不作为主轴标签，而作为接收端机制证据进行汇总。"
    )
  } else {
    paste0(
      "本分析将候选通讯轴定义为 Metabolite–Receiver 二元组，",
      "用于概括特定代谢物作用于接收细胞的总体通讯趋势；",
      "不同 Sender 的贡献在同一 Metabolite–Sensor–Receiver 组合内进行平均，",
      "并在结果表中保留 main_sender 作为主要发送细胞来源注释。"
    )
  }

  sensor_text <- if (identical(axis_level, "sender_metabolite_receiver") ||
      all(c("Sender", "Metabolite_Name", "Receiver") %in% group_by)) {
    meboFuns$glue_ex(
      "在三元轴模式下，同一 Sender–Metabolite–Sensor–Receiver 组合内的 Score、|Log2FC| 和 -log10(FDR) 作为该 sensor 对该候选轴的证据；其中 -log10(FDR) 的上限设为 ⟦significance_cap⟧，以避免极小 FDR 对综合排序产生过度影响。"
    )
  } else {
    meboFuns$glue_ex(
      "在二元轴模式下，同一 Metabolite–Sensor–Receiver 组合内先对不同 Sender 的 Score、|Log2FC| 和 -log10(FDR) 求平均；其中 -log10(FDR) 的上限设为 ⟦significance_cap⟧，以避免极小 FDR 对综合排序产生过度影响。"
    )
  }

  formula_text <- meboFuns$glue_ex(
    "\n\n",
    "$$\n",
    "P_i = \\min\\left\\{-\\log_{10}(FDR_i),\\; ⟦significance_cap⟧\\right\\}\n",
    "$$\n\n",
    "$$\n",
    "\\mathrm{AxisStrength}(A) = \\operatorname{mean}_{k \\in K_A}\\left[\\operatorname{mean}_{i \\in A_k}(Score_i)\\right]\n",
    "$$\n\n",
    "$$\n",
    "\\mathrm{AxisEffect}(A) = \\operatorname{mean}_{k \\in K_A}\\left[\\operatorname{mean}_{i \\in A_k}(|Log2FC_i|)\\right]\n",
    "$$\n\n",
    "$$\n",
    "\\mathrm{AxisSignificance}(A) = \\operatorname{mean}_{k \\in K_A}\\left[\\operatorname{mean}_{i \\in A_k}(P_i)\\right]\n",
    "$$\n\n",
    "$$\n",
    "\\mathrm{DirectionConsistency}(A) = \\frac{\\left|\\sum_{i \\in A}\\operatorname{sign}(Log2FC_i)\\times |Log2FC_i|\\right|}{\\sum_{i \\in A}|Log2FC_i|}\n",
    "$$\n\n",
    "$$\n",
    "\\mathrm{SensorSupport}(A) = \\log(1 + n_{sensor,A})\n",
    "$$\n\n",
    "$$\n",
    "\\mathrm{AxisSpecificity}(m) = \\log\\left(\\frac{1 + N_{cell\\ pair}}{1 + n_{cell\\ pair}(m)}\\right)\n",
    "$$\n\n",
    "$$\n",
    "\\mathrm{overall\\_score}(A) = \\frac{1}{6}\\left(\n",
    "R_{strength}(A) + R_{effect}(A) + R_{significance}(A) + R_{direction}(A) + R_{sensor}(A) + R_{specificity}(A)\n",
    "\\right)\n",
    "$$\n\n"
  )

  meboFuns$glue_ex(
    "具体计算中，每条显著通讯事件首先被定义为 Sender–Metabolite–Sensor–Receiver 四元组。",
    "其中，$Score_i$ 表示第 $i$ 条通讯事件的 MEBOCOST 通讯强度指标（⟦score_label⟧）；$Log2FC_i$ 表示 ⟦model⟧ 组相对于对照组的通讯强度对数倍数变化；$FDR_i$ 表示该差异通讯事件的校正显著性。",
    "⟦axis_text⟧",
    "⟦sensor_text⟧",
    "记候选轴为 $A$，其支持的 sensor 集合为 $K_A$，同一候选轴和同一 sensor 下的通讯事件集合为 $A_k$。overall_score 的计算公式如下：",
    "⟦formula_text⟧",
    "其中，$n_{sensor,A}$ 表示候选轴 $A$ 中支持该轴的 distinct Sensor 数量；$N_{cell\\ pair}$ 表示所有显著候选通讯事件中可观察到的 Sender–Receiver 细胞通讯对总数；$n_{cell\\ pair}(m)$ 表示代谢物 $m$ 覆盖的 Sender–Receiver 细胞通讯对数量。",
    "$R_{strength}$、$R_{effect}$、$R_{significance}$、$R_{direction}$、$R_{sensor}$ 和 $R_{specificity}$ 分别表示上述六个证据维度在所有候选轴中的 0–1 百分位秩。",
    "因此，overall_score 越高，表示该候选轴同时具有更高的通讯强度、更大的组间差异、更强的统计显著性、更一致的变化方向、更多的接收端 sensor 支持，以及更强的细胞通讯对特异性。",
    "该 rank-based 策略保留 MEBOCOST 原始差异通讯证据，同时避免单纯由连接数量、Sender 数量、累积通讯强度或广泛基础代谢物背景决定最终优先级。",
    "最终以 overall_score 作为唯一综合优先级评分，并按 overall_score 从高到低筛选候选代谢通讯轴。\n\n"
  )
}

glue_ex <- meboFuns$glue_ex <- function(..., envir = parent.frame(1))
{
  if (!is.environment(envir)) {
    stop("`envir` must be an environment.")
  }

  glue::glue(
    ...,
    .open = "⟦",
    .close = "⟧",
    .envir = envir
  )
}

setMethod("vis", signature = c(x = "job_mebocost"),
  function(x, mode = c("eventnum_bar", "diff_flow", "commu_dotmap"), ...){
    mode <- match.arg(mode)
    if (mode == "eventnum_bar") {
      plot <- .plot_mebocost_eventnum_bar(x, ...)
    } else if (mode == "diff_flow") {
      plot <- .plot_mebocost_diff_flow(x, ...)
    } else if (mode == "commu_dotmap") {
      plot <- .plot_mebocost_commu_dotmap(x, ...)
    }
    file <- tempfile(mode, fileext = ".pdf")
    plot$savefig(file, bbox_inches = "tight")
    as_data_binary(.file_fig(file))
  })

.mutate_get_chain_in_mebocost_table <- function(data) {
  fun_fix <- function(x) s(x, "^[^~]+~ ", "")
  dplyr::mutate(
    data, Chain = paste0(
      fun_fix(Sender), " -> ", Metabolite_Name, " -> ", Sensor, " -> ", fun_fix(Receiver)
    )
  )
}

.plot_mebocost_eventnum_bar <- function(x, ...) {
  celltypes <- unique(x$metadata[[x$group.by]])
  groups <- unique(x$metadata$group)
  orders <- unlist(lapply(celltypes, function(x) paste0(groups, " ~ ", x)))
  object(x)$eventnum_bar(
    sender_focus = c(),
    metabolite_focus = c(),
    sensor_focus = c(),
    receiver_focus = c(),
    ## uncomment and set to focus on one condition
    # conditions  =  ['Primary'],
    xorder = as.list(orders),
    and_or = "and",
    pval_method = x$use.p,
    pval_cutoff = 0.05,
    comm_score_col = "Norm_Commu_Score",
    comm_score_cutoff = 0,
    cutoff_prop = x$.args$step1$cutoff_prop,
    figsize = c(1 + length(celltypes) * .5, 5),
    save = NULL,
    show_plot = FALSE,
    show_num = TRUE,
    include = list("sender-receiver"),
    group_by_cell = TRUE,
    colorcmap = "tab20",
    return_fig = TRUE
  )
}

.plot_mebocost_commu_dotmap <- function(x, group, cut.p = .05, flux_pass = TRUE)
{
  if (flux_pass && is.null(x$is_compass_run)) {
    flux_pass <- FALSE
  }
  object(x)$commu_dotmap(
    sender_focus = c(),
    metabolite_focus = c(),
    sensor_focus = c(),
    receiver_focus = c(),
    conditions = list(group),
    and_or = 'and',
    flux_pass = flux_pass,
    pval_method = x$use.p,
    pval_cutoff = cut.p, 
    cmap = 'Reds',
    cellpair_order = c(),
    met_sensor_order = c(),
    show_plot = FALSE,
    comm_score_col = 'Commu_Score',
    comm_score_range = NULL,
    comm_score_cutoff = NULL,
    cutoff_prop = x$.args$step1$cutoff_prop,
    return_fig = TRUE
  )
}

.plot_mebocost_diff_flow <- function(x, compare, cut.p = .05, cut.fc = .5) {
  object(x)$DiffFlowPlot(
    comp_cond = compare, 
    pval_method = x$use.p,
    pval_cutoff = cut.p,
    Log2FC_threshold = cut.fc,
    sender_focus = c(),
    metabolite_focus = c(),
    sensor_focus = c(),
    receiver_focus = c(),
    remove_unrelevant = TRUE,
    and_or = 'and',
    node_label_size = 8,
    node_alpha = .8,
    figsize = 'auto',
    node_cmap = 'Set1',
    line_color_col = 'Log2FC',
    line_cmap = 'bwr',
    line_cmap_vmin = -2,
    line_cmap_vmax = 2,
    line_cmap_center = 0,
    linewidth = 1.5,
    node_size_norm = c(10, 150),
    node_value_range = NULL,
    save = NULL, 
    save_plot = FALSE, 
    show_plot = FALSE,
    text_outline = FALSE,
    return_fig = TRUE
  )
}

run_check_cutoff_prop <- function(data_original, cutoff_prop, pval_method, pval_cutoff) {
  tibble::tibble(
    cutoff_prop = cutoff_prop,
    pval_method = pval_method,
    pval_cutoff = pval_cutoff,
    n_pass = sum(
      data_original$Commu_Score >= 0 &
        data_original$metabolite_prop_in_sender > cutoff_prop &
        data_original$sensor_prop_in_receiver > cutoff_prop &
        data_original[[pval_method]] < pval_cutoff,
      na.rm = TRUE
    ),
    n_pass_abundance = sum(
      data_original$Commu_Score >= 0 &
        data_original$metabolite_prop_in_sender > cutoff_prop &
        data_original$sensor_prop_in_receiver > cutoff_prop,
      na.rm = TRUE
    ),
    n_pass_sensor = sum(data_original$sensor_prop_in_receiver > cutoff_prop, na.rm = TRUE),
    n_pass_met = sum(data_original$metabolite_prop_in_sender > cutoff_prop, na.rm = TRUE)
  )
}

.mebocost_network_note <- "基于 MEBOCOST 的细胞间差异代谢通讯网络可视化，按照“发送细胞（Sender）–代谢物（Metabolite）–感受器（Sensor）–接收细胞（Receiver）”四层结构展示完整的通讯路径，其中左侧为分泌代谢物的细胞类型，中间依次为参与通讯的小分子及其对应的受体或转运蛋白，右侧为表达感受器并接收信号的细胞类型；连线表示一条代谢通讯关系，颜色根据组间差异分析的 log2FC 显示变化方向与幅度（红色表示上调，蓝色表示下调，颜色越深代表变化越显著），节点大小表示该节点参与的连接数量（即连接度），用于反映其在网络中的参与程度；该网络仅包含经过差异分析筛选后的显著代谢通讯关系，用于整体呈现不同细胞类型之间通过代谢物介导的通讯模式。"

setMethod("set_remote", signature = c(x = "job_mebocost"),
  function(x, wd)
  {
    x$wd <- wd
    return(x)
  })
