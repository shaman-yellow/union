# ========================================================================== 
# workflow of metaInte
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_metaInte <- setClass("job_metaInte", 
  contains = c("job"),
  prototype = prototype(
    pg = "metaInte",
    info = c(""),
    cite = "",
    method = "",
    tag = "metaInte",
    analysis = "整合代谢组与单细胞代谢通讯分析"
    ))

setGeneric("do_metaInte",
  function(x, ref, ...) standardGeneric("do_metaInte"))

setMethod("do_metaInte", signature = c(x = "job_mebocost", ref = "job_metaboDiff"),
  function(x, ref, ...)
  {
    step_message("Prepare MEBOCOST and metabolomics data for integration.")

    job_mebo <- x
    job_metabo <- ref
    x <- .job_metaInte()

    data_commu_res <- tryCatch(
      job_mebo@tables$step3$t.commu_res,
      error = function(e) NULL
    )
    data_diff_commu <- tryCatch(
      job_mebo@tables$step4$ts.diff_commu[[1L]],
      error = function(e) NULL
    )
    data_overall_score <- tryCatch(
      job_mebo@tables$step5$t.overallScore,
      error = function(e) NULL
    )
    data_metabo_diff <- tryCatch(
      job_metabo@tables$step3$data_diff_report,
      error = function(e) NULL
    )

    if (is.null(data_commu_res)) {
      stop("job_mebocost@tables$step3$t.commu_res is NULL. Please run MEBOCOST step3 first.")
    }
    if (is.null(data_diff_commu)) {
      stop("job_mebocost@tables$step4$ts.diff_commu[[1L]] is NULL. Please run MEBOCOST step4 first.")
    }
    if (is.null(data_overall_score)) {
      stop("job_mebocost@tables$step5$t.overallScore is NULL. Please run MEBOCOST step5 first.")
    }
    if (is.null(data_metabo_diff)) {
      stop("job_metaboDiff@tables$step3$data_diff_report is NULL. Please run metaboDiff step3 first.")
    }

    x$data_sources <- list(
      mebocost = list(
        sig = tryCatch(job_mebo@sig, error = function(e) NA_character_),
        levels = tryCatch(as.character(job_mebo$levels), error = function(e) character(0L)),
        data_commu_res = tibble::as_tibble(data_commu_res),
        data_diff_commu = tibble::as_tibble(data_diff_commu),
        data_overall_score = tibble::as_tibble(data_overall_score)
      ),
      metaboDiff = list(
        sig = tryCatch(job_metabo@sig, error = function(e) NA_character_),
        case_group = metaFuns$get_job_value(job_metabo, c("case_group", "case.group", "case")),
        control_group = metaFuns$get_job_value(job_metabo, c("control_group", "control.group", "control")),
        data_diff_report = tibble::as_tibble(data_metabo_diff)
      )
    )

    x$lst_refine <- list()
    x$compare <- list(
      mebocost_levels = x$data_sources$mebocost$levels,
      metabo_case_group = x$data_sources$metaboDiff$case_group,
      metabo_control_group = x$data_sources$metaboDiff$control_group
    )

    return(x)
  })

setMethod("step0", signature = c(x = "job_metaInte"),
  function(x)
  {
    step_message("Prepare your data with function `do_metaInte`.")
    return(x)
  })

setMethod("step1", signature = c(x = "job_metaInte"),
  function(x, p_cutoff = 0.05, vip_cutoff = 1,
    padj_cutoff = NULL, namespace = "name", use_direction = FALSE,
    direction_mode = c("same", "opposite"), fallback_if_empty = TRUE,
    dir_cache = metaFuns$get_metaInte_cache_dir(x, "pubchem"),
    dic_name = NULL, metabo_filter_fun = NULL, check_compare = TRUE,
    manual_cid_mebo = NULL, manual_cid_metabo = NULL, ...)
  {
    step_message("Try exact metabolite-level intersection by PubChem CID.")
    metaFuns$assert_metaInte_data_sources(x)

    direction_mode <- match.arg(direction_mode)

    if (!dir.exists(dir_cache)) {
      dir.create(dir_cache, recursive = TRUE, showWarnings = FALSE)
    }

    metabo_fc_sign <- metaFuns$resolve_metaInte_fc_sign(
      x = x,
      check_compare = check_compare
    )

    lst_refine <- metaFuns$integrate_mebocost_metabolomics(
      data_mebo_score = x$data_sources$mebocost$data_overall_score,
      data_diff_commu = x$data_sources$mebocost$data_diff_commu,
      data_metabo_diff = x$data_sources$metaboDiff$data_diff_report,
      dir_cache = dir_cache,
      p_cutoff = p_cutoff,
      vip_cutoff = vip_cutoff,
      padj_cutoff = padj_cutoff,
      namespace = namespace,
      dic_name = dic_name,
      use_direction = use_direction,
      direction_mode = direction_mode,
      fallback_if_empty = fallback_if_empty,
      metabo_filter_fun = metabo_filter_fun,
      metabo_fc_sign = metabo_fc_sign,
      manual_cid_mebo = manual_cid_mebo,
      manual_cid_metabo = manual_cid_metabo,
      ...
    )

    metaFuns$message_pubchem_mapping_diagnostics(lst_refine)

    lst_refine$param <- list(
      p_cutoff = p_cutoff,
      vip_cutoff = vip_cutoff,
      padj_cutoff = padj_cutoff,
      namespace = namespace,
      use_direction = use_direction,
      direction_mode = direction_mode,
      fallback_if_empty = fallback_if_empty,
      dir_cache = dir_cache,
      metabo_fc_sign = metabo_fc_sign,
      check_compare = check_compare
    )

    x$lst_refine$exact_intersection <- lst_refine

    if (isTRUE(lst_refine$direction_fallback)) {
      warning("Direction consistency screening returned no result; step1 fell back to CID overlap.")
    }

    if (nrow(lst_refine$data_final) == 0L) {
      warning("No exact PubChem CID-overlap candidate was obtained in step1.")
      return(x)
    }

    t.exact_intersection <- dplyr::arrange(
      lst_refine$data_final,
      dplyr::desc(integrated_score)
    )
    t.exact_intersection <- set_lab_legend(
      t.exact_intersection,
      glue::glue("{x@sig} exact metabolite-level integration candidates"),
      glue::glue(
        "MEBOCOST 与代谢组精确代谢物整合表|||通过 PubChem CID 将 MEBOCOST 差异代谢通讯结果",
        "与整体代谢组差异代谢物进行映射整合，用于筛选同时具有代谢通讯意义和整体代谢组差异证据的",
        "候选关键代谢物及其 Receiver 细胞。"
      )
    )
    x <- tablesAdd(x, t.exact_intersection = t.exact_intersection)
    x$lst_refine$exact_intersection$data_final <- t.exact_intersection

 
    meth_text <- glue::glue(
      "为评估单细胞代谢通讯结果与整体代谢组差异结果是否存在精确化合物层面的共同证据，",
      "将 MEBOCOST 差异代谢通讯中的代谢物名称和整体代谢组差异代谢物名称分别映射至 PubChem CID。",
      "代谢组侧默认保留 pvalue &lt; <<p_cutoff>> 且 VIP ≥ <<vip_cutoff>> 的差异代谢物；",
      "若设置调整后 p 值阈值，则进一步要求 padj 符合对应阈值。随后以 PubChem CID 为统一标识进行精确交集分析。",
      "对于每个交集候选，综合评分定义为：\n\n",
      "$$\nIntegratedScore = S_{MEBOCOST} \\times S_{Metabolomics}\n$$\n\n",
      "其中 $S_{MEBOCOST}$ 为 MEBOCOST metabolite–Receiver 网络综合得分的 0–1 标准化值，",
      "$S_{Metabolomics}$ 为整体代谢组差异强度的 0–1 标准化值。",
      "代谢组差异强度由 $|log_2FC|$、VIP 和 $-log_{10}(pvalue)$ 共同构成。",
      if (use_direction) {
        "进一步比较 MEBOCOST 通讯 Log2FC 与整体代谢组 log2FC 的方向。"
      } else {
        ""
      },
      .open = "<<", .close = ">>"
    )
    x <- methodAdd(x, "{meth_text}")

    snap_met <- t.exact_intersection$Metabolite_Name[1L]
    snap_receiver <- t.exact_intersection$Receiver[1L]
    x <- snapAdd(
      x,
      glue::glue(
        "PubChem CID 精确整合共获得 {nrow(t.exact_intersection)} 条候选 metabolite–Receiver 组合；",
        "综合得分最高的候选轴为 {snap_met} → {snap_receiver}。"
      )
    )

    return(x)
  })

setMethod("step2", signature = c(x = "job_metaInte"),
  function(x, mebo_source = c("overall_score", "diff_commu", "commu_res"),
    p_cutoff = 0.05, vip_cutoff = 1L, padj_cutoff = NULL,
    dir_cache = metaFuns$get_metaInte_cache_dir(x, "kegg"),
    manual_kegg_mebo = NULL, manual_kegg_metabo = NULL, ...)
  {
    step_message("Identify KEGG pathways containing evidence from both sources.")
    metaFuns$assert_metaInte_data_sources(x)

    mebo_source <- match.arg(mebo_source)

    data_kegg_bridge <- metaFuns$diagnose_kegg_pathway_overlap_from_data(
      data_commu_res = x$data_sources$mebocost$data_commu_res,
      data_diff_commu = x$data_sources$mebocost$data_diff_commu,
      data_overall_score = x$data_sources$mebocost$data_overall_score,
      data_metabo_diff = x$data_sources$metaboDiff$data_diff_report,
      mebo_source = mebo_source,
      p_cutoff = p_cutoff,
      vip_cutoff = vip_cutoff,
      padj_cutoff = padj_cutoff,
      dir_cache = dir_cache,
      manual_kegg_mebo = manual_kegg_mebo,
      manual_kegg_metabo = manual_kegg_metabo
    )

    x$lst_refine$kegg_bridge <- data_kegg_bridge

    if (nrow(data_kegg_bridge$data_pathway_summary) == 0L) {
      warning("No shared KEGG pathway was identified in step2.")
      return(x)
    }

    t.kegg_bridge <- set_lab_legend(
      data_kegg_bridge$data_pathway_summary,
      glue::glue("{x@sig} KEGG pathway-level bridge candidates"),
      glue::glue(
        "MEBOCOST 与代谢组 KEGG 桥接通路表|||该表展示同时包含 MEBOCOST 代谢通讯代谢物证据",
        "和整体代谢组显著差异代谢物证据的 KEGG 通路，用于进行通路层面的双组学衔接解释。"
      )
    )
    x <- tablesAdd(x, t.kegg_bridge = t.kegg_bridge)
    x$lst_refine$kegg_bridge$data_pathway_summary <- t.kegg_bridge

    p.kegg_bridge <- metaFuns$plot_kegg_bridge_coverage_barplot(
      data_kegg_bridge$data_pathway_summary,
      n_top = 20L,
      min_mebo = 1L,
      min_metabo = 1L,
      rank_by = "balanced_count"
    )

    p.kegg_bridge <- set_lab_legend(
      p.kegg_bridge,
      glue::glue("{x@sig} KEGG pathway-level bridge coverage"),
      glue::glue(
        "KEGG 通路桥接覆盖图|||该图展示 MEBOCOST 代谢通讯代谢物与整体代谢组显著差异代谢物在 KEGG 通路层面的共同覆盖情况。",
        "每一行代表一条 shared KEGG pathway，两个分面分别表示 MEBOCOST 和 Metabolomics 两个证据来源。",
        "横轴及柱上数字表示该数据源映射到对应通路的代谢物数量。",
        "该图用于展示两类组学结果在通路层面的共同定位关系，而非通路富集显著性。"
      )
    )

    x <- plotsAdd(x, p.kegg_bridge)

    x <- methodAdd(
      x,
      glue::glue(
        "由于单细胞转录组推断得到的代谢通讯分子与整体代谢组实际检测到的化合物可能并不完全相同，",
        "此处采用 KEGG Compound 和 KEGG Pathway 进行通路层面整合。",
        "⟦mark$blue('MEBOCOST 侧使用 `{mebo_source}` 统计的代谢物集合代谢组侧保留 pvalue < {p_cutoff} 且 VIP >= {vip_cutoff} 的差异代谢物')⟧。",
        "两类代谢物分别映射至 KEGG Compound 后，再通过 KEGG pathway membership 识别共同通路。",
        "⟦mark$blue('某一通路被定义为 KEGG 桥接通路，需要同时满足：至少包含一个 MEBOCOST 代谢通讯代谢物，",
        "且至少包含一个整体代谢组显著差异代谢物')⟧。因此，该结果代表通路层面的共同代谢背景。"
      )
    )

    snap_path <- t.kegg_bridge$pathway_name[1L]
    snap_n <- nrow(t.kegg_bridge)
    x <- snapAdd(
      x,
      glue::glue(
        "KEGG 通路层面共识别到 {snap_n} 条同时包含两类组学证据的候选桥接通路{aref(p.kegg_bridge)}；",
        "代谢组显著性最靠前的通路为 {snap_path}。"
      )
    )

    return(x)
  })

setMethod("step3", signature = c(x = "job_metaInte"),
  function(x, method = c("hypergeom", "diffusion", "pagerank"),
    compound_source = c("shared_pathway", "both", "mebocost", "metabolomics"),
    organism = "hsa", threshold = 0.05, approx = "normality",
    dir_fella = .prefix(paste0("fella_", organism), "db"),
    rebuild = FALSE, repair = FALSE, compounds_background = NULL,
    niter = 100L, ...)
  {
    step_message("Run FELLA with compounds from the KEGG bridge.")

    method <- match.arg(method)
    compound_source <- match.arg(compound_source)

    data_fella <- metaFuns$run_fella_from_kegg_bridge(
      data_kegg_bridge = x$lst_refine$kegg_bridge,
      dir_fella = dir_fella,
      organism = organism,
      compound_source = compound_source,
      method = method,
      approx = approx,
      threshold = threshold,
      rebuild = rebuild,
      repair = repair,
      compounds_background = compounds_background,
      niter = niter
    )

    x$lst_refine$fella <- data_fella

    if (nrow(data_fella$data_table) == 0L) {
      warning("FELLA returned no enriched KEGG node under the current threshold in step3.")
      return(x)
    }

    data_fella_bridge <- metaFuns$get_fella_supported_kegg_bridge(
      data_kegg_bridge = x$lst_refine$kegg_bridge,
      data_fella_result = data_fella,
      fella_p_cutoff = threshold,
      require_fella_sig = TRUE
    )

    x$lst_refine$fella_bridge_core <- data_fella_bridge

    t.fella <- set_lab_legend(
      data_fella$data_table,
      glue::glue("{x@sig} FELLA KEGG network enrichment results"),
      glue::glue(
        "FELLA KEGG 网络富集结果|||该表展示基于 KEGG 知识图谱得到的原始通路/网络富集节点，",
        "用于后续筛选获得 FELLA-supported KEGG bridge pathways。"
      )
    )

    t.fella_bridge_core <- set_lab_legend(
      data_fella_bridge,
      glue::glue("{x@sig} FELLA-supported KEGG bridge pathways"),
      glue::glue(
        "FELLA 支持的 KEGG 桥接通路|||该表展示同时满足 KEGG 双组学共同定位和 FELLA 富集支持的 KEGG 通路。",
        "这些通路同时包含 MEBOCOST 代谢通讯代谢物和整体代谢组显著差异代谢物，并在 FELLA 知识图谱富集中达到设定阈值。"
      )
    )

    x <- tablesAdd(x, t.fella = t.fella, t.fella_bridge_core = t.fella_bridge_core)

    x$lst_refine$fella$data_table <- t.fella
    x$lst_refine$fella_bridge_core <- t.fella_bridge_core

    if (method == "hypergeom") {
      text_method <- glue::glue(
        "FELLA 超几何富集 (hypergeom) 用于评估输入 KEGG compound 集合是否在特定 KEGG 通路中呈现超过随机期望的聚集。",
        "设背景 KEGG 图谱中共有 $N$ 个可检测化合物，某通路包含 $K$ 个化合物，输入集合包含 $n$ 个化合物，",
        "其中 $k$ 个落入该通路，则富集 p 值为：\n\n",
        "$$\np = \\sum_{i=k}^{\\min(K,n)} \\frac{\\binom{K}{i}\\binom{N-K}{n-i}}{\\binom{N}{n}}\n$$\n\n",
        .open = "<<", .close = ">>"
      )
    } else {
      text_method <- glue::glue(
        "FELLA `{method}` 方法基于 KEGG compound、reaction、enzyme、module 和 pathway 构成的知识图谱，",
        "评估输入代谢物在 KEGG 网络中的传播关联和子网络富集程度。"
      )
    }

    x <- methodAdd(
      x, "进一步采用 FELLA 对 KEGG 桥接相关化合物进行知识图谱层面的通路/网络富集分析。输入 compound 来源为 '{compound_source}'，⟦mark$blue('富集方法为 {method}，阈值为 {threshold} (&lt; {threshold})')⟧。{text_method} FELLA 富集结果随后与 KEGG 双组学共同通路取交集 {aref(x@plots$step2$p.kegg_bridge)}，仅保留同时满足 MEBOCOST 代谢通讯证据、整体代谢组差异代谢物证据和 FELLA 富集支持的通路。因此，本步骤得到的结果定义为 FELLA-supported KEGG bridge pathways，而不是单纯的合并代谢物富集结果。"
    )

    p.fella_sankey <- metaFuns$plot_fella_bridge_sankey(
      data_fella_bridge = data_fella_bridge,
      data_kegg_bridge = x$lst_refine$kegg_bridge,
      n_pathway = 8L,
      n_metabolite = 40L,
      p_cutoff = threshold
    )

    p.fella_sankey <- set_lab_legend(
      p.fella_sankey,
      glue::glue("{x@sig} FELLA-supported pathway-metabolite bridge"),
      glue::glue(
        "FELLA 支持的通路-代谢物桥接图|||该图展示 FELLA 富集通路与其支持代谢物之间的连接关系。",
        "左侧气泡表示 FELLA 富集通路的富集特征，气泡大小表示通路命中化合物数量，颜色表示 FELLA p value 的 -log10 转换值；",
        "中间为 FELLA-supported KEGG pathway，右侧为映射到这些通路中的支持代谢物。",
        "连线颜色表示代谢物证据来源，其中 MEBOCOST 表示单细胞代谢通讯侧代谢物，Metabolomics 表示整体代谢组显著差异代谢物。",
        "该图用于展示两类组学代谢物如何在 FELLA 支持的 KEGG 通路中形成桥接关系。"
      )
    )

    x <- plotsAdd(x, p.fella_sankey)

    snap_fella <- t.fella_bridge_core$pathway_name[1L]
    snap_p <- signif(t.fella_bridge_core$fella_pvalue[1L], 3L)

    x <- snapAdd(
      x,
      glue::glue(
        "通过 FELLA 富集结果与 KEGG 双组学共同通路取交集，共获得 {nrow(t.fella_bridge_core)} 条 FELLA-supported KEGG bridge pathways {aref(p.fella_sankey)}；",
        "其中 FELLA 支持最强的桥接通路为 {snap_fella}，FELLA pvalue = {snap_p}。"
      )
    )

    return(x)
  })

setMethod("step4", signature = c(x = "job_metaInte"),
  function(x, fella_p_cutoff = 0.05, require_fella_sig = TRUE,
    weights = c(fella = 0.35, metabo = 0.25, coverage = 0.20,
      specificity = 0.10, balance = 0.10),
    exclude_pathway_id = character(0L), exclude_pathway_name = character(0L),
    bridge_score_col = "bridge_integrated_score", n_top_per_pathway = 5L)
  {
    step_message("Score FELLA-supported KEGG bridge pathways and candidate metabolite-receiver axes.")

    data_fella_bridge <- tibble::as_tibble(x$lst_refine$fella_bridge_core)

    if (nrow(data_fella_bridge) == 0L) {
      warning("No FELLA-supported KEGG bridge pathway was retained in step4.")
      return(x)
    }

    data_fella_bridge_ranked <- metaFuns$score_fella_supported_kegg_bridge(
      data_fella_bridge,
      weights = weights,
      exclude_pathway_id = exclude_pathway_id,
      exclude_pathway_name = exclude_pathway_name
    )

    t.fella_bridge_ranked <- set_lab_legend(
      data_fella_bridge_ranked,
      glue::glue("{x@sig} ranked FELLA-supported KEGG bridge pathways"),
      glue::glue("FELLA 支持的 KEGG 桥接通路综合排序表|||根据 FELLA 富集显著性、代谢组显著性、双组学覆盖度、通路特异性及双组学证据平衡性对候选桥接通路进行排序。")
    )

    x <- tablesAdd(x, t.fella_bridge_ranked)

    if (nrow(data_fella_bridge_ranked) == 0L) {
      warning("All FELLA-supported bridge pathways were excluded or removed in step4.")
      return(x)
    }

    data_bridge_candidate <- metaFuns$get_bridge_metabolite_receiver_candidates(
      data_bridge = data_fella_bridge_ranked,
      data_overall_score = x$data_sources$mebocost$data_overall_score,
      bridge_score_col = bridge_score_col,
      n_top_per_pathway = n_top_per_pathway
    )

    t.bridge_candidate <- set_lab_legend(
      data_bridge_candidate,
      glue::glue("{x@sig} candidate metabolite-receiver axes from bridge pathways data"),
      glue::glue("候选代谢物-Receiver 细胞通讯轴|||在 FELLA-supported KEGG bridge pathway 背景下，回到 MEBOCOST 网络综合得分筛选得到的候选 metabolite–receiver 通讯轴。")
    )

    x <- tablesAdd(x, t.bridge_candidate)

    if (nrow(data_bridge_candidate) == 0L) {
      warning("No metabolite-receiver candidate axis was obtained from FELLA-supported pathways in step4.")
      return(x)
    }

    x$lst_refine$fella_bridge_ranked <- t.fella_bridge_ranked
    x$lst_refine$bridge_candidates <- t.bridge_candidate

    text_weights <- glue::glue(
      "FELLA = {weights[['fella']]}, 代谢组 = {weights[['metabo']]}, 覆盖度 = {weights[['coverage']]}, ",
      "通路特异性 = {weights[['specificity']]}, 平衡性 = {weights[['balance']]}"
    )

    meth_text <- glue::glue(
      "在 KEGG pathway-level bridge 和 FELLA 富集结果的基础上，进一步筛选候选关键 metabolite–Receiver 通讯轴。",
      "⟦mark$blue('首先保留同时满足 KEGG 双组学共同定位且 FELLA pvalue &lt; <<fella_p_cutoff>> 的桥接通路')⟧；",
      "随后根据通路综合桥接得分进行排序。通路综合桥接得分由五个标准化分量构成：\n\n",
      "$$\nBridgeScore = w_1F + w_2M + w_3C + w_4S + w_5B\n$$\n\n",
      "其中，$F=-log_{10}(p_{FELLA})$ 表示 FELLA 富集显著性，",
      "$M=-log_{10}(p_{Metabolomics})$ 表示该通路中整体代谢组差异代谢物的最佳显著性，",
      "$C=\\sqrt{log(1+n_{MEBOCOST})\\times log(1+n_{Metabolomics})}$ 表示双组学覆盖度，",
      "$S=CompoundHits/CompoundsInPathway$ 表示通路特异性，",
      "$B=\\min(log(1+n_{MEBOCOST}),log(1+n_{Metabolomics}))/\\max(log(1+n_{MEBOCOST}),log(1+n_{Metabolomics}))$ 表示双组学证据平衡性。",
      "各分量经 0–1 标准化后加权求和，当前权重为：<<text_weights>>。\n\n",
      "最后，在每条 FELLA-supported KEGG bridge pathway 内，将通路桥接得分与 MEBOCOST metabolite–Receiver overall_score 结合，",
      "计算候选通讯轴得分：\n\n",
      "$$\nCandidateScore = 0.65 \\times S_{overall} + 0.35 \\times S_{bridge}\n$$\n\n",
      "其中 $S_{overall}$ 为 MEBOCOST 网络综合得分的 0–1 标准化值，",
      "$S_{bridge}$ 为对应通路桥接得分的 0–1 标准化值。",
      "该评分用于在通路层面证据支持下筛选关键代谢物及关键 Receiver 细胞。",
      .open = "<<", .close = ">>"
    )
    x <- methodAdd(x, "{meth_text}")

    p.bridge_candidate <- metaFuns$plot_bridge_candidate_axis_grouped(
      data_bridge_candidate,
      n_top_pathway = 6L,
      n_top_per_pathway = n_top_per_pathway,
      pathway_rank_by = "bridge_score"
    )

    p.bridge_candidate <- set_lab_legend(
      p.bridge_candidate,
      glue::glue("{x@sig} candidate metabolite-receiver axes from bridge pathways"),
      glue::glue(
        "候选代谢物-Receiver 通讯轴综合图|||该图展示 FELLA-supported KEGG bridge pathway 背景下筛选得到的候选 MEBOCOST 代谢通讯轴。",
        "左侧为 KEGG 桥接通路及其综合桥接得分 BridgeScore，条形长度和数值均表示通路层面的桥接优先级；",
        "右侧为对应通路内的 metabolite–Receiver 通讯轴，横轴为 Receiver 细胞类型，纵向标签为 MEBOCOST 代谢物。",
        "气泡大小和颜色均表示候选通讯轴得分 CandidateScore。",
        "该图同时展示通路层面的桥接优先级和通路内候选通讯轴的优先级，用于筛选关键代谢物及关键 Receiver 细胞。"
      )
    )

    x <- plotsAdd(x, p.bridge_candidate)

    data_axis_ranked <- tibble::as_tibble(t.bridge_candidate)

    data_axis_ranked <- dplyr::arrange(
      data_axis_ranked,
      dplyr::desc(candidate_score),
      dplyr::desc(bridge_score),
      dplyr::desc(overall_score),
      fella_pvalue,
      best_metabo_pvalue
    )

    snap_met <- data_axis_ranked$Metabolite_Name[1L]
    snap_receiver <- data_axis_ranked$Receiver[1L]
    snap_path <- data_axis_ranked$pathway_name[1L]
    snap_axis_score <- signif(data_axis_ranked$candidate_score[1L], 3L)
    snap_bridge_score <- signif(data_axis_ranked$bridge_score[1L], 3L)

    x <- snapAdd(
      x,
      glue::glue(
        "综合 KEGG 双组学通路桥接、FELLA 网络富集支持和 MEBOCOST 通讯轴得分后{aref(p.bridge_candidate)}，",
        "⟦mark$red('CandidateScore 最高的候选 metabolite–Receiver 通讯轴为 {snap_met} → {snap_receiver}')⟧，",
        "候选得分为 {snap_axis_score}。该轴所属桥接通路为 {snap_path}，",
        "其通路层面 BridgeScore 为 {snap_bridge_score}。"
      )
    )
    return(x)
  })

setMethod("step5", signature = c(x = "job_metaInte"),
  function(x, pathway_id = NULL, select_by = c("candidate_score", "bridge_score"),
    species = "hsa", out_dir = create_job_cache_dir(x),
    clean_old = FALSE, ...)
  {
    step_message("Visualize selected KEGG pathway by pathview.")

    select_by <- match.arg(select_by)

    lst_pathview <- metaFuns$run_pathview_for_bridge_pathway(
      x = x,
      pathway_id = pathway_id,
      select_by = select_by,
      species = species,
      out_dir = out_dir,
      clean_old = clean_old,
      ...
    )

    t.pathview_evidence <- set_lab_legend(
      tibble::as_tibble(lst_pathview$data_evidence),
      glue::glue("{x@sig} pathview compound evidence table"),
      glue::glue(
        "Pathview 化合物映射证据表|||该表展示选定 KEGG 通路中用于 pathview 可视化的 compound 映射结果。",
        "source 表示证据来源，log2FC 表示对应组学在该 compound 节点上的变化方向和幅度。",
        "MEBOCOST_Log2FC 来源于单细胞代谢通讯差异结果，Metabolomics_Log2FC 来源于整体代谢组差异结果。"
      )
    )

    x <- tablesAdd(x, t.pathview_evidence = t.pathview_evidence)
    x$lst_refine$pathview <- lst_pathview
    x$lst_refine$pathview$data_evidence <- t.pathview_evidence

    if (!is.null(lst_pathview$plot)) {
      p.pathview <- set_lab_legend(
        lst_pathview$plot,
        glue::glue("{x@sig} pathview visualization of {lst_pathview$hsa_id}"),
        glue::glue(
          "KEGG 通路 pathview 可视化图|||该图展示选定 KEGG 通路中 MEBOCOST 和整体代谢组数据在 compound 节点上的映射情况。",
          "节点颜色表示 log2FC 方向和幅度；当同一 compound 同时存在 MEBOCOST_Log2FC 与 Metabolomics_Log2FC 时，pathview 以 multi-state 形式在同一节点中显示两类数据。",
          "该图用于展示候选桥接通路中两类组学信号在 KEGG 通路图上的定位。"
        )
      )

      x <- plotsAdd(x, p.pathview)
    } else {
      warning("Pathview did not return a PNG file that could be imported as a plot.")
    }

    x <- methodAdd(
      x,
      glue::glue(
        "在候选 metabolite–Receiver 通讯轴筛选结果的基础上，进一步采用 pathview 对选定 KEGG 通路进行可视化。",
        "该分析将 MEBOCOST 差异通讯结果和整体代谢组差异结果分别整理为 MEBOCOST_Log2FC 与 Metabolomics_Log2FC 两个 compound-level 状态，",
        "并映射至 KEGG compound 节点，以展示两类组学信号在同一通路背景下的定位和变化方向。"
      )
    )

    x <- snapAdd(x,
      glue::glue(
        "对选定 KEGG 通路 {lst_pathview$hsa_id} 进行了 pathview 可视化{aref(p.pathview)}。",
        "该图用于展示 MEBOCOST 与整体代谢组 log2FC 信号在 KEGG compound 节点上的共同定位。"
      )
    )

    return(x)
  })


# ==========================================================================

metaFuns <- new.env(parent = emptyenv())

# ------------------------------------------------------------------------------
# metaFuns helper modules for MEBOCOST-metabolomics integration
# ------------------------------------------------------------------------------
#
# Purpose
# -------
# This environment stores helper functions used to connect single-cell MEBOCOST
# metabolic communication results with bulk metabolomics differential results.
# These functions are intentionally isolated in `metaFuns` to avoid polluting the
# package/global namespace.
#
# Overall integration logic
# -------------------------
# The integration is not based only on exact metabolite-name matching, because
# MEBOCOST metabolites and bulk metabolomics features often use different naming
# systems and may not represent identical detected compounds. The workflow is:
#
#   1. Try exact metabolite-level mapping using PubChem CID or KEGG Compound ID.
#   2. If exact overlap is absent, map both data sources to KEGG pathways.
#   3. Use KEGG/FELLA to identify pathway-level bridges between:
#        - MEBOCOST metabolic communication metabolites, and
#        - bulk metabolomics significant differential metabolites.
#   4. Interpret shared or FELLA-supported pathways as candidate bridge
#      mechanisms rather than direct compound-level validation.
#
# Main input classes
# ------------------
# `kegg_results_metaboDiff_and_mebocost`
#   A structured object returned by KEGG bridge diagnosis. It contains:
#     - data_pathway_summary: shared KEGG pathway summary
#     - data_mebo_map: MEBOCOST metabolite -> KEGG compound mapping
#     - data_metabo_map: metabolomics feature -> KEGG compound mapping
#     - data_mebo_path: MEBOCOST compound -> KEGG pathway mapping
#     - data_metabo_path: metabolomics compound -> KEGG pathway mapping
#     - data_metabo_sig: significant metabolomics feature table
#     - diag_mebo / diag_metabo: mapping diagnostics
#
# `fella_results_metaboDiff_and_mebocost`
#   A structured object returned by combined FELLA analysis. This is useful for
#   checking whether KEGG-shared pathways are also supported by FELLA enrichment.
#
# `fella_dual_bridge_metaboDiff_and_mebocost`
#   A structured object returned by dual-source FELLA bridge analysis. This runs
#   FELLA separately for MEBOCOST and metabolomics compounds, then compares the
#   resulting KEGG nodes/pathways. This is stricter and may return few or no
#   shared results when the metabolomics input is small.
#
# ------------------------------------------------------------------------------
# Module 1: Name normalization and ID mapping
# ------------------------------------------------------------------------------
#
# normalize_metabolite_name()
#   Normalizes metabolite names and applies optional manual synonym replacement.
#
# map_pubchem_cids()
#   Maps metabolite names to PubChem CID using PubChemR. Supports manual CID
#   supplements and returns mapping diagnostics.
#
# map_names_to_kegg_compound()
#   Maps metabolite names to KEGG Compound IDs using the KEGG compound name table.
#   Supports manual KEGG ID supplements.
#
# diagnose_kegg_mapping()
#   Reports KEGG mapping quality, including unmapped and multi-mapped compounds.
#
# ------------------------------------------------------------------------------
# Module 2: MEBOCOST and metabolomics coverage diagnosis
# ------------------------------------------------------------------------------
#
# get_mebocost_metabolite_universe()
#   Extracts MEBOCOST metabolite names from step3 communication results, step4
#   differential communication results, or step5 overall score results.
#
# filter_metabo_significant()
#   Filters metabolomics features by pvalue, VIP, and/or adjusted pvalue.
#
# diagnose_mebocost_metabo_coverage()
#   Tests whether significant metabolomics features overlap with MEBOCOST
#   metabolites at the exact CID/name level. This is a diagnostic step and does
#   not produce a report-level conclusion by itself.
#
# ------------------------------------------------------------------------------
# Module 3: KEGG pathway bridge diagnosis
# ------------------------------------------------------------------------------
#
# get_kegg_pathways_for_compounds()
#   Retrieves KEGG pathway memberships for KEGG Compound IDs.
#
# diagnose_kegg_pathway_overlap()
#   Main entry function for KEGG pathway-level bridge diagnosis. It maps MEBOCOST
#   metabolites and significant metabolomics features to KEGG compounds/pathways,
#   then identifies pathways containing evidence from both sources.
#
# Expected output class:
#   "kegg_results_metaboDiff_and_mebocost"
#
# Interpretation:
#   A shared KEGG pathway is not a direct metabolite-level validation. It means
#   that MEBOCOST communication metabolites and metabolomics differential
#   metabolites converge on the same KEGG pathway.
#
# ------------------------------------------------------------------------------
# Module 4: FELLA-based pathway/network support
# ------------------------------------------------------------------------------
#
# assert_kegg_results_metaboDiff_and_mebocost()
#   Validates that the KEGG bridge object has the expected class and fields.
#
# get_or_build_fella_data()
#   Loads an existing FELLA KEGG graph/data object or builds it if needed. The
#   FELLA database should usually be stored in a shared upper-level directory
#   because it depends on organism and KEGG/FELLA settings, not on a specific
#   project.
#
# run_fella_from_kegg_bridge()
#   Combined FELLA enrichment. It uses MEBOCOST and metabolomics compounds
#   together as one input set. This is useful for broad network support, but it
#   should not be interpreted as independent two-omics bridge evidence.
#
# run_fella_dual_bridge_from_kegg_bridge()
#   Dual-source FELLA bridge analysis. It runs FELLA separately for MEBOCOST and
#   metabolomics compounds, then compares their FELLA results. This is stricter
#   and may be underpowered if the metabolomics significant compound set is small.
#
# get_fella_supported_kegg_bridge()
#   Recommended practical bridge function. It first uses KEGG overlap to define
#   pathways supported by both data sources, then checks whether those pathways
#   are supported in combined FELLA enrichment.
#
# ------------------------------------------------------------------------------
# Module 5: Optional visualization support
# ------------------------------------------------------------------------------
#
# run_pathview_safe()
#   A safety wrapper for pathview. It forces all pathview output files into a
#   specified directory, because pathview may otherwise write KEGG files and
#   rendered pathway images to the current working directory.
#
# Recommended usage:
#   Use pathview only for one or two final candidate pathways, not as the main
#   integration method.
#
# ------------------------------------------------------------------------------
# Recommended entry points
# ------------------------------------------------------------------------------
#
# 1. First check exact metabolite overlap:
#      diagnose_mebocost_metabo_coverage()
#
# 2. If exact overlap is absent, run KEGG pathway bridge diagnosis:
#      diagnose_kegg_pathway_overlap()
#
# 3. Use FELLA to support KEGG bridge pathways:
#      run_fella_from_kegg_bridge()
#      get_fella_supported_kegg_bridge()
#
# 4. Use dual FELLA only as a stricter sensitivity check:
#      run_fella_dual_bridge_from_kegg_bridge()
#
# 5. Use pathview only for final selected pathways:
#      run_pathview_safe()
#
# ------------------------------------------------------------------------------
# Report-level caution
# ------------------------------------------------------------------------------
#
# Exact compound overlap supports direct metabolite-level validation.
# KEGG pathway overlap supports pathway-level convergence.
# FELLA-supported KEGG bridge supports network/pathway-level mechanistic linkage.
# These three evidence levels should be reported separately and not mixed.
# ------------------------------------------------------------------------------



# ==========================================================================
# metaInte data access helpers

metaFuns$get_metaInte_cache_dir <- function(x, name = "meta")
{
  dir_cache <- tryCatch(x$dir_cache, error = function(e) NULL)

  if (is.null(dir_cache) || !nzchar(as.character(dir_cache))) {
    dir_cache <- "tmp"
  }

  file.path(dir_cache, name)
}

metaFuns$assert_metaInte_data_sources <- function(x)
{
  if (!is(x, "job_metaInte")) {
    stop('!is(x, "job_metaInte").')
  }

  if (is.null(x$data_sources)) {
    stop("x$data_sources is NULL. Please create the object with `do_metaInte()` first.")
  }

  vec_required_mebo <- c(
    "data_commu_res", "data_diff_commu", "data_overall_score"
  )
  vec_required_metabo <- c("data_diff_report")

  vec_missing_mebo <- setdiff(vec_required_mebo, names(x$data_sources$mebocost))
  vec_missing_metabo <- setdiff(vec_required_metabo, names(x$data_sources$metaboDiff))

  if (length(vec_missing_mebo) > 0L) {
    stop(glue::glue(
      "Missing MEBOCOST data source(s): {paste(vec_missing_mebo, collapse = ', ')}."
    ))
  }

  if (length(vec_missing_metabo) > 0L) {
    stop(glue::glue(
      "Missing metabolomics data source(s): {paste(vec_missing_metabo, collapse = ', ')}."
    ))
  }

  invisible(TRUE)
}

metaFuns$resolve_metaInte_fc_sign <- function(x, check_compare = TRUE)
{
  if (!check_compare) {
    return(1)
  }

  metaFuns$assert_metaInte_data_sources(x)

  vec_mebo_level <- as.character(x$data_sources$mebocost$levels)
  case_group <- x$data_sources$metaboDiff$case_group
  control_group <- x$data_sources$metaboDiff$control_group

  if (length(vec_mebo_level) < 2L || is.null(case_group) || is.null(control_group)) {
    message("Skip comparison-group direction check: group labels were not fully detected.")
    return(1)
  }

  vec_ref_level <- c(as.character(case_group), as.character(control_group))

  if (identical(vec_mebo_level[1L:2L], vec_ref_level)) {
    return(1)
  }

  if (identical(vec_mebo_level[1L:2L], rev(vec_ref_level))) {
    message("Metabolomics contrast appears reversed relative to MEBOCOST; log2FC will be flipped.")
    return(-1)
  }

  message(glue::glue(
    "Comparison groups may be inconsistent: MEBOCOST = {paste(vec_mebo_level[1L:2L], collapse = ' vs ')}, ",
    "metabolomics = {paste(vec_ref_level, collapse = ' vs ')}."
  ))

  1
}

metaFuns$get_mebocost_metabolite_names_from_data <- function(data_commu_res = NULL,
  data_diff_commu = NULL, data_overall_score = NULL,
  source = c("commu_res", "diff_commu", "overall_score"))
{
  source <- match.arg(source)

  if (source == "commu_res") {
    data_mebo <- tibble::as_tibble(data_commu_res)
  } else if (source == "diff_commu") {
    data_mebo <- tibble::as_tibble(data_diff_commu)
  } else {
    data_mebo <- tibble::as_tibble(data_overall_score)
  }

  if (!"Metabolite_Name" %in% colnames(data_mebo)) {
    stop('"Metabolite_Name" was not found in MEBOCOST table.')
  }

  unique(data_mebo$Metabolite_Name)
}

metaFuns$get_mebocost_metabolite_universe_from_data <- function(data_commu_res = NULL,
  data_diff_commu = NULL, data_overall_score = NULL,
  sources = c("commu_res", "diff_commu", "overall_score"))
{
  sources <- match.arg(sources, several.ok = TRUE)
  lst_data <- list()

  if ("commu_res" %in% sources && !is.null(data_commu_res)) {
    data_commu <- tibble::as_tibble(data_commu_res)

    if ("Metabolite_Name" %in% colnames(data_commu)) {
      data_tmp <- dplyr::group_by(data_commu, Metabolite_Name)
      data_tmp <- dplyr::summarise(
        data_tmp,
        mebo_source = "step3_commu_res",
        n_commu_event = dplyr::n(),
        n_sender = dplyr::n_distinct(Sender),
        n_receiver = dplyr::n_distinct(Receiver),
        n_sensor = dplyr::n_distinct(Sensor),
        .groups = "drop"
      )
      lst_data[["step3_commu_res"]] <- data_tmp
    }
  }

  if ("diff_commu" %in% sources && !is.null(data_diff_commu)) {
    data_diff <- tibble::as_tibble(data_diff_commu)

    if ("Metabolite_Name" %in% colnames(data_diff)) {
      data_tmp <- dplyr::group_by(data_diff, Metabolite_Name)
      data_tmp <- dplyr::summarise(
        data_tmp,
        mebo_source = "step4_diff_commu",
        n_commu_event = dplyr::n(),
        n_sender = dplyr::n_distinct(Sender),
        n_receiver = dplyr::n_distinct(Receiver),
        n_sensor = dplyr::n_distinct(Sensor),
        .groups = "drop"
      )
      lst_data[["step4_diff_commu"]] <- data_tmp
    }
  }

  if ("overall_score" %in% sources && !is.null(data_overall_score)) {
    data_score <- tibble::as_tibble(data_overall_score)

    if ("Metabolite_Name" %in% colnames(data_score)) {
      data_tmp <- dplyr::group_by(data_score, Metabolite_Name)
      data_tmp <- dplyr::summarise(
        data_tmp,
        mebo_source = "step5_overall_score",
        n_commu_event = dplyr::n(),
        n_sender = NA_integer_,
        n_receiver = dplyr::n_distinct(Receiver),
        n_sensor = NA_integer_,
        max_overall_score = max(overall_score, na.rm = TRUE),
        .groups = "drop"
      )
      lst_data[["step5_overall_score"]] <- data_tmp
    }
  }

  data_out <- dplyr::bind_rows(lst_data)

  if (nrow(data_out) == 0L) {
    stop("No MEBOCOST metabolite universe was found from prepared data.")
  }

  data_out
}

metaFuns$diagnose_mebocost_metabo_coverage_from_data <- function(data_commu_res = NULL,
  data_diff_commu = NULL, data_overall_score = NULL, data_metabo_diff,
  sources = c("commu_res", "diff_commu", "overall_score"),
  p_cutoff = 0.05, vip_cutoff = NULL, padj_cutoff = NULL,
  dir_cache = "tmp", namespace = "name", dic_name = NULL,
  manual_cid_mebo = NULL, manual_cid_metabo = NULL)
{
  if (!dir.exists(dir_cache)) {
    dir.create(dir_cache, recursive = TRUE, showWarnings = FALSE)
  }

  data_mebo_universe <- metaFuns$get_mebocost_metabolite_universe_from_data(
    data_commu_res = data_commu_res,
    data_diff_commu = data_diff_commu,
    data_overall_score = data_overall_score,
    sources = sources
  )

  data_metabo_sig <- metaFuns$filter_metabo_significant(
    data_metabo_diff = data_metabo_diff,
    p_cutoff = p_cutoff,
    vip_cutoff = vip_cutoff,
    padj_cutoff = padj_cutoff
  )

  message(glue::glue(
    "MEBOCOST metabolite universe: {dplyr::n_distinct(data_mebo_universe$Metabolite_Name)} unique metabolite name(s)."
  ))
  message(glue::glue(
    "Metabolomics significant metabolites: {dplyr::n_distinct(data_metabo_sig$feature_name)} unique feature name(s)."
  ))

  data_mebo_cid <- metaFuns$map_pubchem_cids(
    unique(data_mebo_universe$Metabolite_Name),
    dir_cache = dir_cache,
    namespace = namespace,
    dic_name = dic_name,
    cache_name = "pubchemr_cids_mebocost_universe",
    manual_cid = manual_cid_mebo
  )

  data_metabo_cid <- metaFuns$map_pubchem_cids(
    unique(data_metabo_sig$feature_name),
    dir_cache = dir_cache,
    namespace = namespace,
    dic_name = dic_name,
    cache_name = "pubchemr_cids_metabolomics_significant",
    manual_cid = manual_cid_metabo
  )

  data_mebo_mapped <- dplyr::left_join(
    data_mebo_universe,
    data_mebo_cid,
    by = c("Metabolite_Name" = "name_original")
  )

  data_metabo_mapped <- dplyr::left_join(
    data_metabo_sig,
    data_metabo_cid,
    by = c("feature_name" = "name_original")
  )

  data_overlap_cid <- dplyr::inner_join(
    data_mebo_mapped,
    data_metabo_mapped,
    by = "cid",
    suffix = c("_mebocost", "_metabo")
  )

  data_mebo_name <- unique(metaFuns$normalize_metabolite_name(
    data_mebo_universe$Metabolite_Name,
    dic_name = dic_name
  ))

  data_metabo_name <- unique(metaFuns$normalize_metabolite_name(
    data_metabo_sig$feature_name,
    dic_name = dic_name
  ))

  vec_overlap_name <- intersect(data_mebo_name, data_metabo_name)

  diag_mebo <- metaFuns$mapping_diagnostics(data_mebo_cid, "MEBOCOST universe")
  diag_metabo <- metaFuns$mapping_diagnostics(data_metabo_cid, "significant metabolomics")

  data_summary_source <- dplyr::group_by(data_mebo_mapped, mebo_source)
  data_summary_source <- dplyr::summarise(
    data_summary_source,
    n_mebo_metabolite = dplyr::n_distinct(Metabolite_Name),
    n_mebo_mapped = dplyr::n_distinct(Metabolite_Name[!is.na(cid) & nzchar(cid)]),
    n_cid_overlap = dplyr::n_distinct(
      Metabolite_Name[cid %in% data_metabo_mapped$cid[!is.na(data_metabo_mapped$cid)]]
    ),
    .groups = "drop"
  )

  data_summary <- data.frame(
    n_mebo_metabolite = dplyr::n_distinct(data_mebo_universe$Metabolite_Name),
    n_metabo_sig = dplyr::n_distinct(data_metabo_sig$feature_name),
    n_mebo_mapped = diag_mebo$summary$n_mapped,
    n_metabo_mapped = diag_metabo$summary$n_mapped,
    n_overlap_cid = dplyr::n_distinct(data_overlap_cid$cid),
    n_overlap_name = length(vec_overlap_name),
    stringsAsFactors = FALSE
  )

  message(glue::glue(
    "Coverage diagnosis:\n",
    "  MEBOCOST metabolites: {data_summary$n_mebo_metabolite}; mapped: {data_summary$n_mebo_mapped}.\n",
    "  Significant metabolomics metabolites: {data_summary$n_metabo_sig}; mapped: {data_summary$n_metabo_mapped}.\n",
    "  CID overlap: {data_summary$n_overlap_cid}; normalized-name overlap: {data_summary$n_overlap_name}."
  ))

  list(
    summary = data_summary,
    summary_by_mebo_source = data_summary_source,
    data_mebo_universe = data_mebo_universe,
    data_metabo_sig = data_metabo_sig,
    data_mebo_mapped = data_mebo_mapped,
    data_metabo_mapped = data_metabo_mapped,
    data_overlap_cid = data_overlap_cid,
    overlap_name = vec_overlap_name,
    diag_mebo = diag_mebo,
    diag_metabo = diag_metabo,
    manual_diag_mebo = attr(data_mebo_cid, "manual_diag"),
    manual_diag_metabo = attr(data_metabo_cid, "manual_diag")
  )
}

metaFuns$diagnose_kegg_pathway_overlap_from_data <- function(data_commu_res = NULL,
  data_diff_commu = NULL, data_overall_score = NULL, data_metabo_diff,
  mebo_source = c("commu_res", "diff_commu", "overall_score"),
  p_cutoff = 0.05, vip_cutoff = NULL, padj_cutoff = NULL,
  dir_cache = "tmp", manual_kegg_mebo = NULL, manual_kegg_metabo = NULL)
{
  mebo_source <- match.arg(mebo_source)

  if (!dir.exists(dir_cache)) {
    dir.create(dir_cache, recursive = TRUE, showWarnings = FALSE)
  }

  vec_mebo_name <- metaFuns$get_mebocost_metabolite_names_from_data(
    data_commu_res = data_commu_res,
    data_diff_commu = data_diff_commu,
    data_overall_score = data_overall_score,
    source = mebo_source
  )

  data_metabo_sig <- metaFuns$filter_metabo_significant(
    data_metabo_diff = data_metabo_diff,
    p_cutoff = p_cutoff,
    vip_cutoff = vip_cutoff,
    padj_cutoff = padj_cutoff
  )

  vec_metabo_name <- unique(data_metabo_sig$feature_name)

  data_mebo_map <- metaFuns$map_names_to_kegg_compound(
    vec_mebo_name,
    dir_cache = dir_cache,
    manual_kegg = manual_kegg_mebo,
    cache_name = "kegg_compound_name_table"
  )

  data_metabo_map <- metaFuns$map_names_to_kegg_compound(
    vec_metabo_name,
    dir_cache = dir_cache,
    manual_kegg = manual_kegg_metabo,
    cache_name = "kegg_compound_name_table"
  )

  diag_mebo <- metaFuns$diagnose_kegg_mapping(data_mebo_map, "MEBOCOST metabolites")
  diag_metabo <- metaFuns$diagnose_kegg_mapping(data_metabo_map, "significant metabolomics metabolites")

  vec_mebo_kegg_id <- unique(as.character(data_mebo_map$kegg_id))
  vec_mebo_kegg_id <- vec_mebo_kegg_id[!is.na(vec_mebo_kegg_id) & nzchar(vec_mebo_kegg_id)]

  vec_metabo_kegg_id <- unique(as.character(data_metabo_map$kegg_id))
  vec_metabo_kegg_id <- vec_metabo_kegg_id[!is.na(vec_metabo_kegg_id) & nzchar(vec_metabo_kegg_id)]

  data_mebo_path <- expect_local_data(
    dir_cache,
    glue::glue("kegg_pathway_links_mebocost_{mebo_source}"),
    metaFuns$get_kegg_pathways_for_compounds,
    list(vec_kegg_id = vec_mebo_kegg_id)
  )

  data_metabo_path <- expect_local_data(
    dir_cache,
    "kegg_pathway_links_metabolomics_sig",
    metaFuns$get_kegg_pathways_for_compounds,
    list(vec_kegg_id = vec_metabo_kegg_id)
  )

  data_mebo_path <- dplyr::left_join(
    data_mebo_map,
    data_mebo_path,
    by = "kegg_id"
  )

  data_metabo_path <- dplyr::left_join(
    data_metabo_map,
    data_metabo_path,
    by = "kegg_id"
  )

  data_metabo_path <- dplyr::left_join(
    data_metabo_path,
    data_metabo_sig,
    by = c("name_original" = "feature_name")
  )

  vec_shared_pathway <- intersect(
    unique(data_mebo_path$pathway_id[!is.na(data_mebo_path$pathway_id)]),
    unique(data_metabo_path$pathway_id[!is.na(data_metabo_path$pathway_id)])
  )

  data_pathway_summary <- lapply(vec_shared_pathway, function(pathway) {
    data_mebo_sub <- data_mebo_path[data_mebo_path$pathway_id == pathway, , drop = FALSE]
    data_metabo_sub <- data_metabo_path[data_metabo_path$pathway_id == pathway, , drop = FALSE]

    data.frame(
      pathway_id = pathway,
      pathway_name = data_mebo_sub$pathway_name[which(!is.na(data_mebo_sub$pathway_name))[1L]],
      n_mebo_metabolite = metaFuns$n_distinct_non_na(data_mebo_sub$name_original),
      mebo_metabolites = metaFuns$collapse_non_na(data_mebo_sub$name_original),
      n_metabo_sig = metaFuns$n_distinct_non_na(data_metabo_sub$name_original),
      metabo_features = metaFuns$collapse_non_na(data_metabo_sub$name_original),
      best_metabo_pvalue = min(data_metabo_sub$pvalue, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  data_pathway_summary <- dplyr::bind_rows(data_pathway_summary)

  if (nrow(data_pathway_summary) > 0L) {
    data_pathway_summary <- dplyr::arrange(data_pathway_summary, best_metabo_pvalue)
  }

  message(glue::glue(
    "KEGG pathway overlap diagnosis: ",
    "MEBOCOST source = {mebo_source}; ",
    "shared KEGG pathways = {nrow(data_pathway_summary)}."
  ))

  structure(list(
    data_pathway_summary = data_pathway_summary,
    data_mebo_map = data_mebo_map,
    data_metabo_map = data_metabo_map,
    data_mebo_path = data_mebo_path,
    data_metabo_path = data_metabo_path,
    data_metabo_sig = data_metabo_sig,
    diag_mebo = diag_mebo,
    diag_metabo = diag_metabo
  ), class = "kegg_results_metaboDiff_and_mebocost")
}

# ==========================================================================
# pubchem

metaFuns$message_pubchem_mapping_diagnostics <- function(lst_refine, n_show = 20L)
{
  data_mebo_summary <- lst_refine$diag_mebo$summary
  data_metabo_summary <- lst_refine$diag_metabo$summary

  message(glue::glue(
    "\n================================================================\n",
    "PubChem CID mapping diagnostics:\n",
    "  MEBOCOST metabolites: {data_mebo_summary$n_name}; ",
    "mapped: {data_mebo_summary$n_mapped}; ",
    "unmapped: {data_mebo_summary$n_unmapped}; ",
    "multi-CID matched: {data_mebo_summary$n_multi_cid}.\n",
    "  Metabolomics metabolites: {data_metabo_summary$n_name}; ",
    "mapped: {data_metabo_summary$n_mapped}; ",
    "unmapped: {data_metabo_summary$n_unmapped}; ",
    "multi-CID matched: {data_metabo_summary$n_multi_cid}."
  ))

  data_mebo_unmapped <- lst_refine$diag_mebo$unmapped
  data_metabo_unmapped <- lst_refine$diag_metabo$unmapped
  data_mebo_multi <- lst_refine$diag_mebo$multi_cid
  data_metabo_multi <- lst_refine$diag_metabo$multi_cid

  if (nrow(data_mebo_unmapped) > 0L) {
    message(glue::glue(
      "Unmapped MEBOCOST metabolites, showing up to {n_show}: ",
      paste(utils::head(data_mebo_unmapped$name_original, n_show), collapse = "; ")
    ))
  }

  if (nrow(data_metabo_unmapped) > 0L) {
    message(glue::glue(
      "Unmapped metabolomics metabolites, showing up to {n_show}: ",
      paste(utils::head(data_metabo_unmapped$name_original, n_show), collapse = "; ")
    ))
  }

  if (nrow(data_mebo_multi) > 0L) {
    message(glue::glue(
      "MEBOCOST metabolites with multiple PubChem CIDs, showing up to {n_show}: ",
      paste(utils::head(data_mebo_multi$name_original, n_show), collapse = "; ")
    ))
  }

  if (nrow(data_metabo_multi) > 0L) {
    message(glue::glue(
      "Metabolomics metabolites with multiple PubChem CIDs, showing up to {n_show}: ",
      paste(utils::head(data_metabo_multi$name_original, n_show), collapse = "; ")
    ))
  }
  message("================================================================")
  invisible(lst_refine)
}

metaFuns$normalize_metabolite_name <- function(vec_name, dic_name = NULL)
{
  vec_name <- trimws(as.character(vec_name))

  if (!is.null(dic_name)) {
    vec_name <- ifelse(vec_name %in% names(dic_name), unname(dic_name[vec_name]), vec_name)
  }

  vec_name
}


metaFuns$extract_pubchem_cids <- function(object_cid, identifier = NULL)
{
  data_cid <- tryCatch(
    as.data.frame(PubChemR::retrieve(object_cid), stringsAsFactors = FALSE),
    error = function(e) NULL
  )

  if (is.null(data_cid) || nrow(data_cid) == 0L) {
    data_cid <- tryCatch(
      as.data.frame(PubChemR::CIDs(object_cid, .to.data.frame = TRUE),
        stringsAsFactors = FALSE),
      error = function(e) NULL
    )
  }

  if (is.null(data_cid) || nrow(data_cid) == 0L) {
    return(data.frame(
      name_query = character(0L),
      cid = character(0L),
      stringsAsFactors = FALSE
    ))
  }

  vec_col <- colnames(data_cid)
  vec_col_lower <- tolower(vec_col)

  idx_id <- grep("identifier|query|input|name", vec_col_lower)
  idx_cid <- grep("^cid$|cids|compound", vec_col_lower)

  if (length(idx_cid) == 0L) {
    stop("Cannot identify CID column from PubChemR result.")
  }

  col_cid <- vec_col[idx_cid[1L]]

  if (length(idx_id) > 0L) {
    col_id <- vec_col[idx_id[1L]]
    data_out <- data.frame(
      name_query = as.character(data_cid[[col_id]]),
      cid = data_cid[[col_cid]],
      stringsAsFactors = FALSE
    )
  } else {
    if (!is.null(identifier) && nrow(data_cid) == length(identifier)) {
      data_out <- data.frame(
        name_query = as.character(identifier),
        cid = data_cid[[col_cid]],
        stringsAsFactors = FALSE
      )
    } else {
      stop("Cannot identify identifier column from PubChemR result.")
    }
  }

  if (is.list(data_out$cid)) {
    data_out <- tidyr::unnest_longer(data_out, cid, keep_empty = TRUE)
  }

  data_out$cid <- as.character(data_out$cid)
  data_out$cid <- gsub("^CID:", "", data_out$cid)
  data_out$cid <- trimws(data_out$cid)

  data_out <- tidyr::separate_rows(data_out, cid, sep = "[,; ]+")
  data_out <- unique(data_out)

  data_out
}

metaFuns$get_cids_pubchem_batch <- function(identifier, namespace = "name",
  domain = "compound", searchtype = NULL, options = NULL)
{
  PubChemR::get_cids(
    identifier = identifier,
    namespace = namespace,
    domain = domain,
    searchtype = searchtype,
    options = options
  )
}

metaFuns$as_manual_cid_table <- function(manual_cid, dic_name = NULL)
{
  if (is.null(manual_cid)) {
    return(data.frame(
      name_original = character(0L),
      name_query = character(0L),
      cid = character(0L),
      stringsAsFactors = FALSE
    ))
  }

  if (is.vector(manual_cid) && !is.null(names(manual_cid))) {
    data_manual <- data.frame(
      name_original = names(manual_cid),
      cid = as.character(unname(manual_cid)),
      stringsAsFactors = FALSE
    )
  } else {
    data_manual <- as.data.frame(manual_cid, stringsAsFactors = FALSE)

    if (!all(c("name_original", "cid") %in% colnames(data_manual))) {
      stop('manual_cid must be a named vector or contain columns: "name_original", "cid".')
    }

    data_manual$name_original <- as.character(data_manual$name_original)
    data_manual$cid <- as.character(data_manual$cid)
  }

  data_manual$name_original <- trimws(data_manual$name_original)
  data_manual$cid <- trimws(data_manual$cid)

  data_manual <- tidyr::separate_rows(data_manual, cid, sep = "[,; ]+")
  data_manual <- data_manual[!is.na(data_manual$name_original) & nzchar(data_manual$name_original), ]
  data_manual <- data_manual[!is.na(data_manual$cid) & nzchar(data_manual$cid), ]

  data_manual$name_query <- metaFuns$normalize_metabolite_name(
    data_manual$name_original,
    dic_name = dic_name
  )

  unique(data_manual[, c("name_original", "name_query", "cid")])
}

metaFuns$apply_manual_cids <- function(data_map, manual_cid = NULL, dic_name = NULL)
{
  data_manual <- metaFuns$as_manual_cid_table(manual_cid, dic_name = dic_name)

  if (nrow(data_manual) == 0L) {
    attr(data_map, "manual_diag") <- data.frame(
      name_original = character(0L),
      cid = character(0L),
      matched_by = character(0L),
      stringsAsFactors = FALSE
    )

    return(metaFuns$recount_cid_n(data_map))
  }

  data_map$name_original <- as.character(data_map$name_original)
  data_map$name_query <- as.character(data_map$name_query)
  data_map$cid <- as.character(data_map$cid)

  lst_add <- lapply(seq_len(nrow(data_manual)), function(i) {
    name_i <- data_manual$name_original[i]
    query_i <- data_manual$name_query[i]
    cid_i <- data_manual$cid[i]

    idx_exact <- data_map$name_original == name_i
    idx_query <- data_map$name_query == query_i

    if (any(idx_exact, na.rm = TRUE)) {
      data_hit <- unique(data_map[idx_exact, c("name_original", "name_query")])
      data_hit$cid <- cid_i
      data_hit$cid_source <- "manual"
      data_hit$manual_name <- name_i
      data_hit$matched_by <- "name_original"
      return(data_hit)
    }

    if (any(idx_query, na.rm = TRUE)) {
      data_hit <- unique(data_map[idx_query, c("name_original", "name_query")])
      data_hit$cid <- cid_i
      data_hit$cid_source <- "manual"
      data_hit$manual_name <- name_i
      data_hit$matched_by <- "name_query"
      return(data_hit)
    }

    data.frame(
      name_original = name_i,
      name_query = query_i,
      cid = cid_i,
      cid_source = "manual",
      manual_name = name_i,
      matched_by = "not_found_in_input",
      stringsAsFactors = FALSE
    )
  })

  data_add <- dplyr::bind_rows(lst_add)

  data_map$cid_source <- ifelse(
    is.na(data_map$cid) | !nzchar(data_map$cid),
    NA_character_,
    "pubchem"
  )

  data_out <- dplyr::bind_rows(
    data_map[, c("name_original", "name_query", "cid", "cid_source")],
    data_add[, c("name_original", "name_query", "cid", "cid_source")]
  )

  data_out <- unique(data_out)
  data_out <- metaFuns$recount_cid_n(data_out)

  attr(data_out, "manual_diag") <- data_add[, c(
    "manual_name", "name_original", "name_query", "cid", "matched_by"
  )]

  data_out
}

metaFuns$message_manual_cid_diagnostics <- function(data_map, label = "metabolites")
{
  data_manual <- attr(data_map, "manual_diag")

  if (is.null(data_manual) || nrow(data_manual) == 0L) {
    message(glue::glue("No manual CID supplement was provided for {label}."))
    return(invisible(NULL))
  }

  message("++++++++++++++++++++++++++++++++++++++++++")
  message(glue::glue(
    "Manual CID supplement diagnostics for {label}: ",
    "{nrow(data_manual)} input CID record(s)."
  ))

  print(data_manual)

  n_not_found <- sum(data_manual$matched_by == "not_found_in_input", na.rm = TRUE)

  if (n_not_found > 0L) {
    message(glue::glue(
      "Warning: {n_not_found} manual CID record(s) did not match any input metabolite name. ",
      "Please check spelling, filtering, or whether the metabolite was removed before CID mapping."
    ))
  }
  message("++++++++++++++++++++++++++++++++++++++++++")
  invisible(data_manual)
}

metaFuns$recount_cid_n <- function(data_map)
{
  data_count <- dplyr::group_by(data_map, name_original)
  data_count <- dplyr::summarise(
    data_count,
    cid_n = sum(!is.na(cid) & nzchar(cid)),
    .groups = "drop"
  )

  data_map$cid_n <- NULL

  dplyr::left_join(data_map, data_count, by = "name_original")
}

metaFuns$mapping_diagnostics <- function(data_map, label = "metabolites")
{
  data_by_name <- dplyr::group_by(data_map, name_original)
  data_by_name <- dplyr::summarise(
    data_by_name,
    name_query = dplyr::first(name_query),
    cid_n = sum(!is.na(cid) & nzchar(cid)),
    cid_all = paste(unique(cid[!is.na(cid) & nzchar(cid)]), collapse = ";"),
    .groups = "drop"
  )

  data_unmapped <- dplyr::filter(data_by_name, cid_n == 0L)
  data_multi <- dplyr::filter(data_by_name, cid_n > 1L)

  data_summary <- data.frame(
    label = label,
    n_name = nrow(data_by_name),
    n_mapped = sum(data_by_name$cid_n > 0L),
    n_unmapped = nrow(data_unmapped),
    n_multi_cid = nrow(data_multi),
    stringsAsFactors = FALSE
  )

  list(
    summary = data_summary,
    unmapped = data_unmapped,
    multi_cid = data_multi
  )
}

metaFuns$map_pubchem_cids <- function(vec_name, dir_cache = "tmp",
  namespace = "name", domain = "compound", searchtype = NULL,
  options = NULL, dic_name = NULL, cache_name = NULL,
  manual_cid = NULL)
{
  expect_package("PubChemR", "3.0.0")

  vec_name <- unique(trimws(as.character(vec_name)))
  vec_name <- vec_name[!is.na(vec_name) & nzchar(vec_name)]

  data_query <- data.frame(
    name_original = vec_name,
    name_query = metaFuns$normalize_metabolite_name(vec_name, dic_name = dic_name),
    stringsAsFactors = FALSE
  )

  data_query$name_query <- trimws(data_query$name_query)
  data_query <- data_query[!is.na(data_query$name_query) & nzchar(data_query$name_query), ]

  vec_query <- unique(data_query$name_query)

  if (is.null(cache_name)) {
    cache_name <- glue::glue("pubchemr_cids_{namespace}_{domain}_{length(vec_query)}")
  }

  object_cid <- expect_local_data(
    dir_cache, cache_name, metaFuns$get_cids_pubchem_batch,
    list(
      identifier = vec_query,
      namespace = namespace,
      domain = domain,
      searchtype = searchtype,
      options = options
    )
  )

  data_cid <- metaFuns$extract_pubchem_cids(
    object_cid = object_cid,
    identifier = vec_query
  )

  data_out <- dplyr::left_join(
    data_query,
    data_cid,
    by = "name_query"
  )

  data_out$cid <- as.character(data_out$cid)
  data_out <- metaFuns$apply_manual_cids(
    data_out,
    manual_cid = manual_cid,
    dic_name = dic_name
  )
  data_out
}

# ==========================================================================
# intersect and integrated

metaFuns$scale01 <- function(x)
{
  x <- as.numeric(x)

  if (all(is.na(x))) {
    return(rep(NA_real_, length(x)))
  }

  min_x <- min(x, na.rm = TRUE)
  max_x <- max(x, na.rm = TRUE)

  if (isTRUE(all.equal(min_x, max_x))) {
    return(rep(1, length(x)))
  }

  (x - min_x) / (max_x - min_x)
}

metaFuns$get_job_value <- function(obj, names)
{
  for (name in names) {
    out <- tryCatch(obj[[name]], error = function(e) NULL)

    if (!is.null(out)) {
      return(out)
    }
  }

  NULL
}

metaFuns$resolve_metabo_fc_sign <- function(x, ref, check_compare = TRUE)
{
  if (!check_compare) {
    return(1)
  }

  vec_mebo_level <- tryCatch(as.character(x$levels), error = function(e) character(0L))
  case_group <- metaFuns$get_job_value(ref, c("case_group", "case.group", "case"))
  control_group <- metaFuns$get_job_value(ref, c("control_group", "control.group", "control"))

  if (length(vec_mebo_level) < 2L || is.null(case_group) || is.null(control_group)) {
    message("Skip comparison-group direction check: group labels were not fully detected.")
    return(1)
  }

  vec_ref_level <- c(as.character(case_group), as.character(control_group))

  if (identical(vec_mebo_level[1L:2L], vec_ref_level)) {
    return(1)
  }

  if (identical(vec_mebo_level[1L:2L], rev(vec_ref_level))) {
    message("Metabolomics contrast appears reversed relative to MEBOCOST; log2FC will be flipped.")
    return(-1)
  }

  message(glue::glue(
    "Comparison groups may be inconsistent: MEBOCOST = {paste(vec_mebo_level[1L:2L], collapse = ' vs ')}, ",
    "metabolomics = {paste(vec_ref_level, collapse = ' vs ')}."
  ))

  1
}

metaFuns$summarise_mebocost_direction <- function(data_diff_commu,
  group_by = c("Metabolite_Name", "Receiver"), score_col = NULL)
{
  data_diff_commu <- tibble::as_tibble(data_diff_commu)

  if (!"Log2FC" %in% colnames(data_diff_commu)) {
    stop('"Log2FC" is not in data_diff_commu.')
  }

  if (!is.null(score_col) && score_col %in% colnames(data_diff_commu)) {
    data_diff_commu$direction_weight <- abs(as.numeric(data_diff_commu[[score_col]]))
  } else {
    data_diff_commu$direction_weight <- abs(as.numeric(data_diff_commu$Log2FC))
  }

  data_diff_commu$direction_weight[is.na(data_diff_commu$direction_weight)] <- 0
  data_diff_commu$direction_weight[data_diff_commu$direction_weight == 0] <- 1

  data_direction <- dplyr::group_by(data_diff_commu, !!!rlang::syms(group_by))
  data_direction <- dplyr::summarise(
    data_direction,
    mebo_log2FC = stats::weighted.mean(Log2FC, direction_weight, na.rm = TRUE),
    mebo_direction = sign(mebo_log2FC),
    mebo_n_chain = dplyr::n(),
    .groups = "drop"
  )

  data_direction
}

metaFuns$default_metabo_filter <- function(data_metabo_diff,
  p_cutoff = 0.05, vip_cutoff = 1, padj_cutoff = NULL)
{
  data_metabo_filter <- tibble::as_tibble(data_metabo_diff)

  if ("pvalue" %in% colnames(data_metabo_filter)) {
    data_metabo_filter <- data_metabo_filter[
      !is.na(data_metabo_filter$pvalue) &
        data_metabo_filter$pvalue < p_cutoff,
    ]
  }

  if ("VIP" %in% colnames(data_metabo_filter)) {
    data_metabo_filter <- data_metabo_filter[
      !is.na(data_metabo_filter$VIP) &
        data_metabo_filter$VIP >= vip_cutoff,
    ]
  }

  if (!is.null(padj_cutoff) && "padj" %in% colnames(data_metabo_filter)) {
    data_metabo_filter <- data_metabo_filter[
      !is.na(data_metabo_filter$padj) &
        data_metabo_filter$padj < padj_cutoff,
    ]
  }

  data_metabo_filter
}

metaFuns$integrate_mebocost_metabolomics <- function(data_mebo_score,
  data_diff_commu, data_metabo_diff, dir_cache = "tmp",
  p_cutoff = 0.05, vip_cutoff = 1, padj_cutoff = NULL,
  namespace = "name", dic_name = NULL, use_direction = FALSE,
  direction_mode = c("same", "opposite"), fallback_if_empty = TRUE,
  metabo_filter_fun = NULL, metabo_fc_sign = 1,
  manual_cid_mebo = NULL, manual_cid_metabo = NULL
)
{
  direction_mode <- match.arg(direction_mode)

  data_mebo_score <- tibble::as_tibble(data_mebo_score)
  data_metabo_diff <- tibble::as_tibble(data_metabo_diff)

  if (!all(c("Metabolite_Name", "Receiver", "overall_score") %in% colnames(data_mebo_score))) {
    stop('data_mebo_score must contain "Metabolite_Name", "Receiver", and "overall_score".')
  }

  if (!all(c("feature_name", "log2FC") %in% colnames(data_metabo_diff))) {
    stop('data_metabo_diff must contain "feature_name" and "log2FC".')
  }

  if (is.null(metabo_filter_fun)) {
    data_metabo_filter <- metaFuns$default_metabo_filter(
      data_metabo_diff,
      p_cutoff = p_cutoff,
      vip_cutoff = vip_cutoff,
      padj_cutoff = padj_cutoff
    )
  } else {
    data_metabo_filter <- metabo_filter_fun(data_metabo_diff)
  }

  data_metabo_filter$metabo_log2FC <- as.numeric(data_metabo_filter$log2FC) * metabo_fc_sign
  data_metabo_filter$metabo_direction <- sign(data_metabo_filter$metabo_log2FC)

  vec_mebo_name <- unique(data_mebo_score$Metabolite_Name)
  vec_metabo_name <- unique(data_metabo_filter$feature_name)

  message(glue::glue("Map MEBOCOST metabolites: {length(vec_mebo_name)} names."))
  data_mebo_cid <- metaFuns$map_pubchem_cids(
    vec_mebo_name,
    dir_cache = dir_cache,
    namespace = namespace,
    dic_name = dic_name,
    cache_name = "pubchemr_cids_mebocost",
    manual_cid = manual_cid_mebo
  )
  diag_mebo <- metaFuns$mapping_diagnostics(data_mebo_cid, "MEBOCOST")
  metaFuns$message_manual_cid_diagnostics(data_mebo_cid, "MEBOCOST metabolites")

  message(glue::glue("Map metabolomics features: {length(vec_metabo_name)} names."))
  data_metabo_cid <- metaFuns$map_pubchem_cids(
    vec_metabo_name,
    dir_cache = dir_cache,
    namespace = namespace,
    dic_name = dic_name,
    cache_name = "pubchemr_cids_metabolomics",
    manual_cid = manual_cid_metabo
  )
  diag_metabo <- metaFuns$mapping_diagnostics(data_metabo_cid, "metabolomics")
  metaFuns$message_manual_cid_diagnostics(data_metabo_cid, "metabolomics metabolites")


  data_direction <- metaFuns$summarise_mebocost_direction(
    data_diff_commu,
    group_by = c("Metabolite_Name", "Receiver")
  )

  data_mebo <- dplyr::left_join(
    data_mebo_score,
    data_direction,
    by = c("Metabolite_Name", "Receiver")
  )

  data_mebo <- dplyr::left_join(
    data_mebo,
    data_mebo_cid,
    by = c("Metabolite_Name" = "name_original")
  )

  data_metabo <- dplyr::left_join(
    data_metabo_filter,
    data_metabo_cid,
    by = c("feature_name" = "name_original")
  )

  data_overlap <- dplyr::inner_join(
    data_mebo,
    data_metabo,
    by = "cid",
    suffix = c("_mebocost", "_metabo")
  )

  if (nrow(data_overlap) == 0L) {
    return(list(
      data_overlap = data_overlap,
      data_final = data_overlap,
      data_mebo_mapped = data_mebo,
      data_metabo_mapped = data_metabo,
      data_metabo_filter = data_metabo_filter,
      direction_used = FALSE,
      diag_mebo = diag_mebo,
      diag_metabo = diag_metabo,
      direction_fallback = FALSE
    ))
  }

  data_overlap$direction_match <- NA

  idx_valid_direction <- !is.na(data_overlap$mebo_direction) &
    !is.na(data_overlap$metabo_direction) &
    data_overlap$mebo_direction != 0 &
    data_overlap$metabo_direction != 0

  if (direction_mode == "same") {
    data_overlap$direction_match[idx_valid_direction] <-
      data_overlap$mebo_direction[idx_valid_direction] ==
        data_overlap$metabo_direction[idx_valid_direction]
  } else {
    data_overlap$direction_match[idx_valid_direction] <-
      data_overlap$mebo_direction[idx_valid_direction] !=
        data_overlap$metabo_direction[idx_valid_direction]
  }

  data_overlap$metabo_strength <- abs(data_overlap$metabo_log2FC)

  if ("VIP" %in% colnames(data_overlap)) {
    data_overlap$metabo_strength <- data_overlap$metabo_strength * data_overlap$VIP
  }

  if ("pvalue" %in% colnames(data_overlap)) {
    data_overlap$metabo_strength <- data_overlap$metabo_strength *
      (-log10(pmax(data_overlap$pvalue, .Machine$double.xmin)))
  }

  data_overlap$score_mebo_scaled <- metaFuns$scale01(data_overlap$overall_score)
  data_overlap$score_metabo_scaled <- metaFuns$scale01(data_overlap$metabo_strength)
  data_overlap$integrated_score <- data_overlap$score_mebo_scaled *
    data_overlap$score_metabo_scaled

  data_overlap <- dplyr::arrange(
    data_overlap,
    dplyr::desc(integrated_score)
  )

  if (use_direction) {
    data_direction_match <- data_overlap[
      !is.na(data_overlap$direction_match) & data_overlap$direction_match,
    ]

    if (nrow(data_direction_match) > 0L) {
      data_final <- data_direction_match
      direction_used <- TRUE
      direction_fallback <- FALSE
    } else if (fallback_if_empty) {
      data_final <- data_overlap
      direction_used <- FALSE
      direction_fallback <- TRUE
    } else {
      data_final <- data_direction_match
      direction_used <- TRUE
      direction_fallback <- FALSE
    }
  } else {
    data_final <- data_overlap
    direction_used <- FALSE
    direction_fallback <- FALSE
  }

  data_final$refine_rule <- if (direction_used) {
    paste0("CID overlap + direction ", direction_mode)
  } else if (direction_fallback) {
    "CID overlap fallback because no direction-matched result"
  } else {
    "CID overlap"
  }

  list(
    data_overlap = data_overlap,
    data_final = data_final,
    data_mebo_mapped = data_mebo,
    data_metabo_mapped = data_metabo,
    data_metabo_filter = data_metabo_filter,
    diag_mebo = diag_mebo,
    diag_metabo = diag_metabo,
    direction_used = direction_used,
    direction_fallback = direction_fallback
  )
}

# ==========================================================================
# check whether intersect

metaFuns$get_mebocost_metabolite_universe <- function(x,
  sources = c("commu_res", "diff_commu", "overall_score"))
{
  sources <- match.arg(sources, several.ok = TRUE)

  lst_data <- list()

  if ("commu_res" %in% sources) {
    data_commu <- tryCatch(
      tibble::as_tibble(x@tables$step3$t.commu_res),
      error = function(e) NULL
    )

    if (!is.null(data_commu) && "Metabolite_Name" %in% colnames(data_commu)) {
      data_tmp <- dplyr::group_by(data_commu, Metabolite_Name)
      data_tmp <- dplyr::summarise(
        data_tmp,
        mebo_source = "step3_commu_res",
        n_commu_event = dplyr::n(),
        n_sender = dplyr::n_distinct(Sender),
        n_receiver = dplyr::n_distinct(Receiver),
        n_sensor = dplyr::n_distinct(Sensor),
        .groups = "drop"
      )
      lst_data[["step3_commu_res"]] <- data_tmp
    }
  }

  if ("diff_commu" %in% sources) {
    data_diff <- tryCatch(
      tibble::as_tibble(x@tables$step4$ts.diff_commu[[1L]]),
      error = function(e) NULL
    )

    if (!is.null(data_diff) && "Metabolite_Name" %in% colnames(data_diff)) {
      data_tmp <- dplyr::group_by(data_diff, Metabolite_Name)
      data_tmp <- dplyr::summarise(
        data_tmp,
        mebo_source = "step4_diff_commu",
        n_commu_event = dplyr::n(),
        n_sender = dplyr::n_distinct(Sender),
        n_receiver = dplyr::n_distinct(Receiver),
        n_sensor = dplyr::n_distinct(Sensor),
        .groups = "drop"
      )
      lst_data[["step4_diff_commu"]] <- data_tmp
    }
  }

  if ("overall_score" %in% sources) {
    data_score <- tryCatch(
      tibble::as_tibble(x@tables$step5$t.overallScore),
      error = function(e) NULL
    )

    if (!is.null(data_score) && "Metabolite_Name" %in% colnames(data_score)) {
      data_tmp <- dplyr::group_by(data_score, Metabolite_Name)
      data_tmp <- dplyr::summarise(
        data_tmp,
        mebo_source = "step5_overall_score",
        n_commu_event = dplyr::n(),
        n_sender = NA_integer_,
        n_receiver = dplyr::n_distinct(Receiver),
        n_sensor = NA_integer_,
        max_overall_score = max(overall_score, na.rm = TRUE),
        .groups = "drop"
      )
      lst_data[["step5_overall_score"]] <- data_tmp
    }
  }

  data_out <- dplyr::bind_rows(lst_data)

  if (nrow(data_out) == 0L) {
    stop("No MEBOCOST metabolite universe was found.")
  }

  data_out
}

metaFuns$filter_metabo_significant <- function(data_metabo_diff,
  p_cutoff = 0.05, vip_cutoff = NULL, padj_cutoff = NULL)
{
  data_metabo <- tibble::as_tibble(data_metabo_diff)

  if (!"feature_name" %in% colnames(data_metabo)) {
    stop('"feature_name" is not in metabolomics differential table.')
  }

  if (!is.null(p_cutoff)) {
    if (!"pvalue" %in% colnames(data_metabo)) {
      stop('"pvalue" is not in metabolomics differential table.')
    }

    data_metabo <- data_metabo[
      !is.na(data_metabo$pvalue) &
        data_metabo$pvalue < p_cutoff,
    ]
  }

  if (!is.null(vip_cutoff)) {
    if (!"VIP" %in% colnames(data_metabo)) {
      stop('"VIP" is not in metabolomics differential table.')
    }

    data_metabo <- data_metabo[
      !is.na(data_metabo$VIP) &
        data_metabo$VIP >= vip_cutoff,
    ]
  }

  if (!is.null(padj_cutoff)) {
    if (!"padj" %in% colnames(data_metabo)) {
      stop('"padj" is not in metabolomics differential table.')
    }

    data_metabo <- data_metabo[
      !is.na(data_metabo$padj) &
        data_metabo$padj < padj_cutoff,
    ]
  }

  data_metabo
}

metaFuns$diagnose_mebocost_metabo_coverage <- function(x, ref,
  sources = c("commu_res", "diff_commu", "overall_score"),
  p_cutoff = 0.05, vip_cutoff = NULL, padj_cutoff = NULL,
  dir_cache = NULL, namespace = "name", dic_name = NULL,
  manual_cid_mebo = NULL, manual_cid_metabo = NULL)
{
  if (is.null(dir_cache)) {
    dir_cache <- file.path(x$dir_cache, "pubchem")
  }

  if (!dir.exists(dir_cache)) {
    dir.create(dir_cache, recursive = TRUE, showWarnings = FALSE)
  }

  data_mebo_universe <- metaFuns$get_mebocost_metabolite_universe(
    x = x,
    sources = sources
  )

  data_metabo_sig <- metaFuns$filter_metabo_significant(
    data_metabo_diff = ref@tables$step3$data_diff_report,
    p_cutoff = p_cutoff,
    vip_cutoff = vip_cutoff,
    padj_cutoff = padj_cutoff
  )

  message(glue::glue(
    "MEBOCOST metabolite universe: {dplyr::n_distinct(data_mebo_universe$Metabolite_Name)} unique metabolite name(s)."
  ))

  message(glue::glue(
    "Metabolomics significant metabolites: {dplyr::n_distinct(data_metabo_sig$feature_name)} unique feature name(s)."
  ))

  data_mebo_cid <- metaFuns$map_pubchem_cids(
    unique(data_mebo_universe$Metabolite_Name),
    dir_cache = dir_cache,
    namespace = namespace,
    dic_name = dic_name,
    cache_name = "pubchemr_cids_mebocost_universe",
    manual_cid = manual_cid_mebo
  )

  data_metabo_cid <- metaFuns$map_pubchem_cids(
    unique(data_metabo_sig$feature_name),
    dir_cache = dir_cache,
    namespace = namespace,
    dic_name = dic_name,
    cache_name = "pubchemr_cids_metabolomics_significant",
    manual_cid = manual_cid_metabo
  )

  data_mebo_mapped <- dplyr::left_join(
    data_mebo_universe,
    data_mebo_cid,
    by = c("Metabolite_Name" = "name_original")
  )

  data_metabo_mapped <- dplyr::left_join(
    data_metabo_sig,
    data_metabo_cid,
    by = c("feature_name" = "name_original")
  )

  data_overlap_cid <- dplyr::inner_join(
    data_mebo_mapped,
    data_metabo_mapped,
    by = "cid",
    suffix = c("_mebocost", "_metabo")
  )

  data_mebo_name <- unique(metaFuns$normalize_metabolite_name(
    data_mebo_universe$Metabolite_Name,
    dic_name = dic_name
  ))

  data_metabo_name <- unique(metaFuns$normalize_metabolite_name(
    data_metabo_sig$feature_name,
    dic_name = dic_name
  ))

  vec_overlap_name <- intersect(data_mebo_name, data_metabo_name)

  diag_mebo <- metaFuns$mapping_diagnostics(data_mebo_cid, "MEBOCOST universe")
  diag_metabo <- metaFuns$mapping_diagnostics(data_metabo_cid, "significant metabolomics")

  data_summary_source <- dplyr::group_by(data_mebo_mapped, mebo_source)
  data_summary_source <- dplyr::summarise(
    data_summary_source,
    n_mebo_metabolite = dplyr::n_distinct(Metabolite_Name),
    n_mebo_mapped = dplyr::n_distinct(Metabolite_Name[!is.na(cid) & nzchar(cid)]),
    n_cid_overlap = dplyr::n_distinct(
      Metabolite_Name[cid %in% data_metabo_mapped$cid[!is.na(data_metabo_mapped$cid)]]
    ),
    .groups = "drop"
  )

  data_summary <- data.frame(
    n_mebo_metabolite = dplyr::n_distinct(data_mebo_universe$Metabolite_Name),
    n_metabo_sig = dplyr::n_distinct(data_metabo_sig$feature_name),
    n_mebo_mapped = diag_mebo$summary$n_mapped,
    n_metabo_mapped = diag_metabo$summary$n_mapped,
    n_overlap_cid = dplyr::n_distinct(data_overlap_cid$cid),
    n_overlap_name = length(vec_overlap_name),
    stringsAsFactors = FALSE
  )

  message(glue::glue(
    "Coverage diagnosis:\n",
    "  MEBOCOST metabolites: {data_summary$n_mebo_metabolite}; mapped: {data_summary$n_mebo_mapped}.\n",
    "  Significant metabolomics metabolites: {data_summary$n_metabo_sig}; mapped: {data_summary$n_metabo_mapped}.\n",
    "  CID overlap: {data_summary$n_overlap_cid}; normalized-name overlap: {data_summary$n_overlap_name}."
  ))

  if (nrow(data_overlap_cid) == 0L) {
    message("No CID overlap was found between significant metabolomics metabolites and the selected MEBOCOST metabolite universe.")
  } else {
    message("CID overlap was found. You can inspect `$data_overlap_cid`.")
  }

  list(
    summary = data_summary,
    summary_by_mebo_source = data_summary_source,
    data_mebo_universe = data_mebo_universe,
    data_metabo_sig = data_metabo_sig,
    data_mebo_mapped = data_mebo_mapped,
    data_metabo_mapped = data_metabo_mapped,
    data_overlap_cid = data_overlap_cid,
    overlap_name = vec_overlap_name,
    diag_mebo = diag_mebo,
    diag_metabo = diag_metabo,
    manual_diag_mebo = attr(data_mebo_cid, "manual_diag"),
    manual_diag_metabo = attr(data_metabo_cid, "manual_diag")
  )
}

# ==========================================================================
# axis

metaFuns$default_mebocost_axis_dictionary <- function()
{
  list(
    Eicosanoid_Arachidonic_Acid_Axis = c(
      "arachidonic", "prostaglandin", "leukotriene", "thromboxane",
      "hydroxyeicos", "eicos", "hpete", "hete", "epoxy", "epoxyeicos"
    ),
    Glutamine_Amino_Acid_Axis = c(
      "glutamine", "glutamate", "glutamic", "aspartate", "alanine",
      "arginine", "ornithine", "citrulline", "amino"
    ),
    Purine_Nucleotide_Axis = c(
      "adenosine", "adenine", "amp", "adp", "atp", "inosine",
      "hypoxanthine", "xanthine", "uric", "purine"
    ),
    Cholesterol_Oxysterol_Sterol_Axis = c(
      "cholesterol", "hydroxycholesterol", "oxysterol", "sterol",
      "bile acid", "chenodeoxycholic", "cholic", "deoxycholic"
    ),
    Heme_Iron_Porphyrin_Axis = c(
      "heme", "iron", "porphyrin", "aminolevulinic", "biliverdin",
      "bilirubin", "protoporphyrin"
    ),
    Sphingolipid_Axis = c(
      "sphingosine", "sphingomyelin", "ceramide", "sphingolipid"
    ),
    Choline_Lipid_Axis = c(
      "choline", "phosphocholine", "glycerophosphocholine",
      "phosphatidylcholine", "lysopc", "lysope", "phosphatidylethanolamine"
    )
  )
}

metaFuns$match_metabolite_axis <- function(vec_name,
  axis_dictionary = metaFuns$default_mebocost_axis_dictionary())
{
  vec_name <- unique(trimws(as.character(vec_name)))
  vec_name <- vec_name[!is.na(vec_name) & nzchar(vec_name)]
  vec_lower <- tolower(vec_name)

  lst_hit <- lapply(names(axis_dictionary), function(axis) {
    pattern <- paste0(axis_dictionary[[axis]], collapse = "|")
    hit <- grepl(pattern, vec_lower, ignore.case = TRUE)

    if (!any(hit)) {
      return(data.frame(
        feature_name = character(0L),
        metabolic_axis = character(0L),
        stringsAsFactors = FALSE
      ))
    }

    data.frame(
      feature_name = vec_name[hit],
      metabolic_axis = rep(axis, sum(hit)),
      stringsAsFactors = FALSE
    )
  })

  data_out <- dplyr::bind_rows(lst_hit)
  dplyr::distinct(data_out)
}

metaFuns$diagnose_axis_level_support <- function(x, ref,
  p_cutoff = 0.05, vip_cutoff = NULL, padj_cutoff = NULL,
  axis_dictionary = metaFuns$default_mebocost_axis_dictionary())
{
  data_metabo_sig <- metaFuns$filter_metabo_significant(
    data_metabo_diff = ref@tables$step3$data_diff_report,
    p_cutoff = p_cutoff,
    vip_cutoff = vip_cutoff,
    padj_cutoff = padj_cutoff
  )

  data_mebo_step5 <- tibble::as_tibble(x@tables$step5$t.overallScore)

  data_mebo_axis <- metaFuns$match_metabolite_axis(
    data_mebo_step5$Metabolite_Name,
    axis_dictionary = axis_dictionary
  )

  data_metabo_axis <- metaFuns$match_metabolite_axis(
    data_metabo_sig$feature_name,
    axis_dictionary = axis_dictionary
  )

  data_mebo_axis <- dplyr::left_join(
    data_mebo_axis,
    data_mebo_step5,
    by = c("feature_name" = "Metabolite_Name")
  )

  data_metabo_axis <- dplyr::left_join(
    data_metabo_axis,
    data_metabo_sig,
    by = "feature_name"
  )

  if (nrow(data_mebo_axis) == 0L || nrow(data_metabo_axis) == 0L) {
    message("Axis-level support diagnosis: no metabolic-axis match was found.")

    return(list(
      data_axis_summary = data.frame(),
      data_mebo_axis = data_mebo_axis,
      data_metabo_axis = data_metabo_axis,
      data_metabo_sig = data_metabo_sig
    ))
  }

  vec_axis_shared <- intersect(
    unique(data_mebo_axis$metabolic_axis),
    unique(data_metabo_axis$metabolic_axis)
  )

  if (length(vec_axis_shared) == 0L) {
    message("Axis-level support diagnosis: no shared metabolic axis was found.")

    return(list(
      data_axis_summary = data.frame(),
      data_mebo_axis = data_mebo_axis,
      data_metabo_axis = data_metabo_axis,
      data_metabo_sig = data_metabo_sig
    ))
  }

  data_axis_summary <- lapply(vec_axis_shared, function(axis) {
    data_mebo_sub <- data_mebo_axis[data_mebo_axis$metabolic_axis == axis, , drop = FALSE]
    data_metabo_sub <- data_metabo_axis[data_metabo_axis$metabolic_axis == axis, , drop = FALSE]

    data.frame(
      metabolic_axis = axis,
      n_mebo_metabolite = dplyr::n_distinct(data_mebo_sub$feature_name),
      n_mebo_receiver_axis = nrow(data_mebo_sub),
      n_metabo_sig = dplyr::n_distinct(data_metabo_sub$feature_name),
      best_metabo_pvalue = min(data_metabo_sub$pvalue, na.rm = TRUE),
      max_mebo_overall_score = max(data_mebo_sub$overall_score, na.rm = TRUE),
      top_mebo_metabolite = data_mebo_sub$feature_name[
        which.max(data_mebo_sub$overall_score)
      ],
      top_metabo_feature = data_metabo_sub$feature_name[
        which.min(data_metabo_sub$pvalue)
      ],
      stringsAsFactors = FALSE
    )
  })

  data_axis_summary <- dplyr::bind_rows(data_axis_summary)
  data_axis_summary <- dplyr::arrange(
    data_axis_summary,
    best_metabo_pvalue,
    dplyr::desc(max_mebo_overall_score)
  )

  message(glue::glue(
    "Axis-level support diagnosis: {nrow(data_axis_summary)} shared metabolic axis/axes were found."
  ))

  list(
    data_axis_summary = data_axis_summary,
    data_mebo_axis = data_mebo_axis,
    data_metabo_axis = data_metabo_axis,
    data_metabo_sig = data_metabo_sig
  )
}

# ==========================================================================
# kegg pathway

metaFuns$normalize_kegg_name <- function(vec_name)
{
  vec_name <- trimws(as.character(vec_name))
  vec_name <- tolower(vec_name)
  vec_name <- gsub("α", "alpha", vec_name, fixed = TRUE)
  vec_name <- gsub("β", "beta", vec_name, fixed = TRUE)
  vec_name <- gsub("γ", "gamma", vec_name, fixed = TRUE)
  vec_name <- gsub("\\s+", " ", vec_name)
  vec_name <- gsub("[‘’`]", "'", vec_name)
  vec_name
}

metaFuns$get_kegg_compound_name_table <- function()
{
  vec_cpd <- KEGGREST::keggList("compound")

  data_cpd <- data.frame(
    kegg_id = names(vec_cpd),
    kegg_names = as.character(vec_cpd),
    stringsAsFactors = FALSE
  )

  data_cpd$kegg_id <- sub("^cpd:", "", data_cpd$kegg_id)

  data_cpd <- tidyr::separate_rows(data_cpd, kegg_names, sep = ";")
  data_cpd$kegg_name <- trimws(data_cpd$kegg_names)
  data_cpd$kegg_names <- NULL
  data_cpd$kegg_name_norm <- metaFuns$normalize_kegg_name(data_cpd$kegg_name)

  data_cpd <- data_cpd[!is.na(data_cpd$kegg_name_norm) & nzchar(data_cpd$kegg_name_norm), ]
  dplyr::distinct(data_cpd)
}

metaFuns$as_manual_kegg_table <- function(manual_kegg = NULL)
{
  if (is.null(manual_kegg)) {
    return(data.frame(
      name_original = character(0L),
      kegg_id = character(0L),
      stringsAsFactors = FALSE
    ))
  }

  if (is.vector(manual_kegg) && !is.null(names(manual_kegg))) {
    data_manual <- data.frame(
      name_original = names(manual_kegg),
      kegg_id = as.character(unname(manual_kegg)),
      stringsAsFactors = FALSE
    )
  } else {
    data_manual <- as.data.frame(manual_kegg, stringsAsFactors = FALSE)

    if (!all(c("name_original", "kegg_id") %in% colnames(data_manual))) {
      stop('manual_kegg must be a named vector or contain columns: "name_original", "kegg_id".')
    }
  }

  data_manual$name_original <- trimws(as.character(data_manual$name_original))
  data_manual$kegg_id <- trimws(as.character(data_manual$kegg_id))
  data_manual$kegg_id <- sub("^cpd:", "", data_manual$kegg_id)

  data_manual <- data_manual[!is.na(data_manual$name_original) & nzchar(data_manual$name_original), ]
  data_manual <- data_manual[!is.na(data_manual$kegg_id) & nzchar(data_manual$kegg_id), ]

  dplyr::distinct(data_manual)
}

metaFuns$map_names_to_kegg_compound <- function(vec_name, dir_cache = "tmp",
  manual_kegg = NULL, cache_name = "kegg_compound_name_table")
{
  vec_name <- unique(trimws(as.character(vec_name)))
  vec_name <- vec_name[!is.na(vec_name) & nzchar(vec_name)]

  data_input <- data.frame(
    name_original = vec_name,
    name_norm = metaFuns$normalize_kegg_name(vec_name),
    stringsAsFactors = FALSE
  )

  data_kegg <- expect_local_data(
    dir_cache, cache_name, metaFuns$get_kegg_compound_name_table, list()
  )

  data_map <- dplyr::left_join(
    data_input,
    data_kegg,
    by = c("name_norm" = "kegg_name_norm")
  )

  data_map$map_source <- ifelse(is.na(data_map$kegg_id), NA_character_, "kegg_name_exact")

  data_manual <- metaFuns$as_manual_kegg_table(manual_kegg)

  if (nrow(data_manual) > 0L) {
    data_manual$name_norm <- metaFuns$normalize_kegg_name(data_manual$name_original)
    data_manual$kegg_name <- data_manual$name_original
    data_manual$map_source <- "manual"

    data_map <- data_map[!data_map$name_original %in% data_manual$name_original, ]
    data_map <- dplyr::bind_rows(
      data_map,
      data_manual[, c("name_original", "name_norm", "kegg_id", "kegg_name", "map_source")]
    )
  }

  data_count <- dplyr::group_by(data_map, name_original)
  data_count <- dplyr::summarise(
    data_count,
    kegg_n = dplyr::n_distinct(kegg_id[!is.na(kegg_id) & nzchar(kegg_id)]),
    .groups = "drop"
  )

  data_map$kegg_n <- NULL
  data_map <- dplyr::left_join(data_map, data_count, by = "name_original")
  dplyr::distinct(data_map)
}

metaFuns$diagnose_kegg_mapping <- function(data_map, label = "metabolites", n_show = 20L)
{
  data_by_name <- dplyr::group_by(data_map, name_original)
  data_by_name <- dplyr::summarise(
    data_by_name,
    kegg_n = dplyr::n_distinct(kegg_id[!is.na(kegg_id) & nzchar(kegg_id)]),
    kegg_all = paste(unique(kegg_id[!is.na(kegg_id) & nzchar(kegg_id)]), collapse = ";"),
    .groups = "drop"
  )

  data_unmapped <- dplyr::filter(data_by_name, kegg_n == 0L)
  data_multi <- dplyr::filter(data_by_name, kegg_n > 1L)

  message(glue::glue(
    "KEGG compound mapping diagnostics for {label}: ",
    "input names = {nrow(data_by_name)}, ",
    "mapped = {sum(data_by_name$kegg_n > 0L)}, ",
    "unmapped = {nrow(data_unmapped)}, ",
    "multi-mapped = {nrow(data_multi)}."
  ))

  if (nrow(data_unmapped) > 0L) {
    message(glue::glue(
      "Unmapped {label}, showing up to {n_show}: ",
      paste(utils::head(data_unmapped$name_original, n_show), collapse = "; ")
    ))
  }

  if (nrow(data_multi) > 0L) {
    message(glue::glue(
      "Multi-mapped {label}, showing up to {n_show}: ",
      paste(utils::head(data_multi$name_original, n_show), collapse = "; ")
    ))
  }

  list(
    summary = data.frame(
      label = label,
      n_name = nrow(data_by_name),
      n_mapped = sum(data_by_name$kegg_n > 0L),
      n_unmapped = nrow(data_unmapped),
      n_multi = nrow(data_multi),
      stringsAsFactors = FALSE
    ),
    unmapped = data_unmapped,
    multi = data_multi
  )
}

metaFuns$get_kegg_pathways_for_compounds <- function(vec_kegg_id)
{
  vec_kegg_id <- unique(trimws(as.character(vec_kegg_id)))
  vec_kegg_id <- sub("^cpd:", "", vec_kegg_id)
  vec_kegg_id <- vec_kegg_id[!is.na(vec_kegg_id) & nzchar(vec_kegg_id)]

  if (length(vec_kegg_id) == 0L) {
    return(data.frame(
      kegg_id = character(0L),
      pathway_id = character(0L),
      pathway_name = character(0L),
      stringsAsFactors = FALSE
    ))
  }

  vec_query <- paste0("cpd:", vec_kegg_id)

  vec_link <- tryCatch(
    KEGGREST::keggLink("pathway", vec_query),
    error = function(e) {
      message(glue::glue("KEGG pathway link failed: {conditionMessage(e)}"))

      character(0L)
    }
  )

  if (length(vec_link) == 0L) {
    return(data.frame(
      kegg_id = character(0L),
      pathway_id = character(0L),
      pathway_name = character(0L),
      stringsAsFactors = FALSE
    ))
  }

  data_link <- data.frame(
    kegg_id = sub("^cpd:", "", names(vec_link)),
    pathway_id = sub("^path:", "", as.character(vec_link)),
    stringsAsFactors = FALSE
  )

  vec_path <- tryCatch(
    KEGGREST::keggList("pathway"),
    error = function(e) {
      message(glue::glue("KEGG pathway list failed: {conditionMessage(e)}"))
      return(NULL)
    }
  )

  if (is.null(vec_path)) {
    return(NULL)
  }

  if (length(vec_path) == 0L) {
    data_link$pathway_name <- NA_character_

    return(dplyr::distinct(data_link))
  }

  data_path <- data.frame(
    pathway_id = sub("^path:", "", names(vec_path)),
    pathway_name = as.character(vec_path),
    stringsAsFactors = FALSE
  )

  data_link <- dplyr::left_join(data_link, data_path, by = "pathway_id")
  data_link <- dplyr::distinct(data_link)

  data_link
}

metaFuns$get_mebocost_metabolite_names <- function(x, source = c("commu_res", "diff_commu", "overall_score"))
{
  source <- match.arg(source)

  if (source == "commu_res") {
    data_mebo <- tibble::as_tibble(x@tables$step3$t.commu_res)
  } else if (source == "diff_commu") {
    data_mebo <- tibble::as_tibble(x@tables$step4$ts.diff_commu[[1L]])
  } else {
    data_mebo <- tibble::as_tibble(x@tables$step5$t.overallScore)
  }

  if (!"Metabolite_Name" %in% colnames(data_mebo)) {
    stop('"Metabolite_Name" was not found in MEBOCOST table.')
  }

  unique(data_mebo$Metabolite_Name)
}

metaFuns$diagnose_kegg_pathway_overlap <- function(x, ref,
  mebo_source = c("commu_res", "diff_commu", "overall_score"),
  p_cutoff = 0.05, vip_cutoff = NULL, padj_cutoff = NULL,
  dir_cache = file.path(x$dir_cache, "kegg"),
  manual_kegg_mebo = NULL, manual_kegg_metabo = NULL)
{
  mebo_source <- match.arg(mebo_source)

  if (!dir.exists(dir_cache)) {
    dir.create(dir_cache, recursive = TRUE, showWarnings = FALSE)
  }

  vec_mebo_name <- metaFuns$get_mebocost_metabolite_names(x, source = mebo_source)

  data_metabo_sig <- metaFuns$filter_metabo_significant(
    data_metabo_diff = ref@tables$step3$data_diff_report,
    p_cutoff = p_cutoff,
    vip_cutoff = vip_cutoff,
    padj_cutoff = padj_cutoff
  )

  vec_metabo_name <- unique(data_metabo_sig$feature_name)

  data_mebo_map <- metaFuns$map_names_to_kegg_compound(
    vec_mebo_name,
    dir_cache = dir_cache,
    manual_kegg = manual_kegg_mebo,
    cache_name = "kegg_compound_name_table"
  )

  data_metabo_map <- metaFuns$map_names_to_kegg_compound(
    vec_metabo_name,
    dir_cache = dir_cache,
    manual_kegg = manual_kegg_metabo,
    cache_name = "kegg_compound_name_table"
  )

  diag_mebo <- metaFuns$diagnose_kegg_mapping(data_mebo_map, "MEBOCOST metabolites")
  diag_metabo <- metaFuns$diagnose_kegg_mapping(data_metabo_map, "significant metabolomics metabolites")

  vec_mebo_kegg_id <- unique(as.character(data_mebo_map$kegg_id))
  vec_mebo_kegg_id <- vec_mebo_kegg_id[!is.na(vec_mebo_kegg_id) & nzchar(vec_mebo_kegg_id)]

  vec_metabo_kegg_id <- unique(as.character(data_metabo_map$kegg_id))
  vec_metabo_kegg_id <- vec_metabo_kegg_id[!is.na(vec_metabo_kegg_id) & nzchar(vec_metabo_kegg_id)]

  data_mebo_path <- expect_local_data(
    dir_cache,
    glue::glue("kegg_pathway_links_mebocost_{mebo_source}"),
    metaFuns$get_kegg_pathways_for_compounds,
    list(vec_kegg_id = vec_mebo_kegg_id)
  )

  data_metabo_path <- expect_local_data(
    dir_cache,
    "kegg_pathway_links_metabolomics_sig",
    metaFuns$get_kegg_pathways_for_compounds,
    list(vec_kegg_id = vec_metabo_kegg_id)
  )

  data_mebo_path <- dplyr::left_join(
    data_mebo_map,
    data_mebo_path,
    by = "kegg_id"
  )

  data_metabo_path <- dplyr::left_join(
    data_metabo_map,
    data_metabo_path,
    by = "kegg_id"
  )

  data_metabo_path <- dplyr::left_join(
    data_metabo_path,
    data_metabo_sig,
    by = c("name_original" = "feature_name")
  )

  vec_shared_pathway <- intersect(
    unique(data_mebo_path$pathway_id[!is.na(data_mebo_path$pathway_id)]),
    unique(data_metabo_path$pathway_id[!is.na(data_metabo_path$pathway_id)])
  )

  data_pathway_summary <- lapply(vec_shared_pathway, function(pathway) {
    data_mebo_sub <- data_mebo_path[data_mebo_path$pathway_id == pathway, , drop = FALSE]
    data_metabo_sub <- data_metabo_path[data_metabo_path$pathway_id == pathway, , drop = FALSE]

    data.frame(
      pathway_id = pathway,
      pathway_name = data_mebo_sub$pathway_name[which(!is.na(data_mebo_sub$pathway_name))[1L]],
      n_mebo_metabolite = metaFuns$n_distinct_non_na(data_mebo_sub$name_original),
      mebo_metabolites = metaFuns$collapse_non_na(data_mebo_sub$name_original),
      n_metabo_sig = metaFuns$n_distinct_non_na(data_metabo_sub$name_original),
      metabo_features = metaFuns$collapse_non_na(data_metabo_sub$name_original),
      best_metabo_pvalue = min(data_metabo_sub$pvalue, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  data_pathway_summary <- dplyr::bind_rows(data_pathway_summary)

  if (nrow(data_pathway_summary) > 0L) {
    data_pathway_summary <- dplyr::arrange(data_pathway_summary, best_metabo_pvalue)
  }

  message(glue::glue(
    "KEGG pathway overlap diagnosis: ",
    "MEBOCOST source = {mebo_source}; ",
    "shared KEGG pathways = {nrow(data_pathway_summary)}."
  ))

  structure(list(
    data_pathway_summary = data_pathway_summary,
    data_mebo_map = data_mebo_map,
    data_metabo_map = data_metabo_map,
    data_mebo_path = data_mebo_path,
    data_metabo_path = data_metabo_path,
    data_metabo_sig = data_metabo_sig,
    diag_mebo = diag_mebo,
    diag_metabo = diag_metabo
  ), class = "kegg_results_metaboDiff_and_mebocost")
}

metaFuns$filter_kegg_pathway_summary <- function(data_pathway_summary)
{
  data_pathway_summary <- tibble::as_tibble(data_pathway_summary)

  vec_remove <- c(
    "Metabolic pathways",
    "Biosynthesis of secondary metabolites",
    "Microbial metabolism in diverse environments",
    "Biosynthesis of cofactors",
    "Carbon metabolism"
  )

  data_pathway_summary <- data_pathway_summary[
    !data_pathway_summary$pathway_name %in% vec_remove,
  ]

  data_pathway_summary <- dplyr::arrange(
    data_pathway_summary,
    dplyr::desc(n_metabo_sig),
    dplyr::desc(n_mebo_metabolite),
    best_metabo_pvalue
  )

  data_pathway_summary
}

# ==========================================================================
# FELLA

metaFuns$assert_kegg_results_metaboDiff_and_mebocost <- function(data_kegg_bridge)
{
  if (!inherits(data_kegg_bridge, "kegg_results_metaboDiff_and_mebocost")) {
    stop('!inherits(data_kegg_bridge, "kegg_results_metaboDiff_and_mebocost").')
  }

  vec_required <- c(
    "data_pathway_summary",
    "data_mebo_map",
    "data_metabo_map",
    "data_mebo_path",
    "data_metabo_path",
    "data_metabo_sig",
    "diag_mebo",
    "diag_metabo"
  )

  vec_missing <- setdiff(vec_required, names(data_kegg_bridge))

  if (length(vec_missing) > 0L) {
    stop(glue::glue(
      "Missing required field(s) in data_kegg_bridge: {paste(vec_missing, collapse = ', ')}."
    ))
  }

  invisible(TRUE)
}

metaFuns$get_or_build_fella_data <- function(dir_fella, organism = "hsa",
  method = c("diffusion", "pagerank", "hypergeom"), rebuild = FALSE,
  repair = FALSE, filter.path = c("01100", "01200", "01210", "01212", "01230"),
  niter = 100L)
{
  expect_package("FELLA", "1.22.0")

  method <- match.arg(method)

  load_matrix <- if (method %in% c("diffusion", "pagerank")) {
    method
  } else {
    NULL
  }

  if (!isTRUE(rebuild)) {
    data_fella <- tryCatch(
      FELLA::loadKEGGdata(
        databaseDir = dir_fella,
        internalDir = FALSE,
        loadMatrix = load_matrix
      ),
      error = function(e) {
        message(glue::glue(
          "FELLA database loading failed: {conditionMessage(e)}"
        ))

        NULL
      }
    )

    if (!is.null(data_fella)) {
      message(glue::glue("Loaded existing FELLA database: {dir_fella}"))

      return(data_fella)
    }
  }

  if (dir.exists(dir_fella)) {
    if (isTRUE(rebuild) || isTRUE(repair)) {
      message(glue::glue("Remove existing FELLA database directory: {dir_fella}"))
      unlink(dir_fella, recursive = TRUE, force = TRUE)
    } else {
      stop(glue::glue(
        "FELLA database directory already exists but could not be loaded: {dir_fella}. ",
        "Use `rebuild = TRUE` or `repair = TRUE` to remove and rebuild it."
      ))
    }
  }

  parent_dir <- dirname(dir_fella)

  if (!dir.exists(parent_dir)) {
    dir.create(parent_dir, recursive = TRUE, showWarnings = FALSE)
  }

  message(glue::glue(
    "Build FELLA KEGG knowledge model in: {dir_fella}. ",
    "This may take time and requires network access."
  ))

  graph_kegg <- FELLA::buildGraphFromKEGGREST(
    organism = organism,
    filter.path = filter.path
  )

  matrices <- if (method %in% c("diffusion", "pagerank")) {
    method
  } else {
    "hypergeom"
  }

  normality <- if (method %in% c("diffusion", "pagerank")) {
    method
  } else {
    character(0L)
  }

  FELLA::buildDataFromGraph(
    keggdata.graph = graph_kegg,
    databaseDir = dir_fella,
    internalDir = FALSE,
    matrices = matrices,
    normality = normality,
    niter = as.integer(niter)
  )

  FELLA::loadKEGGdata(
    databaseDir = dir_fella,
    internalDir = FALSE,
    loadMatrix = load_matrix
  )
}

metaFuns$get_fella_compounds_by_source <- function(data_kegg_bridge,
  source = c("mebocost", "metabolomics"), shared_pathway_only = FALSE)
{
  metaFuns$assert_kegg_results_metaboDiff_and_mebocost(data_kegg_bridge)

  source <- match.arg(source)

  if (isTRUE(shared_pathway_only)) {
    data_pathway <- tibble::as_tibble(data_kegg_bridge$data_pathway_summary)
    vec_pathway <- unique(data_pathway$pathway_id)
    vec_pathway <- vec_pathway[!is.na(vec_pathway) & nzchar(vec_pathway)]

    if (source == "mebocost") {
      data_path <- tibble::as_tibble(data_kegg_bridge$data_mebo_path)
    } else {
      data_path <- tibble::as_tibble(data_kegg_bridge$data_metabo_path)
    }

    data_path <- data_path[data_path$pathway_id %in% vec_pathway, ]
    vec_cpd <- data_path$kegg_id
  } else {
    if (source == "mebocost") {
      data_map <- tibble::as_tibble(data_kegg_bridge$data_mebo_map)
    } else {
      data_map <- tibble::as_tibble(data_kegg_bridge$data_metabo_map)
    }

    vec_cpd <- data_map$kegg_id
  }

  vec_cpd <- unique(trimws(as.character(vec_cpd)))
  vec_cpd <- sub("^cpd:", "", vec_cpd)
  vec_cpd <- vec_cpd[!is.na(vec_cpd) & nzchar(vec_cpd)]

  vec_cpd
}

metaFuns$run_fella_enrichment_one <- function(vec_cpd, data_fella,
  label, method = "hypergeom", approx = "normality",
  threshold = 1, compounds_background = NULL, niter = 100L)
{
  vec_cpd <- unique(as.character(vec_cpd))
  vec_cpd <- sub("^cpd:", "", vec_cpd)
  vec_cpd <- vec_cpd[!is.na(vec_cpd) & nzchar(vec_cpd)]

  vec_valid_cpd <- tryCatch(
    FELLA::getCom(data_fella, level = "compound", format = "name"),
    error = function(e) character(0L)
  )

  if (length(vec_valid_cpd) > 0L) {
    vec_excluded <- setdiff(vec_cpd, vec_valid_cpd)
    vec_cpd <- intersect(vec_cpd, vec_valid_cpd)
  } else {
    vec_excluded <- character(0L)
  }

  if (length(vec_cpd) == 0L) {
    stop(glue::glue("No valid FELLA compound was found for {label}."))
  }

  message(glue::glue(
    "Run FELLA for {label}: method = {method}; input compounds = {length(vec_cpd)}; ",
    "excluded compounds = {length(vec_excluded)}."
  ))

  if (is.null(compounds_background)) {
    obj_fella <- FELLA::enrich(
      compounds = vec_cpd,
      method = method,
      approx = approx,
      niter = as.integer(niter),
      data = data_fella
    )
  } else {
    compounds_background <- unique(as.character(compounds_background))
    compounds_background <- sub("^cpd:", "", compounds_background)
    compounds_background <- compounds_background[
      !is.na(compounds_background) & nzchar(compounds_background)
    ]

    obj_fella <- FELLA::enrich(
      compounds = vec_cpd,
      compoundsBackground = compounds_background,
      method = method,
      approx = approx,
      niter = as.integer(niter),
      data = data_fella
    )
  }

  data_table <- tryCatch(
    FELLA::generateResultsTable(
      object = obj_fella,
      data = data_fella,
      method = method,
      threshold = threshold
    ),
    error = function(e) {
      message(glue::glue("No FELLA result table for {label}: {conditionMessage(e)}"))
      data.frame()
    }
  )

  data_table <- metaFuns$standardize_fella_table(data_table)

  list(
    obj_fella = obj_fella,
    data_table = data_table,
    compounds = vec_cpd,
    excluded_compounds = vec_excluded,
    label = label
  )
}

metaFuns$run_fella_dual_bridge_from_kegg_bridge <- function(data_kegg_bridge,
  dir_fella = .prefix(paste0("fella_", organism), "db"),
  organism = "hsa", method = c("hypergeom", "diffusion", "pagerank"),
  approx = "normality", threshold = 1, bridge_p_cutoff = 0.1,
  shared_pathway_only = FALSE, rebuild = FALSE, repair = FALSE,
  niter = 100L)
{
  metaFuns$assert_kegg_results_metaboDiff_and_mebocost(data_kegg_bridge)

  method <- match.arg(method)

  data_fella <- metaFuns$get_or_build_fella_data(
    dir_fella = dir_fella,
    organism = organism,
    method = method,
    rebuild = rebuild,
    repair = repair,
    niter = niter
  )

  vec_mebo_cpd <- metaFuns$get_fella_compounds_by_source(
    data_kegg_bridge = data_kegg_bridge,
    source = "mebocost",
    shared_pathway_only = shared_pathway_only
  )

  vec_metabo_cpd <- metaFuns$get_fella_compounds_by_source(
    data_kegg_bridge = data_kegg_bridge,
    source = "metabolomics",
    shared_pathway_only = shared_pathway_only
  )

  res_mebo <- metaFuns$run_fella_enrichment_one(
    vec_cpd = vec_mebo_cpd,
    data_fella = data_fella,
    label = "MEBOCOST",
    method = method,
    approx = approx,
    threshold = threshold,
    niter = niter
  )

  res_metabo <- metaFuns$run_fella_enrichment_one(
    vec_cpd = vec_metabo_cpd,
    data_fella = data_fella,
    label = "metabolomics",
    method = method,
    approx = approx,
    threshold = threshold,
    niter = niter
  )

  data_mebo <- metaFuns$standardize_fella_table(res_mebo$data_table)
  data_metabo <- metaFuns$standardize_fella_table(res_metabo$data_table)

  data_bridge <- dplyr::inner_join(
    data_mebo,
    data_metabo,
    by = "KEGG.id",
    suffix = c("_mebocost", "_metabolomics")
  )

  if (nrow(data_bridge) > 0L) {
    data_bridge$bridge_stat <- -2 * (
      log(pmax(data_bridge$p.value_mebocost, .Machine$double.xmin)) +
        log(pmax(data_bridge$p.value_metabolomics, .Machine$double.xmin))
    )

    data_bridge$bridge_p_combined <- stats::pchisq(
      data_bridge$bridge_stat,
      df = 4L,
      lower.tail = FALSE
    )

    data_bridge$bridge_score <- -log10(
      pmax(data_bridge$bridge_p_combined, .Machine$double.xmin)
    )

    data_bridge$bridge_pass <- data_bridge$p.value_mebocost < bridge_p_cutoff &
      data_bridge$p.value_metabolomics < bridge_p_cutoff

    data_bridge <- dplyr::arrange(
      data_bridge,
      dplyr::desc(bridge_pass),
      bridge_p_combined,
      p.value_mebocost,
      p.value_metabolomics
    )
  }
  n_bridge_pass <- 0L

  if (nrow(data_bridge) > 0L && "bridge_pass" %in% colnames(data_bridge)) {
    n_bridge_pass <- sum(data_bridge$bridge_pass, na.rm = TRUE)
  }

  message(glue::glue(
      "FELLA dual bridge finished: MEBOCOST results = {nrow(data_mebo)}, ",
      "metabolomics results = {nrow(data_metabo)}, shared KEGG nodes = {nrow(data_bridge)}, ",
      "strict bridge-pass nodes = {n_bridge_pass}."
      ))

  out <- list(
    data_bridge = data_bridge,
    data_mebocost_fella = data_mebo,
    data_metabolomics_fella = data_metabo,
    res_mebocost = res_mebo,
    res_metabolomics = res_metabo,
    data_fella = data_fella,
    method = method,
    threshold = threshold,
    bridge_p_cutoff = bridge_p_cutoff,
    shared_pathway_only = shared_pathway_only,
    data_kegg_bridge = data_kegg_bridge
  )

  class(out) <- c("fella_dual_bridge_metaboDiff_and_mebocost", class(out))

  out
}

metaFuns$standardize_fella_table <- function(data_table)
{
  data_table <- tibble::as_tibble(data_table)

  vec_need <- c("KEGG.id", "KEGG.name", "p.value")
  vec_missing <- setdiff(vec_need, colnames(data_table))

  if (length(vec_missing) > 0L) {
    data_table <- tibble::tibble(
      KEGG.id = character(0L),
      KEGG.name = character(0L),
      p.value = numeric(0L)
    )

    return(data_table)
  }

  if (!"padj" %in% colnames(data_table)) {
    data_table$padj <- stats::p.adjust(data_table$p.value, method = "BH")
  }

  data_table
}

metaFuns$get_fella_compounds_from_kegg_bridge <- function(data_kegg_bridge,
  compound_source = c("shared_pathway", "both", "mebocost", "metabolomics"))
{
  metaFuns$assert_kegg_results_metaboDiff_and_mebocost(data_kegg_bridge)

  compound_source <- match.arg(compound_source)

  data_mebo_map <- tibble::as_tibble(data_kegg_bridge$data_mebo_map)
  data_metabo_map <- tibble::as_tibble(data_kegg_bridge$data_metabo_map)

  if (compound_source == "mebocost") {
    vec_cpd <- data_mebo_map$kegg_id
  } else if (compound_source == "metabolomics") {
    vec_cpd <- data_metabo_map$kegg_id
  } else if (compound_source == "both") {
    vec_cpd <- c(data_mebo_map$kegg_id, data_metabo_map$kegg_id)
  } else {
    data_pathway_summary <- tibble::as_tibble(data_kegg_bridge$data_pathway_summary)
    data_mebo_path <- tibble::as_tibble(data_kegg_bridge$data_mebo_path)
    data_metabo_path <- tibble::as_tibble(data_kegg_bridge$data_metabo_path)

    vec_pathway <- unique(data_pathway_summary$pathway_id)
    vec_pathway <- vec_pathway[!is.na(vec_pathway) & nzchar(vec_pathway)]

    data_mebo_path <- data_mebo_path[
      data_mebo_path$pathway_id %in% vec_pathway,
    ]

    data_metabo_path <- data_metabo_path[
      data_metabo_path$pathway_id %in% vec_pathway,
    ]

    vec_cpd <- c(data_mebo_path$kegg_id, data_metabo_path$kegg_id)
  }

  vec_cpd <- unique(trimws(as.character(vec_cpd)))
  vec_cpd <- sub("^cpd:", "", vec_cpd)
  vec_cpd <- vec_cpd[!is.na(vec_cpd) & nzchar(vec_cpd)]

  vec_cpd
}

metaFuns$run_fella_from_kegg_bridge <- function(data_kegg_bridge,
  dir_fella = .prefix(paste0("fella_", organism), "db"), organism = "hsa",
  compound_source = c("shared_pathway", "both", "mebocost", "metabolomics"),
  method = c("diffusion", "pagerank", "hypergeom"),
  approx = "normality", threshold = 0.05, rebuild = FALSE,
  repair = FALSE, compounds_background = NULL, niter = 100L)
{
  metaFuns$assert_kegg_results_metaboDiff_and_mebocost(data_kegg_bridge)

  compound_source <- match.arg(compound_source)
  method <- match.arg(method)

  data_fella <- metaFuns$get_or_build_fella_data(
    dir_fella = dir_fella,
    organism = organism,
    method = method,
    rebuild = rebuild,
    repair = repair,
    niter = niter
  )

  vec_cpd <- metaFuns$get_fella_compounds_from_kegg_bridge(
    data_kegg_bridge = data_kegg_bridge,
    compound_source = compound_source
  )

  if (length(vec_cpd) == 0L) {
    stop("No KEGG compound IDs were available for FELLA.")
  }

  vec_valid_cpd <- tryCatch(
    FELLA::getCom(data_fella, level = "compound", format = "name"),
    error = function(e) character(0L)
  )

  if (length(vec_valid_cpd) > 0L) {
    vec_excluded <- setdiff(vec_cpd, vec_valid_cpd)
    vec_cpd <- intersect(vec_cpd, vec_valid_cpd)
  } else {
    vec_excluded <- character(0L)
  }

  if (length(vec_cpd) == 0L) {
    stop("No input KEGG compound IDs were found in the FELLA graph.")
  }

  message(glue::glue(
    "Run FELLA: method = {method}; input compounds = {length(vec_cpd)}; ",
    "excluded compounds = {length(vec_excluded)}."
  ))

  if (is.null(compounds_background)) {
    obj_fella <- FELLA::enrich(
      compounds = vec_cpd,
      method = method,
      approx = approx,
      niter = as.integer(niter),
      data = data_fella
    )
  } else {
    compounds_background <- unique(as.character(compounds_background))
    compounds_background <- sub("^cpd:", "", compounds_background)
    compounds_background <- compounds_background[
      !is.na(compounds_background) & nzchar(compounds_background)
    ]

    obj_fella <- FELLA::enrich(
      compounds = vec_cpd,
      compoundsBackground = compounds_background,
      method = method,
      approx = approx,
      niter = as.integer(niter),
      data = data_fella
    )
  }

  data_table <- FELLA::generateResultsTable(
    object = obj_fella,
    data = data_fella,
    method = method,
    threshold = threshold
  )

  graph_fella <- FELLA::generateResultsGraph(
    object = obj_fella,
    data = data_fella,
    method = method,
    threshold = threshold
  )

  out <- list(
    obj_fella = obj_fella,
    data_fella = data_fella,
    data_table = tibble::as_tibble(data_table),
    graph_fella = graph_fella,
    compounds = vec_cpd,
    excluded_compounds = vec_excluded,
    compound_source = compound_source,
    method = method,
    threshold = threshold,
    data_kegg_bridge = data_kegg_bridge
  )

  class(out) <- c("fella_results_metaboDiff_and_mebocost", class(out))

  out
}

metaFuns$assert_fella_results_metaboDiff_and_mebocost <- function(data_fella_result)
{
  if (!inherits(data_fella_result, "fella_results_metaboDiff_and_mebocost")) {
    stop('!inherits(data_fella_result, "fella_results_metaboDiff_and_mebocost").')
  }

  if (is.null(data_fella_result$data_table)) {
    stop("data_fella_result$data_table is NULL.")
  }

  invisible(TRUE)
}

metaFuns$get_fella_supported_kegg_bridge <- function(data_kegg_bridge,
  data_fella_result, fella_p_cutoff = 0.1, require_fella_sig = TRUE)
{
  metaFuns$assert_kegg_results_metaboDiff_and_mebocost(data_kegg_bridge)
  metaFuns$assert_fella_results_metaboDiff_and_mebocost(data_fella_result)

  data_pathway <- tibble::as_tibble(data_kegg_bridge$data_pathway_summary)
  data_fella <- tibble::as_tibble(data_fella_result$data_table)

  if (!all(c("KEGG.id", "KEGG.name", "p.value") %in% colnames(data_fella))) {
    stop('data_fella_result$data_table must contain "KEGG.id", "KEGG.name", and "p.value".')
  }

  data_fella$pathway_id <- sub("^hsa", "map", data_fella$KEGG.id)
  data_fella$fella_pathway_name <- sub(" - Homo sapiens \\(human\\)$", "", data_fella$KEGG.name)
  data_fella$fella_pvalue <- data_fella$p.value
  data_fella$fella_padj <- stats::p.adjust(data_fella$fella_pvalue, method = "BH")

  data_fella <- dplyr::select(
    data_fella,
    pathway_id,
    fella_pathway_name,
    fella_pvalue,
    fella_padj,
    dplyr::everything()
  )

  data_bridge <- dplyr::inner_join(
    data_pathway,
    data_fella,
    by = "pathway_id"
  )

  if (isTRUE(require_fella_sig)) {
    data_bridge <- data_bridge[
      !is.na(data_bridge$fella_pvalue) &
        data_bridge$fella_pvalue < fella_p_cutoff,
    ]
  }

  if (nrow(data_bridge) > 0L) {
    data_bridge$bridge_score <- -log10(
      pmax(data_bridge$fella_pvalue, .Machine$double.xmin)
    ) * log1p(data_bridge$n_mebo_metabolite) *
      log1p(data_bridge$n_metabo_sig)

    data_bridge <- dplyr::arrange(
      data_bridge,
      dplyr::desc(bridge_score),
      fella_pvalue,
      best_metabo_pvalue
    )
  }

  message(glue::glue(
    "FELLA-supported KEGG bridge: {nrow(data_bridge)} pathway(s) retained."
  ))

  data_bridge
}

metaFuns$score_fella_supported_kegg_bridge <- function(data_bridge,
  weights = c(
    fella = 0.35,
    metabo = 0.25,
    coverage = 0.20,
    specificity = 0.10,
    balance = 0.10
  ),
  exclude_pathway_id = character(0L),
  exclude_pathway_name = character(0L))
{
  data_bridge <- tibble::as_tibble(data_bridge)

  vec_required <- c(
    "pathway_id", "pathway_name", "mebo_metabolites", "metabo_features",
    "best_metabo_pvalue", "fella_pvalue", "CompoundHits",
    "CompoundsInPathway"
  )

  vec_missing <- setdiff(vec_required, colnames(data_bridge))

  if (length(vec_missing) > 0L) {
    stop(glue::glue(
      "Missing required column(s): {paste(vec_missing, collapse = ', ')}."
    ))
  }

  weights <- weights / sum(weights)

  lst_mebo <- strsplit(as.character(data_bridge$mebo_metabolites), ";\\s*")
  lst_metabo <- strsplit(as.character(data_bridge$metabo_features), ";\\s*")

  data_bridge$n_mebo_metabolite_clean <- vapply(
    lst_mebo, metaFuns$n_distinct_non_na, integer(1L)
  )

  data_bridge$n_metabo_sig_clean <- vapply(
    lst_metabo, metaFuns$n_distinct_non_na, integer(1L)
  )

  data_bridge$mebo_metabolites_clean <- vapply(
    lst_mebo, metaFuns$collapse_non_na, character(1L)
  )

  data_bridge$metabo_features_clean <- vapply(
    lst_metabo, metaFuns$collapse_non_na, character(1L)
  )

  data_bridge$score_fella <- -log10(
    pmax(data_bridge$fella_pvalue, .Machine$double.xmin)
  )

  data_bridge$score_metabo <- -log10(
    pmax(data_bridge$best_metabo_pvalue, .Machine$double.xmin)
  )

  data_bridge$score_mebo_coverage <- log1p(data_bridge$n_mebo_metabolite_clean)
  data_bridge$score_metabo_coverage <- log1p(data_bridge$n_metabo_sig_clean)

  data_bridge$score_bridge_balance <- pmin(
    data_bridge$score_mebo_coverage,
    data_bridge$score_metabo_coverage
  ) / pmax(
    data_bridge$score_mebo_coverage,
    data_bridge$score_metabo_coverage
  )

  data_bridge$score_pathway_specificity <- data_bridge$CompoundHits /
    pmax(data_bridge$CompoundsInPathway, 1L)

  data_bridge$score_fella_scaled <- metaFuns$scale01(data_bridge$score_fella)
  data_bridge$score_metabo_scaled <- metaFuns$scale01(data_bridge$score_metabo)

  data_bridge$score_coverage_scaled <- metaFuns$scale01(
    sqrt(data_bridge$score_mebo_coverage * data_bridge$score_metabo_coverage)
  )

  data_bridge$score_specificity_scaled <- metaFuns$scale01(
    data_bridge$score_pathway_specificity
  )

  data_bridge$bridge_integrated_score <- 
    weights[["fella"]] * data_bridge$score_fella_scaled +
    weights[["metabo"]] * data_bridge$score_metabo_scaled +
    weights[["coverage"]] * data_bridge$score_coverage_scaled +
    weights[["specificity"]] * data_bridge$score_specificity_scaled +
    weights[["balance"]] * data_bridge$score_bridge_balance

  data_bridge$is_excluded_by_user <- data_bridge$pathway_id %in% exclude_pathway_id |
    data_bridge$pathway_name %in% exclude_pathway_name

  if (length(exclude_pathway_id) > 0L || length(exclude_pathway_name) > 0L) {
    data_bridge <- data_bridge[!data_bridge$is_excluded_by_user, ]
  }

  data_bridge <- dplyr::arrange(
    data_bridge,
    dplyr::desc(bridge_integrated_score),
    fella_pvalue,
    best_metabo_pvalue
  )

  data_bridge
}


metaFuns$split_semicolon_names <- function(x)
{
  x <- unlist(strsplit(as.character(x), ";\\s*"), use.names = FALSE)
  x <- trimws(x)
  x <- x[!is.na(x) & nzchar(x) & x != "NA"]

  unique(x)
}

metaFuns$get_bridge_metabolite_receiver_candidates <- function(data_bridge,
  data_overall_score, bridge_score_col = "bridge_integrated_score",
  n_top_per_pathway = Inf)
{
  data_bridge <- tibble::as_tibble(data_bridge)
  data_overall_score <- tibble::as_tibble(data_overall_score)

  if (!all(c("pathway_id", "pathway_name", "mebo_metabolites") %in% colnames(data_bridge))) {
    stop('data_bridge must contain "pathway_id", "pathway_name", and "mebo_metabolites".')
  }

  if (!all(c("Metabolite_Name", "Receiver", "overall_score") %in% colnames(data_overall_score))) {
    stop('data_overall_score must contain "Metabolite_Name", "Receiver", and "overall_score".')
  }

  if (!bridge_score_col %in% colnames(data_bridge)) {
    data_bridge[[bridge_score_col]] <- 1
  }

  lst_candidate <- lapply(seq_len(nrow(data_bridge)), function(i) {
    vec_metabolite <- metaFuns$split_semicolon_names(data_bridge$mebo_metabolites[i])

    if (length(vec_metabolite) == 0L) {
      return(data.frame())
    }

    data_sub <- data_overall_score[
      data_overall_score$Metabolite_Name %in% vec_metabolite,
    ]

    if (nrow(data_sub) == 0L) {
      return(data.frame())
    }

    data_sub$pathway_id <- data_bridge$pathway_id[i]
    data_sub$pathway_name <- data_bridge$pathway_name[i]
    data_sub$fella_pvalue <- data_bridge$fella_pvalue[i]
    data_sub$best_metabo_pvalue <- data_bridge$best_metabo_pvalue[i]
    data_sub$bridge_score <- data_bridge[[bridge_score_col]][i]
    data_sub$metabo_features <- data_bridge$metabo_features[i]

    data_sub
  })

  data_candidate <- dplyr::bind_rows(lst_candidate)

  if (nrow(data_candidate) == 0L) {
    return(data_candidate)
  }

  data_candidate$score_overall_scaled <- metaFuns$scale01(data_candidate$overall_score)
  data_candidate$score_bridge_scaled <- metaFuns$scale01(data_candidate$bridge_score)

  data_candidate$candidate_score <- 0.65 * data_candidate$score_overall_scaled +
    0.35 * data_candidate$score_bridge_scaled

  data_candidate <- dplyr::group_by(data_candidate, pathway_id)

  data_candidate <- dplyr::arrange(
    data_candidate,
    dplyr::desc(candidate_score),
    dplyr::desc(overall_score),
    fella_pvalue,
    best_metabo_pvalue,
    .by_group = TRUE
  )

  data_candidate <- dplyr::mutate(
    data_candidate,
    rank_in_pathway = dplyr::row_number()
  )

  if (!is.null(n_top_per_pathway) && is.finite(n_top_per_pathway)) {
    data_candidate <- dplyr::filter(
      data_candidate,
      rank_in_pathway <= as.integer(n_top_per_pathway)
    )
  }

  data_candidate <- dplyr::ungroup(data_candidate)

  data_candidate <- dplyr::arrange(
    data_candidate,
    dplyr::desc(candidate_score),
    dplyr::desc(overall_score),
    fella_pvalue,
    best_metabo_pvalue,
    rank_in_pathway
  )

  data_candidate
}

metaFuns$summarise_bridge_metabolite_receiver <- function(data_candidate)
{
  data_candidate <- tibble::as_tibble(data_candidate)

  vec_required <- c(
    "Metabolite_Name", "Receiver", "overall_score", "pathway_id",
    "pathway_name", "fella_pvalue", "best_metabo_pvalue",
    "bridge_score", "metabo_features"
  )

  vec_missing <- setdiff(vec_required, colnames(data_candidate))

  if (length(vec_missing) > 0L) {
    stop(glue::glue(
      "Missing required column(s): {paste(vec_missing, collapse = ', ')}."
    ))
  }

  data_pair <- dplyr::group_by(data_candidate, Metabolite_Name, Receiver)

  data_pair <- dplyr::summarise(
    data_pair,
    overall_score = max(overall_score, na.rm = TRUE),
    n_bridge_pathway = dplyr::n_distinct(pathway_id),
    pathway_ids = metaFuns$collapse_non_na(pathway_id),
    pathway_names = metaFuns$collapse_non_na(pathway_name),
    metabo_features = metaFuns$collapse_non_na(metabo_features),
    best_fella_pvalue = min(fella_pvalue, na.rm = TRUE),
    best_metabo_pvalue = min(best_metabo_pvalue, na.rm = TRUE),
    max_bridge_score = max(bridge_score, na.rm = TRUE),
    .groups = "drop"
  )

  data_pair$score_overall_scaled <- metaFuns$scale01(data_pair$overall_score)
  data_pair$score_bridge_scaled <- metaFuns$scale01(data_pair$max_bridge_score)
  data_pair$score_pathway_count_scaled <- metaFuns$scale01(log1p(data_pair$n_bridge_pathway))

  data_pair$key_axis_score <- 0.70 * data_pair$score_overall_scaled +
    0.20 * data_pair$score_bridge_scaled +
    0.10 * data_pair$score_pathway_count_scaled

  data_pair <- dplyr::arrange(
    data_pair,
    dplyr::desc(key_axis_score),
    dplyr::desc(overall_score),
    best_fella_pvalue,
    best_metabo_pvalue
  )

  data_pair
}

metaFuns$rank_bridge_candidate_by_pathway <- function(data_candidate,
  n_top_per_pathway = 3L)
{
  data_candidate <- tibble::as_tibble(data_candidate)

  vec_required <- c(
    "pathway_id", "pathway_name", "Metabolite_Name", "Receiver",
    "overall_score", "fella_pvalue", "best_metabo_pvalue",
    "bridge_score", "metabo_features"
  )

  vec_missing <- setdiff(vec_required, colnames(data_candidate))

  if (length(vec_missing) > 0L) {
    stop(glue::glue(
      "Missing required column(s): {paste(vec_missing, collapse = ', ')}."
    ))
  }

  data_candidate$score_overall_scaled <- metaFuns$scale01(data_candidate$overall_score)
  data_candidate$score_bridge_scaled <- metaFuns$scale01(data_candidate$bridge_score)

  data_candidate$candidate_score <- 0.75 * data_candidate$score_overall_scaled +
    0.25 * data_candidate$score_bridge_scaled

  data_candidate <- dplyr::group_by(data_candidate, pathway_id, pathway_name)
  data_candidate <- dplyr::arrange(
    data_candidate,
    dplyr::desc(candidate_score),
    dplyr::desc(overall_score),
    fella_pvalue,
    best_metabo_pvalue,
    .by_group = TRUE
  )
  data_candidate <- dplyr::mutate(
    data_candidate,
    rank_in_pathway = dplyr::row_number()
  )

  if (!is.null(n_top_per_pathway)) {
    data_candidate <- dplyr::filter(
      data_candidate,
      rank_in_pathway <= as.integer(n_top_per_pathway)
    )
  }

  data_candidate <- dplyr::ungroup(data_candidate)

  dplyr::arrange(
    data_candidate,
    fella_pvalue,
    best_metabo_pvalue,
    rank_in_pathway
  )
}

metaFuns$collapse_non_na <- function(x)
{
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x) & x != "NA"]

  if (length(x) == 0L) {
    return(NA_character_)
  }

  paste(x, collapse = "; ")
}

metaFuns$n_distinct_non_na <- function(x)
{
  x <- unique(as.character(x))
  x <- x[!is.na(x) & nzchar(x) & x != "NA"]

  length(x)
}

# ==========================================================================

# metaFuns$run_pathview_safe <- function(cpd_data, pathway_id, out_dir,
#   species = "hsa", out_suffix = NULL, kegg_native = TRUE,
#   clean_old = FALSE, ...)
# {
#   expect_package("pathview", "1.40.0")
#
#   pathway_id <- sub("^hsa", "", as.character(pathway_id))
#   pathway_id <- sub("^map", "", pathway_id)
#
#   if (is.null(out_suffix)) {
#     out_suffix <- glue::glue("pathview_{pathway_id}")
#   }
#
#   out_dir <- normalizePath(out_dir, mustWork = FALSE)
#
#   if (!dir.exists(out_dir)) {
#     dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#   }
#
#   if (isTRUE(clean_old)) {
#     vec_old <- list.files(
#       out_dir,
#       pattern = glue::glue("{species}{pathway_id}.*{out_suffix}"),
#       full.names = TRUE
#     )
#
#     if (length(vec_old) > 0L) {
#       unlink(vec_old, force = TRUE)
#     }
#   }
#
#   vec_before <- list.files(out_dir, full.names = TRUE, all.files = FALSE)
#
#   old_wd <- getwd()
#   on.exit(setwd(old_wd), add = TRUE)
#
#   setwd(out_dir)
#
#   out <- pathview::pathview(
#     cpd.data = cpd_data,
#     pathway.id = pathway_id,
#     species = species,
#     out.suffix = out_suffix,
#     kegg.dir = out_dir,
#     kegg.native = kegg_native,
#     ...
#   )
#
#   vec_after <- list.files(out_dir, full.names = TRUE, all.files = FALSE)
#   vec_created <- setdiff(vec_after, vec_before)
#
#   list(
#     pathview_result = out,
#     pathway_id = pathway_id,
#     species = species,
#     out_suffix = out_suffix,
#     out_dir = out_dir,
#     created_files = vec_created
#   )
# }

# ==========================================================================
# plot

metaFuns$plot_kegg_bridge_coverage_barplot <- function(data_pathway_summary,
  n_top = 20L, min_mebo = 1L, min_metabo = 1L,
  rank_by = c("total_count", "balanced_count"),
  exclude_pathway_id = character(0L),
  exclude_pathway_name = character(0L), label_width = 45L)
{
  rank_by <- match.arg(rank_by)
  data_plot <- tibble::as_tibble(data_pathway_summary)

  vec_required <- c(
    "pathway_id", "pathway_name", "n_mebo_metabolite",
    "n_metabo_sig"
  )

  vec_missing <- setdiff(vec_required, colnames(data_plot))

  if (length(vec_missing) > 0L) {
    stop(glue::glue(
      "Missing required column(s): {paste(vec_missing, collapse = ', ')}."
    ))
  }

  .count_non_na_names <- function(x)
  {
    x <- unlist(strsplit(as.character(x), ";\\s*"), use.names = FALSE)
    x <- unique(trimws(x))
    x <- x[!is.na(x) & nzchar(x) & x != "NA"]

    length(x)
  }

  .wrap_text <- function(x, width = 45L)
  {
    vapply(as.character(x), function(s) {
      paste(strwrap(s, width = as.integer(width)), collapse = "\n")
    }, character(1L))
  }

  if ("mebo_metabolites" %in% colnames(data_plot)) {
    data_plot$n_mebo_plot <- vapply(
      data_plot$mebo_metabolites,
      .count_non_na_names,
      integer(1L)
    )
  } else {
    data_plot$n_mebo_plot <- data_plot$n_mebo_metabolite
  }

  if ("metabo_features" %in% colnames(data_plot)) {
    data_plot$n_metabo_plot <- vapply(
      data_plot$metabo_features,
      .count_non_na_names,
      integer(1L)
    )
  } else {
    data_plot$n_metabo_plot <- data_plot$n_metabo_sig
  }

  data_plot <- data_plot[
    !is.na(data_plot$pathway_id) &
      !is.na(data_plot$pathway_name) &
      data_plot$n_mebo_plot >= min_mebo &
      data_plot$n_metabo_plot >= min_metabo,
  ]

  if (length(exclude_pathway_id) > 0L) {
    data_plot <- data_plot[!data_plot$pathway_id %in% exclude_pathway_id, ]
  }

  if (length(exclude_pathway_name) > 0L) {
    data_plot <- data_plot[!data_plot$pathway_name %in% exclude_pathway_name, ]
  }

  if (nrow(data_plot) == 0L) {
    stop("No KEGG bridge pathway remained for plotting.")
  }

  data_plot$total_count <- data_plot$n_mebo_plot + data_plot$n_metabo_plot
  data_plot$balanced_count <- sqrt(data_plot$n_mebo_plot * data_plot$n_metabo_plot)

  if (rank_by == "balanced_count") {
    data_plot <- dplyr::arrange(
      data_plot,
      dplyr::desc(balanced_count),
      dplyr::desc(total_count)
    )
  } else {
    data_plot <- dplyr::arrange(
      data_plot,
      dplyr::desc(total_count),
      dplyr::desc(balanced_count)
    )
  }

  data_plot <- utils::head(data_plot, as.integer(n_top))

  data_plot$pathway_label <- glue::glue(
    "{data_plot$pathway_id}  {data_plot$pathway_name}"
  )
  data_plot$pathway_label <- .wrap_text(data_plot$pathway_label, label_width)

  data_long <- data.frame(
    pathway_label = rep(data_plot$pathway_label, each = 2L),
    evidence_source = rep(c("MEBOCOST", "Metabolomics"), times = nrow(data_plot)),
    n_mapped_metabolite = as.integer(c(rbind(
      data_plot$n_mebo_plot,
      data_plot$n_metabo_plot
    ))),
    stringsAsFactors = FALSE
  )

  data_long$pathway_label <- factor(
    data_long$pathway_label,
    levels = rev(data_plot$pathway_label)
  )

  data_long$evidence_source <- factor(
    data_long$evidence_source,
    levels = c("MEBOCOST", "Metabolomics")
  )

  p <- ggplot(
    data_long,
    aes(x = n_mapped_metabolite, y = pathway_label)
  ) +
    geom_col(aes(fill = evidence_source), width = 0.72) +
    geom_text(
      aes(label = n_mapped_metabolite),
      hjust = -0.25,
      size = 3
    ) +
    facet_wrap(~ evidence_source, nrow = 1L) +
    scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    labs(
      x = "Mapped metabolites in pathway",
      y = "Shared KEGG pathway",
      fill = "Evidence source"
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 9),
      axis.title = element_text(size = 10),
      strip.text = element_text(size = 10),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )

  p
}

metaFuns$plot_fella_bridge_sankey <- function(data_fella_bridge, data_kegg_bridge,
  n_pathway = 10L,
  n_metabolite = 60L,
  p_cutoff = 0.05,
  metabolite_min_pathway = 1L,
  point_pals = c("grey92", "#b22222"),
  path_pals = c("#ce7e73", "#f5a596", "#f8d5ce", "#38546d",
    "#cba3b2", "#fbbe85", "#8390ca", "#c0e0db", "#f7a7a6"),
  x_point_left = 0.05,
  x_point_right = 0.24,
  x_path_min = 0.33,
  x_path_max = 0.55,
  x_path_anchor = 0.56,
  x_met_anchor = 0.86,
  x_met_label = 0.89,
  x_met_block_min = 0.86,
  x_met_block_max = 0.99,
  x_met_text = 0.868,
  met_block_fill = "grey80",
  met_block_alt_fill = "grey91",
  met_block_alpha = 0.92,
  met_block_colour = "white",
  met_block_linewidth = 0.18,
  met_text_colour = "black",
  path_label_wrap = 28L,
  met_label_wrap = 28L,
  path_label_size = 2.7,
  met_label_size = 2.2,
  flow_alpha = 0.25,
  flow_width = 0.025,
  path_alpha = 0.88,
  point_size_range = c(2.2, 5.2),
  show_met_node = TRUE,
  met_node_width = 0.04,
  met_node_fill = "grey94",
  met_node_colour = "white",
  met_node_alpha = 0.85,
  show_bottom_label = TRUE,
  show_point_box = TRUE,
  point_box_pad_x = 0.02,
  point_box_pad_y = 0.12,
  point_box_colour = "grey45",
  point_box_linewidth = 0.3,
  y_axis = -0.20,
  y_tick_label = -0.50,
  y_axis_title = -2,
  y_relation_title = -2,
  y_bottom_pad = 0.35,
  legend_position = "left",
  theme = ggplot2::geom_blank())
{
  if (!requireNamespace("ggalluvial", quietly = TRUE)) {
    stop("Package `ggalluvial` is required.", call. = FALSE)
  }

  metaFuns$assert_kegg_results_metaboDiff_and_mebocost(data_kegg_bridge)

  .split_names <- function(x)
  {
    x <- unlist(strsplit(as.character(x), ";\\s*"), use.names = FALSE)
    x <- unique(trimws(x))
    x <- x[!is.na(x) & nzchar(x) & x != "NA"]

    x
  }

  .wrap_text <- function(x, width)
  {
    vapply(as.character(x), function(s) {
      paste(strwrap(s, width = as.integer(width)), collapse = "\n")
    }, character(1L))
  }

  .get_path_color <- function(n)
  {
    if (length(path_pals) >= n) {
      return(path_pals[seq_len(n)])
    }

    grDevices::colorRampPalette(path_pals)(n)
  }

  .make_stack_layout <- function(vec_name, vec_height)
  {
    data_bottom <- data.frame(
      name = rev(vec_name),
      height = rev(vec_height),
      stringsAsFactors = FALSE
    )

    data_bottom$ymin <- c(0, cumsum(utils::head(data_bottom$height, -1L)))
    data_bottom$ymax <- cumsum(data_bottom$height)
    data_bottom$y <- (data_bottom$ymin + data_bottom$ymax) / 2

    data_bottom[match(vec_name, data_bottom$name), , drop = FALSE]
  }

  .scale_point_x <- function(vec_x, vec_range)
  {
    n_span <- diff(vec_range)

    if (!is.finite(n_span) || n_span == 0) {
      return(rep(mean(c(x_point_left, x_point_right)), length(vec_x)))
    }

    x_point_left + (vec_x - vec_range[1L]) / n_span *
      (x_point_right - x_point_left)
  }

  data_bridge <- tibble::as_tibble(data_fella_bridge)

  vec_required <- c(
    "pathway_id", "pathway_name", "mebo_metabolites",
    "metabo_features", "fella_pvalue", "CompoundHits",
    "CompoundsInPathway"
  )

  vec_missing <- setdiff(vec_required, colnames(data_bridge))

  if (length(vec_missing) > 0L) {
    stop(glue::glue(
      "Missing required column(s): {paste(vec_missing, collapse = ', ')}."
    ))
  }

  data_bridge <- data_bridge[
    !is.na(data_bridge$fella_pvalue) &
      data_bridge$fella_pvalue < p_cutoff,
    ,
    drop = FALSE
  ]

  data_bridge$rich_factor <- data_bridge$CompoundHits /
    pmax(data_bridge$CompoundsInPathway, 1L)

  data_bridge$fill_value <- -log10(
    pmax(data_bridge$fella_pvalue, .Machine$double.xmin)
  )

  data_bridge <- dplyr::arrange(
    data_bridge,
    fella_pvalue,
    dplyr::desc(CompoundHits),
    dplyr::desc(rich_factor)
  )

  data_bridge <- utils::head(data_bridge, as.integer(n_pathway))

  if (nrow(data_bridge) == 0L) {
    stop("No FELLA-supported bridge pathway remained for plotting.")
  }

  data_link <- dplyr::bind_rows(lapply(seq_len(nrow(data_bridge)), function(i) {
    vec_mebo <- .split_names(data_bridge$mebo_metabolites[i])
    vec_metabo <- .split_names(data_bridge$metabo_features[i])

    data_mebo <- data.frame(
      pathway_id = data_bridge$pathway_id[i],
      metabolite = vec_mebo,
      evidence_source = "MEBOCOST",
      stringsAsFactors = FALSE
    )

    data_metabo <- data.frame(
      pathway_id = data_bridge$pathway_id[i],
      metabolite = vec_metabo,
      evidence_source = "Metabolomics",
      stringsAsFactors = FALSE
    )

    dplyr::bind_rows(data_mebo, data_metabo)
  }))

  data_link <- data_link[
    !is.na(data_link$metabolite) &
      nzchar(data_link$metabolite) &
      data_link$metabolite != "NA",
    ,
    drop = FALSE
  ]

  data_link <- dplyr::distinct(
    data_link,
    pathway_id,
    metabolite,
    evidence_source,
    .keep_all = TRUE
  )

  data_source_check <- dplyr::group_by(data_link, pathway_id)
  data_source_check <- dplyr::summarise(
    data_source_check,
    has_mebocost = any(evidence_source == "MEBOCOST"),
    has_metabolomics = any(evidence_source == "Metabolomics"),
    .groups = "drop"
  )

  vec_keep_pathway <- data_source_check$pathway_id[
    data_source_check$has_mebocost &
      data_source_check$has_metabolomics
  ]

  data_link <- data_link[data_link$pathway_id %in% vec_keep_pathway, , drop = FALSE]
  data_bridge <- data_bridge[data_bridge$pathway_id %in% vec_keep_pathway, , drop = FALSE]

  if (nrow(data_bridge) == 0L || nrow(data_link) == 0L) {
    stop("No pathway retained both MEBOCOST and Metabolomics support after cleaning.")
  }

  data_bridge$pathway_order <- seq_len(nrow(data_bridge))

  data_link <- dplyr::left_join(
    data_link,
    data_bridge[, c(
      "pathway_id", "pathway_order", "pathway_name",
      "fella_pvalue", "fill_value", "CompoundHits", "rich_factor"
    )],
    by = "pathway_id"
  )

  data_metabolite_rank <- dplyr::group_by(data_link, metabolite, evidence_source)
  data_metabolite_rank <- dplyr::summarise(
    data_metabolite_rank,
    n_pathway = dplyr::n_distinct(pathway_id),
    path_center = mean(pathway_order, na.rm = TRUE),
    best_fella_pvalue = min(fella_pvalue, na.rm = TRUE),
    .groups = "drop"
  )

  data_metabolite_rank <- data_metabolite_rank[
    data_metabolite_rank$n_pathway >= metabolite_min_pathway,
    ,
    drop = FALSE
  ]

  data_metabolite_rank <- dplyr::arrange(
    data_metabolite_rank,
    path_center,
    dplyr::desc(n_pathway),
    best_fella_pvalue,
    evidence_source,
    metabolite
  )

  data_metabolite_rank <- utils::head(
    data_metabolite_rank,
    as.integer(n_metabolite)
  )

  data_link <- dplyr::inner_join(
    data_link,
    data_metabolite_rank[, c("metabolite", "evidence_source", "path_center")],
    by = c("metabolite", "evidence_source")
  )

  if (nrow(data_link) == 0L) {
    stop("No links remained after metabolite filtering.")
  }

  data_path_height <- stats::aggregate(
    list(height = rep(1, nrow(data_link))),
    list(pathway_id = data_link$pathway_id),
    sum
  )

  data_bridge <- data_bridge[data_bridge$pathway_id %in% data_link$pathway_id, , drop = FALSE]
  data_bridge <- dplyr::left_join(data_bridge, data_path_height, by = "pathway_id")
  data_bridge$height[is.na(data_bridge$height)] <- 1

  data_bridge$path_fill <- .get_path_color(nrow(data_bridge))

  data_bridge$kegg_label <- if ("KEGG.id" %in% colnames(data_bridge)) {
    as.character(data_bridge$KEGG.id)
  } else {
    sub("^map", "hsa", data_bridge$pathway_id)
  }

  data_bridge$path_label <- .wrap_text(
    glue::glue("{data_bridge$kegg_label} {data_bridge$pathway_name}"),
    path_label_wrap
  )

  data_link$path_fill <- data_bridge$path_fill[
    match(data_link$pathway_id, data_bridge$pathway_id)
  ]
  data_link$path_label <- data_bridge$path_label[
    match(data_link$pathway_id, data_bridge$pathway_id)
  ]

  data_link$metabolite_label <- glue::glue(
    "[{data_link$evidence_source}] {data_link$metabolite}"
  )
  data_link$metabolite_label <- .wrap_text(
    data_link$metabolite_label,
    met_label_wrap
  )

  data_met_height <- stats::aggregate(
    list(height = rep(1, nrow(data_link))),
    list(
      metabolite_label = data_link$metabolite_label,
      path_center = data_link$path_center
    ),
    sum
  )

  data_met_height <- dplyr::arrange(
    data_met_height,
    path_center,
    metabolite_label
  )

  data_path_layout <- .make_stack_layout(
    data_bridge$path_label,
    data_bridge$height
  )

  data_met_layout <- .make_stack_layout(
    data_met_height$metabolite_label,
    data_met_height$height
  )

  data_bridge$ymin <- data_path_layout$ymin
  data_bridge$ymax <- data_path_layout$ymax
  data_bridge$y <- data_path_layout$y

  data_met_height$ymin <- data_met_layout$ymin
  data_met_height$ymax <- data_met_layout$ymax
  data_met_height$y <- data_met_layout$y
  data_met_height$met_fill <- rep(
    c(met_block_fill, met_block_alt_fill),
    length.out = nrow(data_met_height)
  )

  data_link$alluvium <- paste(
    data_link$path_label,
    data_link$metabolite_label,
    sep = "___"
  )

  data_alluvium <- rbind(
    data.frame(
      alluvium = data_link$alluvium,
      x = x_path_anchor,
      stratum = data_link$path_label,
      value = 1,
      path_fill = data_link$path_fill,
      stringsAsFactors = FALSE
    ),
    data.frame(
      alluvium = data_link$alluvium,
      x = x_met_anchor,
      stratum = data_link$metabolite_label,
      value = 1,
      path_fill = data_link$path_fill,
      stringsAsFactors = FALSE
    )
  )

  vec_stratum_level <- c(
    rev(data_bridge$path_label),
    rev(data_met_height$metabolite_label)
  )

  data_alluvium$stratum <- factor(
    data_alluvium$stratum,
    levels = vec_stratum_level
  )

  vec_point_range <- range(data_bridge$rich_factor, na.rm = TRUE)
  vec_point_break <- pretty(vec_point_range, n = 4L)
  vec_point_break <- vec_point_break[
    vec_point_break >= vec_point_range[1L] &
      vec_point_break <= vec_point_range[2L]
  ]

  if (!length(vec_point_break)) {
    vec_point_break <- vec_point_range
  }

  data_bridge$x_point <- .scale_point_x(data_bridge$rich_factor, vec_point_range)
  vec_point_break_x <- .scale_point_x(vec_point_break, vec_point_range)

  n_y_max <- max(c(data_bridge$ymax, data_met_height$ymax), na.rm = TRUE)

  data_point_box <- data.frame(
    xmin = x_point_left - point_box_pad_x,
    xmax = x_point_right + point_box_pad_x,
    ymin = min(data_bridge$ymin, na.rm = TRUE) - point_box_pad_y,
    ymax = max(data_bridge$ymax, na.rm = TRUE) + point_box_pad_y
  )

  data_point_axis <- data.frame(
    x = vec_point_break_x,
    label = signif(vec_point_break, 3L)
  )

  y_lower <- min(y_axis, y_tick_label, y_axis_title, y_relation_title) -
    y_bottom_pad

  p <- ggplot2::ggplot()

  if (isTRUE(show_point_box)) {
    p <- p +
      ggplot2::geom_rect(
        data = data_point_box,
        ggplot2::aes(
          xmin = xmin,
          xmax = xmax,
          ymin = ymin,
          ymax = ymax
        ),
        inherit.aes = FALSE,
        fill = NA,
        colour = point_box_colour,
        linewidth = point_box_linewidth
      )
  }

  p <- p +
    ggalluvial::geom_alluvium(
      data = data_alluvium,
      ggplot2::aes(
        x = x,
        stratum = stratum,
        alluvium = alluvium,
        y = value,
        fill = path_fill
      ),
      width = flow_width,
      alpha = flow_alpha,
      colour = NA,
      knot.pos = 0.42,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    ggplot2::geom_rect(
      data = data_bridge,
      ggplot2::aes(
        xmin = x_path_min,
        xmax = x_path_max,
        ymin = ymin,
        ymax = ymax,
        fill = path_fill
      ),
      alpha = path_alpha,
      colour = "white",
      linewidth = 0.25,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = data_bridge,
      ggplot2::aes(
        x = x_path_min + 0.012,
        y = y,
        label = path_label
      ),
      hjust = 0,
      vjust = 0.5,
      lineheight = 0.82,
      size = path_label_size,
      colour = "black",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_point(
      data = data_bridge,
      ggplot2::aes(
        x = x_point,
        y = y,
        size = CompoundHits,
        colour = fill_value
      ),
      inherit.aes = FALSE
    )

  if (isTRUE(show_met_node)) {
    p <- p +
      ggplot2::geom_rect(
        data = data_met_height,
        ggplot2::aes(
          xmin = x_met_anchor - met_node_width / 2,
          xmax = x_met_anchor + met_node_width / 2,
          ymin = ymin,
          ymax = ymax
        ),
        inherit.aes = FALSE,
        fill = met_node_fill,
        colour = met_node_colour,
        alpha = met_node_alpha,
        linewidth = 0.18
      )
  }

  p <- p +
    ggplot2::geom_rect(
      data = data_met_height,
      ggplot2::aes(
        xmin = x_met_block_min,
        xmax = x_met_block_max,
        ymin = ymin,
        ymax = ymax,
        fill = met_fill
      ),
      alpha = met_block_alpha,
      colour = met_block_colour,
      linewidth = met_block_linewidth,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = data_met_height,
      ggplot2::aes(
        x = x_met_text,
        y = y,
        label = metabolite_label
      ),
      hjust = 0,
      vjust = 0.5,
      size = met_label_size,
      colour = met_text_colour,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_colour_gradient(
      low = point_pals[1L],
      high = point_pals[2L],
      name = "-log10(FELLA p value)"
    ) +
    ggplot2::scale_size(
      range = point_size_range,
      name = "Compound hits"
    ) +
    ggplot2::scale_x_continuous(
      breaks = NULL,
      labels = NULL,
      limits = c(0.02, max(0.99, x_met_block_max + 0.015)),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(y_lower, n_y_max + 0.6),
      expand = c(0, 0)
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_bw() +
    theme +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      legend.position = legend_position,
      plot.margin = ggplot2::margin(8, 12, 36, 12)
    )

  if (isTRUE(show_bottom_label)) {
    p <- p +
      ggplot2::geom_segment(
        ggplot2::aes(
          x = x_point_left,
          xend = x_point_right,
          y = y_axis,
          yend = y_axis
        ),
        inherit.aes = FALSE,
        linewidth = 0.25,
        colour = "grey45"
      ) +
      ggplot2::geom_segment(
        data = data_point_axis,
        ggplot2::aes(
          x = x,
          xend = x,
          y = y_axis,
          yend = y_axis - 0.055
        ),
        inherit.aes = FALSE,
        linewidth = 0.25,
        colour = "grey45"
      ) +
      ggplot2::geom_text(
        data = data_point_axis,
        ggplot2::aes(
          x = x,
          y = y_tick_label,
          label = label
        ),
        inherit.aes = FALSE,
        size = 2.7,
        colour = "grey35"
      ) +
      ggplot2::annotate(
        "text",
        x = mean(c(x_point_left, x_point_right)),
        y = y_axis_title,
        label = "Rich factor",
        size = 3.2
      ) +
      ggplot2::annotate(
        "text",
        x = mean(c(x_path_min, x_met_label)),
        y = y_relation_title,
        label = "FELLA pathway-metabolite relationship",
        size = 3.0,
        fontface = "bold"
      )
  }

  p
}

metaFuns$plot_bridge_candidate_axis_grouped <- function(data_bridge_candidate,
  n_top_pathway = 6L, n_top_per_pathway = 5L,
  pathway_rank_by = c("bridge_score", "candidate_score", "fella_pvalue"),
  point_pals = c("grey85", "#b22222"),
  point_size_range = c(2.2, 7.0),
  pathway_label_width = 30L,
  metabolite_label_width = 22L,
  receiver_label_angle = 45,
  x_bridge_min = 0,
  x_bridge_max = 1,
  x_receiver_start = 1.8,
  row_gap = 0.7,
  group_gap = 0.9)
{
  pathway_rank_by <- match.arg(pathway_rank_by)
  data_axis <- tibble::as_tibble(data_bridge_candidate)

  vec_required <- c(
    "Metabolite_Name", "Receiver", "candidate_score",
    "pathway_id", "pathway_name", "fella_pvalue",
    "best_metabo_pvalue", "bridge_score"
  )

  vec_missing <- setdiff(vec_required, colnames(data_axis))

  if (length(vec_missing) > 0L) {
    stop(glue::glue(
      "Missing required column(s): {paste(vec_missing, collapse = ', ')}."
    ))
  }

  .wrap_text <- function(x, width)
  {
    vapply(as.character(x), function(s) {
      paste(strwrap(s, width = as.integer(width)), collapse = "\n")
    }, character(1L))
  }

  .scale01 <- function(x)
  {
    x <- suppressWarnings(as.numeric(x))
    vec_range <- range(x, na.rm = TRUE)

    if (!all(is.finite(vec_range)) || diff(vec_range) == 0) {
      return(rep(1, length(x)))
    }

    (x - vec_range[1L]) / diff(vec_range)
  }

  data_axis <- data_axis[
    !is.na(data_axis$Metabolite_Name) &
      nzchar(data_axis$Metabolite_Name) &
      !is.na(data_axis$Receiver) &
      nzchar(data_axis$Receiver) &
      !is.na(data_axis$candidate_score) &
      !is.na(data_axis$bridge_score),
    ,
    drop = FALSE
  ]

  if (nrow(data_axis) == 0L) {
    stop("No valid bridge candidate axis remained for plotting.")
  }

  data_pathway <- dplyr::group_by(data_axis, pathway_id, pathway_name)
  data_pathway <- dplyr::summarise(
    data_pathway,
    bridge_score = max(bridge_score, na.rm = TRUE),
    max_candidate_score = max(candidate_score, na.rm = TRUE),
    min_fella_pvalue = min(fella_pvalue, na.rm = TRUE),
    n_axis = dplyr::n(),
    .groups = "drop"
  )

  if (pathway_rank_by == "bridge_score") {
    data_pathway <- dplyr::arrange(
      data_pathway,
      dplyr::desc(bridge_score),
      min_fella_pvalue,
      dplyr::desc(max_candidate_score)
    )
  } else if (pathway_rank_by == "candidate_score") {
    data_pathway <- dplyr::arrange(
      data_pathway,
      dplyr::desc(max_candidate_score),
      dplyr::desc(bridge_score),
      min_fella_pvalue
    )
  } else {
    data_pathway <- dplyr::arrange(
      data_pathway,
      min_fella_pvalue,
      dplyr::desc(bridge_score),
      dplyr::desc(max_candidate_score)
    )
  }

  if (!is.null(n_top_pathway) && is.finite(n_top_pathway)) {
    data_pathway <- utils::head(data_pathway, as.integer(n_top_pathway))
  }

  data_pathway$pathway_order <- seq_len(nrow(data_pathway))
  data_pathway$pathway_label <- glue::glue(
    "{data_pathway$pathway_id} {data_pathway$pathway_name}"
  )
  data_pathway$pathway_label <- .wrap_text(
    data_pathway$pathway_label,
    pathway_label_width
  )

  data_axis <- data_axis[
    data_axis$pathway_id %in% data_pathway$pathway_id,
    ,
    drop = FALSE
  ]

  data_axis <- dplyr::left_join(
    data_axis,
    data_pathway[, c("pathway_id", "pathway_order", "pathway_label")],
    by = "pathway_id"
  )

  data_axis <- dplyr::group_by(data_axis, pathway_id)
  data_axis <- dplyr::arrange(
    data_axis,
    dplyr::desc(candidate_score),
    fella_pvalue,
    best_metabo_pvalue,
    .by_group = TRUE
  )
  data_axis <- dplyr::mutate(
    data_axis,
    rank_in_pathway = dplyr::row_number()
  )

  if (!is.null(n_top_per_pathway) && is.finite(n_top_per_pathway)) {
    data_axis <- dplyr::filter(
      data_axis,
      rank_in_pathway <= as.integer(n_top_per_pathway)
    )
  }

  data_axis <- dplyr::ungroup(data_axis)

  data_axis <- dplyr::arrange(
    data_axis,
    pathway_order,
    rank_in_pathway
  )

  data_receiver <- stats::aggregate(
    candidate_score ~ Receiver,
    data = data_axis,
    FUN = max
  )
  data_receiver <- dplyr::arrange(
    data_receiver,
    dplyr::desc(candidate_score),
    Receiver
  )
  data_receiver$x_receiver <- seq_len(nrow(data_receiver)) + x_receiver_start

  data_axis <- dplyr::left_join(data_axis, data_receiver[, c("Receiver", "x_receiver")],
    by = "Receiver")

  data_row <- dplyr::group_by(data_axis, pathway_id, pathway_order,
    pathway_label, Metabolite_Name)
  data_row <- dplyr::summarise(
    data_row,
    row_score = max(candidate_score, na.rm = TRUE),
    .groups = "drop"
  )
  data_row <- dplyr::arrange(
    data_row,
    pathway_order,
    dplyr::desc(row_score),
    Metabolite_Name
  )

  data_row$y <- NA_real_
  y_current <- 0

  lst_group <- lapply(seq_len(nrow(data_pathway)), function(i) {
    pathway_id <- data_pathway$pathway_id[i]
    idx <- which(data_row$pathway_id == pathway_id)

    if (length(idx) == 0L) {
      return(NULL)
    }

    y_start <- y_current + group_gap
    y_value <- y_start + seq_along(idx) * row_gap

    y_current <<- max(y_value)
    data_row$y[idx] <<- y_value

    data.frame(
      pathway_id = pathway_id,
      ymin = min(y_value) - row_gap * 0.45,
      ymax = max(y_value) + row_gap * 0.45,
      ymid = mean(range(y_value)),
      stringsAsFactors = FALSE
    )
  })

  data_group <- dplyr::bind_rows(lst_group)

  data_group <- dplyr::left_join(
    data_group,
    data_pathway,
    by = "pathway_id"
  )

  data_group$bridge_score_scaled <- .scale01(data_group$bridge_score)
  data_group$x_bridge_end <- x_bridge_min +
    data_group$bridge_score_scaled * (x_bridge_max - x_bridge_min)

  data_group$group_fill <- rep(
    c("grey96", "grey91"),
    length.out = nrow(data_group)
  )

  data_row$metabolite_label <- .wrap_text(
    data_row$Metabolite_Name,
    metabolite_label_width
  )

  data_axis <- dplyr::left_join(
    data_axis,
    data_row[, c("pathway_id", "Metabolite_Name", "y", "metabolite_label")],
    by = c("pathway_id", "Metabolite_Name")
  )

  x_limit_left <- -1.75
  x_limit_right <- max(data_receiver$x_receiver) + 0.7
  x_path_label <- -1.65
  x_bridge_text <- x_bridge_max + 0.08
  x_metabolite_text <- x_receiver_start + 0.55

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = data_group,
      ggplot2::aes(
        xmin = x_limit_left,
        xmax = x_limit_right,
        ymin = ymin,
        ymax = ymax,
        fill = group_fill
      ),
      colour = NA,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = data_group,
      ggplot2::aes(
        x = x_path_label,
        y = ymid,
        label = pathway_label
      ),
      hjust = 0,
      vjust = 0.5,
      size = 2.8,
      lineheight = 0.84,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_segment(
      data = data_group,
      ggplot2::aes(
        x = x_bridge_min,
        xend = x_bridge_end,
        y = ymid,
        yend = ymid
      ),
      linewidth = 4.2,
      lineend = "round",
      colour = "grey35",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = data_group,
      ggplot2::aes(
        x = x_bridge_text,
        y = ymid,
        label = signif(bridge_score, 3L)
      ),
      hjust = 0,
      vjust = 0.5,
      size = 2.7,
      colour = "grey25",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = data_row,
      ggplot2::aes(
        x = x_metabolite_text,
        y = y,
        label = metabolite_label
      ),
      hjust = 1,
      vjust = 0.5,
      size = 2.7,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_point(
      data = data_axis,
      ggplot2::aes(
        x = x_receiver,
        y = y,
        size = candidate_score,
        colour = candidate_score
      ),
      alpha = 0.9,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_colour_gradient(
      low = point_pals[1L],
      high = point_pals[2L],
      name = "Axis\nCandidateScore"
    ) +
    ggplot2::scale_size_continuous(
      range = point_size_range,
      name = "Axis\nCandidateScore",
      guide = "none"
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(
        mean(c(x_bridge_min, x_bridge_max)),
        data_receiver$x_receiver
      ),
      labels = c(
        "Pathway\nBridgeScore",
        as.character(data_receiver$Receiver)
      ),
      limits = c(x_limit_left, x_limit_right),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_reverse(
      breaks = NULL,
      labels = NULL,
      expand = ggplot2::expansion(mult = c(0.02, 0.04))
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = receiver_label_angle,
        hjust = 1,
        size = 8
      ),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = 9),
      legend.text = ggplot2::element_text(size = 8),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(8, 12, 8, 16)
    )

  p
}

# ==========================================================================
# pathview

metaFuns$resolve_pathview_pathway_id <- function(pathway_id)
{
  pathway_id <- as.character(pathway_id[1L])
  pathway_id <- sub("^hsa", "", pathway_id)
  pathway_id <- sub("^map", "", pathway_id)
  pathway_id <- sub("^path:", "", pathway_id)
  pathway_id <- sub("^hsa", "", pathway_id)
  pathway_id <- sub("^map", "", pathway_id)

  if (!grepl("^[0-9]{5}$", pathway_id)) {
    stop(glue::glue("Invalid KEGG pathway ID: {pathway_id}."))
  }

  list(
    pathview_id = pathway_id,
    map_id = paste0("map", pathway_id),
    hsa_id = paste0("hsa", pathway_id)
  )
}

metaFuns$resolve_pathview_target_pathway <- function(data_bridge_candidate,
  pathway_id = NULL, select_by = c("candidate_score", "bridge_score"))
{
  select_by <- match.arg(select_by)
  data_axis <- tibble::as_tibble(data_bridge_candidate)

  if (!is.null(pathway_id)) {
    return(metaFuns$resolve_pathview_pathway_id(pathway_id))
  }

  vec_required <- c("pathway_id", "candidate_score", "bridge_score")

  vec_missing <- setdiff(vec_required, colnames(data_axis))

  if (length(vec_missing) > 0L) {
    stop(glue::glue(
      "Missing required column(s): {paste(vec_missing, collapse = ', ')}."
    ))
  }

  if (select_by == "candidate_score") {
    data_axis <- dplyr::arrange(
      data_axis,
      dplyr::desc(candidate_score),
      dplyr::desc(bridge_score)
    )
  } else {
    data_axis <- dplyr::arrange(
      data_axis,
      dplyr::desc(bridge_score),
      dplyr::desc(candidate_score)
    )
  }

  metaFuns$resolve_pathview_pathway_id(data_axis$pathway_id[1L])
}

metaFuns$collapse_log2fc_by_kegg <- function(data_map, data_value,
  by_map = "name_original", by_value, value_col,
  pvalue_col = NULL, method = c("weighted_mean_abs", "mean", "max_abs"))
{
  method <- match.arg(method)

  data_map <- tibble::as_tibble(data_map)
  data_value <- tibble::as_tibble(data_value)

  vec_map_need <- c("kegg_id", by_map)
  vec_value_need <- c(by_value, value_col)

  vec_map_missing <- setdiff(vec_map_need, colnames(data_map))
  vec_value_missing <- setdiff(vec_value_need, colnames(data_value))

  if (length(vec_map_missing) > 0L) {
    stop(glue::glue(
      "data_map lacks required column(s): {paste(vec_map_missing, collapse = ', ')}."
    ))
  }

  if (length(vec_value_missing) > 0L) {
    stop(glue::glue(
      "data_value lacks required column(s): {paste(vec_value_missing, collapse = ', ')}."
    ))
  }

  data_map_use <- dplyr::select(
    data_map,
    kegg_id = "kegg_id",
    name_original = !!rlang::sym(by_map)
  )

  if (!is.null(pvalue_col) && pvalue_col %in% colnames(data_value)) {
    data_value_use <- dplyr::select(
      data_value,
      name_original = !!rlang::sym(by_value),
      .value_log2fc = !!rlang::sym(value_col),
      .value_pvalue = !!rlang::sym(pvalue_col)
    )
  } else {
    data_value_use <- dplyr::select(
      data_value,
      name_original = !!rlang::sym(by_value),
      .value_log2fc = !!rlang::sym(value_col)
    )
    data_value_use$.value_pvalue <- NA_real_
  }

  data_map_use$name_original <- as.character(data_map_use$name_original)
  data_value_use$name_original <- as.character(data_value_use$name_original)

  data_use <- dplyr::left_join(
    data_map_use,
    data_value_use,
    by = "name_original"
  )

  data_use$kegg_id <- sub("^cpd:", "", as.character(data_use$kegg_id))
  data_use$value <- suppressWarnings(as.numeric(data_use$.value_log2fc))
  data_use$pvalue_for_rank <- suppressWarnings(as.numeric(data_use$.value_pvalue))

  data_use <- data_use[
    !is.na(data_use$kegg_id) &
      nzchar(data_use$kegg_id) &
      !is.na(data_use$value),
    ,
    drop = FALSE
  ]

  if (nrow(data_use) == 0L) {
    return(tibble::tibble(
      kegg_id = character(0L),
      log2FC = numeric(0L),
      names = character(0L),
      pvalue = numeric(0L)
    ))
  }

  data_out <- dplyr::group_by(data_use, kegg_id)

  data_out <- dplyr::summarise(
    data_out,
    log2FC = {
      vec_value <- value

      if (method == "mean") {
        mean(vec_value, na.rm = TRUE)
      } else if (method == "max_abs") {
        vec_value[which.max(abs(vec_value))]
      } else {
        vec_weight <- abs(vec_value)

        if (all(is.na(vec_weight)) || sum(vec_weight, na.rm = TRUE) == 0) {
          mean(vec_value, na.rm = TRUE)
        } else {
          stats::weighted.mean(vec_value, vec_weight, na.rm = TRUE)
        }
      }
    },
    names = paste(unique(as.character(name_original)), collapse = "; "),
    pvalue = suppressWarnings(min(pvalue_for_rank, na.rm = TRUE)),
    .groups = "drop"
  )

  data_out$pvalue[is.infinite(data_out$pvalue)] <- NA_real_

  data_out
}


metaFuns$resolve_first_existing_col <- function(data, candidates, label)
{
  vec_hit <- candidates[candidates %in% colnames(data)]

  if (length(vec_hit) == 0L) {
    stop(glue::glue(
        "{label} lacks expected column. Candidates: {paste(candidates, collapse = ', ')}."
        ))
  }

  vec_hit[1L]
}


metaFuns$prepare_pathview_cpd_matrix_from_metaInte <- function(x,
  pathway_id = NULL, select_by = c("candidate_score", "bridge_score"),
  mebocost_method = c("weighted_mean_abs", "mean", "max_abs"),
  metabolomics_method = c("max_abs", "weighted_mean_abs", "mean"),
  missing_value = NA_real_)
{
  select_by <- match.arg(select_by)
  mebocost_method <- match.arg(mebocost_method)
  metabolomics_method <- match.arg(metabolomics_method)

  data_candidate <- tibble::as_tibble(x$lst_refine$bridge_candidates)
  ids <- metaFuns$resolve_pathview_target_pathway(
    data_bridge_candidate = data_candidate,
    pathway_id = pathway_id,
    select_by = select_by
  )

  data_kegg_bridge <- x$lst_refine$kegg_bridge
  data_mebo_path <- tibble::as_tibble(data_kegg_bridge$data_mebo_path)
  data_metabo_path <- tibble::as_tibble(data_kegg_bridge$data_metabo_path)

  data_mebo_path <- data_mebo_path[data_mebo_path$pathway_id == ids$map_id, , drop = FALSE]
  data_metabo_path <- data_metabo_path[data_metabo_path$pathway_id == ids$map_id, , drop = FALSE]

  data_mebo_value <- tibble::as_tibble(x$data_sources$mebocost$data_diff_commu)
  data_metabo_value <- tibble::as_tibble(data_kegg_bridge$data_metabo_sig)

  mebo_value_col <- metaFuns$resolve_first_existing_col(
    data_mebo_value,
    c("Log2FC", "log2FC", "log2fc", "logFC", "FC"),
    "MEBOCOST differential communication table"
  )

  metabo_value_col <- metaFuns$resolve_first_existing_col(
    data_metabo_value,
    c("log2FC", "Log2FC", "log2fc", "logFC", "log2FoldChange", "FC"),
    "Metabolomics differential table"
  )

  metabo_name_col <- metaFuns$resolve_first_existing_col(
    data_metabo_value,
    c("feature_name", "Metabolite_Name", "Metabolite", "name", "Name"),
    "Metabolomics differential table"
  )

  metabo_pvalue_col <- NULL

  if (any(c("pvalue", "p.value", "P.Value", "p_val") %in% colnames(data_metabo_value))) {
    metabo_pvalue_col <- metaFuns$resolve_first_existing_col(
      data_metabo_value,
      c("pvalue", "p.value", "P.Value", "p_val"),
      "Metabolomics differential table"
    )
  }

  data_mebo <- metaFuns$collapse_log2fc_by_kegg(
    data_map = data_mebo_path,
    data_value = data_mebo_value,
    by_map = "name_original",
    by_value = "Metabolite_Name",
    value_col = mebo_value_col,
    pvalue_col = NULL,
    method = mebocost_method
  )

  data_metabo <- metaFuns$collapse_log2fc_by_kegg(
    data_map = data_metabo_path,
    data_value = data_metabo_value,
    by_map = "name_original",
    by_value = metabo_name_col,
    value_col = metabo_value_col,
    pvalue_col = metabo_pvalue_col,
    method = metabolomics_method
  )

  vec_kegg <- unique(c(data_mebo$kegg_id, data_metabo$kegg_id))
  vec_kegg <- vec_kegg[!is.na(vec_kegg) & nzchar(vec_kegg)]

  if (length(vec_kegg) == 0L) {
    stop("No KEGG compound ID was available for pathview.")
  }

  mat_cpd <- matrix(
    missing_value,
    nrow = length(vec_kegg),
    ncol = 2L,
    dimnames = list(vec_kegg, c("MEBOCOST_Log2FC", "Metabolomics_Log2FC"))
  )

  mat_cpd[data_mebo$kegg_id, "MEBOCOST_Log2FC"] <- data_mebo$log2FC
  mat_cpd[data_metabo$kegg_id, "Metabolomics_Log2FC"] <- data_metabo$log2FC

  data_evidence <- dplyr::bind_rows(
    dplyr::mutate(data_mebo, source = "MEBOCOST"),
    dplyr::mutate(data_metabo, source = "Metabolomics")
  )

  data_evidence$pathway_id <- ids$map_id
  data_evidence$pathview_id <- ids$pathview_id
  data_evidence$hsa_id <- ids$hsa_id

  list(
    cpd_data = mat_cpd,
    data_evidence = data_evidence,
    pathway_id = ids$map_id,
    pathview_id = ids$pathview_id,
    hsa_id = ids$hsa_id
  )
}

metaFuns$plot_png_file <- function(file_png)
{
  expect_package("png", "0.1.8")

  mat_img <- png::readPNG(file_png)
  grob_img <- grid::rasterGrob(mat_img, interpolate = TRUE)

  ggplot2::ggplot() +
    ggplot2::annotation_custom(
      grob_img,
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf
    ) +
    ggplot2::theme_void()
}

metaFuns$run_pathview_safe <- function(cpd_data, pathway_id,
  out_dir, species = "hsa", out_suffix = NULL,
  kegg_native = TRUE, multi_state = TRUE, clean_old = FALSE, ...)
{
  expect_package("pathview", "1.40.0")

  ids <- metaFuns$resolve_pathview_pathway_id(pathway_id)

  if (is.null(out_suffix)) {
    out_suffix <- glue::glue("pathview_{ids$hsa_id}")
  }

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  out_dir <- normalizePath(out_dir, mustWork = TRUE)

  if (isTRUE(clean_old)) {
    vec_old <- list.files(
      out_dir,
      pattern = glue::glue("{species}{ids$pathview_id}.*{out_suffix}"),
      full.names = TRUE
    )

    if (length(vec_old) > 0L) {
      unlink(vec_old, force = TRUE)
    }
  }

  vec_before <- list.files(out_dir, full.names = TRUE)

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)

  setwd(out_dir)

  out <- withr::with_package("pathview", {
    pathview::pathview(
      gene.data = NULL,
      cpd.data = cpd_data,
      pathway.id = ids$pathview_id,
      species = species,
      gene.idtype = "KEGG",
      cpd.idtype = "kegg",
      out.suffix = out_suffix,
      kegg.dir = out_dir,
      kegg.native = kegg_native,
      multi.state = multi_state,
      same.layer = TRUE,
      ...
    )
  })

  vec_after <- list.files(out_dir, full.names = TRUE)
  vec_created <- setdiff(vec_after, vec_before)

  file_rendered_png <- metaFuns$resolve_pathview_rendered_png(
    out_dir = out_dir,
    species = species,
    pathview_id = ids$pathview_id,
    out_suffix = out_suffix
  )

  p <- NULL

  if (length(file_rendered_png) > 0L) {
    message(glue::glue("Use rendered pathview PNG: {file_rendered_png}"))
    p <- metaFuns$plot_png_file(file_rendered_png)
  } else {
    warning("No rendered pathview PNG with out_suffix was found.")
  }

  list(
    pathview_result = out,
    plot = p,
    png_files = file_rendered_png,
    created_files = vec_created,
    out_dir = out_dir,
    out_suffix = out_suffix,
    pathway_id = ids$map_id,
    pathview_id = ids$pathview_id,
    hsa_id = ids$hsa_id
  )
}

metaFuns$resolve_pathview_rendered_png <- function(out_dir, species,
  pathview_id, out_suffix)
{
  vec_png <- list.files(
    out_dir,
    pattern = "\\.png$",
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(vec_png) == 0L) {
    return(character(0L))
  }

  base_prefix <- glue::glue("{species}{pathview_id}")
  vec_base <- basename(vec_png)

  idx_rendered_multi <- grepl(base_prefix, vec_base, fixed = TRUE) &
    grepl(out_suffix, vec_base, fixed = TRUE) &
    grepl("\\.multi\\.png$", vec_base, ignore.case = TRUE)

  if (any(idx_rendered_multi)) {
    return(vec_png[idx_rendered_multi][1L])
  }

  idx_rendered <- grepl(base_prefix, vec_base, fixed = TRUE) &
    grepl(out_suffix, vec_base, fixed = TRUE) &
    grepl("\\.png$", vec_base, ignore.case = TRUE)

  if (any(idx_rendered)) {
    return(vec_png[idx_rendered][1L])
  }

  idx_pathview <- grepl(base_prefix, vec_base, fixed = TRUE) &
    !grepl(glue::glue("^{base_prefix}\\.png$"), vec_base)

  if (any(idx_pathview)) {
    return(vec_png[idx_pathview][1L])
  }

  character(0L)
}

metaFuns$run_pathview_for_bridge_pathway <- function(x,
  pathway_id = NULL, select_by = c("candidate_score", "bridge_score"),
  out_dir = file.path(x$dir_cache, "pathview"),
  species = "hsa", out_suffix = NULL,
  mebocost_method = c("weighted_mean_abs", "mean", "max_abs"),
  metabolomics_method = c("max_abs", "weighted_mean_abs", "mean"),
  clean_old = FALSE, ...)
{
  select_by <- match.arg(select_by)
  mebocost_method <- match.arg(mebocost_method)
  metabolomics_method <- match.arg(metabolomics_method)

  lst_cpd <- metaFuns$prepare_pathview_cpd_matrix_from_metaInte(
    x = x,
    pathway_id = pathway_id,
    select_by = select_by,
    mebocost_method = mebocost_method,
    metabolomics_method = metabolomics_method
  )

  if (is.null(out_suffix)) {
    out_suffix <- glue::glue(
      "{x@sig}_{lst_cpd$hsa_id}_{select_by}"
    )
    out_suffix <- gsub("[^A-Za-z0-9_\\-]+", "_", out_suffix)
  }

  lst_pathview <- metaFuns$run_pathview_safe(
    cpd_data = lst_cpd$cpd_data,
    pathway_id = lst_cpd$pathview_id,
    out_dir = out_dir,
    species = species,
    out_suffix = out_suffix,
    clean_old = clean_old,
    ...
  )

  lst_pathview$data_evidence <- lst_cpd$data_evidence
  lst_pathview$cpd_data <- lst_cpd$cpd_data

  lst_pathview
}


