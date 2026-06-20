# ==========================================================================
# workflow of metaboDiff
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_metaboDiff <- setClass("job_metaboDiff", 
  contains = c("job"),
  prototype = prototype(
    pg = "metaboDiff",
    info = c(""),
    cite = "",
    method = "",
    tag = "metaboDiff",
    analysis = "代谢物差异分析"
    ))

job_metaboDiff <- function(mt_metabolomics)
{
  if (!is(mt_metabolomics, "mt_metabolomics")) {
    stop('!is(mt_metabolomics, "mt_metabolomics").')
  }
  x <- .job_metaboDiff()
  x$mtObject <- mt_metabolomics
  return(x)
}

setMethod("step0", signature = c(x = "job_metaboDiff"),
  function(x){
    step_message("Prepare your data with function `job_metaboDiff`.")
  })

setMethod("step1", signature = c(x = "job_metaboDiff"),
  function(x,
    expr_name = "data_expr",
    feature_name = "data_feature",
    expr_name_out = "data_expr_processed",
    feature_name_out = "data_feature_processed",
    missing_rate_cutoff = 0.5,
    impute_method = c("half_min", "median", "zero", "none"),
    log_transform = TRUE,
    pseudo_count = NULL,
    filter_zero_var = TRUE,
    scale_method = c("pareto", "auto", "center", "none"),
    force = FALSE,
    verbose = TRUE,
    ...)
  {
    step_message("Preprocess metabolomics data.")

    impute_method <- match.arg(impute_method)
    scale_method <- match.arg(scale_method)

    obj <- x$mtObject
    if (!is(obj, "mt_metabolomics")) {
      stop('!is(x$mtObject, "mt_metabolomics").')
    }

    if (expr_name_out %in% names(obj) && !isTRUE(force)) {
      message(glue::glue(
        "`{expr_name_out}` already exists in `x$mtObject`; use `force = TRUE` to rerun preprocessing."
      ))
    } else {
      obj <- mtFuns$preprocess_object(
        obj,
        expr_name = expr_name,
        feature_name = feature_name,
        expr_name_out = expr_name_out,
        feature_name_out = feature_name_out,
        missing_rate_cutoff = missing_rate_cutoff,
        impute_method = impute_method,
        log_transform = log_transform,
        pseudo_count = pseudo_count,
        filter_zero_var = filter_zero_var,
        scale_method = scale_method,
        verbose = verbose,
        ...
      )
    }

    x$mtObject <- obj

    data_preprocess <- obj$preprocess_log$summary
    n_sample <- data_preprocess$n_sample
    n_feature_raw <- data_preprocess$n_feature_raw
    n_feature_processed <- data_preprocess$n_feature_processed
    n_removed <- n_feature_raw - n_feature_processed
    n_missing_raw <- data_preprocess$n_missing_raw
    n_missing_processed <- data_preprocess$n_missing_processed
    pseudo_count_used <- data_preprocess$pseudo_count

    text_log <- if (isTRUE(log_transform)) {
      glue::glue(
        "随后对峰强度数据进行 log2 转换，转换时加入伪计数 {signif(pseudo_count_used, 4)}，以降低极端高丰度代谢物对后续多变量分析的影响。"
      )
    } else {
      "本步骤未进行 log2 转换。"
    }

    text_scale <- switch(
      scale_method,
      pareto = "最后采用 Pareto scaling 对代谢物强度矩阵进行尺度校正，即在中心化后除以标准差的平方根，以在保留中高丰度代谢物贡献的同时减弱量纲差异的影响。",
      auto = "最后采用 auto scaling 对代谢物强度矩阵进行尺度校正，即对每个代谢物进行中心化并除以标准差，使不同代谢物在后续多变量分析中具有可比尺度。",
      center = "最后对代谢物强度矩阵进行中心化处理，以消除不同代谢物整体水平差异对后续分析的影响。",
      none = "本步骤未进行额外尺度校正。"
    )

    text_impute <- switch(
      impute_method,
      half_min = "对保留代谢物中仍存在的缺失值，采用该代谢物非零最小值的一半进行填补。",
      median = "对保留代谢物中仍存在的缺失值，采用该代谢物在全部样本中的中位数进行填补。",
      zero = "对保留代谢物中仍存在的缺失值，采用 0 进行填补。",
      none = "本步骤未进行缺失值填补。"
    )

    x <- methodAdd(x,
      "首先对代谢组峰强度矩阵进行预处理 (主要使用 R 包 `stats` ⟦pkgInfo('stats')⟧)。以样本为行、代谢物为列构建定量矩阵后，计算每个代谢物在全部样本中的缺失比例，并剔除缺失比例高于 {missing_rate_cutoff} 的代谢物。{text_impute} 同时，剔除方差为 0 的代谢物，以避免其影响后续统计建模。{text_log} {text_scale} 预处理前共纳入 {n_sample} 个样本和 {n_feature_raw} 个代谢物特征，预处理后保留 {n_feature_processed} 个代谢物特征，共移除 {n_removed} 个特征；原始矩阵中缺失值数量为 {n_missing_raw}，预处理后缺失值数量为 {n_missing_processed}。"
    )

    return(x)
  })

setMethod("step2", signature = c(x = "job_metaboDiff"),
  function(x,
    expr_name = "data_expr_processed",
    feature_name = "data_feature_processed",
    group_col = NULL,
    pca_color_col = NULL,
    plsda_color_col = NULL,
    oplsda_color_col = NULL,
    add_ellipse = TRUE,
    plsda_permI = 100L,
    plsda_crossvalI = 7L,
    plsda_predI = 2L,
    oplsda = TRUE,
    oplsda_permI = 100L,
    oplsda_crossvalI = 7L,
    oplsda_predI = 1L,
    oplsda_orthoI = 1L,
    skip_failed_oplsda = TRUE,
    scaleC = "none",
    force = FALSE,
    verbose = TRUE,
    ...)
  {
    step_message("Run PCA, PLS-DA and OPLS-DA.")

    obj <- x$mtObject
    if (!is(obj, "mt_metabolomics")) {
      stop('!is(x$mtObject, "mt_metabolomics").')
    }

    if (is.null(group_col)) {
      group_col <- obj$params$group_col
    }
    if (is.null(group_col)) {
      stop("`group_col` is NULL.")
    }
    if (!group_col %in% names(obj$data_sample)) {
      stop(glue::glue("`{group_col}` was not found in `obj$data_sample`."))
    }

    if (is.null(pca_color_col)) {
      pca_color_col <- group_col
    }
    if (is.null(plsda_color_col)) {
      plsda_color_col <- group_col
    }
    if (is.null(oplsda_color_col)) {
      oplsda_color_col <- group_col
    }

    if (is.null(obj$res_pca) || isTRUE(force)) {
      obj <- mtFuns$run_pca_object(
        obj,
        expr_name = expr_name,
        feature_name = feature_name,
        group_col = group_col,
        center = FALSE,
        scale. = FALSE,
        verbose = verbose
      )
    }

    obj <- mtFuns$add_pca_plot_object(
      obj,
      color_col = pca_color_col,
      shape_col = NULL,
      add_ellipse = add_ellipse
    )

    data_pca_var <- obj$res_pca$data_variance
    pc1_var <- round(data_pca_var$variance_percent[1L], 2L)
    pc2_var <- round(data_pca_var$variance_percent[2L], 2L)

    x <- methodAdd(x,
      "基于预处理后的代谢物定量矩阵，以 R 包 `stats` ⟦pkgInfo('stats')⟧ 采用主成分分析（Principal Component Analysis，PCA）对样本整体代谢谱分布进行无监督降维分析。PCA 不引入分组信息，主要用于观察样本间整体变异趋势、潜在离群样本及不同分组样本在主要变异轴上的分布情况。本分析使用 `{expr_name}` 矩阵进行 PCA，由于数据已在预处理步骤中完成 log2 转换和尺度校正，本步骤未额外进行中心化或标准化处理。前两个主成分 PC1 和 PC2 分别解释 {pc1_var}% 和 {pc2_var}% 的总体变异。"
    )

    p.pca <- obj$plot_pca
    p.pca <- set_lab_legend(
      wrap(p.pca, 6, 4),
      glue::glue("{x@sig} PCA score plot"),
      glue::glue("PCA 得分图|||基于预处理后的代谢物定量矩阵进行主成分分析。每个点代表一个样本，颜色表示 {pca_color_col} 分组，虚线椭圆表示对应分组样本在 PCA 空间中的分布范围。横轴和纵轴分别为 PC1 和 PC2，括号内为对应主成分解释的总体变异比例。")
    )

    if (is.null(obj$res_plsda) || isTRUE(force)) {
      obj <- mtFuns$run_plsda_object(
        obj,
        expr_name = expr_name,
        feature_name = feature_name,
        group_col = group_col,
        predI = plsda_predI,
        permI = plsda_permI,
        crossvalI = plsda_crossvalI,
        scaleC = scaleC,
        verbose = verbose,
        ...
      )
    }

    obj <- mtFuns$add_ropls_plot_object(
      obj,
      model_type = "plsda",
      color_col = plsda_color_col,
      add_ellipse = add_ellipse
    )

    data_plsda_sum <- obj$res_plsda$data_summary
    plsda_r2x <- round(data_plsda_sum$`R2X(cum)`[1L], 4L)
    plsda_r2y <- round(data_plsda_sum$`R2Y(cum)`[1L], 4L)
    plsda_q2 <- round(data_plsda_sum$`Q2(cum)`[1L], 4L)
    plsda_pr2y <- round(data_plsda_sum$pR2Y[1L], 4L)
    plsda_pq2 <- round(data_plsda_sum$pQ2[1L], 4L)

    x <- methodAdd(x,
      "以 R 包 `ropls` ⟦pkgInfo('ropls')⟧ 进一步采用偏最小二乘判别分析（Partial Least Squares Discriminant Analysis，PLS-DA）评估不同 {group_col} 分组样本的代谢谱判别趋势。PLS-DA 为监督降维方法，可在引入分组信息的基础上提取与组间差异相关的潜变量。本分析使用 `{expr_name}` 矩阵构建 PLS-DA 模型，预测成分数设为 {plsda_predI}，交叉验证折数设为 {plsda_crossvalI}，并进行 {plsda_permI} 次随机置换检验以评估模型稳定性及过拟合风险。模型累计 R2X、R2Y 和 Q2 分别为 {plsda_r2x}、{plsda_r2y} 和 {plsda_q2}，置换检验 pR2Y 和 pQ2 分别为 {plsda_pr2y} 和 {plsda_pq2}。"
    )

    p.plsda <- obj$plot_plsda
    p.plsda <- set_lab_legend(
      wrap(p.plsda, 6, 4),
      glue::glue("{x@sig} PLS-DA score plot"),
      glue::glue("PLS-DA 得分图|||基于预处理后的代谢物定量矩阵进行偏最小二乘判别分析。每个点代表一个样本，颜色表示 {plsda_color_col} 分组，虚线椭圆表示对应分组样本在监督判别空间中的分布范围。该图用于展示不同分组样本在 PLS-DA 模型中的分离趋势。")
    )

    p.oplsda <- NULL

    if (isTRUE(oplsda)) {
      res_oplsda <- tryCatch(
        {
          if (is.null(obj$res_oplsda) || isTRUE(force)) {
            obj <- mtFuns$run_oplsda_object(
              obj,
              expr_name = expr_name,
              feature_name = feature_name,
              group_col = group_col,
              predI = oplsda_predI,
              orthoI = oplsda_orthoI,
              permI = oplsda_permI,
              crossvalI = oplsda_crossvalI,
              scaleC = scaleC,
              verbose = verbose,
              ...
            )
          }

          obj <- mtFuns$add_ropls_plot_object(
            obj,
            model_type = "oplsda",
            color_col = oplsda_color_col,
            add_ellipse = add_ellipse
          )

          TRUE
        },
        error = function(e) {
          message(glue::glue("Skip OPLS-DA: {conditionMessage(e)}"))
          FALSE
        }
      )

      if (isTRUE(res_oplsda)) {
        data_oplsda_sum <- obj$res_oplsda$data_summary
        oplsda_r2x <- round(data_oplsda_sum$`R2X(cum)`[1L], 4L)
        oplsda_r2y <- round(data_oplsda_sum$`R2Y(cum)`[1L], 4L)
        oplsda_q2 <- round(data_oplsda_sum$`Q2(cum)`[1L], 4L)
        oplsda_pr2y <- round(data_oplsda_sum$pR2Y[1L], 4L)
        oplsda_pq2 <- round(data_oplsda_sum$pQ2[1L], 4L)

        x <- methodAdd(x, 
          "同时采用正交偏最小二乘判别分析（Orthogonal Partial Least Squares Discriminant Analysis，OPLS-DA）进一步评估 {group_col} 分组相关的代谢谱差异 (使用 R 包 `ropls` ⟦pkgInfo('ropls')⟧)。OPLS-DA 在 PLS-DA 的基础上分离与分组相关的预测成分和与分组无关的正交成分，有助于突出组间判别信息。本分析设置预测成分数为 {oplsda_predI}，正交成分数为 {oplsda_orthoI}，交叉验证折数为 {oplsda_crossvalI}，并进行 {oplsda_permI} 次随机置换检验。模型累计 R2X、R2Y 和 Q2 分别为 {oplsda_r2x}、{oplsda_r2y} 和 {oplsda_q2}，置换检验 pR2Y 和 pQ2 分别为 {oplsda_pr2y} 和 {oplsda_pq2}。"
        )

        p.oplsda <- obj$plot_oplsda
        p.oplsda <- set_lab_legend(
          wrap(p.oplsda, 6, 4),
          glue::glue("{x@sig} OPLS-DA score plot"),
          glue::glue("OPLS-DA 得分图|||基于预处理后的代谢物定量矩阵进行正交偏最小二乘判别分析。每个点代表一个样本，颜色表示 {oplsda_color_col} 分组，虚线椭圆表示对应分组样本在 OPLS-DA 判别空间中的分布范围。横轴通常表示与分组相关的预测成分，纵轴表示与分组无关的正交成分。")
        )
      } else if (!isTRUE(skip_failed_oplsda)) {
        stop("OPLS-DA failed and `skip_failed_oplsda = FALSE`.")
      } else {
        x <- methodAdd(x, 
          "本步骤尝试构建 OPLS-DA 模型，但模型未能稳定建立或未获得有效得分矩阵，因此未将 OPLS-DA 结果纳入后续报告。后续差异代谢物筛选主要基于 PLS-DA 的 VIP 值及单变量统计结果。"
        )
      }
    }

    x$mtObject <- obj

    if (is.null(p.oplsda)) {
      x <- plotsAdd(
        x,
        p.pca = p.pca,
        p.plsda = p.plsda
      )
    } else {
      x <- plotsAdd(
        x,
        p.pca = p.pca,
        p.plsda = p.plsda,
        p.oplsda = p.oplsda
      )
    }

    return(x)
  })


setMethod("step3", signature = c(x = "job_metaboDiff"),
  function(x,
    case_group,
    control_group,
    group_col = NULL,
    expr_name_raw = "data_expr",
    expr_name_log = "data_expr_log2",
    feature_name = "data_feature",
    test_method = c("welch", "student", "wilcox"),
    p_adjust_method = "BH",
    vip_source = c("plsda", "oplsda", "none"),
    vip_cutoff = 1,
    p_cutoff = 0.05,
    padj_cutoff = NULL,
    log2fc_cutoff = 0,
    pseudo_count = NULL,
    volcano_top_n = 10L,
    volcano_top_by = c("pvalue", "VIP", "abs_log2FC"),
    volcano_log2fc_cutoff = NULL,
    force = FALSE,
    verbose = TRUE,
    ...)
  {
    step_message("Differential analysis.")

    test_method <- match.arg(test_method)
    vip_source <- match.arg(vip_source)
    volcano_top_by <- match.arg(volcano_top_by)

    obj <- x$mtObject
    if (!is(obj, "mt_metabolomics")) {
      stop('!is(x$mtObject, "mt_metabolomics").')
    }

    if (is.null(group_col)) {
      group_col <- obj$params$group_col
    }
    if (is.null(group_col)) {
      stop("`group_col` is NULL.")
    }
    if (!group_col %in% names(obj$data_sample)) {
      stop(glue::glue("`{group_col}` was not found in `obj$data_sample`."))
    }

    if (!expr_name_log %in% names(obj) || isTRUE(force)) {
      obj <- mtFuns$add_log_expr_object(
        obj,
        expr_name = expr_name_raw,
        expr_name_out = expr_name_log,
        pseudo_count = pseudo_count,
        verbose = verbose
      )
    }

    obj <- mtFuns$run_diff_object(
      obj,
      expr_name = expr_name_log,
      feature_name = feature_name,
      group_col = group_col,
      case_group = case_group,
      control_group = control_group,
      test_method = test_method,
      p_adjust_method = p_adjust_method,
      vip_source = vip_source,
      vip_cutoff = vip_cutoff,
      p_cutoff = p_cutoff,
      padj_cutoff = padj_cutoff,
      log2fc_cutoff = log2fc_cutoff,
      verbose = verbose
    )

    data_diff <- obj$res_diff$data_diff

    n_feature <- nrow(data_diff)
    n_up <- sum(data_diff$change == "Up", na.rm = TRUE)
    n_down <- sum(data_diff$change == "Down", na.rm = TRUE)
    n_sig <- n_up + n_down
    n_stable <- sum(data_diff$change == "Not significant", na.rm = TRUE)

    text_test <- switch(
      test_method,
      welch = "Welch's t-test",
      student = "Student's t-test",
      wilcox = "Wilcoxon rank-sum test"
    )

    text_padj <- if (is.null(padj_cutoff)) {
      # glue::glue("同时采用 {p_adjust_method} 方法计算多重检验校正后的 adjusted P value，但不作为本步骤主要筛选阈值。")
      ""
    } else {
      glue::glue("同时采用 {p_adjust_method} 方法计算 adjusted P value，并以 adjusted P value < {padj_cutoff} 作为显著性筛选阈值。")
    }

    text_fc <- if (is.null(log2fc_cutoff) || log2fc_cutoff == 0) {
      "log2FC 用于表示差异方向和效应量大小，但不设置额外 fold change 截断值。"
    } else {
      glue::glue("同时要求 |log2FC| ≥ {abs(log2fc_cutoff)}。")
    }

    text_vip <- switch(
      vip_source,
      plsda = "变量重要性投影值（Variable Importance in Projection，VIP）来源于 PLS-DA 模型。",
      oplsda = "变量重要性投影值（Variable Importance in Projection，VIP）来源于 OPLS-DA 模型。",
      none = "本步骤未纳入 VIP 值作为筛选条件。"
    )

    x <- methodAdd(x,
      "在完成代谢物矩阵预处理后，进一步对 {case_group} 与 {control_group} 两组样本进行单变量差异分析 (使用 R 包 `stats` ⟦pkgInfo('stats')⟧)。为避免尺度校正影响 fold change 的计算，本步骤基于原始峰强度矩阵经 log2 转换后的矩阵进行统计分析，并以 {case_group} 相对于 {control_group} 计算 log2FC。组间差异检验采用 {text_test}，并计算每个代谢物的 P value、adjusted P value、FC 和 log2FC。{text_padj} {text_vip} 本研究以 VIP > {vip_cutoff} 且 P < {p_cutoff} 作为差异代谢物的主要探索性筛选标准，{text_fc}"
    )

    vec_front <- c(
      "mt_feature_id", "feature_name", "duplicate_index", "is_duplicate_name",
      "VIP", "mean_case", "mean_control", "median_case", "median_control",
      "log2FC", "FC", "statistic", "pvalue", "padj", "change"
    )
    vec_front <- intersect(vec_front, names(data_diff))
    vec_rest <- setdiff(names(data_diff), vec_front)
    data_diff_report <- data_diff[, c(vec_front, vec_rest), drop = FALSE]

    data_diff_report <- set_lab_legend(
      data_diff_report,
      glue::glue("{x@sig} differential metabolites"),
      glue::glue(
        "差异代谢物分析结果|||该表展示 {case_group} 与 {control_group} 之间的代谢物差异分析结果。log2FC 表示 {case_group} 相对于 {control_group} 的平均 log2 强度差异，FC = 2^log2FC；P value 由 {text_test} 计算获得，adjusted P value 采用 {p_adjust_method} 方法校正；VIP 来源于 {vip_source} 模型。差异代谢物筛选标准为 VIP > {vip_cutoff} 且 P < {p_cutoff}。"
      )
    )

    x <- tablesAdd(x, data_diff_report)

    .collapse_metabolites <- function(data_x, n_top = 10L)
    {
      if (nrow(data_x) == 0L) {
        return("无")
      }

      data_x <- data_x[order(data_x$pvalue, -abs(data_x$log2FC)), , drop = FALSE]
      data_x <- data_x[seq_len(min(n_top, nrow(data_x))), , drop = FALSE]

      vec_name <- as.character(data_x$feature_name)
      vec_name[is.na(vec_name) | vec_name == ""] <- as.character(data_x$mt_feature_id[is.na(vec_name) | vec_name == ""])

      paste(vec_name, collapse = "、")
    }

    data_up <- data_diff[data_diff$change == "Up", , drop = FALSE]
    data_down <- data_diff[data_diff$change == "Down", , drop = FALSE]

    text_up_top <- .collapse_metabolites(data_up, n_top = 10L)
    text_down_top <- .collapse_metabolites(data_down, n_top = 10L)

    .get_feature_names <- function(data_x)
    {
      if (nrow(data_x) == 0L) {
        return(character())
      }
      vec_name <- as.character(data_x$feature_name)
      idx_bad <- is.na(vec_name) | vec_name == ""
      vec_name[idx_bad] <- as.character(data_x$mt_feature_id[idx_bad])
      unique(vec_name[!is.na(vec_name) & vec_name != ""])
    }

    lst_feature <- list(
      Up = .get_feature_names(data_up),
      Down = .get_feature_names(data_down)
    )
    lst_feature <- lst_feature[vapply(lst_feature, length, integer(1L)) > 0L]

    if (length(lst_feature) > 0L) {
      x$.feature <- as_feature(
        lst_feature,
        glue::glue("{case_group} vs {control_group} 差异代谢物"),
        nature = "compounds",
        type = "差异代谢物"
      )
    }

    x <- snapAdd(x, glue::glue(
      "本步骤共纳入 {n_feature} 个代谢物特征进行差异分析。按照 ⟦mark$blue('VIP > {vip_cutoff} 且 P < {p_cutoff}')⟧ 的筛选标准，⟦mark$red('共识别到 {n_sig} 个差异代谢物，其中 {n_up} 个在 {case_group} 组中上调，{n_down} 个在 {case_group} 组中下调')⟧，另有 {n_stable} 个代谢物未达到差异筛选标准。按 P value 排序，{case_group} 组上调代谢物中排名靠前的代谢物包括：{text_up_top}；下调代谢物中排名靠前的代谢物包括：{text_down_top}。"
    ))

    obj$plot_volcano <- mtFuns$plot_diff_volcano(
      obj$res_diff$data_diff,
      label_col = "feature_name",
      log2fc_col = "log2FC",
      p_col = "pvalue",
      vip_col = "VIP",
      vip_cutoff = vip_cutoff,
      p_cutoff = p_cutoff,
      log2fc_cutoff = volcano_log2fc_cutoff,
      top_n = volcano_top_n,
      top_by = volcano_top_by,
      title = glue::glue("{case_group} vs {control_group}"),
      ...
    )

    p.volcano <- obj$plot_volcano
    p.volcano <- set_lab_legend(
      wrap(p.volcano, 7.5, 7.5),
      glue::glue("{x@sig} volcano plot"),
      glue::glue(
        "差异代谢物火山图|||火山图展示 {case_group} 与 {control_group} 之间代谢物差异分布。横轴为 log2FC，表示 {case_group} 相对于 {control_group} 的变化方向和幅度；纵轴为 -log10(P value)，表示统计显著性；点大小表示 VIP 值；颜色表示差异筛选结果。红色点表示在 {case_group} 组中上调的差异代谢物，蓝色点表示在 {case_group} 组中下调的差异代谢物，灰色点表示未达到筛选标准的代谢物。差异筛选标准为 VIP > {vip_cutoff} 且 P < {p_cutoff}。图中分别标注上调和下调代谢物中按 {volcano_top_by} 排序靠前的代谢物，若某一方向差异代谢物不足 {volcano_top_n} 个，则标注实际数量。"
      )
    )

    x$mtObject <- obj

    x <- plotsAdd(x, p.volcano = p.volcano)

    return(x)
  })


# ==========================================================================
# differential metabolites

mtFuns <- new.env(parent = emptyenv())
mtFuns$.stop <- function(text)
{
  stop(as.character(glue::glue(text, .envir = parent.frame())), call. = FALSE)
}

mtFuns$.msg <- function(text, verbose = TRUE)
{
  if (isTRUE(verbose)) {
    message(as.character(glue::glue(text, .envir = parent.frame())))
  }
  invisible(NULL)
}

mtFuns$.as_data_frame <- function(x, arg)
{
  if (!is.data.frame(x)) {
    mtFuns$.stop("`{arg}` must be a data.frame or tibble.")
  }
  data_x <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  names(data_x) <- names(x)
  data_x
}

mtFuns$.as_data_list <- function(x, arg)
{
  if (is.data.frame(x)) {
    return(list(mtFuns$.as_data_frame(x, arg)))
  }
  if (!is.list(x) || length(x) == 0L) {
    mtFuns$.stop("`{arg}` must be a data.frame or a non-empty list of data.frames.")
  }
  vec_ok <- vapply(x, is.data.frame, logical(1L))
  if (!all(vec_ok)) {
    mtFuns$.stop("All elements in `{arg}` must be data.frames.")
  }
  lapply(seq_along(x), function(i) {
    mtFuns$.as_data_frame(x[[ i ]], glue::glue("{arg}[[{i}]]"))
  })
}

mtFuns$.clean_names <- function(vec_x, trim_names = TRUE)
{
  vec_x <- as.character(vec_x)
  if (isTRUE(trim_names)) {
    vec_x <- trimws(vec_x)
  }
  vec_x
}

mtFuns$.clean_id <- function(vec_x, trim_ids = TRUE)
{
  vec_x <- as.character(vec_x)
  if (isTRUE(trim_ids)) {
    vec_x <- trimws(vec_x)
  }
  vec_x[vec_x %in% c("", "NA", "NaN", "NULL", "null")] <- NA_character_
  vec_x
}

mtFuns$.match_one_col <- function(data_x, col, arg)
{
  idx <- which(names(data_x) == col)
  if (length(idx) == 0L) {
    mtFuns$.stop("Column `{col}` was not found in `{arg}`.")
  }
  if (length(idx) > 1L) {
    mtFuns$.stop("Column `{col}` appears more than once in `{arg}`.")
  }
  idx
}

mtFuns$.to_numeric <- function(vec_x, feature_name, batch, numeric_clean = TRUE)
{
  if (is.numeric(vec_x)) {
    return(as.numeric(vec_x))
  }
  vec_chr <- trimws(as.character(vec_x))
  vec_chr[vec_chr %in% c("", "NA", "NaN", "NULL", "null", "-")] <- NA_character_
  if (isTRUE(numeric_clean)) {
    vec_chr <- gsub(",", "", vec_chr, fixed = TRUE)
  }
  vec_num <- suppressWarnings(as.numeric(vec_chr))
  vec_bad <- !is.na(vec_chr) & is.na(vec_num)
  if (any(vec_bad)) {
    bad_value <- unique(vec_chr[vec_bad])[1L]
    mtFuns$.stop(
      "Non-numeric value `{bad_value}` was found in feature `{feature_name}` of batch `{batch}`."
    )
  }
  vec_num
}

mtFuns$.make_batch_names <- function(lst_data, batch_names = NULL)
{
  n_batch <- length(lst_data)
  if (is.null(batch_names)) {
    batch_names <- names(lst_data)
    if (is.null(batch_names) || any(batch_names == "")) {
      batch_names <- sprintf("batch_%03d", seq_len(n_batch))
    }
  }
  batch_names <- as.character(batch_names)
  if (length(batch_names) != n_batch) {
    mtFuns$.stop("`batch_names` must have the same length as `lst_data`.")
  }
  if (any(is.na(batch_names)) || any(batch_names == "")) {
    mtFuns$.stop("`batch_names` contains empty or NA values.")
  }
  if (any(duplicated(batch_names))) {
    mtFuns$.stop("`batch_names` must be unique.")
  }
  batch_names
}

mtFuns$.prepare_peak_table <- function(data_x, id_col, batch,
  trim_names = TRUE, trim_ids = TRUE, numeric_clean = TRUE)
{
  names(data_x) <- mtFuns$.clean_names(names(data_x), trim_names = trim_names)
  idx_id <- mtFuns$.match_one_col(data_x, id_col, glue::glue("peak table of {batch}"))
  vec_id <- mtFuns$.clean_id(data_x[[ idx_id ]], trim_ids = trim_ids)
  if (any(is.na(vec_id))) {
    mtFuns$.stop("Empty sample IDs were found in peak table of batch `{batch}`.")
  }
  if (any(duplicated(vec_id))) {
    dup_id <- unique(vec_id[duplicated(vec_id)])[1L]
    mtFuns$.stop("Duplicated sample ID `{dup_id}` was found in peak table of batch `{batch}`.")
  }
  idx_feature <- setdiff(seq_along(data_x), idx_id)
  if (length(idx_feature) == 0L) {
    mtFuns$.stop("No metabolite columns were found in peak table of batch `{batch}`.")
  }
  vec_feature <- mtFuns$.clean_names(names(data_x)[idx_feature], trim_names = trim_names)
  if (any(is.na(vec_feature)) || any(vec_feature == "")) {
    mtFuns$.stop("Empty metabolite names were found in peak table of batch `{batch}`.")
  }
  mat_x <- matrix(
    NA_real_,
    nrow = nrow(data_x),
    ncol = length(idx_feature),
    dimnames = list(NULL, vec_feature)
  )
  for (j in seq_along(idx_feature)) {
    mat_x[, j] <- mtFuns$.to_numeric(
      data_x[[ idx_feature[j] ]],
      feature_name = vec_feature[j],
      batch = batch,
      numeric_clean = numeric_clean
    )
  }
  vec_dup_index <- as.integer(stats::ave(
    seq_along(vec_feature),
    vec_feature,
    FUN = seq_along
  ))
  n_dup_name <- length(unique(vec_feature[duplicated(vec_feature)]))
  n_dup_col <- sum(vec_feature %in% vec_feature[duplicated(vec_feature)])
  list(
    batch = batch,
    sample_id_original = vec_id,
    sample_id_internal = vec_id,
    mat_raw = mat_x,
    feature_name = vec_feature,
    duplicate_index = vec_dup_index,
    n_features_raw = length(vec_feature),
    n_duplicate_feature_names = n_dup_name,
    n_duplicate_feature_columns = n_dup_col
  )
}

mtFuns$.collapse_matrix <- function(mat_x, vec_group, fun = c("mean", "median", "sum"))
{
  fun <- match.arg(fun)
  vec_unique <- unique(vec_group)
  lst_col <- lapply(vec_unique, function(group_i) {
    mat_i <- mat_x[, vec_group == group_i, drop = FALSE]
    if (ncol(mat_i) == 1L) {
      return(mat_i[, 1L])
    }
    vec_all_na <- rowSums(!is.na(mat_i)) == 0L
    if (fun == "mean") {
      vec_y <- rowMeans(mat_i, na.rm = TRUE)
    } else if (fun == "median") {
      vec_y <- apply(mat_i, 1L, stats::median, na.rm = TRUE)
    } else {
      vec_y <- rowSums(mat_i, na.rm = TRUE)
    }
    vec_y[vec_all_na] <- NA_real_
    vec_y
  })
  mat_out <- do.call(cbind, lst_col)
  colnames(mat_out) <- vec_unique
  rownames(mat_out) <- rownames(mat_x)
  mat_out
}

mtFuns$.resolve_features <- function(lst_peak, duplicate_strategy, feature_id_prefix, verbose = TRUE)
{
  data_col_raw <- do.call(rbind, lapply(seq_along(lst_peak), function(i) {
    data.frame(
      batch_index = i,
      batch = lst_peak[[ i ]]$batch,
      feature_name = lst_peak[[ i ]]$feature_name,
      duplicate_index = lst_peak[[ i ]]$duplicate_index,
      original_col = seq_along(lst_peak[[ i ]]$feature_name),
      stringsAsFactors = FALSE
    )
  }))
  vec_dup_max <- tapply(
    data_col_raw$duplicate_index,
    data_col_raw$feature_name,
    max
  )
  n_dup_feature <- sum(vec_dup_max > 1L)
  if (n_dup_feature > 0L) {
    mtFuns$.msg(
      "{n_dup_feature} duplicated metabolite name(s) were detected. `duplicate_strategy = \"{duplicate_strategy}\"` will be used.",
      verbose = verbose
    )
  }
  if (duplicate_strategy == "keep") {
    data_col <- data_col_raw
    data_col$feature_key <- ifelse(
      vec_dup_max[data_col$feature_name] > 1L,
      paste0(data_col$feature_name, "\rdup", data_col$duplicate_index),
      data_col$feature_name
    )
  } else {
    for (i in seq_along(lst_peak)) {
      mat_i <- mtFuns$.collapse_matrix(
        lst_peak[[ i ]]$mat_raw,
        lst_peak[[ i ]]$feature_name,
        fun = duplicate_strategy
      )
      lst_peak[[ i ]]$mat_raw <- mat_i
      lst_peak[[ i ]]$feature_name <- colnames(mat_i)
      lst_peak[[ i ]]$duplicate_index <- rep(1L, ncol(mat_i))
    }
    data_col <- do.call(rbind, lapply(seq_along(lst_peak), function(i) {
      data.frame(
        batch_index = i,
        batch = lst_peak[[ i ]]$batch,
        feature_name = lst_peak[[ i ]]$feature_name,
        duplicate_index = lst_peak[[ i ]]$duplicate_index,
        original_col = seq_along(lst_peak[[ i ]]$feature_name),
        stringsAsFactors = FALSE
      )
    }))
    data_col$feature_key <- data_col$feature_name
  }
  vec_key <- unique(data_col$feature_key)
  idx_first <- match(vec_key, data_col$feature_key)
  data_first <- data_col[idx_first, , drop = FALSE]
  vec_feature_id <- sprintf("%s%06d", feature_id_prefix, seq_along(vec_key))
  vec_n_batch <- tapply(data_col$batch, data_col$feature_key, function(x) {
    length(unique(x))
  })
  data_feature <- data.frame(
    mt_feature_id = vec_feature_id,
    feature_name = data_first$feature_name,
    duplicate_index = data_first$duplicate_index,
    duplicate_n_global = as.integer(vec_dup_max[data_first$feature_name]),
    is_duplicate_name = as.integer(vec_dup_max[data_first$feature_name]) > 1L,
    n_batches_detected = as.integer(vec_n_batch[vec_key]),
    has_missing_batch = as.integer(vec_n_batch[vec_key]) < length(lst_peak),
    feature_key = vec_key,
    stringsAsFactors = FALSE
  )
  for (i in seq_along(lst_peak)) {
    if (duplicate_strategy == "keep") {
      vec_key_i <- ifelse(
        vec_dup_max[lst_peak[[ i ]]$feature_name] > 1L,
        paste0(lst_peak[[ i ]]$feature_name, "\rdup", lst_peak[[ i ]]$duplicate_index),
        lst_peak[[ i ]]$feature_name
      )
    } else {
      vec_key_i <- lst_peak[[ i ]]$feature_name
    }
    vec_id_i <- data_feature$mt_feature_id[match(vec_key_i, data_feature$feature_key)]
    if (any(is.na(vec_id_i))) {
      mtFuns$.stop("Feature ID mapping failed in batch `{lst_peak[[ i ]]$batch}`.")
    }
    if (any(duplicated(vec_id_i))) {
      mtFuns$.stop("Duplicated internal feature IDs were generated in batch `{lst_peak[[ i ]]$batch}`.")
    }
    colnames(lst_peak[[ i ]]$mat_raw) <- vec_id_i
    lst_peak[[ i ]]$mat <- lst_peak[[ i ]]$mat_raw
  }
  list(
    lst_peak = lst_peak,
    data_feature = data_feature
  )
}

mtFuns$.prepare_meta_table <- function(meta_x, id_col, group_col, batch,
  sample_id_original, sample_id_internal, require_meta = TRUE,
  trim_names = TRUE, trim_ids = TRUE)
{
  names(meta_x) <- mtFuns$.clean_names(names(meta_x), trim_names = trim_names)
  vec_reserved <- c("mt_sample_id", "mt_batch", "mt_sample_index")
  if (any(names(meta_x) %in% vec_reserved)) {
    bad_col <- names(meta_x)[names(meta_x) %in% vec_reserved][1L]
    mtFuns$.stop("Metadata column `{bad_col}` is reserved. Please rename it before creating the object.")
  }
  if (any(duplicated(names(meta_x)))) {
    bad_col <- names(meta_x)[duplicated(names(meta_x))][1L]
    mtFuns$.stop("Duplicated metadata column `{bad_col}` was found in batch `{batch}`.")
  }
  idx_id <- mtFuns$.match_one_col(meta_x, id_col, glue::glue("metadata of {batch}"))
  if (!is.null(group_col)) {
    mtFuns$.match_one_col(meta_x, group_col, glue::glue("metadata of {batch}"))
  }
  vec_id <- mtFuns$.clean_id(meta_x[[ idx_id ]], trim_ids = trim_ids)
  meta_x[[ idx_id ]] <- vec_id
  if (any(is.na(vec_id))) {
    mtFuns$.stop("Empty sample IDs were found in metadata of batch `{batch}`.")
  }
  if (any(duplicated(vec_id))) {
    dup_id <- unique(vec_id[duplicated(vec_id)])[1L]
    mtFuns$.stop("Duplicated sample ID `{dup_id}` was found in metadata of batch `{batch}`.")
  }
  vec_missing <- setdiff(sample_id_original, vec_id)
  vec_extra <- setdiff(vec_id, sample_id_original)
  if (length(vec_missing) > 0L && isTRUE(require_meta)) {
    miss_id <- vec_missing[1L]
    mtFuns$.stop("Sample `{miss_id}` in peak table of batch `{batch}` was not found in metadata.")
  }
  idx_match <- match(sample_id_original, vec_id)
  data_meta <- meta_x[idx_match, , drop = FALSE]
  data_meta[[ idx_id ]] <- sample_id_original
  data_out <- cbind(
    data.frame(
      mt_sample_id = sample_id_internal,
      mt_batch = batch,
      mt_sample_index = seq_along(sample_id_original),
      stringsAsFactors = FALSE
    ),
    data_meta,
    stringsAsFactors = FALSE
  )
  rownames(data_out) <- sample_id_internal
  list(
    data_sample = data_out,
    n_meta_samples = nrow(meta_x),
    n_missing_meta = length(vec_missing),
    n_extra_meta = length(vec_extra)
  )
}

mtFuns$.rbind_fill <- function(lst_data)
{
  vec_col <- unique(unlist(lapply(lst_data, names)))
  lst_out <- lapply(lst_data, function(data_x) {
    vec_missing <- setdiff(vec_col, names(data_x))
    if (length(vec_missing) > 0L) {
      for (col_i in vec_missing) {
        data_x[[ col_i ]] <- NA
      }
    }
    data_x[, vec_col, drop = FALSE]
  })
  do.call(rbind, lst_out)
}

mtFuns$.combine_matrices <- function(lst_peak, data_feature)
{
  vec_feature_id <- data_feature$mt_feature_id
  lst_mat <- lapply(seq_along(lst_peak), function(i) {
    mat_i <- lst_peak[[ i ]]$mat
    mat_full <- matrix(
      NA_real_,
      nrow = nrow(mat_i),
      ncol = length(vec_feature_id),
      dimnames = list(lst_peak[[ i ]]$sample_id_internal, vec_feature_id)
    )
    idx_col <- match(colnames(mat_i), vec_feature_id)
    mat_full[, idx_col] <- mat_i
    mat_full
  })
  do.call(rbind, lst_mat)
}

mtFuns$is_object <- function(x)
{
  inherits(x, "mt_metabolomics")
}

mtFuns$check_object <- function(x, strict = TRUE)
{
  if (!mtFuns$is_object(x)) {
    if (isTRUE(strict)) {
      mtFuns$.stop("Input is not a `mt_metabolomics` object.")
    }
    return(FALSE)
  }
  vec_need <- c("data_expr", "data_sample", "data_feature", "params", "input_log")
  vec_missing <- setdiff(vec_need, names(x))
  if (length(vec_missing) > 0L) {
    if (isTRUE(strict)) {
      mtFuns$.stop("Object is missing required element `{vec_missing[1L]}`.")
    }
    return(FALSE)
  }
  if (!is.matrix(x$data_expr) || !is.numeric(x$data_expr)) {
    if (isTRUE(strict)) {
      mtFuns$.stop("`data_expr` must be a numeric matrix.")
    }
    return(FALSE)
  }
  if (!identical(rownames(x$data_expr), x$data_sample$mt_sample_id)) {
    if (isTRUE(strict)) {
      mtFuns$.stop("Rows of `data_expr` do not match `data_sample$mt_sample_id`.")
    }
    return(FALSE)
  }
  if (!identical(colnames(x$data_expr), x$data_feature$mt_feature_id)) {
    if (isTRUE(strict)) {
      mtFuns$.stop("Columns of `data_expr` do not match `data_feature$mt_feature_id`.")
    }
    return(FALSE)
  }
  invisible(TRUE)
}

mtFuns$create_object <- function(lst_data, lst_meta,
  id_col = "PatientID", group_col = NULL, group_levels = NULL,
  batch_names = NULL, duplicate_strategy = c("keep", "mean", "median", "sum"),
  sample_id_strategy = c("stop", "make_unique"),
  feature_id_prefix = "M", require_meta = TRUE,
  trim_names = TRUE, trim_ids = TRUE, numeric_clean = TRUE,
  verbose = TRUE)
{
  duplicate_strategy <- match.arg(duplicate_strategy)
  sample_id_strategy <- match.arg(sample_id_strategy)
  lst_data <- mtFuns$.as_data_list(lst_data, "lst_data")
  lst_meta <- mtFuns$.as_data_list(lst_meta, "lst_meta")
  if (length(lst_data) != length(lst_meta)) {
    mtFuns$.stop("`lst_data` and `lst_meta` must have the same length.")
  }
  batch_names <- mtFuns$.make_batch_names(lst_data, batch_names = batch_names)
  lst_peak <- lapply(seq_along(lst_data), function(i) {
    mtFuns$.prepare_peak_table(
      lst_data[[ i ]],
      id_col = id_col,
      batch = batch_names[i],
      trim_names = trim_names,
      trim_ids = trim_ids,
      numeric_clean = numeric_clean
    )
  })
  vec_all_id <- unlist(lapply(lst_peak, function(x) x$sample_id_original), use.names = FALSE)
  vec_all_batch <- unlist(lapply(lst_peak, function(x) {
    rep(x$batch, length(x$sample_id_original))
  }), use.names = FALSE)
  if (any(duplicated(vec_all_id))) {
    if (sample_id_strategy == "stop") {
      dup_id <- unique(vec_all_id[duplicated(vec_all_id)])[1L]
      mtFuns$.stop(
        "Duplicated sample ID `{dup_id}` was found across batches. Use `sample_id_strategy = \"make_unique\"` only if this is expected."
      )
    }
    vec_internal <- paste0(vec_all_id, "__", vec_all_batch)
  } else {
    vec_internal <- vec_all_id
  }
  idx_start <- cumsum(c(1L, vapply(lst_peak, function(x) {
    length(x$sample_id_original)
  }, integer(1L))[-length(lst_peak)]))
  for (i in seq_along(lst_peak)) {
    n_i <- length(lst_peak[[ i ]]$sample_id_original)
    idx_i <- seq.int(idx_start[i], length.out = n_i)
    lst_peak[[ i ]]$sample_id_internal <- vec_internal[idx_i]
    rownames(lst_peak[[ i ]]$mat_raw) <- vec_internal[idx_i]
  }
  res_feature <- mtFuns$.resolve_features(
    lst_peak,
    duplicate_strategy = duplicate_strategy,
    feature_id_prefix = feature_id_prefix,
    verbose = verbose
  )
  lst_peak <- res_feature$lst_peak
  data_feature <- res_feature$data_feature
  lst_meta_res <- lapply(seq_along(lst_meta), function(i) {
    mtFuns$.prepare_meta_table(
      lst_meta[[ i ]],
      id_col = id_col,
      group_col = group_col,
      batch = batch_names[i],
      sample_id_original = lst_peak[[ i ]]$sample_id_original,
      sample_id_internal = lst_peak[[ i ]]$sample_id_internal,
      require_meta = require_meta,
      trim_names = trim_names,
      trim_ids = trim_ids
    )
  })
  data_sample <- mtFuns$.rbind_fill(lapply(lst_meta_res, function(x) x$data_sample))
  mat_expr <- mtFuns$.combine_matrices(lst_peak, data_feature)
  if (!identical(rownames(mat_expr), data_sample$mt_sample_id)) {
    mtFuns$.stop("Internal error: sample order is inconsistent after merging.")
  }
  if (!is.null(group_col)) {
    if (!group_col %in% names(data_sample)) {
      mtFuns$.stop("Group column `{group_col}` was not found after metadata merging.")
    }
    vec_group <- as.character(data_sample[[ group_col ]])
    if (any(is.na(vec_group)) || any(vec_group == "")) {
      mtFuns$.stop("Group column `{group_col}` contains empty or NA values.")
    }
    if (!is.null(group_levels)) {
      vec_unknown <- setdiff(unique(vec_group), group_levels)
      if (length(vec_unknown) > 0L) {
        mtFuns$.stop("Group value `{vec_unknown[1L]}` is not included in `group_levels`.")
      }
      data_sample[[ group_col ]] <- factor(vec_group, levels = group_levels)
    } else {
      data_sample[[ group_col ]] <- factor(vec_group)
    }
    n_group <- length(unique(data_sample[[ group_col ]]))
    if (n_group < 2L) {
      mtFuns$.stop("Group column `{group_col}` must contain at least 2 groups.")
    }
  }
  data_log <- data.frame(
    batch = batch_names,
    n_peak_samples = vapply(lst_peak, function(x) length(x$sample_id_original), integer(1L)),
    n_meta_samples = vapply(lst_meta_res, function(x) x$n_meta_samples, integer(1L)),
    n_features_raw = vapply(lst_peak, function(x) x$n_features_raw, integer(1L)),
    n_features_used = vapply(lst_peak, function(x) ncol(x$mat), integer(1L)),
    n_duplicate_feature_names = vapply(lst_peak, function(x) x$n_duplicate_feature_names, integer(1L)),
    n_duplicate_feature_columns = vapply(lst_peak, function(x) x$n_duplicate_feature_columns, integer(1L)),
    n_missing_meta = vapply(lst_meta_res, function(x) x$n_missing_meta, integer(1L)),
    n_extra_meta = vapply(lst_meta_res, function(x) x$n_extra_meta, integer(1L)),
    stringsAsFactors = FALSE
  )
  obj <- list(
    data_expr = mat_expr,
    data_sample = data_sample,
    data_feature = data_feature,
    params = list(
      id_col = id_col,
      group_col = group_col,
      group_levels = group_levels,
      duplicate_strategy = duplicate_strategy,
      sample_id_strategy = sample_id_strategy,
      feature_id_prefix = feature_id_prefix
    ),
    input_log = data_log
  )
  class(obj) <- c("mt_metabolomics", "list")
  mtFuns$check_object(obj)
  mtFuns$.msg(
    "Created mt_metabolomics object: {nrow(mat_expr)} samples, {ncol(mat_expr)} features, {length(batch_names)} batch(es).",
    verbose = verbose
  )
  obj
}

mtFuns$get_expr <- function(x, label = c("id", "name", "name_id"))
{
  label <- match.arg(label)
  mtFuns$check_object(x)
  mat_x <- x$data_expr
  if (label == "id") {
    return(mat_x)
  }
  if (label == "name") {
    colnames(mat_x) <- make.unique(x$data_feature$feature_name, sep = "__dup")
  } else {
    colnames(mat_x) <- paste0(
      x$data_feature$feature_name,
      "|",
      x$data_feature$mt_feature_id
    )
  }
  mat_x
}

mtFuns$get_sample <- function(x)
{
  mtFuns$check_object(x)
  x$data_sample
}

mtFuns$get_feature <- function(x)
{
  mtFuns$check_object(x)
  x$data_feature
}

mtFuns$summarize_object <- function(x)
{
  mtFuns$check_object(x)
  data_summary <- data.frame(
    n_samples = nrow(x$data_expr),
    n_features = ncol(x$data_expr),
    n_batches = length(unique(x$data_sample$mt_batch)),
    n_duplicate_feature_names = sum(x$data_feature$is_duplicate_name),
    n_missing_values = sum(is.na(x$data_expr)),
    missing_rate = sum(is.na(x$data_expr)) / length(x$data_expr),
    stringsAsFactors = FALSE
  )
  print(data_summary)
  if (!is.null(x$params$group_col)) {
    print(table(x$data_sample[[ x$params$group_col ]], useNA = "ifany"))
  }
  invisible(data_summary)
}

mtFuns$.get_expr_matrix <- function(x, expr_name = "data_expr")
{
  mtFuns$check_object(x)

  if (!expr_name %in% names(x)) {
    mtFuns$.stop("Expression matrix `{expr_name}` was not found in the object.")
  }

  mat_expr <- x[[ expr_name ]]

  if (!is.matrix(mat_expr)) {
    mtFuns$.stop("`{expr_name}` must be a matrix.")
  }
  if (!is.numeric(mat_expr)) {
    mtFuns$.stop("`{expr_name}` must be numeric.")
  }
  if (!identical(rownames(mat_expr), x$data_sample$mt_sample_id)) {
    mtFuns$.stop("Rows of `{expr_name}` do not match `data_sample$mt_sample_id`.")
  }

  mat_expr
}

mtFuns$.get_feature_table <- function(x, feature_name = "data_feature", mat_expr = NULL)
{
  if (!feature_name %in% names(x)) {
    mtFuns$.stop("Feature table `{feature_name}` was not found in the object.")
  }

  data_feature <- x[[ feature_name ]]

  if (!is.data.frame(data_feature)) {
    mtFuns$.stop("`{feature_name}` must be a data.frame.")
  }

  if (!is.null(mat_expr)) {
    if (nrow(data_feature) != ncol(mat_expr)) {
      mtFuns$.stop("Rows of `{feature_name}` do not match columns of expression matrix.")
    }

    if ("mt_feature_id" %in% names(data_feature)) {
      if (!identical(as.character(data_feature$mt_feature_id), colnames(mat_expr))) {
        mtFuns$.stop("`{feature_name}$mt_feature_id` does not match expression matrix columns.")
      }
    }
  }

  data_feature
}

mtFuns$.filter_feature_missing <- function(mat_expr, data_feature = NULL,
  missing_rate_cutoff = 0.5, verbose = TRUE)
{
  if (!is.numeric(missing_rate_cutoff) || length(missing_rate_cutoff) != 1L ||
      is.na(missing_rate_cutoff) || missing_rate_cutoff < 0 || missing_rate_cutoff > 1) {
    mtFuns$.stop("`missing_rate_cutoff` must be a number between 0 and 1.")
  }

  vec_missing_rate <- colMeans(is.na(mat_expr))
  idx_keep <- vec_missing_rate <= missing_rate_cutoff

  if (!any(idx_keep)) {
    mtFuns$.stop("No feature remained after missing-rate filtering.")
  }

  data_missing <- data.frame(
    mt_feature_id = colnames(mat_expr),
    missing_rate = vec_missing_rate,
    keep = idx_keep,
    stringsAsFactors = FALSE
  )

  mat_out <- mat_expr[, idx_keep, drop = FALSE]

  if (!is.null(data_feature)) {
    data_feature <- data_feature[idx_keep, , drop = FALSE]
  }

  n_removed <- sum(!idx_keep)
  mtFuns$.msg(
    "Removed {n_removed} feature(s) with missing rate > {missing_rate_cutoff}.",
    verbose = verbose
  )

  list(
    mat_expr = mat_out,
    data_feature = data_feature,
    data_missing = data_missing
  )
}

mtFuns$.impute_matrix <- function(mat_expr,
  impute_method = c("half_min", "median", "zero", "none"),
  verbose = TRUE)
{
  impute_method <- match.arg(impute_method)
  n_missing <- sum(is.na(mat_expr))

  if (n_missing == 0L) {
    data_impute <- data.frame(
      mt_feature_id = colnames(mat_expr),
      n_missing = 0L,
      impute_value = NA_real_,
      stringsAsFactors = FALSE
    )

    mtFuns$.msg("No missing value needs imputation.", verbose = verbose)

    return(list(
      mat_expr = mat_expr,
      data_impute = data_impute
    ))
  }

  if (impute_method == "none") {
    mtFuns$.stop("Missing values remain but `impute_method = \"none\"`.")
  }

  mat_out <- mat_expr
  vec_n_missing <- colSums(is.na(mat_out))
  vec_impute_value <- rep(NA_real_, ncol(mat_out))

  for (j in seq_len(ncol(mat_out))) {
    if (vec_n_missing[j] == 0L) {
      next
    }

    vec_x <- mat_out[, j]
    vec_valid <- vec_x[!is.na(vec_x)]

    if (length(vec_valid) == 0L) {
      mtFuns$.stop("Feature `{colnames(mat_out)[j]}` has no valid value for imputation.")
    }

    if (impute_method == "half_min") {
      vec_pos <- vec_valid[vec_valid > 0]
      if (length(vec_pos) > 0L) {
        value <- min(vec_pos) / 2
      } else if (all(vec_valid == 0)) {
        value <- 0
      } else {
        mtFuns$.stop(
          "Feature `{colnames(mat_out)[j]}` has no positive value for half-min imputation."
        )
      }
    } else if (impute_method == "median") {
      value <- stats::median(vec_valid)
    } else {
      value <- 0
    }

    mat_out[is.na(vec_x), j] <- value
    vec_impute_value[j] <- value
  }

  data_impute <- data.frame(
    mt_feature_id = colnames(mat_out),
    n_missing = vec_n_missing,
    impute_value = vec_impute_value,
    stringsAsFactors = FALSE
  )

  mtFuns$.msg(
    "Imputed {n_missing} missing value(s) using `{impute_method}`.",
    verbose = verbose
  )

  list(
    mat_expr = mat_out,
    data_impute = data_impute
  )
}

mtFuns$.log2_transform_matrix <- function(mat_expr,
  log_transform = TRUE, pseudo_count = NULL, verbose = TRUE)
{
  if (!isTRUE(log_transform)) {
    return(list(
      mat_expr = mat_expr,
      pseudo_count = NA_real_
    ))
  }

  if (any(mat_expr < 0, na.rm = TRUE)) {
    mtFuns$.stop("Negative values were found before log2 transformation.")
  }

  if (is.null(pseudo_count)) {
    vec_positive <- mat_expr[mat_expr > 0 & !is.na(mat_expr)]

    if (length(vec_positive) == 0L) {
      mtFuns$.stop("No positive value was found for automatic pseudo-count calculation.")
    }

    pseudo_count <- min(vec_positive) / 2
  }

  if (!is.numeric(pseudo_count) || length(pseudo_count) != 1L ||
      is.na(pseudo_count) || pseudo_count < 0) {
    mtFuns$.stop("`pseudo_count` must be a non-negative single numeric value.")
  }

  if (any(mat_expr + pseudo_count <= 0, na.rm = TRUE)) {
    mtFuns$.stop("Non-positive values were found after adding pseudo-count.")
  }

  mat_out <- log2(mat_expr + pseudo_count)

  mtFuns$.msg(
    "Applied log2 transformation with pseudo_count = {signif(pseudo_count, 4)}.",
    verbose = verbose
  )

  list(
    mat_expr = mat_out,
    pseudo_count = pseudo_count
  )
}

mtFuns$.filter_feature_sd <- function(mat_expr, data_feature = NULL,
  sd_cutoff = 0, verbose = TRUE)
{
  vec_sd <- apply(mat_expr, 2L, stats::sd, na.rm = TRUE)
  idx_keep <- !is.na(vec_sd) & vec_sd > sd_cutoff

  if (!any(idx_keep)) {
    mtFuns$.stop("No feature remained after SD filtering.")
  }

  data_sd <- data.frame(
    mt_feature_id = colnames(mat_expr),
    sd = vec_sd,
    keep = idx_keep,
    stringsAsFactors = FALSE
  )

  mat_out <- mat_expr[, idx_keep, drop = FALSE]

  if (!is.null(data_feature)) {
    data_feature <- data_feature[idx_keep, , drop = FALSE]
  }

  n_removed <- sum(!idx_keep)

  mtFuns$.msg(
    "Removed {n_removed} low-variance feature(s) with SD <= {sd_cutoff}.",
    verbose = verbose
  )

  list(
    mat_expr = mat_out,
    data_feature = data_feature,
    data_sd = data_sd
  )
}

mtFuns$.scale_matrix <- function(mat_expr,
  scale_method = c("pareto", "auto", "center", "none"),
  verbose = TRUE)
{
  scale_method <- match.arg(scale_method)

  if (scale_method == "none") {
    mtFuns$.msg("No scaling was applied.", verbose = verbose)

    return(list(
      mat_expr = mat_expr,
      center = NA,
      scale = NA
    ))
  }

  vec_center <- colMeans(mat_expr, na.rm = TRUE)
  mat_out <- sweep(mat_expr, 2L, vec_center, FUN = "-")

  if (scale_method == "center") {
    mtFuns$.msg("Applied center scaling.", verbose = verbose)

    return(list(
      mat_expr = mat_out,
      center = vec_center,
      scale = NA
    ))
  }

  vec_sd <- apply(mat_expr, 2L, stats::sd, na.rm = TRUE)
  vec_bad <- is.na(vec_sd) | vec_sd <= 0

  if (any(vec_bad)) {
    bad_feature <- colnames(mat_expr)[which(vec_bad)[1L]]
    mtFuns$.stop("Feature `{bad_feature}` has invalid SD during scaling.")
  }

  if (scale_method == "auto") {
    vec_scale <- vec_sd
  } else {
    vec_scale <- sqrt(vec_sd)
  }

  mat_out <- sweep(mat_out, 2L, vec_scale, FUN = "/")

  mtFuns$.msg("Applied `{scale_method}` scaling.", verbose = verbose)

  list(
    mat_expr = mat_out,
    center = vec_center,
    scale = vec_scale
  )
}

mtFuns$preprocess_object <- function(x,
  expr_name = "data_expr",
  feature_name = "data_feature",
  expr_name_out = "data_expr_processed",
  feature_name_out = "data_feature_processed",
  missing_rate_cutoff = 0.5,
  impute_method = c("half_min", "median", "zero", "none"),
  log_transform = TRUE,
  pseudo_count = NULL,
  filter_zero_var = TRUE,
  scale_method = c("pareto", "auto", "center", "none"),
  verbose = TRUE)
{
  mtFuns$check_object(x)

  impute_method <- match.arg(impute_method)
  scale_method <- match.arg(scale_method)

  mat_raw <- mtFuns$.get_expr_matrix(x, expr_name = expr_name)
  data_feature <- mtFuns$.get_feature_table(
    x,
    feature_name = feature_name,
    mat_expr = mat_raw
  )

  n_sample_raw <- nrow(mat_raw)
  n_feature_raw <- ncol(mat_raw)
  n_missing_raw <- sum(is.na(mat_raw))

  mtFuns$.msg(
    "Start preprocessing `{expr_name}`: {n_sample_raw} samples, {n_feature_raw} features.",
    verbose = verbose
  )

  res_missing <- mtFuns$.filter_feature_missing(
    mat_raw,
    data_feature = data_feature,
    missing_rate_cutoff = missing_rate_cutoff,
    verbose = verbose
  )

  mat_x <- res_missing$mat_expr
  data_feature <- res_missing$data_feature

  res_impute <- mtFuns$.impute_matrix(
    mat_x,
    impute_method = impute_method,
    verbose = verbose
  )

  mat_x <- res_impute$mat_expr

  res_log <- mtFuns$.log2_transform_matrix(
    mat_x,
    log_transform = log_transform,
    pseudo_count = pseudo_count,
    verbose = verbose
  )

  mat_x <- res_log$mat_expr

  if (isTRUE(filter_zero_var)) {
    res_sd <- mtFuns$.filter_feature_sd(
      mat_x,
      data_feature = data_feature,
      sd_cutoff = 0,
      verbose = verbose
    )

    mat_x <- res_sd$mat_expr
    data_feature <- res_sd$data_feature
    data_sd <- res_sd$data_sd
  } else {
    data_sd <- data.frame(
      mt_feature_id = colnames(mat_x),
      sd = apply(mat_x, 2L, stats::sd, na.rm = TRUE),
      keep = TRUE,
      stringsAsFactors = FALSE
    )
  }

  res_scale <- mtFuns$.scale_matrix(
    mat_x,
    scale_method = scale_method,
    verbose = verbose
  )

  mat_x <- res_scale$mat_expr

  if (!identical(colnames(mat_x), as.character(data_feature$mt_feature_id))) {
    mtFuns$.stop("Processed matrix columns do not match processed feature table.")
  }

  x[[ expr_name_out ]] <- mat_x
  x[[ feature_name_out ]] <- data_feature

  data_summary <- data.frame(
    expr_name_in = expr_name,
    expr_name_out = expr_name_out,
    n_sample = nrow(mat_x),
    n_feature_raw = n_feature_raw,
    n_feature_processed = ncol(mat_x),
    n_missing_raw = n_missing_raw,
    n_missing_processed = sum(is.na(mat_x)),
    missing_rate_cutoff = missing_rate_cutoff,
    impute_method = impute_method,
    log_transform = log_transform,
    pseudo_count = res_log$pseudo_count,
    filter_zero_var = filter_zero_var,
    scale_method = scale_method,
    stringsAsFactors = FALSE
  )

  x$preprocess_log <- list(
    summary = data_summary,
    missing_filter = res_missing$data_missing,
    imputation = res_impute$data_impute,
    sd_filter = data_sd,
    scale_center = res_scale$center,
    scale_factor = res_scale$scale
  )

  if (is.null(x$preprocess_history)) {
    x$preprocess_history <- list()
  }

  x$preprocess_history[[ length(x$preprocess_history) + 1L ]] <- x$preprocess_log

  class(x) <- c("mt_metabolomics", "list")

  mtFuns$.msg(
    "Finished preprocessing: {nrow(mat_x)} samples, {ncol(mat_x)} features.",
    verbose = verbose
  )

  x
}

mtFuns$.resolve_group_col <- function(x, group_col = NULL)
{
  mtFuns$check_object(x)

  if (is.null(group_col)) {
    group_col <- x$params$group_col
  }
  if (is.null(group_col)) {
    mtFuns$.stop("`group_col` is NULL. Please provide a grouping column.")
  }
  if (!group_col %in% names(x$data_sample)) {
    mtFuns$.stop("Group column `{group_col}` was not found in `data_sample`.")
  }

  group_col
}

mtFuns$.get_pc_variance <- function(res_pca)
{
  vec_var <- res_pca$sdev ^ 2
  vec_var_ratio <- vec_var / sum(vec_var)

  data.frame(
    pc = paste0("PC", seq_along(vec_var)),
    variance = vec_var,
    variance_ratio = vec_var_ratio,
    variance_percent = vec_var_ratio * 100,
    cumulative_percent = cumsum(vec_var_ratio) * 100,
    stringsAsFactors = FALSE
  )
}

mtFuns$run_pca_object <- function(x,
  expr_name = "data_expr_processed",
  feature_name = "data_feature_processed",
  group_col = NULL,
  center = FALSE,
  scale. = FALSE,
  verbose = TRUE)
{
  mtFuns$check_object(x)

  group_col <- mtFuns$.resolve_group_col(x, group_col = group_col)

  mat_expr <- mtFuns$.get_expr_matrix(x, expr_name = expr_name)
  data_feature <- mtFuns$.get_feature_table(
    x,
    feature_name = feature_name,
    mat_expr = mat_expr
  )

  if (any(is.na(mat_expr))) {
    mtFuns$.stop("`{expr_name}` contains NA values. Please preprocess or impute first.")
  }

  mtFuns$.msg(
    "Run PCA using `{expr_name}`: {nrow(mat_expr)} samples, {ncol(mat_expr)} features.",
    verbose = verbose
  )

  res_pca <- stats::prcomp(
    mat_expr,
    center = center,
    scale. = scale.
  )

  data_variance <- mtFuns$.get_pc_variance(res_pca)

  data_score <- as.data.frame(res_pca$x, stringsAsFactors = FALSE)
  data_score$mt_sample_id <- rownames(data_score)
  data_score <- merge(
    x$data_sample,
    data_score,
    by = "mt_sample_id",
    all.x = TRUE,
    sort = FALSE
  )
  data_score <- data_score[match(rownames(mat_expr), data_score$mt_sample_id), , drop = FALSE]

  data_loading <- as.data.frame(res_pca$rotation, stringsAsFactors = FALSE)
  data_loading$mt_feature_id <- rownames(data_loading)
  data_loading <- merge(
    data_feature,
    data_loading,
    by = "mt_feature_id",
    all.x = TRUE,
    sort = FALSE
  )
  data_loading <- data_loading[match(colnames(mat_expr), data_loading$mt_feature_id), , drop = FALSE]

  x$res_pca <- list(
    model = res_pca,
    data_score = data_score,
    data_loading = data_loading,
    data_variance = data_variance,
    expr_name = expr_name,
    feature_name = feature_name,
    group_col = group_col,
    center = center,
    scale. = scale.
  )

  mtFuns$.msg(
    "Finished PCA: PC1 = {round(data_variance$variance_percent[1L], 2)}%, PC2 = {round(data_variance$variance_percent[2L], 2)}%.",
    verbose = verbose
  )

  x
}

mtFuns$plot_pca_score <- function(x,
  pc_x = "PC1",
  pc_y = "PC2",
  color_col = NULL,
  shape_col = NULL,
  add_ellipse = TRUE,
  point_size = 3,
  alpha = 0.85,
  title = NULL)
{
  mtFuns$check_object(x)

  if (is.null(x$res_pca)) {
    mtFuns$.stop("`x$res_pca` was not found. Please run `mtFuns$run_pca_object()` first.")
  }

  data_score <- x$res_pca$data_score
  data_variance <- x$res_pca$data_variance

  if (!pc_x %in% names(data_score)) {
    mtFuns$.stop("`{pc_x}` was not found in PCA score table.")
  }
  if (!pc_y %in% names(data_score)) {
    mtFuns$.stop("`{pc_y}` was not found in PCA score table.")
  }

  if (is.null(color_col)) {
    color_col <- x$res_pca$group_col
  }
  if (!color_col %in% names(data_score)) {
    mtFuns$.stop("Color column `{color_col}` was not found in PCA score table.")
  }
  if (!is.null(shape_col) && !shape_col %in% names(data_score)) {
    mtFuns$.stop("Shape column `{shape_col}` was not found in PCA score table.")
  }

  idx_x <- match(pc_x, data_variance$pc)
  idx_y <- match(pc_y, data_variance$pc)

  xlab <- glue::glue("{pc_x} ({round(data_variance$variance_percent[idx_x], 2)}%)")
  ylab <- glue::glue("{pc_y} ({round(data_variance$variance_percent[idx_y], 2)}%)")

  if (is.null(title)) {
    title <- "PCA score plot"
  }

  if (is.null(shape_col)) {
    p <- ggplot2::ggplot(
      data_score,
      ggplot2::aes(
        x = !!rlang::sym(pc_x),
        y = !!rlang::sym(pc_y),
        color = !!rlang::sym(color_col)
      )
    )
  } else {
    p <- ggplot2::ggplot(
      data_score,
      ggplot2::aes(
        x = !!rlang::sym(pc_x),
        y = !!rlang::sym(pc_y),
        color = !!rlang::sym(color_col),
        shape = !!rlang::sym(shape_col)
      )
    )
  }

  p <- p +
    ggplot2::geom_point(size = point_size, alpha = alpha) +
    ggplot2::labs(
      title = title,
      x = xlab,
      y = ylab,
      color = color_col,
      shape = shape_col
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 9)
    )

  if (isTRUE(add_ellipse)) {
    data_group_count <- as.data.frame(table(data_score[[ color_col ]]), stringsAsFactors = FALSE)

    if (all(data_group_count$Freq >= 3L)) {
      p <- p + ggplot2::stat_ellipse(
        ggplot2::aes(group = !!rlang::sym(color_col)),
        type = "norm",
        linewidth = 0.5,
        linetype = 2,
        show.legend = FALSE
      )
    }
  }

  p
}

mtFuns$add_pca_plot_object <- function(x,
  pc_x = "PC1",
  pc_y = "PC2",
  color_col = NULL,
  shape_col = NULL,
  add_ellipse = TRUE)
{
  mtFuns$check_object(x)

  x$plot_pca <- mtFuns$plot_pca_score(
    x,
    pc_x = pc_x,
    pc_y = pc_y,
    color_col = color_col,
    shape_col = shape_col,
    add_ellipse = add_ellipse
  )

  x
}

mtFuns$.check_ropls <- function()
{
  if (!requireNamespace("ropls", quietly = TRUE)) {
    mtFuns$.stop(
      "Package `ropls` is required. Please install it with `BiocManager::install(\"ropls\")`."
    )
  }

  invisible(TRUE)
}

mtFuns$.extract_ropls_summary <- function(model)
{
  data_summary <- tryCatch(
    ropls::getSummaryDF(model),
    error = function(e) NULL
  )

  if (is.null(data_summary)) {
    data_summary <- data.frame(stringsAsFactors = FALSE)
  }

  data_summary
}

mtFuns$run_ropls_da_object <- function(x,
  model_type = c("plsda", "oplsda"),
  expr_name = "data_expr_processed",
  feature_name = "data_feature_processed",
  group_col = NULL,
  predI = NULL,
  orthoI = NULL,
  permI = 100L,
  crossvalI = 7L,
  scaleC = "none",
  verbose = TRUE)
{
  mtFuns$check_object(x)
  mtFuns$.check_ropls()

  model_type <- match.arg(model_type)
  group_col <- mtFuns$.resolve_group_col(x, group_col = group_col)

  mat_expr <- mtFuns$.get_expr_matrix(x, expr_name = expr_name)
  data_feature <- mtFuns$.get_feature_table(
    x,
    feature_name = feature_name,
    mat_expr = mat_expr
  )

  if (any(is.na(mat_expr))) {
    mtFuns$.stop("`{expr_name}` contains NA values. Please preprocess or impute first.")
  }

  vec_group <- x$data_sample[[ group_col ]]

  if (any(is.na(vec_group))) {
    mtFuns$.stop("Group column `{group_col}` contains NA values.")
  }

  vec_group <- factor(vec_group)

  if (nlevels(vec_group) != 2L) {
    mtFuns$.stop("Currently `{model_type}` expects exactly 2 groups.")
  }

  if (is.null(predI)) {
    if (model_type == "plsda") {
      predI <- 2L
    } else {
      predI <- 1L
    }
  }

  if (is.null(orthoI)) {
    if (model_type == "plsda") {
      orthoI <- 0L
    } else {
      orthoI <- 1L
    }
  }

  mtFuns$.msg(
    "Run {toupper(model_type)} using `{expr_name}` and group `{group_col}`.",
    verbose = verbose
  )

  model <- ropls::opls(
    mat_expr,
    vec_group,
    predI = predI,
    orthoI = orthoI,
    permI = permI,
    crossvalI = crossvalI,
    scaleC = scaleC,
    fig.pdfC = "none",
    info.txtC = "none"
  )

  data_summary <- mtFuns$.extract_ropls_summary(model)

  mat_check <- tryCatch(
    ropls::getScoreMN(model, orthoL = FALSE),
    error = function(e) NULL
  )

  if (mtFuns$.is_empty_ropls_matrix(mat_check)) {
    mtFuns$.stop(
      "ropls did not build a valid {toupper(model_type)} model. Try reducing model complexity or use PLS-DA only."
    )
  }

  data_score <- mtFuns$.extract_ropls_score(
    model,
    data_sample = x$data_sample,
    model_type = model_type
  )

  data_loading <- mtFuns$.extract_ropls_loading(
    model,
    data_feature = data_feature
  )

  data_vip_raw <- mtFuns$.extract_ropls_vip(model)

  if (nrow(data_vip_raw) > 0L) {
    if (nrow(data_vip_raw) == nrow(data_feature) &&
        (any(is.na(data_vip_raw$mt_feature_id)) ||
         any(data_vip_raw$mt_feature_id == "") ||
         any(!data_vip_raw$mt_feature_id %in% data_feature$mt_feature_id))) {
      data_vip_raw$mt_feature_id <- data_feature$mt_feature_id
    }

    data_vip <- data_feature
    data_vip$VIP <- data_vip_raw$VIP[
      match(data_feature$mt_feature_id, data_vip_raw$mt_feature_id)
    ]
  } else {
    data_vip <- data_feature
    data_vip$VIP <- NA_real_
  }

  res_name <- paste0("res_", model_type)

  x[[ res_name ]] <- list(
    model = model,
    data_summary = data_summary,
    data_score = data_score,
    data_loading = data_loading,
    data_vip = data_vip,
    expr_name = expr_name,
    feature_name = feature_name,
    group_col = group_col,
    model_type = model_type,
    predI = predI,
    orthoI = orthoI,
    permI = permI,
    crossvalI = crossvalI,
    scaleC = scaleC
  )

  mtFuns$.msg(
    "Finished {toupper(model_type)}.",
    verbose = verbose
  )

  x
}

mtFuns$run_plsda_object <- function(x, ...)
{
  mtFuns$run_ropls_da_object(
    x,
    model_type = "plsda",
    ...
  )
}

mtFuns$run_oplsda_object <- function(x, ...)
{
  mtFuns$run_ropls_da_object(
    x,
    model_type = "oplsda",
    ...
  )
}

mtFuns$plot_ropls_score <- function(x,
  model_type = c("plsda", "oplsda"),
  comp_x = NULL,
  comp_y = NULL,
  color_col = NULL,
  add_ellipse = TRUE,
  point_size = 3,
  alpha = 0.85,
  title = NULL)
{
  mtFuns$check_object(x)

  model_type <- match.arg(model_type)
  res_name <- paste0("res_", model_type)

  if (is.null(x[[ res_name ]])) {
    mtFuns$.stop("`x${res_name}` was not found. Please run `{model_type}` first.")
  }

  res_model <- x[[ res_name ]]
  data_score <- res_model$data_score

  if (is.null(color_col)) {
    color_col <- res_model$group_col
  }

  if (!color_col %in% names(data_score)) {
    mtFuns$.stop("Color column `{color_col}` was not found in score table.")
  }

  if (is.null(comp_x)) {
    comp_x <- "t1"
  }

  if (is.null(comp_y)) {
    if (model_type == "oplsda" && "to1" %in% names(data_score)) {
      comp_y <- "to1"
    } else {
      comp_y <- "t2"
    }
  }

  if (!comp_x %in% names(data_score)) {
    mtFuns$.stop("Component `{comp_x}` was not found in score table.")
  }

  if (!comp_y %in% names(data_score)) {
    mtFuns$.stop("Component `{comp_y}` was not found in score table.")
  }

  if (is.null(title)) {
    title <- paste0(toupper(model_type), " score plot")
  }

  p <- ggplot2::ggplot(
    data_score,
    ggplot2::aes(
      x = !!rlang::sym(comp_x),
      y = !!rlang::sym(comp_y),
      color = !!rlang::sym(color_col)
    )
  ) +
    ggplot2::geom_point(size = point_size, alpha = alpha) +
    ggplot2::labs(
      title = title,
      x = comp_x,
      y = comp_y,
      color = color_col
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 9)
    )

  if (isTRUE(add_ellipse)) {
    data_group_count <- as.data.frame(table(data_score[[ color_col ]]), stringsAsFactors = FALSE)

    if (all(data_group_count$Freq >= 3L)) {
      p <- p + ggplot2::stat_ellipse(
        ggplot2::aes(group = !!rlang::sym(color_col)),
        type = "norm",
        linewidth = 0.5,
        linetype = 2,
        show.legend = FALSE
      )
    }
  }

  p
}

mtFuns$add_ropls_plot_object <- function(x,
  model_type = c("plsda", "oplsda"),
  comp_x = NULL,
  comp_y = NULL,
  color_col = NULL,
  add_ellipse = TRUE)
{
  mtFuns$check_object(x)

  model_type <- match.arg(model_type)
  plot_name <- paste0("plot_", model_type)

  x[[ plot_name ]] <- mtFuns$plot_ropls_score(
    x,
    model_type = model_type,
    comp_x = comp_x,
    comp_y = comp_y,
    color_col = color_col,
    add_ellipse = add_ellipse
  )

  x
}

# ==========================================================================

mtFuns$.as_ropls_matrix <- function(x, value_name)
{
  if (is.null(x)) {
    return(NULL)
  }

  if (is.vector(x) && !is.list(x)) {
    mat_x <- matrix(as.numeric(x), ncol = 1L)
    rownames(mat_x) <- names(x)
    colnames(mat_x) <- value_name
    return(mat_x)
  }

  if (is.data.frame(x)) {
    mat_x <- as.matrix(x)
    return(mat_x)
  }

  if (is.matrix(x)) {
    return(x)
  }

  mtFuns$.stop("Cannot convert `{value_name}` from ropls result to matrix.")
}

mtFuns$.has_valid_ids <- function(vec_id, vec_expected_id)
{
  !is.null(vec_id) &&
    !any(is.na(vec_id)) &&
    !any(vec_id == "") &&
    !any(duplicated(vec_id)) &&
    all(vec_id %in% vec_expected_id)
}

mtFuns$.resolve_ropls_row_id <- function(mat_x, vec_expected_id, value_name)
{
  if (is.null(dim(mat_x))) {
    mat_x <- matrix(as.numeric(mat_x), ncol = 1L)
  }

  vec_id <- rownames(mat_x)

  if (mtFuns$.has_valid_ids(vec_id, vec_expected_id)) {
    return(list(
      mat = mat_x,
      id = vec_id
    ))
  }

  if (nrow(mat_x) == length(vec_expected_id)) {
    rownames(mat_x) <- vec_expected_id

    return(list(
      mat = mat_x,
      id = vec_expected_id
    ))
  }

  if (ncol(mat_x) == length(vec_expected_id)) {
    mat_x <- t(mat_x)
    rownames(mat_x) <- vec_expected_id

    return(list(
      mat = mat_x,
      id = vec_expected_id
    ))
  }

  actual_dim <- paste(dim(mat_x), collapse = " x ")

  mtFuns$.stop(
    "Cannot resolve row IDs for `{value_name}`: matrix dimension is {actual_dim}, expected sample number is {length(vec_expected_id)}."
  )
}

mtFuns$.resolve_ropls_feature_id <- function(mat_x, vec_expected_id, value_name)
{
  if (is.null(dim(mat_x))) {
    mat_x <- matrix(as.numeric(mat_x), ncol = 1L)
  }

  vec_id <- rownames(mat_x)

  if (mtFuns$.has_valid_ids(vec_id, vec_expected_id)) {
    return(list(
      mat = mat_x,
      id = vec_id
    ))
  }

  if (nrow(mat_x) == length(vec_expected_id)) {
    rownames(mat_x) <- vec_expected_id

    return(list(
      mat = mat_x,
      id = vec_expected_id
    ))
  }

  if (ncol(mat_x) == length(vec_expected_id)) {
    mat_x <- t(mat_x)
    rownames(mat_x) <- vec_expected_id

    return(list(
      mat = mat_x,
      id = vec_expected_id
    ))
  }

  actual_dim <- paste(dim(mat_x), collapse = " x ")

  mtFuns$.stop(
    "Cannot resolve feature IDs for `{value_name}`: matrix dimension is {actual_dim}, expected feature number is {length(vec_expected_id)}."
  )
}

mtFuns$.is_empty_ropls_matrix <- function(mat_x)
{
  is.null(mat_x) || length(mat_x) == 0L || any(dim(mat_x) == 0L)
}

mtFuns$.extract_ropls_score <- function(model, data_sample, model_type)
{
  mat_score <- tryCatch(
    ropls::getScoreMN(model, orthoL = FALSE),
    error = function(e) NULL
  )

  if (mtFuns$.is_empty_ropls_matrix(mat_score)) {
    mtFuns$.stop(
      "No predictive score was extracted from ropls model. The model was probably not successfully built."
    )
  }

  mat_score <- as.matrix(mat_score)

  if (nrow(mat_score) != nrow(data_sample)) {
    mtFuns$.stop(
      "Predictive score row number is {nrow(mat_score)}, but sample number is {nrow(data_sample)}."
    )
  }

  rownames(mat_score) <- data_sample$mt_sample_id
  colnames(mat_score) <- paste0("t", seq_len(ncol(mat_score)))

  data_score <- cbind(
    data_sample,
    as.data.frame(mat_score, stringsAsFactors = FALSE),
    stringsAsFactors = FALSE
  )

  mat_ortho <- tryCatch(
    ropls::getScoreMN(model, orthoL = TRUE),
    error = function(e) NULL
  )

  if (!mtFuns$.is_empty_ropls_matrix(mat_ortho)) {
    mat_ortho <- as.matrix(mat_ortho)

    if (nrow(mat_ortho) == nrow(data_sample)) {
      rownames(mat_ortho) <- data_sample$mt_sample_id
      colnames(mat_ortho) <- paste0("to", seq_len(ncol(mat_ortho)))

      data_score <- cbind(
        data_score,
        as.data.frame(mat_ortho, stringsAsFactors = FALSE),
        stringsAsFactors = FALSE
      )
    }
  }

  rownames(data_score) <- data_score$mt_sample_id
  data_score
}

mtFuns$.extract_ropls_loading <- function(model, data_feature)
{
  mat_loading <- tryCatch(
    ropls::getLoadingMN(model, orthoL = FALSE),
    error = function(e) NULL
  )

  if (mtFuns$.is_empty_ropls_matrix(mat_loading)) {
    return(data.frame(
      mt_feature_id = data_feature$mt_feature_id,
      stringsAsFactors = FALSE
    ))
  }

  mat_loading <- as.matrix(mat_loading)

  if (nrow(mat_loading) != nrow(data_feature)) {
    mtFuns$.stop(
      "Predictive loading row number is {nrow(mat_loading)}, but feature number is {nrow(data_feature)}."
    )
  }

  rownames(mat_loading) <- data_feature$mt_feature_id
  colnames(mat_loading) <- paste0("p", seq_len(ncol(mat_loading)))

  data_loading <- cbind(
    data_feature,
    as.data.frame(mat_loading, stringsAsFactors = FALSE),
    stringsAsFactors = FALSE
  )

  mat_ortho_loading <- tryCatch(
    ropls::getLoadingMN(model, orthoL = TRUE),
    error = function(e) NULL
  )

  if (!mtFuns$.is_empty_ropls_matrix(mat_ortho_loading)) {
    mat_ortho_loading <- as.matrix(mat_ortho_loading)

    if (nrow(mat_ortho_loading) == nrow(data_feature)) {
      rownames(mat_ortho_loading) <- data_feature$mt_feature_id
      colnames(mat_ortho_loading) <- paste0("po", seq_len(ncol(mat_ortho_loading)))

      data_loading <- cbind(
        data_loading,
        as.data.frame(mat_ortho_loading, stringsAsFactors = FALSE),
        stringsAsFactors = FALSE
      )
    }
  }

  rownames(data_loading) <- data_loading$mt_feature_id
  data_loading
}

mtFuns$.extract_ropls_vip <- function(model)
{
  vec_vip <- tryCatch(
    ropls::getVipVn(model, orthoL = FALSE),
    error = function(e) NULL
  )

  if (is.null(vec_vip) || length(vec_vip) == 0L) {
    return(data.frame(
      mt_feature_id = character(),
      VIP = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    mt_feature_id = names(vec_vip),
    VIP = as.numeric(vec_vip),
    stringsAsFactors = FALSE
  )
}

mtFuns$add_log_expr_object <- function(x,
  expr_name = "data_expr",
  expr_name_out = "data_expr_log2",
  pseudo_count = NULL,
  verbose = TRUE)
{
  mtFuns$check_object(x)

  mat_expr <- mtFuns$.get_expr_matrix(x, expr_name = expr_name)

  res_log <- mtFuns$.log2_transform_matrix(
    mat_expr,
    log_transform = TRUE,
    pseudo_count = pseudo_count,
    verbose = verbose
  )

  x[[ expr_name_out ]] <- res_log$mat_expr

  x$log_expr_info <- list(
    expr_name = expr_name,
    expr_name_out = expr_name_out,
    pseudo_count = res_log$pseudo_count
  )

  class(x) <- c("mt_metabolomics", "list")
  x
}

mtFuns$run_diff_object <- function(x,
  expr_name = "data_expr_log2",
  feature_name = "data_feature",
  group_col = NULL,
  case_group = "TP53mut",
  control_group = "TP53wt",
  test_method = c("welch", "student", "wilcox"),
  p_adjust_method = "BH",
  vip_source = c("plsda", "oplsda", "none"),
  vip_cutoff = 1,
  p_cutoff = 0.05,
  padj_cutoff = NULL,
  log2fc_cutoff = 0,
  verbose = TRUE)
{
  mtFuns$check_object(x)

  test_method <- match.arg(test_method)
  vip_source <- match.arg(vip_source)
  group_col <- mtFuns$.resolve_group_col(x, group_col = group_col)

  mat_expr <- mtFuns$.get_expr_matrix(x, expr_name = expr_name)
  data_feature <- mtFuns$.get_feature_table(
    x,
    feature_name = feature_name,
    mat_expr = mat_expr
  )

  vec_group <- as.character(x$data_sample[[ group_col ]])

  if (!all(c(case_group, control_group) %in% vec_group)) {
    mtFuns$.stop("Both `case_group` and `control_group` must exist in `{group_col}`.")
  }

  idx_case <- vec_group == case_group
  idx_control <- vec_group == control_group

  mtFuns$.msg(
    "Run differential analysis: {case_group} (n = {sum(idx_case)}) vs {control_group} (n = {sum(idx_control)}).",
    verbose = verbose
  )

  .test_one <- function(vec_x)
  {
    vec_case <- vec_x[idx_case]
    vec_control <- vec_x[idx_control]

    mean_case <- mean(vec_case, na.rm = TRUE)
    mean_control <- mean(vec_control, na.rm = TRUE)
    median_case <- stats::median(vec_case, na.rm = TRUE)
    median_control <- stats::median(vec_control, na.rm = TRUE)

    if (test_method == "welch") {
      res_test <- tryCatch(
        stats::t.test(vec_case, vec_control, var.equal = FALSE),
        error = function(e) NULL
      )
    } else if (test_method == "student") {
      res_test <- tryCatch(
        stats::t.test(vec_case, vec_control, var.equal = TRUE),
        error = function(e) NULL
      )
    } else {
      res_test <- tryCatch(
        stats::wilcox.test(vec_case, vec_control, exact = FALSE),
        error = function(e) NULL
      )
    }

    if (is.null(res_test)) {
      pvalue <- NA_real_
      statistic <- NA_real_
    } else {
      pvalue <- res_test$p.value
      statistic <- as.numeric(res_test$statistic[1L])
    }

    c(
      mean_case = mean_case,
      mean_control = mean_control,
      median_case = median_case,
      median_control = median_control,
      log2FC = mean_case - mean_control,
      statistic = statistic,
      pvalue = pvalue
    )
  }

  mat_stat <- t(apply(mat_expr, 2L, .test_one))

  data_diff <- cbind(
    data_feature,
    as.data.frame(mat_stat, stringsAsFactors = FALSE),
    stringsAsFactors = FALSE
  )

  data_diff$FC <- 2 ^ data_diff$log2FC
  data_diff$padj <- stats::p.adjust(data_diff$pvalue, method = p_adjust_method)

  if (vip_source == "plsda") {
    if (is.null(x$res_plsda) || is.null(x$res_plsda$data_vip)) {
      mtFuns$.stop("PLS-DA VIP was not found. Please run PLS-DA first.")
    }
    data_vip <- x$res_plsda$data_vip
  } else if (vip_source == "oplsda") {
    if (is.null(x$res_oplsda) || is.null(x$res_oplsda$data_vip)) {
      mtFuns$.stop("OPLS-DA VIP was not found. Please run OPLS-DA first.")
    }
    data_vip <- x$res_oplsda$data_vip
  } else {
    data_vip <- data.frame(
      mt_feature_id = data_diff$mt_feature_id,
      VIP = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  data_diff$VIP <- data_vip$VIP[
    match(data_diff$mt_feature_id, data_vip$mt_feature_id)
  ]

  if (is.null(padj_cutoff)) {
    idx_sig <- !is.na(data_diff$VIP) &
      data_diff$VIP > vip_cutoff &
      !is.na(data_diff$pvalue) &
      data_diff$pvalue < p_cutoff &
      abs(data_diff$log2FC) >= log2fc_cutoff
  } else {
    idx_sig <- !is.na(data_diff$VIP) &
      data_diff$VIP > vip_cutoff &
      !is.na(data_diff$padj) &
      data_diff$padj < padj_cutoff &
      abs(data_diff$log2FC) >= log2fc_cutoff
  }

  data_diff$change <- "Not significant"
  data_diff$change[idx_sig & data_diff$log2FC > 0] <- "Up"
  data_diff$change[idx_sig & data_diff$log2FC < 0] <- "Down"

  data_diff <- data_diff[order(data_diff$pvalue, -abs(data_diff$log2FC)), , drop = FALSE]

  x$res_diff <- list(
    data_diff = data_diff,
    expr_name = expr_name,
    feature_name = feature_name,
    group_col = group_col,
    case_group = case_group,
    control_group = control_group,
    test_method = test_method,
    p_adjust_method = p_adjust_method,
    vip_source = vip_source,
    vip_cutoff = vip_cutoff,
    p_cutoff = p_cutoff,
    padj_cutoff = padj_cutoff,
    log2fc_cutoff = log2fc_cutoff
  )

  mtFuns$.msg(
    "Finished differential analysis: {sum(data_diff$change == 'Up')} up, {sum(data_diff$change == 'Down')} down.",
    verbose = verbose
  )

  x
}

mtFuns$plot_diff_volcano <- function(data_diff,
  label_col = "feature_name",
  log2fc_col = "log2FC",
  p_col = "pvalue",
  vip_col = "VIP",
  vip_cutoff = 1,
  p_cutoff = 0.05,
  log2fc_cutoff = NULL,
  top_n = 10L,
  top_by = c("pvalue", "VIP", "abs_log2FC"),
  label_significant_only = TRUE,
  show_count = TRUE,
  show_vip_size = TRUE,
  show_threshold = TRUE,
  point_size = 1.8,
  size_range = c(1.2, 4.5),
  alpha = 0.8,
  seed = 2L,
  f_nudge = 0.6,
  title = "Volcano plot")
{
  top_by <- match.arg(top_by)
  set.seed(seed)

  vec_need <- c(label_col, log2fc_col, p_col, vip_col)
  vec_missing <- setdiff(vec_need, names(data_diff))
  if (length(vec_missing) > 0L) {
    mtFuns$.stop("Column `{vec_missing[1L]}` was not found in `data_diff`.")
  }

  data_plot <- data_diff

  data_plot$plot_label <- as.character(data_plot[[ label_col ]])
  data_plot$plot_log2FC <- as.numeric(data_plot[[ log2fc_col ]])
  data_plot$plot_pvalue <- as.numeric(data_plot[[ p_col ]])
  data_plot$plot_VIP <- as.numeric(data_plot[[ vip_col ]])

  data_plot <- data_plot[
    !is.na(data_plot$plot_log2FC) &
      !is.na(data_plot$plot_pvalue) &
      data_plot$plot_pvalue > 0,
    ,
    drop = FALSE
  ]

  data_plot$plot_neg_log10_p <- -log10(data_plot$plot_pvalue)

  idx_sig <- !is.na(data_plot$plot_VIP) &
    data_plot$plot_VIP > vip_cutoff &
    data_plot$plot_pvalue < p_cutoff

  if (!is.null(log2fc_cutoff)) {
    idx_sig <- idx_sig & abs(data_plot$plot_log2FC) >= abs(log2fc_cutoff)
  }

  data_plot$plot_change <- "Not significant"
  data_plot$plot_change[idx_sig & data_plot$plot_log2FC > 0] <- "Up"
  data_plot$plot_change[idx_sig & data_plot$plot_log2FC < 0] <- "Down"

  data_plot$plot_change <- factor(
    data_plot$plot_change,
    levels = c("Down", "Not significant", "Up")
  )

  n_up <- sum(data_plot$plot_change == "Up", na.rm = TRUE)
  n_down <- sum(data_plot$plot_change == "Down", na.rm = TRUE)

  text_criteria <- glue::glue("Criteria: VIP > {vip_cutoff}, P < {p_cutoff}")
  if (!is.null(log2fc_cutoff)) {
    text_criteria <- glue::glue(
      "{text_criteria}, |log2FC| ≥ {abs(log2fc_cutoff)}"
    )
  }

  if (isTRUE(show_count)) {
    title <- glue::glue("{title} (Up: {n_up}; Down: {n_down})")
  }

  if (isTRUE(label_significant_only)) {
    data_lab <- data_plot[data_plot$plot_change %in% c("Up", "Down"), , drop = FALSE]
  } else {
    data_lab <- data_plot
  }

  .slice_top_side <- function(data_x, side, top_n, top_by)
  {
    data_x <- data_x[data_x$plot_change == side, , drop = FALSE]

    if (nrow(data_x) == 0L) {
      return(data_x)
    }

    if (top_by == "pvalue") {
      data_x <- data_x[order(data_x$plot_pvalue, -abs(data_x$plot_log2FC)), , drop = FALSE]
    } else if (top_by == "VIP") {
      data_x <- data_x[order(-data_x$plot_VIP, data_x$plot_pvalue), , drop = FALSE]
    } else {
      data_x <- data_x[order(-abs(data_x$plot_log2FC), data_x$plot_pvalue), , drop = FALSE]
    }

    data_x[seq_len(min(top_n, nrow(data_x))), , drop = FALSE]
  }

  if (nrow(data_lab) > 0L && top_n > 0L) {
    data_lab_up <- .slice_top_side(data_lab, "Up", top_n = top_n, top_by = top_by)
    data_lab_down <- .slice_top_side(data_lab, "Down", top_n = top_n, top_by = top_by)
    data_lab <- rbind(data_lab_up, data_lab_down)
  } else {
    data_lab <- data_lab[0L, , drop = FALSE]
  }

  if (nrow(data_lab) > 0L) {
    vec_abs_fc <- abs(data_lab$plot_log2FC)
    med_fc <- stats::median(vec_abs_fc[is.finite(vec_abs_fc) & vec_abs_fc > 0], na.rm = TRUE)

    if (!is.finite(med_fc) || is.na(med_fc) || med_fc == 0) {
      med_fc <- 1
    }

    data_lab$nudge_x <- med_fc * sign(data_lab$plot_log2FC) * f_nudge
    data_lab$label_side <- ifelse(data_lab$plot_log2FC >= 0, "right", "left")
  }

  if (isTRUE(show_vip_size)) {
    p <- ggplot2::ggplot(
      data_plot,
      ggplot2::aes(
        x = plot_log2FC,
        y = plot_neg_log10_p,
        color = plot_change,
        size = plot_VIP
      )
    ) +
      ggplot2::geom_point(alpha = alpha, stroke = 0)
  } else {
    p <- ggplot2::ggplot(
      data_plot,
      ggplot2::aes(
        x = plot_log2FC,
        y = plot_neg_log10_p,
        color = plot_change
      )
    ) +
      ggplot2::geom_point(alpha = alpha, stroke = 0, size = point_size)
  }

  p <- p +
    ggplot2::scale_color_manual(
      values = c(
        "Down" = "#053061FF",
        "Not significant" = "grey85",
        "Up" = "#67001FFF"
      ),
      drop = FALSE
    ) +
    ggplot2::labs(
      title = title,
      subtitle = text_criteria,
      x = "log2(FC)",
      y = glue::glue("-log10({p_col})"),
      color = "Change",
      size = "VIP"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 9)
    )

  if (isTRUE(show_vip_size)) {
    p <- p + ggplot2::scale_size_continuous(range = size_range)
  }

  if (isTRUE(show_threshold)) {
    p <- p +
      ggplot2::geom_hline(
        yintercept = -log10(p_cutoff),
        linetype = 4L,
        linewidth = 0.5
      ) +
      ggplot2::geom_vline(
        xintercept = 0,
        linetype = 3L,
        linewidth = 0.4
      )

    if (!is.null(log2fc_cutoff) && abs(log2fc_cutoff) > 0) {
      p <- p +
        ggplot2::geom_vline(
          xintercept = c(-abs(log2fc_cutoff), abs(log2fc_cutoff)),
          linetype = 4L,
          linewidth = 0.5
        )
    }
  }

  if (nrow(data_lab) > 0L) {
    p <- p +
      ggrepel::geom_label_repel(
        data = data_lab,
        ggplot2::aes(
          x = plot_log2FC,
          y = plot_neg_log10_p,
          label = plot_label
        ),
        nudge_x = data_lab$nudge_x,
        direction = "y",
        seed = seed,
        size = 3,
        box.padding = 0.35,
        point.padding = 0.2,
        label.padding = grid::unit(0.12, "lines"),
        label.r = grid::unit(0.12, "lines"),
        min.segment.length = 0,
        segment.alpha = 0.6,
        force = 2,
        force_pull = 0.5,
        max.overlaps = Inf,
        show.legend = FALSE
      )
  }

  p
}

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


