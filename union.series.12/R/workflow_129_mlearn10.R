# ==========================================================================
# workflow of mlearn10
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_mlearn10 <- setClass("job_mlearn10", 
  contains = c("job"),
  prototype = prototype(
    pg = "mlearn10",
    info = c(""),
    cite = "",
    method = "",
    tag = "mlearn10",
    analysis = ""
    ))

setMethod_traceable("step0", signature = c(x = "job_mlearn10"),
  function(x){
    step_message("Prepare your data with function `job_mlearn10`.")
  })

setGeneric("asjob_mlearn10",
  function(x, ...) standardGeneric("asjob_mlearn10"))

setMethod_traceable("asjob_mlearn10", signature = c(x = "job_mlearn"),
  function(x, ...)
  {
    .job_mlearn10(copy_job(x))
  })

setMethod_traceable("step1", signature = c(x = "job_mlearn10"),
  function(x, n = 10L, nthread = 1L, validate = NULL, 
    use_rfe = TRUE, nkeep = 10L, debug = FALSE, repeats = 10L,
    exclude = c("XGBoost"), sizes = NULL)
  {
    step_message("Run caret")
    if (!is.null(validate)) {
      if (!is(validate, "job_mlearn10") && !is(validate, "job_mlearn")) {
        stop('Not valid input of validation dataset.')
      }
    }
    x <- methodAdd(x, "采用多种机器学习算法，从候选基因集中筛选关键基因。")
    if (!debug) {
      res_ml_train10 <- x$res_ml_train10 <- run_ml_train10(
        object(x), x$target, x$levels[2], cv_folds = n,
        nthread = nthread, seed = x$seed, use_rfe = use_rfe, 
        repeats = repeats, exclude = exclude, rfe_sizes = sizes
      )
    } else {
      res_ml_train10 <- x$res_ml_train10
    }
    if (use_rfe) {
      x <- methodAdd(x, "在机器学习模型构建前，为减少冗余特征、降低过拟合风险并提升模型稳定性，本研究采用递归特征消除法（recursive feature elimination，RFE）进行变量筛选。RFE 通过迭代训练模型并依据变量贡献度逐步剔除低重要性特征，从而获得最优特征子集。针对不同分类算法，分别调用对应的特征筛选函数进行独立分析，包括随机森林（RF）、逻辑回归（LR）、决策树（DT）、朴素贝叶斯（NB）及线性判别分析（LDA）等模型。最终选择交叉验证表现最佳的特征组合用于后续模型训练与验证。")
      t.rfe_imp <- .get_rfe_importance_ml_train10(res_ml_train10)
      t.rfe_imp <- set_lab_legend(
        t.rfe_imp,
        glue::glue("{x@sig} rfe feature importance"),
        glue::glue("RFE feature 重要性数据")
      )
      x <- tablesAdd(x, t.rfe_imp)
      layout <- wrap_layout(NULL, length(res_ml_train10$models))
      p.rfe_imp <- .plot_importance_ml_train10(t.rfe_imp, ncol = layout$ncol)
      x$rfe_avg_imp <- .get_ggplot_content(p.rfe_imp)
      p.rfe_imp <- set_lab_legend(
        add(layout, p.rfe_imp),
        glue::glue("{x@sig} RFE feature importance"),
        glue::glue("RFE 特征重要性|||每个分面（facet）代表一种独立的建模算法，不同颜色对应不同算法类别。横坐标表示特征名称，纵坐标表示平均特征重要性（Average Importance），数值越高说明该变量在模型判别过程中贡献越大。柱状图长度反映特征的重要性强弱，并按照重要性从高到低排序显示。")
      )
      t.rfe_perf <- .get_rfe_performance_ml_train10(res_ml_train10)
      t.rfe_perf <- set_lab_legend(
        t.rfe_perf,
        glue::glue("{x@sig} rfe performance"),
        glue::glue("RFE 性能数据")
      )
      lst_rfe_perf <- .plot_rfe_performance_ml_train10(
        t.rfe_perf, ncol = layout$ncol
      )
      p.rfe_kappa <- lst_rfe_perf$p.kappa
      p.rfe_kappa <- set_lab_legend(
        add(layout, p.rfe_kappa),
        glue::glue("{x@sig} RFE performance Kappa curve"),
        glue::glue("Kappa 值随特征数量变化曲线|||代表随着输入特征数量变化所对应的交叉验证 Kappa 系数。每个分面代表一种独立建模算法。横坐标表示纳入模型的特征数量（Number of Variables），纵坐标表示交叉验证 Kappa 值（Cross-validated Kappa）。折线及散点表示不同特征子集下模型的平均 Kappa 值，误差线表示标准差（SD）。红色虚线表示该模型达到最佳分类一致性时对应的最优特征数量。Kappa 系数综合考虑随机一致性的影响，可更客观评估模型预测结果与真实分类之间的一致程度。该图可用于辅助判断模型在不同特征维度下的稳健性，尤其适用于类别分布不均衡数据的性能评估。")
      )
      p.rfe_acc <- lst_rfe_perf$p.acc
      p.rfe_acc <- set_lab_legend(
        add(layout, p.rfe_acc),
        glue::glue("{x@sig} RFE performance accuracy curve"),
        glue::glue("准确率随特征数量变化曲线|||各机器学习模型在递归特征消除（recursive feature elimination，RFE）过程中，随着输入特征数量变化所对应的交叉验证准确率（Accuracy）表现。每个分面代表一种独立建模算法。横坐标表示纳入模型的特征数量（Number of Variables），纵坐标表示交叉验证准确率（Cross-validated Accuracy）。折线及散点表示不同特征子集下模型的平均分类准确率，误差线表示标准差（SD），用于反映模型在重复验证中的稳定性。红色虚线表示该模型获得最高准确率时对应的最优特征数量。")
      )
      x <- plotsAdd(x, p.rfe_imp, p.rfe_kappa, p.rfe_acc)
    }
    x <- methodAdd(
      x, "为提高模型评估的稳定性并减少由于样本划分带来的随机性影响，本研究采用重复k折交叉验证（repeated k-fold cross-validation）进行模型性能评估。具体而言，数据集被随机划分为 k 个互斥子集，每次选择其中一个子集作为验证集，其余作为训练集，完成一次 k 折交叉验证。该过程在不同随机划分下重复进行 n 次（repeats, n = {repeats}），最终模型性能取所有重复实验结果的平均值。该方法能够有效降低单次数据划分带来的偏差，尤其适用于小样本数据集的模型评估。"
    )
    text_methods <- .description_ml_train10(
      object(x), cv_folds = n, mlt10 = res_ml_train10
    )
    text_methods <- paste0(seq_along(text_methods), ". ", text_methods)
    x <- methodAdd(
      x, paste0(
        "\n\n\n", bind(text_methods, co = "\n\n\n"), "\n\n\n"
      )
    )
    if (is.null(validate)) {
      # case, no external dataset, thus extract k-fold training validation data from caret.
      res_ml_train10$is_external_data <- FALSE
      x <- methodAdd(x, "对各模型以混淆矩阵、ROC 以及预测性能评估。")
    } else {
      Class <- res_ml_train10$data$Class
      inputs_validate <- .as_input_for_ml_train10(
        object(validate), validate$target, validate$levels[2]
      )
      res_ml_train10$data <- inputs_validate$dat
      res_ml_train10$data$Class <- map_factor(
        res_ml_train10$data$Class, Class
      )
      res_ml_train10$is_external_data <- TRUE
      x$validate_ml_train10 <- res_ml_train10
      x$validate_ml_train10$project <- validate$project
      x <- methodAdd(x, "以 ⟦mark$blue('{validate$project}')⟧ 作为外部验证集，对各模型以混淆矩阵和 ROC 评估。")
    }
    x$res_evaluation <- .evaluation_ml_train10(res_ml_train10)
    x$res_roc <- .roc_ml_train10(res_ml_train10)
    p.rocs <- .plot_roc_ml_train10(x$res_roc)
    p.rocs <- set_lab_legend(
      wrap(p.rocs, 7, 6.5),
      glue::glue("{x@sig} ROC evaluation of models"),
      glue::glue("各模型 ROC|||受试者工作特征曲线（ROC）评价各模型对高负荷组与低负荷组样本的区分能力。横坐标表示假阳性率（False Positive Rate, FPR），纵坐标表示真阳性率（True Positive Rate, TPR）。对角虚线表示随机分类水平（AUC = 0.5）。不同颜色曲线代表不同模型，曲线越接近左上角表示模型分类性能越优，曲线下面积（Area Under the Curve, AUC）越大说明模型判别能力越强。")
    )
    layout <- wrap_layout(NULL, length(res_ml_train10$models))
    p.confusions <- .plot_confusion_ml_train10(
      x$res_evaluation, ncol = layout$ncol
    )
    p.confusions <- set_lab_legend(
      add(layout, p.confusions),
      glue::glue("{x@sig} confusion matrix of models"),
      glue::glue("各模型混淆矩阵|||子图对应一种模型，横坐标表示真实类别（Reference），纵坐标表示预测类别（Prediction）。颜色深浅代表该单元格中的样本数量，颜色越深表示数量越多。矩阵左上角与右下角分别表示正确分类的阴性样本和阳性样本，右上角与左下角表示误分类样本。混淆矩阵可直观反映模型的分类准确率、敏感性及特异性，用于综合比较不同模型的预测性能。")
    )
    layout <- wrap_layout(NULL, 5L)
    p.evaluation <- .plot_evaluation_ml_train10(
      x$res_evaluation, ncol = layout$ncol
    )
    p.evaluation <- set_lab_legend(
      add(layout, p.evaluation),
      glue::glue("{x@sig} prediction comprehensive comparison"),
      glue::glue("机器学习模型分类性能指标综合比较图|||每个分面代表一种评价指标，包括准确率（Accuracy）、Kappa 系数（Kappa）、灵敏度（Sensitivity）、特异度（Specificity）及 F1 值（F1-score）。横坐标为模型名称，纵坐标为对应指标得分，柱状图高度表示模型在该指标上的性能水平，图中数字为具体数值。其中，Accuracy 反映总体预测正确率；Kappa 用于评价模型预测结果与真实分类的一致性，并校正随机因素影响；Sensitivity 表示识别阳性样本的能力；Specificity 表示识别阴性样本的能力；F1-score 为精确率与召回率的综合指标，适用于类别不平衡数据的评估。数值越高通常表示模型表现越优。图中出现 “NA” 表示该指标在当前预测结果下无法计算，提示模型可能存在类别预测失衡或判别能力不足。")
    )
    x <- plotsAdd(x, p.rocs, p.confusions, p.evaluation)
    x$.feature_rfe <- as_feature(
      res_ml_train10$feature_rfe, "RFE 筛选后用于各模型训练的基因"
    )
    return(x)
  })

setMethod_traceable("step2", signature = c(x = "job_mlearn10"),
  function(x, which = 1L)
  {
    step_message("Select best model.")
    x <- snapAdd(x, "为客观筛选综合性能最优的机器学习模型，本研究结合判别能力与分类表现，对各候选模型进行多指标综合评价。首先整合受试者工作特征曲线下面积（area under the curve，AUC）、准确率（Accuracy）、Kappa 系数（Kappa）及 F1-score 等指标。随后构建加权综合评分体系：AUC、Kappa、F1-score 与 Accuracy 分别赋予 0.40、0.25、0.20 和 0.15 的权重，以突出模型判别能力与稳定性的重要性。对于个别模型无法计算的指标（NA 值），按 0 处理，以避免高估其性能。最终依据综合得分对所有模型进行降序排序，选择得分最高者作为最优预测模型，并用于后续诊断效能验证及关键特征解析。")
    lst_best <- .select_best_model_ml_train10(x$res_roc, x$res_evaluation)
    p.score <- .plot_best_model_score_ml_train10(lst_best$summary)
    p.score <- set_lab_legend(
      wrap_scale(p.score, 12, nrow(lst_best$summary), h.size = .15),
      glue::glue("{x@sig} comprehensive score evaluation cross models"),
      glue::glue("各模型综合评分|||综合评分由受试者工作特征曲线下面积（AUC）、准确率（Accuracy）、Kappa 系数（Kappa）及 F1-score 加权计算获得。")
    )
    x <- plotsAdd(x, p.score)
    dat <- lst_best$summary[which, , drop = FALSE]
    fea <- feature(x, "rfe")[[ dat$Model ]]
    x$.feature_best <- as_feature(fea@.Data, "Features of Best Model")
    rfeGenesGest <- dplyr::filter(x$rfe_avg_imp, Algorithm == dat$Model)$Feature
    x <- snapAdd(x, "基于多指标综合评分体系对各机器学习模型进行排序后{aref(p.score)}，⟦mark$red('{dat$Model} 模型获得最高综合得分（Score = {round(dat$Score, 3)}），被确定为最优预测模型')⟧。该模型的受试者工作特征曲线下面积（AUC）为 {round(dat$AUC, 3)}，准确率（Accuracy）为 {round(dat$Accuracy, 3)}，Kappa 系数为 {round(dat$Kappa, 3)}，F1-score 为 {round(dat$F1, 3)}。其 RFE 排名前 {length(rfeGenesGest)} 重要的基因为 {bind(rfeGenesGest)}。该模型对应特征基因为: {bind(fea)}。")
    x$comprehensive_summary <- lst_best$summary
    return(x)
  })

setMethod_traceable("step3", signature = c(x = "job_mlearn10"),
  function(x, use = c("validate", "train"), model = x$comprehensive_summary$Model[1])
  {
    use <- match.arg(use)
    dalex_explain <- .generate_dalex_explain_ml_train10(x$res_ml_train10, model)
    x <- snapAdd(x, "为进一步解释最优机器学习模型的预测机制并识别关键特征基因，在模型构建完成后进行 SHAP（Shapley Additive Explanations）分析。SHAP 基于博弈论思想，通过计算各特征对单一样本预测结果的边际贡献值，量化变量对模型输出的影响方向与作用强度。")
    line <- .get_description("shap.md")
    x <- snapAdd(x, "{line}")
    x <- snapAdd(x, "本研究采用训练完成的最优分类模型，计算各候选基因在全部样本中的 SHAP 值，并进一步汇总平均绝对 SHAP 值评估全局特征重要性，同时结合样本层面的 SHAP 分布展示不同基因在个体预测中的贡献差异。\n\n")
    if (use == "validate" && !is.null(x$validate_ml_train10)) {
      lst <- .as_input_for_dalex_explain_ml_train10(
        x$res_ml_train10, model, data = x$validate_ml_train10$data
      )
      data <- lst$data
      x <- snapAdd(x, "以外部验证集 {x$validate_ml_train10$project} 作为评估数据，以 R 包 `DALEX` ⟦pkgInfo('DALEX')⟧ 使用 SHAP 方法来解释模型的预测结果。")
    } else {
      data <- dalex_explain$data
      x <- snapAdd(x, "以 R 包 `DALEX` ⟦pkgInfo('DALEX')⟧ 使用 SHAP 方法对数据集 {x$project} 来解释模型的预测结果。")
    }
    lst_shap <- .shap_analysis_dalex_explain(dalex_explain, data)
    ps.shap <- .plot_shap_analysis(lst_shap)
    p.importance <- ps.shap$p.importance
    p.importance <- set_lab_legend(
      wrap_scale(p.importance, 12, ncol(data), h.size = .2),
      glue::glue("{x@sig} Global SHAP Feature Importance"),
      glue::glue("全局 SHAP 特征重要性条形图|||各变量平均绝对 SHAP 值（Mean |SHAP value|），反映该特征在全部样本中的总体影响强度。SHAP 绝对值越大，说明该变量对模型预测结果贡献越大，重要性越高")
    )
    p.summary <- ps.shap$p.summary
    p.summary <- set_lab_legend(
      wrap_scale(p.summary, 12, ncol(data), h.size = .2),
      glue::glue("{x@sig} SHAP Summary Distribution"),
      glue::glue(
        "SHAP 分布汇总图|||每个散点代表一个样本在该特征上的 SHAP 值贡献。纵坐标为特征名称，横坐标为 SHAP 值。SHAP 值大于 0 表示该特征推动模型预测向阳性组/疾病组方向，SHAP 值小于 0 表示推动预测向阴性组/正常组方向。绝对值越大说明该特征对预测结果的影响越强。"
      )
    )
    snap <- .stat_shap_analysis(lst_shap, p.importance, p.summary)
    x <- snapAdd(x, "{snap}")
    x <- plotsAdd(x, p.importance, p.summary)
    return(x)
  })


.stat_shap_analysis <- function(x, p1 = NULL, p2 = NULL, top_n = 5L, digits = 3L)
{
  if (is.null(x$summary) || is.null(x$long)) {
    stop("Invalid input object.")
  }

  top_n <- min(top_n, nrow(x$summary))

  imp <- dplyr::slice_head(
    x$summary,
    n = top_n
  )

  vars <- imp$variable_name
  vals <- round(imp$mean_abs_shap, digits)

  rank_txt <- paste(
    glue::glue("{vars}（{vals}）"),
    collapse = "、"
  )

  top1 <- vars[1]
  top1v <- vals[1]

  if (length(vars) >= 2L) {
    remain_txt <- paste(vars[-1], collapse = "、")
  } else {
    remain_txt <- ""
  }

  dat_long <- dplyr::filter(
    x$long,
    variable_name %in% vars
  )

  dir_dat <- dplyr::summarise(
    dplyr::group_by(dat_long, variable_name),
    pos_ratio = mean(contribution > 0, na.rm = TRUE),
    neg_ratio = mean(contribution < 0, na.rm = TRUE),
    span = diff(range(contribution, na.rm = TRUE)),
    .groups = "drop"
  )

  widest <- dir_dat$variable_name[which.max(dir_dat$span)]

  glue::glue(
    "SHAP分析结果显示{aref(p1)}，前{top_n}个重要特征按平均绝对SHAP值排序依次为：{rank_txt}。其中，{top1}的重要性最高（mean |SHAP| = {top1v}），提示其对模型预测结果贡献最大，是最核心的判别特征。其余重要变量包括{remain_txt}，共同参与模型分类决策。

    SHAP分布图进一步显示{aref(p2)}，不同特征在各样本中的贡献方向和强度存在差异。SHAP值大于0表示该特征推动模型预测趋向目标类别，SHAP值小于0表示推动模型预测趋向对照类别。{widest}的SHAP值分布范围最广，提示其在不同个体中具有较强判别能力及一定异质性。

    综合分析表明，上述关键基因在模型识别目标表型过程中发挥重要作用，其中{top1}可能为最具潜力的核心生物标志物。"
  )

}

.plot_shap_analysis <- function(x, n_top = 10L) {

  # ------------------------------------------------------------------
  # Check input object
  # ------------------------------------------------------------------
  if (is.null(x$summary) || is.null(x$long)) {
    stop("Invalid input object.")
  }

  n_top <- min(n_top, nrow(x$summary))

  # ------------------------------------------------------------------
  # Slice top features by mean absolute SHAP
  # ------------------------------------------------------------------
  data_imp <- dplyr::slice_head(x$summary, n = n_top)
  lst_vars <- data_imp$variable_name

  data_long <- dplyr::filter(
    x$long,
    variable_name %in% lst_vars
  )

  data_imp$variable_name <- factor(
    data_imp$variable_name,
    levels = rev(lst_vars)
  )

  data_long$variable_name <- factor(
    data_long$variable_name,
    levels = rev(lst_vars)
  )

  message(glue::glue("Plotting top {length(lst_vars)} SHAP features"))

  # ------------------------------------------------------------------
  # Global feature importance bar plot
  # ------------------------------------------------------------------
  p_importance <- ggplot2::ggplot(
    data_imp,
    ggplot2::aes(
      x = variable_name,
      y = mean_abs_shap
    )
    ) +
  ggplot2::geom_col(width = 0.75) +
  ggplot2::geom_text(
    ggplot2::aes(label = round(mean_abs_shap, 3)),
    hjust = -0.1,
    size = 3.5
    ) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    title = "Global SHAP Feature Importance",
    x = "Feature",
    y = "Mean |SHAP value|"
    ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),
    panel.grid.minor = ggplot2::element_blank()
    ) +
  ggplot2::expand_limits(
    y = max(data_imp$mean_abs_shap, na.rm = TRUE) * 1.08
  )

  # ------------------------------------------------------------------
  # SHAP summary distribution (Beeswarm / Violin-like)
  # ------------------------------------------------------------------
  if (!"feature_value" %in% colnames(data_long)) {
    stop("Column 'feature_value' is required in x$long for Beeswarm coloring.")
  }

  p_summary <- ggplot2::ggplot(
    data_long,
    ggplot2::aes(
      y = contribution,
      x = variable_name,
      color = feature_value
    )
    ) +
  ggbeeswarm::geom_quasirandom(
    ggplot2::aes(y = contribution),
    width = 0.2,
    alpha = 0.7,
    size = 1.8,
    groupOnX = FALSE
    ) +
  ggplot2::scale_color_gradient(
    low = "blue",
    high = "red",
    name = "Feature value"
    ) +
  ggplot2::labs(
    title = "SHAP Summary Distribution (Beeswarm)",
    x = "Feature",
    y = "SHAP value"
    ) +
  ggplot2::theme_minimal() +
  coord_flip() +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),
    panel.grid.minor = ggplot2::element_blank()
  )

  # ------------------------------------------------------------------
  # Return results
  # ------------------------------------------------------------------
  list(
    p.importance = p_importance,
    p.summary = p_summary
  )
}

.shap_analysis_dalex_explain <- function(exp, data, B = 100L)
{
  # Check input
  if (missing(data) || is.null(data)) {
    stop("`data` is required.")
  }

  data <- as.data.frame(data)

  # ------------------------------------------------------------
  # Run local SHAP for each sample
  # ------------------------------------------------------------
  message("DALEX::predict_parts (type = 'shap')")

  fun_shap <- function(...) {
    lapply(seq_len(nrow(data)),
      function(i) {
        DALEX::predict_parts(
          explainer = exp,
          new_observation = data[i, , drop = FALSE],
          type = "shap",
          B = B
        )
      }
    )
  }

  id_args <- list(
    if (inherits(exp$model, "train")) exp$model$finalModel else exp$model,
    exp$data,
    data,
    B
  )

  lst_res <- expect_local_data(
    "tmp",
    "dalex_shap",
    fun_shap,
    id_args
  )

  message("Finished SHAP calculation.")

  # ------------------------------------------------------------
  # Bind results
  # ------------------------------------------------------------
  data_long <- dplyr::bind_rows(
    lapply(lst_res, function(x) tibble::as_tibble(x)),
    .id = "id"
  )

  # Keep only feature rows
  data_long <- dplyr::filter(
    data_long,
    variable_name != "",
    variable_name != "_baseline_"
  )

  # Attach original feature value
  # For each row, look up the value in original data
  data_long$feature_value <- sapply(seq_len(nrow(data_long)),
    function(i) {
      row_idx <- as.integer(data_long$id[i])
      feat_name <- data_long$variable_name[i]
      data[row_idx, feat_name]
    })

  # ------------------------------------------------------------
  # Global summary from local SHAP
  # ------------------------------------------------------------
  sum_dat <- dplyr::summarise(
    dplyr::group_by(data_long, variable_name),
    mean_abs_shap = mean(abs(contribution), na.rm = TRUE),
    mean_shap = mean(contribution, na.rm = TRUE),
    sd_shap = stats::sd(contribution, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  )

  sum_dat <- dplyr::arrange(
    sum_dat,
    dplyr::desc(mean_abs_shap)
  )

  message("Finished SHAP aggregation.")

  list(
    long = data_long,
    summary = sum_dat
  )
}

.as_input_for_dalex_explain_ml_train10 <- function(mlt10, nm, data = NULL)
{
  # Validate input
  if (!is(mlt10, "ml_train10")) {
    stop('!is(mlt10, "ml_train10").')
  }

  if (is.null(mlt10$models[[nm]])) {
    stop("Model not found.")
  }

  # Use internal data by default
  if (is.null(data)) {
    dat <- mlt10$data
  } else {
    dat <- data
  }

  lv <- levels(mlt10$data$Class)
  vars <- mlt10$feature_rfe[[nm]]

  # Keep selected variables only
  if (!is.null(vars)) {
    keep <- intersect(colnames(dat), c("Class", vars))
    dat <- dat[, keep, drop = FALSE]
  }

  # Prepare x and y for DALEX::explain
  if ("Class" %in% colnames(dat)) {
    x <- dplyr::select(dat, -Class)
    y <- as.integer(dat$Class == lv[1])
  } else {
    x <- dat
    y <- NULL
  }
  list(data = x, y = y, levels = lv)
}

.generate_dalex_explain_ml_train10 <- function(mlt10, nm)
{
  obj <- .as_input_for_dalex_explain_ml_train10(mlt10 = mlt10, nm = nm)
  fit <- mlt10$models[[nm]]
  lv <- obj$levels
  fun_pred <- function(model, newdata)
  {
    if (inherits(model, "xgb.Booster")) {
      pred <- predict(
        model, newdata = as.matrix(newdata)
      )
      return(as.numeric(pred))
    }
    pred <- .caret_predict(
      fit = model, dat = newdata, type = "prob"
    )
    if (is.data.frame(pred) && lv[1] %in% colnames(pred)) {
      return(as.numeric(pred[, lv[1]]))
    }
    rep(NA_real_, nrow(newdata))
  }
  DALEX::explain(model = fit, data = obj$data,
    y = obj$y, predict_function = fun_pred, label = nm, verbose = TRUE
  )
}

.plot_best_model_score_ml_train10 <- function(dat)
{
  dat$Model <- factor(
    dat$Model,
    levels = rev(dat$Model)
  )
  ggplot(dat, aes(x = Model, y = Score)) +
    geom_col(width = 0.75) +
    geom_text(aes(label = round(Score, 3)), hjust = -0.1, size = 3.5) +
    coord_flip() +
    labs(
      title = "Comprehensive Model Score Ranking",
      x = "Model",
      y = "Score"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.minor = element_blank()
    ) +
    expand_limits(y = max(dat$Score, na.rm = TRUE) * 1.08)
}

.select_best_model_ml_train10 <- function(res_roc, res_evaluation)
{
  if (!is(res_roc, "ml_roc")) {
    stop('!is(res_roc, "ml_roc").')
  }
  if (!is(res_evaluation, "ml_evaluation")) {
    stop('!is(res_evaluation, "ml_evaluation").')
  }
  roc_df <- res_roc$summary
  eval_df <- res_evaluation$summary

  dat <- dplyr::left_join(
    dplyr::select(roc_df, Model, AUC),
    dplyr::select(
      eval_df,
      Model, Accuracy, Kappa, F1
    ),
    by = "Model"
  )

  dat <- dplyr::mutate(
    tibble::as_tibble(dat),
    dplyr::across(
      dplyr::where(is.numeric), function(x) dplyr::coalesce(x, 0)
    ),
    # Weighted comprehensive score
    Score =
      0.40 * AUC +
      0.25 * Kappa +
      0.20 * F1 +
      0.15 * Accuracy
  )

  dat <- dplyr::arrange(
    dat,
    dplyr::desc(Score),
    dplyr::desc(AUC),
    dplyr::desc(Kappa)
  )

  dat$Rank <- seq_len(nrow(dat))

  dat <- dplyr::select(
    dat,
    Rank, Model,
    AUC, Accuracy, Kappa, F1,
    Score
  )
  list(
    best_model = dat$Model[1],
    summary = dat
  )
}


map_factor <- function(x, template) {
  factor(
    levels(template)[as.integer(x)],
    levels = levels(template)
  )
}

setMethod_traceable("clear", signature = c(x = "job_mlearn10"),
  function(x, ..., name = substitute(x, parent.frame(1)))
  {
    eval(name)
    x <- callNextMethod(
      x, ..., name = name,
      expr_lite = expression({
        x$res_ml_train10$models <- NULL
      })
    )
    return(x)
  })

.get_rfe_importance_ml_train10 <- function(mlt10) {
  data <- lapply(mlt10$rfes, function(x) x$fit_rfe$variables)
  data <- dplyr::bind_rows(data, .id = "model")
  dplyr::select(
    data, Algorithm = model,
    Feature = var,
    Importance = Overall,
    Resample
  )
}

.get_rfe_performance_ml_train10 <- function(mlt10) {
  data <- lapply(mlt10$rfes, function(x) x$fit_rfe$results)
  data <- dplyr::bind_rows(data, .id = "model")
}

.plot_rfe_performance_ml_train10 <- function(data, ...) {
  opt_df <- dplyr::slice_max(
    dplyr::group_by(data, model),
    order_by = Accuracy,
    n = 1L,
    with_ties = FALSE
  )

  # Accuracy plot
  p.acc <- ggplot(data, aes(x = Variables, y = Accuracy)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_errorbar(
      aes(ymin = Accuracy - AccuracySD, ymax = Accuracy + AccuracySD),
      width = 0.2) +
    geom_vline(
      data = opt_df,
      aes(xintercept = Variables),
      color = "red",
      linetype = "dashed",
      linewidth = 0.8
      ) +
    facet_wrap(~ model, scales = "free_x", ...) +
    scale_x_continuous(
      breaks = scales::pretty_breaks(n = 5L)
      ) +
    scale_y_continuous(
      labels = scales::percent
      ) +
    labs(
      title = "Accuracy vs Feature Count",
      x = "Number of Variables",
      y = "Cross-validated Accuracy"
      ) +
    theme_minimal() +
    theme(legend.position = "none")

  # Kappa plot
  p.kappa <- ggplot(data, aes(x = Variables, y = Kappa)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = Kappa - KappaSD, ymax = Kappa + KappaSD),
      width = 0.2) +
    geom_vline(
      data = opt_df,
      aes(xintercept = Variables),
      color = "red",
      linetype = "dashed",
      linewidth = 0.8
      ) +
    facet_wrap(~ model, scales = "free_x", ...) +
    scale_x_continuous(
      breaks = scales::pretty_breaks(n = 5L)
      ) +
    labs(
      title = "Kappa vs Feature Count",
      x = "Number of Variables",
      y = "Cross-validated Kappa"
      ) +
    theme_minimal() +
    theme(legend.position = "none")
  namel(p.kappa, p.acc)
}

.plot_importance_ml_train10 <- function(data, ...) {
  avg_imp <- dplyr::summarise(
    dplyr::group_by(data, Algorithm, Feature),
    Avg_Importance = mean(Importance, na.rm = TRUE),
    .groups = "drop"
  )
  # scale the importance
  avg_imp <- dplyr::mutate(
    dplyr::group_by(avg_imp, Algorithm),
    Avg_Importance = (
      Avg_Importance - min(Avg_Importance)
      ) / (
      max(Avg_Importance) - min(Avg_Importance)
    )
  )
  avg_imp <- dplyr::group_modify(
    dplyr::group_by(avg_imp, Algorithm),
    ~ dplyr::slice_max(.x, order_by = Avg_Importance, n = 10L, with_ties = FALSE)
  )
  avg_imp <- dplyr::ungroup(avg_imp)
  ggplot(avg_imp,
    aes(x = tidytext::reorder_within(Feature, Avg_Importance, Algorithm),
      y = Avg_Importance)
    ) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~ Algorithm, scales = "free_y", ...) +
  tidytext::scale_x_reordered() +
  labs(title = "Algorithm Feature Importance Ranking",
    x = "Feature", y = "Average Importance Scaled") +
  theme_minimal()
}

.as_input_for_ml_train10 <- function(data, target, positive_class)
{
  data <- as.data.frame(data)
  target <- as.factor(target)
  if (nrow(data) != length(target)) {
    stop("nrow(data) must equal length(target)")
  }
  if (length(levels(target)) != 2) {
    stop("target must be binary")
  }
  lv <- levels(target)
  if (is.null(positive_class)) {
    positive_class <- lv[2]
    message(glue::glue("Use '{lv[2]}' as case, is sure that?"))
  }
  if (!positive_class %in% lv) {
    stop("positive_class not found in target")
  }
  target <- factor(target, levels = c(positive_class, setdiff(lv, positive_class)))
  dat <- data.frame(Class = target, data)
  colnames(dat) <- make.names(colnames(dat), unique = TRUE)
  namel(data, dat, target, positive_class)
}

run_ml_train10 <- function(
  data, target, positive_class = NULL,
  seed = 123L, cv_folds = 10L,
  tune_length = 5L, nthread = 1L,
  use_rfe = TRUE, repeats = 10L,
  rfe_sizes = NULL, exclude = NULL
)
{
  # ------------------------------------------------------------------
  # Check input
  # ------------------------------------------------------------------
  inputs <- .as_input_for_ml_train10(
    data, target, positive_class
  )

  data <- inputs$data
  dat <- inputs$dat
  target <- inputs$target
  positive_class <- inputs$positive_class

  # ------------------------------------------------------------------
  # Model configuration
  # ------------------------------------------------------------------
  cfg <- make_cfg(data, tune_length)
  if (!is.null(exclude)) {
    cfg <- cfg[ !names(cfg) %in% exclude ]
  }

  ctrl <- caret::trainControl(
    method = "repeatedcv",
    number = cv_folds,
    classProbs = TRUE,
    summaryFunction = caret::twoClassSummary,
    savePredictions = "final",
    allowParallel = TRUE,
    repeats = repeats,
    seeds = make_seeds(
      n_folds = cv_folds,
      repeats = repeats,
      n_candidates = 50L,
      seed = seed
    )
  )

  fun_caret <- function(...) {
    fit <- tryCatch(
      do.call(caret::train, append(base_args, args)),
      error = function(e) {
        print(e)
        NULL
      }
    )
    if (is.null(fit)) {
      stop('is.null(fit).')
    }
    return(fit)
  }

  # ------------------------------------------------------------------
  # Train caret models
  # ------------------------------------------------------------------
  model_list <- list()
  rfe_list <- list()
  feature_rfe <- list()

  for (nm in names(cfg)) {

    args <- cfg[[nm]]
    dat_use <- dat

    if (use_rfe) {
      message(glue::glue("RFE selecting for {nm} ..."))
      rfe_list[[ nm ]] <- perform_rfe(
        dat, nm, cv_folds, seed = seed, rfe_sizes = rfe_sizes
      )
      vars <- rfe_list[[ nm ]]$vars
      feature_rfe[[ nm ]] <- vars
      dat_use <- dat[, c("Class", vars), drop = FALSE]
    }

    base_args <- list(
      form = Class ~ .,
      data = dat_use,
      metric = "ROC",
      trControl = ctrl
    )

    message(glue::glue("Training {nm} ..."))

    fit <- expect_local_data(
      "tmp",
      glue::glue("ml_{nm}"),
      fun_caret,
      list(colnames(dat_use), rownames(dat_use), levels(target), positive_class,
        seed, cv_folds, tune_length, use_rfe, repeats)
    )

    if (is.null(fit)) {
      message(glue::glue("Failed in training {nm}"))
    }

    model_list[[nm]] <- fit
  }

  # ------------------------------------------------------------------
  # XGBoost (manual training)
  # ------------------------------------------------------------------
  xgb_fun <- function(...) {

    params <- list(
      objective = "binary:logistic",
      eval_metric = "auc",
      eta = 0.05,
      max_depth = 3,
      subsample = 0.8,
      colsample_bytree = 0.8
    )

    cv <- tryCatch(
      xgboost::xgb.cv(
        params = params,
        data = dtrain,
        nrounds = 5,
        nfold = cv_folds,
        early_stopping_rounds = 3,
        verbose = 0,
        prediction = FALSE
      ),
      error = function(e) NULL
    )

    if (is.null(cv)) {
      return(NULL)
    }

    eval_log <- cv$evaluation_log

    best_nrounds <- which.max(
      eval_log$test_auc_mean
    )

    model <- tryCatch(
      xgboost::xgb.train(
        params = params,
        data = dtrain,
        nrounds = best_nrounds,
        verbose = 0
      ),
      error = function(e) NULL
    )

    attr(model, "cv_auc") <- max(
      eval_log$test_auc_mean
    )
    attr(model, "best_nrounds") <- best_nrounds
    return(model)
  }

  if (is.null(exclude) || !any(exclude == "XGBoost")) {
    dat_xgb <- dat
    if (use_rfe) {
      message("RFE selecting for XGBoost ...")
      rfe_list[[ "XGBoost" ]] <- perform_rfe(
        dat, "XGBoost", cv_folds, seed = seed, rfe_sizes = rfe_sizes
      )
      vars <- rfe_list[[ "XGBoost" ]]$vars
      feature_rfe[[ "XGBoost" ]] <- vars
      dat_xgb <- dat[, c("Class", vars), drop = FALSE]
    }
    x_mat <- as.matrix(dat_xgb[, -1, drop = FALSE])
    y_vec <- ifelse(dat_xgb$Class == levels(dat_xgb$Class)[1], 1L, 0L)
    dtrain <- xgboost::xgb.DMatrix(data = x_mat, label = y_vec)
    message("Training XGBoost ...")
    model_list$XGBoost <- expect_local_data(
      "tmp", "ml_XGBoost", xgb_fun,
      list(
        colnames(x_mat), rownames(x_mat), levels(target), positive_class,
        seed, cv_folds, tune_length, use_rfe
      )
    )
  }

  # ------------------------------------------------------------------
  # Summary table
  # ------------------------------------------------------------------
  res <- data.frame(
    Model = names(model_list),
    CV_AUC = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(model_list)) {
    fit <- model_list[[i]]
    if (is.null(fit)) {
      next
    }
    if (inherits(fit, "train")) {
      if ("ROC" %in% colnames(fit$results)) {
        res$CV_AUC[i] <- max(
          fit$results$ROC,
          na.rm = TRUE
        )
      }
    } else if (inherits(fit, "xgb.Booster")) {
      res$CV_AUC[i] <- attr(fit, "cv_auc")
    }
  }

  # ------------------------------------------------------------------
  # Feature genes
  # ------------------------------------------------------------------
  feat_res <- list()
  feat_score <- list()

  for (nm in names(model_list)) {
    fit <- model_list[[nm]]
    if (is.null(fit)) {
      next
    }
    if (inherits(fit, "train")) {
      vi <- tryCatch(
        caret::varImp(fit)$importance,
        error = function(e) NULL
      )
      if (!is.null(vi)) {
        score <- vi[, 1]
        names(score) <- rownames(vi)
        score <- sort(
          score,
          decreasing = TRUE
        )
        feat_res[[nm]] <- names(score)
        feat_score[[nm]] <- score
      }
    }
  }

  if (!is.null(model_list$XGBoost)) {
    imp <- tryCatch(
      xgboost::xgb.importance(
        model = model_list$XGBoost
      ),
      error = function(e) NULL
    )
    if (!is.null(imp)) {
      feat_res$XGBoost <- imp$Feature
    }
  }

  # ------------------------------------------------------------------
  # Order result
  # ------------------------------------------------------------------
  res <- res[order(-res$CV_AUC), ]
  rownames(res) <- NULL

  # ------------------------------------------------------------------
  # Return
  # ------------------------------------------------------------------
  out <- list(
    models = model_list, rfes = rfe_list, summary = res,
    feature_genes = feat_res, feature_score = feat_score,
    feature_rfe = feature_rfe, data = dat,
    positive_class = positive_class,
    levels = levels(target), use_rfe = use_rfe
  )
  class(out) <- "ml_train10"
  return(out)
}

make_seeds <- function(
  n_folds, n_candidates = 50L, seed = 12345L, repeats = 1L
)
{
  set.seed(seed)
  n_folds <- n_folds * repeats
  seeds <- vector(
    mode = "list",
    length = n_folds + 1L
  )
  for (i in seq_len(n_folds)) {
    seeds[[i]] <- sample.int(1000000L, n_candidates)
  }
  seeds[[n_folds + 1]] <- sample.int(1000000L, 1L)
  return(seeds)
}

perform_rfe <- function(dat, alg_name, cv_folds, rfe_sizes = NULL, seed = 12345L)
{
  fun_rfe <- function(...) {
    message(glue::glue("Run RFE: {alg_name}"))
    if (ncol(dat) <= 2) {
      return(colnames(dat)[colnames(dat) != "Class"])
    }
    rfe_funcs <- switch(
      alg_name,
      "RF"      = caret::rfFuncs,
      "LR"      = caret::lrFuncs,
      "DT"      = caret::treebagFuncs,
      "NB"      = caret::nbFuncs,
      "LDA"     = caret::ldaFuncs,
      "BT"      = caret::rfFuncs,
      caret::caretFuncs
    )
    p <- ncol(dat) - 1L
    sizes <- rfe_sizes
    if (is.null(sizes)) {
      sizes <- unique(
        round(seq(
            from = min(5L, p),
            to   = min(p, 10L),
            length.out = min(10L, p)
            ))
      )
    }
    rfe_funcs$selectSize <- function(x, metric, maximize) {
      x <- x[x$Variables %in% sizes, , drop = FALSE]
      if (!nrow(x)) {
        stop("No valid subset sizes found in resampling results.")
      }
      if (maximize) {
        best <- which.max(x[[metric]])
      } else {
        best <- which.min(x[[metric]])
      }
      x$Variables[best]
    }
    sizes <- sizes[sizes <= p]
    seeds <- make_seeds(n_folds = cv_folds, n_candidates = 50L, seed = seed)
    ctrl_rfe <- caret::rfeControl(
      functions = rfe_funcs,
      method = "cv",
      number = cv_folds,
      verbose = FALSE,
      returnResamp = "final",
      seeds = seeds
    )
    x <- dat[, setdiff(colnames(dat), "Class"), drop = FALSE]
    y <- dat$Class
    fit_rfe <- tryCatch(
      caret::rfe(x = x, y = y, sizes = sizes, rfeControl = ctrl_rfe),
      error = function(e) {
        print(e)
        NULL
      }
    )
    if (is.null(fit_rfe)) {
      stop(glue::glue("Failed to perform RFE: {alg_name}"))
      # return(colnames(x))
    }
    vars <- predictors(fit_rfe)
    if (length(vars) == 0) {
      vars <- colnames(x)
    }
    list(fit_rfe = fit_rfe, vars = vars)
  }
  res <- expect_local_data(
    "tmp", glue::glue("rfe_{alg_name}"), fun_rfe,
    list(
      colnames(dat), rownames(dat),
      dat$Class, cv_folds, seed,
      rfe_sizes
    )
  )
  return(res)
}

# ==========================================================
# Dynamic method text generator for ml_train10
# Return: named character vector / list
# Each model = one independent paragraph
# ==========================================================

.description_ml_train10 <- function(data, cv_folds = 10L, 
  tune_length = 5L, mlt10 = NULL)
{
  cfg <- make_cfg(data, tune_length)

  n_sample <- nrow(data)
  p_feature <- ncol(data)

  txt <- list(

    LR = glue::glue(
      "采用逻辑回归（Logistic Regression，LR）模型进行分类分析，基于 R 包 `stats` ⟦pkgInfo('stats')⟧ 中广义线性模型函数（caret method = glm）进行拟合。以二分类结局变量为因变量，纳入{n_sample}例样本及{p_feature}个特征变量，并采用{cv_folds}折交叉验证评估模型稳定性与预测性能。"
    ),

    DT = glue::glue(
      "采用决策树（Decision Tree，DT）模型进行分类分析，基于 R 包 `rpart` ⟦pkgInfo('rpart')⟧（caret method = rpart）构建树模型，并通过复杂度参数（cp）进行剪枝优化。采用{cv_folds}折交叉验证筛选最优参数组合。"
    ),

    SVM = glue::glue(
      "采用支持向量机（Support Vector Machine，SVM）模型进行分类分析，基于 R 包 `kernlab` ⟦pkgInfo('kernlab')⟧ 的线性支持向量机算法（caret method = svmLinear）构建分类模型。采用{cv_folds}折交叉验证确定最优参数。"
    ),

    RF = glue::glue(
      "采用随机森林（Random Forest，RF）模型进行分类分析，基于 R 包 `randomForest` ⟦pkgInfo('randomForest')⟧（caret method = rf）构建{cfg$RF$ntree}棵决策树进行集成学习。采用{cv_folds}折交叉验证优化模型参数。"
    ),

    KNN = glue::glue(
      "采用 K 近邻（K-Nearest Neighbor，KNN）模型进行分类分析，基于 R 包 `class` ⟦pkgInfo('class')⟧（caret method = knn）根据邻近样本投票完成分类。采用{cv_folds}折交叉验证筛选最优参数。"
    ),

    BT = glue::glue(
      "采用 Bagging Tree（BT）模型进行分类分析，基于 R 包 `randomForest` ⟦pkgInfo('randomForest')⟧（caret method = rf）实现装袋集成学习。构建{cfg$BT$ntree}棵树，并设定 mtry = {cfg$BT$tuneGrid$mtry}，最终通过多数投票获得预测结果。采用{cv_folds}折交叉验证优化模型参数。"
    ),

    LDA = glue::glue(
      "采用线性判别分析（Linear Discriminant Analysis，LDA）模型进行分类分析，基于 R 包 `MASS` ⟦pkgInfo('MASS')⟧（caret method = lda）构建线性判别函数，并采用{cv_folds}折交叉验证评估分类性能。"
    ),

    NNET = glue::glue(
      "采用人工神经网络（Neural Network，NNET）模型进行分类分析，基于 R 包 `nnet` ⟦pkgInfo('nnet')⟧（caret method = nnet）构建单隐层前馈神经网络。采用{cv_folds}折交叉验证筛选最优结构。"
    ),

    NB = glue::glue(
      "采用朴素贝叶斯（Naive Bayes，NB）模型进行分类分析，基于 R 包 `klaR` ⟦pkgInfo('klaR')⟧（caret method = nb）估计各特征的后验概率，并采用{cv_folds}折交叉验证评估模型性能。参数搜索强度设定为 tuneLength = {tune_length}。"
    ),

    EN = glue::glue(
      "采用弹性网络（Elastic Net，EN）模型进行分类分析，基于 R 包 `glmnet` ⟦pkgInfo('glmnet')⟧（caret method = glmnet）进行建模。采用{cv_folds}折交叉验证确定最优参数组合。"
    ),

    Lasso = glue::glue(
      "采用 LASSO 回归模型进行分类分析，基于 R 包 `glmnet` ⟦pkgInfo('glmnet')⟧（caret method = glmnet）进行建模，并通过 L1 正则化筛选关键变量。采用{cv_folds}折交叉验证确定最优参数。"
    ),

    RR = glue::glue(
      "采用岭回归（Ridge Regression，RR）模型进行分类分析，基于 R 包 `glmnet` ⟦pkgInfo('glmnet')⟧（caret method = glmnet）进行建模，并通过 L2 正则化提高模型稳定性。采用{cv_folds}折交叉验证筛选最优模型。"
    ),

    GBM = glue::glue(
      "采用梯度提升机（Gradient Boosting Machine，GBM）模型进行分类分析，基于 R 包 `gbm` ⟦pkgInfo('gbm')⟧（caret method = gbm）逐步迭代构建弱学习器。参数搜索范围包括树深度（{paste(unique(cfg$GBM$tuneGrid$interaction.depth), collapse = '、')}）、树数量（{paste(unique(cfg$GBM$tuneGrid$n.trees), collapse = '、')}）、学习率（{paste(unique(cfg$GBM$tuneGrid$shrinkage), collapse = '、')}）及最小叶节点样本数（{paste(unique(cfg$GBM$tuneGrid$n.minobsinnode), collapse = '、')}），采用{cv_folds}折交叉验证优化模型。"
    ),

    XGBoost = glue::glue(
      "采用极端梯度提升算法（Extreme Gradient Boosting，XGBoost）进行分类分析，基于 R 包 `xgboost` ⟦pkgInfo('xgboost')⟧ 构建模型。设置学习率 eta = 0.05、最大树深 max_depth = 3、行采样比例 subsample = 0.8、列采样比例 colsample_bytree = 0.8，并采用{cv_folds}折交叉验证与早停策略确定最优迭代轮数。"
    )
  )
  if (!is.null(mlt10)) {
    use <- names(mlt10$models)
    use <- use[ !vapply(mlt10, is.null, logical(1)) ]
    txt <- txt[ use ]
  }
  txt
}

.caret_predict <- function(fit, dat, type, ...) {
  if (!is.null(fit$coefnames)) {
    dat <- dat[, colnames(dat) %in% fit$coefnames ]
  }
  predict(fit, newdata = dat, type = type) 
}

.predict_or_extract_ml_train10 <- function(mlt10, type = c("raw", "prob"), cutoff = 0.5,
  predict = !is.null(mlt10$is_external_data) && mlt10$is_external_data
)
{
  type <- match.arg(type)

  if (!is(mlt10, "ml_train10")) {
    stop('!is(mlt10, "ml_train10").')
  }
  if (is.null(mlt10$models)) {
    stop("Invalid input object")
  }

  dat <- mlt10$data
  truth <- dat$Class
  lv <- levels(truth)

  pred_list <- sapply(names(mlt10$models), simplify = FALSE,
    function(nm) {
      fit <- mlt10$models[[nm]]
      if (is.null(fit)) {
        return(NULL)
      }
      vars <- mlt10$feature_rfe[[nm]]
      dat_use <- if (is.null(vars)) {
        dat[, colnames(dat) != "Class" ]
      } else {
        dat[, colnames(dat) %in% vars, drop = FALSE]
      }
      pred <- NULL
      if (inherits(fit, "train")) {
        if (predict) {
          ## Predict by external/new data
          pred <- tryCatch(
            .caret_predict(fit = fit, dat = dat_use, type = type),
            error = function(e) NULL
          )
          if (!is.null(pred) && type == "raw") {
            pred <- factor(pred, levels = lv)
          }
          if (!is.null(pred) && type == "prob") {
            if (is.data.frame(pred) && lv[1L] %in% colnames(pred)) {
              pred <- pred[, lv[1L]]
            } else {
              pred <- NULL
            }
          }
        } else {

          ## Extract repeatedcv out-of-fold predictions
          pred_df <- fit$pred
          if (is.null(pred_df)) {
            rlang::abort(
              'is.null(pred_df), try extract prediction results failed.'
            )
            return(NULL)
          }
          ## Keep only best tuning parameters
          if (!is.null(fit$bestTune)) {
            for (kk in colnames(fit$bestTune)) {
              pred_df <- pred_df[
                pred_df[[kk]] == fit$bestTune[[kk]],
                ,
                drop = FALSE
              ]
            }
          }
          ## Average repeated predictions by rowIndex
          if ("rowIndex" %in% colnames(pred_df)) {
            idx_split <- split(seq_len(nrow(pred_df)), pred_df$rowIndex)
            pred_df <- do.call(rbind,
              lapply(idx_split,
                function(ii) {
                  z <- pred_df[ii, , drop = FALSE]
                  out <- z[1L, , drop = FALSE]
                  out$obs <- z$obs[1L]

                  if (lv[1L] %in% colnames(z)) {
                    out[[lv[1L]]] <- mean(z[[lv[1L]]], na.rm = TRUE)
                  }

                  if (lv[2L] %in% colnames(z)) {
                    out[[lv[2L]]] <- mean(z[[lv[2L]]], na.rm = TRUE)
                  }

                  if (lv[1L] %in% colnames(out)) {
                    out$pred <- ifelse(out[[lv[1L]]] >= cutoff, lv[1L], lv[2L])
                  }
                  out$pred <- factor(out$pred, levels = lv)
                  out$Resample <- NA_character_
                  out
                })
            )
            rownames(pred_df) <- NULL
          }

          if (type == "raw") {
            pred <- factor(pred_df$pred, levels = lv)
          }

          if (type == "prob") {
            if (lv[1L] %in% colnames(pred_df)) {
              pred <- pred_df[[lv[1L]]]
            }
          }

          truth_sub <- factor(pred_df$obs, levels = lv)
          attr(pred, "truth") <- truth_sub
        }

      } else if (inherits(fit, "xgb.Booster")) {

        pred <- tryCatch(
          predict(fit, newdata = as.matrix(dat_use[, -1L, drop = FALSE])),
          error = function(e) NULL
        )

        if (!is.null(pred) && type == "raw") {
          pred <- ifelse(pred >= cutoff, lv[1L], lv[2L])
          pred <- factor(pred, levels = lv)
        }
      }

      pred
    })

  pred_list <- pred_list[!vapply(pred_list, is.null, logical(1L))]

  truth_out <- if (predict) truth else NULL

  if (!predict) {
    truth_out <- lapply(pred_list, attr, which = "truth")
  }
  if (all(duplicated(truth_out)[-1])) {
    truth_out <- truth_out[[ 1 ]]
  }

  list(truth = truth_out, levels = lv, pred = pred_list)
}

# ===============================================================
# Confusion Matrix Evaluation + Plot
# Support:
#   caret::train
#   xgboost::xgb.Booster
# Input:
#   result from ml_train10()
# ===============================================================

.evaluation_ml_train10 <- function(mlt10, cutoff = 0.5)
{
  res_pred <- .predict_or_extract_ml_train10(
    mlt10 = mlt10,
    type = "raw",
    cutoff = cutoff
  )
  truth <- res_pred$truth
  lv <- res_pred$levels
  cm_list <- list()
  stat_list <- list()
  for (nm in names(res_pred$pred)) {
    pred <- res_pred$pred[[nm]]
    cm <- caret::confusionMatrix(
      data = pred,
      reference = truth,
      positive = lv[1]
    )
    cm_list[[nm]] <- cm
    stat_list[[nm]] <- data.frame(
      Model = nm,
      Accuracy = unname(cm$overall["Accuracy"]),
      Kappa = unname(cm$overall["Kappa"]),
      Sensitivity = unname(cm$byClass["Sensitivity"]),
      Specificity = unname(cm$byClass["Specificity"]),
      F1 = unname(cm$byClass["F1"]),
      stringsAsFactors = FALSE
    )
  }
  stat_df <- do.call(rbind, stat_list)
  rownames(stat_df) <- NULL
  stat_df <- stat_df[order(-stat_df$Accuracy), ]
  out <- list(
    matrix = cm_list,
    summary = stat_df
  )
  class(out) <- "ml_evaluation"
  out
}

# ===============================================================
# ROC Evaluation + Plot
# Support:
#   caret::train
#   xgboost::xgb.Booster
# ===============================================================

.roc_ml_train10 <- function(mlt10)
{
  res_pred <- .predict_or_extract_ml_train10(
    mlt10 = mlt10,
    type = "prob"
  )

  truth <- res_pred$truth
  lv <- res_pred$levels

  roc_list <- list()
  auc_list <- list()

  for (nm in names(res_pred$pred)) {
    prob <- res_pred$pred[[nm]]
    roc_obj <- pROC::roc(
      response = truth,
      predictor = prob,
      levels = rev(lv),
      quiet = TRUE
    )

    roc_list[[nm]] <- roc_obj
    auc_list[[nm]] <- data.frame(
      Model = nm,
      AUC = as.numeric(pROC::auc(roc_obj)),
      stringsAsFactors = FALSE
    )
  }

  auc_df <- do.call(rbind, auc_list)
  rownames(auc_df) <- NULL
  auc_df <- auc_df[order(-auc_df$AUC), ]

  out <- list(roc = roc_list, summary = auc_df)
  class(out) <- "ml_roc"
  out
}

# ===============================================================
# Plot single confusion matrix
# ===============================================================

.plot_confusion_ml_train10 <- function(el_obj, ...)
{
  if (!is(el_obj, "ml_evaluation")) {
    stop('!is(el_obj, "ml_evaluation").')
  }
  df_all <- list()
  for (nm in names(el_obj$matrix)) {
    cm <- el_obj$matrix[[nm]]
    if (is.null(cm)) {
      next
    }
    df <- as.data.frame(cm$table)
    colnames(df) <- c("Prediction", "Reference", "Freq")
    df$Model <- nm
    df_all[[nm]] <- df
  }
  dat_plot <- do.call(rbind, df_all)
  rownames(dat_plot) <- NULL
  set.seed(123L)
  ggplot2::ggplot(dat_plot, ggplot2::aes(x = Reference, y = Prediction, fill = Freq)) +
    ggplot2::geom_tile() +
    scale_fill_gradient(low = "lightblue", high = "Blue4") +
    ggplot2::geom_text(
      ggplot2::aes(label = Freq), color = "white", size = 4
    ) +
    ggplot2::facet_wrap(~ Model, ...) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Confusion Matrices",
      x = "Reference",
      y = "Prediction"
    )
}

.plot_evaluation_ml_train10 <- function(el_obj, ...)
{
  if (!is(el_obj, "ml_evaluation")) {
    stop('!is(el_obj, "ml_evaluation").')
  }
  dat <- el_obj$summary
  ord <- dat$Model[order(dat$Accuracy, decreasing = TRUE)]

  dat_long <- tidyr::pivot_longer(
    dat,
    cols = c(
      Accuracy,
      Kappa,
      Sensitivity,
      Specificity,
      F1
    ),
    names_to = "Metric",
    values_to = "Value"
  )

  dat_long <- dplyr::mutate(
    dat_long,
    Model = factor(Model, levels = unique(ord)),
    Is_NA = is.na(Value),
    Label = dplyr::if_else(
      Is_NA,
      "NA",
      as.character(round(Value, 3))
    ),
    Value_plot = dplyr::if_else(
      Is_NA,
      0,
      Value
    )
  )

  ggplot(dat_long, aes(x = Model, y = Value_plot)) +
    geom_col(
      data = dplyr::filter(dat_long, !Is_NA),
      width = 0.75
    ) +
    geom_col(
      data = dplyr::filter(dat_long, Is_NA),
      width = 0.75,
      alpha = 0.25
    ) +
    geom_text(
      aes(label = Label, hjust = ifelse(Value_plot > 0, -.1, 1.1)),
      size = 3
    ) +
    coord_flip(clip = "off") +
    facet_wrap(~ Metric, scales = "free_x", ...) +
    labs(
      title = "Classification Performance Across Models",
      x = "Model",
      y = "Score"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold")
      ) +
    expand_limits(y = 1.08)
}


# ===============================================================
# Plot all ROC curves
# ===============================================================

.plot_roc_ml_train10 <- function(roc_obj)
{
  if (!is(roc_obj, "ml_roc")) {
    stop("!is(roc_obj, 'ml_roc').")
  }

  df_all  <- list()
  auc_lab <- character(0)

  for (nm in names(roc_obj$roc)) {

    r <- roc_obj$roc[[nm]]

    df_all[[nm]] <- data.frame(
      FPR   = 1 - r$specificities,
      TPR   = r$sensitivities,
      Model = nm,
      stringsAsFactors = FALSE
    )

    auc_val <- as.numeric(r$auc)

    auc_lab[nm] <- sprintf(
      "%s (AUC = %.3f)",
      nm,
      auc_val
    )
  }

  df_plot <- do.call(rbind, df_all)
  rownames(df_plot) <- NULL

  df_plot$Model <- factor(
    df_plot$Model,
    levels = names(auc_lab),
    labels = auc_lab
  )

  df_plot <- df_plot[
    order(df_plot$Model, df_plot$FPR, df_plot$TPR),
  ]

  cols <- color_set()

  if (length(cols) < length(auc_lab)) {
    cols <- rep(cols, length.out = length(auc_lab))
  }

  names(cols) <- auc_lab

  ggplot2::ggplot(
    df_plot,
    ggplot2::aes(
      x = FPR,
      y = TPR,
      color = Model
    )
  ) +
    ggplot2::geom_step(
      direction = "vh",
      linewidth = 1.15
    ) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = 2,
      linewidth = 0.7,
      color = "grey40"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      expand = c(0, 0)
    ) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::labs(
      title = "ROC Curves",
      x = "False Positive Rate",
      y = "True Positive Rate",
      color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0
      ),
      panel.grid.minor = ggplot2::element_blank(),

      legend.position = c(0.98, 0.02),
      legend.justification = c(1, 0),

      legend.background = ggplot2::element_rect(
        fill = grDevices::adjustcolor("white", alpha.f = 0.75),
        color = "grey80"
      ),
      legend.key = ggplot2::element_blank(),

      legend.text = ggplot2::element_text(size = 10)
    )
}



# ------------------------------------------------------------------
# Dynamic cfg (compact version)
# ------------------------------------------------------------------
make_cfg <- function(data, tune_length = 5L) {

  n_sample  <- nrow(data)
  p_feature <- ncol(data)

  # --------------------------------------------------------------
  # Sample-size scaling factor
  # smaller sample -> stronger regularization
  # --------------------------------------------------------------
  small_mode <- n_sample <= 40
  mid_mode   <- n_sample > 40 & n_sample <= 80

  # --------------------------------------------------------------
  # Shared dynamic parameters
  # --------------------------------------------------------------
  ntree <- if (small_mode) 300L else if (mid_mode) 500L else 1000L

  svm_C <- if (small_mode) {
    c(0.01, 0.1, 0.5)
  } else if (mid_mode) {
    c(0.01, 0.1, 1, 5)
  } else {
    c(0.01, 0.1, 1, 10)
  }

  knn_k <- if (small_mode) {
    c(5, 7, 9)
  } else {
    c(3, 5, 7, 9)
  }

  cp_seq <- if (small_mode) {
    c(0.05, 0.1, 0.2)
  } else {
    c(0.01, 0.05, 0.1)
  }

  lambda_seq <- if (small_mode) {
    seq(0.01, 1, length.out = 15)
  } else if (mid_mode) {
    seq(0.005, 0.5, length.out = 15)
  } else {
    seq(0.001, 0.2, length.out = 20)
  }

  alpha_en <- if (small_mode) {
    c(0.1, 0.3, 0.5)
  } else {
    c(0.2, 0.5, 0.8)
  }

  nnet_size <- if (small_mode) c(1, 2) else c(1, 2, 3)

  nnet_decay <- if (small_mode) {
    c(0.5, 1, 2)
  } else {
    c(0.1, 0.5, 1)
  }

  gbm_depth <- if (small_mode) c(1) else if (mid_mode) c(1, 2) else c(1, 2, 3)
  gbm_tree <- if (small_mode) c(20, 50) else if (mid_mode) c(50, 100) else c(100, 200)
  gbm_lr <- if (small_mode) c(0.03) else c(0.03, 0.05)
  gbm_node <- if (small_mode) c(5, 8) else c(5, 10)

  # --------------------------------------------------------------
  # Final cfg
  # --------------------------------------------------------------
  cfg <- list(

    LR = list(
      method = "glm"
    ),

    DT = list(
      method = "rpart",
      tuneGrid = expand.grid(
        cp = cp_seq
      )
    ),

    SVM = list(
      method = "svmLinear",
      tuneGrid = expand.grid(
        C = svm_C
      )
    ),

    RF = list(
      method = "rf",
      ntree = ntree,
      tuneGrid = expand.grid(
        mtry = unique(
          pmax(1, c(2, floor(sqrt(p_feature))))
        )
      )
    ),

    KNN = list(
      method = "knn",
      tuneGrid = expand.grid(
        k = knn_k
      )
    ),

    BT = list(
      method = "rf",
      ntree = ntree,
      tuneGrid = data.frame(
        mtry = p_feature
      )
    ),

    LDA = list(
      method = "lda"
    ),

    NNET = list(
      method = "nnet",
      trace = FALSE,
      MaxNWts = 5000,
      tuneGrid = expand.grid(
        size = nnet_size,
        decay = nnet_decay
      )
    ),

    NB = list(
      method = "nb",
      tuneLength = tune_length
    ),

    EN = list(
      method = "glmnet",
      tuneGrid = expand.grid(
        alpha = alpha_en,
        lambda = lambda_seq
      )
    ),

    Lasso = list(
      method = "glmnet",
      tuneGrid = expand.grid(
        alpha = 1,
        lambda = lambda_seq
      )
    ),

    RR = list(
      method = "glmnet",
      tuneGrid = expand.grid(
        alpha = 0,
        lambda = lambda_seq
      )
    ),

    GBM = list(
      method = "gbm",
      verbose = FALSE,
      tuneGrid = expand.grid(
        interaction.depth = gbm_depth,
        n.trees = gbm_tree,
        shrinkage = gbm_lr,
        n.minobsinnode = gbm_node
      )
    )
  )

  return(cfg)
}


