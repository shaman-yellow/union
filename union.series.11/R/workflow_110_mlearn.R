# ==========================================================================
# workflow of mlearn
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_mlearn <- setClass("job_mlearn", 
  contains = c("job"),
  prototype = prototype(
    pg = "mlearn",
    info = c(""),
    cite = "",
    method = "",
    tag = "mlearn",
    analysis = "Machine-Learning 机器学习"
    ))

setGeneric("asjob_mlearn",
  function(x, ...) standardGeneric("asjob_mlearn"))

setMethod("asjob_mlearn", signature = c(x = "job_deseq2"),
  function(x, ref, group = "group", levels = rev(.guess_compare_deseq2(x, 1L)), seed = 987456L)
  {
    if (x@step < 1L) {
      stop('x@step < 1L.')
    }
    object <- x$vst
    if (is.null(object)) {
      stop('is.null(object).')
    }
    data <- SummarizedExperiment::assay(object)
    if (is(ref, "feature")) {
      snap <- snap(ref)
      ref <- resolve_feature(ref)
      if (length(ref) <= 2) {
        stop('length(ref) <= 2, too few genes.')
      }
      snapAdd_onExit("x", "以 {x$project} 为训练集，将{snap}用于机器学习筛选关键基因。")
    }
    if (any(!ref %in% rownames(data))) {
      stop('any(!ref %in% rownames(data)).')
    }
    data <- t(data[ rownames(data) %in% ref, ])
    metadata <- data.frame(object@colData)
    levels <- eval(levels)
    project <- x$project
    x <- .job_mlearn(object = data)
    x$metadata <- metadata
    x$project <- project
    x$levels <- levels
    x$target <- factor(metadata[[ group ]], levels = levels)
    x$seed <- seed
    return(x)
  })

setMethod("asjob_mlearn", signature = c(x = "job_limma"),
  function(x, ref, group = "group", levels = rev(.guess_compare_limma(x, 1L)), seed = 987456L)
  {
    if (x@step < 1L) {
      stop('x@step < 1L.')
    }
    object <- x$normed_data
    if (is.null(object)) {
      stop('is.null(object).')
    }
    data <- object$E
    if (is(ref, "feature")) {
      snap <- snap(ref)
      ref <- resolve_feature(ref)
      if (length(ref) <= 2) {
        stop('length(ref) <= 2, too few genes.')
      }
      snapAdd_onExit("x", "以 {x$project} 为训练集，将{snap}用于机器学习筛选关键基因。")
    }
    if (any(!ref %in% rownames(data))) {
      stop('any(!ref %in% rownames(data)).')
    }
    data <- t(data[ rownames(data) %in% ref, ])
    metadata <- data.frame(object$targets)
    levels <- eval(levels)
    project <- x$project
    x <- .job_mlearn(object = data)
    x$project <- project
    x$metadata <- metadata
    x$levels <- levels
    x$target <- factor(metadata[[ group ]], levels = levels)
    x$seed <- seed
    return(x)
  })

setMethod("step0", signature = c(x = "job_mlearn"),
  function(x){
    step_message("Prepare your data with function `job_mlearn`.")
  })

setMethod("step1", signature = c(x = "job_mlearn"),
  function(x, subset_sizes = 15:50, n = 10, method = "cv", kernel = "linear", seed = x$seed,
    workers = NULL, ..., rerun = FALSE, skip = FALSE)
  {
    step_message("SVM-RFE.")
    if (skip) {
      return(x)
    }
    data <- object(x)
    target <- x$target
    args <- as.list(environment())
    args$rerun <- args$x <- NULL
    dir.create("tmp", FALSE)
    svm_rfe <- expect_local_data(
      "tmp", "svm_rfe", .run_svm_rfe, args, rerun = rerun
    )
    x$svm_rfe <- svm_rfe
    data_svm <- as.data.frame(svm_rfe$results)
    data_svm <- data_svm[data_svm$Variables %in% subset_sizes, , drop = FALSE]
    data_svm <- data_svm[order(data_svm$Variables), , drop = FALSE]
    data_svm$Error <- 1 - data_svm$Accuracy
    p.acc <- .plot_svm_rfe_metric(
      data_svm,
      y = "Accuracy",
      ylab = glue::glue("{n} × CV Accuracy"),
      best = "max", subset_sizes = subset_sizes
    )
    p.error <- .plot_svm_rfe_metric(
      data_svm,
      y = "Error",
      ylab = glue::glue("{n} × CV Error Rate"),
      best = "min", subset_sizes = subset_sizes
    )
    p.svm <- patchwork::wrap_plots(p.acc, p.error, ncol = 2L)
    p.svm <- set_lab_legend(
      wrap(p.svm, 7, 4),
      glue::glue("{x@sig} SVM-RFE candidate subset sizes evaluation"),
      glue::glue("SVM-RFE候选子集正确率与错误率曲线|||{n}折交叉验证准确率与错误率（{n}x CV Accuracy）随特征数量变化的趋势。")
    )
    svm_rfe_res <- list(best_size = svm_rfe$bestSubset, features = caret::predictors(svm_rfe))
    x$svm_rfe_res <- svm_rfe_res
    x$t.svm_rfe_accuracy <- dplyr::arrange(
      as_tibble(svm_rfe$results), dplyr::desc(Accuracy)
    )
    x <- plotsAdd(x, p.svm)
    x <- methodAdd(x, "以 R 包 `e1071` ⟦pkgInfo('e1071')⟧ 构建支持向量机递归特征消除模型 (SVM-RFE)。核函数设定为线性核 (kernel = {kernel}) ；以 R 包 `caret` ⟦pkgInfo('caret')⟧ 实现递归特征消除流程，采用 {n} 折交叉验证评估模型分类性能，依据分类准确率 (Accuracy) 从高到低排序筛选最优特征基因子集。")
    x <- snapAdd(
      x, "SVM-RFE 最佳子集数为 {svm_rfe_res$best_size}{aref(p.svm)}，准确率 (Accuracy) 为 {round(x$t.svm_rfe_accuracy$Accuracy[1], 3)}，误差值 (AccuracySD) 为 {round(x$t.svm_rfe_accuracy$AccuracySD[1], 3)}，对应 feature 为：{bind(svm_rfe_res$features)}。\n\n\n\n"
    )
    return(x)
  })

setMethod("step2", signature = c(x = "job_mlearn"),
  function(x, n = 10, lambda.type = c("1se", "min"), alpha = 1,
    style_plot = c("internal", "ggplot"), seed = x$seed, ...)
  {
    step_message("Lasso")
    lambda.type <- match.arg(lambda.type)
    lambda.type <- paste0("lambda.", lambda.type)
    data <- object(x)
    target <- x$target
    set.seed(seed)
    cv_lasso <- e(glmnet::cv.glmnet(
        x = data, y = target,
        family = "binomial", alpha = alpha, nfolds = n,
        # type.measure = "deviance",
        standardize = TRUE, parallel = FALSE
        ))
    lambda <- cv_lasso[[ lambda.type ]]
    coefs <- coef(cv_lasso, s = lambda)
    coefs_matrix <- as.matrix(coefs)
    whichCoefs <- which(coefs_matrix[, 1] != 0 & rownames(coefs_matrix) != "(Intercept)")
    selected <- rownames(coefs_matrix)[ whichCoefs ]
    # coef_values <- coefs_matrix[selected, 1]
    x$lasso_res <- list(
      cv_lasso = cv_lasso, coefs = coefs_matrix, features = selected, type = lambda.type
    )
    style_plot <- match.arg(style_plot)
    if (style_plot == "internal") {
      expr <- expression({
        fun <- function() {
          cv <- cv_lasso
          requireNamespace("glmnet")
          suffix <- c("1se", "min")
          types <- paste0("lambda.", suffix)
          y <- max(cv$cvm)
          lambdas <- vapply(types, function(x) cv[[x]], double(1))
          x <- log(lambdas)
          labels <- glue::glue("log(λ) ({suffix})\n = {signif(log(lambdas), 2)}")
          plot(cv, sign.lambda = 1)
          text(x, y, labels, adj = 1)
        }
        fun()
      })
      p.lasso_cv <- as_grob(expr, environment())
      expr <- expression({
        fun <- function() {
          cv <- cv_lasso
          requireNamespace("glmnet")
          suffix <- c("1se", "min")
          types <- paste0("lambda.", suffix)
          lambdas <- vapply(types, function(x) cv[[x]], double(1))
          x <- log(lambdas)
          if (any(formalArgs(glmnet:::plot.glmnet) == "sign.lambda")) {
            plot(cv$glmnet.fit, sign.lambda = 1, label = FALSE, xvar = "lambda")
          } else {
            plot(cv$glmnet.fit, label = FALSE, xvar = "lambda")
          }
          abline(v = x, lty = 2)
          labels <- glue::glue("log(λ) ({suffix})\n = {signif(log(lambdas), 2)}")
          text(x, par("usr")[4] * 0.7, labels, adj = 1)
        }
        fun()
      })
      p.coefs_path <- as_grob(expr, environment())
    } else {
      ps.lasso <- .plot_cv_glmnet_with_ggStyle(cv_lasso)
      p.lasso_cv <- ps.lasso$p_cv
      p.coefs_path <- ps.lasso$p_coef
    }
    p.lasso_cv <- set_lab_legend(
      wrap(p.lasso_cv, 5.5, 4, showtext = TRUE),
      glue::glue("{x@sig} LASSO Cross Validation"),
      glue::glue("LASSO 交叉验证误差|||Lasso 回归模型的交叉验证图，用于选择正则化参数 λ。图中展示了不同 λ 值下的 {cv_lasso$name}。横坐标是log(λ)，即正则化参数 λ 的对数值。随着 λ 值的增加，模型的复杂度降低，正则化强度增加。纵坐标是二项式偏差。")
    )
    p.coefs_path <- set_lab_legend(
      wrap(p.coefs_path, if (style_plot == "internal") 5.5 else 8, 4, showtext = TRUE),
      glue::glue("{x@sig} Lasso Coefficient path"),
      glue::glue("LASSO 系数路径|||Lasso 回归系数路径图，展示了不同特征的系数随正则化参数 log(λ) 变化的情况。横坐标是 log(λ)，纵坐标是模型中各个特征的系数值。随着 λ 值的增加（从右到左），更多的特征系数被压缩至零，这是Lasso回归的特征选择过程。")
    )
    x <- plotsAdd(x, p.lasso_cv, p.coefs_path)
    prin <- if (lambda.type == "lambda.1se") "1-SE" else "最小误差"
    if (alpha == 1) {
      x <- methodAdd(
        x, "以 R 包 glmnet ⟦pkgInfo('glmnet')⟧ 开展 LASSO 逻辑回归分析。设置 α = 1 实现 L1 正则化，通过 {n} 折交叉验证结合{prin}准则确定最优 λ 值 (λ = {fmt(cv_lasso[[ lambda.type ]])}, 
        Log(λ) = {signif(log(cv_lasso[[ lambda.type ]]), 2)}) 。"
)
    } else if (alpha < 1) {
      x <- methodAdd(
        x, "以 R 包 glmnet ⟦pkgInfo('glmnet')⟧ 开展 LASSO (Elastic Net) 回归分析。设置 α = {alpha} 实现 L1 与 L2 正则化的加权组合，通过 {n} 折交叉验证结合{prin}准则确定最优 λ 值 (λ = {fmt(cv_lasso[[ lambda.type ]])}, Log(λ) = {signif(log(cv_lasso[[ lambda.type ]]), 2)})。"
      )
    } else {
      stop('alpha?')
    }
    x <- snapAdd(x, "LASSO 筛选的核心 feature（非零系数）数量为 {length(selected)}{aref(p.lasso_cv)}，对应为：{bind(selected)}。\n\n\n\n")
    return(x)
  })

setMethod("step3", signature = c(x = "job_mlearn"),
  function(x, ntree = 1000, top = 10, seed = x$seed, ...)
  {
    step_message("Random Forest.")
    data <- object(x)
    target <- x$target
    mtry = floor(sqrt(ncol(data)))
    set.seed(seed)
    rf_model <- e(randomForest::randomForest(
      x = data, y = target, ntree = ntree, mtry = mtry,
      importance = TRUE, proximity = FALSE,
      oob.prox = FALSE, keep.forest = TRUE
    ))
    error_data <- as_tibble(rf_model$err.rate)
    error_data <- dplyr::mutate(error_data, trees = seq_len(nrow(error_data)))
    error_data <- tidyr::pivot_longer(error_data, -trees, names_to = "Error_Type", values_to = "Error_Rate")
    p.error <- ggplot(error_data, aes(x = trees, y = Error_Rate, color = Error_Type)) +
      geom_line() +
      labs(x = "Number of trees", y = "Error Rate") +
      theme_minimal() + theme(legend.title = element_blank())
    p.error <- set_lab_legend(
      wrap(p.error, 5, 3),
      glue::glue("{x@sig} Trend of random forest error rate"),
      glue::glue("随机森林误差率随树数量变化趋势图|||在训练过程中模型对不同组别识别的错误概。OOB 为总体袋外误差（OOB error），即所有类别的平均误差率。随着树的数量增加，总体误差率逐渐趋于稳定。")
    )
    importance_df <- as_tibble(
      randomForest::importance(rf_model), idcol = "feature"
    )
    importance_df <- dplyr::arrange(importance_df, dplyr::desc(MeanDecreaseGini))
    t.tops <- head(importance_df, top)
    t.tops <- set_lab_legend(
      t.tops,
      glue::glue("{x@sig} top importance feature"),
      glue::glue("按 MeanDecreaseGini 降低排序的Top Feature。")
    )
    p.tops <- .plot_rf_importance(importance_df, top)
    p.tops <- set_lab_legend(
      wrap(p.tops, 5, top * .2 + 2),
      glue::glue("{x@sig} RF top importance feature"),
      glue::glue("随机森林（Random Forest, RF）特征重要性分析图|||基于随机森林模型计算各特征的重要性评分，并按 MeanDecreaseGini 指标降序排列展示 Top 特征。其中，MeanDecreaseAccuracy 表示变量对模型预测准确率的贡献程度，数值越高说明该变量对分类性能影响越显著；MeanDecreaseGini 表示变量对节点纯度提升的贡献程度，数值越高说明该变量在随机森林决策过程中具有更强的区分能力。")
    )
    x$rf_res <- list(rf_model = rf_model, features = t.tops$feature)
    x <- tablesAdd(x, t.tops)
    x <- plotsAdd(x, p.error, p.tops)
    x <- methodAdd(x, "以 R 包 `randomForest` ⟦pkgInfo('randomForest')⟧ 构建随机森林分类模型，设定决策树数量（ntree）为 {ntree}，特征选择数 (mtry) 为基因总数的平方根，通过袋外数据 (OOB) 评估模型误差率，计算 Feature 重要性评分，筛选相对重要性 top {top}；同时分析分类树数量与误差率的关联趋势，确定模型最优复杂度。")
    x <- snapAdd(x, "随机森林特征重要性 Top {top} 基因 (PMID: 37065165; PMID: 16398926; PMID: 41243474)：{bind(t.tops$feature)}{aref(p.error)}。\n\n\n\n")
    return(x)
  })

setMethod("step4", signature = c(x = "job_mlearn"),
  function(x, n = 10, seed = x$seed,
    early = 5L, eta = .05, nrounds = 100L, rerun = FALSE)
  {
    step_message("XGBoost")
    data <- object(x)
    target <- x$target
    fun_xgb <- function(...) {
      .mlearn_alter_xgboost(
        data, target, nfold = n, seed = seed, early_stopping_rounds = early, 
        eta = eta, nrounds = nrounds
      )
    }
    res <- expect_local_data(
      "tmp", "xgboost", fun_xgb,
      list(n, seed, target, nrounds, eta, early), rerun = rerun
    )
    eval_log <- res$cv$evaluation_log
    p.importance <- ggplot(res$importance, aes(x = reorder(Feature, Gain), y = Gain)) +
      geom_col() +
      coord_flip() +
      theme_bw() +
      xlab("Gene") +
      ylab("Importance (Gain)") +
      ggtitle("Feature Importance")
    p.importance <- set_lab_legend(
      wrap_scale(p.importance, 20, nrow(res$importance), size = .1),
      glue::glue("{x@sig} XGBoost Feature Importance"),
      glue::glue("XGBoost 模型特征重要性排序|||展示对分类任务贡献度最高的基因（按 Gain 值衡量），横轴表示特征在模型中的相对重要性，纵轴为基因名称，重要性越高表示该基因在模型决策中贡献越大。")
    )
    data <- dplyr::select(
      eval_log, iter, Validate = test_auc_mean, Train = train_auc_mean
    )
    data <- tidyr::pivot_longer(data, -iter, names_to = "type", values_to = "value")
    p.auc <- ggplot(data, aes(x = iter, y = value, color = type)) +
      geom_line() +
      theme_bw() +
      labs(x = "Iteration", y = "AUC", color = "Type") +
      ggtitle("Cross-validation AUC")
    p.auc <- set_lab_legend(
      wrap(p.auc, 5, 3.5),
      glue::glue("{x@sig} XGBoost Cross-validation AUC"),
      glue::glue("交叉验证模型性能迭代曲线图|||横轴为迭代轮数（Number of Trees），纵轴为模型在验证集上的 AUC 值，用于评估模型随训练过程的收敛趋势及最优迭代轮数的选择。")
    )
    x <- methodAdd(
      x, "以 R 包 `xgboost` ⟦pkgInfo('xgboost')⟧ 构建梯度提升树二分类模型，设定最大迭代轮数（nrounds）为 {nrounds}，学习率（eta）为 {eta}，最大树深（max_depth）为 4，并结合 {n} 折交叉验证与早停策略（early stopping）确定最优迭代轮数；基于模型计算特征重要性（Gain），筛选重要基因；同时分析迭代轮数与模型性能（AUC）变化趋势，以评估模型收敛过程与复杂度。"
    )
    features <- res$selected_genes
    x <- snapAdd(x, "XGBoost 所有重要基因 (n = {length(features)})：{bind(features)}{aref(p.importance)}。\n\n\n\n")
    x$xgb_res <- list(
      rf_model = res$model, rf_cv = res$cv, features = features
    )
    x <- plotsAdd(x, p.auc, p.importance)
    return(x)
  })

.plot_cv_glmnet_with_ggStyle <- function(cv, n_label_genes = 10L)
{

  if (!inherits(cv, "cv.glmnet")) {
    stop("Input must be 'cv.glmnet' object.")
  }

  message(glue::glue(
    "Preparing plotting data for cv.glmnet object..."
  ))

  vec_lambda <- cv$lambda
  vec_log_lambda <- log(vec_lambda)

  data_cv <- data.frame(
    lambda = vec_lambda,
    log_lambda = vec_log_lambda,
    cvm = cv$cvm,
    cvsd = cv$cvsd,
    nzero = cv$nzero
  )

  data_cv$ymin <- data_cv$cvm - data_cv$cvsd
  data_cv$ymax <- data_cv$cvm + data_cv$cvsd

  vec_lambda_selected <- c(
    cv$lambda.min,
    cv$lambda.1se
  )

  vec_log_lambda_selected <- log(vec_lambda_selected)

  vec_suffix <- c("min", "1se")

  data_vline <- data.frame(
    lambda = vec_lambda_selected,
    log_lambda = vec_log_lambda_selected,
    suffix = vec_suffix
  )

  data_vline$label <- glue::glue(
    "lambda.{data_vline$suffix}\nlog(λ) = {signif(data_vline$log_lambda, 3L)}"
  )

  message(glue::glue(
    "Extracting coefficient matrix..."
  ))

  mat_beta <- as.matrix(
    stats::coef(cv$glmnet.fit)
  )

  vec_features <- rownames(mat_beta)

  mat_beta <- mat_beta[-1L, , drop = FALSE]
  vec_features <- vec_features[-1L]

  data_coef <- reshape2::melt(
    mat_beta,
    varnames = c("feature_index", "lambda_index"),
    value.name = "coef"
  )

  data_coef$feature <- vec_features[
    data_coef$feature_index
  ]

  data_coef$lambda <- vec_lambda[
    data_coef$lambda_index
  ]

  data_coef$log_lambda <- log(
    data_coef$lambda
  )

  message(glue::glue(
    "Detected {length(unique(data_coef$feature))} features."
  ))

  data_abs <- stats::aggregate(
    abs(coef) ~ feature,
    data = data_coef,
    FUN = max
  )

  data_abs <- data_abs[
    order(
      data_abs[, 2L],
      decreasing = TRUE
    ),
  ]

  vec_label_features <- head(
    data_abs$feature,
    n_label_genes
  )

  data_label <- data_coef[
    data_coef$feature %in% vec_label_features,
  ]

  vec_idx_label <- unlist(
    tapply(
      seq_len(nrow(data_label)),
      data_label$feature,
      function(x) {
        x[which.max(
          data_label$log_lambda[x]
        )]
      }
    )
  )

  data_label <- data_label[
    vec_idx_label,
  ]

  data_cv$nzero_factor <- factor(
    data_cv$nzero,
    levels = unique(data_cv$nzero)
  )

  p_cv <- ggplot(
    data_cv,
    aes(
      x = log_lambda,
      y = cvm
    )
  ) +
    geom_errorbar(
      aes(
        ymin = ymin,
        ymax = ymax,
        color = nzero_factor
      ),
      width = 0.03
    ) +
    geom_point(
      aes(color = nzero_factor),
      size = 1.8
    ) +
    geom_line(
      color = "grey40",
      linewidth = 0.5
    ) +
    geom_vline(
      data = data_vline,
      aes(xintercept = log_lambda),
      linetype = 2L
    ) +
    geom_text(
      data = data_vline,
      aes(
        x = log_lambda,
        y = max(data_cv$ymax) * 0.98,
        label = label
      ),
      hjust = -0.05,
      vjust = 1,
      size = 3
    ) +
    labs(
      x = "Log(λ)",
      y = cv$name,
      color = "Variables"
    ) +
    theme_bw()

  p_coef <- ggplot(
    data_coef,
    aes(
      x = log_lambda,
      y = coef,
      color = feature,
      group = feature
    )
  ) +
    geom_line(
      linewidth = 0.8
    ) +
    geom_vline(
      data = data_vline,
      aes(xintercept = log_lambda),
      linetype = 2L,
      color = "grey40"
    ) +
    geom_text(
      data = data_vline,
      aes(
        x = log_lambda,
        y = max(data_coef$coef) * 0.95,
        label = label
      ),
      hjust = -0.05,
      vjust = 1,
      inherit.aes = FALSE,
      size = 3
    ) +
    ggrepel::geom_text_repel(
      data = data_label,
      aes(label = feature),
      size = 3,
      show.legend = FALSE
    ) +
    labs(
      x = "Log(λ)",
      y = "Coefficients",
      color = "Feature"
    ) +
    theme_bw()

  list(p_cv = p_cv, p_coef = p_coef)
}

.plot_rf_importance <- function(data_importance, n_top) {
  data_importance <- dplyr::arrange(
    data_importance, dplyr::desc(MeanDecreaseGini)
  )
  data_importance <- head(data_importance, n_top)

  data_importance <- dplyr::arrange(
    data_importance, MeanDecreaseGini
  )

  data_plot <- tidyr::pivot_longer(
    data_importance,
    cols = c(
      MeanDecreaseAccuracy,
      MeanDecreaseGini
      ),
    names_to = "metric",
    values_to = "importance"
  )

  data_plot$feature <- factor(
    data_plot$feature,
    levels = data_importance$feature
  )

  message(
    glue::glue(
      "Prepared plotting table with {nrow(data_plot)} rows."
    )
  )

  p_rf <- ggplot2::ggplot(
    data_plot, ggplot2::aes(x = feature, y = importance, fill = metric)
    ) +
  ggplot2::geom_col(
    position = ggplot2::position_dodge(width = 0.75),
    width = 0.65
    ) +
  ggplot2::geom_text(
    ggplot2::aes(label = round(importance, 2L)),
    position = ggplot2::position_dodge(width = 0.75),
    hjust = -0.1,
    size = 3.5
    ) +
  ggplot2::scale_fill_manual(
    values = c("MeanDecreaseAccuracy" = "#3C77C4",
      "MeanDecreaseGini" = "#F28E2B"
    )) +
  ggplot2::labs(
    title = "RF Feature Importance",
    x = "Feature",
    y = "Importance Score",
    fill = "Type"
    ) +
  ggplot2::coord_flip() +
  ggplot2::theme_classic()
}

.mlearn_alter_xgboost <- function(
  data, target,
  nrounds = 100, nfold = 10,
  early_stopping_rounds = 5,
  eta = .05,
  seed = 123)
{
  set.seed(seed)
  # -----------------------------
  # Prepare data
  # -----------------------------
  X <- as.matrix(data)
  if (is.factor(target)) {
    y <- as.integer(target) - 1
  } else {
    y <- target
  }

  dtrain <- e(xgboost::xgb.DMatrix(data = X, label = y))

  # -----------------------------
  # Parameters (robust defaults)
  # -----------------------------
  params <- list(
    objective = "binary:logistic",
    eval_metric = "auc",
    max_depth = 4,
    eta = eta,
    subsample = 0.7,
    colsample_bytree = 0.4,
    lambda = 1,
    alpha = 0
  )

  # -----------------------------
  # Cross-validation
  # -----------------------------
  cv <- e(xgboost::xgb.cv(
    params = params,
    data = dtrain,
    nrounds = nrounds,
    nfold = nfold,
    early_stopping_rounds = early_stopping_rounds,
    verbose = 1
  ))

  eval_log <- cv$evaluation_log
  best_nrounds <- which.max(eval_log$test_auc_mean)

  # -----------------------------
  # Train final model
  # -----------------------------
  model <- e(xgboost::xgb.train(
    params = params,
    data = dtrain,
    nrounds = best_nrounds,
    verbose = 1
  ))

  # -----------------------------
  # Feature importance
  # -----------------------------
  importance <- e(xgboost::xgb.importance(
    model = model,
    feature_names = colnames(X)
  ))
  # -----------------------------
  # Output
  # -----------------------------
  list(
    model = model, cv = cv, best_nrounds = best_nrounds,
    importance = importance,
    selected_genes = importance$Feature
  )
}

run_mlean_with_seeds <- function(x, ..., ntry = 100, expect = 2)
{
  x@sig <- "test"
  x <- copy_job(x)
  lapply(seq_len(ntry),
    function(n) {
      seeds <- sample(1:100000, 3)
      capture.output({
        x <- suppressMessages(step1(x, seed = seeds[1], ...))
        x <- suppressMessages(step2(x, seed = seeds[2], ...))
        x <- suppressMessages(step3(x, seed = seeds[3], ...))
      })
      alls <- list(
        SVM_RFE = x$svm_rfe_res$features,
        LASSO = x$lasso_res$features,
        Random_Forest = x$rf_res$features
      )
      res <- ins(lst = alls)
      message(glue::glue("N = {n}, use seeds: {bind(seeds)}, got: {bind(res)}"))
      if (length(res) >= expect) {
        return(seeds)
      } else {
        NULL
      }
    })
}

setMethod("asjob_venn", signature = c(x = "job_mlearn"),
  function(x){
    job_venn(lst = feature(x), mode = "ck")
  })

setMethod("feature", signature = c(x = "job_mlearn"),
  function(x){
    lst <- list(SVM_RFE = x$svm_rfe_res$features,
      LASSO = x$lasso_res$features,
      Random_Forest = x$rf_res$features,
      XGBoost = x$xgb_res$features
    )
    lst <- lst[ !vapply(lst, is.null, logical(1)) ]
    as_feature(lst, "Machine Learning")
  })

.run_svm_rfe <- function(data, target, n,
  method = "cv", kernel = "linear", cost = 10,
  subset_sizes, workers = NULL, seed = 123, ...)
{
  # -----------------------------
  # Basic checks & preprocessing
  # -----------------------------
  set.seed(seed)

  if (!is.factor(target)) {
    stop('!is.factor(target).')
  }

  # Filter valid subset sizes
  subset_sizes <- subset_sizes[subset_sizes <= ncol(data)]
  if (length(subset_sizes) == 0) {
    stop("No valid subset_sizes after filtering.")
  }

  # -----------------------------
  # Seeds for reproducibility
  # -----------------------------
  seeds <- vector("list", length = n + 1)
  size <- length(subset_sizes) + 1
  for (i in seq_len(n)) {
    seeds[[i]] <- sample.int(100000L, size)
  }
  seeds[[n + 1]] <- sample.int(100000L, 1)

  # -----------------------------
  # Custom SVM functions for RFE
  # -----------------------------
  svm_funcs <- caret::caretFuncs

  svm_funcs$fit <- function(x, y, first, last, ...) {
    e1071::svm(x = x, y = y,
      kernel = kernel, scale = TRUE, cost = cost,
      probability = FALSE, ...
    )
  }

  svm_funcs$pred <- function(object, x) {
    predict(object, x)
  }

  svm_funcs$selectSize <- function(x, metric, maximize) {
    x <- x[x$Variables %in% subset_sizes, , drop = FALSE]
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

  svm_funcs$rank <- function(object, x, y) {
    # Feature importance only valid for linear kernel
    if (kernel != "linear" || is.null(object$coefs)) {
      return(data.frame(
        var = colnames(x),
        Overall = rep(NA_real_, ncol(x))
      ))
    }

    # Compute weight vector: w = t(coefs) %*% SV
    w <- t(object$coefs) %*% object$SV
    importance <- abs(as.vector(w))

    # Safety check
    if (length(importance) != ncol(x)) {
      stop(
        sprintf(
          "Invalid importance length: got %d, expected %d. Usually caused by multiclass SVM or non-linear kernel.",
          length(importance), ncol(x)
        )
      )
    }

    data.frame(var = colnames(x), Overall = importance)
  }

  # -----------------------------
  # RFE control
  # -----------------------------
  ctrl <- caret::rfeControl(
    functions = svm_funcs,
    method = method,
    number = n,
    seeds = seeds,
    verbose = TRUE,
    allowParallel = TRUE
  )

  # -----------------------------
  # Parallel setup
  # -----------------------------
  if (!is.null(workers)) {
    workers <- min(workers, parallel::detectCores() - 1)

    cl <- e(parallel::makeCluster(workers))
    lib_paths <- .libPaths()
    parallel::clusterExport(cl, "lib_paths", envir = environment())

    # Ensure workers have required packages
    e(parallel::clusterEvalQ(cl, {
      .libPaths(lib_paths)
      require(caret)
      require(e1071)
    }))

    e(doParallel::registerDoParallel(cl))
    on.exit({
      parallel::stopCluster(cl)
      foreach::registerDoSEQ()
    }, add = TRUE)
  }

  # -----------------------------
  # Run RFE
  # -----------------------------
  res <- e(caret::rfe(x = data, y = target,
    sizes = subset_sizes, rfeControl = ctrl,
    metric = "Accuracy", maximize = TRUE
  ))

  return(res)
}

.plot_svm_rfe_metric <- function(data_metric, y, ylab, best = c("max", "min"),
  line_colour = "#2b9fd8", mark_colour = "#b83b3b", subset_sizes)
{
  best <- match.arg(best)
  vec_y <- data_metric[[y]]

  n_best <- if (best == "max") {
    which.max(vec_y)
  } else {
    which.min(vec_y)
  }

  data_best <- data_metric[n_best, , drop = FALSE]
  data_best$label <- paste0(
    data_best$Variables,
    " – ",
    formatC(data_best[[y]], format = "f", digits = 3L)
  )

  n_y_range <- diff(range(data_metric[[y]], na.rm = TRUE))
  if (!is.finite(n_y_range) || n_y_range == 0) {
    n_y_range <- abs(data_best[[y]]) * 0.05
  }
  if (!is.finite(n_y_range) || n_y_range == 0) {
    n_y_range <- 0.05
  }

  vec_x_range <- range(data_metric$Variables, na.rm = TRUE)
  n_x_mid <- mean(vec_x_range)

  n_push_x <- if (data_best$Variables <= n_x_mid) 1.6 else -1.6
  n_push_y <- n_y_range * 0.08
  n_hjust <- if (n_push_x > 0) 0 else 1

  data_best$x_text <- data_best$Variables + n_push_x
  data_best$y_text <- data_best[[y]] + n_push_y
  data_best$x_seg <- data_best$Variables + n_push_x * 0.7
  data_best$y_seg <- data_best[[y]] + n_push_y * 0.7

  ggplot2::ggplot(data_metric, ggplot2::aes(x = Variables, y = .data[[y]])) +
    ggplot2::geom_line(colour = line_colour, linewidth = 0.7) +
    ggplot2::geom_point(
      data = data_best,
      ggplot2::aes(x = Variables, y = .data[[y]]),
      inherit.aes = FALSE,
      shape = 21,
      size = 2.8,
      stroke = 0.9,
      colour = mark_colour,
      fill = "white"
    ) +
    ggplot2::geom_segment(
      data = data_best,
      ggplot2::aes(
        x = Variables,
        y = .data[[y]],
        xend = x_seg,
        yend = y_seg
      ),
      inherit.aes = FALSE,
      linewidth = 0.5,
      colour = mark_colour
    ) +
    ggplot2::geom_text(
      data = data_best,
      ggplot2::aes(
        x = x_text,
        y = y_text,
        label = label
      ),
      inherit.aes = FALSE,
      hjust = n_hjust,
      vjust = 0,
      colour = mark_colour,
      size = 3.6,
      fontface = "bold"
    ) +
    ggplot2::scale_x_continuous(
      limits = range(subset_sizes),
      breaks = pretty(range(subset_sizes), n = 5L)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(
        min(data_metric[[y]], na.rm = TRUE) - n_y_range * 0.12,
        max(data_metric[[y]], na.rm = TRUE) + n_y_range * 0.22
      )
    ) +
    ggplot2::labs(x = NULL, y = ylab) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(
        fill = NA,
        colour = "grey35",
        linewidth = 0.5
      ),
      axis.line = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(colour = "black"),
      axis.title = ggplot2::element_text(colour = "black"),
      plot.margin = ggplot2::margin(6, 8, 6, 6)
    )
}

# .run_svm_rfe <- function(data, target, n, method, kernel, 
#   subset_sizes, workers, seed)
# {
#   set.seed(seed)
#   seeds <- vector(mode = "list", length = n + 1)
#   size <- length(subset_sizes) + 1
#   for (i in seq_len(n)) {
#     seeds[[i]] <- sample.int(100000L, size)
#   }
#   seeds[[ n + 1 ]] <- sample.int(100000L, 1)
#   ctrl <- e(caret::rfeControl(functions = caret::caretFuncs,
#       method = method, number = n, seeds = seeds, verbose = TRUE, allowParallel = TRUE))
#   svm_funcs <- caret::caretFuncs
#   svm_funcs$fit <- function(x, y, first, last, ...) {
#     e1071::svm(x, y, kernel = kernel, scale = TRUE, probability = TRUE, ...)
#   }
#   svm_funcs$pred <- function(object, x) {
#     predict(object, x)
#   }
#   svm_funcs$rank <- function(object, x, y) {
#     # Calculate feature weights (coefficients of linear kernels)
#     if (is.null(object$coefs)) {
#       # If there are no coefficients, return a random ranking
#       data.frame(var = colnames(x), Overall = runif(ncol(x)))
#     } else {
#       # Calculate feature weights: w = t(x) %*% coefs
#       w <- t(object$coefs) %*% object$SV
#       if (ncol(w) == ncol(x)) {
#         importance <- abs(as.numeric(w))
#       } else {
#         importance <- rep(0, ncol(x))
#       }
#       data.frame(var = colnames(x), Overall = importance)
#     }
#   }
#   subset_sizes <- subset_sizes[subset_sizes <= (ncol(data) - 1)]
#   # run
#   # if (missing(workers) || is.null(workers)) {
#   #   n_cores <- e(parallel::detectCores()) - 1
#   # } else {
#   #   n_cores <- workers
#   # }
#   # if (n_cores > 10) {
#   #   stop('n_cores > 10, too many cores set to run.')
#   # }
#   if (!is.null(workers)) {
#     cl <- e(parallel::makeCluster(workers))
#     e(doParallel::registerDoParallel(cl))
#   }
#   res <- e(caret::rfe(x = data, y = target,
#       sizes = subset_sizes, rfeControl = ctrl,
#       metric = "Accuracy", maximize = TRUE, funcs = svm_funcs
#       ))
#   if (!is.null(workers)) {
#     e(parallel::stopCluster(cl))
#     e(foreach::registerDoSEQ())
#   }
#   res
# }

setMethod("set_remote", signature = c(x = "job_mlearn"),
  function(x, wd)
  {
    x$wd <- wd
    return(x)
  })
