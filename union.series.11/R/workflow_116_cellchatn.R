# ==========================================================================
# workflow of cellchatn
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_cellchatn <- setClass("job_cellchatn", 
  contains = c("job"),
  prototype = prototype(
    info = c("https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html"),
    cite = "[@InferenceAndAJinS2021]",
    method = "R package `CellChat` used for cell communication analysis",
    tag = "scrna:cellchat",
    analysis = "CellChat 细胞通讯分析"
    ))

setGeneric("asjob_cellchatn",
  function(x, ...) standardGeneric("asjob_cellchatn"))

setMethod("asjob_cellchatn", signature = c(x = "list"),
  function(x, nms = names(x))
  {
    if (any(!vapply(x, is, logical(1), "job_cellchat"))) {
      message("All elements of list should be 'job_cellchat'.")
    }
    if (any(vapply(x, function(x) x@step < 2, logical(1)))) {
      stop('any(vapply(x, function(x) x@step < 2, logical(1))).')
    }
    meth <- bind(meth(x[[1]])$step0, meth(x[[1]])$step1, co = "\n")
    p.lr_comm_bubbles <- lapply(x, 
      function(x) {
        x@plots$step2$lr_comm_bubble
      })
    t.lr_comm_bubble <- lapply(x, 
      function(x) {
        x@tables$step2$t.lr_comm_bubble
      })
    group.by <- vapply(x, function(x) x$group.by, character(1))
    if (!length(group.by <- unique(group.by))) {
      stop('!length(unique(group.by)).')
    }
    args_inters <- sapply(c("count", "weight"), simplify = FALSE,
      function(type) {
        weight.max <- e(CellChat::getMaxWeight(
            lapply(x, function(x) object(x)), attribute = c("idents", type)
            ))
        nets <- lapply(x, function(x) object(x)@net[[type]])
        list(nets = nets, wmax = weight.max, type = type)
      })
    object <- e(CellChat::mergeCellChat(lapply(x, function(x) object(x)), add.names = nms))
    x <- .job_cellchatn(object = object)
    x$.pre_meth <- meth
    x$each_lr_comm <- namel(p.lr_comm_bubbles, t.lr_comm_bubble)
    x$args_inters <- args_inters
    x$group.by <- group.by
    return(x)
  })


setMethod("step0", signature = c(x = "job_cellchatn"),
  function(x){
    step_message("Prepare your data with function `job_cellchatn`.")
  })

setMethod("step1", signature = c(x = "job_cellchatn"),
  function(x){
    step_message("Plot group comparison.")
    p.inters_counts <- funPlot(
      .plot_interactions_across_datasets, x$args_inters$count
    )
    p.inters_counts <- set_lab_legend(
      wrap(p.inters_counts, 14, 7),
      glue::glue("{x@sig} Comparison Number communication network"),
      glue::glue("数量通讯网络|||不同细胞的连线表示潜在的通讯关系，通讯强度数量越多，连线越粗。")
    )
    p.inters_weights <- funPlot(
      .plot_interactions_across_datasets, x$args_inters$weight
    )
    p.inters_weights <- set_lab_legend(
      wrap(p.inters_weights, 14, 7),
      glue::glue("{x@sig} Comparison Strength communication network"),
      glue::glue("强度通讯网络|||不同细胞的连线表示潜在的通讯关系，通讯强度数量越多，连线越粗。")
    )
    x <- methodAdd(x, x$.pre_meth)
    x <- methodAdd(x, "参照 <https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html> 对多数据集单细胞数据进行组间比较分析。")
    x <- plotsAdd(x, p.inters_counts, p.inters_weights)
    return(x)
  })

setMethod("step2", signature = c(x = "job_cellchatn"),
  function(x, cell){
    step_message("Cell as source or target interactions with other cells.")
    allCells <- levels(object(x)@meta[[ x$group.by ]])
    if (!is(cell, "feature_char")) {
      stop('!is(cell, "feature_char").')
    }
    if (length(cell) > 1) {
      stop('length(cell) > 1.')
    }
    if (!any(isThat <- grpl(allCells, cell))) {
      stop('!any(isThat <- grpl(allCells, cell)), can not found.')
    }
    keyCell <- allCells[ isThat ]
    otherCells <- allCells[ !isThat ]
    p.keyCellAsSource <- e(
      CellChat::netVisual_bubble(
        object(x), comparison = c(1, 2),
        sources.use = keyCell, targets.use = otherCells, 
        angle.x = 45
      )
    )
    data.keyCellAsSource <- .get_ggplot_content(p.keyCellAsSource)
    p.keyCellAsSource <- set_lab_legend(
      p.keyCellAsSource,
      glue::glue("{x@sig} {cell} as source comparison"),
      glue::glue("{cell} 作为通讯发送方与其他细胞之间的通讯在组间的比较|||颜色表示细胞通讯概率（Communication Probability）大小，颜色由蓝色逐渐过渡至红色，分别代表由低到高的通讯强度。点的大小表示统计学显著性水平，点越大代表显著性越高；小点表示 p > 0.05，中等点表示 0.01 < p ≤ 0.05，大点表示 p < 0.01。")
    )
    p.keyCellAsTarget <- e(
      CellChat::netVisual_bubble(
        object(x), comparison = c(1, 2),
        sources.use = otherCells, targets.use = keyCell, 
        angle.x = 45
      )
    )
    x <- methodAdd(x, "利用 netVisual_bubble 函数绘制气泡图对各个细胞类型中配体受体介导的相互作用。")
    data.keyCellAsTarget <- .get_ggplot_content(p.keyCellAsTarget)
    p.keyCellAsTarget <- set_lab_legend(
      p.keyCellAsTarget,
      glue::glue("{x@sig} key cell as target comparison"),
      glue::glue("{cell} 作为通讯接收方与其他细胞之间的通讯在组间的比较|||颜色表示细胞通讯概率（Communication Probability）大小，颜色由蓝色逐渐过渡至红色，分别代表由低到高的通讯强度。点的大小表示统计学显著性水平，点越大代表显著性越高；小点表示 p > 0.05，中等点表示 0.01 < p ≤ 0.05，大点表示 p < 0.01。")
    )
    p.allCells_LP_comm_each_group <- x$each_lr_comm$p.lr_comm_bubbles
    p.allCells_LP_comm_each_group <- set_lab_legend(
      p.allCells_LP_comm_each_group,
      glue::glue("{x@sig} {names(p.allCells_LP_comm_each_group)} ligand receptor interactions bubble plot"),
      glue::glue("Group {names(p.allCells_LP_comm_each_group)}: 不同细胞类型之间的配体-受体对相互作用气泡图。|||纵坐标为配体-受体对，横坐标为细胞-细胞相互作用方向，颜色代表交互的可能性，颜色越红代表通讯可能越高，气泡的点大小代表显著性。")
    )
    x <- methodAdd(x, "以 `CellChat::netVisual_bubble` 比较 {snap(cell)} 的组间差异配体与受体相互作用 (⟦mark$blue('P value &lt; 0.05')⟧)。")
    lr_comm <- x$each_lr_comm$t.lr_comm_bubble
    s.com <- glue::glue(
      "{names(lr_comm)} 组包含 {vapply(lr_comm, nrow, integer(1))} 对唯一互作"
    )
    x <- snapAdd(x, "各组细胞的配体、受体通讯统计{aref(p.allCells_LP_comm_each_group)}，{bind(s.com)} (P value &lt; 0.05)。")
    x$keyCell_data <- namel(data.keyCellAsSource, data.keyCellAsTarget)
    snap_compare <- .stat_cellchat_keycell_summary(
      x$keyCell_data, p.keyCellAsSource, p.keyCellAsTarget
    )
    x <- snapAdd(x, "\n\n\n{snap_compare}")
    # x <- snapAdd(x, "其中，{keyCell} 作为互作发送方。")
    x <- plotsAdd(x, p.keyCellAsSource, p.keyCellAsTarget, p.allCells_LP_comm_each_group)
    return(x)
  })

.plot_interactions_across_datasets <- function(nets, 
  wmax, type, nms = names(nets))
{
  if (is.null(nms)) {
    stop('is.null(nms).')
  }
  par(mfrow = c(1, length(nets)), xpd = TRUE)
  type <- switch(type, count = "Number", weight = "Strength")
  lapply(seq_along(nets), 
    function(i) {
      CellChat::netVisual_circle(
        nets[[i]], weight.scale = TRUE,
        label.edge= FALSE, edge.weight.max = wmax[2],
        edge.width.max = 12, title.name = paste0(
          type, " of interactions - ", nms[i]
        )
      )
    })
}

.stat_cellchat_keycell_summary <- function(
  keyCell_data,
  p.source = NULL,
  p.target = NULL,
  top_n = 3L,
  digits = 3L
)
{
  .make_one <- function(dat, slot_name, top_n, digits, p = NULL) {

    if (is.null(dat)) {
      return(NULL)
    }

    dat <- tibble::as_tibble(dat)

    if (nrow(dat) == 0L) {
      return(NULL)
    }

    dat <- dat[!is.na(dat$prob), , drop = FALSE]

    if (nrow(dat) == 0L) {
      return(NULL)
    }

    is_target <- identical(slot_name, "data.keyCellAsTarget")

    key_cell <- if (is_target) {
      as.character(stats::na.omit(unique(dat$target)))[1L]
    } else {
      as.character(stats::na.omit(unique(dat$source)))[1L]
    }

    role_txt <- if (is_target) {
      "受体端"
    } else {
      "发送端"
    }

    partner_col <- if (is_target) "source" else "target"

    n_links <- nrow(dat)
    n_partner <- length(unique(dat[[partner_col]]))
    n_dataset <- length(unique(dat$dataset))
    mean_prob <- mean(dat$prob, na.rm = TRUE)
    max_prob <- max(dat$prob, na.rm = TRUE)
    strong_n <- sum(dat$pval >= 3L, na.rm = TRUE)

    partner_stat <- dplyr::summarise(
      dplyr::group_by(dat, .data[[partner_col]]),
      n_links = dplyr::n(),
      mean_prob = mean(.data$prob, na.rm = TRUE)
    )

    partner_stat <- dplyr::arrange(
      partner_stat,
      dplyr::desc(.data$n_links),
      dplyr::desc(.data$mean_prob)
    )

    top_partner <- utils::head(partner_stat, top_n)

    partner_txt <- paste(
      vapply(
        seq_len(nrow(top_partner)),
        function(i) {
          glue::glue(
            "{top_partner[[1L]][i]}（{top_partner$n_links[i]}条）"
          )
        },
        character(1L)
      ),
      collapse = "、"
    )

    lr_stat <- dplyr::mutate(
      dat,
      lr_pair = paste(.data$ligand, .data$receptor, sep = " - ")
    )

    lr_stat <- dplyr::summarise(
      dplyr::group_by(lr_stat, .data$lr_pair),
      n_links = dplyr::n(),
      mean_prob = mean(.data$prob, na.rm = TRUE)
    )

    lr_stat <- dplyr::arrange(
      lr_stat,
      dplyr::desc(.data$n_links),
      dplyr::desc(.data$mean_prob)
    )

    top_lr <- utils::head(lr_stat, top_n)

    lr_txt <- paste(
      vapply(
        seq_len(nrow(top_lr)),
        function(i) {
          glue::glue(
            "{top_lr$lr_pair[i]}（{top_lr$n_links[i]}次）"
          )
        },
        character(1L)
      ),
      collapse = "、"
    )

    ds_stat <- dplyr::summarise(
      dplyr::group_by(dat, .data$dataset),
      n_links = dplyr::n()
    )

    ds_stat <- dplyr::arrange(ds_stat, dplyr::desc(.data$n_links))

    if (nrow(ds_stat) >= 2L) {
      ds_txt <- glue::glue(
        "{ds_stat$dataset[1L]} 通讯事件更多（{ds_stat$n_links[1L]} vs {ds_stat$n_links[nrow(ds_stat)]}）"
      )
    } else {
      ds_txt <- "仅检测到单一分组结果"
    }

    glue::glue(
      "以 {key_cell} 为核心进行{role_txt}通讯分析{aref(p)}，共识别到 {n_links} 条显著互作，涉及 {n_partner} 类互作细胞，平均通讯概率为 {round(mean_prob, digits)}，最高为 {round(max_prob, digits)}。显著性较强（p < 0.01）的互作共 {strong_n} 条。主要互作对象集中于 {partner_txt}。优势配体-受体轴主要包括 {lr_txt}。组间比较显示，{ds_txt}。"
    )
  }

  txt_target <- .make_one(
    keyCell_data[["data.keyCellAsTarget"]],
    "data.keyCellAsTarget",
    top_n,
    digits, p = p.target
  )

  txt_source <- .make_one(
    keyCell_data[["data.keyCellAsSource"]],
    "data.keyCellAsSource",
    top_n,
    digits, p = p.source
  )

  out <- c(txt_target, txt_source)
  out <- out[!vapply(out, is.null, logical(1L))]

  paste(out, collapse = "\n\n")
}



