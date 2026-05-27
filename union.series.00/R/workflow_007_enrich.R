# ==========================================================================
# workflow of enrich
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_enrich <- setClass("job_enrich", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("..."),
    cite = "[@ClusterprofilerWuTi2021]",
    method = "R package `clusterProfiler` used for gene enrichment analysis",
    tag = "enrich:clusterProfiler",
    analysis = "ClusterProfiler 富集分析"
    ))

setGeneric("asjob_enrich",
  function(x, ...) standardGeneric("asjob_enrich"))

setMethod("asjob_enrich", signature = c(x = "feature"),
  function(x, unlist = TRUE, ...){
    names <- names(x)
    x <- resolve_feature_snapAdd_onExit("x", x)
    if (unlist) {
      x <- job_enrich(unlist(x), ...)
    } else {
      names(x) <- names
      x <- job_enrich(x, ...)
    }
    return(x)
  })

job_enrich <- function(ids, annotation, from = "hgnc_symbol", to = "entrezgene_id")
{
  if (!is(ids, "list")) {
    ids <- list(ids = gname(rm.no(ids)))
  } else {
    ids <- lapply(ids, function(x) gname(rm.no(x)))
  }
  if (is.null(names(ids))) {
    stop("is.null(names(ids))")
  }
  if (missing(annotation)) {
    if (from != "hgnc_symbol") {
      stop('from != "hgnc_symbol".')
    }
    if (TRUE) {
      ids <- lapply(ids, gname)
      annotation <- e(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
          keys = unique(unlist(ids)), 
          keytype = "SYMBOL", columns = c("SYMBOL", "ENTREZID")))
      annotation <- annotation[!is.na(annotation[, 2]), ]
      annotation <- dplyr::rename(
        annotation, !!!nl(c(from, to), c("SYMBOL", "ENTREZID"))
      )
    } else {
      mart <- new_biomart()
      annotation <- filter_biomart(
        mart, c(from, "entrezgene_id"), from, unique(unlist(ids))
      )
    }
  }
  maps <- lapply(ids,
    function(id) {
      res <- unique(filter(annotation, !!rlang::sym(from) %in% !!id)[[ to ]])
      res[!is.na(res)]
    })
  en <- .job_enrich(object = maps)
  en@params$raw <- ids
  en@params$annotation <- annotation
  en$from <- from
  en
}

setMethod("step0", signature = c(x = "job_enrich"),
  function(x){
    step_message("Prepare your data with function `job_enrich`.
      "
    )
  })

# Biological Process, Molecular Function, and Cellular Component groups
setMethod("step1", signature = c(x = "job_enrich"),
  function(x, organism = c("hsa", "mmu", "rno"),
    orgDb = switch(organism,
      hsa = "org.Hs.eg.db", mmu = "org.Mm.eg.db", rno = "org.Rn.eg.db"),
    cl = 4, maxShow.kegg = 20, exclude_disease = TRUE,
    maxShow.go = 10, use = c("p.adjust", "pvalue"))
  {
    step_message("Use clusterProfiler for enrichment.")
    cli::cli_alert_info("clusterProfiler::enrichKEGG")
    organism <- match.arg(organism)
    orgDb <- match.arg(orgDb)
    use <- match.arg(use)
    res.kegg <- expect_local_data(
      "tmp", "kegg", multi_enrichKEGG, 
      list(lst.entrez_id = object(x), organism = organism)
    )
    if (exclude_disease) {
      methodAdd_onExit("x", "对于 KEGG 富集分析，富集后根据 Category 去除 Human Diseases 条目，随后根据阈值统计结果。")
      res.kegg <- lapply(res.kegg, 
        function(x) {
          dplyr::filter(x, category != "Human Diseases")
        })
    }
    fun <- function(sets) {
      lapply(sets,
        function(set) {
          from_ids <- x@params$annotation$entrezgene_id
          to_names <- x@params$annotation[[ x$from ]]
          to_names[ match(set, from_ids) ]
        })
    }
    res.kegg <- lapply(res.kegg, mutate, geneName_list = fun(geneID_list))
    p.kegg <- vis_enrich.kegg(res.kegg, maxShow = maxShow.kegg, use = use)
    use.p <- attr(p.kegg, "use.p")
    p.kegg <- set_lab_legend(
      p.kegg,
      glue::glue("{x@sig} {names(p.kegg)} KEGG enrichment"),
      glue::glue("{.enNames(p.kegg)} KEGG 富集分析气泡图|||横坐标为 GeneRatio (目标基因在基因集中的比例与目标基因的总数的比值)；纵坐标代表富集到通路的名称；点的大小代表富集到的属于该通路的基因的数量；颜色代表-log10({use.p})，值越大代表富集到的通路越显著。")
    )
    res.kegg <- set_lab_legend(
      res.kegg,
      glue::glue("{x@sig} {names(res.kegg)} KEGG enrichment data"),
      glue::glue("为 {.enNames(res.kegg)} KEGG 富集分析统计表。")
    )
    cli::cli_alert_info("clusterProfiler::enrichGO")
    res.go <- expect_local_data(
      "tmp", "go", multi_enrichGO, list(lst.entrez_id = object(x), orgDb = orgDb, cl = cl),
      ignore = "cl"
    )
    p.go <- vis_enrich.go(res.go, maxShow = maxShow.go, use = use.p)
    p.go <- set_lab_legend(
      p.go,
      glue::glue("{x@sig} {names(p.go)} GO enrichment"),
      glue::glue("为 {.enNames(p.go)} GO 富集分析气泡图|||横坐标为 GeneRatio (目标基因在基因集中的比例与目标基因的总数的比值)；纵坐标代表富集到通路的名称；点的大小代表富集到的属于该通路的基因的数量；颜色代表-log10({use.p})，值越大代表富集到的通路越显著。GO 富集囊括了三个子集 Cellular Component (CC), the Molecular Function (MF) and the Biological Process (BP)。")
    )
    res.go <- lapply(res.go,
      function(data) {
        if (all(vapply(data, is.data.frame, logical(1)))) {
          data <- as_tibble(data.table::rbindlist(data, idcol = TRUE))
          data <- dplyr::mutate(data, geneName_list = fun(geneID_list))
          dplyr::relocate(data, ont = .id)
        }
      })
    res.go <- set_lab_legend(
      res.go,
      glue::glue("{x@sig} {names(res.go)} GO enrichment data"),
      glue::glue(" {.enNames(res.go)} GO 富集分析统计表。")
    )
    if (length(res.kegg) == 1) {

      x <- snapAdd(
        x, "KEGG 一共富集到 {.stat_table_by_pvalue(res.kegg[[1]], 10, use.p = use.p, enumeration = FALSE)}"
      )
      x <- snapAdd(
        x, "GO 一共富集到 {.stat_table_by_pvalue(res.go[[1]], 5, 'ont', use.p = use.p, enumeration = FALSE)}"
      )
    }
    x <- tablesAdd(x, res.kegg, res.go)
    x <- plotsAdd(x, p.kegg, p.go)
    x@params$check_go <- check_enrichGO(res.go)
    x$organism <- organism
    x <- methodAdd(x, "以 R 包 `clusterProfiler` ⟦pkgInfo('clusterProfiler')⟧进行 KEGG 和 GO 富集分析。以 {use.p} 表示显著水平 (且按 {use.p} 对结果排序)。富集筛选条件为 {use.p} &lt; 0.05。")
    return(x)
  })


.stat_table_by_pvalue <- function(data, n = 5, 
  split = NULL, use.p = "p.adjust", colName = "Description", target = "通路", by = "富集到",
  needSum = TRUE, enumeration = TRUE)
{
  if (is.null(split)) {
    data <- dplyr::arrange(data, !!rlang::sym(use.p))
    data <- dplyr::filter(data, !!rlang::sym(use.p) < .05)
    paths <- head(data[[colName]], n = n)
    if (enumeration) {
      glue::glue(
        "{nrow(data)} 个{target}，按 {use.p} 值从低到高排序的前 {length(paths)} 个{target}分别为：{atrans(paths)}。"
      )
    } else {
      glue::glue("{nrow(data)} 个{target}。")
    }
  } else {
    data <- dplyr::filter(data, !!rlang::sym(use.p) < .05)
    ele <- split(data, data[[split]])
    ele <- vapply(
      ele, .stat_table_by_pvalue, character(1), 
      n = n, use.p = use.p, target = target, by = by, 
      colName = colName, enumeration = enumeration
    )
    if (enumeration) {
      ele.snap <- bind(glue::glue("在 {names(ele)} {by} {ele}"), co = "\n\n")
    } else {
      ele.snap <- bind(glue::glue("在 {names(ele)} {by} {ele}"), co = "\n")
    }
    if (needSum) {
      glue::glue("{nrow(data)} 个{target}。{ele.snap}")
    } else {
      glue::glue("{ele.snap}")
    }
  }
}

.enNames <- function(x) {
  if (length(x) == 1 && names(x) == "ids") {
    return("")
  } else {
    return(names(x))
  }
}

setMethod("step2", signature = c(x = "job_enrich"),
  function(x, pathways, which.lst = 1, species = x$organism,
    name = paste0("pathview", gs(Sys.time(), " |:", "_")),
    search = "pathview",
    external = FALSE, gene.level = NULL, gene.level.name = "hgnc_symbol")
  {
    stop("deprecated, use 'asjob_pathview' instead.")
    return(x)
  })

setMethod("res", signature = c(x = "job_enrich", ref = "character"),
  function(x, ref = c("id", "des", "cate", "sub", "p", "adj"),
    which = 1, key = 1, from = c("kegg", "go"))
  {
    type <- match.arg(ref)
    type <- switch(
      type, id = "ID", des = "Description", cate = "category", 
      sub = "subcategory", p = "pvalue", adj = "p.adjust"
    )
    from <- match.arg(from)
    data <- x@tables$step1[[ paste0("res.", from) ]][[ key ]]
    data[[ type ]][ which ]
  })

setMethod("asjob_enrich", signature = c(x = "job_seurat"),
  function(x, exclude.pattern = NULL, exclude.use = NULL,
    ignore.case = TRUE, marker.list = x@params$contrasts, geneType = "hgnc_symbol")
  {
    if (is.null(marker.list)) {
      if (x@step < 5) {
        stop("x@step < 5")
      }
      data <- x@tables$step5[[ "all_markers" ]]
      if (!is.null(exclude.pattern)) {
        exclude.cluster <- dplyr::filter(object(x)@meta.data,
          grepl(!!exclude.pattern, !!rlang::sym(exclude.use), ignore.case))$seurat_clusters
        exclude.cluster <- unique(exclude.cluster)
        message("Exclude clasters:\n  ", paste0(exclude.cluster, collapse = ", "))
        data <- dplyr::filter(data, !cluster %in% exclude.cluster)
      }
    } else {
      data <- marker.list
    }
    data <- dplyr::mutate(data, gene = gs(gene, "\\.[0-9]*$", ""))
    mart <- new_biomart()
    anno <- filter_biomart(mart, general_attrs(), geneType, unique(data$gene))
    data <- dplyr::filter(data, gene %in% anno[[ !!geneType ]])
    if (is.null(data$contrast)) {
      ids <- split(data$gene, data$cluster)
    } else {
      ids <- split(data$gene, data$contrast)
    }
    ids <- lst_clear0(ids)
    job_enrich(ids, anno)
  })

setMethod("focus", signature = c(x = "job_enrich"),
  function(x, symbols, data = x@tables$step1$res.kegg[[1]])
  {
    if (x@step < 1L) {
      stop("x@step < 1L")
    }
    isThat <- vapply(data$geneName_list, FUN.VALUE = logical(1),
      function(genes) {
        any(genes %in% symbols)
      })
    data <- dplyr::filter(data, !!isThat)
    data
  })

multi_enrichKEGG <- function(lst.entrez_id, organism = 'hsa')
{
  res <- pbapply::pblapply(lst.entrez_id,
    function(ids) {
      res.kegg <- clusterProfiler::enrichKEGG(ids, organism = organism)
      res.path <- tibble::as_tibble(res.kegg@result)
      res.path <- dplyr::mutate(res.path, geneID_list = lapply(strsplit(geneID, "/"), as.integer))
      res.path
    })
  res
}

multi_enrichGO <- function(lst.entrez_id, orgDb = 'org.Hs.eg.db', cl = NULL)
{
  res <- pbapply::pblapply(lst.entrez_id, cl = cl,
    function(ids) {
      onts <- c("BP", "CC", "MF")
      res <- sapply(onts, simplify = FALSE,
        function(ont) {
          res.go <- try(clusterProfiler::enrichGO(ids, orgDb, ont = ont), TRUE)
          if (inherits(res.go, "try-error")) {
            return("try-error of enrichment")
          }
          res.res <- try(res.go@result, TRUE)
          if (inherits(res.res, "try-error")) {
            value <- "try-error of enrichment"
            attr(value, "data") <- res.res
            return(value)
          }
          res.path <- tibble::as_tibble(res.res)
          res.path <- dplyr::mutate(res.path, geneID_list = lapply(strsplit(geneID, "/"), as.integer))
          res.path
        })
    })
}

check_enrichGO <- function(res.go) {
  isthat <- lapply(res.go,
    function(res) {
      !vapply(res, FUN.VALUE = logical(1), is.character)
    })
  isthat
}

vis_enrich.kegg <- function(lst, cutoff = .05, maxShow = 10,
  use = c("p.adjust", "pvalue"), least = 3L, sankey = TRUE)
{
  use <- match.arg(use)
  use.p <- use
  res <- lapply(lst,
    function(data) {
      data <- dplyr::filter(raw <- data, !!rlang::sym(use) < !!cutoff)
      if (!nrow(data) | nrow(data) < least) {
        message("\n", "Too few of the results (", nrow(data), ")")
        if (use == "p.adjust") {
          message("Switch to use `pvalue`.")
          use <- "pvalue"
          use.p <<- use
          data <- dplyr::filter(raw, !!rlang::sym(use) < !!cutoff)
        }
      }
      data <- dplyr::arrange(data, !!rlang::sym(use))
      data <- head(data, n = maxShow)
      data <- dplyr::mutate(data, GeneRatio = as_double.ratioCh(GeneRatio))
      if (sankey) {
        p <- wrap(
          .plot_kegg_sankey(data, use, n_pathway = maxShow), 8, 3 * (maxShow / 7)
        )
      } else {
        p <- .plot_kegg(data, use)
        p <- wrap(p, 7, 3 * (maxShow / 10))
      }
      p <- setLegend(p, "KEGG 富集图展示了以 {use} 排序，前 {maxShow} 的富集通路。")
      p
    })
  attr(res, "use.p") <- use.p
  res
}

# for external use
plot_kegg <- function(data, cutoff = .05, maxShow = 10,
  use = c("p.adjust", "pvalue"), 
  use.log = FALSE, pattern = NULL, pals = c("grey90", "darkred"), ...)
{
  use <- match.arg(use)
  data <- .format_enrich(data, use, cutoff, maxShow, pattern)
  if (use.log) {
    data <- dplyr::mutate(data, logUse = -log10(!!rlang::sym(use)))
    use <- glue::glue("-log10({use})")
    data <- dplyr::rename(data, !!!setNames(list("logUse"), use))
    pals <- rev(pals)
  }
  p <- .plot_kegg(data, use, pals = pals, ...)
  firstlines <- vapply(
    strsplit(data$Description, "\n"), function(x) x[1], character(1)
  )
  p <- wrap_scale(
    p, max(nchar(firstlines)) / 4, nrow(data)
  )
  p <- .set_lab(p, "KEGG-enrichment")
  p
}

.plot_kegg <- function(data, use, ratio = "GeneRatio", count = "Count", order.by = ratio,
  order.desc = FALSE, pals = c("grey90", "darkred"), theme = geom_blank())
{
  p <- ggplot(data) +
    geom_point(
      aes(x = reordern(Description, !!rlang::sym(order.by), decreasing = order.desc),
        y = !!rlang::sym(ratio), size = !!rlang::sym(count), fill = !!rlang::sym(use)),
      shape = 21, stroke = 0, color = "transparent"
    ) +
    scale_fill_gradient(high = pals[1], low = pals[2]) +
    scale_size(range = c(4, 6)) +
    labs(x = "", y = "Hits Ratio") +
    guides(size = guide_legend(override.aes = list(color = "grey70", stroke = 1))) +
    coord_flip() +
    ylim(zoRange(data[[ ratio ]], 1.3)) +
    rstyle("theme") +
    theme
  p
}

.plot_kegg_sankey <- function(data, use = "pvalue",
  ratio = "GeneRatio", count = "Count", gene = "geneName_list",
  order.by = NULL, order.desc = FALSE,
  n_pathway = 10L, n_gene = Inf,
  gene_min_pathway = 1L,
  gene_value_mode = c("count", "sqrt", "equal"),
  point_pals = c("grey92", "#b22222"),
  path_pals = c("#ce7e73", "#f5a596", "#f8d5ce", "#38546d",
    "#cba3b2", "#fbbe85", "#8390ca", "#c0e0db", "#f7a7a6"),
  x_ratio_left = 0.05,
  x_ratio_right = 0.24,
  x_path_min = 0.33,
  x_path_max = 0.55,
  x_path_anchor = 0.56,
  x_gene_anchor = 0.86,
  x_gene_label = 0.89,
  x_gene_block_min = 0.86,
  x_gene_block_max = 0.99,
  x_gene_text = 0.868,
  gene_block_fill = "grey80",
  gene_block_alt_fill = "grey91",
  gene_block_alpha = 0.92,
  gene_block_colour = "white",
  gene_block_linewidth = 0.18,
  gene_text_colour = "black",
  path_label_wrap = 28L,
  path_label_size = 2.7,
  gene_label_size = 2.2,
  flow_alpha = 0.23,
  flow_width = 0.04,
  path_alpha = 0.85,
  point_size_range = c(2.2, 5.2),
  show_gene_node = TRUE,
  gene_node_width = 0.018,
  gene_node_fill = "grey94",
  gene_node_colour = "white",
  gene_node_alpha = 0.85,
  show_bottom_label = TRUE,
  show_ratio_box = TRUE,
  ratio_box_pad_x = 0.02,
  ratio_box_pad_y = 0.12,
  ratio_box_colour = "grey45",
  ratio_box_linewidth = 0.3,
  y_axis = -0.20,
  y_tick_label = -0.50,
  y_axis_title = -2,
  y_relation_title = -2,
  y_bottom_pad = 0.35,
  legend_position = "left",
  theme = ggplot2::geom_blank())
{
  if (!requireNamespace("ggalluvial", quietly = TRUE)) {
    stop("Package `ggalluvial` is required. Please install it with install.packages('ggalluvial').",
      call. = FALSE)
  }

  gene_value_mode <- match.arg(gene_value_mode)

  .parse_ratio <- function(vec_x)
  {
    if (is.numeric(vec_x)) {
      return(as.numeric(vec_x))
    }

    vec_chr <- as.character(vec_x)
    vec_out <- suppressWarnings(as.numeric(vec_chr))
    idx_fraction <- grepl("/", vec_chr, fixed = TRUE)

    if (any(idx_fraction)) {
      lst_fraction <- strsplit(vec_chr[idx_fraction], "/", fixed = TRUE)

      vec_out[idx_fraction] <- vapply(lst_fraction, function(vec_item) {
        vec_num <- suppressWarnings(as.numeric(vec_item))

        if (length(vec_num) != 2L || any(is.na(vec_num)) || vec_num[2L] == 0) {
          return(NA_real_)
        }

        vec_num[1L] / vec_num[2L]
      }, numeric(1L))
    }

    vec_out
  }

  .safe_neg_log10 <- function(vec_x)
  {
    vec_x <- suppressWarnings(as.numeric(vec_x))
    vec_x[is.na(vec_x) | vec_x <= 0] <- NA_real_

    -log10(vec_x)
  }

  .resolve_gene_list <- function(vec_gene)
  {
    if (is.list(vec_gene)) {
      return(lapply(vec_gene, function(vec_x) {
        vec_x <- trimws(as.character(vec_x))
        unique(vec_x[nzchar(vec_x)])
      }))
    }

    lapply(as.character(vec_gene), function(x) {
      vec_x <- unlist(strsplit(x, "[/;,|]", perl = TRUE))
      vec_x <- trimws(vec_x)
      unique(vec_x[nzchar(vec_x)])
    })
  }

  .get_path_color <- function(n)
  {
    if (length(path_pals) >= n) {
      return(path_pals[seq_len(n)])
    }

    grDevices::colorRampPalette(path_pals)(n)
  }

  .scale_ratio_x <- function(vec_x, vec_range)
  {
    n_span <- diff(vec_range)

    if (!is.finite(n_span) || n_span == 0) {
      return(rep(mean(c(x_ratio_left, x_ratio_right)), length(vec_x)))
    }

    x_ratio_left + (vec_x - vec_range[1L]) / n_span * (x_ratio_right - x_ratio_left)
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

  data_path <- as.data.frame(data)
  vec_required <- c("Description", ratio, count, use, gene)

  if (!all(vec_required %in% names(data_path))) {
    stop("Required columns were not found in `data`.", call. = FALSE)
  }

  data_path$ratio_value <- .parse_ratio(data_path[[ratio]])
  data_path$count_value <- suppressWarnings(as.numeric(data_path[[count]]))

  if (use %in% c("pvalue", "p.adjust", "qvalue")) {
    data_path$fill_value <- .safe_neg_log10(data_path[[use]])
    fill_title <- glue::glue("-log10({use})")
  } else {
    data_path$fill_value <- suppressWarnings(as.numeric(data_path[[use]]))
    fill_title <- use
  }

  data_path <- data_path[
    is.finite(data_path$ratio_value) &
      is.finite(data_path$count_value) &
      is.finite(data_path$fill_value),
    ,
    drop = FALSE
  ]

  if (!is.null(order.by)) {
    if (!order.by %in% names(data_path)) {
      stop(glue::glue("Column `{order.by}` was not found."), call. = FALSE)
    }

    data_path$order_value <- .parse_ratio(data_path[[order.by]])
    data_path <- data_path[
      order(data_path$order_value, decreasing = order.desc, na.last = TRUE),
      ,
      drop = FALSE
    ]
  }

  if (is.finite(n_pathway)) {
    data_path <- data_path[seq_len(min(n_pathway, nrow(data_path))), , drop = FALSE]
  }

  if (!nrow(data_path)) {
    stop("No valid KEGG terms were available for plotting.", call. = FALSE)
  }

  vec_path_order <- data_path$Description
  data_path$path_fill <- .get_path_color(nrow(data_path))
  lst_gene <- .resolve_gene_list(data_path[[gene]])

  data_link <- do.call(rbind, lapply(seq_len(nrow(data_path)), function(i) {
    vec_gene <- lst_gene[[i]]

    if (!length(vec_gene)) {
      return(NULL)
    }

    data.frame(
      pathway = data_path$Description[i],
      gene = vec_gene,
      path_index = i,
      value = 1,
      stringsAsFactors = FALSE
    )
  }))

  if (is.null(data_link) || !nrow(data_link)) {
    stop("No valid genes were found in the gene list column.", call. = FALSE)
  }

  data_gene_n <- stats::aggregate(
    pathway ~ gene,
    data = data_link,
    FUN = function(vec_x) length(unique(vec_x))
  )
  names(data_gene_n)[2L] <- "n_pathway"

  data_gene_center <- stats::aggregate(
    path_index ~ gene,
    data = data_link,
    FUN = mean
  )
  names(data_gene_center)[2L] <- "path_center"

  data_gene <- merge(data_gene_n, data_gene_center, by = "gene", all = TRUE)
  data_gene <- data_gene[data_gene$n_pathway >= gene_min_pathway, , drop = FALSE]
  data_gene <- data_gene[
    order(data_gene$path_center, -data_gene$n_pathway, data_gene$gene),
    ,
    drop = FALSE
  ]

  if (is.finite(n_gene)) {
    data_gene <- data_gene[seq_len(min(n_gene, nrow(data_gene))), , drop = FALSE]
  }

  if (!nrow(data_gene)) {
    stop("No valid genes remained after filtering.", call. = FALSE)
  }

  data_link <- data_link[data_link$gene %in% data_gene$gene, , drop = FALSE]

  data_gene_degree <- stats::aggregate(
    pathway ~ gene,
    data = data_link,
    FUN = function(vec_x) length(unique(vec_x))
  )
  names(data_gene_degree)[2L] <- "gene_degree"

  data_link$gene_degree <- data_gene_degree$gene_degree[
    match(data_link$gene, data_gene_degree$gene)
  ]

  if (gene_value_mode == "equal") {
    data_link$value <- 1 / data_link$gene_degree
  } else if (gene_value_mode == "sqrt") {
    data_link$value <- 1 / sqrt(data_link$gene_degree)
  } else {
    data_link$value <- 1
  }

  data_path <- data_path[data_path$Description %in% data_link$pathway, , drop = FALSE]
  data_path <- data_path[
    match(vec_path_order[vec_path_order %in% data_path$Description], data_path$Description),
    ,
    drop = FALSE
  ]

  data_path_height <- stats::aggregate(
    value ~ pathway,
    data = data_link,
    FUN = sum
  )

  data_path$height <- data_path_height$value[
    match(data_path$Description, data_path_height$pathway)
  ]
  data_path$height[is.na(data_path$height)] <- 0

  data_gene_height <- stats::aggregate(
    value ~ gene,
    data = data_link,
    FUN = sum
  )

  data_gene$height <- data_gene_height$value[match(data_gene$gene, data_gene_height$gene)]
  data_gene <- data_gene[!is.na(data_gene$height) & data_gene$height > 0, , drop = FALSE]

  data_path_layout <- .make_stack_layout(data_path$Description, data_path$height)
  data_gene_layout <- .make_stack_layout(data_gene$gene, data_gene$height)

  data_path$ymin <- data_path_layout$ymin
  data_path$ymax <- data_path_layout$ymax
  data_path$y <- data_path_layout$y

  data_gene$ymin <- data_gene_layout$ymin
  data_gene$ymax <- data_gene_layout$ymax
  data_gene$y <- data_gene_layout$y

  data_gene$gene_fill <- rep(
    c(gene_block_fill, gene_block_alt_fill),
    length.out = nrow(data_gene)
  )

  data_link$path_fill <- data_path$path_fill[match(data_link$pathway, data_path$Description)]
  data_link$alluvium <- paste(data_link$pathway, data_link$gene, sep = "___")

  vec_stratum_level <- c(rev(data_path$Description), rev(data_gene$gene))

  data_alluvium <- rbind(
    data.frame(
      alluvium = data_link$alluvium,
      x = x_path_anchor,
      stratum = data_link$pathway,
      value = data_link$value,
      path_fill = data_link$path_fill,
      stringsAsFactors = FALSE
    ),
    data.frame(
      alluvium = data_link$alluvium,
      x = x_gene_anchor,
      stratum = data_link$gene,
      value = data_link$value,
      path_fill = data_link$path_fill,
      stringsAsFactors = FALSE
    )
  )

  data_alluvium$stratum <- factor(data_alluvium$stratum, levels = vec_stratum_level)

  vec_ratio_range <- range(data_path$ratio_value, na.rm = TRUE)
  vec_ratio_break <- pretty(vec_ratio_range, n = 4L)
  vec_ratio_break <- vec_ratio_break[
    vec_ratio_break >= vec_ratio_range[1L] &
      vec_ratio_break <= vec_ratio_range[2L]
  ]

  if (!length(vec_ratio_break)) {
    vec_ratio_break <- vec_ratio_range
  }

  data_path$x_ratio <- .scale_ratio_x(data_path$ratio_value, vec_ratio_range)
  vec_ratio_break_x <- .scale_ratio_x(vec_ratio_break, vec_ratio_range)

  n_y_max <- max(c(data_path$ymax, data_gene$ymax), na.rm = TRUE)

  data_ratio_box <- data.frame(
    xmin = x_ratio_left - ratio_box_pad_x,
    xmax = x_ratio_right + ratio_box_pad_x,
    ymin = min(data_path$ymin, na.rm = TRUE) - ratio_box_pad_y,
    ymax = max(data_path$ymax, na.rm = TRUE) + ratio_box_pad_y
  )

  data_ratio_axis <- data.frame(
    x = vec_ratio_break_x,
    label = signif(vec_ratio_break, 3L)
  )

  y_lower <- min(y_axis, y_tick_label, y_axis_title, y_relation_title) - y_bottom_pad

  p <- ggplot2::ggplot()

  if (isTRUE(show_ratio_box)) {
    p <- p + ggplot2::geom_rect(
      data = data_ratio_box,
      ggplot2::aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax
      ),
      inherit.aes = FALSE,
      fill = NA,
      colour = ratio_box_colour,
      linewidth = ratio_box_linewidth
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
      data = data_path,
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
      data = data_path,
      ggplot2::aes(
        x = x_path_min + 0.012,
        y = y,
        label = stringr::str_wrap(Description, path_label_wrap)
      ),
      hjust = 0,
      vjust = 0.5,
      lineheight = 0.82,
      size = path_label_size,
      colour = "black",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_point(
      data = data_path,
      ggplot2::aes(
        x = x_ratio,
        y = y,
        size = count_value,
        colour = fill_value
      ),
      inherit.aes = FALSE
    )

  if (isTRUE(show_gene_node)) {
    p <- p + ggplot2::geom_rect(
      data = data_gene,
      ggplot2::aes(
        xmin = x_gene_anchor - gene_node_width / 2,
        xmax = x_gene_anchor + gene_node_width / 2,
        ymin = ymin,
        ymax = ymax
      ),
      inherit.aes = FALSE,
      fill = gene_node_fill,
      colour = gene_node_colour,
      alpha = gene_node_alpha,
      linewidth = 0.18
    )
  }

  p <- p +
    ggplot2::geom_rect(
      data = data_gene,
      ggplot2::aes(
        xmin = x_gene_block_min,
        xmax = x_gene_block_max,
        ymin = ymin,
        ymax = ymax,
        fill = gene_fill
      ),
      alpha = gene_block_alpha,
      colour = gene_block_colour,
      linewidth = gene_block_linewidth,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = data_gene,
      ggplot2::aes(
        x = x_gene_text,
        y = y,
        label = gene
      ),
      hjust = 0,
      vjust = 0.5,
      size = gene_label_size,
      colour = gene_text_colour,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_colour_gradient(
      low = point_pals[1L],
      high = point_pals[2L],
      name = fill_title
    ) +
    ggplot2::scale_size(
      range = point_size_range,
      name = count
    ) +
    ggplot2::scale_x_continuous(
      breaks = NULL,
      labels = NULL,
      limits = c(0.02, max(0.99, x_gene_block_max + 0.015)),
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
          x = x_ratio_left,
          xend = x_ratio_right,
          y = y_axis,
          yend = y_axis
        ),
        inherit.aes = FALSE,
        linewidth = 0.25,
        colour = "grey45"
      ) +
      ggplot2::geom_segment(
        data = data_ratio_axis,
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
        data = data_ratio_axis,
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
        x = mean(c(x_ratio_left, x_ratio_right)),
        y = y_axis_title,
        label = ratio,
        size = 3.2
      ) +
      ggplot2::annotate(
        "text",
        x = mean(c(x_path_min, x_gene_label)),
        y = y_relation_title,
        label = "Gene-Pathway relationship",
        size = 3.0,
        fontface = "bold"
      )
  }

  p
}


reordern <- function(x, by, ...) {
  if (!is.numeric(by)) {
    if (is.character(by)) {
      by <- as.integer(as.factor(by))
    } else if (is.factor(by)) {
      by <- as.integer(by)
    } else {
      stop('`by` should be either numeric or factor or character.')
    }
  }
  reorder(x, by, ...)
}

vis_enrich.go <- function(lst, cutoff = .05, maxShow = 10,
  use = c("p.adjust", "pvalue"), least = 3L)
{
  use <- match.arg(use)
  fun <- function(data) {
    data <- lapply(data,
      function(data) {
        if (is.character(data))
          return()
        data <- dplyr::filter(data, !!rlang::sym(use) < cutoff)
        data <- dplyr::arrange(data, !!rlang::sym(use))
        data <- head(data, n = maxShow)
        data
      })
    data <- data.table::rbindlist(data, idcol = TRUE)
    if (!nrow(data)) {
      return()
    }
    data <- dplyr::mutate(
      data, GeneRatio = as_double.ratioCh(GeneRatio),
      stringr::str_wrap(Description, width = 30)
    )
    p <- .plot_go_use_ggplot(data, use)
    p <- wrap(p, 7.5, 1 + nrow(data) * .2)
    p <- setLegend(p, "GO 富集图展示了基因集在 GO 的 BP (Biological Process), MF (Molecular Function), CC (Cellular Component) 组中的富集结果 (以 {use} 排序，各自展示前 {maxShow} 的富集通路) 。")
    p
  }
  res <- lapply(lst,
    function(x) {
      try(fun(x), silent = TRUE)
    })
  res
}

# this function is for external use. 
plot_go <- function(data, cutoff = .05, maxShow = 10,
  use = c("p.adjust", "pvalue"), facet = ".id", pattern = NULL)
{
  use <- match.arg(use)
  data <- lapply(split(data, data[[ facet ]]),
    function(data) {
      .format_enrich(data, use, cutoff, maxShow, pattern)
    })
  data <- frbind(data)
  p <- .plot_go_use_ggplot(data, use, facet)
  p <- wrap(p, 7.5, 1 + nrow(data) * .2)
  p <- .set_lab(p, "GO-enrichment")
  p <- setLegend(p, "GO 富集图展示了基因集在 GO 的 BP (Biological Process), MF (Molecular Function), CC (Cellular Component) 组中的富集结果 (以 {use} 排序，各自展示前 {maxShow} 的富集通路) 。")
  p
}

plot_go_polor <- function(data,
  pr = c(inner_blank = 70, layer1 = 10, layer2 = 30, layer3 = 30),
  space. = 1 / 30)
{
  space <- sum(pr) * space.
  data <- dplyr::select(data, ONTOLOGY, ID, pvalue, Count, GeneRatio)
  data <- dplyr::group_by(data, ONTOLOGY)
  data <- dplyr::arrange(data, dplyr::desc(Count), .by_group = TRUE)
  data <- dplyr::ungroup(data)
  data <- dplyr::mutate(
    data, xmax = cumsum(Count), xmin = xmax - Count + sum(Count) / 500,
    xmid = xmin + Count / 2,
    angle_pi = -(xmid / max(xmax) * 2 * pi) - pi / 2,
    angle = angle_pi / (2 * pi) * 360,
    angle = angle + as.integer(cos(angle_pi) < 0) * 180,
    logP = -log10(pvalue),
    logG = -log10(GeneRatio),
    norm_logG = logG / max(logG)
  )
  pals <- c(BP ="#8FC93A", CC = "#1982C4", MF = "#FF595E")
  fp <- function(a, scale = 1) {
    if (a >= length(pr) + 1) {
      stop('a >= length(pr)')
    }
    a <- a - 1L
    if (a == -1L) {
      a <- length(pr) - 1
    }
    sum(pr[seq_len(a)]) + pr[a + 1L] * scale
  }
  p <- ggplot(data, aes(xmin = xmin, xmax = xmax)) +
    # p-value
    geom_rect(aes(ymin = fp(1), ymax = fp(2) - space, fill = logP)) +
    scale_fill_gradient(low = "lightyellow", high = "darkred") +
    labs(fill = "-Log10(Pvalue)") +
    ggnewscale::new_scale_fill() +
    # Ontology
    geom_rect(aes(ymin = fp(2), ymax = fp(3)- space, fill = ONTOLOGY)) +
    # other
    geom_rect(aes(ymin = fp(3), ymax = fp(4, norm_logG) - space, fill = ONTOLOGY)) +
    scale_fill_manual(values = pals) +
    labs(fill = "Ontology") +
    guides(fill = guide_legend(order = 1)) +
    # text
    geom_text(
      aes(x = xmid, y = fp(3, 0), label = ID, angle = angle,
        hjust = as.integer(cos(angle_pi) > 0)), 
      size = 4.5
      ) +
    coord_polar() +
    lims(y = c(0, fp(0)), x = c(0, sum(data$Count))) +
    theme_void() +
    theme(
      plot.margin = unit(rep(-.5, 4), "in"), legend.position = c(.5, .5)
      ) +
    geom_blank()
  p <- .set_lab(p, "GO-enrichment-coord-polar")
  p <- setLegend(p, "GO 富集图展示了基因集在 GO 的 BP (Biological Process), MF (Molecular Function), CC (Cellular Component) 组中的富集结果。")
  return(p)
}

.format_enrich <- function(data, use, cutoff, maxShow, pattern = NULL)
{
  data <- dplyr::filter(data, !!rlang::sym(use) < cutoff)
  data <- dplyr::arrange(data, !!rlang::sym(use))
  less <- head(data, n = maxShow)
  if (!is.null(pattern)) {
    extra <- dplyr::filter(data, grpl(Description, pattern, TRUE))
    if (nrow(extra)) {
      less <- dplyr::bind_rows(less, extra)
      less <- dplyr::distinct(less)
    }
  }
  if (!is.null(data$GeneRatio) && !is.double(data$GeneRatio)) {
    less <- dplyr::mutate(less, GeneRatio = as_double.ratioCh(GeneRatio))
    less <- dplyr::arrange(less, GeneRatio)
  }
  less
}

.plot_go_use_ggplot <- function(data, use = "p.adjust", facet = ".id",
  x = "Count", y = "Description", wrap_width = 42L,
  size_range = c(2.5, 6.5), point_stroke = 0.2,
  point_colour = "grey35", strip_text_colour = "white",
  theme = ggplot2::geom_blank())
{
  if (!requireNamespace("ggnewscale", quietly = TRUE)) {
    stop("Package `ggnewscale` is required.", call. = FALSE)
  }

  if (!requireNamespace("ggh4x", quietly = TRUE)) {
    stop("Package `ggh4x` is required.", call. = FALSE)
  }

  if (!requireNamespace("ggtext", quietly = TRUE)) {
    stop("Package `ggtext` is required.", call. = FALSE)
  }

  .lighten_color <- function(col, factor = 0.75)
  {
    mat_rgb <- grDevices::col2rgb(col) / 255
    mat_new <- 1 - (1 - mat_rgb) * (1 - factor)

    grDevices::rgb(mat_new[1L], mat_new[2L], mat_new[3L])
  }

  .build_palette <- function(vec_id)
  {
    vec_dark <- c(
      BP = "#eb6d82",
      CC = "#b8d8ad",
      MF = "#526fb3",
      KEGG = "#3ea88a"
    )

    vec_missing <- setdiff(vec_id, names(vec_dark))

    if (length(vec_missing)) {
      vec_extra <- grDevices::hcl.colors(length(vec_missing), palette = "Dark 3")
      names(vec_extra) <- vec_missing
      vec_dark <- c(vec_dark, vec_extra)
    }

    vec_dark <- vec_dark[vec_id]
    vec_light <- stats::setNames(
      vapply(vec_dark, .lighten_color, character(1L)),
      vec_id
    )

    list(dark = vec_dark, light = vec_light)
  }

  data_plot <- as.data.frame(data)
  vec_need <- c(facet, x, y, use, "Count")

  if (!all(vec_need %in% colnames(data_plot))) {
    stop("Missing required columns in `data`.", call. = FALSE)
  }

  data_plot[[facet]] <- as.character(data_plot[[facet]])
  data_plot[[x]] <- suppressWarnings(as.numeric(data_plot[[x]]))
  data_plot[[use]] <- suppressWarnings(as.numeric(data_plot[[use]]))
  data_plot[["Count"]] <- suppressWarnings(as.numeric(data_plot[["Count"]]))
  data_plot[[y]] <- stringr::str_wrap(as.character(data_plot[[y]]), width = wrap_width)

  data_plot <- data_plot[
    !is.na(data_plot[[facet]]) &
      is.finite(data_plot[[x]]) &
      is.finite(data_plot[[use]]) &
      is.finite(data_plot[["Count"]]),
    ,
    drop = FALSE
  ]

  if (!nrow(data_plot)) {
    stop("No valid rows remained after filtering.", call. = FALSE)
  }

  vec_preferred <- c("BP", "CC", "MF", "KEGG")
  vec_id <- unique(data_plot[[facet]])
  vec_id <- unique(c(intersect(vec_preferred, vec_id), setdiff(vec_id, vec_preferred)))

  data_plot$.facet <- factor(data_plot[[facet]], levels = vec_id)

  lst_level <- lapply(vec_id, function(id) {
    data_sub <- data_plot[data_plot[[facet]] == id, , drop = FALSE]
    data_sub <- data_sub[order(data_sub[[x]], decreasing = FALSE), , drop = FALSE]

    paste(id, data_sub[[y]], sep = "|||")
  })

  vec_term_level <- unlist(lst_level, use.names = FALSE)

  data_plot$.term <- paste(data_plot[[facet]], data_plot[[y]], sep = "|||")
  data_plot$.term <- factor(data_plot$.term, levels = vec_term_level)

  lst_pal <- .build_palette(vec_id)
  vec_dark <- lst_pal$dark
  vec_light <- lst_pal$light

  data_term_label <- unique(
    data_plot[, c(".term", facet, y), drop = FALSE]
  )

  data_term_label$.term_chr <- as.character(data_term_label$.term)
  data_term_label$.label_html <- sprintf(
    "<span style='color:%s;'>%s</span>",
    vec_dark[data_term_label[[facet]]],
    gsub("\n", "<br>", data_term_label[[y]], fixed = TRUE)
  )

  vec_term_label <- stats::setNames(
    data_term_label$.label_html,
    data_term_label$.term_chr
  )

  obj_theme <- if (exists("rstyle", mode = "function")) {
    rstyle("theme")
  } else {
    ggplot2::theme_bw()
  }

  p <- ggplot2::ggplot()

  for (i in seq_along(vec_id)) {
    id <- vec_id[i]
    data_sub <- data_plot[data_plot[[facet]] == id, , drop = FALSE]

    p <- p + ggplot2::geom_point(
      data = data_sub,
      ggplot2::aes(
        x = .data[[x]],
        y = .term,
        size = Count,
        fill = .data[[use]]
      ),
      shape = 21,
      stroke = point_stroke,
      colour = point_colour,
      alpha = 0.95,
      inherit.aes = FALSE
    )

    p <- p + ggplot2::scale_fill_gradient(
      low = vec_light[id],
      high = vec_dark[id],
      trans = "reverse",
      name = paste0(id, " ", use),
      guide = ggplot2::guide_colourbar(
        order = i,
        title.position = "top",
        barheight = grid::unit(16, "mm"),
        barwidth = grid::unit(4, "mm")
      )
    )

    if (i < length(vec_id)) {
      p <- p + ggnewscale::new_scale_fill()
    }
  }

  p <- p +
    ggplot2::scale_size(
      range = size_range,
      name = "Count",
      guide = ggplot2::guide_legend(
        order = length(vec_id) + 1L,
        override.aes = list(fill = "grey70", colour = "grey40", stroke = 0.2)
      )
    ) +
    ggh4x::facet_grid2(
      rows = ggplot2::vars(.facet),
      scales = "free_y",
      space = "free_y",
      switch = "y",
      strip = ggh4x::strip_themed(
        background_y = ggh4x::elem_list_rect(
          fill = unname(vec_dark),
          colour = unname(vec_dark)
        ),
        text_y = ggh4x::elem_list_text(
          colour = strip_text_colour,
          face = "bold",
          angle = 0,
          size = 9
        )
      )
    ) +
    ggplot2::scale_y_discrete(
      position = "right",
      labels = function(x) {
        vec_term_label[as.character(x)]
      }
    ) +
    ggplot2::labs(x = x, y = NULL) +
    obj_theme +
    ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      axis.text.y.left = ggplot2::element_blank(),
      axis.ticks.y.left = ggplot2::element_blank(),
      axis.text.y.right = ggtext::element_markdown(
        size = 7.5,
        hjust = 0,
        lineheight = 0.95
      ),
      axis.ticks.y.right = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.placement = "outside",
      legend.position = "right",
      legend.box = "vertical",
      legend.title = ggplot2::element_text(size = 8, face = "bold"),
      legend.text = ggplot2::element_text(size = 7),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.margin = ggplot2::margin(8, 18, 8, 8)
    ) +
    theme

  p
}


setMethod("map", signature = c(x = "job_enrich", ref = "job_enrich"),
  function(x, ref, use = c("kegg", "go"), key = 1, cutoff = .05, use.cutoff = c("p.adjust", "pvalue"))
  {
    message("Find intersection pathways across two 'job_enrich'")
    use <- match.arg(use)
    use.cutoff <- match.arg(use.cutoff)
    fun_extract <- function(x) {
      dplyr::filter(x@tables$step1[[ paste0("res.", use) ]][[ key ]],
        !!rlang::sym(use.cutoff) < cutoff)
    }
    lst <- lapply(list(x, ref), fun_extract)
    data <- dplyr::filter(lst[[1]], ID %in% !!lst[[2]]$ID)
    data <- .set_lab(data, sig(x), "pathways intersection")
    x$intersect_paths <- data
    return(x)
  })

map_gene <- function(data, col,
  from = "ENTREZID", to = "SYMBOL", split = ",",
  from_bm = "entrezgene_id", to_bm = "hgnc_symbol", try_bm = FALSE, force_bm = FALSE,
  get = to, OrgDb = org.Hs.eg.db::org.Hs.eg.db, gname = TRUE)
{
  if (gname && is.character(data[[ col ]])) {
    data[[ get ]] <- gname(data[[ col ]])
  } else {
    data[[ get ]] <- data[[ col ]]
  }
  data <- dplyr::relocate(data, !!rlang::sym(get), .after = !!rlang::sym(col))
  ids <- data[[ get ]]
  mixSets <- NULL
  if (!is.null(split) && is.character(ids)) {
    mixSets <- strsplit(ids, split)
    maybeMergedIds <- ids
    numMix <- sum(lengths(mixSets) > 1)
    message(glue::glue("Try split id colomn, number of {numMix} merged id rows found."))
    ids <- unlist(mixSets)
    mixSets <- data.frame(
      id = rep(maybeMergedIds, lengths(mixSets)),
      splits = unlist(mixSets)
    )
  }
  if (any(duplicated(ids))) {
    ids <- unique(ids)
  }
  backup <- data[[ get ]]
  funStat <- function(which) {
    misFreq <- sum(is.na(data[[ get ]])) / nrow(data)
    message(glue::glue("Missing '{which}': {misFreq} of data."))
    return(misFreq)
  }
  if (!force_bm) {
    annotation <- e(AnnotationDbi::select(OrgDb, keys = ids, 
        keytype = from, columns = c(from, to)))
    annotation <- annotation[!is.na(annotation[, 2]), ]
    if (!is.null(mixSets)) {
      annotation <- map(
        annotation, from, mixSets, "splits", "id", col = from
      )
    }
    data <- map(data, col, annotation, from, to, col = get)
    misFreq <- funStat(to)
  }
  if (force_bm || (try_bm && misFreq > .3)) {
    message(glue::glue("Try mapping via biomaRt."))
    mart <- new_biomart()
    annotation <- filter_biomart(
      mart, c(from_bm, to_bm), from_bm, ids
    )
    if (!is.null(mixSets)) {
      annotation <- map(
        annotation, from_bm, mixSets, "splits", "id", col = from_bm
      )
    }
    data[[ get ]] <- backup
    data <- map(data, get, annotation, from_bm, to_bm, col = get)
    funStat(to)
  }
  return(data)
}

get_genes.keggPath <- function(name) {
  if (!is(name, "character")) {
    stop("is(name, 'character')")
  }
  lst <- e(KEGGREST::keggGet(name))
  x <- strx(lst[[1]]$GENE, "^[A-Za-z][^;]+")
  x[ !is.na(x) ]
}

as_double.ratioCh <- function(ch) {
  values <- stringr::str_extract_all(ch, "[0-9]{1,}")
  vapply(values, FUN.VALUE = double(1),
    function(values) {
      values <- as.double(values)
      values[ 1 ] / values[ 2 ]
    })
}
