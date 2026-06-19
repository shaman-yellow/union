# ==========================================================================
# workflow of stringdb
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_stringdb <- setClass("job_stringdb", 
  contains = c("job"),
  representation = representation(
    object = "ANY",
    params = "list",
    plots = "list",
    tables = "list",
    others = "ANY"),
  prototype = prototype(
    info = c("Tutorial: https://www.bioconductor.org/packages/release/bioc/html/STRINGdb.html"),
    cite = "[@TheStringDataSzklar2021; @CytohubbaIdenChin2014]",
    method = "R package `STEINGdb` used for PPI network construction",
    tag = "ppi:stringdb",
    analysis = "STRINGdb PPI 分析"
    ))

job_stringdb <- function(data)
{
  if (is.character(data)) {
    data <- data.frame(Symbol = rm.no(data))
  }
  data <- dplyr::distinct(data, Symbol)
  data$Symbol <- genes(data$Symbol)
  .job_stringdb(object = data)
}

setGeneric("asjob_stringdb",
  function(x, ...) standardGeneric("asjob_stringdb"))

setMethod("asjob_stringdb", signature = c(x = "job_herb"),
  function(x){
    job_stringdb(data = x$ppi_used)
  })

setMethod("asjob_stringdb", signature = c(x = "character"),
  function(x){
    message("`x` should be gene Symbols.")
    job_stringdb(data.frame(Symbol = x))
  })

setMethod("asjob_stringdb", signature = c(x = "feature"),
  function(x, extra = NULL){
    x <- resolve_feature_snapAdd_onExit("x", x)
    if (!is.null(extra)) {
      x <- c(x, extra)
    }
    x <- job_stringdb(unlist(x))
    return(x)
  })

setMethod("step0", signature = c(x = "job_stringdb"),
  function(x){
    step_message("Prepare your data with function `job_stringdb`.
      "
    )
  })

setMethod("step1", signature = c(x = "job_stringdb"),
  function(x, tops = 30, layout = "spiral", species = 9606, score_threshold = 400,
    network_type = "phy", input_directory = .prefix("stringdb_physical_v12.0", name = "db"),
    version = "12.0", label = FALSE, HLs = NULL, use.anno = TRUE, link_data = "detailed",
    file_anno = .prefix("stringdb_physical_v12.0/9606.protein.physical.links.full.v12.0.txt.gz", name = "db"),
    filter.exp = 0, filter.text = 0, MCC = TRUE, args_spiral = list())
  {
    step_message("Create PPI network.")
    set.seed(x$seed)
    require(ggraph)
    x$network_type <- network_type <- match.arg(network_type, c("physical", "full"))
    if (!dir.exists(input_directory)) {
      dir.create(input_directory)
    }
    if (is.null(x$sdb)) {
      message("Use STRINGdb network type of '", network_type, "'")
      sdb <- new_stringdb(
        score_threshold = score_threshold, species = species, network_type = network_type,
        input_directory = input_directory, version = version, link_data = link_data
      )
      x$sdb <- sdb
    } else {
      sdb <- x$sdb 
    }
    if (is.null(x$res.str)) {
      res.str <- create_interGraph(sdb, data.frame(object(x)), col = "Symbol")
      x$res.str <- res.str
    } else {
      res.str <- x$res.str
    }
    if (is.null(x$graph)) {
      graph <- fast_layout(
        x$res.str$graph, layout = layout, args_spiral = args_spiral
      )
      graph$name <- x$res.str$mapped$Symbol[match(graph$name, x$res.str$mapped$STRING_id)]
      # igraph <- dedup.edges(igraph)
      x$graph <- graph
    } else {
      graph <- x$graph 
    }
    edges <- as_tibble(igraph::as_data_frame(res.str$graph))
    edges <- dplyr::distinct(edges, from, to, .keep_all = TRUE)
    if (species == 9606) {
      message("`use.anno` not available for non hsa.")
      use.anno <- FALSE
    }
    if (use.anno) {
      message("Get PPI annotation from:\n\t", file_anno)
      anno <- ftibble(file_anno)
      edges <- tbmerge(edges[, 1:2], anno, by.x = c("from", "to"), by.y = paste0("protein", 1:2),
        all.x = TRUE, sort = FALSE)
      if (filter.exp || filter.text) {
        edges <- dplyr::filter(edges, experiments >= !!filter.exp, textmining >= !!filter.text)
      }
    }
    edges <- map(edges, "from", res.str$mapped, "STRING_id", "Symbol", rename = FALSE)
    edges <- map(edges, "to", res.str$mapped, "STRING_id", "Symbol", rename = FALSE)
    des_edges <- list("STRINGdb network type:" = match.arg(network_type, c("physical", "full")))
    if (use.anno) {
      des_edges <- c(des_edges,
        list("Filter experiments score:" = paste0("At least score ", filter.exp),
          "Filter textmining score:" = paste0("At least score ", filter.text)
        ))
    }
    edges <- .set_lab(edges, sig(x), "PPI annotation")
    attr(edges, "lich") <- new_lich(des_edges)
    x$edges <- edges
    if (FALSE && network_type == "full") {
      p.ppi <- NULL
    } else {
      p.ppi <- ppiFuns$plot_network_str(graph, label = label)
      h.ppi <- nrow(x$res.str$mapped) %/% 10
      p.ppi <- set_lab_legend(
        wrap(p.ppi, (h.ppi + 1.5) * 1.5, h.ppi * 1.5),
        glue::glue("{x@sig} PPI network"),
        glue::glue("PPI 网络图|||每个节点表示一个蛋白 (基因)，连线表示可能存在的相互作用。")
      )
    }
    ## hub genes
    if (MCC) {
      message("Calculate MCC score.")
    }
    hub_genes <- cal_mcc.str(res.str, "Symbol", FALSE, MCC = MCC)
    graph_mcc <- get_subgraph.mcc(res.str$graph, hub_genes, top = tops)
    x$graph_mcc <- graph_mcc <- fast_layout(graph_mcc, layout = "linear", circular = TRUE)
    x$graph_mcc <- .set_lab(x$graph_mcc, sig(x), "graph MCC layout data")
    snap.mcc <- if (MCC) "(带有 Cytohubba {cite_show('CytohubbaIdenChin2014')} MCC 得分)" else ""
    x$graph_mcc <- setLegend(x$graph_mcc, "PPI {snap.mcc}附表")
    if (!is.null(tops)) {
      feature(x) <- head(hub_genes$Symbol, n = tops)
    } else {
      feature(x) <- hub_genes$Symbol
    }
    p.mcc <- plot_networkFill.str(
      graph_mcc, label = "Symbol", HLs = HLs, netType = network_type
    )
    p.mcc <- .set_lab(wrap(p.mcc), sig(x), paste0("Top", tops, " MCC score"))
    p.mcc <- setLegend(p.mcc, "PPI {snap.mcc}网络图")
    x <- plotsAdd(x, p.ppi, p.mcc)
    x <- tablesAdd(x, hub_genes, mapped = dplyr::relocate(res.str$mapped, Symbol, STRING_id))
    x$tops <- tops
    if (MCC) {
      ex <- glue::glue("以 Cytohubba {cite_show('CytohubbaIdenChin2014')} 的算法在 R 中计算 MCC (Maximal Clique Centrality) 。")
    } else {
      ex <- ""
    }
    nAll <- nrow(x$res.str$mapped)
    nIsolate <- nAll - length(unique(c(x$edges$from, x$edges$to)))
    x <- methodAdd(x, "STRING database 蛋白–蛋白相互作用（PPI）网络分析是一种基于已知与预测相互作用信息构建分子互作网络的方法，其主要目的是从系统层面解析基因或蛋白之间的功能关联关系。通过将候选基因映射至 STRING 数据库 (<https://string-db.org/>) ，构建 PPI 网络并分析其拓扑结构（如节点连接度、聚类系数等），可以识别在网络中处于核心地位的关键蛋白（hub genes）及其参与的功能模块。进一步结合功能富集分析，可揭示这些关键节点在特定生物学过程或疾病机制中的潜在作用，从而为筛选重要调控分子及后续实验验证提供依据。")
    x <- methodAdd(x, "以 R 包 `STEINGdb` ⟦pkgInfo('STRINGdb')⟧ {cite_show('TheStringDataSzklar2021')} 构建 PPI 网络。数据版本为 {version}，互作类型为 {network_type}。置信评分 (confidence score) 阈值为 {score_threshold / 1000} (&gt; {score_threshold / 1000})。{ex}随后，以 R 包 `ggraph` ⟦pkgInfo('ggraph')⟧ 可视化网络。")
    x <- snapAdd(x, "PPI 网络图 {aref(p.ppi)} 共包含 {nAll} 个蛋白 (基因)，存在 {nrow(edges)} 对相互作用，孤立蛋白数量为{nIsolate}。")
    return(x)
  })

setMethod("step2", signature = c(x = "job_stringdb"),
  function(x, n_top = 10L, component = "largest", as_key = TRUE, ...)
  {
    step_message("Screening hub genes by PPI topology.")

    res_cytohubba <- ppiFuns$run_ppi_cytohubba_like(
      edges = x$edges,
      col_from = "from",
      col_to = "to",
      n_top = n_top,
      directed = FALSE,
      component = component
    )

    x$res_cytohubba <- res_cytohubba

    meth_cytohubba <- ppiFuns$get_ppi_topology_method_text(
      n_top = n_top,
      component = component
    )
    x <- methodAdd(x, "{meth_cytohubba}")

    p.venn <- new_venn(
      lst = res_cytohubba$lst_top,
      force_upset = FALSE,
      ...
    )

    text_intersect <- if (length(p.venn$ins) > 0L) {
      less(p.venn$ins, 20L)
    } else {
      "None"
    }

    p.venn <- set_lab_legend(
      p.venn,
      glue::glue("{x@sig} intersection of PPI hub genes"),
      glue::glue(
        "PPI 核心基因交集图|||",
        "该图展示基于 PPI 网络拓扑中心性筛选得到的核心基因交集。分别按照 Degree、Betweenness Centrality ",
        "和 Closeness Centrality 对网络节点进行降序排序，并选取每种算法排名前 {n_top} 的基因进行交集分析。",
        "图中 {length(p.venn$ins)} 个交集基因为：{text_intersect}。"
      )
    )

    x <- plotsAdd(x, p.venn)

    x <- snapAdd(
      x,
      glue::glue(
        "基于 STRING PPI 网络，分别采用 Degree、Betweenness Centrality 和 Closeness Centrality ",
        "筛选排名前 {n_top} 的核心节点，并对三种拓扑算法的候选基因取交集，",
        "最终获得 {length(p.venn$ins)} 个基因{aref(p.venn)}。"
      )
    )
    if (as_key) {
      x <- snapAdd(x, "⟦mark$red('将交集基因 ({text_intersect}) 定义为关键基因')⟧。")
      x$.feature <- as_feature(p.venn$ins, "关键基因")
    }

    return(x)
  })



setMethod("filter", signature = c(x = "job_stringdb"),
  function(x, ref.x, ref.y, lab.x = "Source", lab.y = "Target",
    use = "preferred_name", data = x$graph, level.x = NULL,
    lab.fill = "log2FC",
    ## this top is used for 'from' or 'to'
    top = 10, use.top = c("from", "to"),
    top_in = NULL, keep.ref = if (is.null(top)) TRUE else FALSE,
    keep_extra_link = TRUE, show.mcc = FALSE,
    arrow = TRUE, ...)
  {
    message("Search and filter: ref.x in from, ref.y in to; or, reverse.")
    ref.x <- resolve_feature(ref.x)
    ref.y <- resolve_feature(ref.y)
    use.top <- match.arg(use.top)
    data <- tibble::as_tibble(get_edges()(data), .name_repair = "minimal")
    data <- dplyr::select(data, dplyr::ends_with(use))
    data <- dplyr::rename(data, from = 1, to = 2)
    if (keep.ref) {
      data <- dplyr::filter(data, from %in% c(ref.x, ref.y), to %in% c(ref.x, ref.y))
    } else {
      data <- dplyr::filter(data,
        (from %in% ref.x & to %in% ref.y) | (from %in% ref.y & to %in% ref.x)
      )  
    }
    data <- dplyr::mutate(data, needRev = ifelse(from %in% ref.x, FALSE, TRUE))
    data <- apply(data, 1,
      function(x) {
        if (x[[ "needRev" ]]) {
          x[2:1]
        } else {
          x[1:2]
        }
      })
    data <- tibble::as_tibble(t(data))
    edges <- data <- dplyr::rename(data, from = 1, to = 2)
    nodes <- tibble::tibble(name = unique(c(data$from, data$to)))
    nodes <- dplyr::mutate(nodes,
      type = ifelse(name %in% ref.x, "from", "to")
    )
    data <- dplyr::mutate(data, id = "pseudo")
    data <- dplyr::relocate(data, id)
    p.ppi <- plot_network.pharm(data, edge_width = 1, ax2 = lab.x, ax3 = lab.y,
      ax2.level = level.x, lab.fill = lab.fill, ...)
    p.ppi <- .set_lab(p.ppi, sig(x), "filtered and formated PPI network")
    mcc <- cal_mcc(edges)
    nodes <- map(nodes, "name", mcc, "name", "MCC_score", col = "MCC_score")
    fun_tops <- function() {
      ## tops from ref.x or ref.y
      tops <- dplyr::filter(nodes, type == !!use.top)
      if (!is.null(top_in)) {
        fun <- function() {
          name <- switch(use.top, from = lab.x, to = lab.y)
          if (!is(top_in, "list")) {
            message("`top_in` is not 'list' with names, converted as 'list'.")
            top_in <- nl("Set", list(top_in))
          }
          ## venn plot
          new_venn(lst = c(nl(name, list(tops$name)), top_in))
        }
        p.top_in <- fun()
        p.top_in <- .set_lab(p.top_in, sig(x), "intersection with pre-filter data")
        if (length(p.top_in$ins)) {
          tops <- dplyr::filter(tops, name %in% !!unlist(top_in))
        } else {
          stop("length(p.top_in$ins) == 0, no features in the `top_in`.")
        }
      } else {
        p.top_in <- NULL
      }
      all_edges <- edges
      if (!is.null(top) && !keep.ref) {
        tops <- dplyr::slice_max(tops, MCC_score, n = top)
        edges <- dplyr::filter(edges, !!rlang::sym(use.top) %in% tops$name)
      }
      nodes <- dplyr::slice(nodes, c(which(name %in% ref.x), which(name %in% ref.y)))
      nodes <- dplyr::distinct(nodes, name, .keep_all = TRUE)
      if (keep_extra_link) {
        extras <- dplyr::filter(all_edges, from %in% nodes$name & to %in% nodes$name)
        edges <- dplyr::bind_rows(edges, extras)
      }
      graph <- igraph::graph_from_data_frame(edges, vertices = nodes)
      graph <- fast_layout(graph, layout = "linear", circular = TRUE)
      p.mcc <- plot_networkFill.str(graph, label = "name",
        arrow = arrow, shape = TRUE,
        levels = level.x,
        lab.fill = if (is.null(level.x)) "MCC score" else "Log2(FC)",
        netType = x$network_type, ...
      )
      if (!show.mcc) {
        p.mcc <- p.mcc + guides(fill = "none")
      }
      p.mcc <- .set_lab(p.mcc, sig(x), "Top MCC score")
      colnames(edges) <- c(lab.x, lab.y)
      namel(p.mcc, nodes, edges, p.top_in)
    }
    p.mcc <- fun_tops()
    x$filter_results <- namel(p.ppi, nodes, edges = edges, p.top_in = p.mcc$p.top_in,
      p.mcc = p.mcc$p.mcc, nodes_mcc = p.mcc$nodes, edges_mcc = p.mcc$edges
    )
    return(x)
  })

setMethod("asjob_enrich", signature = c(x = "job_stringdb"),
  function(x, tops = x$tops){
    ids <- head(x@tables$step1$hub_genes$Symbol, tops)
    job_enrich(list(hub_genes = ids), x@tables$step1$hub_genes)
  })

setMethod("vis", signature = c(x = "job_stringdb"),
  function(x, HLs){
    p <- ggraph(x$graph_mcc) +
      geom_edge_arc(
        aes(color = ifelse(node1.Symbol == !!HLs | node2.Symbol == !!HLs,
            paste0(node1.Symbol, " <-> ", node2.Symbol), "...Others"))) +
      geom_node_point(aes(x = x, y = y, fill = ifelse(is.na(MCC_score), 0, MCC_score)),
        size = 12, shape = 21, stroke = .3) +
      geom_node_text(aes(x = x * 1.2, y = y * 1.2, label = Symbol,
          angle = -((-node_angle(x,  y) + 90) %% 180) + 90), size = 4) +
      scale_fill_gradient(low = "lightyellow", high = "red") +
      scale_edge_color_manual(values = c("grey92", color_set())) +
      labs(edge_color = "Link", fill = "MCC score") +
      theme_void()
    p <- wrap(p, 12, 9)
    p <- .set_lab(p, sig(x), "MCC score of PPI top feature")
    data <- get_edges()(x$graph_mcc)
    data <- dplyr::filter(data, node1.Symbol %in% !!HLs | node2.Symbol %in% !!HLs)
    data <- dplyr::select(data, node1.Symbol, node2.Symbol)
    namel(p.mcc = p, data)
  })

new_stringdb <- function(
  score_threshold = 200,
  species = 9606,
  network_type = c("physical", "full"),
  link_data = c("detailed", "full", "combined_only"),
  input_directory = .prefix("stringdb_physical_v12.0", name = "db"),
  version = "12.0")
{
  if (packageVersion("STRINGdb") < "2.22.0") {
    e(STRINGdb::STRINGdb$new(score_threshold = score_threshold,
        species = species, network_type = match.arg(network_type), 
        input_directory = input_directory, version = version
        ))
  } else {
    e(STRINGdb::STRINGdb$new(score_threshold = score_threshold,
        species = species, network_type = match.arg(network_type), 
        link_data = match.arg(link_data),
        input_directory = input_directory, version = version
        ))
  }
}

create_interGraph <- function(sdb, data, col = "name", rm.na = TRUE) {
  cli::cli_alert_info("sdb$map")
  mapped <- sdb$map(data.frame(data), col, removeUnmappedRows = rm.na)
  message()
  cli::cli_alert_info("sdb$get_subnetwork")
  graph <- sdb$get_subnetwork(mapped$STRING_id)
  list(mapped = mapped, graph = graph)
}

fast_layout.str <- function(res.str, sdb, layout = "fr", seed = runif(1, max = 1000000000)) {
  ## get annotation
  target <- file.path(sdb$input_directory, paste0(sdb$species, ".protein.info.v", sdb$version, ".txt"))
  dir.create(dir <- file.path(sdb$input_directory, "temp"), FALSE)
  if (!file.exists(file <- file.path(dir, basename(target)))) {
    file.copy(gfile <- paste0(target, ".gz"), dir)
    R.utils::gunzip(paste0(file, ".gz"))
  }
  anno <- data.table::fread(file)
  ## merge with annotation
  igraph <- add_attr.igraph(res.str$graph, anno, by.y = "#string_protein_id")
  set.seed(seed)
  graph <- fast_layout(igraph, layout)
  attr(graph, "igraph") <- igraph 
  graph
}

cal_pagerank <- function(igraph) {
  res <- igraph::page_rank(igraph)$vector
  data <- data.frame(name = names(res), weight = unname(res))
  igraph <- add_attr.igraph(igraph, data, by.y = "name")
  igraph
}

get_nodes <- function(igraph, from = "vertices") {
  tibble::as_tibble(igraph::as_data_frame(igraph, from))
} 

dedup.edges <- function(igraph){
  add_attr.igraph(igraph)
}

add_attr.igraph <- function(igraph, data, by.x = "name", by.y, dedup.edges = FALSE)
{
  comps <- igraph::as_data_frame(igraph, "both")
  if (!missing(data)) {
    nodes <- merge(comps$vertices, data, by.x = by.x, by.y = by.y, all.x = TRUE)
  } else {
    nodes <- comps$vertices
  }
  if (dedup.edges) {
    edges <- dplyr::distinct(comps$edges)
  } else {
    edges <- comps$edges
  }
  igraph <- igraph::graph_from_data_frame(edges, vertices = nodes)
  igraph
}

output_graph <- function(igraph, file, format = "graphml", toCyDir = TRUE) {
  igraph::write_graph(igraph, file, format = format)
}

plot_network.str <- function(graph, scale.x = 1.1, scale.y = 1.1,
  label.size = 4, sc = 5, ec = 5, 
  arr.len = 2, edge.color = 'grey80', edge.width = .4, label = FALSE)
{
  if (label) {
    layer.nodes <- geom_node_label(aes(label = name), size = label.size)
  } else {
    layer.nodes <- geom_node_point(
      aes(x = x, y = y, color = centrality_degree),
      stroke = .3, alpha = .8, size = 8)
  }
  p <- ggraph(graph) +
    geom_edge_fan(aes(x = x, y = y),
      color = edge.color, width = edge.width) +
    layer.nodes +
    scale_x_continuous(limits = zoRange(graph$x, scale.x)) +
    scale_y_continuous(limits = zoRange(graph$y, scale.y)) +
    theme_minimal() +
    theme(axis.text = element_blank(), axis.title = element_blank())
  if (!label) {
    pal <- c("#7d1339", "#fba8ad")
    p <- p + geom_node_text(aes(label = name), size = 4) +
      ggplot2::scale_color_gradient(low = pal[2], high = pal[1])
  }
  p
} 

plot_networkFill.str <- function(graph, scale.x = 1.1, scale.y = 1.1,
  label.size = 4, node.size = 12, sc = 5, ec = 5, 
  arr.len = 2, edge.color = NULL, edge.width = 1,
  lab.fill = if (is.null(levels)) "MCC score" else "Levels",
  label = "genes", HLs = NULL, arrow = FALSE, shape = FALSE, levels = NULL,
  label.shape = c(from = "from", to = "to"), netType = c("physical", "full"), ...)
{
  dataNodes <- ggraph::get_nodes()(graph)
  if (is.null(levels)) {
    fill <- "MCC_score"
    dataNodes <- dplyr::mutate(dataNodes, Levels = ifelse(is.na(MCC_score), 0, MCC_score))
    pal <- rev(color_set2())
    scale_fill_gradient <- scale_fill_gradient(low = pal[1], high = pal[2])
  } else {
    fill <- "Levels"
    if (!is.data.frame(levels)) {
      stop("!is.data.frame(levels) == FALSE")
    }
    message("The first and second columns of `levels` were used as name and levels.")
    dataNodes <- map(tibble::tibble(dataNodes), "name",
      levels, colnames(levels)[1], colnames(levels)[2], col = "Levels")
    dataNodes <- dplyr::mutate(dataNodes, Levels = ifelse(is.na(Levels), 0, Levels))
    pal <- color_set2()
    scale_fill_gradient <- ggplot2::scale_fill_gradient2(low = pal[2], high = pal[1])
  }
  if (shape) {
    geom_node_point <- geom_node_point(
      data = dataNodes,
      aes(x = x, y = y, fill = !!rlang::sym(fill), shape = type),
      size = node.size, stroke = .3, alpha = .7)
  } else {
    geom_node_point <- geom_node_point(
      data = dataNodes,
      aes(x = x, y = y, fill = !!rlang::sym(fill)),
      size = node.size, shape = 21, stroke = .3, alpha = .7)
  }
  if (is.null(edge.color)) {
    edge.color <- sample(color_set()[1:10], 1)
  }
  netType <- match.arg(netType)
  p <- ggraph(graph) +
    ggraph::geom_edge_arc(
      aes(x = x, y = y, edge_linetype = !!netType),
      start_cap = circle(sc, 'mm'),
      end_cap = circle(ec, 'mm'),
      arrow = if (arrow) arrow(length = unit(arr.len, 'mm')) else NULL,
      color = edge.color, width = edge.width, alpha = .5) +
    geom_node_point +
    geom_node_text(aes(label = !!rlang::sym(label)), size = label.size) +
    scale_fill_gradient +
    scale_x_continuous(limits = zoRange(graph$x, scale.x)) +
    scale_y_continuous(limits = zoRange(graph$y, scale.y)) +
    labs(fill = lab.fill, shape = "Type", edge_linetype = "Interaction") +
    theme_void() +
    theme(plot.margin = margin(r = .05, unit = "npc")) +
    geom_blank()
  if (shape) {
    p <- p + scale_shape_manual(values = c(24, 21, 22, 23), labels = label.shape)
  }
  if (!is.null(HLs)) {
    data <- dplyr::filter(dataNodes, name %in% !!HLs)
    p <- p + geom_point(data = data, aes(x = x, y = y), shape = 21, color = "red", size = 20)
  }
  if (fill == "MCC_score" && all(dataNodes[[ fill ]] %in% c(NA, 0))) {
    p <- p + guides(fill = "none")
  }
  p
} 

cal_mcc.str <- function(res.str, name = "name", rename = TRUE, ...){
  hubs_score <- cal_mcc(res.str$graph, ...)
  hubs_score <- tbmerge(res.str$mapped, hubs_score, by.x = "STRING_id", by.y = "name", all.x = TRUE)
  hubs_score <- dplyr::relocate(hubs_score, !!rlang::sym(name), MCC_score)
  hubs_score <- dplyr::arrange(hubs_score, dplyr::desc(MCC_score))
  if (rename)
    hubs_score <- dplyr::rename(hubs_score, genes = name)
  hubs_score
}

cal_mcc <- function(edges, MCC = TRUE)
{
  if (is(edges, "igraph")) {
    igraph <- edges
    edges <- igraph::as_data_frame(edges, "edges")
  } else if (is(edges, "data.frame")) {
    igraph <- igraph::graph_from_data_frame(edges, FALSE)
  }
  nodes <- unique(unlist(c(edges[, 1], edges[ , 2])))
  if (MCC) {
    maxCliques <- igraph::max_cliques(igraph)
    scores <- vapply(nodes, FUN.VALUE = double(1), USE.NAMES = FALSE,
      function(node) {
        if.contains <- vapply(maxCliques, FUN.VALUE = logical(1), USE.NAMES = FALSE,
          function(clique) {
            members <- attributes(clique)$names
            if (any(members == node)) TRUE else FALSE
          })
        in.cliques <- maxCliques[ if.contains ]
        scores <- vapply(in.cliques, FUN.VALUE = double(1),
          function(clique) {
            num <- length(attributes(clique)$names)
            factorial(num - 1)
          })
        sum(scores)
      })
  } else {
    message("Escape from calculating MCC score.")
    scores <- 0L
  }
  res <- data.frame(name = nodes, MCC_score = scores)
  res <- tibble::as_tibble(dplyr::arrange(res, dplyr::desc(MCC_score)))
  if (exists(".add_internal_job") && MCC) {
    .add_internal_job(
      .job(method = "The MCC score was calculated referring to algorithm of `CytoHubba`",
        cite = "[@CytohubbaIdenChin2014]"))
  }
  res
}

get_subgraph.mcc <- function(igraph, resMcc, top = 10)
{
  tops <- dplyr::arrange(resMcc, dplyr::desc(MCC_score))
  if (!is.null(top)) {
    tops <- head(tops$STRING_id, n = top)
  } else {
    tops <- tops$STRING_id
  }
  data <- igraph::as_data_frame(igraph, "both")
  nodes <- dplyr::filter(data$vertices, name %in% !!tops)
  nodes <- merge(nodes, resMcc, by.x = "name", by.y = "STRING_id", all.x = TRUE)
  edges <- dplyr::filter(data$edges, (from %in% !!tops) & (to %in% !!tops))
  nodes <- dplyr::distinct(nodes, name, .keep_all = TRUE)
  igraph <- igraph::graph_from_data_frame(edges, FALSE, nodes)
  igraph
}

sortDup_edges <- function(edges) {
  edges.sort <- apply(dplyr::select(edges, 1:2), 1,
    function(vec) {
      sort(vec)
    })
  edges.sort <- tibble::as_tibble(data.frame(t(edges.sort)))
  edges.sort <- dplyr::distinct(edges.sort)
  edges.sort
}

getBelong_edges <- function(edges) {
  nodes <- unique(unlist(c(edges[, 1], edges[, 2])))
  edges.rev <- edges
  edges.rev[, 1:2] <- edges.rev[, 2:1]
  links.db <- rbind(edges, edges.rev)
  lst.belong <- split(links.db, unlist(links.db[, 1]))
  lst.belong <- lapply(lst.belong,
    function(data) unlist(data[, 2], use.names = FALSE))
  lst.belong
}

# ==========================================================================
# functions

ppiFuns <- new.env(parent = emptyenv())

ppiFuns$resolve_col <- function(data, col, arg_name = "col")
{
  if (is.numeric(col) && length(col) == 1L) {
    if (col < 1L || col > ncol(data)) {
      stop(sprintf("%s is out of range.", arg_name))
    }

    return(colnames(data)[col])
  }

  if (is.character(col) && length(col) == 1L) {
    if (!col %in% colnames(data)) {
      stop(sprintf("%s column was not found: %s.", arg_name, col))
    }

    return(col)
  }

  stop(sprintf("%s should be one column name or one column index.", arg_name))
}

ppiFuns$prepare_ppi_edges <- function(edges,
  col_from,
  col_to,
  remove_self_loop = TRUE
)
{
  data_edges <- as.data.frame(edges, stringsAsFactors = FALSE)

  col_from <- ppiFuns$resolve_col(data_edges, col_from, "col_from")
  col_to <- ppiFuns$resolve_col(data_edges, col_to, "col_to")

  data_edges <- data.frame(
    from = trimws(as.character(data_edges[[col_from]])),
    to = trimws(as.character(data_edges[[col_to]])),
    stringsAsFactors = FALSE
  )

  data_edges <- data_edges[
    !is.na(data_edges$from) &
      !is.na(data_edges$to) &
      data_edges$from != "" &
      data_edges$to != "",
    ,
    drop = FALSE
  ]

  if (isTRUE(remove_self_loop)) {
    data_edges <- data_edges[
      data_edges$from != data_edges$to,
      ,
      drop = FALSE
    ]
  }

  data_edges
}

ppiFuns$prepare_ppi_graph <- function(edges,
  col_from,
  col_to,
  directed = FALSE,
  remove_self_loop = TRUE,
  simplify_graph = TRUE,
  component = c("all", "largest")
)
{
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop('Package "igraph" is required.')
  }

  component <- match.arg(component)

  data_edges <- ppiFuns$prepare_ppi_edges(
    edges = edges,
    col_from = col_from,
    col_to = col_to,
    remove_self_loop = remove_self_loop
  )

  if (nrow(data_edges) == 0L) {
    stop("No valid PPI edge was found.")
  }

  graph <- igraph::graph_from_data_frame(
    d = data_edges,
    directed = directed
  )

  if (isTRUE(simplify_graph)) {
    graph <- igraph::simplify(
      graph,
      remove.multiple = TRUE,
      remove.loops = remove_self_loop
    )
  }

  if (component == "largest") {
    vec_comp <- igraph::components(graph)$membership
    tab_comp <- table(vec_comp)
    id_comp <- names(tab_comp)[which.max(tab_comp)]

    graph <- igraph::induced_subgraph(
      graph,
      vids = igraph::V(graph)[vec_comp == id_comp]
    )
  }

  graph
}

ppiFuns$run_ppi_topology_score <- function(edges,
  col_from,
  col_to,
  n_top = 10L,
  directed = FALSE,
  normalized = FALSE,
  remove_self_loop = TRUE,
  simplify_graph = TRUE,
  component = c("all", "largest"),
  rank_ties_method = "min"
)
{
  graph <- ppiFuns$prepare_ppi_graph(
    edges = edges,
    col_from = col_from,
    col_to = col_to,
    directed = directed,
    remove_self_loop = remove_self_loop,
    simplify_graph = simplify_graph,
    component = component
  )

  vec_gene <- igraph::V(graph)$name

  vec_degree <- igraph::degree(
    graph = graph,
    v = igraph::V(graph),
    mode = "all",
    loops = FALSE,
    normalized = normalized
  )

  vec_betweenness <- igraph::betweenness(
    graph = graph,
    v = igraph::V(graph),
    directed = directed,
    weights = NULL,
    normalized = normalized,
    cutoff = -1L
  )

  vec_closeness <- igraph::closeness(
    graph = graph,
    vids = igraph::V(graph),
    mode = "all",
    weights = NULL,
    normalized = normalized,
    cutoff = -1L
  )

  data_score <- data.frame(
    gene = vec_gene,
    degree = as.numeric(vec_degree),
    betweenness = as.numeric(vec_betweenness),
    closeness = as.numeric(vec_closeness),
    stringsAsFactors = FALSE
  )

  vec_score_col <- c("degree", "betweenness", "closeness")

  data_score[vec_score_col] <- lapply(
    data_score[vec_score_col],
    function(vec_score) {
      vec_score[!is.finite(vec_score)] <- 0
      vec_score
    }
  )

  data_score$degree_rank <- rank(
    -data_score$degree,
    ties.method = rank_ties_method,
    na.last = "keep"
  )

  data_score$betweenness_rank <- rank(
    -data_score$betweenness,
    ties.method = rank_ties_method,
    na.last = "keep"
  )

  data_score$closeness_rank <- rank(
    -data_score$closeness,
    ties.method = rank_ties_method,
    na.last = "keep"
  )

  data_score$is_top_degree <- data_score$degree_rank <= n_top
  data_score$is_top_betweenness <- data_score$betweenness_rank <= n_top
  data_score$is_top_closeness <- data_score$closeness_rank <= n_top

  data_score$n_top_method <- rowSums(
    data_score[
      c("is_top_degree", "is_top_betweenness", "is_top_closeness")
    ]
  )

  data_score$rank_sum <- data_score$degree_rank +
    data_score$betweenness_rank +
    data_score$closeness_rank

  data_score <- data_score[
    order(data_score$rank_sum, data_score$gene),
    ,
    drop = FALSE
  ]

  rownames(data_score) <- NULL

  data_score
}

ppiFuns$get_top_table <- function(data_score,
  score_col,
  method_name,
  n_top = 10L
)
{
  data_top <- data_score[
    order(-data_score[[score_col]], data_score$gene),
    ,
    drop = FALSE
  ]

  data_top <- data_top[
    seq_len(min(n_top, nrow(data_top))),
    ,
    drop = FALSE
  ]

  data.frame(
    method = method_name,
    rank = seq_len(nrow(data_top)),
    gene = data_top$gene,
    score = data_top[[score_col]],
    stringsAsFactors = FALSE
  )
}

ppiFuns$run_ppi_cytohubba_like <- function(edges,
  col_from,
  col_to,
  n_top = 10L,
  directed = FALSE,
  normalized = FALSE,
  remove_self_loop = TRUE,
  simplify_graph = TRUE,
  component = c("all", "largest"),
  rank_ties_method = "min"
)
{
  data_score <- ppiFuns$run_ppi_topology_score(
    edges = edges,
    col_from = col_from,
    col_to = col_to,
    n_top = n_top,
    directed = directed,
    normalized = normalized,
    remove_self_loop = remove_self_loop,
    simplify_graph = simplify_graph,
    component = component,
    rank_ties_method = rank_ties_method
  )

  lst_method <- list(
    Degree = "degree",
    Betweenness = "betweenness",
    Closeness = "closeness"
  )

  data_top <- do.call(
    rbind,
    lapply(
      names(lst_method),
      function(method_name) {
        ppiFuns$get_top_table(
          data_score = data_score,
          score_col = lst_method[[method_name]],
          method_name = method_name,
          n_top = n_top
        )
      }
    )
  )

  rownames(data_top) <- NULL

  lst_top <- lapply(
    names(lst_method),
    function(method_name) {
      data_top$gene[data_top$method == method_name]
    }
  )

  names(lst_top) <- names(lst_method)

  vec_intersect <- Reduce(intersect, lst_top)

  data_key <- data_score[
    data_score$gene %in% vec_intersect,
    ,
    drop = FALSE
  ]

  data_key <- data_key[
    order(data_key$rank_sum, data_key$gene),
    ,
    drop = FALSE
  ]

  rownames(data_key) <- NULL

  list(
    data_score = data_score,
    data_top = data_top,
    lst_top = lst_top,
    vec_intersect = vec_intersect,
    data_key = data_key
  )
}

ppiFuns$get_ppi_topology_method_text <- function(n_top = 10L,
  component = c("all", "largest")
)
{
  component <- match.arg(component)

  text_component <- if (component == "largest") {
    "为减少离散小网络对中心性排序的影响，本分析进一步提取 PPI 网络中的最大连通子图用于拓扑中心性计算。"
  } else {
    "本分析保留过滤后的全部非游离 PPI 网络节点用于拓扑中心性计算。"
  }

  glue::glue(
    "为识别 PPI 网络中的核心调控节点，本研究基于 STRING 数据库获得的蛋白互作边表构建无向、无权重 PPI 网络，",
    "并在 R 环境中采用 igraph 包计算节点拓扑中心性指标。<<text_component>>",
    "设 PPI 网络为 $G=(V,E)$，其中 $V$ 表示蛋白节点集合，$E$ 表示蛋白互作边集合。",
    "对任意节点 $v \\in V$，度值中心性用于衡量该节点直接连接的邻居数量，定义为：",
    "\n\n",
    "$$\n",
    "Degree(v)=|N(v)|\n",
    "$$\n",
    "\n\n",
    "其中 $N(v)$ 表示与节点 $v$ 直接相连的邻居节点集合。介数中心性用于衡量节点位于其他节点最短路径上的程度，定义为：",
    "\n\n",
    "$$\n",
    "Betweenness(v)=\\sum_{s \\ne v \\ne t}\\frac{\\sigma_{st}(v)}{\\sigma_{st}}\n",
    "$$\n",
    "\n\n",
    "其中 $\\sigma_{st}$ 表示节点 $s$ 与节点 $t$ 之间的最短路径总数，$\\sigma_{st}(v)$ 表示这些最短路径中经过节点 $v$ 的路径数。",
    "接近中心性用于衡量节点到网络中其他节点的整体距离，定义为：",
    "\n\n",
    "$$\n",
    "Closeness(v)=\\frac{1}{\\sum_{u \\in V, u \\ne v}d(v,u)}\n",
    "$$\n",
    "\n\n",
    "其中 $d(v,u)$ 表示节点 $v$ 与节点 $u$ 之间的最短路径长度。",
    "分别按 Degree、Betweenness Centrality 和 Closeness Centrality 对所有节点进行降序排序，",
    "选取每种指标排名前 <<n_top>> 的基因，并取三种指标 Top <<n_top>> 基因的交集作为候选关键基因。",
    "该流程在 R 语言中完成 (主要使用 R 包 `igraph` ⟦pkgInfo('igraph')⟧)，与 Cytoscape cytoHubba 中 Degree、Betweenness 和 Closeness 三类拓扑算法的核心计算逻辑一致，",
    "避免了手动软件操作带来的重复性和效率问题。",
    .open = "<<",
    .close = ">>"
  )
}


ppiFuns$prepare_radial_label_layout <- function(graph,
  center_x = NULL,
  center_y = NULL,
  label_nudge = .01,
  nudge_mode = c("relative", "absolute"),
  keep_upright = TRUE
)
{
  nudge_mode <- match.arg(nudge_mode)

  if (!all(c("x", "y") %in% colnames(graph))) {
    stop('The graph layout should contain columns "x" and "y".')
  }

  data_graph <- graph

  vec_x <- data_graph$x
  vec_y <- data_graph$y

  if (is.null(center_x)) {
    center_x <- mean(range(vec_x, na.rm = TRUE))
  }

  if (is.null(center_y)) {
    center_y <- mean(range(vec_y, na.rm = TRUE))
  }

  vec_dx <- vec_x - center_x
  vec_dy <- vec_y - center_y
  vec_radius <- sqrt(vec_dx^2 + vec_dy^2)

  vec_unit_x <- ifelse(vec_radius > 0, vec_dx / vec_radius, 1)
  vec_unit_y <- ifelse(vec_radius > 0, vec_dy / vec_radius, 0)

  range_x <- range(vec_x, na.rm = TRUE)
  range_y <- range(vec_y, na.rm = TRUE)
  span_xy <- max(diff(range_x), diff(range_y), na.rm = TRUE)

  if (!is.finite(span_xy) || span_xy <= 0) {
    span_xy <- 1
  }

  if (nudge_mode == "relative") {
    label_nudge <- label_nudge * span_xy
  }

  vec_angle <- atan2(vec_dy, vec_dx) * 180 / pi
  vec_left <- vec_angle > 90 | vec_angle < -90

  vec_label_angle <- vec_angle

  if (isTRUE(keep_upright)) {
    vec_label_angle <- ifelse(
      vec_left,
      vec_label_angle + 180,
      vec_label_angle
    )
  }

  vec_label_angle <- ifelse(
    vec_label_angle > 180,
    vec_label_angle - 360,
    vec_label_angle
  )

  vec_label_angle <- ifelse(
    vec_label_angle < -180,
    vec_label_angle + 360,
    vec_label_angle
  )

  data_graph$.label_x <- vec_x + vec_unit_x * label_nudge
  data_graph$.label_y <- vec_y + vec_unit_y * label_nudge
  data_graph$.label_angle <- vec_label_angle
  data_graph$.label_hjust <- ifelse(vec_left & keep_upright, 1, 0)
  data_graph$.label_vjust <- .5

  data_graph
}

ppiFuns$plot_network_str <- function(graph,
  scale.x = 1.1,
  scale.y = 1.1,
  label.size = 3,
  sc = 5,
  ec = 5,
  arr.len = 2,
  edge.color = "grey80",
  edge.width = .4,
  label = FALSE,
  radial_label = TRUE,
  radial_nudge = .01,
  radial_nudge_mode = c("relative", "absolute"),
  radial_center_x = NULL,
  radial_center_y = NULL,
  keep_label_upright = TRUE,
  check_overlap = FALSE,
  centrality_col = "centrality_degree",
  node_size = 8,
  node_alpha = .8,
  pal = c("#7d1339", "#fba8ad")
)
{
  radial_nudge_mode <- match.arg(radial_nudge_mode)

  if (isTRUE(radial_label)) {
    graph <- ppiFuns$prepare_radial_label_layout(
      graph = graph,
      center_x = radial_center_x,
      center_y = radial_center_y,
      label_nudge = radial_nudge,
      nudge_mode = radial_nudge_mode,
      keep_upright = keep_label_upright
    )
  }

  data_label <- as.data.frame(graph)

  if (!centrality_col %in% colnames(data_label)) {
    stop(glue::glue("Column was not found: {centrality_col}."))
  }

  if (isTRUE(label)) {
    layer_nodes <- ggraph::geom_node_label(
      ggplot2::aes(label = name),
      size = label.size
    )
  } else {
    layer_nodes <- ggraph::geom_node_point(
      ggplot2::aes(color = .data[[centrality_col]]),
      stroke = .3,
      alpha = node_alpha,
      size = node_size
    )
  }

  p <- ggraph::ggraph(graph) +
    ggraph::geom_edge_fan(
      color = edge.color,
      width = edge.width
    ) +
    layer_nodes +
    ggplot2::scale_x_continuous(
      limits = zoRange(graph$x, scale.x)
    ) +
    ggplot2::scale_y_continuous(
      limits = zoRange(graph$y, scale.y)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank()
    )

  if (!isTRUE(label)) {
    if (isTRUE(radial_label)) {
      p <- p +
        ggplot2::geom_text(
          data = data_label,
          inherit.aes = FALSE,
          check_overlap = check_overlap,
          size = label.size,
          ggplot2::aes(
            x = .label_x,
            y = .label_y,
            label = name,
            angle = .label_angle,
            hjust = .label_hjust,
            vjust = .label_vjust
          )
        )
    } else {
      p <- p +
        ggraph::geom_node_text(
          ggplot2::aes(label = name),
          size = label.size
        )
    }

    p <- p +
      ggplot2::scale_color_gradient(
        low = pal[2L],
        high = pal[1L],
        name = centrality_col
      )
  }

  p
}

