# ==========================================================================
# workflow of regNet
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_regNet <- setClass("job_regNet", 
  contains = c("job"),
  prototype = prototype(
    pg = "regNet",
    info = c(
      targetScanHuman = "https://www.targetscan.org/cgi-bin/targetscan/data_download.vert80.cgi",
      mirdb = "https://mirdb.org/download.html",
      encori = "https://rnasysu.com/encori/tutorialAPI.php",
      miRNetR = "https://github.com/xia-lab/miRNetR/blob/master/vignettes/miRNetR.Rmd",
      tarbase = "https://dianalab.e-ce.uth.gr/tarbasev9/downloads",
      LncBase = "https://diana.e-ce.uth.gr/lncbasev3/home",
      npinter5 = "http://bigdata.ibp.ac.cn/npinter5/download/"
    ),
    cite = "",
    method = "",
    tag = "regNet",
    analysis = "分子调控网络分析"
    ))

setGeneric("asjob_regNet",
  function(x, ...) standardGeneric("asjob_regNet"))

setMethod_traceable("asjob_regNet", signature = c(x = "feature"),
  function(x){
    fea <- resolve_feature_snapAdd_onExit("x", x)
    x <- .job_regNet(object = fea)
    return(x)
  })

setMethod_traceable("step0", signature = c(x = "job_regNet"),
  function(x){
    step_message("Prepare your data with function `job_regNet`.")
  })

setMethod_traceable("step1", signature = c(x = "job_regNet"),
  function(x, recode = NULL, tsh = TRUE, tbs = TRUE, 
    miRNetR = FALSE, mdb = FALSE, miRNetR_cache = NULL)
  {
    step_message("Get miRNA data.")
    targets <- object(x)
    if (!is.null(recode)) {
      targets <- unique(c(targets, names(recode)))
      recode <- as.list(recode)
    }
    fun_recode <- function(data, col, recode) {
      if (!is.null(recode)) {
        data[[ col ]] <- dplyr::recode(
          data[[ col ]], !!!recode, .default = data[[ col ]]
        )
      }
      data
    }
    x$all_miRNA <- list()
    if (tsh) {
      targetScanHuman <- ftibble(get_url_data(
          "Conserved_Family_Info.txt",
          # "Predicted_Targets_Info.default_predictions.txt",
          "https://www.targetscan.org/vert_80/vert_80_data_download/Conserved_Family_Info.txt.zip",
          # "https://www.targetscan.org/vert_80/vert_80_data_download/Predicted_Targets_Info.default_predictions.txt.zip",
          "targetScanHuman"
          ))
      colnames(targetScanHuman) <- formal_name(colnames(targetScanHuman))
      targetScanHuman <- dplyr::filter(targetScanHuman, Gene_Symbol %in% targets)
      which_not_in_data(targetScanHuman, "Gene_Symbol", targets)
      targetScanHuman <- fun_recode(
        targetScanHuman, "Gene_Symbol", recode
      )
      targetScanHuman <- dplyr::rename(
        targetScanHuman, mirna_name = miR_Family, gene_name = Gene_Symbol
      )
      targetScanHuman <- reframe_split(
        targetScanHuman, "mirna_name", "/"
      )
      targetScanHuman <- dplyr::mutate(
        targetScanHuman, mirna_name = ifelse(
          grpl(mirna_name, "^miR"), 
          mirna_name, paste0("miR-", mirna_name)
        ),
        mirna_name = paste0("hsa-", mirna_name)
      )
      targetScanHuman <- dplyr::relocate(
        targetScanHuman, mirna_name, gene_name
      )
      targetScanHuman <- set_lab_legend(
        targetScanHuman,
        glue::glue("{x@sig} miRNA-RNA data from targetScanHuman"),
        glue::glue("从 targetScanHuman 获取的 miRNA-RNA 数据。")
      )
      x$all_miRNA$targetScanHuman <- targetScanHuman
      x <- methodAdd(
        x, "TargetScanHuman (<https://www.targetscan.org/cgi-bin/targetscan/data_download.vert80.cgi>)，"
      )
    }
    if (tbs) {
      tarbase <- ftibble(get_url_data(
          "Homo_sapiens_TarBase-v9.tsv.gz",
          "https://dianalab.e-ce.uth.gr/tarbasev9/data/Homo_sapiens_TarBase-v9.tsv.gz",
          "tarbase", fun_decompress = NULL
          ))
      tarbase <- dplyr::filter(tarbase, gene_name %in% targets)
      which_not_in_data(tarbase, "gene_name", targets)
      tarbase <- fun_recode(tarbase, "gene_name", recode)
      tarbase <- dplyr::relocate(tarbase, mirna_name, gene_name)
      tarbase <- set_lab_legend(
        tarbase,
        glue::glue("{x@sig} miRNA-RNA data from TarBase"),
        glue::glue("从 Tarbase 获取的 miRNA-RNA 数据。")
      )
      x$all_miRNA$tarbase <- tarbase
      x <- methodAdd(
        x, "TarBase (<https://dianalab.e-ce.uth.gr/tarbasev9/downloads>)，"
      )
    }
    if (FALSE && mdb) {
      mirdb <- ftibble(get_url_data(
          "miRDB_v6.0_prediction_result.txt.gz",
          "https://mirdb.org/download/miRDB_v6.0_prediction_result.txt.gz",
          "mirdb", fun_decompress = NULL
          ))
      x$all_miRNA$mirdb <- mirdb
    }
    if (miRNetR) {
      if (!is.null(miRNetR_cache)) {
        miRNetR <- ftibble(miRNetR_cache)
      } else {
        miRNetR <- callr::r(
          map_genes_in_miRNetR,
          list(genes = targets, wd = .prefix("miRNetR", "db"),
            extra = "--no-check-certificate"), show = TRUE
          )$mir.res
        miRNetR <- tibble::as_tibble(miRNetR)
      }
      which_not_in_data(miRNetR, "Target", targets)
      miRNetR <- fun_recode(miRNetR, "Target", recode)
      miRNetR <- dplyr::relocate(miRNetR, mirna_name = ID, gene_name = Target)
      miRNetR <- dplyr::mutate(
        miRNetR, mirna_name = gs(mirna_name, "hsa-mir", "hsa-miR")
      )
      miRNetR <- set_lab_legend(
        miRNetR,
        glue::glue("{x@sig} miRNA-RNA data from miRNet"),
        glue::glue("从 miRNet 获取的 miRNA-RNA 数据。")
      )
      x$all_miRNA$miRNetR <- miRNetR
      x <- methodAdd(x, "miRnet (<https://www.mirnet.ca/>)，")
    }
    x <- methodAdd(x, "以上数据库用于检索以基因集 (mRNA) 为靶点 miRNA 数据。")
    return(x)
  })

setMethod_traceable("step2", signature = c(x = "job_regNet"),
  function(x, use = "all"){
    step_message("miRNA intersection.")
    sets <- x$all_miRNA
    if (!identical(use, "all")) {
      sets <- sets[ use ]
    }
    ins <- .merge_list_by_cols(sets, by = c("mirna_name", "gene_name"))
    ins <- dplyr::distinct(ins, mirna_name, gene_name)
    message(glue::glue("Got data: {try_snap(ins, 'gene_name', 'mirna_name')}"))
    x$ins_mirna <- ins
    ins.snap <- try_snap(ins, "gene_name", "mirna_name")
    if (length(sets) > 1) {
      x <- snapAdd(x, "将数据库 {bind(names(sets))} 预测或记录的 mRNA 的上游 miRNA，二者 mRNA-miRNA 关系对取交集，共得到 {nrow(ins)} 对调控关系【{ins.snap}】。")
      x <- methodAdd(x, "取以上数据库的交集，确定高可信度的miRNA-mRNA调控配对。")
    } else {
      x <- snapAdd(x, "将数据库 {bind(names(sets))} 预测或记录的 mRNA 的上游 miRNA 建立 mRNA-miRNA 关系对，共得到 {nrow(ins)} 对调控关系【{ins.snap}】。")
    }
    return(x)
  })

setMethod_traceable("step3", signature = c(x = "job_regNet"),
  function(x, enc = TRUE, npi = TRUE, mi = unique(x$ins_mirna$mirna_name)){
    step_message("Got lncRNA.")
    x$all_lncRNA <- list()
    if (enc) {
      encori <- get_encori_miRNA_lncRNA(mi)
      encori <- dplyr::relocate(encori, mirna_name = miRNAname, lncrna_name = geneName)
      encori <- set_lab_legend(
        encori,
        glue::glue("{x@sig} miRNA-lncRNA data from ENCORI"),
        glue::glue("从 ENCORI 获取的 miRNA-lncRNA 数据。")
      )
      x$all_lncRNA$encori <- encori
      x <- methodAdd(x, "ENCORI (<https://rnasysu.com/encori/>)，")
    }
    if (npi) {
      npinter <- ftibble(get_url_data(
          "miRNA_interaction.txt.gz",
          "http://bigdata.ibp.ac.cn/npinter5/download/file/miRNA_interaction.txt.gz",
          "npinter_miRNA_interaction", fun_decompress = NULL
          ))
      npinter <- dplyr::relocate(npinter, mirna_name = V5, lncrna_name = V2)
      if (!grpl(npinter[1, 1, drop = TRUE], "miR")) {
        stop('!grpl(npinter[1, 1, drop = TRUE], "miR")')
      }
      npinter <- dplyr::filter(npinter, mirna_name %in% mi)
      which_not_in_data(npinter, "mirna_name", mi)
      npinter <- set_lab_legend(
        npinter,
        glue::glue("{x@sig} miRNA-lncRNA data from NPinter"),
        glue::glue("从 NPinter 获取的 miRNA-lncRNA 数据。")
      )
      x$all_lncRNA$npinter <- npinter
      x <- methodAdd(x, "NPInter (<http://bigdata.ibp.ac.cn/npinter5>)，")
    }
    x <- methodAdd(x, "数据库用于获取与上述 miRNA 靶向结合的 lncRNA。")
    return(x)
  })

setMethod_traceable("step4", signature = c(x = "job_regNet"),
  function(x, use = "all"){
    step_message("lncRNA intersection")
    sets <- x$all_lncRNA
    if (!identical(use, "all")) {
      sets <- sets[ use ]
    }
    ins <- .merge_list_by_cols(sets, by = c("lncrna_name", "mirna_name"))
    ins <- dplyr::distinct(ins, lncrna_name, mirna_name)
    message(glue::glue("Got data: {try_snap(ins, 'mirna_name', 'lncrna_name')}"))
    ins.snap <- try_snap(ins, "mirna_name", "lncrna_name")
    x <- snapAdd(x, "将数据库 {bind(names(sets))} 预测或记录的 miRNA 的上游 lncRNA，二者 miRNA-lncRNA 关系对取交集，共得到 {nrow(ins)} 对调控关系【{ins.snap}】。")
    x <- methodAdd(x, "取两个数据库结果的交集筛选出高可信度的 lncRNA-miRNA 调控配对。")
    x$ins_lncrna <- ins
    return(x)
  })

setMethod_traceable("step5", signature = c(x = "job_regNet"),
  function(x, layout = "fr"){
    step_message("Network.")
    funName <- function(x) {
      setNames(x, c("from", "to"))
    }
    edges <- rbind(funName(x$ins_mirna), funName(x$ins_lncrna))
    nodes <- list(
      mRNA = x$ins_mirna$gene_name, miRNA = x$ins_mirna$mirna_name,
      lncRNA = x$ins_lncrna$lncrna_name
    )
    nodes <- as_df.lst(lapply(nodes, unique), "type", "name")[, 2:1]
    graph <- igraph::graph_from_data_frame(
      edges, directed = TRUE, vertices = nodes
    )
    set.seed(x$seed)
    layout <- ggraph::create_layout(graph, layout = layout)
    require(ggraph)
    p.regNet <- ggraph(layout) + 
      geom_edge_link(
        edge_width = .5, color = "grey80", show.legend = FALSE,
        end_cap = ggraph::circle(7, 'mm'),
        arrow = arrow(length = unit(1, 'mm'))) + 
      geom_node_point(aes(color = type, size = type, shape = type)) + 
      ggrepel::geom_label_repel(aes(x = x, y = y, label = name), size = 3) +
      guides(
        size = "none", shape = "none",
        color = guide_legend(override.aes = list(size = 4))
      ) +
      scale_size_manual(values = c(mRNA = 10, miRNA = 6, lncRNA = 6)) +
      scale_shape_manual(
        values = c(mRNA = 16, miRNA = 17, lncRNA = 18)
      ) +
      labs(color = "Type") +
      theme_void()
    p.regNet <- set_lab_legend(
      p.regNet,
      glue::glue("{x@sig} expression regulation networking"),
      glue::glue("lncRNA-miRNA-mRNA 表达网络分析|||图中的节点表示对应 RNA 类型，边代表相互作用。")
    )
    x <- snapAdd(x, "构建 lncRNA-miRNA-mRNA 表达网络{aref(p.regNet)}，如图所示，共 {nrow(nodes)} 个节点，{nrow(edges)} 个边。")
    x <- methodAdd(x, "以 R 包 `ggraph` ⟦pkgInfo('ggraph')⟧ 与 `ggplot2` ⟦pkgInfo('ggplot2')⟧ 对 lncRNA-miRNA-mRNA 转录后调控网络进行整合，可视化多层级分子调控网络的整体结构。")
    x <- plotsAdd(x, p.regNet)
    return(x)
  })

setMethod_traceable("step6", signature = c(x = "job_regNet"),
  function(x, trrust = TRUE, chipbase = FALSE,
    layout = "fr", filter_chipbase = TRUE, num_quantile = .9, n_min_support = 3L)
  {
    targets <- object(x)
    if (chipbase) {
      file_chipbase <- get_url_data(
        "hg38_network.bed.gz",
        "https://rnasysu.com/chipbase3/data/download/network/hg38_network.bed.gz",
        "chipbase", fun_decompress = NULL
      )
      fun_read <- function(...) {
        .shFilter_read_table_by_id(
          file_chipbase, targets, "gene_symbol", "\t"
        )
      }
      chipbase <- expect_local_data(
        "tmp", "chipbase_network", fun_read, list(targets)
      )
      chipbase <- dplyr::filter(chipbase, protein_tf_type == "tf")
      which_not_in_data(chipbase, "gene_symbol", targets)
      chipbase <- dplyr::rename(chipbase, TF = protein, Target = gene_symbol)
      if (nrow(chipbase)) {
        if (filter_chipbase) {
          chipbase <- .filter_chipbase_tf(
            chipbase, num_quantile = num_quantile, n_min_support = n_min_support
          )
        }
        chipbase <- dplyr::relocate(chipbase, TF, Target)
        chipbase <- dplyr::distinct(chipbase, TF, Target, .keep_all = TRUE)
      }
      x$all_tf$chipbase <- chipbase
      x <- methodAdd(x, "\n\n基于 ChIPBase v3.0 数据库 <https://rnasysu.com/chipbase3/index.php> 获取转录因子与靶基因之间的潜在调控关系。该数据库整合了大量 ChIP-seq 数据，可用于系统分析转录因子在基因启动子或增强子区域的结合情况，并提供转录调控关系及表达相关性信息。通过将候选基因映射至 ChIPBase v3.0，可识别其潜在上游转录因子，并进一步构建 TF–target 调控网络，从而挖掘关键转录调控轴及其在相关生物学过程中的潜在作用机制。")
      meth <- .description_filter_chipbase(num_quantile, n_min_support)
      x <- methodAdd(x, "{meth}")
    }
    if (trrust) {
      trrust <- ftibble(get_url_data(
          "trrust_rawdata.human.tsv",
          "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv",
          "trrust", fun_decompress = NULL
          ))
      trrust <- dplyr::rename(
        trrust, TF = V1, Target = V2, Regulation = V3, PMID = V4
      )
      trrust <- dplyr::filter(trrust, Target %in% !!targets)
      which_not_in_data(trrust, "Target", targets)
      x$all_tf$trrust <- trrust
      x <- methodAdd(x, "基于 TRRUST 数据库 (<https://www.grnpedia.org>) 获取转录因子（TF）与靶基因之间的调控关系，用于构建转录调控网络并筛选关键调控因子。TRRUST 收录了经文献证据支持的人类和小鼠转录调控关系，包含转录因子、靶基因及其激活或抑制作用等信息。通过将基因映射至该数据库，可识别其上游调控转录因子，并进一步构建 TF–target 调控网络，挖掘核心转录因子及潜在调控轴，为解析基因表达调控机制提供依据。")
    }
    data <- .merge_list_by_cols(x$all_tf, by = c("TF", "Target"))
    nodes <- list(
      TF = data$TF, mRNA = data$Target
    )
    nodes <- as_df.lst(lapply(nodes, unique), "type", "name")[, 2:1]
    graph <- igraph::graph_from_data_frame(
      data, directed = TRUE, vertices = nodes
    )
    set.seed(x$seed)
    require(ggraph)
    layout <- ggraph::create_layout(graph, layout = layout)
    p.regTF <- ggraph(layout) +
      geom_edge_link(
        edge_width = .5, color = "grey80", show.legend = FALSE,
        end_cap = ggraph::circle(7, 'mm'),
        arrow = arrow(length = unit(1, 'mm'))) + 
      geom_node_point(aes(color = type, size = type, shape = type)) + 
      ggrepel::geom_label_repel(aes(x = x, y = y, label = name), size = 3) +
      guides(
        size = "none", shape = "none",
        color = guide_legend(override.aes = list(size = 4))
        ) +
      scale_size_manual(values = c(mRNA = 10, TF = 6)) +
      scale_shape_manual(
        values = c(mRNA = 16, TF = 17)
        ) +
      labs(color = "Type") +
      theme_void()
    p.regTF <- set_lab_legend(
      p.regTF,
      glue::glue("{x@sig} TF regulation networking"),
      glue::glue("TF-mRNA 表达调控网络分析|||图中的节点表示对应 mRNA 或 TF，边代表相互作用。")
    )
    x <- snapAdd(x, "以 TF、mRNA 构建 TF-mRNA 表达调控网络{aref(p.regTF)}，如图所示，共 {nrow(nodes)} 个节点，{nrow(data)} 个边。")
    x <- plotsAdd(x, p.regTF)
    return(x)
  })



.filter_chipbase_tf <- function(
  data_chipbase,
  col_tf = "TF",
  col_target = "Target",
  vec_promoter_col = c("up1", "down1"),
  col_n_sample = "total_samples_num",
  num_quantile = 0.75,
  n_min_support = 1L
)
{
  vec_keep_promoter <- apply(
    data_chipbase[, vec_promoter_col, drop = FALSE],
    MARGIN = 1L,
    FUN = function(x)
    {
      any(x > 0L, na.rm = TRUE)
    }
  )

  data_chipbase <- data_chipbase[
    vec_keep_promoter,
    ,
    drop = FALSE
    ]

  message(
    glue::glue(
      "Promoter proximal records retained: {nrow(data_chipbase)}"
    )
  )

  data_pair <- dplyr::rename(
    data_chipbase, n_sample = !!rlang::sym(col_n_sample)
  )

  data_pair <- dplyr::group_by(
    data_pair,
    .data[[col_target]]
  )

  data_pair <- dplyr::mutate(
    data_pair,
    n_cutoff = max(
      n_min_support,
      stats::quantile(
        n_sample,
        probs = num_quantile,
        na.rm = TRUE
      )
    )
  )

  data_pair <- dplyr::filter(
    data_pair,
    n_sample >= n_cutoff
  )

  data_pair <- dplyr::ungroup(
    data_pair
  )

  message(
    glue::glue(
      "Retained TF-target pairs: {nrow(data_pair)}"
    )
  )

  return(data_pair)
}

.merge_list_by_cols <- function(sets, by, keep = seq_along(by)) {
  sets <- lapply(sets, function(x) dplyr::distinct(x[, keep]))
  if (length(sets) > 1) {
    ins <- sets[[1]]
    for (i in 2:length(sets)) {
      ins <- merge(ins, sets[[ i ]], by = by)
    }
  } else {
    ins <- sets[[1]]
  }
  ins
}


get_encori_miRNA_lncRNA <- function(miRNA, dir_db = .prefix("encori", "db"))
{
  if (any(grpl(miRNA, "^miR-"))) {
    stop('any(grpl(miRNA, "^miR-")), the organisms such as hsa- should be added to prefix.')
  }
  # IDs: your query, col: the ID column, res: results table
  dir.create(dir_db, FALSE)
  db <- new_db(file.path(dir_db, "encori_miRNA_lncRNA.rdata"), ".id")
  db <- not(db, miRNA)
  query <- db@query
  if (length(query)) {
    res <- pbapply::pbsapply(query, .api_encori_miRNA_lncRNA, simplify = FALSE)
    res <- frbind(res, idcol = ".id", fill = TRUE)
    db <- upd(db, res)
  }
  res <- dplyr::filter(db@db, .id %in% !!miRNA)
}

.api_encori_miRNA_lncRNA <- function(miRNA) {
  if (length(miRNA) > 1) {
    stop('length(miRNA) > 1.')
  }
  url <- glue::glue("https://rnasysu.com/encori/api/miRNATarget/?assembly=hg38&geneType=lncRNA&miRNA={miRNA}&interNum=1&expNum=1&cellType=all")
  data.table::fread(text = RCurl::getURL(url))
}

get_url_data <- function(expect_filename, url,
  name = formal_name(tools::file_path_sans_ext(basename(expect_filename))),
  dir = .prefix(name, "db"), fun_decompress = utils::unzip)
{
  expect_file <- file.path(dir, expect_filename)
  if (!file.exists(expect_file)) {
    dir.create(dir, FALSE)
    if (is.null(fun_decompress)) {
      download_file <- expect_file
    } else {
      download_file <- file.path(
        dir, paste0(name, ".", tools::file_ext(url))
      )
    }
    utils::download.file(url, destfile = download_file)
    if (!is.null(fun_decompress)) {
      fun_decompress(download_file, exdir = dir)
    }
  }
  message(glue::glue("The `expect_file` exist: {file.exists(expect_file)}"))
  normalizePath(expect_file)
}

which_not_in_data <- function(data, col, items, prefix = NULL, stop = FALSE)
{
  name <- rlang::expr_text(substitute(data))
  whichNot <- !items %in% data[[ col ]]
  if (any(whichNot)) {
    if (!is.null(prefix)) {
      message(glue::glue("{prefix} {name} not got: {bind(items[ whichNot ])}"))
    } else {
      message(glue::glue("{name} not got: {bind(items[ whichNot ])}"))
    }
    if (stop) {
      stop("See above.")
    }
  }
  return(whichNot)
}

map_genes_in_miRNetR <- function(genes, wd, method = "wget", extra = "--no-check-certificate")
{
  if (length(ls(envir = .GlobalEnv))) {
    stop("`miRNetR` is a terrible package with Poor code standards.
      Please run this function in a completely isolated subprocess!")
  }
  dir.create(wd, FALSE)
  setwd(wd)
  message(glue::glue("Set work directory in {normalizePath(wd)}.
      This is done because miRNetR will download the file to the current workspace."))
  fileWeb <- "https://www.xialab.ca/rest/sqlite/mir2gene.sqlite"
  filename <- basename(fileWeb)
  if (!file.exists(filename)) {
    download.file(fileWeb, filename, method = method, extra = extra)
  }
  miRNetR::Init.Data("mir", "gene")
  miRNetR::SetupIndListData(genes, "hsa", "gene", "symbol", "na", "na")
  nms.vec <<- c("gene")
  miRNetR::SetCurrentDataMulti()
  sqlite.path <<- "https://www.xialab.ca/rest/sqlite/"
  .on.public.web <<- FALSE
  miRNetR::QueryMultiList()
  dataSet
}

get_miRNetR.huibang <- function() {
  pak::pkg_install
}

.description_filter_chipbase <- function(num_quantile, n_min_support)
{
  num_quantile_percent <- num_quantile * 100
  glue::glue('为降低不同目标基因间公共 ChIP-seq 数据覆盖度差异带来的偏倚，本研究未采用统一固定的全局样本数阈值进行筛选。首先，仅保留在转录起始位点（TSS）上下游 ±1 kb 启动子邻近区域存在结合证据的 TF-target 相互作用。随后，统计每个 TF-target 关系对应的支持记录数，并在各 target 基因内部进行分位数筛选。具体筛选标准定义为：

$$
n_{sample} \\geq \\max\\left(
<<n_min_support>>,
Q_{<<num_quantile_percent>>}(n_{sample})
\\right)
$$

其中，$n_{sample}$ 表示对应 TF-target 关系的支持记录数，$Q$ 表示各 target 基因内部支持记录数分布的分位数函数。最终保留满足条件的高可信度 TF-target 相互作用用于后续调控网络构建。', .open = "<<", .close = ">>")
}


#
# chipbase <- dplyr::mutate(
#   chipbase, sample_id = strx(name, "ChIPBase3ID[0-9]+"), .before = 1
# )
# data_motif <- ftibble(get_url_data(
#     "hg38_motif.txt.gz",
#     "https://rnasysu.com/chipbase3/data/download/motif/hg38_motif.txt.gz",
#     "chipbase", fun_decompress = NULL
#     ))
# chipbase <- map(
#   chipbase, "sample_id", data_motif, "sample_id", "best_match", col = "Regulator"
# )
# chipbase <- dplyr::mutate(
#   chipbase, TF = Regulator
# )
# chipbase <- dplyr::relocate(
#   chipbase, TF, Target = gene_symbol
# )
