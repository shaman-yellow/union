# ==========================================================================
# workflow of locate
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_locate <- setClass("job_locate", 
  contains = c("job"),
  prototype = prototype(
    pg = "locate",
    info = c("https://mrcuizhe.github.io/interacCircos_documentation/html/users_from_rcircos.html",
      "http://www.rnalocate.org/"),
    cite = "",
    method = "",
    tag = "locate",
    analysis = "染色体定位和亚细胞定位"
    ))

setGeneric("asjob_locate",
  function(x, ...) standardGeneric("asjob_locate"))

setMethod("asjob_locate", signature = c(x = "feature"),
  function(x, ...){
    fea <- resolve_feature_snapAdd_onExit("x", x)
    x <- .job_locate(object = fea)
    dir.create("tmp", FALSE)
    x$gene_data <- expect_local_data(
      "tmp", "biomart", locateFuns$get_gene_positions, list(genes = fea), ...
    )
    return(x)
  })

setMethod("step0", signature = c(x = "job_locate"),
  function(x){
    step_message("Prepare your data with function `job_locate`.")
  })

setMethod("step1", signature = c(x = "job_locate"),
  function(x, skip = FALSE)
  {
    step_message("Chromosome location")
    if (skip) {
      message("Skip this step.")
      return(x)
    }
    p.locateChr <- funPlot(locateFuns$plot_genes_in_RCircos, list(gene_data = x$gene_data))
    p.locateChr <- set_lab_legend(
      p.locateChr,
      glue::glue("{x@sig} Chromosome localization"),
      glue::glue("基因于染色体定位|||外圈数字表示染色体（1-22表示1-22条人类染色体，XY对应性染色体）")
    )
    snap <- glue::glue(
      "基因 {x$gene_data$Gene} 位于 {s(x$gene_data$Chromosome, 'chr', '')} 染色体上"
    )
    x <- snapAdd(x, "如图所示{aref(p.locateChr)}，{bind(snap)}。")
    x <- methodAdd(x, "染色体定位分析可有效揭示基因在染色体上的分布特征。以 R 包 `RCircos` ⟦pkgInfo('RCircos')⟧ 生成基因染色体定位图谱。")
    x <- plotsAdd(x, p.locateChr)
    return(x)
  })

setMethod("step2", signature = c(x = "job_locate"),
  function(x,
    RNALocate = TRUE,
    COMPARTMENTS = FALSE,
    COMPARTMENTS_channels = c(
      "experiments", "knowledge", "predictions", "textmining"
    ),
    COMPARTMENTS_top_n_per_gene = 5L,
    COMPARTMENTS_drop_broad_terms = TRUE,
    COMPARTMENTS_min_confidence = NULL,
    COMPARTMENTS_dir = .prefix("COMPARTMENTS", "db")
  )
  {
    step_message("Get location data.")

    if (RNALocate) {
      org <- "Homo sapiens"
      data <- locateFuns$get_RNALocate_subcellular_data(org = org)
      data <- dplyr::filter(data, RNA_Symbol %in% object(x))

      whichNot <- !object(x) %in% data$RNA_Symbol

      if (any(whichNot)) {
        message(glue::glue("Not got: {object(x)[ whichNot ]}"))
      }

      data <- dplyr::arrange(data, dplyr::desc(RNALocate_Score))
      data <- dplyr::distinct(
        data, RNA_Symbol, Subcellular_Localization, .keep_all = TRUE
      )

      x$locateData_RNALocate <- data
      x$locateData <- data

      p.locateScore <- wrap_scale(
        locateFuns$.plot_subcellular_scores(data),
        length(unique(data$Subcellular_Localization)),
        length(object(x)[!whichNot]),
        pre_height = 3.5,
        min_width = 2
      )

      x <- methodAdd(
        x,
        "从 RNALocate v3.0 (<http://www.rnalocate.org/>) 获取 mRNA 亚细胞定位数据，并用 R 包 `ggplot2` ⟦pkgInfo('ggplot2')⟧将定位和得分数据可视化。"
      )

      p.locateScore <- set_lab_legend(
        p.locateScore,
        glue::glue("{x@sig} RNA Subcellular Localization Distribution"),
        glue::glue("RNA 亚细胞定位分布|||纵坐标为不同基因，横坐标为的预测的蛋白质亚细胞定位得分：RNA 亚细胞定位关联信息来自不同类型的资源，包括实验证据和预测证据；实验证据对置信度评分的贡献应该比预测证据更大；强有力的实验证据应该比薄弱的实验证据提供更可靠的证据；有更多证据支持的 RNA 亚细胞定位关联应比证据较少支持的关联具有更高的置信度评分 (<http://www.rnalocate.org/help>)。")
      )

      top <- dplyr::distinct(data, RNA_Symbol, .keep_all = TRUE)
      snap <- glue::glue(
        "蛋白 {top$RNA_Symbol} 分布于 {top$Subcellular_Localization}"
      )

      x <- snapAdd(
        x,
        "如图{aref(p.locateScore)}，依据 RNALocate 评分准则 (见图注或 <http://www.rnalocate.org/help>) 最有证据证明 {bind(snap)}。"
      )

      x <- plotsAdd(x, p.locateScore)
    }

    if (COMPARTMENTS) {
      res_compartments <- locateFuns$get_COMPARTMENTS_subcellular_data(
        genes = object(x),
        channels = COMPARTMENTS_channels,
        dir = COMPARTMENTS_dir,
        drop_broad_terms = COMPARTMENTS_drop_broad_terms,
        min_confidence = COMPARTMENTS_min_confidence
      )

      data_compartments <- res_compartments$data
      data_summary_all <- locateFuns$summarize_COMPARTMENTS_data(
        data = data_compartments,
        top_n_per_gene = NULL
      )
      data_report <- locateFuns$summarize_COMPARTMENTS_data(
        data = data_compartments,
        top_n_per_gene = COMPARTMENTS_top_n_per_gene
      )
      data_channel_stat <- locateFuns$summarize_COMPARTMENTS_channels(
        data = data_compartments,
        genes = object(x),
        channels = COMPARTMENTS_channels
      )

      x$locateData_COMPARTMENTS <- data_compartments
      x$locateSummary_COMPARTMENTS <- data_summary_all
      x$locateReport_COMPARTMENTS <- data_report
      x$locateStat_COMPARTMENTS <- data_channel_stat

      if (!RNALocate) {
        x$locateData <- data_report
      }

      t.compartments <- locateFuns$format_COMPARTMENTS_report_table(
        data_report
      )

      t.channel <- locateFuns$format_COMPARTMENTS_channel_table(
        data_channel_stat
      )

      t.compartments <- set_lab_legend(
        t.compartments,
        glue::glue("{x@sig} COMPARTMENTS subcellular localization summary"),
        glue::glue(
          "COMPARTMENTS 亚细胞定位证据汇总表|||",
          "该表基于 COMPARTMENTS 数据库整合目标基因的亚细胞定位证据。",
          "每行对应一个基因-亚细胞定位条目，并分别展示 experiments、knowledge、predictions 和 textmining 四类证据通道的置信度分数。",
          "表格按证据通道数量和最高置信度排序，并默认保留每个基因排名靠前的定位条目，以便在报告中简洁展示。"
        )
      )

      t.channel <- set_lab_legend(
        t.channel,
        glue::glue("{x@sig} COMPARTMENTS evidence channel statistics"),
        glue::glue(
          "COMPARTMENTS 证据通道统计表|||",
          "该表展示 COMPARTMENTS 数据库中 experiments、knowledge、predictions 和 textmining 四类证据通道的检索结果数量、覆盖基因数量、定位条目数量和置信度分布。"
        )
      )

      x <- tablesAdd(
        x,
        t.compartments = t.compartments,
        t.channel = t.channel
      )

      p.compartments <- locateFuns$plot_COMPARTMENTS_scores(
        data = data_report,
        channels = COMPARTMENTS_channels
      )

      p.compartments <- wrap_scale(
        p.compartments,
        length(unique(data_report$compartment)),
        length(unique(data_report$query_gene)),
        pre_height = 3.5,
        min_width = 3
      )

      p.compartments <- set_lab_legend(
        p.compartments,
        glue::glue("{x@sig} COMPARTMENTS subcellular localization evidence distribution"),
        glue::glue(
          "COMPARTMENTS 亚细胞定位证据分布图|||",
          "该图展示目标基因在 COMPARTMENTS 数据库四类证据通道中的亚细胞定位证据分布。",
          "横坐标为亚细胞区室，纵坐标为目标基因；点的大小和颜色表示对应证据通道下的置信度分数。",
          "分面分别对应实验、知识库、计算预测和文献挖掘证据，用于展示不同来源证据对目标基因亚细胞定位的支持情况。"
        )
      )

      x <- plotsAdd(x, p.compartments)

      n_gene_input <- length(unique(object(x)))
      n_gene_found <- length(unique(data_compartments$query_gene))
      n_record <- nrow(data_compartments)
      n_location <- nrow(data_summary_all)

      x <- methodAdd(
        x,
        locateFuns$get_COMPARTMENTS_method_text(
          channels = COMPARTMENTS_channels,
          drop_broad_terms = COMPARTMENTS_drop_broad_terms,
          top_n_per_gene = COMPARTMENTS_top_n_per_gene
        )
      )

      x <- snapAdd(
        x,
        glue::glue(
          "基于 COMPARTMENTS 数据库对 {n_gene_input} 个目标基因进行亚细胞定位证据检索，",
          "共获得 {n_record} 条数据库证据记录，覆盖 {n_gene_found} 个目标基因；",
          "经基因-定位层面汇总后得到 {n_location} 个定位证据条目，",
          "四类证据通道的覆盖情况见表格{aref(t.channel)}。"
        )
      )

      x <- snapAdd(
        x,
        glue::glue(
          "为便于报告展示，COMPARTMENTS 结果进一步整理为每个基因排名靠前的定位汇总表，",
          "该表同时展示 experiments、knowledge、predictions 和 textmining 四类证据通道的置信度分数{aref(t.compartments)}；",
          "四类证据通道下的基因-亚细胞区室定位分布见图{aref(p.compartments)}。"
        )
      )
    }

    return(x)
  })



setMethod("step3", signature = c(x = "job_locate"),
  function(x,
    recode = NULL,
    dir_save = create_job_cache_dir(x)
  )
  {
    step_message("Get protein sequence.")

    genes <- object(x)

    if (!is.null(recode)) {
      if (!all(names(recode) %in% genes)) {
        stop('!all(names(recode) %in% genes).')
      }

      genes <- dplyr::recode(genes, !!!recode)
    }

    fun_getseq <- function(...) {
      mart <- new_biomart("hsa")
      get_seq.pro(genes, mart)
    }

    x$seqs <- expect_local_data(
      "tmp", "seq", fun_getseq, list(genes)
    )

    if (!is.null(recode)) {
      recode <- stats::setNames(names(recode), unname(recode))

      if (!all(names(recode) %in% x$seqs$data$hgnc_symbol)) {
        stop('!all(names(recode) %in% x$seqs$data$hgnc_symbol).')
      }

      x$seqs$data <- dplyr::mutate(
        x$seqs$data,
        hgnc_symbol = dplyr::recode(hgnc_symbol, !!!recode)
      )
    }

    x$file_seqs <- union.publish:::write(
      x$seqs$fasta,
      name = "peptide",
      dir = dir_save,
      max = NULL
    )

    x <- methodAdd(
      x,
      "使用 `biomaRt` ⟦pkgInfo('biomaRt')⟧ 获取目标基因对应的蛋白序列。分析首先根据基因 Symbol 检索 ensembl peptide ID 与 canonical transcript 信息，并优先选取经典转录本对应的多肽序列，用于后续 DeepLoc 亚细胞定位预测。"
    )

    return(x)
  })

setMethod("step4", signature = c(x = "job_locate"),
  function(x,
    url = NULL,
    file_json = NULL,
    dir_save = create_job_cache_dir(x),
    use_cache = TRUE,
    overwrite = FALSE,
    add_detail_plots = TRUE,
    stop_if_missing = TRUE
  )
  {
    step_message("Collect DeepLoc result.")

    if (is.null(url) && !is.null(x$DeepLoc_url)) {
      url <- x$DeepLoc_url
    }

    if (is.null(file_json) && !is.null(x$DeepLoc_files$file_json) &&
        file.exists(x$DeepLoc_files$file_json)) {
      file_json <- x$DeepLoc_files$file_json
    }

    res_file <- locateFuns$cache_DeepLoc_summary(
      url = url,
      file_json = file_json,
      dir = dir_save,
      use_cache = use_cache,
      overwrite = overwrite,
      stop_if_missing = stop_if_missing
    )

    res_data <- locateFuns$parse_DeepLoc_json(
      file_json = res_file$file_json
    )

    asset_manifest <- tibble::tibble(
      entry_id = character(0),
      Gene_display = character(0),
      asset_type = character(0),
      remote_file = character(0),
      url = character(0),
      file = character(0),
      file_exists = logical(0)
    )

    plot_map <- tibble::tibble(
      entry_id = character(0),
      Gene_display = character(0),
      file_tree = character(0),
      file_alpha = character(0)
    )

    if (isTRUE(add_detail_plots)) {
      asset_manifest <- locateFuns$cache_DeepLoc_assets_from_json(
        data_assets = res_data$data_assets,
        dir = dir_save,
        use_cache = use_cache,
        overwrite = overwrite
      )

      plot_map <- locateFuns$prepare_DeepLoc_plot_map(
        data_summary = res_data$data_summary,
        asset_manifest = asset_manifest
      )
    }

    x$DeepLoc_url <- url
    x$DeepLoc_files <- res_file
    x$DeepLoc_asset_manifest <- asset_manifest
    x$DeepLoc_plot_map <- plot_map
    x$locateData_DeepLoc <- res_data$data_summary
    x$locateData_DeepLoc_long <- res_data$data_long

    t.deeploc <- locateFuns$format_DeepLoc_report_table(
      res_data$data_summary
    )

    t.deeploc <- set_lab_legend(
      t.deeploc,
      glue::glue("{x@sig} DeepLoc subcellular localization prediction"),
      glue::glue(
        "DeepLoc 蛋白亚细胞定位预测结果表|||",
        "该表展示基于目标蛋白氨基酸序列获得的 DeepLoc 亚细胞定位预测结果。",
        "表中列出各目标蛋白的主要预测定位、定位概率以及可溶性/膜蛋白类型判断，",
        "用于从序列特征角度补充亚细胞定位证据。"
      )
    )

    x <- tablesAdd(x, t.deeploc = t.deeploc)

    p.deeploc <- locateFuns$plot_DeepLoc_scores(
      res_data$data_long
    )

    p.deeploc <- set_lab_legend(
      p.deeploc,
      glue::glue("{x@sig} DeepLoc localization likelihood distribution"),
      glue::glue(
        "DeepLoc 亚细胞定位预测概率分布图|||",
        "该图展示目标蛋白在 DeepLoc 各亚细胞定位类别中的预测概率分布。",
        "横坐标为亚细胞定位类别，纵坐标为目标蛋白；点的大小和颜色表示对应定位类别的预测概率，",
        "用于展示各目标蛋白主要定位及其他候选定位类别之间的概率差异。"
      )
    )

    x <- plotsAdd(x, p.deeploc)

    ps.detail <- list()

    if (isTRUE(add_detail_plots) && nrow(plot_map) > 0L) {
      ps.detail <- locateFuns$plot_DeepLoc_core_images(plot_map)

      if (length(ps.detail) > 0L) {
        lab_detail <- as.character(glue::glue(
          "{x@sig} {plot_map$Gene_display} DeepLoc prediction details"
        ))

        labs_detail <- as.character(glue::glue(
          "DeepLoc {plot_map$Gene_display} 亚细胞定位预测细节图|||",
          "该图展示目标蛋白 {plot_map$Gene_display} 的 DeepLoc 预测细节结果。",
          "左侧为 hierarchical tree，显示模型层级定位分类路径；",
          "右侧为 position importance 分布，显示蛋白序列不同位置对定位预测的贡献。",
          "该结果用于辅助解释 DeepLoc 序列预测，不作为实验定位证据。"
        ))

        ps.detail <- set_lab_legend(ps.detail, lab_detail, labs_detail)
        x <- plotsAdd(x, ps.detail)
      }
    }

    n_gene <- nrow(res_data$data_summary)

    x <- methodAdd(
      x,
      glue::glue(
        "为从蛋白序列层面补充亚细胞定位证据，本研究采用 DeepLoc 对目标基因对应的经典蛋白序列进行预测。",
        "DeepLoc 根据蛋白氨基酸序列识别定位相关特征，输出各亚细胞定位类别的预测概率及主要定位结果。",
        "该结果作为 RNALocate 与 COMPARTMENTS 数据库定位证据之外的序列预测补充。"
      )
    )

    x <- snapAdd(
      x,
      glue::glue(
        "DeepLoc 共完成 {n_gene} 个目标蛋白的亚细胞定位预测，预测结果见表格{aref(t.deeploc)}，",
        "各定位类别的预测概率分布见图{aref(p.deeploc)}。"
      )
    )

    if (length(ps.detail) > 0L) {
      x <- snapAdd(
        x,
        glue::glue(
          "DeepLoc 的 hierarchical tree 与 position importance 结果用于辅助展示模型预测依据。"
        )
      )
    }

    return(x)
  })

locateFuns <- new.env(parent = emptyenv())

locateFuns$get_url_data <- function(expect_filename, url,
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

locateFuns$get_RNALocate_subcellular_data <- function(org = "Homo sapiens", dir = .prefix("mRNA_subcellular", "db"))
{
  filename <- "mRNA subcellular localization information.txt"
  file <- file.path(dir, filename)
  if (!file.exists(file)) {
    dir.create(dir, FALSE)
    url <- "http://www.rnalocate.org/static/download/mRNA%20subcellular%20localization%20information.zip"
    zipfile <- file.path(dir, "mRNA.zip")
    utils::download.file(url, zipfile)
    unzip(zipfile, exdir = dir)
  }
  dplyr::filter(ftibble(file), Species == !!org)
}

locateFuns$get_COMPARTMENTS_urls <- function()
{
  c(
    experiments = "https://download.jensenlab.org/human_compartment_experiments_full.tsv",
    knowledge = "https://download.jensenlab.org/human_compartment_knowledge_full.tsv",
    predictions = "https://download.jensenlab.org/human_compartment_predictions_full.tsv",
    textmining = "https://download.jensenlab.org/human_compartment_textmining_full.tsv"
  )
}

locateFuns$get_COMPARTMENTS_broad_terms <- function()
{
  c(
    "GO:0005575",
    "GO:0110165",
    "GO:0005622",
    "GO:0043226",
    "GO:0043227",
    "GO:0043228",
    "GO:0043229",
    "GO:0043231",
    "GO:0043232"
  )
}

locateFuns$collapse_unique <- function(x, sep = "; ")
{
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]
  paste(unique(x), collapse = sep)
}

locateFuns$read_COMPARTMENTS_table <- function(file)
{
  if (requireNamespace("data.table", quietly = TRUE)) {
    data <- data.table::fread(
      file,
      sep = "\t",
      header = FALSE,
      data.table = FALSE,
      showProgress = FALSE,
      quote = ""
    )
  } else {
    data <- as.data.frame(ftibble(file), stringsAsFactors = FALSE)
  }

  if (ncol(data) < 7L) {
    stop("COMPARTMENTS table should contain at least 7 columns.")
  }

  data <- data[, seq_len(7L), drop = FALSE]
  colnames(data) <- paste0("V", seq_len(7L))
  data
}

locateFuns$standardize_COMPARTMENTS_channel <- function(data,
  channel
)
{
  channel <- match.arg(
    channel,
    c("experiments", "knowledge", "predictions", "textmining")
  )

  if (channel == "textmining") {
    score_raw <- suppressWarnings(as.numeric(data$V5))
    confidence_score <- suppressWarnings(as.numeric(data$V6))

    data_res <- data.frame(
      channel = channel,
      protein_id = as.character(data$V1),
      gene_symbol = as.character(data$V2),
      go_id = as.character(data$V3),
      compartment = as.character(data$V4),
      source = rep("Text mining", nrow(data)),
      evidence = as.character(data$V7),
      score_raw = score_raw,
      confidence_score = confidence_score,
      url = as.character(data$V7),
      stringsAsFactors = FALSE
    )
  } else {
    confidence_score <- suppressWarnings(as.numeric(data$V7))

    data_res <- data.frame(
      channel = channel,
      protein_id = as.character(data$V1),
      gene_symbol = as.character(data$V2),
      go_id = as.character(data$V3),
      compartment = as.character(data$V4),
      source = as.character(data$V5),
      evidence = as.character(data$V6),
      score_raw = confidence_score,
      confidence_score = confidence_score,
      url = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  data_res$gene_symbol[is.na(data_res$gene_symbol)] <- ""
  data_res$compartment[is.na(data_res$compartment)] <- ""
  data_res$go_id[is.na(data_res$go_id)] <- ""
  data_res
}

locateFuns$read_COMPARTMENTS_channel <- function(channel,
  genes,
  dir = .prefix("COMPARTMENTS", "db"),
  drop_broad_terms = TRUE,
  min_confidence = NULL
)
{
  urls <- locateFuns$get_COMPARTMENTS_urls()

  if (!channel %in% names(urls)) {
    stop(glue::glue("Unknown COMPARTMENTS channel: {channel}."))
  }

  file <- locateFuns$get_url_data(
    expect_filename = basename(urls[[channel]]),
    url = urls[[channel]],
    name = tools::file_path_sans_ext(basename(urls[[channel]])),
    dir = dir,
    fun_decompress = NULL
  )

  data <- locateFuns$read_COMPARTMENTS_table(file)
  data <- locateFuns$standardize_COMPARTMENTS_channel(
    data = data,
    channel = channel
  )

  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & genes != ""]

  data_gene <- data.frame(
    query_gene = genes,
    gene_key = toupper(genes),
    stringsAsFactors = FALSE
  )

  data$gene_key <- toupper(data$gene_symbol)
  data <- data[data$gene_key %in% data_gene$gene_key, , drop = FALSE]

  if (nrow(data) == 0L) {
    data$query_gene <- character(0)
    data$gene_key <- NULL
    return(tibble::as_tibble(data))
  }

  data$query_gene <- data_gene$query_gene[
    match(data$gene_key, data_gene$gene_key)
  ]
  data$gene_key <- NULL

  if (isTRUE(drop_broad_terms)) {
    data <- data[
      !data$go_id %in% locateFuns$get_COMPARTMENTS_broad_terms(),
      ,
      drop = FALSE
    ]
  }

  if (!is.null(min_confidence)) {
    data <- data[
      !is.na(data$confidence_score) &
        data$confidence_score >= min_confidence,
      ,
      drop = FALSE
    ]
  }

  tibble::as_tibble(data)
}

locateFuns$get_COMPARTMENTS_subcellular_data <- function(genes,
  channels = c("experiments", "knowledge", "predictions", "textmining"),
  dir = .prefix("COMPARTMENTS", "db"),
  drop_broad_terms = TRUE,
  min_confidence = NULL
)
{
  channels <- unique(as.character(channels))
  channels <- channels[!is.na(channels) & channels != ""]

  lst_data <- lapply(
    channels,
    function(channel) {
      locateFuns$read_COMPARTMENTS_channel(
        channel = channel,
        genes = genes,
        dir = dir,
        drop_broad_terms = drop_broad_terms,
        min_confidence = min_confidence
      )
    }
  )

  names(lst_data) <- channels
  data <- do.call(rbind, lapply(lst_data, as.data.frame))

  if (is.null(data) || nrow(data) == 0L) {
    data <- data.frame(
      channel = character(0),
      protein_id = character(0),
      gene_symbol = character(0),
      go_id = character(0),
      compartment = character(0),
      source = character(0),
      evidence = character(0),
      score_raw = numeric(0),
      confidence_score = numeric(0),
      url = character(0),
      query_gene = character(0),
      stringsAsFactors = FALSE
    )
  }

  rownames(data) <- NULL

  list(
    data = tibble::as_tibble(data),
    data_by_channel = lst_data
  )
}

locateFuns$summarize_COMPARTMENTS_data <- function(data,
  top_n_per_gene = 5L
)
{
  if (is.null(data) || nrow(data) == 0L) {
    return(tibble::tibble(
      query_gene = character(0),
      compartment = character(0),
      best_confidence = numeric(0),
      n_channel = integer(0),
      n_record = integer(0),
      evidence_channels = character(0),
      source_summary = character(0),
      experiments = numeric(0),
      knowledge = numeric(0),
      predictions = numeric(0),
      textmining = numeric(0)
    ))
  }

  data <- as.data.frame(data, stringsAsFactors = FALSE)

  data_channel <- dplyr::group_by(
    data,
    query_gene,
    compartment,
    channel
  )

  data_channel <- dplyr::summarise(
    data_channel,
    confidence = ifelse(
      all(is.na(confidence_score)),
      NA_real_,
      max(confidence_score, na.rm = TRUE)
    ),
    n_record = dplyr::n(),
    source = locateFuns$collapse_unique(source),
    .groups = "drop"
  )

  data_key <- unique(data_channel[, c("query_gene", "compartment")])
  channels <- c("experiments", "knowledge", "predictions", "textmining")

  for (channel in channels) {
    data_tmp <- data_channel[
      data_channel$channel == channel,
      c("query_gene", "compartment", "confidence"),
      drop = FALSE
    ]
    colnames(data_tmp)[3L] <- channel

    data_key <- merge(
      data_key,
      data_tmp,
      by = c("query_gene", "compartment"),
      all.x = TRUE,
      sort = FALSE
    )
  }

  data_record <- dplyr::group_by(
    data_channel,
    query_gene,
    compartment
  )

  data_record <- dplyr::summarise(
    data_record,
    n_record = sum(n_record, na.rm = TRUE),
    evidence_channels = locateFuns$collapse_unique(channel),
    source_summary = locateFuns$collapse_unique(source),
    .groups = "drop"
  )

  data_res <- merge(
    data_key,
    data_record,
    by = c("query_gene", "compartment"),
    all.x = TRUE,
    sort = FALSE
  )

  data_res$n_channel <- rowSums(!is.na(data_res[, channels, drop = FALSE]))

  data_res$best_confidence <- apply(
    data_res[, channels, drop = FALSE],
    1L,
    function(x) {
      x <- as.numeric(x)
      if (all(is.na(x))) {
        return(NA_real_)
      }
      max(x, na.rm = TRUE)
    }
  )

  data_res <- data_res[
    order(
      data_res$query_gene,
      -data_res$n_channel,
      -data_res$best_confidence,
      data_res$compartment
    ),
    ,
    drop = FALSE
  ]

  data_res <- data_res[
    c(
      "query_gene",
      "compartment",
      "best_confidence",
      "n_channel",
      "n_record",
      "evidence_channels",
      "source_summary",
      channels
    )
  ]

  rownames(data_res) <- NULL

  if (!is.null(top_n_per_gene)) {
    data_res <- do.call(
      rbind,
      lapply(
        split(data_res, data_res$query_gene),
        function(data_item) {
          data_item[seq_len(min(top_n_per_gene, nrow(data_item))), , drop = FALSE]
        }
      )
    )
    rownames(data_res) <- NULL
  }

  tibble::as_tibble(data_res)
}

locateFuns$summarize_COMPARTMENTS_channels <- function(data,
  genes,
  channels = c("experiments", "knowledge", "predictions", "textmining")
)
{
  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & genes != ""]
  channels <- unique(as.character(channels))

  data <- as.data.frame(data, stringsAsFactors = FALSE)

  data_res <- lapply(
    channels,
    function(channel) {
      data_channel <- data[data$channel == channel, , drop = FALSE]

      confidence <- data_channel$confidence_score
      confidence <- confidence[!is.na(confidence)]

      data.frame(
        channel = channel,
        n_record = nrow(data_channel),
        n_gene = length(unique(data_channel$query_gene)),
        n_missing_gene = length(setdiff(genes, unique(data_channel$query_gene))),
        n_compartment = length(unique(data_channel$compartment)),
        median_confidence = ifelse(
          length(confidence) > 0L,
          stats::median(confidence),
          NA_real_
        ),
        max_confidence = ifelse(
          length(confidence) > 0L,
          max(confidence),
          NA_real_
        ),
        stringsAsFactors = FALSE
      )
    }
  )

  data_res <- do.call(rbind, data_res)
  rownames(data_res) <- NULL
  tibble::as_tibble(data_res)
}

locateFuns$format_COMPARTMENTS_report_table <- function(data)
{
  data <- as.data.frame(data, stringsAsFactors = FALSE)

  if (nrow(data) == 0L) {
    return(tibble::as_tibble(data))
  }

  data <- data.frame(
    Gene = data$query_gene,
    Compartment = data$compartment,
    `Best confidence` = round(data$best_confidence, 2L),
    `Evidence channels` = data$evidence_channels,
    Experiments = round(data$experiments, 2L),
    Knowledge = round(data$knowledge, 2L),
    Predictions = round(data$predictions, 2L),
    `Text mining` = round(data$textmining, 2L),
    `Record count` = data$n_record,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  tibble::as_tibble(data)
}

locateFuns$format_COMPARTMENTS_channel_table <- function(data)
{
  data <- as.data.frame(data, stringsAsFactors = FALSE)

  if (nrow(data) == 0L) {
    return(tibble::as_tibble(data))
  }

  data$channel <- factor(
    data$channel,
    levels = c("experiments", "knowledge", "predictions", "textmining")
  )

  data <- data[order(data$channel), , drop = FALSE]

  data <- data.frame(
    Channel = as.character(data$channel),
    `Evidence records` = data$n_record,
    `Covered genes` = data$n_gene,
    `Missing genes` = data$n_missing_gene,
    `Compartments` = data$n_compartment,
    `Median confidence` = round(data$median_confidence, 2L),
    `Max confidence` = round(data$max_confidence, 2L),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  tibble::as_tibble(data)
}


locateFuns$prepare_COMPARTMENTS_plot_data <- function(data,
  channels = c("experiments", "knowledge", "predictions", "textmining")
)
{
  data <- as.data.frame(data, stringsAsFactors = FALSE)

  if (is.null(data) || nrow(data) == 0L) {
    return(tibble::tibble(
      query_gene = character(0),
      compartment = character(0),
      channel = character(0),
      channel_label = character(0),
      confidence = numeric(0)
    ))
  }

  channels <- unique(as.character(channels))
  channels <- channels[!is.na(channels) & channels != ""]
  channels <- intersect(
    channels,
    c("experiments", "knowledge", "predictions", "textmining")
  )

  channel_lab <- c(
    experiments = "Experiments",
    knowledge = "Knowledge",
    predictions = "Predictions",
    textmining = "Text mining"
  )

  lst_data <- lapply(
    channels,
    function(channel) {
      data.frame(
        query_gene = data$query_gene,
        compartment = data$compartment,
        channel = channel,
        channel_label = channel_lab[channel],
        confidence = suppressWarnings(as.numeric(data[[channel]])),
        stringsAsFactors = FALSE
      )
    }
  )

  data_plot <- do.call(rbind, lst_data)
  data_plot <- data_plot[
    !is.na(data_plot$confidence) &
      data_plot$confidence > 0,
    ,
    drop = FALSE
  ]

  rownames(data_plot) <- NULL
  tibble::as_tibble(data_plot)
}

locateFuns$wrap_axis_text <- function(x,
  width = 24L
)
{
  x <- as.character(x)
  vapply(
    x,
    function(item) {
      paste(strwrap(item, width = width), collapse = "\n")
    },
    character(1)
  )
}

locateFuns$plot_COMPARTMENTS_scores <- function(data,
  channels = c("experiments", "knowledge", "predictions", "textmining"),
  compartment_wrap_width = 22L,
  point_alpha = .75,
  min_point_size = 2,
  max_point_size = 8
)
{
  data_plot <- locateFuns$prepare_COMPARTMENTS_plot_data(
    data = data,
    channels = channels
  )

  if (nrow(data_plot) == 0L) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::annotate(
          "text",
          x = 0,
          y = 0,
          label = "No COMPARTMENTS localization evidence was available."
        )
    )
  }

  data_order <- stats::aggregate(
    confidence ~ compartment,
    data = data_plot,
    FUN = max
  )

  data_order <- data_order[
    order(-data_order$confidence, data_order$compartment),
    ,
    drop = FALSE
  ]

  data_plot$compartment <- factor(
    data_plot$compartment,
    levels = data_order$compartment
  )

  data_plot$query_gene <- factor(
    data_plot$query_gene,
    levels = rev(unique(data$query_gene))
  )

  data_plot$channel_label <- factor(
    data_plot$channel_label,
    levels = c("Experiments", "Knowledge", "Predictions", "Text mining")
  )

  ggplot2::ggplot(
    data_plot,
    ggplot2::aes(
      x = compartment,
      y = query_gene,
      size = confidence,
      color = confidence
    )
  ) +
    ggplot2::geom_point(alpha = point_alpha) +
    ggplot2::facet_wrap(
      ~channel_label,
      ncol = 2L
    ) +
    ggplot2::scale_x_discrete(
      labels = function(x) {
        locateFuns$wrap_axis_text(x, width = compartment_wrap_width)
      }
    ) +
    ggplot2::scale_size_continuous(
      range = c(min_point_size, max_point_size)
    ) +
    ggplot2::scale_color_gradient(
      low = "#2b8cbe",
      high = "#d7301f"
    ) +
    ggplot2::labs(
      x = "Subcellular compartment",
      y = "Gene",
      size = "Confidence",
      color = "Confidence"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      panel.grid.major = ggplot2::element_line(linewidth = .25),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(
        fill = "grey95",
        color = NA
      ),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

locateFuns$get_COMPARTMENTS_method_text <- function(channels,
  drop_broad_terms = TRUE,
  top_n_per_gene = 5L
)
{
  text_channel <- paste(channels, collapse = ", ")

  text_filter <- if (isTRUE(drop_broad_terms)) {
    "为减少过于泛化的 Gene Ontology cellular component 条目对结果展示的影响，分析时去除了 cellular_component、cellular anatomical entity、intracellular anatomical structure、organelle 等宽泛定位条目。"
  } else {
    "分析保留 COMPARTMENTS 原始定位条目，不额外去除宽泛 Gene Ontology cellular component 条目。"
  }

  glue::glue(
    "基因亚细胞定位证据从 COMPARTMENTS 数据库获取。该数据库整合实验、知识库、计算预测和文献挖掘等多来源证据，",
    "并将蛋白与亚细胞区室的关联映射至 Gene Ontology cellular component 条目。",
    "本分析使用的证据通道包括 {text_channel}。",
    "首先按目标基因符号筛选数据库记录，然后在基因-亚细胞区室-证据通道层面保留最高置信度分数，",
    "再将 experiments、knowledge、predictions 和 textmining 四类证据整理为同一张汇总表。",
    "{text_filter}",
    "为便于报告展示，每个基因默认展示排名靠前的 {top_n_per_gene} 个定位条目，排序依据为支持该定位的证据通道数量及最高置信度分数。"
  )
}

locateFuns$.plot_subcellular_scores <- function(data) {
  p <- ggplot(data, aes(x = Subcellular_Localization, 
      y = reorder(RNA_Symbol, RNALocate_Score), 
      size = RNALocate_Score,
      color = RNALocate_Score)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(range = c(3, 10)) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Subcellular Localization", 
    y = "RNA",
    size = "Score",
    color = "Score")
  p
}

locateFuns$plot_genes_in_RCircos <- function(gene_data) {
  require(RCircos)
  data(UCSC.HG38.Human.CytoBandIdeogram, package = "RCircos")
  cyto.info <- UCSC.HG38.Human.CytoBandIdeogram
  RCircos::RCircos.Set.Core.Components(
    cyto.info, 
    chr.exclude = NULL,
    tracks.inside = 2,
    tracks.outside = 0
  )
  params <- RCircos::RCircos.Get.Plot.Parameters()
  params$text.size <- 1
  params$point.size <- 1.2
  RCircos::RCircos.Reset.Plot.Parameters(params)
  RCircos::RCircos.Set.Plot.Area()
  RCircos::RCircos.Chromosome.Ideogram.Plot()
  RCircos::RCircos.Gene.Connector.Plot(gene_data, track.num = 1, side = "in")
  RCircos::RCircos.Gene.Name.Plot(gene_data, name.col = 4, track.num = 2, side = "in")
}


locateFuns$get_gene_positions <- function(genes) {
  ensembl <- new_biomart()
  gene_positions <- biomaRt::getBM(
    attributes = c("chromosome_name", "start_position", "end_position", "hgnc_symbol"),
    filters = "hgnc_symbol",
    values = genes,
    mart = ensembl
  )
  valid_chrs <- c(as.character(1:22), "X", "Y")
  gene_positions <- gene_positions[ gene_positions$chromosome_name %in% valid_chrs, ]
  gene_positions$chromosome <- paste0("chr", gene_positions$chromosome_name)
  result <- data.frame(
    Chromosome = gene_positions$chromosome,
    Start = gene_positions$start_position,
    End = gene_positions$end_position,
    Gene = gene_positions$hgnc_symbol,
    stringsAsFactors = FALSE
  )
  chr_order <- order(gene_positions$chromosome_name)
  result <- result[chr_order, ]
  if (nrow(result) < length(genes)) {
    not_found <- setdiff(genes, result$Gene)
    message(glue::glue("Not found: {bind(not_found)}"))
  }
  return(result)
}

locateFuns$get_DeepLoc_job_id <- function(url)
{
  if (is.null(url) || is.na(url) || url == "") {
    return(NA_character_)
  }

  job_id <- strx(url, "(?<=jobid=)[A-Za-z0-9]+")

  if (is.na(job_id)) {
    job_id <- strx(url, "[A-Fa-f0-9]{20,}")
  }

  job_id
}

locateFuns$get_DeepLoc_tmp_url <- function(url)
{
  job_id <- locateFuns$get_DeepLoc_job_id(url)

  if (is.na(job_id)) {
    return(NA_character_)
  }

  paste0(
    "https://services.healthtech.dtu.dk/services/DeepLoc-1.0/tmp/",
    job_id,
    "/"
  )
}

locateFuns$as_absolute_DeepLoc_url <- function(path)
{
  path <- as.character(path)
  path[is.na(path)] <- ""

  ifelse(
    grepl("^https?://", path),
    path,
    paste0("https://services.healthtech.dtu.dk", path)
  )
}

locateFuns$download_cached_file <- function(url,
  file,
  use_cache = TRUE,
  overwrite = FALSE,
  mode = "wb"
)
{
  if (isTRUE(use_cache) && !isTRUE(overwrite) &&
      file.exists(file) && file.info(file)$size > 0) {
    return(TRUE)
  }

  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)

  ok <- tryCatch(
    {
      utils::download.file(
        url,
        destfile = file,
        mode = mode,
        quiet = TRUE
      )
      file.exists(file) && file.info(file)$size > 0
    },
    error = function(e) {
      FALSE
    },
    warning = function(w) {
      FALSE
    }
  )

  if (!isTRUE(ok) && file.exists(file)) {
    unlink(file)
  }

  isTRUE(ok)
}

locateFuns$is_DeepLoc_json <- function(file)
{
  if (is.na(file) || !file.exists(file) || file.info(file)$size == 0) {
    return(FALSE)
  }

  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop('Package "jsonlite" is required.')
  }

  ok <- tryCatch(
    {
      jsonlite::fromJSON(file, simplifyVector = FALSE)
      TRUE
    },
    error = function(e) {
      FALSE
    }
  )

  isTRUE(ok)
}

locateFuns$cache_DeepLoc_summary <- function(url = NULL,
  file_json = NULL,
  dir = ".",
  use_cache = TRUE,
  overwrite = FALSE,
  stop_if_missing = TRUE
)
{
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  if (!is.null(file_json)) {
    file_json <- normalizePath(file_json, mustWork = FALSE)
  }

  if (!is.null(file_json) && locateFuns$is_DeepLoc_json(file_json) &&
      !isTRUE(overwrite)) {
    return(list(
      file_json = file_json,
      job_id = locateFuns$get_DeepLoc_job_id(url)
    ))
  }

  if (is.null(url) || is.na(url) || url == "") {
    if (isTRUE(stop_if_missing)) {
      stop("DeepLoc JSON file was not found. Please provide `url` or `file_json`.")
    }

    return(list(file_json = NA_character_, job_id = NA_character_))
  }

  job_id <- locateFuns$get_DeepLoc_job_id(url)
  tmp_url <- locateFuns$get_DeepLoc_tmp_url(url)

  if (is.na(job_id) || is.na(tmp_url)) {
    stop("Could not resolve DeepLoc job ID from URL.")
  }

  if (is.null(file_json)) {
    file_json <- file.path(dir, paste0("deeploc_", job_id, "_summary.json"))
  }

  urls_json <- c(
    paste0(tmp_url, "results_", job_id, ".json"),
    paste0(tmp_url, "results.json"),
    paste0(tmp_url, "summary.json"),
    paste0(tmp_url, "output.json")
  )

  for (url_json in urls_json) {
    locateFuns$download_cached_file(
      url = url_json,
      file = file_json,
      use_cache = use_cache,
      overwrite = overwrite,
      mode = "wb"
    )

    if (locateFuns$is_DeepLoc_json(file_json)) {
      break
    }
  }

  if (!locateFuns$is_DeepLoc_json(file_json)) {
    if (isTRUE(stop_if_missing)) {
      stop("DeepLoc JSON summary was not found. Please download JSON Summary manually and pass it as `file_json`.")
    }

    file_json <- NA_character_
  }

  list(
    file_json = file_json,
    job_id = job_id
  )
}

locateFuns$get_DeepLoc_gene_label <- function(entry_id)
{
  entry_id <- as.character(entry_id)
  entry_id <- sub("_ENSP[0-9]+.*$", "", entry_id)
  entry_id <- sub("\\|.*$", "", entry_id)
  entry_id
}

locateFuns$pick_DeepLoc_png <- function(x)
{
  x <- as.character(unlist(x))
  x <- x[!is.na(x) & x != ""]
  x <- x[grepl("[.]png($|[?])", x, ignore.case = TRUE)]

  if (!length(x)) {
    return(NA_character_)
  }

  x[1L]
}

locateFuns$parse_DeepLoc_json <- function(file_json)
{
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop('Package "jsonlite" is required.')
  }

  obj <- jsonlite::fromJSON(file_json, simplifyVector = FALSE)
  sequences <- obj$sequences

  if (is.null(sequences) || !length(sequences)) {
    stop("Could not find DeepLoc sequence predictions in JSON output.")
  }

  entry_ids <- names(sequences)

  if (is.null(entry_ids) || any(entry_ids == "")) {
    entry_ids <- vapply(
      sequences,
      function(sequence) {
        name <- sequence$Name
        if (is.null(name) || is.na(name) || name == "") NA_character_ else name
      },
      character(1)
    )
  }

  lst_long <- lapply(
    seq_along(sequences),
    function(i) {
      sequence <- sequences[[i]]
      entry_id <- entry_ids[i]
      localization <- as.character(unlist(sequence$Localization))
      likelihood <- suppressWarnings(as.numeric(unlist(sequence$Likelihood)))

      n <- min(length(localization), length(likelihood))

      if (n == 0L) {
        return(NULL)
      }

      data.frame(
        entry_id = entry_id,
        Gene_display = locateFuns$get_DeepLoc_gene_label(entry_id),
        localization = localization[seq_len(n)],
        likelihood = likelihood[seq_len(n)],
        rank = seq_len(n),
        stringsAsFactors = FALSE
      )
    }
  )

  lst_long <- Filter(Negate(is.null), lst_long)
  data_long <- do.call(rbind, lst_long)
  rownames(data_long) <- NULL

  lst_summary <- lapply(
    seq_along(sequences),
    function(i) {
      sequence <- sequences[[i]]
      entry_id <- entry_ids[i]
      localization <- as.character(unlist(sequence$Localization))
      likelihood <- suppressWarnings(as.numeric(unlist(sequence$Likelihood)))
      membrane <- as.character(sequence$Membrane)
      membrane_likelihood <- suppressWarnings(as.numeric(unlist(sequence$Membrane_likelihood)))

      data.frame(
        entry_id = entry_id,
        Gene = locateFuns$get_DeepLoc_gene_label(entry_id),
        Prediction = if (length(localization)) localization[1L] else NA_character_,
        Likelihood = if (length(likelihood)) likelihood[1L] else NA_real_,
        Protein_type = if (length(membrane)) membrane[1L] else NA_character_,
        Soluble_likelihood = if (length(membrane_likelihood) >= 1L) membrane_likelihood[1L] else NA_real_,
        Membrane_likelihood = if (length(membrane_likelihood) >= 2L) membrane_likelihood[2L] else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  )

  data_summary <- do.call(rbind, lst_summary)
  rownames(data_summary) <- NULL

  lst_assets <- lapply(
    seq_along(sequences),
    function(i) {
      sequence <- sequences[[i]]
      entry_id <- entry_ids[i]

      data.frame(
        entry_id = entry_id,
        Gene_display = locateFuns$get_DeepLoc_gene_label(entry_id),
        asset_type = c("tree", "alpha"),
        remote_file = c(
          locateFuns$pick_DeepLoc_png(sequence$Tree),
          locateFuns$pick_DeepLoc_png(sequence$Attention)
        ),
        stringsAsFactors = FALSE
      )
    }
  )

  data_assets <- do.call(rbind, lst_assets)
  data_assets <- data_assets[
    !is.na(data_assets$remote_file) & data_assets$remote_file != "",
    ,
    drop = FALSE
  ]
  rownames(data_assets) <- NULL

  list(
    data_summary = tibble::as_tibble(data_summary),
    data_long = tibble::as_tibble(data_long),
    data_assets = tibble::as_tibble(data_assets)
  )
}

locateFuns$cache_DeepLoc_assets_from_json <- function(data_assets,
  dir = ".",
  use_cache = TRUE,
  overwrite = FALSE
)
{
  data_assets <- as.data.frame(data_assets, stringsAsFactors = FALSE)

  if (nrow(data_assets) == 0L) {
    return(tibble::as_tibble(data_assets))
  }

  dir_assets <- file.path(dir, "deeploc_assets")
  dir.create(dir_assets, recursive = TRUE, showWarnings = FALSE)

  data_assets$url <- locateFuns$as_absolute_DeepLoc_url(data_assets$remote_file)

  data_assets$file <- file.path(
    dir_assets,
    paste0(
      formal_name(data_assets$Gene_display),
      "_",
      data_assets$asset_type,
      ".png"
    )
  )

  data_assets$file_exists <- vapply(
    seq_len(nrow(data_assets)),
    function(i) {
      locateFuns$download_cached_file(
        url = data_assets$url[i],
        file = data_assets$file[i],
        use_cache = use_cache,
        overwrite = overwrite,
        mode = "wb"
      )
    },
    logical(1)
  )

  data_assets$file <- normalizePath(data_assets$file, mustWork = FALSE)
  tibble::as_tibble(data_assets)
}

locateFuns$prepare_DeepLoc_plot_map <- function(data_summary,
  asset_manifest
)
{
  data_summary <- as.data.frame(data_summary, stringsAsFactors = FALSE)
  asset_manifest <- as.data.frame(asset_manifest, stringsAsFactors = FALSE)

  if (nrow(data_summary) == 0L || nrow(asset_manifest) == 0L) {
    return(tibble::tibble(
      entry_id = character(0),
      Gene_display = character(0),
      file_tree = character(0),
      file_alpha = character(0)
    ))
  }

  data_tree <- asset_manifest[
    asset_manifest$asset_type == "tree" & asset_manifest$file_exists,
    c("entry_id", "file"),
    drop = FALSE
  ]

  data_alpha <- asset_manifest[
    asset_manifest$asset_type == "alpha" & asset_manifest$file_exists,
    c("entry_id", "file"),
    drop = FALSE
  ]

  colnames(data_tree)[2L] <- "file_tree"
  colnames(data_alpha)[2L] <- "file_alpha"

  data_map <- merge(
    data_summary[, c("entry_id", "Gene"), drop = FALSE],
    data_tree,
    by = "entry_id",
    all.x = TRUE,
    sort = FALSE
  )

  data_map <- merge(
    data_map,
    data_alpha,
    by = "entry_id",
    all.x = TRUE,
    sort = FALSE
  )

  data_map$Gene_display <- data_map$Gene
  data_map <- data_map[, c("entry_id", "Gene_display", "file_tree", "file_alpha"), drop = FALSE]
  data_map <- data_map[
    file.exists(data_map$file_tree) | file.exists(data_map$file_alpha),
    ,
    drop = FALSE
  ]
  rownames(data_map) <- NULL

  tibble::as_tibble(data_map)
}

locateFuns$format_DeepLoc_report_table <- function(data)
{
  data <- as.data.frame(data, stringsAsFactors = FALSE)

  if (nrow(data) == 0L) {
    return(tibble::tibble())
  }

  data <- data[, c(
    "Gene", "Prediction", "Likelihood", "Protein_type",
    "Soluble_likelihood", "Membrane_likelihood"
  ), drop = FALSE]

  data$Likelihood <- round(data$Likelihood, 4L)
  data$Soluble_likelihood <- round(data$Soluble_likelihood, 4L)
  data$Membrane_likelihood <- round(data$Membrane_likelihood, 4L)
  tibble::as_tibble(data)
}

locateFuns$plot_DeepLoc_scores <- function(data_long,
  label_wrap_width = 18L
)
{
  data_long <- as.data.frame(data_long, stringsAsFactors = FALSE)

  if (nrow(data_long) == 0L) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate(
          "text",
          x = 0,
          y = 0,
          label = "No DeepLoc prediction result was available."
        ) +
        ggplot2::theme_void()
    )
  }

  data_long$localization_plot <- vapply(
    data_long$localization,
    function(x) paste(strwrap(x, width = label_wrap_width), collapse = "\n"),
    character(1)
  )

  ggplot2::ggplot(
    data_long,
    ggplot2::aes(
      x = localization_plot,
      y = Gene_display,
      size = likelihood,
      color = likelihood
    )
  ) +
    ggplot2::geom_point(alpha = .9) +
    ggplot2::scale_size_continuous(range = c(1.5, 7)) +
    ggplot2::scale_color_gradient(low = "#bdd7e7", high = "#08519c") +
    ggplot2::labs(
      x = "Subcellular localization",
      y = "Gene",
      size = "Likelihood",
      color = "Likelihood"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.major = ggplot2::element_line(color = "grey90")
    )
}

locateFuns$read_png_grob <- function(file)
{
  if (is.na(file) || !file.exists(file)) {
    return(NULL)
  }

  mat_img <- png::readPNG(file)
  grid::rasterGrob(mat_img, interpolate = TRUE)
}

locateFuns$plot_DeepLoc_image_pair <- function(file_tree,
  file_alpha
)
{
  grob_tree <- locateFuns$read_png_grob(file_tree)
  grob_alpha <- locateFuns$read_png_grob(file_alpha)

  p <- ggplot2::ggplot() +
    ggplot2::coord_cartesian(xlim = c(0, 2), ylim = c(0, 1), expand = FALSE) +
    ggplot2::theme_void()

  if (!is.null(grob_tree)) {
    p <- p + ggplot2::annotation_custom(
      grob_tree,
      xmin = 0,
      xmax = 1,
      ymin = 0,
      ymax = 1
    )
  }

  if (!is.null(grob_alpha)) {
    p <- p + ggplot2::annotation_custom(
      grob_alpha,
      xmin = 1,
      xmax = 2,
      ymin = 0,
      ymax = 1
    )
  }

  wrap(p, 10, 5)
}

locateFuns$plot_DeepLoc_core_images <- function(plot_map)
{
  plot_map <- as.data.frame(plot_map, stringsAsFactors = FALSE)

  if (nrow(plot_map) == 0L) {
    return(list())
  }

  ps <- lapply(
    seq_len(nrow(plot_map)),
    function(i) {
      locateFuns$plot_DeepLoc_image_pair(
        file_tree = plot_map$file_tree[i],
        file_alpha = plot_map$file_alpha[i]
      )
    }
  )

  names(ps) <- formal_name(plot_map$Gene_display)
  ps
}

