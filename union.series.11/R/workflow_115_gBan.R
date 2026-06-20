# ==========================================================================
# workflow of gBan
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_gBan <- setClass("job_gBan", 
  contains = c("job"),
  prototype = prototype(
    pg = "gBan",
    info = c("https://github.com/HamidHadipour/GraphBAN"),
    cite = "",
    method = "",
    tag = "gBan",
    analysis = "GraphBAN 药物预测"
    ))

job_gBan <- function()
{
  .job_gBan()
}

setGeneric("asjob_gBan",
  function(x, ...) standardGeneric("asjob_gBan"))

setMethod("asjob_gBan", signature = c(x = "feature"),
  function(x){
    fea <- resolve_feature_snapAdd_onExit("x", x)
    x <- .job_gBan(object = fea)
    x <- methodAdd(x, "为进一步挖掘筛选得到的关键基因在临床转化中的潜在应用价值，基于 GraphBAN 模型开展药物预测分析。该分析旨在从分子靶点层面连接疾病相关基因与可干预药物之间的桥梁，识别可能影响疾病发生发展的药物分子，为后续机制研究及药物重定位提供理论依据，并为精准治疗策略的制定提供潜在靶点与候选干预方案。")
    return(x)
  })

setMethod("step0", signature = c(x = "job_gBan"),
  function(x){
    step_message("Prepare your data with function `job_gBan`.")
  })

setMethod("step1", signature = c(x = "job_gBan"),
  function(x, dir_save = paste0("GraphBAN_", x@sig),
    db = c("batman", "zinc", "cmnpd", "dgidb", "drugbank", "custom"),
    batman = FALSE, zinc = FALSE, cmnpd = FALSE, dgidb = FALSE,
    drugbank = FALSE, custom = FALSE, file_dgidb = NULL,
    file_batman_compounds_info = getOption("file_batman_compounds_info"),
    db_custom = NULL, db_custom_name = "Custom compound set", recode = NULL)
  {
    step_message("Got amino acid sequence and drug smiles.")
    x$dir_save <- dir_save
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
      recode <- setNames(names(recode), unname(recode))
      if (!all(names(recode) %in% x$seqs$data$hgnc_symbol)) {
        stop('!all(names(recode) %in% x$seqs$data$hgnc_symbol).')
      }
      x$seqs$data <- dplyr::mutate(
        x$seqs$data, hgnc_symbol = dplyr::recode(hgnc_symbol, !!!recode)
      )
    }
    x$file_seqs <- union.publish:::write(
      x$seqs$fasta, name = "peptide", dir = dir_save, max = NULL
    )
    x <- methodAdd(x, "以 `biomaRt` ⟦pkgInfo('biomaRt')⟧ 获取多肽序列 (先以 Symbol 获取 'ensembl_peptide_id' 以及 'transcript_is_canonical'，从而选择经典转录本的 'ensembl_peptide_id' 获取多肽序列) 。")
    x$smiles_compounds <- list()
    db <- match.arg(db)
    assign(db, TRUE)
    if (batman) {
      if (is.null(file_batman_compounds_info)) {
        stop('is.null(file_batman_compounds_info).')
      }
      if (!file.exists(file_batman_compounds_info)) {
        stop('!file.exists(file_batman_compounds_info).')
      }
      data_batman <- ftibble(file_batman_compounds_info)
      if (!any(colnames(data_batman) == "SMILES")) {
        stop('!any(colnames(data_batman) == "SMILES").')
      }
      x$smiles_compounds$data_batman <- dplyr::distinct(data_batman, smiles = SMILES)
      x <- methodAdd(x, "获取 BATMAN-TCM (<http://bionet.ncpsb.org.cn/batman-tcm>) 化合物 SMILES 结构式。")
    }
    if (zinc) {
      data_zinc <- ftibble(get_url_data(
        "zinc.csv",
        "https://raw.githubusercontent.com/aspuru-guzik-group/chemical_vae/master/models/zinc_properties/250k_rndm_zinc_drugs_clean_3.csv",
        "zinc", fun_decompress = NULL
      ))
      data_zinc <- dplyr::mutate(
        data_zinc, smiles = sub("\n$", "", smiles)
      )
      x$smiles_compounds$data_zinc <- data_zinc
      x <- methodAdd(x, "获取 ZINC-250k 小分子库 (<https://www.kaggle.com/datasets/basu369victor/zinc250k>) 的化合物 SMILES 结构式。")
    }
    if (cmnpd) {
      data_cmnpd <- ftibble(get_url_data(
          "cmnpd.tsv",
          "https://www.cmnpd.org/cmnpd/supplement/Downloads/CMNPD_1.0_calc_prop.tsv",
          "cmnpd", fun_decompress = NULL
          ))
      x$smiles_compounds$data_cmnpd <- dplyr::distinct(data_cmnpd, smiles = SMILES)
      x <- methodAdd(x, "获取 CMNPD 数据库 (<https://www.cmnpd.org>) 的化合物 SMILES 结构式。")
    }
    if (drugbank) {
      data_drugbank <- ftibble(pg("db_drugbank"))
      x$smiles_compounds$data_drugbank <- dplyr::distinct(data_drugbank, smiles = SMILES)
      x <- methodAdd(x, "获取 Drugbank 数据库 (<https://go.drugbank.com/>) 的化合物 SMILES 结构式。")
    }
    if (dgidb) {
      if (!is.null(file_dgidb)) {
        data_dgidb <- ftibble(file_dgidb)
        data_dgidb <- dplyr::select(data_dgidb, gene, drug)
      } else {
        data_dgidb <- ftibble(get_url_data(
            "interactions.tsv",
            "https://dgidb.org/data/2024-Dec/interactions.tsv",
            "dgidb", fun_decompress = NULL
            ))
        data_dgidb <- dplyr::select(
          data_dgidb, gene = gene_claim_name, drug = drug_claim_name
        )
      }
      expect_package("PubChemR", "3.0.0")
      drugs <- s(unique(data_dgidb$drug), "^CHEMBL:", "")
      ndrugs <- length(unique(drugs))
      cli::cli_alert_info("PubChemR::get_cids")
      cids <- expect_local_data(
        "tmp", "pubchemr_cids", PubChemR::get_cids, list(
          identifier = drugs,
          namespace = "name"
        )
      )
      cids <- PubChemR::CIDs(cids)
      data_dgidb <- map(
        data_dgidb, "drug", cids, "Name", "CID", col = "CID"
      )
      data_dgidb <- dplyr::filter(data_dgidb, !is.na(data_dgidb$CID))
      smiles <- get_smiles_batch(unique(data_dgidb$CID))
      data_dgidb <- map(
        data_dgidb, "CID", smiles, "CID", "SMILES", col = "SMILES"
      )
      x$data_dgidb <- data_dgidb
      x$smiles_compounds$data_dgidb <- dplyr::distinct(
        data_dgidb, smiles = SMILES
      )
      x <- methodAdd(
        x, "以 DGIdb 数据库 (<https://www.cmnpd.org>) 初步预测与输入基因存在相互作用的候选药物化合物。共计得到 {nrow(ndrugs)} 条记录。剔除无法从 PubChemR 搜索到对应化合物信息记录的条目。余下 {nrow(data_dgidb)} 条记录，按各基因统计为: {try_snap(data_dgidb, 'gene', 'drug')}。"
      )
    }
    if (custom) {
      data_custom <- gBanFuns$resolve_custom_compound_db(
        db_custom = db_custom,
        db_custom_name = db_custom_name,
        require_smiles = TRUE
      )
      x$custom_compound_db <- data_custom
      x$smiles_compounds$data_custom <- dplyr::distinct(
        data_custom,
        smiles
      )
      x <- methodAdd(
        x,
        glue::glue(
          "采用自定义候选化合物库“{db_custom_name}”作为 GraphBAN 输入化合物来源。",
          "该化合物库仅保留具有 SMILES 结构式的化合物，用于后续 GraphBAN 药物–靶点相互作用预测；",
          "原始化合物名称、PubChem CID、InChIKey 及来源证据等信息保留为候选化合物注释字段。"
        )
      )
    }
    return(x)
  })

setMethod("step2", signature = c(x = "job_gBan"),
  function(x, cl = 10, mem = 1000, w.cutoff = 1000, 
    rerun = FALSE, filter_by = c("rcdk", "rdkit"))
  {
    step_message("Filter drugs.")
    compounds <- unique(unlist(lapply(x$smiles_compounds, function(x) x$smiles)))
    filter_by <- match.arg(filter_by)
    if (filter_by == "rdkit") {
      x$info_compounds <- expect_local_data(
        "tmp", "compounds_weight_rdkit", inBatches_get_compounds_weight.rdkit,
        list(smiles_list = compounds), rerun = rerun
      )
      x <- methodAdd(
        x, "以 RDKit 过滤无法解析的分子结构，剔除分子量 &gt; {w.cutoff}，含重金属 (以原子序号大于钙为条件过滤) 的化合物，生成标准化 SMILES。"
      )
    } else if (filter_by == "rcdk") {
      x$info_compounds <- expect_local_data(
        "tmp", "compounds_weight", inBatches_get_compounds_weight.rcdk,
        list(smiles_list = compounds, mem = mem, cl = cl), ignore = "cl", rerun = rerun
      )
      x <- methodAdd(
        x, "以 R 包 `rcdk` ⟦pkgInfo('rcdk')⟧ 过滤无法解析的分子结构，剔除分子量 &gt; {w.cutoff}，含重金属 (以原子序号大于钙为条件过滤) 的化合物。"
      )
    }
    input_compounds <- dplyr::filter(
      x$info_compounds, MolecularWeight < w.cutoff, !HasHeavyMetal
    )
    input_compounds <- dplyr::distinct(input_compounds, smiles)
    x$input_compounds <- dplyr::mutate(input_compounds, id = seq_len(nrow(input_compounds)), .before = 1)
    combn <- expand.grid(x$input_compounds$id, x$seqs$data$hgnc_symbol)
    combn <- dplyr::rename(combn, id = 1, hgnc_symbol = 2)
    combn <- Reduce(merge, list(combn, x$input_compounds, x$seqs$data))
    combn <- dplyr::mutate(combn, Y = 1)
    combn <- dplyr::relocate(
      combn, hgnc_symbol, id, SMILES = smiles, Protein = peptide, Y
    )
    if (!is.null(x$data_dgidb)) {
      layout <- dplyr::select(x$data_dgidb, gene, SMILES)
      combn <- merge(
        combn, layout, by.x = c("hgnc_symbol", "SMILES"), 
        by.y = c("gene", "SMILES")
      )
    }
    x$file_combn <- file.path(
      x$dir_save, "graphBan_input.csv"
    )
    write.csv(combn, x$file_combn, row.names = FALSE)
    x$combn <- combn

    if (!is.null(x$custom_compound_db)) {
      data_custom_keep <- merge(
        x$input_compounds,
        as.data.frame(x$custom_compound_db, stringsAsFactors = FALSE),
        by = "smiles",
        all.x = TRUE
      )
      x$custom_input_compound_db <- tibble::as_tibble(data_custom_keep)
    }

    x <- snapAdd(x, "从数据库获取到的化合物共 {length(compounds)} 个。经过滤后得到 {nrow(input_compounds)} 个唯一化合物。")
    return(x)
  })

setMethod("step3", signature = c(x = "job_gBan"),
  function(x){
    step_message("Do nothing")
    return(x)
  })

setMethod("step4", signature = c(x = "job_gBan"),
  function(x,
    pattern = "graphBan_res_",
    cutoff = .95,
    min_cutoff = .5,
    max_cutoff = 1,
    cutoff_step = .05,
    target_min = 100L,
    target_max = 1000L,
    target_center = NULL,
    reRead = FALSE,
    method_model = c("auto", "intersection", "union"),
    method_keep = c("auto", "all", "respective"),
    plot_venn = TRUE,
    plot_model_venn_max_gene = 6L,
    plot_network = TRUE,
    network_layout = "fr",
    network_label_top_n = 40L
  )
  {
    step_message("Collate GraphBAN prediction results.")

    method_model <- match.arg(method_model)
    method_keep <- match.arg(method_keep)

    if (is.null(x$combn) || nrow(x$combn) == 0L) {
      stop("GraphBAN input combinations were not found in x$combn.")
    }

    files_res <- list.files(x$dir_save, pattern, full.names = TRUE)

    if (!length(files_res)) {
      stop("No GraphBAN result file was found.")
    }

    cutoff_grid <- gBanFuns$resolve_cutoff_grid(
      cutoff = cutoff,
      min_cutoff = min_cutoff,
      max_cutoff = max_cutoff,
      cutoff_step = cutoff_step
    )

    data_pred <- gBanFuns$read_graphBan_predictions(
      files_res = files_res,
      combn = x$combn,
      pattern = pattern,
      reRead = reRead
    )

    data_plan <- gBanFuns$evaluate_graphBan_plans(
      data_pred = data_pred,
      cutoff_grid = cutoff_grid,
      target_min = target_min,
      target_max = target_max,
      target_center = target_center,
      method_model = method_model,
      method_keep = method_keep
    )

    res_decision <- gBanFuns$select_graphBan_plan(
      data_pred = data_pred,
      data_plan = data_plan
    )

    data_decision <- res_decision$data_decision
    data_pair_model <- res_decision$data_pair_model
    data_selected <- res_decision$data_selected

    if (nrow(data_selected) == 0L) {
      warning("No candidate compound was retained from GraphBAN predictions.")
      return(x)
    }

    selected_cutoff <- data_decision$cutoff[1L]
    selected_model_strategy <- data_decision$model_strategy[1L]
    selected_gene_strategy <- data_decision$gene_strategy[1L]
    n_model_total <- length(unique(data_pred$model))

    data_model_stat <- gBanFuns$summarize_graphBan_model_counts(
      data_pred = data_pred,
      cutoff = selected_cutoff
    )

    data_gene_stat <- gBanFuns$summarize_graphBan_gene_selection(
      data_pair_model = data_pair_model,
      data_selected = data_selected,
      n_model_total = n_model_total
    )

    data_flow <- gBanFuns$summarize_graphBan_selection_flow(
      data_pred = data_pred,
      data_pair_model = data_pair_model,
      data_selected = data_selected,
      cutoff = selected_cutoff,
      model_strategy = selected_model_strategy,
      gene_strategy = selected_gene_strategy
    )

    x$res_graphBan_all <- data_pred
    x$res_graphBan_pair <- data_pair_model
    x$res_graphBan <- data_selected
    x$graphBan_plan <- data_plan
    x$graphBan_decision <- data_decision
    x$graphBan_selection_flow <- data_flow
    x$graphBan_model_stat <- data_model_stat
    x$graphBan_gene_stat <- data_gene_stat

    if (is.null(x$lst_refine)) {
      x$lst_refine <- list()
    }

    x$lst_refine$graphBan_selection <- list(
      decision = data_decision,
      plan = data_plan,
      flow = data_flow,
      selected = data_selected,
      model_integrated = data_pair_model,
      model_stat = data_model_stat,
      gene_stat = data_gene_stat
    )

    x$smiles_keep <- unique(data_selected$SMILES)
    x$smiles_from_gban <- x$smiles_keep
    x$admet_skipped <- TRUE
    x$swiss_skipped <- TRUE
    x$smiles_from_admet <- NULL
    x$smiles_from_admet_toxicity <- NULL
    x$smiles_for_swiss <- NULL
    x$smiles_from_swiss <- NULL
    x$t_candidate_admet <- NULL
    x$t_candidate_swissAdme <- NULL
    x$split_by_genes <- split(data_selected, data_selected$hgnc_symbol)

    x$file_smiles_for_admet <- file.path(x$dir_save, "smiles_for_admet.txt")
    writeLines(x$smiles_keep, x$file_smiles_for_admet)

    t.decision <- set_lab_legend(
      data_decision,
      glue::glue("{x@sig} GraphBAN candidate selection strategy"),
      glue::glue(
        "GraphBAN 候选化合物筛选策略表|||",
        "该表汇总 GraphBAN 预测结果采用的互作概率阈值、模型整合方式、靶点整合方式，",
        "以及最终获得的候选化合物数量。"
      )
    )

    t.flow <- set_lab_legend(
      data_flow,
      glue::glue("{x@sig} GraphBAN candidate filtering summary"),
      glue::glue(
        "GraphBAN 候选化合物筛选流程统计表|||",
        "该表按筛选阶段汇总候选化合物-靶点互作对和唯一候选化合物数量，",
        "用于说明从原始预测、概率阈值筛选、模型整合到多靶点整合后的候选集合变化。"
      )
    )

    t.model_stat <- set_lab_legend(
      data_model_stat,
      glue::glue("{x@sig} GraphBAN model-level prediction statistics"),
      glue::glue(
        "GraphBAN 模型层面预测统计表|||",
        "该表按靶点基因和 GraphBAN 模型汇总达到互作概率阈值的候选化合物数量，",
        "用于展示不同模型对候选化合物的预测支持情况。"
      )
    )

    t.gene_stat <- set_lab_legend(
      data_gene_stat,
      glue::glue("{x@sig} GraphBAN target-level candidate statistics"),
      glue::glue(
        "GraphBAN 靶点层面候选化合物统计表|||",
        "该表按靶点基因汇总模型整合后及最终候选集合中的化合物数量，",
        "并列出候选互作概率分数的分布情况。"
      )
    )

    x <- tablesAdd(
      x,
      t.decision = t.decision,
      t.flow = t.flow,
      t.model_stat = t.model_stat,
      t.gene_stat = t.gene_stat
    )

    if (nrow(data_model_stat) > 0L) {
      p.model_stat <- gBanFuns$plot_graphBan_model_counts(data_model_stat)
      p.model_stat <- set_lab_legend(
        p.model_stat,
        glue::glue("{x@sig} GraphBAN candidate number by target and model"),
        glue::glue(
          "GraphBAN 候选化合物数量统计图|||",
          "该图按靶点基因和 GraphBAN 模型展示达到互作概率阈值的候选化合物数量，",
          "用于比较不同靶点及不同模型下的候选预测规模。"
        )
      )
      x <- plotsAdd(x, p.model_stat)
    }

    if (isTRUE(plot_venn) && selected_model_strategy == "intersection") {
      genes <- names(x$split_by_genes)

      if (length(genes) <= plot_model_venn_max_gene) {
        data_pass <- dplyr::filter(data_pred, pred >= selected_cutoff)
        ps.model_venn <- list()

        for (gene in genes) {
          data_gene <- dplyr::filter(data_pass, hgnc_symbol == gene)
          lst_model <- lapply(
            split(data_gene, data_gene$model),
            function(data_item) {
              unique(data_item$SMILES)
            }
          )

          if (length(lst_model) > 1L) {
            ps.model_venn[[gene]] <- new_venn(
              lst = lst_model,
              force_upset = FALSE
            )
          }
        }

        if (length(ps.model_venn) > 0L) {
          genes_venn <- names(ps.model_venn)
          lab_model_venn <- as.character(glue::glue(
            "{x@sig} {genes_venn} GraphBAN model intersection"
          ))
          labs_model_venn <- as.character(glue::glue(
            "GraphBAN {genes_venn} 模型交集图|||",
            "该图展示靶点 {genes_venn} 在不同 GraphBAN 模型中达到互作概率阈值的候选化合物重叠情况。",
            "交集区域表示同时获得多个模型支持的候选化合物。"
          ))
          ps.model_venn <- set_lab_legend(
            ps.model_venn,
            lab_model_venn,
            labs_model_venn
          )
          x <- plotsAdd(x, ps.model_venn = ps.model_venn)
        }
      }
    }

    text_network <- ""

    if (isTRUE(plot_network) && selected_gene_strategy != "intersection" &&
        length(unique(data_selected$hgnc_symbol)) > 1L) {
      res_network <- gBanFuns$prepare_graphBan_target_compound_network(
        data_selected = data_selected
      )

      if (nrow(res_network$edges) > 0L) {
        x$graphBan_network_nodes <- res_network$nodes
        x$graphBan_network_edges <- res_network$edges
        x$graphBan_network_graph <- res_network$graph
        x$graphBan_compound_map <- res_network$compound_map

        p.graphBan_network <- gBanFuns$plot_graphBan_target_compound_network(
          res_network = res_network,
          layout = network_layout,
          label_top_n = network_label_top_n,
          seed = x$seed
        )

        p.graphBan_network <- set_lab_legend(
          p.graphBan_network,
          glue::glue("{x@sig} GraphBAN target-compound network"),
          glue::glue(
            "GraphBAN 靶点-候选化合物网络图|||",
            "该图展示 GraphBAN 筛选后候选化合物与靶点蛋白之间的预测关联。",
            "图中靶点节点表示输入蛋白，化合物节点表示候选化合物；边表示对应化合物-靶点组合通过互作概率筛选。",
            "若某一候选化合物同时连接多个靶点，则提示该化合物可能具有多靶点预测关联。"
          )
        )

        x <- plotsAdd(x, p.graphBan_network)
        text_network <- glue::glue(
          " 候选化合物与靶点之间的预测关联见网络图{aref(p.graphBan_network)}。"
        )
      }
    }

    if (isTRUE(plot_venn) && selected_gene_strategy == "intersection" &&
        length(x$split_by_genes) > 1L) {
      lst_gene <- lapply(
        split(data_pair_model, data_pair_model$hgnc_symbol),
        function(data_item) {
          unique(data_item$SMILES)
        }
      )

      p.gene_venn <- new_venn(
        lst = lst_gene,
        force_upset = length(lst_gene) > 5L
      )
      p.gene_venn <- set_lab_legend(
        p.gene_venn,
        glue::glue("{x@sig} intersection of predicted compounds across targets"),
        glue::glue(
          "GraphBAN 多靶点候选化合物交集图|||",
          "该图展示不同靶点基因对应候选化合物的重叠情况，",
          "交集部分表示在多个靶点预测结果中共同出现的候选化合物。"
        )
      )
      x <- plotsAdd(x, p.gene_venn)
    }

    text_cutoff <- gBanFuns$format_cutoff(selected_cutoff)
    text_model_strategy <- gBanFuns$get_model_strategy_text(
      selected_model_strategy,
      n_model_total = n_model_total
    )
    text_gene_strategy <- gBanFuns$get_gene_strategy_text(
      selected_gene_strategy
    )
    text_flow <- gBanFuns$format_graphBan_flow_text(data_flow)
    text_gene <- gBanFuns$format_graphBan_gene_text(data_gene_stat)

    x <- methodAdd(
      x,
      glue::glue(
        "基于 GraphBAN（Graph-Based Attention Network）模型预测候选化合物与靶点蛋白之间的相互作用概率分数（Interaction Probability Score，范围 0-1）。",
        "按照 GraphBAN 预测结果的判定标准，互作概率分数大于 0.5 的化合物-靶点组合视为潜在互作对；",
        "本分析最终以 {text_cutoff} 作为互作概率筛选阈值。",
        "在模型层面，{text_model_strategy}；在多靶点层面，{text_gene_strategy}。"
      )
    )

    x <- snapAdd(
      x,
      glue::glue(
        "GraphBAN 预测结果经互作概率阈值筛选、模型层面整合和靶点层面整合后，",
        "最终保留 {length(x$smiles_keep)} 个唯一候选化合物。{text_flow}{text_gene}{text_network}"
      )
    )

    return(x)
  })



setMethod("step5", signature = c(x = "job_gBan"),
  function(x,
    file_admet = NULL,
    admet_tool = c("admetlab", "admetsar"),
    route = c("systemic", "topical_nasal"),
    cutoff = .7,
    warning_cutoff = .5,
    toxicity_cols = NULL,
    optional_toxicity_cols = NULL,
    include_optional_toxicity = FALSE,
    max_toxicity_flags = 0L,
    max_warning_flags = Inf,
    prepare_swiss_input = TRUE,
    filter_long_smiles_for_swiss = TRUE,
    swiss_max_smiles_chars = 200L,
    add_candidate_table = FALSE,
    add_stat_table = TRUE,
    skip = FALSE
  )
  {
    step_message("ADMET safety evaluation.")
    admet_tool <- match.arg(admet_tool)
    route <- match.arg(route)

    smiles <- x$smiles_from_gban

    if (is.null(smiles) || !length(smiles)) {
      stop("No GraphBAN candidate SMILES was found in x$smiles_from_gban.")
    }

    if (isTRUE(skip)) {
      x$admet_skipped <- TRUE
      x$swiss_skipped <- TRUE
      x$admet <- NULL
      x$admet_eval <- NULL
      x$admet_filter_toxicity <- NULL
      x$admet_filter <- NULL
      x$admet_toxicity_cols <- NULL
      x$admet_filter_stat <- NULL
      x$admet_selection_flow <- NULL
      x$swissAdme <- NULL
      x$swissAdme_eval <- NULL
      x$swiss_selection_flow <- NULL
      x$swiss_selection_decision <- NULL
      x$t_candidate_admet <- NULL
      x$t_candidate_swissAdme <- NULL
      x$smiles_from_swiss <- NULL
      x$smiles_keep <- x$smiles_from_admet <- smiles
      x$smiles_from_admet_toxicity <- smiles

      if (isTRUE(prepare_swiss_input)) {
        x$smiles_for_swiss <- smiles
        x$file_smiles_for_swiss <- file.path(x$dir_save, "smiles_for_swiss.txt")
        writeLines(x$smiles_for_swiss, x$file_smiles_for_swiss)
      } else {
        x$smiles_for_swiss <- NULL
        x$file_smiles_for_swiss <- NULL
      }

      message("ADMET safety evaluation was skipped; GraphBAN-retained SMILES were passed to the next layer.")
      return(x)
    }

    x$admet_skipped <- FALSE

    if (is.null(file_admet)) {
      stop("is.null(file_admet).")
    }

    data_admet <- ftibble(file_admet)
    data_admet <- as.data.frame(data_admet, stringsAsFactors = FALSE)

    data_admet <- gBanFuns$resolve_admet_smiles_identity(
      data_admet = data_admet,
      smiles = smiles,
      admet_tool = admet_tool
    )
    x$admet_smiles_diagnostics <- attr(data_admet, "smiles_diagnostics")

    if (nrow(data_admet) == 0L) {
      stop("No ADMET record matched x$smiles_from_gban.")
    }

    admet_config <- gBanFuns$get_admet_filter_config(
      data = data_admet,
      admet_tool = admet_tool,
      route = route,
      toxicity_cols = toxicity_cols,
      optional_toxicity_cols = optional_toxicity_cols,
      include_optional_toxicity = include_optional_toxicity
    )
    toxicity_cols <- admet_config$hard_cols
    warning_cols <- admet_config$warning_cols

    if (!length(toxicity_cols)) {
      stop("No valid core toxicity endpoint column was found in ADMET table.")
    }

    data_eval <- gBanFuns$evaluate_admet_toxicity(
      data = data_admet,
      toxicity_cols = toxicity_cols,
      warning_cols = warning_cols,
      cutoff = cutoff,
      warning_cutoff = warning_cutoff,
      max_toxicity_flags = max_toxicity_flags,
      max_warning_flags = max_warning_flags
    )

    data_pass_toxicity <- data_eval[
      data_eval$pass_admet_toxicity,
      ,
      drop = FALSE
    ]

    data_pass <- data_pass_toxicity
    data_pass$smiles_nchar <- nchar(as.character(data_pass$raw_smiles))

    n_long_smiles <- 0L

    if (isTRUE(prepare_swiss_input) && isTRUE(filter_long_smiles_for_swiss)) {
      data_pass <- data_pass[
        !is.na(data_pass$smiles_nchar) &
          data_pass$smiles_nchar <= swiss_max_smiles_chars,
        ,
        drop = FALSE
      ]

      n_long_smiles <- nrow(data_pass_toxicity) - nrow(data_pass)
    }

    x$admet <- tibble::as_tibble(data_admet)
    x$admet_eval <- tibble::as_tibble(data_eval)
    x$admet_filter_toxicity <- tibble::as_tibble(data_pass_toxicity)
    x$admet_filter <- tibble::as_tibble(data_pass)
    x$admet_tool <- admet_tool
    x$admet_route <- route
    x$admet_toxicity_cols <- toxicity_cols
    x$admet_warning_cols <- warning_cols
    x$admet_filter_config <- admet_config
    stat_hard <- gBanFuns$summarize_admet_toxicity_filter(
      data_eval = data_eval,
      toxicity_cols = toxicity_cols,
      cutoff = cutoff
    )
    stat_hard$endpoint_role <- "core filter"
    stat_warning <- if (length(warning_cols)) {
      gBanFuns$summarize_admet_toxicity_filter(
        data_eval = data_eval,
        toxicity_cols = warning_cols,
        cutoff = warning_cutoff
      )
    } else {
      stat_hard[0L, , drop = FALSE]
    }
    stat_warning$endpoint_role <- "risk annotation"
    x$admet_filter_stat <- dplyr::bind_rows(stat_hard, stat_warning)

    x$admet_selection_flow <- data.frame(
      stage = c(
        "GraphBAN candidates",
        admet_config$stage_label_en,
        "Structure-scope control"
      ),
      n_compound = c(
        length(unique(smiles)),
        length(unique(data_pass_toxicity$raw_smiles)),
        length(unique(data_pass$raw_smiles))
      ),
      stringsAsFactors = FALSE
    )

    x$swiss_input_stat <- data.frame(
      prepare_swiss_input = isTRUE(prepare_swiss_input),
      filter_long_smiles_for_swiss = isTRUE(filter_long_smiles_for_swiss),
      swiss_max_smiles_chars = swiss_max_smiles_chars,
      n_after_toxicity = length(unique(data_pass_toxicity$raw_smiles)),
      n_long_smiles_removed = n_long_smiles,
      n_for_next_step = length(unique(data_pass$raw_smiles)),
      stringsAsFactors = FALSE
    )

    cols_report <- gBanFuns$resolve_admet_report_cols(
      data = data_pass,
      toxicity_cols = toxicity_cols,
      warning_cols = warning_cols,
      admet_tool = admet_tool,
      route = route
    )

    t.candidate_admet <- data_pass[, cols_report, drop = FALSE]
    t.candidate_admet <- dplyr::rename(
      t.candidate_admet,
      SMILES = raw_smiles,
      Safety_decision = admet_decision,
      Hard_risk_n = n_toxicity_risk,
      Max_hard_risk = max_toxicity_risk,
      Warning_risk_n = n_warning_risk,
      Max_warning_risk = max_warning_risk,
      Major_risk_endpoints = major_risk_endpoints
    )

    t.candidate_admet <- dplyr::mutate(
      t.candidate_admet,
      dplyr::across(
        dplyr::where(is.numeric),
        function(x) round(x, 3L)
      )
    )

    x$t_candidate_admet <- tibble::as_tibble(t.candidate_admet)

    if (isTRUE(add_candidate_table)) {
      t.candidate_admet <- set_lab_legend(
        t.candidate_admet,
        glue::glue("{x@sig} candidate compounds retained by {admet_config$tool_label} safety evaluation"),
        glue::glue(
          "候选化合物安全性初筛保留表|||",
          "该表展示经 {admet_config$tool_label} 安全性初筛后保留的候选化合物。",
          "筛选时重点考察 {admet_config$route_label} 相关核心风险终点，",
          "并保留核心风险指标未达到显著风险阈值的化合物。"
        )
      )

      x <- tablesAdd(
        x,
        t.candidate_admet = t.candidate_admet
      )
    }

    if (isTRUE(add_stat_table)) {
      t.admet_stat <- x$admet_filter_stat
      t.admet_stat <- set_lab_legend(
        t.admet_stat,
        glue::glue("{x@sig} {admet_config$tool_label} safety endpoint summary"),
        glue::glue(
          "候选化合物安全性指标统计表|||",
          "该表统计 {admet_config$tool_label} 各安全性终点在候选化合物中的风险分布情况，",
          "包括达到相应风险阈值的化合物数量及比例。"
        )
      )

      x <- tablesAdd(
        x,
        t.admet_stat = t.admet_stat
      )
    }

    x$smiles_keep <- x$smiles_from_admet <- data_pass$raw_smiles
    x$smiles_from_admet_toxicity <- data_pass_toxicity$raw_smiles
    x$swiss_skipped <- TRUE
    x$swissAdme <- NULL
    x$swissAdme_eval <- NULL
    x$swiss_selection_flow <- NULL
    x$swiss_selection_decision <- NULL
    x$t_candidate_swissAdme <- NULL
    x$smiles_from_swiss <- NULL

    if (isTRUE(prepare_swiss_input)) {
      x$smiles_for_swiss <- data_pass$raw_smiles
      x$file_smiles_for_swiss <- file.path(x$dir_save, "smiles_for_swiss.txt")
      writeLines(x$smiles_for_swiss, x$file_smiles_for_swiss)
    } else {
      x$smiles_for_swiss <- NULL
      x$file_smiles_for_swiss <- NULL
    }

    n_input <- length(unique(smiles))
    n_pass_toxicity <- length(unique(data_pass_toxicity$raw_smiles))
    n_pass <- length(unique(data_pass$raw_smiles))
    text_endpoint <- gBanFuns$format_admet_endpoint_text(
      stat = stat_hard,
      cutoff = cutoff
    )
    text_cutoff <- gBanFuns$format_cutoff(cutoff)
    text_risk_rule <- if (max_toxicity_flags == 0L) {
      glue::glue("全部核心毒性指标均低于 {text_cutoff}")
    } else {
      glue::glue("达到显著风险阈值的核心毒性指标数量不超过 {max_toxicity_flags} 个")
    }

    text_structure_method <- ""
    text_structure_snap <- ""

    if (isTRUE(prepare_swiss_input) && isTRUE(filter_long_smiles_for_swiss) &&
        n_long_smiles > 0L) {
      text_structure_method <- glue::glue(
        "同时对结构式异常冗长的条目进行控制，剔除 SMILES 字符长度超过 {swiss_max_smiles_chars} 的化合物，",
        "以减少超大分子或复杂结构对小分子候选物后续评价的影响。"
      )

      text_structure_snap <- glue::glue(
        "经结构式长度控制后，保留 {n_pass} 个候选化合物用于后续评价。"
      )
    }

    x <- methodAdd(
      x,
      glue::glue(
        "将 GraphBAN 筛选得到的候选化合物 SMILES 输入 {admet_config$tool_label} 平台（{admet_config$tool_url}）进行安全性风险评估。",
        "{admet_config$method_scope}",
        "对于分类风险模型，输出值可作为化合物属于相应风险类别的概率或风险分数进行解释。",
        "本分析依据较高置信风险优先排除的原则，以核心风险终点分数大于等于 {text_cutoff} 作为显著风险判定阈值，",
        "并保留{text_risk_rule}的候选化合物；其他局部接触或系统暴露相关风险终点作为风险注释保留。{text_structure_method}"
      )
    )

    x <- snapAdd(
      x,
      glue::glue(
        "{admet_config$tool_label} 安全性评估前共有 {n_input} 个候选化合物；",
        "经核心风险终点筛选后保留 {n_pass_toxicity} 个候选化合物。",
        "{text_structure_snap}",
        "{text_endpoint}"
      )
    )

    return(x)
  })

setMethod("step6", signature = c(x = "job_gBan"),
  function(x,
    file_swiss = NULL,
    method = c("adaptive", "strict", "core", "lipinski"),
    min_keep = 5L,
    lipinski_min_pass = 3L,
    bioavailability_min = .1,
    tpsa_max = 140,
    mw_max = 500,
    mw_max_relaxed = 800,
    logp_max = 5,
    logp_col = c("Consensus_Log_P", "MLOGP", "XLOGP3", "WLOGP"),
    add_candidate_table = TRUE,
    skip = FALSE
  )
  {
    step_message("SwissADME drug-likeness evaluation.")

    method <- match.arg(method)

    if (isTRUE(skip)) {
      smiles <- if (!is.null(x$smiles_from_admet)) {
        x$smiles_from_admet
      } else {
        x$smiles_from_gban
      }

      if (is.null(smiles) || !length(smiles)) {
        stop("No candidate SMILES was found from GraphBAN or ADMETlab results.")
      }

      x$swiss_skipped <- TRUE
      x$swissAdme <- NULL
      x$swissAdme_eval <- NULL
      x$swiss_selection_flow <- NULL
      x$swiss_selection_decision <- NULL
      x$t_candidate_swissAdme <- NULL
      x$smiles_keep <- x$smiles_from_swiss <- smiles
      message("SwissADME drug-likeness evaluation was skipped; inherited candidate SMILES were passed to PubChem annotation.")
      return(x)
    }

    x$swiss_skipped <- TRUE
    x$swissAdme <- NULL
    x$swissAdme_eval <- NULL
    x$swiss_selection_flow <- NULL
    x$swiss_selection_decision <- NULL
    x$t_candidate_swissAdme <- NULL
    x$smiles_from_swiss <- NULL
    context_input <- gBanFuns$resolve_candidate_filter_context(x)

    if (is.null(file_swiss)) {
      stop("is.null(file_swiss).")
    }

    data_swiss <- ftibble(file_swiss)
    data_swiss <- as.data.frame(data_swiss, stringsAsFactors = FALSE)

    smiles_for_swiss <- if (isTRUE(context_input$has_admet) &&
        !is.null(x$smiles_for_swiss)) {
      x$smiles_for_swiss
    } else if (isTRUE(context_input$has_admet)) {
      x$smiles_from_admet
    } else {
      x$smiles_from_gban
    }

    if (is.null(smiles_for_swiss)) {
      stop("No SMILES vector was found from previous steps.")
    }

    if (nrow(data_swiss) != length(smiles_for_swiss)) {
      stop("nrow(data_swiss) != length(smiles_for_swiss).")
    }

    data_swiss <- dplyr::mutate(
      data_swiss,
      raw_smiles = smiles_for_swiss,
      .after = 1
    )

    colnames(data_swiss) <- formal_name(colnames(data_swiss))

    data_eval <- gBanFuns$evaluate_swiss_druglikeness(
      data = data_swiss,
      lipinski_min_pass = lipinski_min_pass,
      bioavailability_min = bioavailability_min,
      tpsa_max = tpsa_max,
      mw_max = mw_max,
      mw_max_relaxed = mw_max_relaxed,
      logp_max = logp_max,
      logp_col = logp_col
    )

    data_flow <- gBanFuns$summarize_swiss_selection_flow(data_eval)

    data_decision <- gBanFuns$choose_swiss_selection_rule(
      data_eval = data_eval,
      method = method,
      min_keep = min_keep
    )

    rule_col <- data_decision$rule_col[1L]
    data_keep <- data_eval[data_eval[[rule_col]], , drop = FALSE]

    x$swiss_skipped <- FALSE
    x$swissAdme <- tibble::as_tibble(data_swiss)
    x$swissAdme_eval <- tibble::as_tibble(data_eval)
    x$swiss_selection_flow <- tibble::as_tibble(data_flow)
    x$swiss_selection_decision <- tibble::as_tibble(data_decision)
    x$smiles_keep <- x$smiles_from_swiss <- data_keep$raw_smiles

    t.candidate_swissAdme <- gBanFuns$make_swiss_candidate_table(
      data_eval = data_keep,
      include_smiles = TRUE
    )

    x$t_candidate_swissAdme <- tibble::as_tibble(t.candidate_swissAdme)

    t.swiss_flow <- set_lab_legend(
      data_flow,
      glue::glue("{x@sig} SwissADME drug-likeness screening statistics"),
      glue::glue(
        "SwissADME 成药性筛选统计表|||",
        "该表汇总候选化合物在 SwissADME 成药性评价中的逐层保留情况，",
        "包括 Lipinski 规则、生物利用度评分及 TPSA 条件下的候选化合物数量变化。"
      )
    )

    x <- tablesAdd(
      x,
      t.swiss_flow = t.swiss_flow
    )

    if (isTRUE(add_candidate_table)) {
      t.candidate_report <- gBanFuns$make_swiss_candidate_table(
        data_eval = data_keep,
        include_smiles = FALSE
      )

      t.candidate_report <- set_lab_legend(
        t.candidate_report,
        glue::glue("{x@sig} candidates retained by SwissADME drug-likeness evaluation"),
        glue::glue(
          "SwissADME 成药性评价后保留的候选化合物|||",
          "该表展示经 SwissADME 成药性评价后保留的候选化合物及其主要分子性质，",
          "包括分子量、氢键受体数、氢键供体数、脂水分配系数、TPSA、Lipinski 通过项数、",
          "生物利用度评分及可合成性评分。"
        )
      )

      x <- tablesAdd(
        x,
        t.candidate_swissAdme = t.candidate_report
      )
    }

    text_rule <- gBanFuns$format_swiss_rule_text(data_decision)
    text_flow <- gBanFuns$format_swiss_flow_text(data_flow)
    n_input <- nrow(data_eval)
    n_keep <- nrow(data_keep)
    x <- methodAdd(
      x,
      glue::glue(
        "将 {context_input$stage_label}保留的候选化合物输入 SwissADME 平台（<https://www.swissadme.ch/>）进行成药性评价。",
        "根据 SwissADME 输出的分子量、氢键受体数、氢键供体数、脂水分配系数、拓扑极性表面积（TPSA）",
        "及生物利用度评分等参数，评估候选化合物的理化性质和类药性。",
        "Lipinski 规则按 MW≤{mw_max}、HBA≤10、HBD≤5 和 LogP≤{logp_max} 统计通过项，",
        "保留至少 {lipinski_min_pass} 项满足 Lipinski 条件且生物利用度评分大于 {bioavailability_min} 的候选化合物。",
        "{text_rule}"
      )
    )

    x <- snapAdd(
      x,
      glue::glue(
        "SwissADME 成药性评价前共有 {n_input} 个候选化合物。",
        "{text_flow}",
        "本分析最终保留 {n_keep} 个候选化合物。"
      )
    )

    return(x)
  })

setMethod("step7", signature = c(x = "job_gBan"),
  function(x,
    require_pubchem = TRUE,
    annotation_source = c("auto", "custom", "pubchem"),
    identity_types = c(
      "same_stereo_isotope",
      "same_stereo",
      "same_isotope",
      "same_connectivity"
    ),
    use_similarity = FALSE,
    similarity_threshold = 99L,
    use_cache = TRUE,
    overwrite = FALSE,
    dir_save = create_job_cache_dir(x),
    sleep = .2
  )
  {
    step_message("PubChem compound annotation.")
    expect_package("jsonlite")
    annotation_source <- match.arg(annotation_source)

    if (is.null(x$smiles_keep)) {
      stop("No SMILES vector was found from previous steps.")
    }

    smiles <- unique(as.character(x$smiles_keep))
    smiles <- smiles[!is.na(smiles) & smiles != ""]

    if (!length(smiles)) {
      stop("No valid SMILES was found from previous steps.")
    }

    context <- gBanFuns$resolve_candidate_filter_context(x)

    data_input <- data.frame(
      SMILES = smiles,
      stringsAsFactors = FALSE
    )

    if (!is.null(x$t_candidate_swissAdme) &&
        all(c("SMILES", "Compound") %in% colnames(x$t_candidate_swissAdme))) {
      data_name <- as.data.frame(x$t_candidate_swissAdme, stringsAsFactors = FALSE)
      data_name <- data_name[!duplicated(data_name$SMILES), , drop = FALSE]
      data_input <- merge(
        data_input,
        data_name[, c("SMILES", "Compound"), drop = FALSE],
        by = "SMILES",
        all.x = TRUE
      )
    } else {
      data_input$Compound <- paste0("Candidate_", seq_len(nrow(data_input)))
    }

    data_custom_lookup <- gBanFuns$make_custom_compound_annotation(
      x = x,
      smiles = data_input$SMILES,
      compound = data_input$Compound
    )
    has_custom_annotation <- any(data_custom_lookup$Custom_annotation, na.rm = TRUE)

    annotation_source_used <- annotation_source
    if (identical(annotation_source_used, "auto")) {
      annotation_source_used <- if (isTRUE(has_custom_annotation)) "custom" else "pubchem"
    }

    if (identical(annotation_source_used, "custom") && !isTRUE(has_custom_annotation)) {
      stop("No custom compound annotation was found. Use `annotation_source = 'pubchem'` or provide custom compound metadata.")
    }

    if (identical(annotation_source_used, "custom")) {
      data_lookup <- data_custom_lookup
    } else {
      dir.create(dir_save, recursive = TRUE, showWarnings = FALSE)
      file_cache <- file.path(dir_save, "pubchem_smiles_annotation.rds")

      if (isTRUE(use_cache) && !isTRUE(overwrite) && file.exists(file_cache)) {
        data_lookup <- readRDS(file_cache)
      } else {
        data_lookup <- gBanFuns$annotate_pubchem_smiles(
          smiles = data_input$SMILES,
          compound = data_input$Compound,
          identity_types = identity_types,
          use_similarity = use_similarity,
          similarity_threshold = similarity_threshold,
          sleep = sleep
        )
        saveRDS(data_lookup, file_cache)
      }
    }

    data_lookup <- as.data.frame(data_lookup, stringsAsFactors = FALSE)

    if ("CID" %in% colnames(data_lookup)) {
      data_lookup$CID <- gBanFuns$normalize_pubchem_cid(data_lookup$CID)
    }

    data_lookup$PubChem_match <- !is.na(data_lookup$CID) & data_lookup$CID != ""

    if (!identical(annotation_source_used, "custom")) {
      file_syn_cache <- file.path(dir_save, "pubchem_selected_synonyms.rds")

      data_synonym <- gBanFuns$get_pubchem_selected_synonyms(
        data_lookup = data_lookup,
        use_cache = use_cache,
        overwrite = overwrite,
        file_cache = file_syn_cache,
        sleep = sleep
      )

      if (nrow(data_synonym) > 0L) {
        data_lookup <- merge(
          data_lookup,
          data_synonym,
          by = "CID",
          all.x = TRUE,
          suffixes = c("", "_selected")
        )

        data_lookup$Synonym <- ifelse(
          !is.na(data_lookup$Synonym_selected) & data_lookup$Synonym_selected != "",
          data_lookup$Synonym_selected,
          data_lookup$Synonym
        )

        data_lookup$Synonym_selected <- NULL
      }
    }

    if (identical(annotation_source_used, "custom")) {
      data_lookup_final <- data_lookup
    } else if (isTRUE(require_pubchem)) {
      data_lookup_final <- data_lookup[data_lookup$PubChem_match, , drop = FALSE]
    } else {
      data_lookup_final <- data_lookup
    }

    x$pubchem_lookup <- tibble::as_tibble(data_lookup)
    x$pubchem_lookup_final <- tibble::as_tibble(data_lookup_final)
    x$smiles_from_pubchem <- data_lookup_final$SMILES
    x$smiles_keep <- data_lookup_final$SMILES

    x$cids <- tibble::as_tibble(
      data_lookup[data_lookup$PubChem_match, c("SMILES", "CID"), drop = FALSE]
    )

    x$synos <- tibble::as_tibble(
      data_lookup_final[, c(
        "SMILES",
        "CID",
        "Synonym",
        "PubChem_match",
        "PubChem_strategy",
        "PubChem_n_cid"
      ), drop = FALSE]
    )

    if (identical(annotation_source_used, "custom")) {
      data_stat <- data.frame(
        Category = c(
          context$stage_label_en,
          "Annotated from custom compound database",
          "With PubChem CID in custom annotation",
          "Final retained candidates"
        ),
        Count = c(
          length(smiles),
          sum(data_lookup$Custom_annotation, na.rm = TRUE),
          sum(data_lookup$PubChem_match, na.rm = TRUE),
          nrow(data_lookup_final)
        ),
        stringsAsFactors = FALSE
      )
    } else {
      data_stat <- data.frame(
        Category = c(
          context$stage_label_en,
          "Matched with PubChem CID",
          "No PubChem CID matched",
          "Final retained candidates"
        ),
        Count = c(
          length(smiles),
          sum(data_lookup$PubChem_match, na.rm = TRUE),
          sum(!data_lookup$PubChem_match, na.rm = TRUE),
          nrow(data_lookup_final)
        ),
        stringsAsFactors = FALSE
      )
    }

    x$pubchem_annotation_stat <- tibble::as_tibble(data_stat)

    text_annotation_title <- if (identical(annotation_source_used, "custom")) {
      "自定义化合物库注释统计表"
    } else {
      "PubChem 化合物注释统计表"
    }
    text_annotation_caption <- if (identical(annotation_source_used, "custom")) {
      glue::glue(
        "该表汇总 {context$stage_label}候选化合物从自定义化合物库继承名称、PubChem CID 和来源注释的情况，",
        "并统计最终保留的候选化合物数量。"
      )
    } else {
      glue::glue(
        "该表汇总 {context$stage_label}候选化合物在 PubChem 数据库中的结构身份检索情况，",
        "包括成功匹配 PubChem CID 的化合物数量及最终保留的候选化合物数量。"
      )
    }

    t.pubchem_stat <- set_lab_legend(
      data_stat,
      glue::glue("{x@sig} compound annotation statistics"),
      glue::glue(
        "{text_annotation_title}|||",
        "{text_annotation_caption}"
      )
    )

    t.swiss_for_final <- if (isTRUE(context$has_swiss)) {
      if (!is.null(x$t_candidate_swissAdme)) {
        x$t_candidate_swissAdme
      } else if (!is.null(x@tables$step6$t.candidate_swissAdme)) {
        x@tables$step6$t.candidate_swissAdme
      } else {
        NULL
      }
    } else {
      NULL
    }

    t.admet_for_final <- if (isTRUE(context$has_admet)) {
      if (!is.null(x$t_candidate_admet)) {
        x$t_candidate_admet
      } else if (!is.null(x@tables$step5$t.candidate_admet)) {
        x@tables$step5$t.candidate_admet
      } else {
        NULL
      }
    } else {
      NULL
    }

    t.final_candidates <- gBanFuns$merge_final_candidate_tables(
      data_pubchem = data_lookup_final,
      data_swiss = t.swiss_for_final,
      data_admet = t.admet_for_final
    )

    t.final_report <- gBanFuns$make_final_candidate_report_table(
      data = t.final_candidates,
      include_smiles = FALSE
    )

    text_final_caption <- if (identical(annotation_source_used, "custom")) {
      glue::glue(
        "该表展示经 {context$pipeline_label}后保留的候选化合物，",
        "并继承自定义化合物库中的化合物名称、PubChem CID、InChIKey 及来源证据等注释。"
      )
    } else {
      glue::glue(
        "该表展示经 {context$pipeline_label}后，",
        "进一步通过 PubChem 结构身份检索确认并保留的候选化合物，",
        "表中补充 PubChem CID 及对应化合物名称注释。"
      )
    }

    t.final_report <- set_lab_legend(
      t.final_report,
      glue::glue("{x@sig} candidate compounds retained after {context$pipeline_label_en}"),
      glue::glue(
        "最终候选化合物汇总表|||",
        "{text_final_caption}"
      )
    )

    max_transpose_candidates <- 6L
    t.final_candidates_mutate <- gBanFuns$make_final_candidate_transposed_table(
      data = t.final_report,
      name_col = "Synonym",
      id_col = "CID",
      max_candidates_per_block = max_transpose_candidates
    )

    text_mutate_caption <- if (nrow(t.final_report) > max_transpose_candidates) {
      glue::glue(
        "该表以分组转置形式展示最终候选化合物的核心注释及评价指标；",
        "每组最多展示 {max_transpose_candidates} 个候选化合物，",
        "通过 Candidate_ID 行对应完整候选化合物编号，以避免 Word 表格横向列数过多。"
      )
    } else {
      "该表以转置形式展示最终候选化合物的核心注释及评价指标，便于在报告正文中展示。"
    }

    t.final_candidates_mutate <- set_lab_legend(
      t.final_candidates_mutate,
      glue::glue("{x@sig} candidate compounds retained after {context$pipeline_label_en} grouped transposition"),
      glue::glue(
        "最终候选化合物分组转置汇总表|||",
        "{text_mutate_caption}"
      )
    )

    x$t_final_candidates <- tibble::as_tibble(t.final_candidates)
    x$t_final_candidates_report <- tibble::as_tibble(t.final_report)
    x$t_final_candidates_mutate <- tibble::as_tibble(t.final_candidates_mutate)

    x$.feature_genes <- gBanFuns$make_final_gene_feature(
      x = x,
      smiles = data_lookup_final$SMILES
    )

    x <- tablesAdd(
      x,
      t.pubchem_stat = t.pubchem_stat,
      t.final_candidates = t.final_report,
      t.final_candidates_mutate = t.final_candidates_mutate
    )

    n_input <- length(smiles)
    n_match <- sum(data_lookup$PubChem_match, na.rm = TRUE)
    n_unmatch <- n_input - n_match
    n_final <- nrow(data_lookup_final)

    if (identical(annotation_source_used, "custom")) {
      n_custom_annotated <- sum(data_lookup$Custom_annotation, na.rm = TRUE)
      x <- methodAdd(
        x,
        glue::glue(
          "对 {context$stage_label}保留的候选化合物进行注释整理。",
          "由于本流程使用自定义候选化合物库作为输入，优先继承自定义库中已整理的化合物名称、PubChem CID、InChIKey 及来源证据，",
          "不再将 PubChem 重新检索作为最终候选化合物的保留条件。",
          "该策略可避免同一结构因 SMILES 表达形式、立体化学标注或数据库检索差异导致的候选化合物误删。"
        )
      )

      x <- snapAdd(
        x,
        glue::glue(
          "{context$stage_label}共获得 {n_input} 个候选化合物；",
          "其中 {n_custom_annotated} 个化合物继承自定义化合物库注释，",
          "{n_match} 个化合物具有 PubChem CID。",
          "最终候选表保留 {n_final} 个化合物。"
        )
      )
    } else {
      text_pubchem_rule <- if (isTRUE(require_pubchem)) {
        "仅保留成功匹配 PubChem CID 的条目作为最终候选化合物。"
      } else {
        "PubChem 检索用于化合物注释，未匹配 CID 的结构条目在该注释层中仍予保留。"
      }
      text_pubchem_snap <- if (isTRUE(require_pubchem)) {
        glue::glue("PubChem 可注释条目作为最终候选化合物，共保留 {n_final} 个化合物。")
      } else {
        glue::glue("PubChem 注释后保留全部结构条目，共保留 {n_final} 个化合物。")
      }

      x <- methodAdd(
        x,
        glue::glue(
          "对 {context$stage_label}保留的候选化合物进行 PubChem 数据库注释。",
          "检索以 SMILES 结构式为输入，首先进行严格结构身份匹配；对于未直接匹配的条目，",
          "进一步采用同立体化学、同同位素及相同连接关系等不同身份层级进行检索，",
          "以降低 SMILES 表达形式、立体化学标注或同位素标注差异造成的漏检。",
          "成功匹配的条目获取 PubChem CID，并结合 PubChem 记录标题、同义名及 IUPAC 名称筛选代表性化合物名称；当缺少简短通用名称时，以 IUPAC 名称作为候选化合物名称兜底。",
          "{text_pubchem_rule}"
        )
      )

      x <- snapAdd(
        x,
        glue::glue(
          "{context$stage_label}共获得 {n_input} 个候选化合物；",
          "经 PubChem 结构身份检索，{n_match} 个化合物成功匹配 PubChem CID，",
          "{n_unmatch} 个化合物未匹配到 PubChem CID。",
          "{text_pubchem_snap}"
        )
      )
    }

    vec_feature <- x$synos$Synonym
    vec_feature <- ifelse(
      is.na(vec_feature) | vec_feature == "",
      paste0("CID", x$synos$CID),
      vec_feature
    )

    feature(x) <- as_feature(
      vec_feature,
      "候选化合物",
      nature = "compounds"
    )

    return(x)
  })

setMethod("asjob_vina", signature = c(x = "job_gBan"),
  function(x,
    require_pubchem = TRUE
  )
  {
    if (is.null(x$pubchem_lookup_final)) {
      stop("Please run step7 before asjob_vina().")
    }

    cpds <- as.data.frame(x$pubchem_lookup_final, stringsAsFactors = FALSE)

    if (isTRUE(require_pubchem)) {
      cpds <- cpds[!is.na(cpds$CID) & cpds$CID != "", , drop = FALSE]
    }

    if (!nrow(cpds)) {
      stop("No candidate compound with PubChem CID was found for vina workflow.")
    }

    layout <- dplyr::filter(
      x$res_graphBan,
      SMILES %in% cpds$SMILES
    )

    layout <- merge(
      layout,
      cpds[, c("SMILES", "CID", "Synonym"), drop = FALSE],
      by = "SMILES",
      all.x = TRUE
    )

    layout <- dplyr::select(
      layout,
      Synonym,
      hgnc_symbol,
      CID
    )

    layout <- dplyr::mutate(
      layout,
      hgnc_symbol = as.character(hgnc_symbol)
    )

    layout <- dplyr::distinct(layout)
    x <- job_vina(.layout = layout)
    x <- snapAdd(x, "对候选化合物与对应靶点进行 AutoDock Vina 分子对接。")
    return(x)
  })


if (!exists("gBanFuns")) {
  gBanFuns <- new.env(parent = emptyenv())
}


gBanFuns$resolve_custom_compound_db <- function(db_custom,
  db_custom_name = "Custom compound set", require_smiles = TRUE)
{
  if (is.null(db_custom)) {
    stop("`db_custom` must be provided when `db = 'custom'`.")
  }

  as_chr <- function(x) {
    if (is.null(x)) {
      return(rep(NA_character_, 0L))
    }
    if (is.numeric(x)) {
      return(trimws(format(x, scientific = FALSE, trim = TRUE)))
    }
    trimws(as.character(x))
  }

  collapse_unique <- function(x) {
    x <- as_chr(x)
    x <- unique(x[!is.na(x) & nzchar(x)])
    if (!length(x)) {
      return(NA_character_)
    }
    paste(x, collapse = "; ")
  }

  find_col <- function(data, candidates) {
    cols <- colnames(data)
    cols_lower <- tolower(gsub("[ ._-]+", "", cols))
    candidates_lower <- tolower(gsub("[ ._-]+", "", candidates))
    idx <- match(candidates_lower, cols_lower)
    idx <- idx[!is.na(idx)]
    if (!length(idx)) {
      return(NA_character_)
    }
    cols[idx[1L]]
  }

  collection <- NULL
  data_custom <- NULL

  if (is.character(db_custom) && length(db_custom) == 1L && file.exists(db_custom)) {
    data_custom <- ftibble(db_custom)
    data_custom <- as.data.frame(data_custom, stringsAsFactors = FALSE)
  } else if (inherits(db_custom, "job_herbsCollection")) {
    collection <- tryCatch(db_custom$collection, error = function(e) NULL)
    if (is.null(collection)) {
      stop("`db_custom` is a job_herbsCollection object, but `db_custom$collection` was not found.")
    }
  } else if (inherits(db_custom, "herbs_collection") ||
      (is.list(db_custom) && !is.null(db_custom$compound_unique))) {
    collection <- db_custom
  } else if (is.data.frame(db_custom)) {
    data_custom <- as.data.frame(db_custom, stringsAsFactors = FALSE)
  } else {
    stop("`db_custom` must be a job_herbsCollection object, a herbs_collection list, a data.frame, or a readable table path.")
  }

  if (!is.null(collection)) {
    if (is.null(collection$compound_unique)) {
      stop("`db_custom` collection does not contain `compound_unique`.")
    }

    data_custom <- as.data.frame(collection$compound_unique, stringsAsFactors = FALSE)

    if (!is.null(collection$herb_compound) &&
        "compound_key" %in% colnames(data_custom)) {
      data_rel <- as.data.frame(collection$herb_compound, stringsAsFactors = FALSE)

      if ("compound_key" %in% colnames(data_rel)) {
        if (!"query_herb" %in% colnames(data_rel)) {
          data_rel$query_herb <- NA_character_
        }
        if (!"herb_latin_name" %in% colnames(data_rel)) {
          data_rel$herb_latin_name <- NA_character_
        }
        if (!"source" %in% colnames(data_rel)) {
          data_rel$source <- NA_character_
        }

        data_rel_split <- split(data_rel, data_rel$compound_key)
        data_rel_sum <- do.call(rbind, lapply(names(data_rel_split), function(key) {
          data_one <- data_rel_split[[key]]
          data.frame(
            compound_key = key,
            herb = collapse_unique(data_one$query_herb),
            latin_name = collapse_unique(data_one$herb_latin_name),
            evidence_sources = collapse_unique(data_one$source),
            stringsAsFactors = FALSE
          )
        }))

        data_custom <- merge(
          data_custom,
          data_rel_sum,
          by = "compound_key",
          all.x = TRUE
        )
      }
    }
  }

  if (is.null(data_custom) || !nrow(data_custom)) {
    stop("No custom compound record was available.")
  }

  col_smiles <- find_col(data_custom, c(
    "smiles", "SMILES", "canonical_smiles", "isomeric_smiles",
    "Canonical SMILES", "Isomeric SMILES"
  ))

  if (is.na(col_smiles)) {
    stop("The custom compound database must contain a SMILES column.")
  }

  col_name <- find_col(data_custom, c(
    "compound_name", "Compound name", "Compound", "name", "Name",
    "Synonym", "compound"
  ))
  col_cid <- find_col(data_custom, c(
    "pubchem_cid", "PubChem CID", "PubChem_CID", "CID", "cid"
  ))
  col_inchikey <- find_col(data_custom, c(
    "inchikey", "InChIKey", "InChI Key"
  ))
  col_key <- find_col(data_custom, c(
    "compound_key", "Compound key", "compound_id", "Compound ID"
  ))
  col_herb <- find_col(data_custom, c(
    "herb", "Herb", "query_herb", "herb_cn_name", "Chinese herb"
  ))
  col_latin <- find_col(data_custom, c(
    "latin_name", "Latin name", "herb_latin_name", "Latin.Name"
  ))
  col_evidence <- find_col(data_custom, c(
    "evidence_sources", "Evidence sources", "source_list", "source",
    "Source"
  ))

  n <- nrow(data_custom)
  data_out <- data.frame(
    smiles = as_chr(data_custom[[col_smiles]]),
    compound_name = if (!is.na(col_name)) as_chr(data_custom[[col_name]]) else rep(NA_character_, n),
    pubchem_cid = if (!is.na(col_cid)) as_chr(data_custom[[col_cid]]) else rep(NA_character_, n),
    inchikey = if (!is.na(col_inchikey)) as_chr(data_custom[[col_inchikey]]) else rep(NA_character_, n),
    compound_key = if (!is.na(col_key)) as_chr(data_custom[[col_key]]) else rep(NA_character_, n),
    herb = if (!is.na(col_herb)) as_chr(data_custom[[col_herb]]) else rep(NA_character_, n),
    latin_name = if (!is.na(col_latin)) as_chr(data_custom[[col_latin]]) else rep(NA_character_, n),
    evidence_sources = if (!is.na(col_evidence)) as_chr(data_custom[[col_evidence]]) else rep(NA_character_, n),
    stringsAsFactors = FALSE
  )

  data_out$smiles <- trimws(data_out$smiles)
  data_out <- data_out[!is.na(data_out$smiles) & nzchar(data_out$smiles), , drop = FALSE]

  if (isTRUE(require_smiles) && nrow(data_out) == 0L) {
    stop("No compound with a valid SMILES string was found in `db_custom`.")
  }

  data_split <- split(data_out, data_out$smiles)
  data_out <- do.call(rbind, lapply(names(data_split), function(smiles_i) {
    data_one <- data_split[[smiles_i]]
    data.frame(
      smiles = smiles_i,
      compound_name = collapse_unique(data_one$compound_name),
      pubchem_cid = collapse_unique(data_one$pubchem_cid),
      inchikey = collapse_unique(data_one$inchikey),
      compound_key = collapse_unique(data_one$compound_key),
      herb = collapse_unique(data_one$herb),
      latin_name = collapse_unique(data_one$latin_name),
      evidence_sources = collapse_unique(data_one$evidence_sources),
      db_custom_name = db_custom_name,
      stringsAsFactors = FALSE
    )
  }))

  rownames(data_out) <- NULL
  tibble::as_tibble(data_out)
}

gBanFuns$extract_swiss_numeric <- function(x)
{
  if (is.numeric(x)) {
    return(as.numeric(x))
  }

  x <- as.character(x)
  x[x %in% c("", "n/d", "N/D", "NA", "NaN")] <- NA_character_
  x <- gsub(",", "", x, fixed = TRUE)
  x <- gsub("[^0-9eE+\\.\\-]+", "", x, perl = TRUE)
  suppressWarnings(as.numeric(x))
}

gBanFuns$resolve_first_existing_col <- function(data, cols)
{
  cols <- cols[cols %in% colnames(data)]

  if (!length(cols)) {
    return(NA_character_)
  }

  cols[1L]
}

gBanFuns$evaluate_swiss_druglikeness <- function(data,
  lipinski_min_pass = 3L,
  bioavailability_min = .1,
  tpsa_max = 140,
  mw_max = 500,
  mw_max_relaxed = 800,
  logp_max = 5,
  logp_col = c("Consensus_Log_P", "MLOGP", "XLOGP3", "WLOGP")
)
{
  data_eval <- as.data.frame(data, stringsAsFactors = FALSE)

  col_hba <- gBanFuns$resolve_first_existing_col(
    data_eval,
    c("X_H_bond_acceptors", "HBA", "H_bond_acceptors")
  )

  col_hbd <- gBanFuns$resolve_first_existing_col(
    data_eval,
    c("X_H_bond_donors", "HBD", "H_bond_donors")
  )

  col_logp <- gBanFuns$resolve_first_existing_col(
    data_eval,
    logp_col
  )

  col_bio <- gBanFuns$resolve_first_existing_col(
    data_eval,
    c("Bioavailability_Score", "Bioavailability_score")
  )

  data_eval$.MW <- gBanFuns$extract_swiss_numeric(data_eval$MW)
  data_eval$.HBA <- if (!is.na(col_hba)) {
    gBanFuns$extract_swiss_numeric(data_eval[[col_hba]])
  } else {
    NA_real_
  }

  data_eval$.HBD <- if (!is.na(col_hbd)) {
    gBanFuns$extract_swiss_numeric(data_eval[[col_hbd]])
  } else {
    NA_real_
  }

  data_eval$.LogP <- if (!is.na(col_logp)) {
    gBanFuns$extract_swiss_numeric(data_eval[[col_logp]])
  } else {
    NA_real_
  }

  data_eval$.TPSA <- if ("TPSA" %in% colnames(data_eval)) {
    gBanFuns$extract_swiss_numeric(data_eval$TPSA)
  } else {
    NA_real_
  }

  data_eval$.Bioavailability_Score <- if (!is.na(col_bio)) {
    gBanFuns$extract_swiss_numeric(data_eval[[col_bio]])
  } else {
    NA_real_
  }

  data_eval$.Synthetic_Accessibility <- if ("Synthetic_Accessibility" %in% colnames(data_eval)) {
    gBanFuns$extract_swiss_numeric(data_eval$Synthetic_Accessibility)
  } else {
    NA_real_
  }

  data_eval$.pass_mw <- !is.na(data_eval$.MW) & data_eval$.MW <= mw_max
  data_eval$.pass_mw_relaxed <- !is.na(data_eval$.MW) & data_eval$.MW <= mw_max_relaxed
  data_eval$.pass_hba <- !is.na(data_eval$.HBA) & data_eval$.HBA <= 10
  data_eval$.pass_hbd <- !is.na(data_eval$.HBD) & data_eval$.HBD <= 5
  data_eval$.pass_logp <- !is.na(data_eval$.LogP) & data_eval$.LogP <= logp_max
  data_eval$.pass_tpsa <- !is.na(data_eval$.TPSA) & data_eval$.TPSA <= tpsa_max
  data_eval$.pass_bioavailability <- !is.na(data_eval$.Bioavailability_Score) &
    data_eval$.Bioavailability_Score > bioavailability_min

  data_eval$.lipinski_pass_n <- rowSums(
    data_eval[, c(".pass_mw", ".pass_hba", ".pass_hbd", ".pass_logp")],
    na.rm = TRUE
  )

  data_eval$.pass_lipinski <- data_eval$.lipinski_pass_n >= lipinski_min_pass

  data_eval$.pass_strict <- data_eval$.pass_lipinski &
    data_eval$.pass_bioavailability &
    data_eval$.pass_tpsa

  data_eval$.pass_core <- data_eval$.pass_lipinski &
    data_eval$.pass_bioavailability

  data_eval$.pass_lipinski_only <- data_eval$.pass_lipinski

  data_eval
}

gBanFuns$summarize_swiss_selection_flow <- function(data_eval)
{
  data.frame(
    Stage = c(
      "SwissADME input",
      "Lipinski >= 3 criteria",
      "Lipinski + bioavailability score",
      "Lipinski + bioavailability score + TPSA"
    ),
    Criteria = c(
      "Candidate compounds submitted to SwissADME",
      "At least three of MW, HBA, HBD and LogP criteria were satisfied",
      "Lipinski criteria plus bioavailability score > 0.1",
      "Core rule plus TPSA <= 140 A2"
    ),
    Compounds = c(
      nrow(data_eval),
      sum(data_eval$.pass_lipinski, na.rm = TRUE),
      sum(data_eval$.pass_core, na.rm = TRUE),
      sum(data_eval$.pass_strict, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
}

gBanFuns$choose_swiss_selection_rule <- function(data_eval,
  method = c("adaptive", "strict", "core", "lipinski"),
  min_keep = 5L
)
{
  method <- match.arg(method)

  data_rule <- data.frame(
    method = c("strict", "core", "lipinski"),
    rule_col = c(".pass_strict", ".pass_core", ".pass_lipinski_only"),
    rule_label = c(
      "Lipinski + bioavailability score + TPSA",
      "Lipinski + bioavailability score",
      "Lipinski criteria"
    ),
    n_compound = c(
      sum(data_eval$.pass_strict, na.rm = TRUE),
      sum(data_eval$.pass_core, na.rm = TRUE),
      sum(data_eval$.pass_lipinski_only, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )

  if (method == "adaptive") {
    id <- which(data_rule$n_compound >= min_keep)

    if (!length(id)) {
      id <- which.max(data_rule$n_compound)
    } else {
      id <- id[1L]
    }
  } else {
    id <- match(method, data_rule$method)
  }

  data_rule[id, , drop = FALSE]
}

gBanFuns$make_swiss_candidate_table <- function(data_eval,
  include_smiles = TRUE
)
{
  data_table <- data.frame(
    SMILES = data_eval$raw_smiles,
    Compound = if ("Molecule" %in% colnames(data_eval)) {
      as.character(data_eval$Molecule)
    } else {
      paste0("Compound_", seq_len(nrow(data_eval)))
    },
    MW = round(data_eval$.MW, 2),
    HBA = data_eval$.HBA,
    HBD = data_eval$.HBD,
    LogP = round(data_eval$.LogP, 2),
    TPSA = round(data_eval$.TPSA, 2),
    Lipinski_pass = data_eval$.lipinski_pass_n,
    Bioavailability_Score = round(data_eval$.Bioavailability_Score, 2),
    Synthetic_Accessibility = round(data_eval$.Synthetic_Accessibility, 2),
    stringsAsFactors = FALSE
  )

  if (!isTRUE(include_smiles)) {
    data_table <- data_table[, setdiff(colnames(data_table), "SMILES"), drop = FALSE]
  }

  tibble::as_tibble(data_table)
}

gBanFuns$format_swiss_rule_text <- function(data_decision)
{
  method <- data_decision$method[1L]

  if (method == "strict") {
    return("同时将 TPSA≤140 Å² 作为筛选条件，用于保留极性表面积处于推荐范围内的候选化合物。")
  }

  if (method == "core") {
    return("TPSA 作为辅助药物性质指标进行统计与展示，不作为本轮硬性剔除条件。")
  }

  "生物利用度评分和 TPSA 作为辅助药物性质指标进行统计与展示，不作为本轮硬性剔除条件。"
}

gBanFuns$format_swiss_flow_text <- function(data_flow)
{
  n_lipinski <- data_flow$Compounds[data_flow$Stage == "Lipinski >= 3 criteria"]
  n_core <- data_flow$Compounds[data_flow$Stage == "Lipinski + bioavailability score"]
  n_strict <- data_flow$Compounds[data_flow$Stage == "Lipinski + bioavailability score + TPSA"]

  glue::glue(
    "其中 {n_lipinski} 个化合物满足至少 3 项 Lipinski 条件，",
    "{n_core} 个化合物同时满足 Lipinski 条件和生物利用度评分要求，",
    "{n_strict} 个化合物进一步满足 TPSA 条件。"
  )
}

gBanFuns$resolve_candidate_filter_context <- function(x)
{
  has_admet <- !isTRUE(x$admet_skipped) && !is.null(x$smiles_from_admet)
  has_swiss <- !isTRUE(x$swiss_skipped) && !is.null(x$smiles_from_swiss)

  admet_stage_label <- if (!is.null(x$admet_filter_config$stage_label)) {
    x$admet_filter_config$stage_label
  } else if (identical(x$admet_tool, "admetsar")) {
    "admetSAR 安全性初筛后"
  } else {
    "ADMETlab 毒性风险评估后"
  }

  admet_stage_label_en <- if (!is.null(x$admet_filter_config$stage_label_en)) {
    x$admet_filter_config$stage_label_en
  } else if (identical(x$admet_tool, "admetsar")) {
    "admetSAR-retained candidates"
  } else {
    "ADMETlab-retained candidates"
  }

  stage_label <- if (isTRUE(has_swiss)) {
    "SwissADME 成药性评价后"
  } else if (isTRUE(has_admet)) {
    admet_stage_label
  } else {
    "GraphBAN 预测后"
  }

  stage_label_en <- if (isTRUE(has_swiss)) {
    "SwissADME-retained candidates"
  } else if (isTRUE(has_admet)) {
    admet_stage_label_en
  } else {
    "GraphBAN-retained candidates"
  }

  vec_pipeline <- c(
    "GraphBAN 预测",
    if (isTRUE(has_admet)) admet_stage_label else NULL,
    if (isTRUE(has_swiss)) "SwissADME 成药性评价" else NULL
  )

  vec_pipeline_en <- c(
    "GraphBAN prediction",
    if (isTRUE(has_admet)) admet_stage_label_en else NULL,
    if (isTRUE(has_swiss)) "SwissADME drug-likeness evaluation" else NULL
  )

  list(
    has_admet = has_admet,
    has_swiss = has_swiss,
    stage_label = stage_label,
    stage_label_en = stage_label_en,
    pipeline_label = paste(vec_pipeline, collapse = "、"),
    pipeline_label_en = paste(vec_pipeline_en, collapse = " and ")
  )
}


gBanFuns$resolve_admet_smiles_identity <- function(data_admet, smiles,
  admet_tool = c("admetlab", "admetsar"))
{
  admet_tool <- match.arg(admet_tool)
  data_admet <- as.data.frame(data_admet, stringsAsFactors = FALSE)

  clean_smiles <- function(x) {
    x <- trimws(as.character(x))
    x[is.na(x)] <- NA_character_
    x
  }

  vec_upstream <- clean_smiles(smiles)
  vec_upstream <- vec_upstream[!is.na(vec_upstream) & nzchar(vec_upstream)]

  if (!length(vec_upstream)) {
    stop("No upstream SMILES was provided for ADMET identity mapping.")
  }

  col_admet_smiles <- NA_character_
  if ("SMILES" %in% colnames(data_admet)) {
    col_admet_smiles <- "SMILES"
  } else if ("smiles" %in% colnames(data_admet)) {
    col_admet_smiles <- "smiles"
  } else if ("raw_smiles" %in% colnames(data_admet)) {
    col_admet_smiles <- "raw_smiles"
  }

  if (is.na(col_admet_smiles)) {
    stop('Column "raw_smiles", "smiles" or "SMILES" was not found in ADMET table.')
  }

  data_admet$admet_smiles <- clean_smiles(data_admet[[col_admet_smiles]])

  n_upstream <- length(vec_upstream)
  n_admet <- nrow(data_admet)

  if (admet_tool == "admetsar" && n_admet != n_upstream) {
    stop(glue::glue(
      "admetSAR output row number changed: upstream SMILES = {n_upstream}, ",
      "admetSAR rows = {n_admet}. Please provide an ADMET table containing ",
      "the original raw_smiles column for reliable identity mapping."
    ))
  }

  if ("raw_smiles" %in% colnames(data_admet)) {
    data_admet$raw_smiles <- clean_smiles(data_admet$raw_smiles)

    id_match <- match(vec_upstream, data_admet$raw_smiles)

    if (length(id_match) == n_upstream && all(!is.na(id_match))) {
      data_admet <- data_admet[id_match, , drop = FALSE]
      data_admet$raw_smiles <- vec_upstream
      mapping_mode <- "raw_smiles_match"
    } else if (admet_tool == "admetsar") {
      stop(glue::glue(
        "admetSAR raw_smiles could not be fully matched to upstream SMILES: ",
        "matched {sum(!is.na(id_match))}/{n_upstream}."
      ))
    } else {
      data_admet <- data_admet[data_admet$raw_smiles %in% vec_upstream, , drop = FALSE]
      mapping_mode <- "raw_smiles_subset"
    }
  } else {
    if (n_admet == n_upstream) {
      data_admet$raw_smiles <- vec_upstream
      mapping_mode <- "row_order"
    } else {
      id_keep <- data_admet$admet_smiles %in% vec_upstream
      data_admet <- data_admet[id_keep, , drop = FALSE]
      data_admet$raw_smiles <- data_admet$admet_smiles
      mapping_mode <- "admet_smiles_subset"
    }
  }

  vec_admet <- clean_smiles(data_admet$admet_smiles)
  vec_raw <- clean_smiles(data_admet$raw_smiles)

  n_position_same <- sum(vec_admet == vec_raw, na.rm = TRUE)
  n_admet_in_upstream <- sum(unique(vec_admet) %in% unique(vec_upstream))
  n_raw_in_upstream <- sum(unique(vec_raw) %in% unique(vec_upstream))

  data_diag <- data.frame(
    admet_tool = admet_tool,
    mapping_mode = mapping_mode,
    n_upstream_smiles = n_upstream,
    n_admet_rows = n_admet,
    n_output_rows = nrow(data_admet),
    n_admet_smiles_same_position = n_position_same,
    n_admet_smiles_overlap_upstream = n_admet_in_upstream,
    n_raw_smiles_overlap_upstream = n_raw_in_upstream,
    stringsAsFactors = FALSE
  )

  message(glue::glue(
    "ADMET SMILES identity check: upstream = {n_upstream}, ADMET rows = {n_admet}, ",
    "mapping = {mapping_mode}; ADMET SMILES identical to workflow key by row = ",
    "{n_position_same}/{nrow(data_admet)}, ADMET SMILES overlapping upstream key = ",
    "{n_admet_in_upstream}/{length(unique(vec_admet))}. Workflow raw_smiles keeps upstream SMILES."
  ))

  attr(data_admet, "smiles_diagnostics") <- data_diag
  data_admet
}

gBanFuns$extract_admet_numeric <- function(x)
{
  if (is.numeric(x)) {
    return(as.numeric(x))
  }

  x <- as.character(x)
  x[is.na(x)] <- NA_character_

  value <- rep(NA_character_, length(x))

  id_paren <- grepl("\\([-+]?[0-9]*\\.?[0-9]+\\)", x, perl = TRUE)

  if (any(id_paren, na.rm = TRUE)) {
    value[id_paren] <- sub(
      "^.*\\(([-+]?[0-9]*\\.?[0-9]+)\\).*$",
      "\\1",
      x[id_paren],
      perl = TRUE
    )
  }

  id_rest <- is.na(value) & !is.na(x)

  if (any(id_rest)) {
    value[id_rest] <- sub(
      "^.*?([-+]?[0-9]*\\.?[0-9]+).*$",
      "\\1",
      x[id_rest],
      perl = TRUE
    )

    value[id_rest & value == x] <- NA_character_
  }

  suppressWarnings(as.numeric(value))
}

gBanFuns$get_admet_filter_config <- function(data,
  admet_tool = c("admetlab", "admetsar"),
  route = c("systemic", "topical_nasal"),
  toxicity_cols = NULL,
  optional_toxicity_cols = NULL,
  include_optional_toxicity = FALSE
)
{
  admet_tool <- match.arg(admet_tool)
  route <- match.arg(route)

  if (admet_tool == "admetsar" && route == "topical_nasal") {
    hard_default <- c(
      "Respiratory_toxicity",
      "Eye_corrosion",
      "Skin_corrosion",
      "Ames",
      "Micronucleus",
      "Mouse_carcinogenicity_c",
      "Rat_carcinogenicity_c"
    )
    warning_default <- c(
      "Eye_irritation",
      "Skin_irritation",
      "Skin_sensitisation",
      "ADT",
      "Photoinduced_toxicity",
      "Phototoxicity_Photoirritation",
      "Photoallergy",
      "Repeated_dose_toxicity",
      "Reproductive_toxicity",
      "Mitochondrial_toxicity",
      "Hemolytic_toxicity",
      "DILI",
      "Nephrotoxicity",
      "hERG_1uM",
      "hERG_10uM",
      "hERG_30uM"
    )
    method_scope <- paste0(
      "admetSAR 3.0 提供多类 ADMET 与毒性端点预测；针对鼻腔外用给药场景，",
      "本分析将呼吸道毒性、眼/皮肤腐蚀、Ames、微核及啮齿动物致癌性等强安全性终点作为核心筛选指标，",
      "并将眼/皮肤刺激、皮肤致敏、急性经皮毒性、光毒性及系统暴露相关毒性作为风险注释指标。"
    )
    tool_label <- "admetSAR 3.0"
    tool_url <- "<https://lmmd.ecust.edu.cn/admetsar3/index.php>"
    route_label <- "鼻腔外用给药"
    stage_label <- "admetSAR 外用安全性初筛后"
    stage_label_en <- "admetSAR topical safety filter"
  } else {
    hard_default <- c(
      "DILI",
      "H-HT",
      "Carcinogenicity",
      "Ames",
      "Genotoxicity",
      "RPMI-8226"
    )
    warning_default <- c(
      "hERG",
      "hERG-10um",
      "ROA"
    )
    method_scope <- paste0(
      "ADMETlab 3.0 可基于分子结构预测多类 ADMET 相关性质及毒性终点；",
      "本分析重点纳入系统暴露相关毒性和结构毒性风险终点进行筛选。"
    )
    tool_label <- "ADMETlab 3.0"
    tool_url <- "<https://admetlab3.scbdd.com/>"
    route_label <- "系统给药"
    stage_label <- "ADMETlab 毒性风险评估后"
    stage_label_en <- "ADMETlab toxicity filter"
  }

  if (is.null(toxicity_cols)) {
    toxicity_cols <- hard_default
  }
  if (is.null(optional_toxicity_cols)) {
    optional_toxicity_cols <- warning_default
  }
  if (isTRUE(include_optional_toxicity)) {
    toxicity_cols <- unique(c(toxicity_cols, optional_toxicity_cols))
    optional_toxicity_cols <- setdiff(optional_toxicity_cols, toxicity_cols)
  }

  list(
    hard_cols = toxicity_cols[toxicity_cols %in% colnames(data)],
    warning_cols = optional_toxicity_cols[optional_toxicity_cols %in% colnames(data)],
    admet_tool = admet_tool,
    route = route,
    tool_label = tool_label,
    tool_url = tool_url,
    route_label = route_label,
    stage_label = stage_label,
    stage_label_en = stage_label_en,
    method_scope = method_scope
  )
}

gBanFuns$resolve_admet_toxicity_cols <- function(data,
  toxicity_cols = NULL,
  optional_toxicity_cols = NULL,
  include_optional_toxicity = FALSE
)
{
  if (is.null(toxicity_cols)) {
    toxicity_cols <- c(
      "DILI",
      "H-HT",
      "Carcinogenicity",
      "Ames",
      "Genotoxicity",
      "RPMI-8226"
    )
  }

  if (is.null(optional_toxicity_cols)) {
    optional_toxicity_cols <- c(
      "hERG",
      "hERG-10um",
      "ROA"
    )
  }

  if (isTRUE(include_optional_toxicity)) {
    toxicity_cols <- unique(c(toxicity_cols, optional_toxicity_cols))
  }

  toxicity_cols[toxicity_cols %in% colnames(data)]
}

gBanFuns$evaluate_admet_toxicity <- function(data,
  toxicity_cols,
  warning_cols = character(0L),
  cutoff = .7,
  warning_cutoff = .5,
  max_toxicity_flags = 0L,
  max_warning_flags = Inf
)
{
  data_eval <- as.data.frame(data, stringsAsFactors = FALSE)
  warning_cols <- warning_cols[warning_cols %in% colnames(data_eval)]

  for (col in unique(c(toxicity_cols, warning_cols))) {
    data_eval[[col]] <- gBanFuns$extract_admet_numeric(data_eval[[col]])
  }

  mat_toxicity <- as.matrix(data_eval[, toxicity_cols, drop = FALSE])
  mat_flag <- mat_toxicity >= cutoff
  mat_flag[is.na(mat_flag)] <- FALSE

  data_eval$n_toxicity_risk <- rowSums(mat_flag)
  data_eval$max_toxicity_risk <- apply(
    mat_toxicity,
    1L,
    function(x) {
      if (all(is.na(x))) {
        return(NA_real_)
      }

      max(x, na.rm = TRUE)
    }
  )

  data_eval$n_warning_risk <- 0L
  data_eval$max_warning_risk <- NA_real_

  if (length(warning_cols)) {
    mat_warning <- as.matrix(data_eval[, warning_cols, drop = FALSE])
    mat_warning_flag <- mat_warning >= warning_cutoff
    mat_warning_flag[is.na(mat_warning_flag)] <- FALSE
    data_eval$n_warning_risk <- rowSums(mat_warning_flag)
    data_eval$max_warning_risk <- apply(
      mat_warning,
      1L,
      function(x) {
        if (all(is.na(x))) {
          return(NA_real_)
        }

        max(x, na.rm = TRUE)
      }
    )
  }

  data_eval$pass_admet_toxicity <- data_eval$n_toxicity_risk <=
    max_toxicity_flags

  if (is.finite(max_warning_flags)) {
    data_eval$pass_admet_toxicity <- data_eval$pass_admet_toxicity &
      data_eval$n_warning_risk <= max_warning_flags
  }

  data_eval$major_risk_endpoints <- apply(
    mat_flag,
    1L,
    function(x) {
      cols <- toxicity_cols[which(x)]
      if (!length(cols)) {
        return("None")
      }
      paste(cols, collapse = "; ")
    }
  )

  data_eval$admet_decision <- ifelse(
    data_eval$pass_admet_toxicity,
    "Retained",
    "Excluded by core safety risk"
  )

  data_eval
}

gBanFuns$summarize_admet_toxicity_filter <- function(data_eval,
  toxicity_cols,
  cutoff = .7
)
{
  data_stat <- lapply(
    toxicity_cols,
    function(col) {
      value <- gBanFuns$extract_admet_numeric(data_eval[[col]])
      n_valid <- sum(!is.na(value))
      n_risk <- sum(value >= cutoff, na.rm = TRUE)
      percent_risk <- if (n_valid > 0L) {
        round(n_risk / n_valid * 100, 2L)
      } else {
        NA_real_
      }

      data.frame(
        endpoint = col,
        n_valid = n_valid,
        n_risk = n_risk,
        percent_risk = percent_risk,
        median_score = if (n_valid > 0L) stats::median(value, na.rm = TRUE) else NA_real_,
        max_score = if (n_valid > 0L) max(value, na.rm = TRUE) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  )

  data_stat <- do.call(rbind, data_stat)
  rownames(data_stat) <- NULL
  tibble::as_tibble(data_stat)
}

gBanFuns$resolve_admet_report_cols <- function(data,
  toxicity_cols,
  warning_cols = character(0L),
  admet_tool = "admetlab",
  route = "systemic"
)
{
  cols_aux <- c(
    "raw_smiles",
    "admet_decision",
    "pass_admet_toxicity",
    "n_toxicity_risk",
    "max_toxicity_risk",
    "n_warning_risk",
    "max_warning_risk",
    "major_risk_endpoints",
    toxicity_cols,
    warning_cols,
    "MW",
    "TPSA",
    "SlogP",
    "logS",
    "QED",
    "DILI",
    "Nephrotoxicity",
    "hERG",
    "hERG_1uM",
    "hERG_10uM",
    "hERG_30uM",
    "caco2",
    "hia",
    "PPB",
    "ROA",
    "LD50_oral",
    "Synth",
    "Lipinski"
  )

  unique(cols_aux[cols_aux %in% colnames(data)])
}

gBanFuns$format_admet_endpoint_text <- function(stat,
  cutoff = .7,
  max_items = 8L
)
{
  if (is.null(stat) || nrow(stat) == 0L) {
    return("")
  }

  stat <- as.data.frame(stat, stringsAsFactors = FALSE)
  stat <- stat[order(-stat$n_risk, stat$endpoint), , drop = FALSE]
  stat <- stat[seq_len(min(max_items, nrow(stat))), , drop = FALSE]

  text <- paste(
    glue::glue(
      "{stat$endpoint}: {stat$n_risk}/{stat$n_valid}"
    ),
    collapse = "；"
  )

  glue::glue("按终点统计，风险分数达到 {gBanFuns$format_cutoff(cutoff)} 的化合物数量为：{text}。")
}

gBanFuns$format_cutoff <- function(cutoff)
{
  cutoff <- as.numeric(cutoff)[1L]
  text <- formatC(cutoff, format = "f", digits = 2L)
  text <- sub("0+$", "", text)
  text <- sub("\\.$", "", text)
  text
}

gBanFuns$get_model_strategy_text <- function(model_strategy,
  n_model_total = 3L
)
{
  if (identical(model_strategy, "intersection")) {
    return(glue::glue(
      "仅保留全部 {n_model_total} 个 GraphBAN 模型均达到阈值的化合物-靶点互作对"
    ))
  }

  "保留任一 GraphBAN 模型达到阈值的化合物-靶点互作对"
}

gBanFuns$get_gene_strategy_text <- function(gene_strategy)
{
  if (identical(gene_strategy, "intersection")) {
    return("保留在不同靶点候选集合中共同出现的候选化合物")
  }

  "保留各靶点候选集合，并对重复 SMILES 进行去重"
}

gBanFuns$resolve_cutoff_grid <- function(cutoff = .95,
  min_cutoff = .5,
  max_cutoff = 1,
  cutoff_step = .05
)
{
  cutoff <- as.numeric(cutoff)[1L]
  min_cutoff <- as.numeric(min_cutoff)[1L]
  max_cutoff <- as.numeric(max_cutoff)[1L]
  cutoff_step <- abs(as.numeric(cutoff_step)[1L])

  if (!is.finite(cutoff) || !is.finite(min_cutoff) ||
      !is.finite(max_cutoff) || !is.finite(cutoff_step)) {
    stop("Cutoff parameters should be finite numeric values.")
  }

  if (min_cutoff < .5 || max_cutoff > 1 || min_cutoff > max_cutoff) {
    stop("The cutoff range should be within 0.5-1.0.")
  }

  cutoff <- min(max(cutoff, min_cutoff), max_cutoff)

  vec_cutoff <- seq(cutoff, min_cutoff, by = -cutoff_step)
  vec_cutoff <- unique(round(c(vec_cutoff, min_cutoff), 4L))
  vec_cutoff <- vec_cutoff[vec_cutoff >= min_cutoff & vec_cutoff <= max_cutoff]
  sort(vec_cutoff, decreasing = TRUE)
}

gBanFuns$read_graphBan_predictions <- function(files_res,
  combn,
  pattern = "graphBan_res_",
  reRead = FALSE
)
{
  data_all <- pbapply::pblapply(
    files_res,
    function(file) {
      data_pred <- expect_local_data(
        "tmp", "gbanResRead", ftibble,
        list(files = file, select = "pred", fill = TRUE),
        rerun = reRead
      )

      if (nrow(data_pred) != nrow(combn)) {
        stop("nrow(data_pred) != nrow(combn).")
      }

      model <- basename(file)
      model <- sub(pattern, "", model)
      model <- sub("\\..*$", "", model)

      data_pred <- dplyr::bind_cols(
        data.frame(
          model = rep(model, nrow(data_pred)),
          stringsAsFactors = FALSE
        ),
        data_pred,
        combn
      )

      data_pred
    }
  )

  data_all <- dplyr::bind_rows(data_all)
  data_all$model <- as.character(data_all$model)
  data_all$hgnc_symbol <- as.character(data_all$hgnc_symbol)
  data_all$SMILES <- as.character(data_all$SMILES)
  data_all$pred <- as.numeric(data_all$pred)
  data_all
}

gBanFuns$get_graphBan_pair_table <- function(data_pred,
  cutoff,
  model_strategy = c("intersection", "union")
)
{
  model_strategy <- match.arg(model_strategy)

  data_pass <- dplyr::filter(data_pred, pred >= cutoff)

  if (nrow(data_pass) == 0L) {
    data_empty <- data_pred[0L, , drop = FALSE]
    data_empty$model_strategy <- model_strategy
    data_empty$cutoff <- cutoff
    return(data_empty)
  }

  n_model_total <- length(unique(data_pred$model))

  data_pair <- dplyr::group_by(
    data_pass,
    hgnc_symbol,
    SMILES
  )

  data_pair <- dplyr::summarise(
    data_pair,
    pred_max = max(pred, na.rm = TRUE),
    pred_mean = mean(pred, na.rm = TRUE),
    n_models_pass = dplyr::n_distinct(model),
    models_pass = paste(sort(unique(model)), collapse = "; "),
    id = dplyr::first(id),
    Protein = dplyr::first(Protein),
    Y = dplyr::first(Y),
    .groups = "drop"
  )

  if (model_strategy == "intersection") {
    data_pair <- dplyr::filter(data_pair, n_models_pass >= n_model_total)
  }

  data_pair$model_strategy <- model_strategy
  data_pair$cutoff <- cutoff
  data_pair
}

gBanFuns$apply_graphBan_gene_strategy <- function(data_pair,
  genes,
  gene_strategy = c("intersection", "respective")
)
{
  gene_strategy <- match.arg(gene_strategy)

  if (nrow(data_pair) == 0L) {
    data_pair$gene_strategy <- gene_strategy
    return(data_pair)
  }

  if (gene_strategy == "intersection") {
    data_gene <- dplyr::group_by(data_pair, SMILES)
    data_gene <- dplyr::summarise(
      data_gene,
      n_genes_pass = dplyr::n_distinct(hgnc_symbol),
      genes_pass = paste(sort(unique(hgnc_symbol)), collapse = "; "),
      .groups = "drop"
    )

    vec_common <- data_gene$SMILES[data_gene$n_genes_pass >= length(genes)]
    data_pair <- dplyr::filter(data_pair, SMILES %in% vec_common)
  }

  data_pair$gene_strategy <- gene_strategy
  data_pair
}

gBanFuns$evaluate_graphBan_plans <- function(data_pred,
  cutoff_grid,
  target_min = 100L,
  target_max = 1000L,
  target_center = NULL,
  method_model = c("auto", "intersection", "union"),
  method_keep = c("auto", "all", "respective")
)
{
  method_model <- match.arg(method_model)
  method_keep <- match.arg(method_keep)

  if (is.null(target_center)) {
    target_center <- mean(c(target_min, target_max))
  }

  genes <- sort(unique(data_pred$hgnc_symbol))

  vec_model_strategy <- c("intersection", "union")
  vec_gene_strategy <- c("respective", "intersection")

  if (method_model != "auto") {
    vec_model_strategy <- method_model
  }

  if (method_keep == "all") {
    vec_gene_strategy <- "intersection"
  } else if (method_keep == "respective") {
    vec_gene_strategy <- "respective"
  }

  lst_plan <- list()
  i_plan <- 0L

  for (cutoff in cutoff_grid) {
    for (model_strategy in vec_model_strategy) {
      data_pair <- gBanFuns$get_graphBan_pair_table(
        data_pred = data_pred,
        cutoff = cutoff,
        model_strategy = model_strategy
      )

      for (gene_strategy in vec_gene_strategy) {
        data_selected <- gBanFuns$apply_graphBan_gene_strategy(
          data_pair = data_pair,
          genes = genes,
          gene_strategy = gene_strategy
        )

        i_plan <- i_plan + 1L
        lst_plan[[i_plan]] <- data.frame(
          cutoff = cutoff,
          model_strategy = model_strategy,
          gene_strategy = gene_strategy,
          n_row = nrow(data_selected),
          n_smiles = length(unique(data_selected$SMILES)),
          n_gene = length(unique(data_selected$hgnc_symbol)),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  data_plan <- dplyr::bind_rows(lst_plan)

  data_plan$in_target <- data_plan$n_smiles >= target_min &
    data_plan$n_smiles <= target_max

  data_plan$range_distance <- ifelse(
    data_plan$n_smiles < target_min,
    target_min - data_plan$n_smiles,
    ifelse(
      data_plan$n_smiles > target_max,
      data_plan$n_smiles - target_max,
      0
    )
  )

  data_plan$center_distance <- abs(data_plan$n_smiles - target_center)

  data_plan$model_rank <- ifelse(
    data_plan$model_strategy == "intersection",
    1L,
    2L
  )

  data_plan$gene_rank <- ifelse(
    data_plan$gene_strategy == "respective",
    1L,
    2L
  )

  data_in <- data_plan[data_plan$in_target, , drop = FALSE]
  data_out <- data_plan[!data_plan$in_target, , drop = FALSE]

  if (nrow(data_in) > 0L) {
    data_in <- data_in[
      order(
        data_in$model_rank,
        data_in$gene_rank,
        -data_in$cutoff,
        data_in$center_distance
      ),
      ,
      drop = FALSE
    ]
  }

  if (nrow(data_out) > 0L) {
    data_out <- data_out[
      order(
        data_out$range_distance,
        data_out$model_rank,
        data_out$gene_rank,
        -data_out$cutoff,
        data_out$center_distance
      ),
      ,
      drop = FALSE
    ]
  }

  data_plan <- rbind(data_in, data_out)
  rownames(data_plan) <- NULL
  data_plan
}

gBanFuns$select_graphBan_plan <- function(data_pred,
  data_plan
)
{
  data_decision <- data_plan[1L, , drop = FALSE]
  genes <- sort(unique(data_pred$hgnc_symbol))

  data_pair_model <- gBanFuns$get_graphBan_pair_table(
    data_pred = data_pred,
    cutoff = data_decision$cutoff[1L],
    model_strategy = data_decision$model_strategy[1L]
  )

  data_selected <- gBanFuns$apply_graphBan_gene_strategy(
    data_pair = data_pair_model,
    genes = genes,
    gene_strategy = data_decision$gene_strategy[1L]
  )

  data_selected <- dplyr::arrange(
    data_selected,
    hgnc_symbol,
    dplyr::desc(n_models_pass),
    dplyr::desc(pred_max),
    SMILES
  )

  list(
    data_decision = data_decision,
    data_pair_model = data_pair_model,
    data_selected = data_selected
  )
}

gBanFuns$summarize_graphBan_model_counts <- function(data_pred,
  cutoff
)
{
  data_pass <- dplyr::filter(data_pred, pred >= cutoff)

  if (nrow(data_pass) == 0L) {
    return(tibble::tibble(
      hgnc_symbol = character(0),
      model = character(0),
      n_candidate = integer(0),
      pred_median = numeric(0),
      pred_max = numeric(0)
    ))
  }

  data_stat <- dplyr::group_by(data_pass, hgnc_symbol, model)
  data_stat <- dplyr::summarise(
    data_stat,
    n_candidate = dplyr::n_distinct(SMILES),
    pred_median = stats::median(pred, na.rm = TRUE),
    pred_max = max(pred, na.rm = TRUE),
    .groups = "drop"
  )

  dplyr::arrange(data_stat, hgnc_symbol, model)
}

gBanFuns$replace_na_number <- function(x, value = 0)
{
  x[is.na(x)] <- value
  x
}

gBanFuns$summarize_graphBan_gene_selection <- function(data_pair_model,
  data_selected,
  n_model_total = NULL
)
{
  if (is.null(n_model_total) && nrow(data_pair_model) > 0L) {
    n_model_total <- max(data_pair_model$n_models_pass, na.rm = TRUE)
  }

  if (is.null(n_model_total) || !is.finite(n_model_total)) {
    n_model_total <- 0L
  }

  if (nrow(data_pair_model) > 0L) {
    data_model <- dplyr::group_by(data_pair_model, hgnc_symbol)
    data_model <- dplyr::summarise(
      data_model,
      n_candidate_after_model = dplyr::n_distinct(SMILES),
      n_all_model_supported_after_model = sum(
        n_models_pass >= n_model_total,
        na.rm = TRUE
      ),
      pred_max_after_model = max(pred_max, na.rm = TRUE),
      pred_mean_after_model = mean(pred_mean, na.rm = TRUE),
      .groups = "drop"
    )
  } else {
    data_model <- tibble::tibble(
      hgnc_symbol = character(0),
      n_candidate_after_model = integer(0),
      n_all_model_supported_after_model = integer(0),
      pred_max_after_model = numeric(0),
      pred_mean_after_model = numeric(0)
    )
  }

  if (nrow(data_selected) > 0L) {
    data_final <- dplyr::group_by(data_selected, hgnc_symbol)
    data_final <- dplyr::summarise(
      data_final,
      n_candidate_final = dplyr::n_distinct(SMILES),
      n_all_model_supported_final = sum(
        n_models_pass >= n_model_total,
        na.rm = TRUE
      ),
      pred_max_final = max(pred_max, na.rm = TRUE),
      pred_mean_final = mean(pred_mean, na.rm = TRUE),
      .groups = "drop"
    )
  } else {
    data_final <- tibble::tibble(
      hgnc_symbol = character(0),
      n_candidate_final = integer(0),
      n_all_model_supported_final = integer(0),
      pred_max_final = numeric(0),
      pred_mean_final = numeric(0)
    )
  }

  data_stat <- merge(
    as.data.frame(data_model),
    as.data.frame(data_final),
    by = "hgnc_symbol",
    all = TRUE
  )

  if (!nrow(data_stat)) {
    return(tibble::as_tibble(data_stat))
  }

  vec_count <- c(
    "n_candidate_after_model",
    "n_all_model_supported_after_model",
    "n_candidate_final",
    "n_all_model_supported_final"
  )

  for (col in vec_count) {
    if (col %in% colnames(data_stat)) {
      data_stat[[col]] <- gBanFuns$replace_na_number(data_stat[[col]], 0L)
    }
  }

  data_stat <- data_stat[
    order(-data_stat$n_candidate_final, data_stat$hgnc_symbol),
    ,
    drop = FALSE
  ]

  rownames(data_stat) <- NULL
  tibble::as_tibble(data_stat)
}

gBanFuns$summarize_graphBan_selection_flow <- function(data_pred,
  data_pair_model,
  data_selected,
  cutoff,
  model_strategy,
  gene_strategy
)
{
  data_pass <- dplyr::filter(data_pred, pred >= cutoff)

  .make_row <- function(stage,
    strategy,
    data,
    n_model_record = NA_integer_
  )
  {
    data.frame(
      stage = stage,
      strategy = strategy,
      n_model_record = n_model_record,
      n_target_compound_pair = length(unique(paste(data$hgnc_symbol, data$SMILES, sep = "|||"))),
      n_unique_compound = length(unique(data$SMILES)),
      n_target = length(unique(data$hgnc_symbol)),
      stringsAsFactors = FALSE
    )
  }

  data_flow <- rbind(
    .make_row(
      stage = "Input combinations",
      strategy = "All GraphBAN input pairs",
      data = data_pred,
      n_model_record = nrow(data_pred)
    ),
    .make_row(
      stage = "Probability threshold",
      strategy = glue::glue("Interaction probability >= {gBanFuns$format_cutoff(cutoff)}"),
      data = data_pass,
      n_model_record = nrow(data_pass)
    ),
    .make_row(
      stage = "Model-level integration",
      strategy = model_strategy,
      data = data_pair_model
    ),
    .make_row(
      stage = "Target-level integration",
      strategy = gene_strategy,
      data = data_selected
    )
  )

  rownames(data_flow) <- NULL
  tibble::as_tibble(data_flow)
}

gBanFuns$format_graphBan_flow_text <- function(data_flow)
{
  if (nrow(data_flow) < 4L) {
    return("")
  }

  data_threshold <- data_flow[data_flow$stage == "Probability threshold", , drop = FALSE]
  data_model <- data_flow[data_flow$stage == "Model-level integration", , drop = FALSE]
  data_final <- data_flow[data_flow$stage == "Target-level integration", , drop = FALSE]

  glue::glue(
    " 在阈值筛选阶段，共有 {data_threshold$n_model_record} 条模型预测记录达到阈值；",
    "经模型层面整合后得到 {data_model$n_target_compound_pair} 个化合物-靶点候选互作对，",
    "对应 {data_model$n_unique_compound} 个唯一化合物；",
    "经靶点层面整合后保留 {data_final$n_target_compound_pair} 个候选互作对。"
  )
}

gBanFuns$format_graphBan_gene_text <- function(data_gene_stat,
  max_show = 6L
)
{
  if (nrow(data_gene_stat) == 0L) {
    return("")
  }

  data_gene_stat <- data_gene_stat[
    order(data_gene_stat$hgnc_symbol),
    ,
    drop = FALSE
  ]

  data_show <- data_gene_stat[
    seq_len(min(max_show, nrow(data_gene_stat))),
    ,
    drop = FALSE
  ]

  text_item <- paste0(
    data_show$hgnc_symbol,
    "：",
    data_show$n_candidate_final,
    "个"
  )

  text <- paste(text_item, collapse = "；")

  if (nrow(data_gene_stat) > max_show) {
    text <- paste0(text, "；其余靶点见统计表")
  }

  glue::glue(" 各靶点最终候选化合物数量为 {text}。")
}

gBanFuns$plot_graphBan_model_counts <- function(data_stat)
{
  if (nrow(data_stat) == 0L) {
    stop("No data to plot.")
  }

  ggplot2::ggplot(
    data_stat,
    ggplot2::aes(
      x = hgnc_symbol,
      y = n_candidate,
      fill = model
    )
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = .75),
      width = .68
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = NULL,
      y = "Number of predicted compounds",
      fill = "Model"
    ) +
    ggplot2::theme_bw()
}


gBanFuns$prepare_graphBan_target_compound_network <- function(data_selected)
{
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop('Package "igraph" is required.')
  }

  data_selected <- as.data.frame(data_selected, stringsAsFactors = FALSE)

  vec_need <- c("hgnc_symbol", "SMILES", "pred_max", "pred_mean", "n_models_pass")
  vec_miss <- setdiff(vec_need, colnames(data_selected))

  if (length(vec_miss) > 0L) {
    stop(glue::glue(
      "Missing required column(s): {paste(vec_miss, collapse = ', ')}."
    ))
  }

  data_selected <- data_selected[
    !is.na(data_selected$hgnc_symbol) &
      data_selected$hgnc_symbol != "" &
      !is.na(data_selected$SMILES) &
      data_selected$SMILES != "",
    ,
    drop = FALSE
  ]

  if (nrow(data_selected) == 0L) {
    return(list(
      nodes = data.frame(stringsAsFactors = FALSE),
      edges = data.frame(stringsAsFactors = FALSE),
      graph = igraph::make_empty_graph(),
      compound_map = data.frame(stringsAsFactors = FALSE)
    ))
  }

  data_compound <- dplyr::group_by(data_selected, SMILES)
  data_compound <- dplyr::summarise(
    data_compound,
    n_targets = dplyr::n_distinct(hgnc_symbol),
    targets = paste(sort(unique(hgnc_symbol)), collapse = "; "),
    pred_max = max(pred_max, na.rm = TRUE),
    pred_mean = mean(pred_mean, na.rm = TRUE),
    n_models_pass_max = max(n_models_pass, na.rm = TRUE),
    .groups = "drop"
  )

  data_compound <- dplyr::arrange(
    data_compound,
    dplyr::desc(n_targets),
    dplyr::desc(pred_max),
    SMILES
  )

  data_compound$compound_id <- sprintf("C%03d", seq_len(nrow(data_compound)))

  data_edge <- merge(
    data_selected,
    as.data.frame(data_compound[, c("SMILES", "compound_id")]),
    by = "SMILES",
    all.x = TRUE
  )

  data_edges <- data.frame(
    from = data_edge$hgnc_symbol,
    to = data_edge$compound_id,
    edge_type = "Predicted interaction",
    SMILES = data_edge$SMILES,
    pred_max = data_edge$pred_max,
    pred_mean = data_edge$pred_mean,
    n_models_pass = data_edge$n_models_pass,
    stringsAsFactors = FALSE
  )

  data_edges <- unique(data_edges)

  data_target <- dplyr::group_by(data_edges, from)
  data_target <- dplyr::summarise(
    data_target,
    n_link = dplyr::n_distinct(to),
    pred_max = max(pred_max, na.rm = TRUE),
    .groups = "drop"
  )

  nodes_target <- data.frame(
    name = data_target$from,
    type = "Target",
    label = data_target$from,
    n_link = data_target$n_link,
    n_targets = NA_integer_,
    pred_max = data_target$pred_max,
    SMILES = NA_character_,
    stringsAsFactors = FALSE
  )

  nodes_compound <- data.frame(
    name = data_compound$compound_id,
    type = "Compound",
    label = data_compound$compound_id,
    n_link = data_compound$n_targets,
    n_targets = data_compound$n_targets,
    pred_max = data_compound$pred_max,
    SMILES = data_compound$SMILES,
    stringsAsFactors = FALSE
  )

  data_nodes <- rbind(nodes_target, nodes_compound)
  data_nodes$type <- factor(data_nodes$type, levels = c("Target", "Compound"))

  graph <- igraph::graph_from_data_frame(
    d = data_edges[, c("from", "to"), drop = FALSE],
    directed = FALSE,
    vertices = data_nodes
  )

  list(
    nodes = tibble::as_tibble(data_nodes),
    edges = tibble::as_tibble(data_edges),
    graph = graph,
    compound_map = tibble::as_tibble(data_compound)
  )
}


gBanFuns$plot_graphBan_target_compound_network <- function(res_network,
  layout = "fr",
  label_top_n = 40L,
  seed = 1L
)
{
  if (!requireNamespace("ggraph", quietly = TRUE)) {
    stop('Package "ggraph" is required.')
  }

  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop('Package "ggrepel" is required.')
  }

  graph <- res_network$graph

  if (igraph::gorder(graph) == 0L || igraph::gsize(graph) == 0L) {
    stop("No graph edge was found.")
  }

  set.seed(seed)
  data_layout <- ggraph::create_layout(graph, layout = layout)
  data_label <- as.data.frame(data_layout)

  data_label$.label_priority <- ifelse(
    data_label$type == "Target",
    1L,
    ifelse(data_label$type == "Compound" & data_label$n_targets > 1L, 2L, 3L)
  )

  data_label <- data_label[
    order(
      data_label$.label_priority,
      -data_label$n_link,
      -data_label$pred_max,
      data_label$name
    ),
    ,
    drop = FALSE
  ]

  data_label <- data_label[
    data_label$type == "Target" |
      (data_label$type == "Compound" & data_label$n_targets > 1L),
    ,
    drop = FALSE
  ]

  if (!is.null(label_top_n)) {
    data_label <- data_label[
      seq_len(min(label_top_n, nrow(data_label))),
      ,
      drop = FALSE
    ]
  }

  ggraph::ggraph(data_layout) +
    ggraph::geom_edge_link(
      color = "grey80",
      edge_width = .35,
      alpha = .75,
      show.legend = FALSE
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(
        color = type,
        shape = type,
        size = n_link
      ),
      alpha = .85
    ) +
    ggrepel::geom_label_repel(
      data = data_label,
      inherit.aes = FALSE,
      size = 3,
      max.overlaps = Inf,
      min.segment.length = 0,
      ggplot2::aes(
        x = x,
        y = y,
        label = label
      )
    ) +
    ggplot2::scale_shape_manual(
      values = c(Target = 16L, Compound = 17L)
    ) +
    ggplot2::scale_size_continuous(
      range = c(2.5, 8)
    ) +
    ggplot2::labs(
      color = "Node type",
      shape = "Node type",
      size = "Connections"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "right"
    )
}



.stat_rapp_table_by_fun <- function(data, levels, cover = "得到", fun = function(x) length(unique(x)))
{
  # n <- 0L
  # cols <- names(levels)
  # des <- unname(levels)
  # if (length(cols) != length(des)) {
  #   stop('length(cols) != length(des).')
  # }
  # rapp <- function(data, name) {
  #   if (n + 1 > length(cols)) {
  #     return("")
  #   }
  #   lst <- split(data, data[[cols[n + 1]]])
  #   n <<- n + 1L
  #   snap <- vapply(names(lst), FUN.VALUE = character(1), 
  #     function(name) {
  #       rapp(lst[[name]], name)
  #     })
  #   n <<- n - 1L
  #   glue::glue("{des[n + 1]}{name}{cover}{fun(data[[ cols[ n + 2 ] ]])}个唯一{des[ n + 2 ]}。")
  # }
  # rapp(data, "")
}

# <https://admetlab3.scbdd.com/explanation/#/>
# 
# Ames, f30, hia (intestine), BBB, Drug-induced liver injury (DILI)
# 0-0.3: excellent (green); 0.3-0.7: medium (yellow); 0.7-1.0: poor (red)
# Caco-2
# > -5.15: excellent (green); otherwise: poor (red)
# PPB
# ≤ 90%: excellent (green); otherwise: poor (red)
# "CYP1A2-inh"  "CYP2C19-inh" "CYP2C9-inh"  "CYP2D6-inh"  "CYP3A4-inh"  "CYP2B6-inh"  "CYP2C8-inh"
# Category 0: Non-substrate / Non-inhibitor; Category 1: substrate / inhibitor. 


get_remote_graphBan.huibang <- function(x, pattern = "graphBan_res_*",
  remote = "graphBan", remote_to = "remote", expect = 3L)
{
  dir_save <- x$dir_save
  dir.create(dir_save, FALSE)
  existFiles <- list.files(dir_save, pattern)
  if (!length(existFiles) || length(existFiles) < expect) {
    remote_dir <- glue::glue(
      "~/{s(guess_project(), '^[0-9]+_', '')}"
    )
    files <- pattern
    cmd_get <- glue::glue("scp {remote}:{remote_dir}/{files} {x$dir_save}")
    cdRun(cmd_get)
  }
  existFiles <- list.files(dir_save, pattern, full.names = TRUE)
  if (!is_sshfs_mount(remote_to)) {
    stop('!is_sshfs_mount(remote_to).')
  }
  toDir <- file.path(remote_to, x$dir_save)
  toDir_existFiles <- list.files(toDir, pattern)
  if (!length(toDir_existFiles) || length(toDir_existFiles) < expect) {
    message(glue::glue("Send to '{remote_to}'"))
    file.copy(existFiles, toDir)
  }
}

run_remote_graphBan.huibang <- function(x,
  remote = "graphBan", remote_from = "remote")
{
  remote_dir <- glue::glue(
    "~/{s(guess_project(), '^[0-9]+_', '')}"
  )
  if (!file.exists(x$file_combn)) {
    if (!is_sshfs_mount(remote_from)) {
      stop('!is_sshfs_mount(remote_from).')
    }
    message(glue::glue("Get `file_combn` from: '{remote_from}'"))
    dir.create(dirname(x$file_combn), FALSE)
    file_from <- file.path(remote_from, x$file_combn)
    file.copy(file_from, x$file_combn)
  }
  files <- paste(
    x$file_combn, file.path(.expath, "job_templ", "graphBan.sh")
  )
  cmd_prepare <- glue::glue("ssh {remote} 'mkdir {remote_dir}'")
  cmd_send <- glue::glue("scp {files} {remote}:{remote_dir}")
  cmd_run <- glue::glue("ssh {remote} 'cd {remote_dir} && nohup sh graphBan.sh > task.log 2>&1 &'")
  cdRun(cmd_prepare)
  cdRun(cmd_send)
  cdRun(cmd_run, wait = FALSE)
}

inBatches_get_compounds_weight.rdkit <- function(smiles_list, python = getOption("rdkit_python"))
{
  stop(glue::glue("..."))
  # if (!is.null(python)) {
  #   e(base::Sys.setenv(RETICULATE_PYTHON = python))
  #   e(reticulate::use_python(python))
  #   e(reticulate::py_config())
  # }
  # rdkit::import("rdkit.Chem")
  # rdkit::import("rdkit.Chem.Descriptors")
  # rdkit::import("rdkit.Chem.rdMolDescriptors")
}

inBatches_get_compounds_weight.rcdk <- function(smiles_list, 
  cl = NULL, ..., mem = 1000, dir_db = .prefix("smiles_compounds_weight", "db"))
{
  # IDs: your query, col: the ID column, res: results table
  # smiles_list <- head(smiles_list, n = 3000)
  dir.create(dir_db, FALSE)
  db <- new_db(file.path(dir_db, "compounds_weight.rdata"), "smiles")
  db <- not(db, smiles_list)
  query <- db@query
  if (length(query)) {
    groups <- grouping_vec2list(query, mem, TRUE)
    message("Total group: ", length(groups))
    cli::cli_alert_info("rcdk::parse.smiles")
    res <- pbapply::pblapply(groups, cl = cl,
      function(smiles) {
        res <- try(silent = TRUE, callr::r(
          get_compounds_weight, args = c(
            list(smiles_list = smiles, libPaths = .libPaths()), list(...)
          ),
          show = TRUE
          ))
        message(glue::glue("Group {attr(smiles, 'name')} finished."))
        if (!inherits(res, "try-error")) {
          res
        } else NULL
      })
    res <- frbind(res)
    db <- upd(db, res)
  }
  res <- dplyr::filter(db@db, smiles %in% !!smiles_list)
  res
}


gBanFuns$make_final_gene_feature <- function(x, smiles)
{
  smiles <- unique(as.character(smiles))
  smiles <- smiles[!is.na(smiles) & smiles != ""]

  if (!length(smiles)) {
    return(NULL)
  }

  data_gene <- NULL

  if (!is.null(x$res_graphBan)) {
    data_gene <- as.data.frame(x$res_graphBan, stringsAsFactors = FALSE)
  } else if (!is.null(x$res_graphBan_pair)) {
    data_gene <- as.data.frame(x$res_graphBan_pair, stringsAsFactors = FALSE)
  }

  if (is.null(data_gene) || !all(c("SMILES", "hgnc_symbol") %in% colnames(data_gene))) {
    return(NULL)
  }

  data_gene <- data_gene[data_gene$SMILES %in% smiles, , drop = FALSE]
  gene <- unique(as.character(data_gene$hgnc_symbol))
  gene <- gene[!is.na(gene) & gene != ""]

  if (!length(gene)) {
    return(NULL)
  }

  as_feature(
    gene,
    "最终候选化合物对应靶点基因",
    nature = "genes"
  )
}

gBanFuns$read_pubchem_json <- function(url)
{
  res <- tryCatch(
    suppressWarnings(jsonlite::fromJSON(url)),
    warning = function(w) NULL,
    error = function(e) NULL
  )

  res
}

gBanFuns$normalize_pubchem_cid <- function(x)
{
  if (is.null(x)) {
    return(character(0))
  }

  if (is.factor(x)) {
    x <- as.character(x)
  }

  if (is.numeric(x) || is.integer(x)) {
    x <- format(
      x,
      scientific = FALSE,
      trim = TRUE,
      digits = 22L
    )
  } else {
    x <- as.character(x)
  }

  x <- gsub("\\.0+$", "", x, perl = TRUE)
  x <- gsub("\\s+", "", x, perl = TRUE)
  x[x %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  x
}

gBanFuns$get_pubchem_cids_from_json <- function(x)
{
  if (is.null(x)) {
    return(character(0))
  }

  if (!is.null(x$IdentifierList$CID)) {
    cid <- gBanFuns$normalize_pubchem_cid(x$IdentifierList$CID)
    cid <- cid[!is.na(cid) & cid != ""]
    return(unique(cid))
  }

  character(0)
}

gBanFuns$query_pubchem_identity <- function(smiles,
  identity_type = "same_connectivity"
)
{
  smiles_url <- utils::URLencode(smiles, reserved = TRUE)

  url <- paste0(
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastidentity/smiles/",
    smiles_url,
    "/cids/JSON?identity_type=",
    identity_type
  )

  res <- gBanFuns$read_pubchem_json(url)
  gBanFuns$get_pubchem_cids_from_json(res)
}

gBanFuns$query_pubchem_similarity <- function(smiles,
  threshold = 99L
)
{
  smiles_url <- utils::URLencode(smiles, reserved = TRUE)

  url <- paste0(
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/smiles/",
    smiles_url,
    "/cids/JSON?Threshold=",
    threshold
  )

  res <- gBanFuns$read_pubchem_json(url)
  gBanFuns$get_pubchem_cids_from_json(res)
}

gBanFuns$get_pubchem_property_table <- function(cids,
  sleep = .2
)
{
  cids <- unique(gBanFuns$normalize_pubchem_cid(cids))
  cids <- cids[!is.na(cids) & cids != ""]

  if (!length(cids)) {
    return(data.frame(stringsAsFactors = FALSE))
  }

  property <- paste(
    c(
      "MolecularFormula",
      "MolecularWeight",
      "CanonicalSMILES",
      "IsomericSMILES",
      "InChIKey",
      "IUPACName"
    ),
    collapse = ","
  )

  url <- paste0(
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
    paste(cids, collapse = ","),
    "/property/",
    property,
    "/JSON"
  )

  res <- gBanFuns$read_pubchem_json(url)

  if (is.null(res$PropertyTable$Properties)) {
    return(data.frame(
      CID = cids,
      stringsAsFactors = FALSE
    ))
  }

  Sys.sleep(sleep)
  data <- as.data.frame(res$PropertyTable$Properties, stringsAsFactors = FALSE)

  if ("CID" %in% colnames(data)) {
    data$CID <- gBanFuns$normalize_pubchem_cid(data$CID)
  }

  data
}

gBanFuns$read_pubchem_record_title <- function(cid)
{
  cid <- gBanFuns$normalize_pubchem_cid(cid)
  cid <- cid[!is.na(cid) & cid != ""]

  if (!length(cid)) {
    return(NA_character_)
  }

  url <- paste0(
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/",
    cid[1L],
    "/JSON"
  )

  res <- gBanFuns$read_pubchem_json(url)

  title <- NA_character_

  if (!is.null(res$Record$RecordTitle)) {
    title <- as.character(res$Record$RecordTitle)[1L]
  }

  if (is.na(title) || title == "") {
    return(NA_character_)
  }

  title
}

gBanFuns$annotate_one_pubchem_smiles <- function(smiles,
  compound = NA_character_,
  identity_types = c(
    "same_stereo_isotope",
    "same_stereo",
    "same_isotope",
    "same_connectivity"
  ),
  use_similarity = FALSE,
  similarity_threshold = 99L,
  sleep = .2
)
{
  cids <- character(0)
  strategy <- NA_character_

  for (identity_type in identity_types) {
    cids <- gBanFuns$query_pubchem_identity(
      smiles = smiles,
      identity_type = identity_type
    )

    if (length(cids) > 0L) {
      strategy <- paste0("identity:", identity_type)
      break
    }

    Sys.sleep(sleep)
  }

  if (!length(cids) && isTRUE(use_similarity)) {
    cids <- gBanFuns$query_pubchem_similarity(
      smiles = smiles,
      threshold = similarity_threshold
    )

    if (length(cids) > 0L) {
      strategy <- paste0("similarity_2d:", similarity_threshold)
    }

    Sys.sleep(sleep)
  }

  cid <- if (length(cids) > 0L) {
    cids[1L]
  } else {
    NA_character_
  }

  data_property <- gBanFuns$get_pubchem_property_table(
    cids = cid,
    sleep = sleep
  )

  if (!nrow(data_property)) {
    data_property <- data.frame(
      CID = cid,
      stringsAsFactors = FALSE
    )
  }

  for (col in c(
    "MolecularFormula",
    "MolecularWeight",
    "CanonicalSMILES",
    "IsomericSMILES",
    "InChIKey",
    "IUPACName",
    "Title"
  )) {
    if (!col %in% colnames(data_property)) {
      data_property[[col]] <- NA
    }
  }

  data_property$CID <- gBanFuns$normalize_pubchem_cid(data_property$CID)

  title <- gBanFuns$read_pubchem_record_title(cid)

  if (!is.na(title) && title != "") {
    data_property$Title[1L] <- title
  }

  synonym <- if (!is.na(data_property$Title[1L]) && data_property$Title[1L] != "") {
    as.character(data_property$Title[1L])
  } else if (!is.na(data_property$IUPACName[1L]) && data_property$IUPACName[1L] != "") {
    as.character(data_property$IUPACName[1L])
  } else if (!is.na(compound) && compound != "") {
    as.character(compound)
  } else {
    NA_character_
  }

  data.frame(
    SMILES = smiles,
    Compound = compound,
    CID = cid,
    Synonym = synonym,
    PubChem_strategy = strategy,
    PubChem_n_cid = length(cids),
    PubChem_all_cids = paste(cids, collapse = ";"),
    MolecularFormula = data_property$MolecularFormula[1L],
    MolecularWeight_PubChem = data_property$MolecularWeight[1L],
    CanonicalSMILES_PubChem = data_property$CanonicalSMILES[1L],
    IsomericSMILES_PubChem = data_property$IsomericSMILES[1L],
    InChIKey = data_property$InChIKey[1L],
    IUPACName = data_property$IUPACName[1L],
    Title = data_property$Title[1L],
    stringsAsFactors = FALSE
  )
}

gBanFuns$annotate_pubchem_smiles <- function(smiles,
  compound = NULL,
  identity_types = c(
    "same_stereo_isotope",
    "same_stereo",
    "same_isotope",
    "same_connectivity"
  ),
  use_similarity = FALSE,
  similarity_threshold = 99L,
  sleep = .2
)
{
  smiles <- as.character(smiles)
  smiles <- smiles[!is.na(smiles) & smiles != ""]

  if (is.null(compound)) {
    compound <- paste0("Candidate_", seq_along(smiles))
  }

  compound <- as.character(compound)

  lst <- lapply(
    seq_along(smiles),
    function(i) {
      gBanFuns$annotate_one_pubchem_smiles(
        smiles = smiles[i],
        compound = compound[i],
        identity_types = identity_types,
        use_similarity = use_similarity,
        similarity_threshold = similarity_threshold,
        sleep = sleep
      )
    }
  )

  data <- do.call(rbind, lst)
  rownames(data) <- NULL
  data
}


gBanFuns$read_pubchem_synonyms <- function(cid)
{
  cid <- gBanFuns$normalize_pubchem_cid(cid)
  cid <- cid[!is.na(cid) & cid != ""]

  if (!length(cid)) {
    return(character(0))
  }

  url <- paste0(
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
    cid[1L],
    "/synonyms/JSON"
  )

  res <- gBanFuns$read_pubchem_json(url)

  if (is.null(res$InformationList$Information)) {
    return(character(0))
  }

  data_info <- res$InformationList$Information

  if (is.data.frame(data_info) && "Synonym" %in% colnames(data_info)) {
    syn <- data_info$Synonym[[1L]]
  } else if (is.list(data_info) && !is.null(data_info[[1L]]$Synonym)) {
    syn <- data_info[[1L]]$Synonym
  } else {
    syn <- character(0)
  }

  syn <- unique(as.character(syn))
  syn <- syn[!is.na(syn) & syn != ""]
  syn
}

gBanFuns$filter_general_synonyms <- function(data)
{
  if (!nrow(data)) {
    return(data)
  }

  data <- as.data.frame(data, stringsAsFactors = FALSE)

  if (!"syno" %in% colnames(data)) {
    stop('Column "syno" was not found.')
  }

  data <- data[!is.na(data$syno) & data$syno != "", , drop = FALSE]
  data <- data[!grepl("[0-9]{5}", data$syno), , drop = FALSE]
  data <- data[!grepl("^[A-Z-]{1,5}$", data$syno), , drop = FALSE]
  data <- data[!grepl("^[A-Z0-9]{1,}$", data$syno), , drop = FALSE]
  data <- data[!grepl("(?<=-)[A-Z0-9]{5,}$", data$syno, perl = TRUE), , drop = FALSE]
  data <- data[!grepl("^[0-9-]*$", data$syno), , drop = FALSE]
  data <- data[!grepl("^[A-Z]{14}-", data$syno), , drop = FALSE]
  data <- data[nchar(data$syno) <= 120L, , drop = FALSE]

  data
}

gBanFuns$pick_one_synonym <- function(syno,
  title = NA_character_,
  iupac = NA_character_,
  max_display_chars = 80L
)
{
  syno <- unique(as.character(syno))
  syno <- syno[!is.na(syno) & syno != ""]

  title <- as.character(title)[1L]
  iupac <- as.character(iupac)[1L]

  data <- data.frame(
    syno = syno,
    stringsAsFactors = FALSE
  )

  data <- gBanFuns$filter_general_synonyms(data)
  data_short <- data[nchar(data$syno) <= max_display_chars, , drop = FALSE]

  if (nrow(data_short) > 0L && exists("PickGeneral", mode = "function")) {
    res <- tryCatch(
      PickGeneral(data_short$syno),
      error = function(e) NA_character_
    )

    if (!is.na(res) && res != "" && nchar(res) <= max_display_chars) {
      return(as.character(res))
    }
  }

  if (!is.na(title) && title != "") {
    data_title <- data.frame(syno = title, stringsAsFactors = FALSE)
    data_title <- gBanFuns$filter_general_synonyms(data_title)

    if (nrow(data_title) > 0L && nchar(title) <= max_display_chars) {
      return(title)
    }
  }

  if (nrow(data_short) > 0L) {
    data_short$.score <- nchar(data_short$syno)
    data_short$.score <- data_short$.score + ifelse(grepl(";|\\||<|>|\\[|\\]", data_short$syno), 50L, 0L)
    data_short$.score <- data_short$.score + ifelse(grepl("^InChI|^CID|^SCHEMBL", data_short$syno, ignore.case = TRUE), 100L, 0L)
    data_short <- data_short[order(data_short$.score, nchar(data_short$syno), data_short$syno), , drop = FALSE]
    return(data_short$syno[1L])
  }

  if (!is.na(iupac) && iupac != "") {
    return(iupac)
  }

  NA_character_
}

gBanFuns$get_pubchem_selected_synonyms <- function(data_lookup,
  use_cache = TRUE,
  overwrite = FALSE,
  file_cache = NULL,
  sleep = .2,
  max_display_chars = 80L
)
{
  data_lookup <- as.data.frame(data_lookup, stringsAsFactors = FALSE)

  if (!"CID" %in% colnames(data_lookup)) {
    return(data.frame(
      CID = character(0),
      Synonym = character(0),
      stringsAsFactors = FALSE
    ))
  }

  if (!"Title" %in% colnames(data_lookup)) {
    data_lookup$Title <- NA_character_
  }

  if (!"IUPACName" %in% colnames(data_lookup)) {
    data_lookup$IUPACName <- NA_character_
  }

  data_lookup$CID <- gBanFuns$normalize_pubchem_cid(data_lookup$CID)
  data_lookup <- data_lookup[!is.na(data_lookup$CID) & data_lookup$CID != "", , drop = FALSE]

  if (!nrow(data_lookup)) {
    return(data.frame(
      CID = character(0),
      Synonym = character(0),
      stringsAsFactors = FALSE
    ))
  }

  cids <- unique(data_lookup$CID)
  cids <- cids[!is.na(cids) & cids != ""]

  data_cache <- data.frame(
    CID = character(0),
    syno = character(0),
    stringsAsFactors = FALSE
  )

  if (!is.null(file_cache) && isTRUE(use_cache) && !isTRUE(overwrite) &&
      file.exists(file_cache)) {
    data_cache <- readRDS(file_cache)
    data_cache <- as.data.frame(data_cache, stringsAsFactors = FALSE)

    if ("CID" %in% colnames(data_cache)) {
      data_cache$CID <- gBanFuns$normalize_pubchem_cid(data_cache$CID)
    }
  }

  cids_done <- unique(data_cache$CID)
  cids_query <- setdiff(cids, cids_done)

  if (length(cids_query) > 0L) {
    data_new <- lapply(
      cids_query,
      function(cid) {
        syn <- gBanFuns$read_pubchem_synonyms(cid)
        Sys.sleep(sleep)

        if (!length(syn)) {
          return(data.frame(
            CID = cid,
            syno = NA_character_,
            stringsAsFactors = FALSE
          ))
        }

        data.frame(
          CID = rep(cid, length(syn)),
          syno = syn,
          stringsAsFactors = FALSE
        )
      }
    )

    data_new <- Filter(Negate(is.null), data_new)

    if (length(data_new) > 0L) {
      data_new <- do.call(rbind, data_new)
      data_cache <- rbind(data_cache, data_new)
    }

    if (!is.null(file_cache)) {
      dir.create(dirname(file_cache), recursive = TRUE, showWarnings = FALSE)
      saveRDS(data_cache, file_cache)
    }
  }

  data_pick <- lapply(
    cids,
    function(cid) {
      data_one <- data_lookup[data_lookup$CID == cid, , drop = FALSE]
      syn <- data_cache$syno[data_cache$CID == cid]
      syn <- c(data_one$Title[1L], syn, data_one$IUPACName[1L], data_one$Synonym[1L])
      syn <- syn[!is.na(syn) & syn != ""]

      data.frame(
        CID = cid,
        Synonym = gBanFuns$pick_one_synonym(
          syno = syn,
          title = data_one$Title[1L],
          iupac = data_one$IUPACName[1L],
          max_display_chars = max_display_chars
        ),
        stringsAsFactors = FALSE
      )
    }
  )

  data_pick <- do.call(rbind, data_pick)
  rownames(data_pick) <- NULL
  data_pick
}

gBanFuns$make_final_candidate_transposed_table <- function(data,
  name_col = "Synonym", id_col = "CID", max_candidates_per_block = 6L)
{
  data <- as.data.frame(data, stringsAsFactors = FALSE)

  if (!nrow(data)) {
    return(tibble::as_tibble(data.frame(stringsAsFactors = FALSE)))
  }

  max_candidates_per_block <- as.integer(max_candidates_per_block[1L])
  if (is.na(max_candidates_per_block) || max_candidates_per_block < 1L) {
    max_candidates_per_block <- 6L
  }

  if (!name_col %in% colnames(data)) {
    name_col <- colnames(data)[1L]
  }

  data <- data.frame(lapply(data, function(v) {
    if (is.numeric(v)) {
      return(as.character(signif(v, 3L)))
    }
    as.character(v)
  }), stringsAsFactors = FALSE)

  preferred_cols <- c(
    name_col,
    id_col,
    "Safety_decision",
    "Hard_risk_n",
    "Warning_risk_n",
    "Major_risk_endpoints",
    "MW",
    "LogP",
    "TPSA",
    "Lipinski_pass"
  )
  preferred_cols <- unique(preferred_cols[preferred_cols %in% colnames(data)])

  if (!length(preferred_cols)) {
    preferred_cols <- setdiff(colnames(data), c("SMILES", "Canonical SMILES"))
  }

  candidate_id <- paste0("C", seq_len(nrow(data)))
  blocks <- split(seq_len(nrow(data)), ceiling(seq_len(nrow(data)) / max_candidates_per_block))

  data_blocks <- lapply(seq_along(blocks), function(i) {
    idx <- blocks[[i]]
    data_one <- data[idx, preferred_cols, drop = FALSE]
    data_one <- data.frame(lapply(data_one, function(v) {
      v <- as.character(v)
      v[is.na(v)] <- ""
      v
    }), stringsAsFactors = FALSE)

    rows <- c("Candidate_ID", preferred_cols)
    data_t <- data.frame(
      Group = paste0("Group_", i),
      Index = rows,
      stringsAsFactors = FALSE
    )

    for (j in seq_len(max_candidates_per_block)) {
      col_j <- paste0("Candidate_", j)
      if (j <= length(idx)) {
        value_j <- c(candidate_id[idx[j]], as.character(unlist(data_one[j, preferred_cols], use.names = FALSE)))
      } else {
        value_j <- rep("", length(rows))
      }
      data_t[[col_j]] <- value_j
    }

    data_t
  })

  data_out <- dplyr::bind_rows(data_blocks)
  tibble::as_tibble(data_out)
}

gBanFuns$make_candidate_column_notes <- function(data,
  name_col = "Synonym",
  id_col = "CID"
)
{
  data <- as.data.frame(data, stringsAsFactors = FALSE)

  if (!nrow(data)) {
    return("")
  }

  if (!name_col %in% colnames(data)) {
    return("")
  }

  names_cpd <- as.character(data[[name_col]])
  names_cpd[is.na(names_cpd) | names_cpd == ""] <- "Unknown"

  if (id_col %in% colnames(data)) {
    ids <- as.character(data[[id_col]])
    ids[is.na(ids) | ids == ""] <- "NA"
    notes <- glue::glue("C{seq_along(names_cpd)}: {names_cpd} (PubChem CID: {ids})")
  } else {
    notes <- glue::glue("C{seq_along(names_cpd)}: {names_cpd}")
  }

  paste(notes, collapse = "\n")
}


gBanFuns$get_custom_compound_metadata <- function(x)
{
  data_custom <- tryCatch(x$custom_input_compound_db, error = function(e) NULL)

  if (is.null(data_custom)) {
    data_custom <- tryCatch(x$custom_compound_db, error = function(e) NULL)
  }

  if (is.null(data_custom)) {
    data_custom <- tryCatch(x@params$custom_compound_db, error = function(e) NULL)
  }

  if (is.null(data_custom)) {
    return(NULL)
  }

  as.data.frame(data_custom, stringsAsFactors = FALSE)
}

gBanFuns$make_custom_compound_annotation <- function(x, smiles, compound = NULL)
{
  smiles <- as.character(smiles)
  smiles <- smiles[!is.na(smiles) & smiles != ""]

  if (is.null(compound)) {
    compound <- paste0("Candidate_", seq_along(smiles))
  }
  compound <- as.character(compound)
  if (length(compound) != length(smiles)) {
    compound <- rep(NA_character_, length(smiles))
  }

  data_base <- data.frame(
    SMILES = smiles,
    Compound = compound,
    stringsAsFactors = FALSE
  )

  data_custom <- gBanFuns$get_custom_compound_metadata(x)
  if (is.null(data_custom) || !nrow(data_custom)) {
    data_base$CID <- NA_character_
    data_base$Synonym <- data_base$Compound
    data_base$PubChem_strategy <- "no_custom_annotation"
    data_base$PubChem_n_cid <- 0L
    data_base$PubChem_all_cids <- NA_character_
    data_base$MolecularFormula <- NA_character_
    data_base$MolecularWeight_PubChem <- NA_real_
    data_base$CanonicalSMILES_PubChem <- data_base$SMILES
    data_base$IsomericSMILES_PubChem <- data_base$SMILES
    data_base$InChIKey <- NA_character_
    data_base$IUPACName <- NA_character_
    data_base$Title <- data_base$Compound
    data_base$Annotation_source <- NA_character_
    data_base$Custom_annotation <- FALSE
    data_base$PubChem_match <- FALSE
    return(data_base)
  }

  find_col <- function(data, candidates) {
    cols <- colnames(data)
    cols_key <- tolower(gsub("[ ._-]+", "", cols))
    cand_key <- tolower(gsub("[ ._-]+", "", candidates))
    idx <- match(cand_key, cols_key)
    idx <- idx[!is.na(idx)]
    if (!length(idx)) {
      return(NA_character_)
    }
    cols[idx[1L]]
  }

  pick_col <- function(data, col, default = NA_character_) {
    if (is.na(col) || !col %in% colnames(data)) {
      return(rep(default, nrow(data)))
    }
    as.character(data[[col]])
  }

  col_smiles <- find_col(data_custom, c("smiles", "SMILES", "canonical_smiles", "isomeric_smiles"))
  if (is.na(col_smiles)) {
    data_base$CID <- NA_character_
    data_base$Synonym <- data_base$Compound
    data_base$PubChem_strategy <- "custom_annotation_without_smiles"
    data_base$PubChem_n_cid <- 0L
    data_base$PubChem_all_cids <- NA_character_
    data_base$MolecularFormula <- NA_character_
    data_base$MolecularWeight_PubChem <- NA_real_
    data_base$CanonicalSMILES_PubChem <- data_base$SMILES
    data_base$IsomericSMILES_PubChem <- data_base$SMILES
    data_base$InChIKey <- NA_character_
    data_base$IUPACName <- NA_character_
    data_base$Title <- data_base$Compound
    data_base$Annotation_source <- NA_character_
    data_base$Custom_annotation <- FALSE
    data_base$PubChem_match <- FALSE
    return(data_base)
  }

  col_name <- find_col(data_custom, c("compound_name", "Compound name", "Compound", "name", "Synonym"))
  col_cid <- find_col(data_custom, c("pubchem_cid", "PubChem CID", "PubChem_CID", "CID", "cid"))
  col_inchikey <- find_col(data_custom, c("inchikey", "InChIKey", "InChI Key"))
  col_formula <- find_col(data_custom, c("molecular_formula", "MolecularFormula", "Formula"))
  col_weight <- find_col(data_custom, c("molecular_weight", "MolecularWeight", "MW"))
  col_herb <- find_col(data_custom, c("herb", "Herb", "query_herb", "herb_cn_name"))
  col_source <- find_col(data_custom, c("evidence_sources", "Evidence sources", "source_list", "source"))
  col_db_name <- find_col(data_custom, c("db_custom_name", "database", "source_database"))

  data_std <- data.frame(
    SMILES = trimws(as.character(data_custom[[col_smiles]])),
    Synonym_custom = pick_col(data_custom, col_name),
    CID_custom = pick_col(data_custom, col_cid),
    InChIKey_custom = pick_col(data_custom, col_inchikey),
    MolecularFormula_custom = pick_col(data_custom, col_formula),
    MolecularWeight_custom = suppressWarnings(as.numeric(pick_col(data_custom, col_weight))),
    herb = pick_col(data_custom, col_herb),
    evidence_sources = pick_col(data_custom, col_source),
    Annotation_source_custom = pick_col(data_custom, col_db_name, "Custom compound database"),
    stringsAsFactors = FALSE
  )
  data_std <- data_std[!is.na(data_std$SMILES) & data_std$SMILES != "", , drop = FALSE]

  if (!nrow(data_std)) {
    data_base$CID <- NA_character_
    data_base$Synonym <- data_base$Compound
    data_base$PubChem_strategy <- "no_custom_annotation"
    data_base$PubChem_n_cid <- 0L
    data_base$PubChem_all_cids <- NA_character_
    data_base$MolecularFormula <- NA_character_
    data_base$MolecularWeight_PubChem <- NA_real_
    data_base$CanonicalSMILES_PubChem <- data_base$SMILES
    data_base$IsomericSMILES_PubChem <- data_base$SMILES
    data_base$InChIKey <- NA_character_
    data_base$IUPACName <- NA_character_
    data_base$Title <- data_base$Compound
    data_base$Annotation_source <- NA_character_
    data_base$Custom_annotation <- FALSE
    data_base$PubChem_match <- FALSE
    return(data_base)
  }

  collapse_unique <- function(v) {
    v <- as.character(v)
    v <- unique(v[!is.na(v) & v != ""])
    if (!length(v)) {
      return(NA_character_)
    }
    paste(v, collapse = "; ")
  }

  data_split <- split(data_std, data_std$SMILES)
  data_sum <- do.call(rbind, lapply(names(data_split), function(smiles_i) {
    data_one <- data_split[[smiles_i]]
    data.frame(
      SMILES = smiles_i,
      Synonym_custom = collapse_unique(data_one$Synonym_custom),
      CID_custom = collapse_unique(gBanFuns$normalize_pubchem_cid(data_one$CID_custom)),
      InChIKey_custom = collapse_unique(data_one$InChIKey_custom),
      MolecularFormula_custom = collapse_unique(data_one$MolecularFormula_custom),
      MolecularWeight_custom = suppressWarnings(as.numeric(collapse_unique(data_one$MolecularWeight_custom))),
      herb = collapse_unique(data_one$herb),
      evidence_sources = collapse_unique(data_one$evidence_sources),
      Annotation_source_custom = collapse_unique(data_one$Annotation_source_custom),
      stringsAsFactors = FALSE
    )
  }))

  data_out <- merge(data_base, data_sum, by = "SMILES", all.x = TRUE)
  data_out$CID <- gBanFuns$normalize_pubchem_cid(data_out$CID_custom)
  data_out$Synonym <- ifelse(
    !is.na(data_out$Synonym_custom) & data_out$Synonym_custom != "",
    data_out$Synonym_custom,
    data_out$Compound
  )
  data_out$InChIKey <- data_out$InChIKey_custom
  data_out$MolecularFormula <- data_out$MolecularFormula_custom
  data_out$MolecularWeight_PubChem <- data_out$MolecularWeight_custom
  data_out$CanonicalSMILES_PubChem <- data_out$SMILES
  data_out$IsomericSMILES_PubChem <- data_out$SMILES
  data_out$IUPACName <- NA_character_
  data_out$Title <- data_out$Synonym
  data_out$PubChem_match <- !is.na(data_out$CID) & data_out$CID != ""
  data_out$PubChem_strategy <- ifelse(
    data_out$PubChem_match,
    "custom_database_cid",
    "custom_database_annotation"
  )
  data_out$PubChem_n_cid <- ifelse(data_out$PubChem_match, 1L, 0L)
  data_out$PubChem_all_cids <- data_out$CID
  data_out$Annotation_source <- data_out$Annotation_source_custom
  data_out$Custom_annotation <- !is.na(data_out$Synonym_custom) |
    (!is.na(data_out$CID) & data_out$CID != "")

  data_out$Synonym_custom <- NULL
  data_out$CID_custom <- NULL
  data_out$InChIKey_custom <- NULL
  data_out$MolecularFormula_custom <- NULL
  data_out$MolecularWeight_custom <- NULL
  data_out$Annotation_source_custom <- NULL

  data_out
}

gBanFuns$merge_final_candidate_tables <- function(data_pubchem,
  data_swiss = NULL,
  data_admet = NULL
)
{
  data_final <- data_pubchem

  if (!is.null(data_swiss) && "SMILES" %in% colnames(data_swiss)) {
    data_final <- merge(
      data_final,
      as.data.frame(data_swiss, stringsAsFactors = FALSE),
      by = "SMILES",
      all.x = TRUE
    )
  }

  if (!is.null(data_admet) && "SMILES" %in% colnames(data_admet)) {
    data_final <- merge(
      data_final,
      as.data.frame(data_admet, stringsAsFactors = FALSE),
      by = "SMILES",
      all.x = TRUE,
      suffixes = c("", "_ADMETlab")
    )
  }

  data_final
}

gBanFuns$make_final_candidate_report_table <- function(data,
  include_smiles = FALSE)
{
  data <- as.data.frame(data, stringsAsFactors = FALSE)

  if (!"PubChem_match" %in% colnames(data)) {
    data$PubChem_match <- !is.na(data$CID) & data$CID != ""
  }

  data$PubChem_status <- ifelse(
    data$PubChem_match,
    "Matched",
    "Not matched"
  )

  cols_core <- c(
    "SMILES",
    "Synonym",
    "CID",
    "Safety_decision",
    "Hard_risk_n",
    "Warning_risk_n",
    "Major_risk_endpoints"
  )

  cols_swiss <- c(
    "MW",
    "LogP",
    "TPSA",
    "Lipinski_pass",
    "Bioavailability_Score",
    "Synthetic_Accessibility"
  )

  cols_fallback <- c(
    "Risk_count",
    "Max_risk"
  )

  cols <- unique(c(cols_core, cols_swiss, cols_fallback))
  cols <- cols[cols %in% colnames(data)]

  if ("PubChem_status" %in% colnames(data)) {
    status_unique <- unique(as.character(data$PubChem_status[!is.na(data$PubChem_status)]))
    if (length(status_unique) > 1L) {
      cols <- unique(c(cols, "PubChem_status"))
    }
  }

  if (!isTRUE(include_smiles)) {
    cols <- setdiff(cols, "SMILES")
  }

  data_report <- data[, cols, drop = FALSE]

  if ("CID" %in% colnames(data_report)) {
    data_report$CID <- gBanFuns$normalize_pubchem_cid(data_report$CID)
  }

  data_report <- dplyr::mutate(
    data_report,
    dplyr::across(
      dplyr::where(is.numeric),
      function(x) signif(x, 3L)
    )
  )

  tibble::as_tibble(data_report)
}

get_compounds_weight <- function(smiles_list, HEAVY_METAL_THRESHOLD = 20, 
  molecules = NULL, libPaths)
{
  # smiles_list <- head(smiles_list, n = 1000)
  .libPaths(libPaths)
  if (is.null(molecules)) {
    molecules <- rcdk::parse.smiles(smiles_list)
  }
  results <- pbapply::pblapply(smiles_list,
    function(smile) {
      molecule <- molecules[[ smile ]]
      if (is.null(molecule)) {
        warning(paste("Unable to parse SMILES: ", smile))
        return()
      }
      rcdk::convert.implicit.to.explicit(molecule)
      formula <- rcdk::get.mol2formula(molecule, charge = 0)
      mol_weight <- formula@mass
      has_heavy_metal <- FALSE
      atoms <- rcdk::get.atoms(molecule)
      for (atom in atoms) {
        atomic_num <- rcdk::get.atomic.number(atom)
        if (atomic_num > HEAVY_METAL_THRESHOLD) {
          has_heavy_metal <- TRUE
          break
        }
      }
      data.frame(
        smiles = smile,
        MolecularWeight = mol_weight,
        HasHeavyMetal = has_heavy_metal,
        stringsAsFactors = FALSE
      )
    })
  do.call(dplyr::bind_rows, results)
  dplyr::as_tibble(data.table::rbindlist(results))
}

