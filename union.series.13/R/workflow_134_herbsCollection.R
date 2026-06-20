# ==========================================================================
# workflow of herbsCollection
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_herbsCollection <- setClass("job_herbsCollection",
  contains = c("job"),
  prototype = prototype(
    pg = "herbsCollection",
    info = c(""),
    cite = "",
    method = "",
    tag = "herbsCollection",
    analysis = ""
    ))

herbFuns <- new.env(parent = emptyenv())

if (!exists("metaFuns", inherits = TRUE)) {
  metaFuns <- new.env(parent = emptyenv())
}

job_herbsCollection <- function(herbs,
  sources = c("tcmsp", "batman"), data_literature = NULL, ...)
{
  collection <- herbFuns$collect_herbs(
    herbs = herbs,
    sources = sources,
    ...
  )

  if (!is.null(data_literature)) {
    collection <- herbFuns$inject_literature_table(
      collection,
      data_literature = data_literature
    )
  }

  x <- .job_herbsCollection()
  x$collection <- collection

  return(x)
}

setMethod("step0", signature = c(x = "job_herbsCollection"),
  function(x)
  {
    step_message("Prepare herbs collection with `job_herbsCollection()`.")
    return(x)
  })

setMethod("step1", signature = c(x = "job_herbsCollection"),
  function(x)
  {
    step_message("Stat compounds.")

    if (is.null(x$collection)) {
      stop("`x$collection` was not found. Please create the object with `job_herbsCollection()` first.")
    }

    collection <- x$collection
    herbFuns$validate_collection(collection, stop_if_invalid = TRUE)

    if (is.null(x$lst_refine)) {
      x$lst_refine <- list()
    }

    data_source <- herbFuns$stat_source_compounds(collection)
    data_herb <- herbFuns$stat_herb_compounds(collection)
    data_report <- herbFuns$stat_report_herb_summary(collection)
    data_missing <- herbFuns$diagnose_missing_structures(collection)

    x$lst_refine$stat_source_compounds <- data_source
    x$lst_refine$stat_herb_compounds <- data_herb
    x$lst_refine$report_herb_compound_summary <- data_report
    x$lst_refine$missing_structures <- data_missing

    t.report <- set_lab_legend(
      data_report,
      glue::glue("{x@sig} herb-compound evidence summary"),
      glue::glue(
        "药物–化合物证据汇总表|||该表基于 TCMSP、BATMAN-TCM 及文献来源整合各药物的候选化合物信息，并统计不同证据来源下的化合物数量及结构信息覆盖情况，包括 SMILES、InChIKey 和 PubChem CID。"
      )
    )
    x <- tablesAdd(x, t.report)

    t.source <- set_lab_legend(
      data_source,
      glue::glue("{x@sig} compound source summary"),
      glue::glue(
        "化合物来源统计表|||该表展示 TCMSP、BATMAN-TCM 及文献来源支持的药物-化合物关系数量、",
        "去重化合物数量，以及 SMILES、InChIKey 和 PubChem CID 等结构信息覆盖情况。"
      )
    )
    x <- tablesAdd(x, t.source)

    t.herb <- set_lab_legend(
      data_herb,
      glue::glue("{x@sig} herb-level compound summary"),
      glue::glue(
        "药物层面化合物统计表|||该表按药物及来源汇总候选化合物数量和结构信息覆盖情况，",
        "用于评估不同药物及不同证据来源对候选化合物集合的贡献。"
      )
    )
    x <- tablesAdd(x, t.herb)

    if (nrow(data_missing) > 0L) {
      t.missing <- set_lab_legend(
        data_missing,
        glue::glue("{x@sig} compounds requiring structure completion"),
        glue::glue(
          "待补全结构信息化合物表|||该表列出当前缺少 SMILES、InChIKey 或 PubChem CID 的候选化合物，",
          "用于后续基于 PubChem 或人工校正的结构信息补全。"
        )
      )
      x <- tablesAdd(x, t.missing)
    }

    n_source <- length(unique(data_source$source))
    n_relationship <- sum(data_source$n_relationship, na.rm = TRUE)
    n_compound_unique <- length(unique(collection$compound_unique$compound_key))
    n_with_smiles <- sum(collection$compound_unique$has_smiles, na.rm = TRUE)
    n_without_smiles <- sum(!collection$compound_unique$has_smiles, na.rm = TRUE)
    source_text <- herbFuns$format_source_label(unique(data_source$source))

    x <- methodAdd(
      x,
      "本分析整合 {source_text} 来源中的药物-化合物信息，建立候选化合物集合。对于每个药物，分别统计不同来源支持的化合物数量，并进一步评估化合物结构信息，包括 SMILES、InChIKey 和 PubChem CID 的覆盖情况。"
    )

    x <- snapAdd(
      x,
      "共整合 {n_source} 类来源、{n_relationship} 条药物-化合物关系，获得 {n_compound_unique} 个去重候选化合物；其中 {n_with_smiles} 个具有 SMILES 结构信息，{n_without_smiles} 个暂未具有 SMILES 结构信息。"
    )

    return(x)
  })

setMethod("step2", signature = c(x = "job_herbsCollection"),
  function(x, dir_cache = file.path("tmp", "pubchem_herbsCollection"),
    manual_cid = NULL, manual_id = NULL, dic_name = NULL,
    prefer_isomeric_smiles = TRUE, verbose = TRUE)
  {
    step_message("Complete compound structures with PubChem.")

    if (is.null(x$collection)) {
      stop("`x$collection` was not found. Please run `job_herbsCollection()` first.")
    }

    if (!is.null(manual_id)) {
      if (!is.null(manual_cid)) {
        stop("Please provide only one of `manual_cid` or `manual_id`.")
      }
      manual_cid <- manual_id
    }

    collection <- x$collection
    herbFuns$validate_collection(collection, stop_if_invalid = TRUE)

    if (!dir.exists(dir_cache)) {
      dir.create(dir_cache, recursive = TRUE, showWarnings = FALSE)
    }
    if (!dir.exists(dir_cache)) {
      stop("Cannot create PubChem cache directory: ", dir_cache)
    }

    if (is.null(x$lst_refine)) {
      x$lst_refine <- list()
    }

    data_before <- herbFuns$stat_structure_coverage(collection, label = "Before PubChem completion")
    n_manual_input <- herbFuns$count_manual_cid_input(manual_cid)

    res_pubchem <- herbFuns$complete_compound_structures_pubchem(
      collection = collection,
      dir_cache = dir_cache,
      manual_cid = manual_cid,
      dic_name = dic_name,
      prefer_isomeric_smiles = prefer_isomeric_smiles,
      verbose = verbose
    )

    x$collection <- res_pubchem$collection

    data_after <- herbFuns$stat_structure_coverage(
      x$collection,
      label = "After PubChem completion"
    )
    data_coverage <- dplyr::bind_rows(data_before, data_after)
    data_report <- herbFuns$stat_report_pubchem_final_summary(x$collection)
    data_catalog <- herbFuns$make_final_compound_catalog(x$collection)

    x$lst_refine$pubchem_completion <- res_pubchem
    x$lst_refine$pubchem_structure_coverage <- data_coverage
    x$lst_refine$pubchem_completion_summary <- res_pubchem$summary
    x$lst_refine$report_pubchem_final_summary <- data_report
    x$lst_refine$final_compound_catalog <- data_catalog
    x$lst_refine$pubchem_cid_mapping <- res_pubchem$cid_mapping
    x$lst_refine$pubchem_cid_mapping_inchikey <- res_pubchem$cid_mapping_inchikey
    x$lst_refine$pubchem_cid_mapping_smiles <- res_pubchem$cid_mapping_smiles
    x$lst_refine$pubchem_property_mapping <- res_pubchem$property_mapping
    x$lst_refine$pubchem_no_smiles_compounds <- res_pubchem$unresolved

    t.report <- set_lab_legend(
      data_report,
      glue::glue("{x@sig} final herb-compound structure summary"),
      glue::glue(
        "最终药物–化合物结构信息汇总表|||该表展示 PubChem 结构信息映射后各药物候选化合物的最终数量、证据来源以及 SMILES、InChIKey 和 PubChem CID 覆盖情况。"
      )
    )
    x <- tablesAdd(x, t.report)

    t.catalog <- set_lab_legend(
      data_catalog,
      glue::glue("{x@sig} final compound catalog"),
      glue::glue(
        "候选化合物总表|||该表列出各药物对应的候选化合物、证据来源、PubChem CID、SMILES、InChIKey、InChI、分子式和分子量，用于展示本分析纳入的药物–化合物数据基础。"
      )
    )
    x <- tablesAdd(x, t.catalog)

    t.coverage <- set_lab_legend(
      data_coverage,
      glue::glue("{x@sig} PubChem structure mapping coverage"),
      glue::glue(
        "PubChem 结构信息映射覆盖度表|||该表展示 PubChem 结构信息映射前后候选化合物的 SMILES、InChIKey、PubChem CID 及完整结构信息覆盖情况。"
      )
    )
    x <- tablesAdd(x, t.coverage)

    t.summary <- set_lab_legend(
      res_pubchem$summary,
      glue::glue("{x@sig} PubChem structure mapping summary"),
      glue::glue(
        "PubChem 结构信息映射统计表|||该表汇总基于 PubChem CID、化合物名称、InChIKey 和 SMILES 映射候选化合物结构信息的结果，包括新增 PubChem CID、SMILES 和 InChIKey 的数量，以及最终未获得 SMILES 的化合物数量。"
      )
    )
    x <- tablesAdd(x, t.summary)

    n_smiles_before <- data_before$n_with_smiles[[1L]]
    n_smiles_after <- data_after$n_with_smiles[[1L]]
    n_inchikey_before <- data_before$n_with_inchikey[[1L]]
    n_inchikey_after <- data_after$n_with_inchikey[[1L]]
    n_cid_before <- data_before$n_with_pubchem_cid[[1L]]
    n_cid_after <- data_after$n_with_pubchem_cid[[1L]]
    n_no_smiles <- nrow(res_pubchem$unresolved)

    manual_text <- ""
    if (n_manual_input > 0L) {
      manual_text <- paste0(
        "同时，对于经文献信息或 PubChem 页面检索可明确对应的化合物，纳入 ",
        n_manual_input,
        " 条手工 PubChem CID 映射作为补充证据。"
      )
    }

    x <- methodAdd(
      x,
      paste0(
        "为统一候选化合物的结构信息，本分析使用 R 包 `PubChemR` ⟦pkgInfo('PubChemR')⟧基于 PubChem Compound 数据库进行结构信息映射。",
        "对于已具有 PubChem CID 的化合物，按 CID 获取 SMILES、InChIKey、InChI、分子式和分子量；对于缺少 CID 的条目，先对化合物名称中由 PDF 换行、连字符、逗号和多余空格造成的检索噪声进行规范化，并以原始名称、规范化名称及括号内别名作为候选名称映射 PubChem CID，再依次根据 InChIKey 和已有 SMILES 映射 PubChem CID，最终获取对应结构信息。",
        manual_text,
        "对于经上述映射后仍未获得明确 SMILES 的化合物，保留其名称和来源证据；后续涉及结构依赖分析时，以具备 SMILES 结构式的候选化合物为主要分析对象。"
      )
    )

    x <- snapAdd(
      x,
      "经 PubChem 结构信息映射后，具有 SMILES 的候选化合物由 {n_smiles_before} 个增加至 {n_smiles_after} 个，具有 InChIKey 的候选化合物由 {n_inchikey_before} 个增加至 {n_inchikey_after} 个，具有 PubChem CID 的候选化合物由 {n_cid_before} 个增加至 {n_cid_after} 个；最终 {n_no_smiles} 个候选化合物未获得明确 SMILES 结构式，已作为名称来源证据保留。"
    )

    return(x)
  })


setMethod("step3", signature = c(x = "job_herbsCollection"),
  function(x, data_tcmsp_targets = NULL, data_batman_targets = NULL,
    read_local = TRUE, verbose = TRUE)
  {
    step_message("Collect compound targets.")

    if (is.null(x$collection)) {
      stop("`x$collection` was not found. Please run `job_herbsCollection()` first.")
    }

    if (is.null(x$lst_refine) ||
        is.null(x$lst_refine$pubchem_completion) ||
        is.null(x$lst_refine$final_compound_catalog)) {
      stop("Please run `step2()` before `step3()` so compound structures and PubChem CIDs are finalized.")
    }

    collection <- x$collection
    herbFuns$validate_collection(collection, stop_if_invalid = TRUE)

    data_index <- herbFuns$prepare_current_compound_target_index(collection)
    if (nrow(data_index) == 0L) {
      warning("No candidate compound was available for target annotation.")
      return(x)
    }

    if (is.null(data_tcmsp_targets)) {
      if (!isTRUE(read_local)) {
        data_tcmsp_targets <- NULL
      } else {
        data_tcmsp_targets <- herbFuns$read_tcmsp_targets(verbose = verbose)
      }
    }

    if (is.null(data_batman_targets)) {
      if (!isTRUE(read_local)) {
        data_batman_targets <- NULL
      } else {
        data_batman_targets <- herbFuns$read_batman_targets(verbose = verbose)
      }
    }

    data_tcmsp <- herbFuns$collect_tcmsp_targets(
      data_index = data_index,
      data_tcmsp_targets = data_tcmsp_targets,
      verbose = verbose
    )

    data_batman <- herbFuns$collect_batman_targets(
      data_index = data_index,
      data_batman_targets = data_batman_targets,
      verbose = verbose
    )

    data_full <- dplyr::bind_rows(data_tcmsp, data_batman)
    data_full <- dplyr::distinct(data_full)

    if (nrow(data_full) == 0L) {
      warning("No compound-target annotation was obtained from TCMSP or BATMAN-TCM for the current candidate compounds.")
      return(x)
    }

    data_full <- dplyr::arrange(
      data_full,
      query_herb,
      compound_name,
      target_source,
      target_evidence_type,
      target_gene
    )

    data_catalog <- herbFuns$make_compound_target_catalog(data_full)
    data_report <- herbFuns$stat_report_compound_target_summary(data_full)
    data_source <- herbFuns$stat_target_source_summary(data_full)

    x$lst_refine$compound_target_catalog_full <- data_full
    x$lst_refine$compound_target_catalog <- data_catalog
    x$lst_refine$report_compound_target_summary <- data_report
    x$lst_refine$target_source_summary <- data_source

    t.report <- set_lab_legend(
      data_report,
      glue::glue("{x@sig} compound-target summary"),
      glue::glue(
        "候选成分–靶点汇总表|||该表按药物汇总候选成分对应的靶点注释结果，",
        "包括具有靶点注释的化合物数量、成分–靶点关系数量、唯一靶点基因数量及不同靶点证据来源的贡献。"
      )
    )
    x <- tablesAdd(x, t.report)

    t.catalog <- set_lab_legend(
      data_catalog,
      glue::glue("{x@sig} compound-target catalog"),
      glue::glue(
        "候选成分–靶点总表|||该表列出各药物候选化合物对应的靶点基因、靶点名称、PubChem CID、",
        "化合物证据来源、靶点来源及证据类型，用于展示本分析纳入的成分–靶点数据基础。"
      )
    )
    x <- tablesAdd(x, t.catalog)

    t.source <- set_lab_legend(
      data_source,
      glue::glue("{x@sig} target evidence source summary"),
      glue::glue(
        "靶点证据来源统计表|||该表按 TCMSP、BATMAN-TCM 已知靶点和 BATMAN-TCM 预测靶点汇总候选成分–靶点关系数量、",
        "具有靶点注释的化合物数量及唯一靶点基因数量。"
      )
    )
    x <- tablesAdd(x, t.source)

    n_target_links <- nrow(data_full)
    n_compounds_with_targets <- dplyr::n_distinct(data_full$compound_key)
    vec_target_gene <- herbFuns$clean_text(data_full$target_gene)
    n_target_genes <- length(unique(vec_target_gene[!is.na(vec_target_gene) & nzchar(vec_target_gene)]))
    n_tcmsp_links <- sum(data_full$target_source == "TCMSP", na.rm = TRUE)
    n_batman_known_links <- sum(
      data_full$target_source == "BATMAN-TCM" &
        data_full$target_evidence_type == "known target",
      na.rm = TRUE
    )
    n_batman_predicted_links <- sum(
      data_full$target_source == "BATMAN-TCM" &
        data_full$target_evidence_type == "predicted target",
      na.rm = TRUE
    )

    x <- methodAdd(
      x,
      "为补充候选化合物的潜在靶点信息，本分析进一步整合 TCMSP 与 BATMAN-TCM 的成分–靶点注释。TCMSP 靶点基于候选化合物对应的 MOL_ID 进行匹配；BATMAN-TCM 靶点基于 PubChem CID 匹配已知靶点与预测靶点。对于 BATMAN-TCM 预测靶点，原始 Entrez Gene ID 通过 R 包 `AnnotationDbi` ⟦pkgInfo('AnnotationDbi')⟧与 `org.Hs.eg.db` ⟦pkgInfo('org.Hs.eg.db')⟧转换为 HGNC gene symbol。最终保留药物、候选化合物、靶点基因、靶点来源及证据类型，用于形成候选成分–靶点信息表。"
    )

    x <- snapAdd(
      x,
      "基于 TCMSP 与 BATMAN-TCM 成分–靶点注释，共获得 {n_target_links} 条候选成分–靶点关系，涉及 {n_compounds_with_targets} 个候选化合物和 {n_target_genes} 个靶点基因；其中 TCMSP 支持 {n_tcmsp_links} 条关系，BATMAN-TCM 已知靶点支持 {n_batman_known_links} 条关系，BATMAN-TCM 预测靶点支持 {n_batman_predicted_links} 条关系。"
    )

    return(x)
  })

# --------------------------------------------------------------------------
# General utilities
# --------------------------------------------------------------------------

herbFuns$clean_text <- function(x)
{
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "Na", "na", "NULL", "Null", "null", "N/A", "n/a", "-")] <- NA_character_
  x
}

herbFuns$get_col <- function(data, col, default = NA_character_)
{
  if (is.null(col) || is.na(col) || !col %in% colnames(data)) {
    return(rep(default, nrow(data)))
  }

  data[[col]]
}

herbFuns$resolve_col <- function(data, candidates, required = TRUE)
{
  hit <- candidates[candidates %in% colnames(data)]

  if (length(hit) > 0L) {
    return(hit[[1L]])
  }

  if (required) {
    stop("Cannot find required column. Candidates: ", paste(candidates, collapse = ", "))
  }

  NA_character_
}

herbFuns$match_key <- function(x)
{
  x <- herbFuns$clean_text(x)
  x <- tolower(x)
  x <- gsub("\\s+", "", x)
  x
}

herbFuns$collapse_unique <- function(x, sep = ";")
{
  x <- herbFuns$clean_text(x)
  x <- unique(x[!is.na(x) & nzchar(x)])

  if (length(x) == 0L) {
    return(NA_character_)
  }

  paste(x, collapse = sep)
}

herbFuns$first_non_missing <- function(x)
{
  x <- herbFuns$clean_text(x)
  x <- x[!is.na(x) & nzchar(x)]

  if (length(x) == 0L) {
    return(NA_character_)
  }

  x[[1L]]
}

herbFuns$format_n <- function(x)
{
  format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

herbFuns$log_progress <- function(verbose = TRUE, ...)
{
  if (isTRUE(verbose)) {
    message("[herbFuns] ", paste0(..., collapse = ""))
  }
  invisible(NULL)
}

herbFuns$format_source_label <- function(source)
{
  source <- unique(as.character(source))
  source <- source[!is.na(source) & nzchar(source)]

  dic <- c(
    tcmsp = "TCMSP",
    batman = "BATMAN-TCM",
    literature = "Literature"
  )

  out <- ifelse(source %in% names(dic), unname(dic[source]), source)
  paste(out, collapse = ", ")
}


herbFuns$ensure_columns <- function(data, cols, value = NA_character_)
{
  data <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)

  for (col in cols) {
    if (!col %in% colnames(data)) {
      data[[col]] <- rep(value, nrow(data))
    }
  }

  data
}

herbFuns$ensure_compound_columns <- function(data)
{
  herbFuns$ensure_columns(
    data,
    c(
      "source", "compound_source_id", "compound_name", "mol_id",
      "pubchem_cid", "smiles", "inchi", "inchikey",
      "molecular_formula", "molecular_weight", "source_priority",
      "compound_key"
    )
  )
}

herbFuns$ensure_relationship_columns <- function(data)
{
  data <- herbFuns$ensure_columns(
    data,
    c(
      "query_id", "query_herb", "source", "source_herb_id",
      "herb_cn_name", "herb_en_name", "herb_pinyin", "herb_latin_name",
      "compound_source_id", "compound_name", "compound_name_source",
      "relationship_source", "evidence_source", "compound_key"
    )
  )

  # Structure fields are resolved from `compound`; keeping them in relationship
  # tables may create `.x/.y` suffixes after joins.
  dplyr::select(
    data,
    -dplyr::any_of(c(
      "pubchem_cid", "smiles", "inchi", "inchikey",
      "molecular_formula", "molecular_weight", "has_smiles",
      "has_inchikey", "has_pubchem_cid", "structure_status"
    ))
  )
}

herbFuns$prepare_query <- function(herbs)
{
  herbs <- herbFuns$clean_text(herbs)
  herbs <- herbs[!is.na(herbs) & nzchar(herbs)]
  herbs <- unique(herbs)

  tibble::tibble(
    query_id = seq_along(herbs),
    query_herb = herbs,
    query_key = herbFuns$match_key(herbs)
  )
}

herbFuns$make_source_herb_id <- function(source, ...)
{
  value <- unlist(list(...), use.names = FALSE)
  value <- herbFuns$clean_text(value)
  value <- value[!is.na(value) & nzchar(value)]

  if (length(value) == 0L) {
    value <- "unknown"
  }

  paste0(source, ":", paste(value, collapse = "|"))
}

herbFuns$make_compound_key <- function(inchikey = NA_character_, smiles = NA_character_,
  pubchem_cid = NA_character_, source = NA_character_, compound_source_id = NA_character_,
  compound_name = NA_character_)
{
  inchikey <- herbFuns$clean_text(inchikey)
  smiles <- herbFuns$clean_text(smiles)
  pubchem_cid <- herbFuns$clean_text(pubchem_cid)
  source <- herbFuns$clean_text(source)
  compound_source_id <- herbFuns$clean_text(compound_source_id)
  compound_name <- herbFuns$clean_text(compound_name)

  out <- rep(NA_character_, max(
    length(inchikey), length(smiles), length(pubchem_cid),
    length(source), length(compound_source_id), length(compound_name)
  ))

  inchikey <- rep(inchikey, length.out = length(out))
  smiles <- rep(smiles, length.out = length(out))
  pubchem_cid <- rep(pubchem_cid, length.out = length(out))
  source <- rep(source, length.out = length(out))
  compound_source_id <- rep(compound_source_id, length.out = length(out))
  compound_name <- rep(compound_name, length.out = length(out))

  idx <- !is.na(inchikey) & nzchar(inchikey)
  out[idx] <- paste0("inchikey:", inchikey[idx])

  idx <- is.na(out) & !is.na(smiles) & nzchar(smiles)
  out[idx] <- paste0("smiles:", smiles[idx])

  idx <- is.na(out) & !is.na(pubchem_cid) & nzchar(pubchem_cid)
  out[idx] <- paste0("cid:", pubchem_cid[idx])

  idx <- is.na(out) & !is.na(source) & nzchar(source) &
    !is.na(compound_source_id) & nzchar(compound_source_id)
  out[idx] <- paste0(source[idx], ":", compound_source_id[idx])

  idx <- is.na(out) & !is.na(compound_name) & nzchar(compound_name)
  out[idx] <- paste0("name:", tolower(trimws(compound_name[idx])))

  idx <- is.na(out)
  out[idx] <- paste0("unknown:", seq_len(sum(idx)))

  out
}

herbFuns$add_structure_status <- function(data)
{
  data <- herbFuns$ensure_compound_columns(data)
  data$molecular_weight <- suppressWarnings(as.numeric(data$molecular_weight))
  data$has_smiles <- !is.na(herbFuns$clean_text(data$smiles)) & nzchar(herbFuns$clean_text(data$smiles))
  data$has_inchikey <- !is.na(herbFuns$clean_text(data$inchikey)) & nzchar(herbFuns$clean_text(data$inchikey))
  data$has_pubchem_cid <- !is.na(herbFuns$clean_text(data$pubchem_cid)) & nzchar(herbFuns$clean_text(data$pubchem_cid))
  data$has_compound_name <- !is.na(herbFuns$clean_text(data$compound_name)) & nzchar(herbFuns$clean_text(data$compound_name))

  data$structure_status <- ifelse(
    data$has_smiles & data$has_inchikey,
    "complete_structure",
    ifelse(
      data$has_smiles | data$has_inchikey,
      "partial_structure",
      ifelse(
        data$has_pubchem_cid,
        "cid_only",
        ifelse(data$has_compound_name, "name_only", "missing_structure")
      )
    )
  )

  data
}

herbFuns$new_empty_source_result <- function(query = NULL)
{
  if (is.null(query)) {
    query <- tibble::tibble(
      query_id = integer(),
      query_herb = character(),
      query_key = character()
    )
  }

  list(
    query = query,
    herb_map = tibble::tibble(),
    herb_compound = tibble::tibble(),
    compound = tibble::tibble(),
    logs = list()
  )
}

# --------------------------------------------------------------------------
# PubChem structure completion helpers
# --------------------------------------------------------------------------

metaFuns$get_properties_pubchem_batch <- function(identifier,
  properties = c(
    "MolecularFormula", "MolecularWeight", "CanonicalSMILES",
    "IsomericSMILES", "InChI", "InChIKey"
  ), namespace = "cid", searchtype = NULL, options = NULL,
  propertyMatch = list(.ignore.case = FALSE, type = "match"))
{
  PubChemR::get_properties(
    properties = properties,
    identifier = identifier,
    namespace = namespace,
    searchtype = searchtype,
    options = options,
    propertyMatch = propertyMatch
  )
}

metaFuns$extract_pubchem_properties <- function(object_property,
  identifier = NULL, prefer_isomeric_smiles = TRUE)
{
  data_prop <- tryCatch(
    as.data.frame(
      PubChemR::retrieve(
        object_property,
        .combine.all = TRUE,
        .to.data.frame = TRUE
      ),
      stringsAsFactors = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(data_prop) || nrow(data_prop) == 0L) {
    data_prop <- tryCatch(
      as.data.frame(
        PubChemR::retrieve(object_property, .combine.all = TRUE),
        stringsAsFactors = FALSE
      ),
      error = function(e) NULL
    )
  }

  if (is.null(data_prop) || nrow(data_prop) == 0L) {
    data_prop <- tryCatch(
      as.data.frame(PubChemR::retrieve(object_property), stringsAsFactors = FALSE),
      error = function(e) NULL
    )
  }

  if (is.null(data_prop) || nrow(data_prop) == 0L) {
    data_prop <- tryCatch(
      as.data.frame(object_property, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
  }

  if (is.null(data_prop) || nrow(data_prop) == 0L) {
    return(tibble::tibble(
      pubchem_cid = character(0L),
      canonical_smiles = character(0L),
      isomeric_smiles = character(0L),
      smiles = character(0L),
      inchi = character(0L),
      inchikey = character(0L),
      molecular_formula = character(0L),
      molecular_weight = numeric(0L)
    ))
  }

  if (!is.null(rownames(data_prop))) {
    row_id <- rownames(data_prop)
    if (!is.null(row_id) && !all(row_id %in% as.character(seq_len(nrow(data_prop))))) {
      data_prop$pubchem_cid_from_rownames <- row_id
    }
  }

  col_lower <- tolower(colnames(data_prop))
  get_col_by_pattern <- function(pattern, default = NA_character_) {
    idx <- grep(pattern, col_lower)
    if (length(idx) == 0L) {
      return(rep(default, nrow(data_prop)))
    }
    data_prop[[idx[[1L]]]]
  }

  cid <- get_col_by_pattern(
    "^cid$|compoundcid|pubchem.*cid|identifier|input|query|pubchem_cid_from_rownames"
  )

  if (all(is.na(cid)) && !is.null(identifier) && length(identifier) == nrow(data_prop)) {
    cid <- identifier
  }

  canonical_smiles <- get_col_by_pattern("canonical.*smiles|canonicalsmiles")
  isomeric_smiles <- get_col_by_pattern("isomeric.*smiles|isomericsmiles")
  smiles <- ifelse(
    isTRUE(prefer_isomeric_smiles) & !is.na(isomeric_smiles) & nzchar(as.character(isomeric_smiles)),
    isomeric_smiles,
    canonical_smiles
  )
  smiles <- ifelse(
    is.na(smiles) | !nzchar(as.character(smiles)),
    get_col_by_pattern("^smiles$"),
    smiles
  )

  out <- tibble::tibble(
    pubchem_cid = as.character(cid),
    canonical_smiles = herbFuns$clean_text(canonical_smiles),
    isomeric_smiles = herbFuns$clean_text(isomeric_smiles),
    smiles = herbFuns$clean_text(smiles),
    inchi = herbFuns$clean_text(get_col_by_pattern("^inchi$")),
    inchikey = herbFuns$clean_text(get_col_by_pattern("inchikey")),
    molecular_formula = herbFuns$clean_text(get_col_by_pattern("molecular.*formula|molecularformula")),
    molecular_weight = suppressWarnings(as.numeric(get_col_by_pattern("molecular.*weight|molecularweight")))
  )

  out$pubchem_cid <- gsub("^CID:", "", out$pubchem_cid)
  out$pubchem_cid <- trimws(out$pubchem_cid)
  out <- dplyr::filter(out, !is.na(pubchem_cid), nzchar(pubchem_cid))
  dplyr::distinct(out)
}

metaFuns$map_pubchem_cids_by_identifier <- function(identifier,
  namespace = "name", dir_cache = file.path("tmp", "pubchem_herbsCollection"),
  cache_name = NULL, domain = "compound", searchtype = NULL, options = NULL)
{
  if (exists("expect_package", mode = "function")) {
    expect_package("PubChemR", "3.0.0")
  } else {
    if (!requireNamespace("PubChemR", quietly = TRUE)) {
      stop("Package `PubChemR` is required for PubChem CID mapping.")
    }
  }

  identifier <- unique(trimws(as.character(identifier)))
  identifier <- identifier[!is.na(identifier) & nzchar(identifier)]

  if (length(identifier) == 0L) {
    return(tibble::tibble(
      identifier_original = character(0L),
      identifier_query = character(0L),
      namespace = character(0L),
      cid = character(0L),
      cid_n = integer(0L)
    ))
  }

  if (is.null(cache_name)) {
    cache_name <- glue::glue(
      "pubchemr_cids_{namespace}_{domain}_{length(identifier)}"
    )
  }

  if (!dir.exists(dir_cache)) {
    dir.create(dir_cache, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(dir_cache)) {
    stop("Cannot create PubChem cache directory: ", dir_cache)
  }

  object_cid <- tryCatch({
    if (exists("expect_local_data", mode = "function")) {
      expect_local_data(
        dir_cache,
        cache_name,
        metaFuns$get_cids_pubchem_batch,
        list(
          identifier = identifier,
          namespace = namespace,
          domain = domain,
          searchtype = searchtype,
          options = options
        )
      )
    } else {
      metaFuns$get_cids_pubchem_batch(
        identifier = identifier,
        namespace = namespace,
        domain = domain,
        searchtype = searchtype,
        options = options
      )
    }
  }, error = function(e) {
    warning(
      "PubChem CID mapping failed for namespace `", namespace, "`: ",
      conditionMessage(e)
    )
    NULL
  })

  if (is.null(object_cid)) {
    return(tibble::tibble(
      identifier_original = identifier,
      identifier_query = identifier,
      namespace = namespace,
      cid = NA_character_,
      cid_n = 0L
    ))
  }

  data_cid <- tryCatch(
    metaFuns$extract_pubchem_cids(object_cid, identifier = identifier),
    error = function(e) {
      warning(
        "Cannot extract PubChem CID mapping for namespace `", namespace, "`: ",
        conditionMessage(e)
      )
      NULL
    }
  )

  if (is.null(data_cid) || nrow(data_cid) == 0L) {
    data_out <- tibble::tibble(
      identifier_original = identifier,
      identifier_query = identifier,
      namespace = namespace,
      cid = NA_character_
    )
  } else {
    if (!"name_query" %in% colnames(data_cid)) {
      data_cid$name_query <- identifier[seq_len(min(length(identifier), nrow(data_cid)))]
    }

    data_out <- tibble::tibble(
      identifier_query = as.character(data_cid$name_query),
      namespace = namespace,
      cid = as.character(data_cid$cid)
    )

    data_query <- tibble::tibble(
      identifier_original = identifier,
      identifier_query = identifier
    )

    data_out <- dplyr::left_join(
      data_query,
      data_out,
      by = "identifier_query"
    )
    data_out$namespace <- ifelse(is.na(data_out$namespace), namespace, data_out$namespace)
  }

  data_out$cid <- as.character(data_out$cid)
  data_out$cid <- gsub("^CID:", "", data_out$cid)
  data_out$cid <- trimws(data_out$cid)

  data_count <- dplyr::group_by(data_out, identifier_original)
  data_count <- dplyr::summarise(
    data_count,
    cid_n = sum(!is.na(cid) & nzchar(cid)),
    .groups = "drop"
  )

  data_out$cid_n <- NULL
  data_out <- dplyr::left_join(data_out, data_count, by = "identifier_original")
  dplyr::distinct(data_out)
}


metaFuns$clean_pubchem_compound_name <- function(x)
{
  x <- herbFuns$clean_text(x)
  x <- gsub("[\r\n\t]+", " ", x)
  x <- gsub("[\u00A0\u2007\u202F]", " ", x, perl = TRUE)
  x <- gsub("[‐‑‒–—−－]", "-", x)
  x <- gsub("\\s+", " ", x)

  # Remove line-break artifacts around chemical punctuation.
  x <- gsub("\\s*-\\s*", "-", x)
  x <- gsub("\\s*,\\s*", ",", x)
  x <- gsub("\\s*;\\s*", ";", x)
  x <- gsub("\\s*:\\s*", ":", x)
  x <- gsub("\\(\\s+", "(", x)
  x <- gsub("\\s+\\)", ")", x)
  x <- gsub("\\[\\s+", "[", x)
  x <- gsub("\\s+\\]", "]", x)

  # Remove formula fragments occasionally inserted into names during PDF extraction.
  x <- gsub("\\bC[0-9]+H[0-9]+(?:[A-Z][a-z]?[0-9]*)*\\b", "", x, perl = TRUE)

  # Remove category labels accidentally appended to compound names.
  x <- gsub("\\bC21 Steroids\\b", "", x, ignore.case = TRUE)
  x <- gsub("\\bVolatile oils\\b", "", x, ignore.case = TRUE)
  x <- gsub("\\bPhenanthroindolizidine alkaloids\\b", "", x, ignore.case = TRUE)
  x <- gsub("\\bCarbohydrates\\b", "", x, ignore.case = TRUE)
  x <- gsub("\\bOther classes\\b", "", x, ignore.case = TRUE)

  # Remove table-column residues that may be attached after the actual name.
  x <- gsub("\\s+roots?\\s+[0-9/% ]*(water|methanol|ethanol|ethyl[- ]acetate).*extracts?$", "", x, ignore.case = TRUE)
  x <- gsub("\\s+stems?\\s+[0-9/% ]*(water|methanol|ethanol|ethyl[- ]acetate).*extracts?$", "", x, ignore.case = TRUE)
  x <- gsub("\\s+rhizomes?\\s+[0-9/% ]*(water|methanol|ethanol|ethyl[- ]acetate).*extracts?$", "", x, ignore.case = TRUE)
  x <- gsub("\\s+leaves?\\s+oils$", "", x, ignore.case = TRUE)
  x <- gsub("\\s+all herb\\s+oils$", "", x, ignore.case = TRUE)
  x <- gsub("rootsethanol extract$", "", x, ignore.case = TRUE)

  # Correct a small number of common PDF/OCR artifacts before PubChem search.
  x <- gsub("Didodecvl", "Didodecyl", x, ignore.case = FALSE)
  x <- gsub("Octadeceneamid\\b", "Octadecenamide", x, ignore.case = TRUE)
  x <- gsub("myristcin", "myristicin", x, ignore.case = TRUE)
  x <- gsub("Penten1-ol", "Penten-1-ol", x, ignore.case = TRUE)
  x <- gsub("5a-lanost", "5α-lanost", x, fixed = TRUE)
  x <- gsub("11a-hydroxy", "11α-hydroxy", x, fixed = TRUE)

  x <- gsub("\\s+", " ", x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NULL", "-", ".")] <- NA_character_
  x
}

metaFuns$make_pubchem_name_candidates <- function(vec_name)
{
  vec_name <- herbFuns$clean_text(vec_name)
  vec_name <- vec_name[!is.na(vec_name) & nzchar(vec_name)]
  vec_name <- unique(vec_name)

  if (length(vec_name) == 0L) {
    return(tibble::tibble(
      name_original = character(0L),
      name_query = character(0L),
      candidate_rank = integer(0L)
    ))
  }

  lst_candidate <- lapply(vec_name, function(name_i) {
    clean_i <- metaFuns$clean_pubchem_compound_name(name_i)
    out_i <- c(name_i, clean_i)

    if (!is.na(clean_i) && nzchar(clean_i)) {
      inside_i <- regmatches(clean_i, gregexpr("\\([^()]+\\)", clean_i))[[1L]]
      if (length(inside_i) > 0L && !identical(inside_i, character(0L))) {
        inside_i <- gsub("^\\(|\\)$", "", inside_i)
        inside_i <- metaFuns$clean_pubchem_compound_name(inside_i)
        out_i <- c(out_i, inside_i)
      }

      no_parentheses_i <- gsub("\\s*\\([^()]+\\)", "", clean_i)
      no_parentheses_i <- metaFuns$clean_pubchem_compound_name(no_parentheses_i)
      out_i <- c(out_i, no_parentheses_i)
    }

    out_i <- herbFuns$clean_text(out_i)
    out_i <- out_i[!is.na(out_i) & nzchar(out_i)]
    out_i <- unique(out_i)

    if (length(out_i) == 0L) {
      return(tibble::tibble(
        name_original = name_i,
        name_query = NA_character_,
        candidate_rank = NA_integer_
      ))
    }

    tibble::tibble(
      name_original = rep(name_i, length(out_i)),
      name_query = out_i,
      candidate_rank = seq_along(out_i)
    )
  })

  data_candidate <- dplyr::bind_rows(lst_candidate)
  data_candidate <- dplyr::filter(data_candidate, !is.na(name_query), nzchar(name_query))
  dplyr::distinct(data_candidate)
}

metaFuns$map_pubchem_cids_by_name_variants <- function(vec_name,
  dir_cache = file.path("tmp", "pubchem_herbsCollection"),
  namespace = "name", domain = "compound", searchtype = NULL,
  options = NULL, dic_name = NULL, cache_name = NULL,
  manual_cid = NULL)
{
  data_candidate <- metaFuns$make_pubchem_name_candidates(vec_name)

  if (nrow(data_candidate) == 0L) {
    out_empty <- data.frame(
      name_original = character(0L),
      name_query = character(0L),
      cid = character(0L),
      cid_n = integer(0L),
      stringsAsFactors = FALSE
    )
    attr(out_empty, "candidate_table") <- data_candidate
    attr(out_empty, "candidate_n") <- 0L
    return(out_empty)
  }

  vec_candidate <- unique(data_candidate$name_query)
  vec_candidate <- vec_candidate[!is.na(vec_candidate) & nzchar(vec_candidate)]

  if (is.null(cache_name)) {
    cache_name <- glue::glue(
      "pubchemr_cids_name_compound_variants_v1_{length(vec_candidate)}"
    )
  }

  data_cid <- metaFuns$map_pubchem_cids(
    vec_name = vec_candidate,
    dir_cache = dir_cache,
    namespace = namespace,
    domain = domain,
    searchtype = searchtype,
    options = options,
    dic_name = dic_name,
    cache_name = cache_name,
    manual_cid = manual_cid
  )

  if (is.null(data_cid) || nrow(data_cid) == 0L) {
    data_out <- data.frame(
      name_original = unique(data_candidate$name_original),
      name_query = NA_character_,
      cid = NA_character_,
      cid_n = 0L,
      stringsAsFactors = FALSE
    )
    attr(data_out, "candidate_table") <- data_candidate
    attr(data_out, "candidate_n") <- length(vec_candidate)
    return(data_out)
  }

  data_cid$name_original <- as.character(data_cid$name_original)
  data_cid$cid <- as.character(data_cid$cid)

  data_hit <- dplyr::left_join(
    data_candidate,
    dplyr::select(
      data_cid,
      name_query_candidate = name_original,
      cid,
      cid_n_raw = cid_n
    ),
    by = c("name_query" = "name_query_candidate")
  )

  lst_out <- lapply(unique(data_candidate$name_original), function(name_i) {
    data_i <- data_hit[data_hit$name_original == name_i, , drop = FALSE]
    cid_i <- unique(data_i$cid[!is.na(data_i$cid) & nzchar(data_i$cid)])
    query_i <- unique(data_i$name_query[!is.na(data_i$cid) & nzchar(data_i$cid)])

    if (length(cid_i) == 0L) {
      return(data.frame(
        name_original = name_i,
        name_query = metaFuns$clean_pubchem_compound_name(name_i),
        cid = NA_character_,
        cid_n = 0L,
        stringsAsFactors = FALSE
      ))
    }

    data.frame(
      name_original = rep(name_i, length(cid_i)),
      name_query = paste(query_i, collapse = "; "),
      cid = cid_i,
      cid_n = length(cid_i),
      stringsAsFactors = FALSE
    )
  })

  data_out <- dplyr::bind_rows(lst_out)
  data_out <- dplyr::distinct(data_out)
  attr(data_out, "candidate_table") <- data_candidate
  attr(data_out, "candidate_n") <- length(vec_candidate)
  data_out
}

metaFuns$map_pubchem_properties_by_cid <- function(cid,
  dir_cache = file.path("tmp", "pubchem_herbsCollection"),
  properties = c(
    "MolecularFormula", "MolecularWeight", "CanonicalSMILES",
    "IsomericSMILES", "InChI", "InChIKey"
  ), cache_name = NULL, prefer_isomeric_smiles = TRUE)
{
  if (exists("expect_package", mode = "function")) {
    expect_package("PubChemR", "3.0.0")
  } else {
    if (!requireNamespace("PubChemR", quietly = TRUE)) {
      stop("Package `PubChemR` is required for PubChem structure completion.")
    }
  }

  cid <- unique(trimws(as.character(cid)))
  cid <- cid[!is.na(cid) & nzchar(cid)]

  if (length(cid) == 0L) {
    return(tibble::tibble(
      pubchem_cid = character(0L),
      canonical_smiles = character(0L),
      isomeric_smiles = character(0L),
      smiles = character(0L),
      inchi = character(0L),
      inchikey = character(0L),
      molecular_formula = character(0L),
      molecular_weight = numeric(0L)
    ))
  }

  if (is.null(cache_name)) {
    cache_name <- glue::glue("pubchemr_properties_cid_v3_{length(cid)}")
  }

  if (!dir.exists(dir_cache)) {
    dir.create(dir_cache, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(dir_cache)) {
    stop("Cannot create PubChem cache directory: ", dir_cache)
  }

  object_property <- tryCatch({
    if (exists("expect_local_data", mode = "function")) {
      expect_local_data(
        dir_cache,
        cache_name,
        metaFuns$get_properties_pubchem_batch,
        list(
          identifier = cid,
          properties = properties,
          namespace = "cid",
          propertyMatch = list(.ignore.case = FALSE, type = "match")
        )
      )
    } else {
      metaFuns$get_properties_pubchem_batch(
        identifier = cid,
        properties = properties,
        namespace = "cid",
        propertyMatch = list(.ignore.case = FALSE, type = "match")
      )
    }
  }, error = function(e) {
    warning("PubChem property retrieval failed: ", conditionMessage(e))
    NULL
  })

  if (is.null(object_property)) {
    return(tibble::tibble(
      pubchem_cid = character(0L),
      canonical_smiles = character(0L),
      isomeric_smiles = character(0L),
      smiles = character(0L),
      inchi = character(0L),
      inchikey = character(0L),
      molecular_formula = character(0L),
      molecular_weight = numeric(0L)
    ))
  }

  metaFuns$extract_pubchem_properties(
    object_property = object_property,
    identifier = cid,
    prefer_isomeric_smiles = prefer_isomeric_smiles
  )
}

herbFuns$stat_structure_coverage <- function(x, label = "Structure coverage")
{
  herbFuns$validate_collection(x, stop_if_invalid = TRUE)
  data_comp <- herbFuns$add_structure_status(x$compound_unique)

  tibble::tibble(
    stage = label,
    n_unique_compound = nrow(data_comp),
    n_with_smiles = sum(data_comp$has_smiles, na.rm = TRUE),
    n_without_smiles = sum(!data_comp$has_smiles, na.rm = TRUE),
    n_with_inchikey = sum(data_comp$has_inchikey, na.rm = TRUE),
    n_without_inchikey = sum(!data_comp$has_inchikey, na.rm = TRUE),
    n_with_pubchem_cid = sum(data_comp$has_pubchem_cid, na.rm = TRUE),
    n_without_pubchem_cid = sum(!data_comp$has_pubchem_cid, na.rm = TRUE),
    n_complete_structure = sum(data_comp$structure_status == "complete_structure", na.rm = TRUE),
    n_partial_structure = sum(data_comp$structure_status == "partial_structure", na.rm = TRUE),
    n_cid_only = sum(data_comp$structure_status == "cid_only", na.rm = TRUE),
    n_name_only = sum(data_comp$structure_status == "name_only", na.rm = TRUE)
  )
}

herbFuns$fill_missing_text <- function(old, new)
{
  old <- herbFuns$clean_text(old)
  new <- herbFuns$clean_text(new)
  idx <- (is.na(old) | !nzchar(old)) & !is.na(new) & nzchar(new)
  old[idx] <- new[idx]
  old
}

herbFuns$complete_compound_structures_pubchem <- function(collection,
  dir_cache = file.path("tmp", "pubchem_herbsCollection"),
  manual_cid = NULL, dic_name = NULL, prefer_isomeric_smiles = TRUE,
  verbose = TRUE)
{
  herbFuns$validate_collection(collection, stop_if_invalid = TRUE)

  if (is.null(metaFuns$map_pubchem_cids)) {
    stop(
      "`metaFuns$map_pubchem_cids()` was not found. ",
      "Please source the PubChem helper workflow before running step2."
    )
  }

  if (!dir.exists(dir_cache)) {
    dir.create(dir_cache, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(dir_cache)) {
    stop("Cannot create PubChem cache directory: ", dir_cache)
  }

  empty_name_cid <- function() {
    data.frame(
      name_original = character(0L),
      name_query = character(0L),
      cid = character(0L),
      cid_n = integer(0L),
      stringsAsFactors = FALSE
    )
  }

  apply_cid_by_name <- function(data_comp, data_cid, source_label) {
    if (is.null(data_cid) || nrow(data_cid) == 0L) {
      return(list(data = data_comp, n = 0L))
    }
    if (!"cid_n" %in% colnames(data_cid)) {
      data_cid <- metaFuns$recount_cid_n(data_cid)
    }
    data_cid_single <- dplyr::filter(
      data_cid,
      cid_n == 1L,
      !is.na(cid),
      nzchar(cid)
    )
    data_cid_single <- dplyr::select(
      data_cid_single,
      compound_name = name_original,
      pubchem_name_query_mapped = name_query,
      pubchem_cid_mapped = cid
    )
    data_cid_single <- dplyr::distinct(data_cid_single)

    data_comp <- dplyr::left_join(data_comp, data_cid_single, by = "compound_name")
    idx <- (!data_comp$has_pubchem_cid) &
      !is.na(data_comp$pubchem_cid_mapped) & nzchar(data_comp$pubchem_cid_mapped)
    data_comp$pubchem_cid[idx] <- data_comp$pubchem_cid_mapped[idx]
    data_comp$pubchem_cid_source[idx] <- source_label
    if (!"pubchem_name_query" %in% colnames(data_comp)) {
      data_comp$pubchem_name_query <- NA_character_
    }
    data_comp$pubchem_name_query[idx] <- data_comp$pubchem_name_query_mapped[idx]
    data_comp$pubchem_cid_mapped <- NULL
    data_comp$pubchem_name_query_mapped <- NULL
    data_comp <- herbFuns$add_structure_status(data_comp)
    list(data = data_comp, n = sum(idx, na.rm = TRUE))
  }

  apply_cid_by_identifier <- function(data_comp, data_cid, id_col, source_label) {
    if (is.null(data_cid) || nrow(data_cid) == 0L) {
      return(list(data = data_comp, n = 0L))
    }
    data_cid_single <- dplyr::filter(
      data_cid,
      cid_n == 1L,
      !is.na(cid),
      nzchar(cid)
    )
    if (nrow(data_cid_single) == 0L) {
      return(list(data = data_comp, n = 0L))
    }
    data_cid_single <- dplyr::select(
      data_cid_single,
      identifier_original,
      pubchem_cid_mapped = cid
    )
    data_cid_single <- dplyr::distinct(data_cid_single)

    data_comp$join_identifier_value <- data_comp[[id_col]]
    data_comp <- dplyr::left_join(
      data_comp,
      data_cid_single,
      by = c("join_identifier_value" = "identifier_original")
    )
    idx <- (!data_comp$has_pubchem_cid) &
      !is.na(data_comp$pubchem_cid_mapped) & nzchar(data_comp$pubchem_cid_mapped)
    data_comp$pubchem_cid[idx] <- data_comp$pubchem_cid_mapped[idx]
    data_comp$pubchem_cid_source[idx] <- source_label
    data_comp$pubchem_cid_mapped <- NULL
    data_comp$join_identifier_value <- NULL
    data_comp <- herbFuns$add_structure_status(data_comp)
    list(data = data_comp, n = sum(idx, na.rm = TRUE))
  }

  data_comp <- herbFuns$add_structure_status(collection$compound)
  data_hc <- herbFuns$ensure_relationship_columns(collection$herb_compound)

  data_before <- herbFuns$stat_structure_coverage(collection, label = "before")
  data_comp$compound_key_old <- data_comp$compound_key
  data_comp$pubchem_cid_source <- ifelse(data_comp$has_pubchem_cid, "original", NA_character_)
  data_comp$pubchem_name_query <- NA_character_
  data_comp$compound_name_cleaned <- metaFuns$clean_pubchem_compound_name(data_comp$compound_name)
  data_comp$structure_completion_source <- NA_character_

  data_need_cid_name <- dplyr::filter(data_comp, !has_pubchem_cid, has_compound_name)
  vec_name <- unique(data_need_cid_name$compound_name)
  vec_name <- vec_name[!is.na(vec_name) & nzchar(vec_name)]
  data_name_candidate <- metaFuns$make_pubchem_name_candidates(vec_name)
  n_name_candidate <- length(unique(data_name_candidate$name_query))

  herbFuns$log_progress(
    verbose,
    "PubChem CID mapping by compound names: ", herbFuns$format_n(length(vec_name)),
    " names requiring CID; ", herbFuns$format_n(n_name_candidate),
    " normalized name candidates."
  )

  if (length(vec_name) > 0L) {
    data_cid_name <- tryCatch(
      metaFuns$map_pubchem_cids_by_name_variants(
        vec_name = vec_name,
        dir_cache = dir_cache,
        namespace = "name",
        domain = "compound",
        dic_name = dic_name,
        manual_cid = manual_cid
      ),
      error = function(e) {
        warning("PubChem CID mapping by normalized name candidates failed: ", conditionMessage(e))
        empty_name_cid()
      }
    )
  } else {
    data_cid_name <- empty_name_cid()
  }

  data_apply <- apply_cid_by_name(data_comp, data_cid_name, "pubchem_name")
  data_comp <- data_apply$data
  n_cid_filled_name <- data_apply$n

  data_need_cid_inchikey <- dplyr::filter(data_comp, !has_pubchem_cid, has_inchikey)
  vec_inchikey <- unique(data_need_cid_inchikey$inchikey)
  vec_inchikey <- vec_inchikey[!is.na(vec_inchikey) & nzchar(vec_inchikey)]

  herbFuns$log_progress(
    verbose,
    "PubChem CID mapping by InChIKey: ", herbFuns$format_n(length(vec_inchikey)),
    " identifiers requiring CID."
  )

  data_cid_inchikey <- metaFuns$map_pubchem_cids_by_identifier(
    identifier = vec_inchikey,
    namespace = "inchikey",
    dir_cache = dir_cache,
    cache_name = glue::glue("pubchemr_cids_inchikey_compound_{length(vec_inchikey)}"),
    domain = "compound"
  )

  data_apply <- apply_cid_by_identifier(
    data_comp,
    data_cid_inchikey,
    id_col = "inchikey",
    source_label = "pubchem_inchikey"
  )
  data_comp <- data_apply$data
  n_cid_filled_inchikey <- data_apply$n

  data_need_cid_smiles <- dplyr::filter(data_comp, !has_pubchem_cid, has_smiles)
  vec_smiles <- unique(data_need_cid_smiles$smiles)
  vec_smiles <- vec_smiles[!is.na(vec_smiles) & nzchar(vec_smiles)]

  herbFuns$log_progress(
    verbose,
    "PubChem CID mapping by SMILES: ", herbFuns$format_n(length(vec_smiles)),
    " identifiers requiring CID."
  )

  data_cid_smiles <- metaFuns$map_pubchem_cids_by_identifier(
    identifier = vec_smiles,
    namespace = "smiles",
    dir_cache = dir_cache,
    cache_name = glue::glue("pubchemr_cids_smiles_compound_{length(vec_smiles)}"),
    domain = "compound"
  )

  data_apply <- apply_cid_by_identifier(
    data_comp,
    data_cid_smiles,
    id_col = "smiles",
    source_label = "pubchem_smiles"
  )
  data_comp <- data_apply$data
  n_cid_filled_smiles <- data_apply$n

  n_cid_filled <- n_cid_filled_name + n_cid_filled_inchikey + n_cid_filled_smiles

  cid_need <- unique(data_comp$pubchem_cid[
    data_comp$has_pubchem_cid &
      (!data_comp$has_smiles | !data_comp$has_inchikey |
         is.na(herbFuns$clean_text(data_comp$inchi)) |
         is.na(herbFuns$clean_text(data_comp$molecular_formula)))
  ])
  cid_need <- cid_need[!is.na(cid_need) & nzchar(cid_need)]

  herbFuns$log_progress(
    verbose,
    "PubChem property retrieval by CID: ", herbFuns$format_n(length(cid_need)),
    " CIDs requiring structure fields."
  )

  data_prop <- metaFuns$map_pubchem_properties_by_cid(
    cid = cid_need,
    dir_cache = dir_cache,
    prefer_isomeric_smiles = prefer_isomeric_smiles
  )

  if (nrow(data_prop) > 0L) {
    data_prop_join <- dplyr::select(
      data_prop,
      pubchem_cid,
      prop_smiles = smiles,
      prop_inchi = inchi,
      prop_inchikey = inchikey,
      prop_molecular_formula = molecular_formula,
      prop_molecular_weight = molecular_weight
    )
    data_prop_join <- dplyr::distinct(data_prop_join)
    data_comp <- dplyr::left_join(data_comp, data_prop_join, by = "pubchem_cid")
  } else {
    data_comp$prop_smiles <- NA_character_
    data_comp$prop_inchi <- NA_character_
    data_comp$prop_inchikey <- NA_character_
    data_comp$prop_molecular_formula <- NA_character_
    data_comp$prop_molecular_weight <- NA_real_
  }

  old_smiles_missing <- is.na(herbFuns$clean_text(data_comp$smiles)) |
    !nzchar(herbFuns$clean_text(data_comp$smiles))
  old_inchikey_missing <- is.na(herbFuns$clean_text(data_comp$inchikey)) |
    !nzchar(herbFuns$clean_text(data_comp$inchikey))

  data_comp$smiles <- herbFuns$fill_missing_text(data_comp$smiles, data_comp$prop_smiles)
  data_comp$inchi <- herbFuns$fill_missing_text(data_comp$inchi, data_comp$prop_inchi)
  data_comp$inchikey <- herbFuns$fill_missing_text(data_comp$inchikey, data_comp$prop_inchikey)
  data_comp$molecular_formula <- herbFuns$fill_missing_text(
    data_comp$molecular_formula,
    data_comp$prop_molecular_formula
  )

  idx_weight <- (is.na(data_comp$molecular_weight) | !is.finite(data_comp$molecular_weight)) &
    !is.na(data_comp$prop_molecular_weight) & is.finite(data_comp$prop_molecular_weight)
  data_comp$molecular_weight[idx_weight] <- data_comp$prop_molecular_weight[idx_weight]

  new_smiles_present <- !is.na(herbFuns$clean_text(data_comp$smiles)) &
    nzchar(herbFuns$clean_text(data_comp$smiles))
  new_inchikey_present <- !is.na(herbFuns$clean_text(data_comp$inchikey)) &
    nzchar(herbFuns$clean_text(data_comp$inchikey))

  n_smiles_filled <- sum(old_smiles_missing & new_smiles_present, na.rm = TRUE)
  n_inchikey_filled <- sum(old_inchikey_missing & new_inchikey_present, na.rm = TRUE)

  idx_structure <- (!is.na(data_comp$prop_smiles) & nzchar(herbFuns$clean_text(data_comp$prop_smiles))) |
    (!is.na(data_comp$prop_inchikey) & nzchar(herbFuns$clean_text(data_comp$prop_inchikey)))
  data_comp$structure_completion_source[idx_structure] <- "pubchem_cid"

  data_comp$prop_smiles <- NULL
  data_comp$prop_inchi <- NULL
  data_comp$prop_inchikey <- NULL
  data_comp$prop_molecular_formula <- NULL
  data_comp$prop_molecular_weight <- NULL

  data_comp$compound_key <- herbFuns$make_compound_key(
    inchikey = data_comp$inchikey,
    smiles = data_comp$smiles,
    pubchem_cid = data_comp$pubchem_cid,
    source = data_comp$source,
    compound_source_id = data_comp$compound_source_id,
    compound_name = data_comp$compound_name
  )

  data_key <- dplyr::select(
    data_comp,
    source,
    compound_source_id,
    compound_key_old,
    compound_key_new = compound_key
  )
  data_key <- dplyr::distinct(data_key)

  data_hc$compound_key_old <- data_hc$compound_key
  data_hc <- dplyr::left_join(
    data_hc,
    data_key,
    by = c("source", "compound_source_id", "compound_key_old")
  )
  data_hc$compound_key <- ifelse(
    is.na(data_hc$compound_key_new),
    data_hc$compound_key,
    data_hc$compound_key_new
  )
  data_hc$compound_key_old <- NULL
  data_hc$compound_key_new <- NULL

  data_comp$compound_key_old <- NULL
  data_comp <- herbFuns$add_structure_status(data_comp)
  data_comp <- dplyr::distinct(data_comp)
  data_hc <- dplyr::distinct(data_hc)

  compound_unique <- herbFuns$deduplicate_compounds(data_comp)
  compound_unique <- herbFuns$add_structure_status(compound_unique)

  collection_out <- collection
  collection_out$herb_compound <- data_hc
  collection_out$compound <- data_comp
  collection_out$compound_unique <- compound_unique
  collection_out$logs <- herbFuns$make_logs(collection_out)
  class(collection_out) <- "herbs_collection"

  data_after <- herbFuns$stat_structure_coverage(collection_out, label = "after")

  data_unresolved <- dplyr::filter(compound_unique, !has_smiles)
  if (nrow(data_unresolved) > 0L) {
    data_unresolved$unresolved_reason <- ifelse(
      !data_unresolved$has_compound_name & !data_unresolved$has_inchikey & !data_unresolved$has_pubchem_cid,
      "missing_searchable_identifier",
      ifelse(
        !data_unresolved$has_pubchem_cid,
        "missing_pubchem_cid_or_multi_cid",
        "pubchem_property_missing_smiles"
      )
    )
    data_unresolved <- dplyr::select(
      data_unresolved,
      compound_key,
      compound_name,
      compound_name_cleaned,
      pubchem_name_query,
      source_list,
      pubchem_cid,
      smiles,
      inchikey,
      structure_status,
      unresolved_reason
    )
  }

  diag_cid_name <- if (nrow(data_cid_name) > 0L) {
    metaFuns$mapping_diagnostics(data_cid_name, label = "name")
  } else {
    list(
      summary = data.frame(
        label = "name",
        n_name = 0L,
        n_mapped = 0L,
        n_unmapped = 0L,
        n_multi_cid = 0L,
        stringsAsFactors = FALSE
      ),
      unmapped = tibble::tibble(),
      multi_cid = tibble::tibble()
    )
  }

  summarise_identifier_cid <- function(data_cid, label) {
    if (is.null(data_cid) || nrow(data_cid) == 0L) {
      return(data.frame(
        label = label,
        n_name = 0L,
        n_mapped = 0L,
        n_unmapped = 0L,
        n_multi_cid = 0L,
        stringsAsFactors = FALSE
      ))
    }
    data_by_id <- dplyr::group_by(data_cid, identifier_original)
    data_by_id <- dplyr::summarise(
      data_by_id,
      cid_n = sum(!is.na(cid) & nzchar(cid)),
      .groups = "drop"
    )
    data.frame(
      label = label,
      n_name = nrow(data_by_id),
      n_mapped = sum(data_by_id$cid_n > 0L),
      n_unmapped = sum(data_by_id$cid_n == 0L),
      n_multi_cid = sum(data_by_id$cid_n > 1L),
      stringsAsFactors = FALSE
    )
  }

  data_summary <- tibble::tibble(
    n_unique_compound_before = data_before$n_unique_compound[[1L]],
    n_unique_compound_after = data_after$n_unique_compound[[1L]],
    n_pubchem_cid_before = data_before$n_with_pubchem_cid[[1L]],
    n_pubchem_cid_after = data_after$n_with_pubchem_cid[[1L]],
    n_smiles_before = data_before$n_with_smiles[[1L]],
    n_smiles_after = data_after$n_with_smiles[[1L]],
    n_inchikey_before = data_before$n_with_inchikey[[1L]],
    n_inchikey_after = data_after$n_with_inchikey[[1L]],
    n_cid_filled = n_cid_filled,
    n_cid_filled_by_name = n_cid_filled_name,
    n_cid_filled_by_inchikey = n_cid_filled_inchikey,
    n_cid_filled_by_smiles = n_cid_filled_smiles,
    n_smiles_filled = n_smiles_filled,
    n_inchikey_filled = n_inchikey_filled,
    n_cid_query_name = length(vec_name),
    n_cid_query_name_candidate = n_name_candidate,
    n_cid_mapped_name = diag_cid_name$summary$n_mapped[[1L]],
    n_cid_multi_match_name = diag_cid_name$summary$n_multi_cid[[1L]],
    n_cid_unmapped_name = diag_cid_name$summary$n_unmapped[[1L]],
    n_cid_query_inchikey = length(vec_inchikey),
    n_cid_mapped_inchikey = summarise_identifier_cid(data_cid_inchikey, "inchikey")$n_mapped[[1L]],
    n_cid_multi_match_inchikey = summarise_identifier_cid(data_cid_inchikey, "inchikey")$n_multi_cid[[1L]],
    n_cid_unmapped_inchikey = summarise_identifier_cid(data_cid_inchikey, "inchikey")$n_unmapped[[1L]],
    n_cid_query_smiles = length(vec_smiles),
    n_cid_mapped_smiles = summarise_identifier_cid(data_cid_smiles, "smiles")$n_mapped[[1L]],
    n_cid_multi_match_smiles = summarise_identifier_cid(data_cid_smiles, "smiles")$n_multi_cid[[1L]],
    n_cid_unmapped_smiles = summarise_identifier_cid(data_cid_smiles, "smiles")$n_unmapped[[1L]],
    n_property_query_cid = length(cid_need),
    n_property_mapped_cid = dplyr::n_distinct(data_prop$pubchem_cid),
    n_no_smiles_compound = nrow(data_unresolved)
  )

  data_cid_inchikey_out <- dplyr::mutate(data_cid_inchikey, cid_source = "pubchem_inchikey")
  data_cid_smiles_out <- dplyr::mutate(data_cid_smiles, cid_source = "pubchem_smiles")
  data_cid_name_out <- dplyr::mutate(data_cid_name, cid_source = "pubchem_name")

  list(
    collection = collection_out,
    summary = data_summary,
    cid_mapping = data_cid_name_out,
    cid_mapping_inchikey = data_cid_inchikey_out,
    cid_mapping_smiles = data_cid_smiles_out,
    cid_diagnostics = diag_cid_name,
    property_mapping = data_prop,
    unresolved = data_unresolved
  )
}

# --------------------------------------------------------------------------
# File helpers and readers
# --------------------------------------------------------------------------

herbFuns$resolve_db_files <- function(files, db_dir)
{
  if (is.null(names(files)) || any(!nzchar(names(files)))) {
    names(files) <- basename(files)
  }

  vec_file <- as.character(files)
  is_abs <- grepl("^/", vec_file) | grepl("^[A-Za-z]:", vec_file)
  out <- ifelse(is_abs, vec_file, file.path(pg(db_dir), vec_file))
  names(out) <- names(files)

  out
}

herbFuns$files_tcmsp <- function(files)
{
  keyCols <- list(
    "Herb_Ingredients_relationship.xlsx" = c(
      "herb_cn_name", "herb_en_name", "MOL_ID"
    ),
    "Ingredients_Targets_relationship.xlsx" = c(
      "herb_cn_name", "herb_en_name", "MOL_ID",
      "target_name", "gene_name", "target_ID"
    ),
    "cid_info.csv" = c("Mol ID", "Smiles")
  )

  herbFuns$resolve_db_files(files, "db_local_tcmsp")
}

herbFuns$files_batman <- function(files)
{
  keyCols <- list(
    "herb_browse.txt" = c("Chinese.Name", "English.Name", "Ingredients"),
    "cid_info.csv" = c("Compound_CID", "SMILES"),
    "known_browse_by_ingredients.txt.gz" = c(
      "PubChem_CID", "IUPAC_name", "known_target_proteins"
    ),
    "predicted_browse_by_ingredients.txt.gz" = c(
      "PubChem_CID", "IUPAC_name", "predicted_target_proteins"
    )
  )

  herbFuns$resolve_db_files(files, "db_local_batman")
}

herbFuns$read_table <- function(file)
{
  if (!file.exists(file)) {
    stop("File was not found: ", file)
  }

  ext <- tolower(tools::file_ext(file))

  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Package `readxl` is required to read Excel files.")
    }
    data <- readxl::read_excel(file, col_types = "text")
    return(as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE))
  }

  if (ext %in% c("csv")) {
    return(utils::read.csv(file, stringsAsFactors = FALSE, check.names = FALSE))
  }

  if (ext %in% c("tsv", "txt")) {
    return(utils::read.delim(file, stringsAsFactors = FALSE, check.names = FALSE))
  }

  stop("Unsupported file extension: ", ext)
}

herbFuns$read_tcmsp <- function(files = NULL, verbose = TRUE)
{
  if (is.null(files)) {
    files <- c(
      "Herb_Ingredients_relationship.xlsx",
      "cid_info.csv"
    )
    names(files) <- files
    files <- herbFuns$files_tcmsp(files)
  }

  herbFuns$log_progress(verbose, "Reading TCMSP herb-compound records.")
  data_rel <- herbFuns$read_table(files[["Herb_Ingredients_relationship.xlsx"]])

  herbFuns$log_progress(verbose, "Reading TCMSP compound structure records.")
  data_cid <- herbFuns$read_table(files[["cid_info.csv"]])

  herbFuns$log_progress(
    verbose,
    "TCMSP loaded: ", herbFuns$format_n(nrow(data_rel)),
    " herb-compound rows; ", herbFuns$format_n(nrow(data_cid)),
    " compound rows."
  )

  list(
    herb_ingredients = data_rel,
    cid_info = data_cid
  )
}

herbFuns$read_batman <- function(files = NULL, verbose = TRUE)
{
  if (is.null(files)) {
    files <- c("herb_browse.txt", "cid_info.csv")
    names(files) <- files
    files <- herbFuns$files_batman(files)
  }

  herbFuns$log_progress(verbose, "Reading BATMAN-TCM herb-compound records.")
  data_herb <- herbFuns$read_table(files[["herb_browse.txt"]])

  herbFuns$log_progress(verbose, "Reading BATMAN-TCM compound structure records.")
  data_cid <- herbFuns$read_table(files[["cid_info.csv"]])

  herbFuns$log_progress(
    verbose,
    "BATMAN-TCM loaded: ", herbFuns$format_n(nrow(data_herb)),
    " herb rows; ", herbFuns$format_n(nrow(data_cid)),
    " compound rows."
  )

  list(
    herb_browse = data_herb,
    cid_info = data_cid
  )
}


herbFuns$read_tcmsp_targets <- function(files = NULL, verbose = TRUE)
{
  if (is.null(files)) {
    files <- c("Ingredients_Targets_relationship.xlsx")
    names(files) <- files
    files <- herbFuns$files_tcmsp(files)
  }

  herbFuns$log_progress(verbose, "Reading TCMSP compound-target records.")
  data_target <- herbFuns$read_table(files[["Ingredients_Targets_relationship.xlsx"]])

  herbFuns$log_progress(
    verbose,
    "TCMSP target records loaded: ", herbFuns$format_n(nrow(data_target)),
    " compound-target rows."
  )

  data_target
}

herbFuns$read_batman_tsv <- function(file)
{
  if (!file.exists(file)) {
    stop("File was not found: ", file)
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package `data.table` is required to read BATMAN-TCM target files.")
  }

  data.table::fread(
    file,
    sep = "\t",
    data.table = FALSE,
    quote = "",
    fill = TRUE,
    showProgress = FALSE
  )
}

herbFuns$parse_batman_predicted_table <- function(data)
{
  data <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)

  if (ncol(data) == 1L) {
    vec_line <- herbFuns$clean_text(data[[1L]])
    vec_line <- vec_line[!is.na(vec_line) & nzchar(vec_line)]
    pubchem_cid <- sub("^([0-9]+)\\s+.*$", "\\1", vec_line)
    pubchem_cid[!grepl("^[0-9]+$", pubchem_cid)] <- NA_character_
    predicted_target_proteins <- sub("^.*\\s([^\\s]+)$", "\\1", vec_line)
    predicted_target_proteins[predicted_target_proteins == vec_line] <- NA_character_
    iupac_name <- sub("^[^\\s]+\\s+(.*)\\s+[^\\s]+$", "\\1", vec_line)
    iupac_name[iupac_name == vec_line] <- NA_character_

    return(tibble::tibble(
      PubChem_CID = pubchem_cid,
      IUPAC_name = iupac_name,
      predicted_target_proteins = predicted_target_proteins
    ))
  }

  col_cid <- herbFuns$resolve_col(data, c("PubChem_CID", "PubChem CID", "cid"), required = TRUE)
  col_iupac <- herbFuns$resolve_col(data, c("IUPAC_name", "IUPAC name", "IUPAC"), required = FALSE)
  col_target <- herbFuns$resolve_col(
    data,
    c("predicted_target_proteins", "predicted target proteins", "Targets"),
    required = TRUE
  )

  tibble::tibble(
    PubChem_CID = herbFuns$clean_text(data[[col_cid]]),
    IUPAC_name = herbFuns$clean_text(herbFuns$get_col(data, col_iupac)),
    predicted_target_proteins = herbFuns$clean_text(data[[col_target]])
  )
}

herbFuns$read_batman_targets <- function(files = NULL, verbose = TRUE)
{
  if (is.null(files)) {
    files <- c(
      "known_browse_by_ingredients.txt.gz",
      "predicted_browse_by_ingredients.txt.gz"
    )
    names(files) <- files
    files <- herbFuns$files_batman(files)
  }

  herbFuns$log_progress(verbose, "Reading BATMAN-TCM known compound-target records.")
  data_known <- herbFuns$read_batman_tsv(files[["known_browse_by_ingredients.txt.gz"]])

  herbFuns$log_progress(verbose, "Reading BATMAN-TCM predicted compound-target records.")
  data_predicted <- herbFuns$read_batman_tsv(files[["predicted_browse_by_ingredients.txt.gz"]])
  data_predicted <- herbFuns$parse_batman_predicted_table(data_predicted)

  herbFuns$log_progress(
    verbose,
    "BATMAN-TCM target records loaded: ", herbFuns$format_n(nrow(data_known)),
    " known compound rows; ", herbFuns$format_n(nrow(data_predicted)),
    " predicted compound rows."
  )

  list(
    known = data_known,
    predicted = data_predicted
  )
}


# --------------------------------------------------------------------------
# Matching helpers
# --------------------------------------------------------------------------

herbFuns$match_herb_records <- function(data_herb, query, source,
  fields, field_labels = NULL)
{
  if (is.null(field_labels)) {
    field_labels <- fields
  }

  if (nrow(query) == 0L || nrow(data_herb) == 0L) {
    return(tibble::tibble())
  }

  lst_match <- list()
  k <- 0L

  for (i in seq_len(nrow(query))) {
    for (j in seq_along(fields)) {
      field <- fields[[j]]

      if (!field %in% colnames(data_herb)) {
        next
      }

      key <- herbFuns$match_key(data_herb[[field]])
      idx <- which(!is.na(key) & key == query$query_key[[i]])

      if (length(idx) == 0L) {
        next
      }

      k <- k + 1L
      data_hit <- data_herb[idx, , drop = FALSE]
      data_hit$query_id <- query$query_id[[i]]
      data_hit$query_herb <- query$query_herb[[i]]
      data_hit$source <- source
      data_hit$match_field <- field_labels[[j]]
      data_hit$match_type <- "exact"
      data_hit$match_score <- 1
      data_hit$match_priority <- j
      lst_match[[k]] <- data_hit
    }
  }

  if (length(lst_match) == 0L) {
    return(tibble::tibble())
  }

  data_match <- dplyr::bind_rows(lst_match)
  data_match <- dplyr::arrange(data_match, query_id, source_herb_id, match_priority)
  data_match <- dplyr::distinct(
    data_match,
    query_id,
    source,
    source_herb_id,
    .keep_all = TRUE
  )

  data_match
}

herbFuns$filter_herb_candidate_records <- function(data_herb, query, fields)
{
  if (nrow(data_herb) == 0L || nrow(query) == 0L) {
    return(data_herb[0L, , drop = FALSE])
  }

  query_key <- unique(query$query_key)
  query_key <- query_key[!is.na(query_key)]
  keep <- rep(FALSE, nrow(data_herb))

  for (field in fields) {
    if (!field %in% colnames(data_herb)) {
      next
    }

    key <- herbFuns$match_key(data_herb[[field]])
    keep <- keep | (!is.na(key) & key %in% query_key)
  }

  data_herb[keep, , drop = FALSE]
}

# --------------------------------------------------------------------------
# BATMAN-TCM duplicate herb-name handling
# --------------------------------------------------------------------------

herbFuns$diagnose_duplicate_herb_names <- function(data_herb,
  group_col = "herb_cn_name", source = "batman",
  detail_cols = c("herb_pinyin", "herb_en_name", "herb_latin_name"))
{
  if (!group_col %in% colnames(data_herb) || nrow(data_herb) == 0L) {
    return(tibble::tibble())
  }

  data_work <- data_herb
  data_work$duplicate_group_key <- herbFuns$match_key(data_work[[group_col]])
  data_work$duplicate_group_value <- data_work[[group_col]]
  data_work <- dplyr::filter(data_work, !is.na(duplicate_group_key))

  if (nrow(data_work) == 0L) {
    return(tibble::tibble())
  }

  detail_cols <- detail_cols[detail_cols %in% colnames(data_work)]
  data_dup <- dplyr::group_by(data_work, duplicate_group_key)
  data_dup <- dplyr::summarize(
    data_dup,
    source = source,
    herb_name = herbFuns$first_non_missing(duplicate_group_value),
    n_records = dplyr::n(),
    source_herb_ids = herbFuns$collapse_unique(source_herb_id),
    .groups = "drop"
  )

  if (length(detail_cols) > 0L) {
    for (col in detail_cols) {
      data_tmp <- data_work
      data_tmp$detail_value <- data_tmp[[col]]
      data_col <- dplyr::group_by(data_tmp, duplicate_group_key)
      data_col <- dplyr::summarize(
        data_col,
        value = herbFuns$collapse_unique(detail_value),
        n_value = dplyr::n_distinct(detail_value[!is.na(detail_value)]),
        .groups = "drop"
      )
      colnames(data_col)[colnames(data_col) == "value"] <- paste0(col, "_values")
      colnames(data_col)[colnames(data_col) == "n_value"] <- paste0("n_", col)
      data_dup <- dplyr::left_join(data_dup, data_col, by = "duplicate_group_key")
    }
  }

  data_dup <- dplyr::filter(data_dup, n_records > 1L)

  if (nrow(data_dup) == 0L) {
    return(data_dup)
  }

  consistency_cols <- c("herb_en_name", "herb_latin_name")
  data_dup$is_consistent <- TRUE

  for (col in consistency_cols) {
    n_col <- paste0("n_", col)
    if (n_col %in% colnames(data_dup)) {
      data_dup$is_consistent <- data_dup$is_consistent & data_dup[[n_col]] <= 1L
    }
  }

  data_dup
}

herbFuns$print_duplicate_herb_group <- function(data_herb, group_key,
  group_col = "herb_cn_name")
{
  idx <- which(herbFuns$match_key(data_herb[[group_col]]) == group_key)
  data_show <- data_herb[idx, , drop = FALSE]
  cols <- c(
    "source_herb_id", "herb_cn_name", "herb_pinyin",
    "herb_en_name", "herb_latin_name"
  )
  cols <- cols[cols %in% colnames(data_show)]
  print(dplyr::select(data_show, dplyr::all_of(cols)), n = Inf)
  invisible(data_show)
}

herbFuns$resolve_duplicate_herb_names <- function(data_herb,
  source = "batman", group_col = "herb_cn_name",
  duplicate_policy = c("merge_if_consistent", "interactive", "keep_all", "error"))
{
  duplicate_policy <- match.arg(duplicate_policy)

  if (!"source_herb_id" %in% colnames(data_herb)) {
    stop("`data_herb` must contain `source_herb_id` before duplicate resolution.")
  }

  data_herb$source_herb_id_original <- data_herb$source_herb_id
  data_herb$duplicate_group_key <- NA_character_
  data_herb$duplicate_group_n <- 1L
  data_herb$duplicate_resolution <- "unique"

  data_dup <- herbFuns$diagnose_duplicate_herb_names(
    data_herb,
    group_col = group_col,
    source = source
  )

  if (nrow(data_dup) == 0L) {
    return(list(
      data_herb = data_herb,
      duplicate_herb_names = data_dup,
      duplicate_decisions = tibble::tibble()
    ))
  }

  if (identical(duplicate_policy, "error")) {
    print(data_dup)
    stop("Duplicated herb names were found in `", source, "`.")
  }

  if (identical(duplicate_policy, "interactive") && !base::interactive()) {
    warning(
      "Duplicated herb names were found, but the current session is not interactive. ",
      "Use `merge_if_consistent` fallback."
    )
    duplicate_policy <- "merge_if_consistent"
  }

  data_decision <- data_dup
  data_decision$merge_decision <- FALSE
  data_decision$decision_reason <- "kept_separate"

  if (identical(duplicate_policy, "keep_all")) {
    data_decision$decision_reason <- "keep_all_policy"
  }

  if (identical(duplicate_policy, "merge_if_consistent")) {
    data_decision$merge_decision <- data_decision$is_consistent
    data_decision$decision_reason <- ifelse(
      data_decision$is_consistent,
      "merged_consistent_english_latin",
      "kept_separate_inconsistent_english_latin"
    )
  }

  if (identical(duplicate_policy, "interactive")) {
    for (i in seq_len(nrow(data_decision))) {
      group_key <- data_decision$duplicate_group_key[[i]]
      cat("\nDuplicated herb name in ", source, ": ",
        data_decision$herb_name[[i]], "\n", sep = "")
      herbFuns$print_duplicate_herb_group(
        data_herb,
        group_key = group_key,
        group_col = group_col
      )
      ans <- readline("Merge these records into one source herb? [y/N]: ")
      ans <- tolower(trimws(ans))
      data_decision$merge_decision[[i]] <- ans %in% c("y", "yes")
      data_decision$decision_reason[[i]] <- ifelse(
        data_decision$merge_decision[[i]],
        "merged_interactive_confirmation",
        "kept_separate_interactive_confirmation"
      )
    }
  }

  for (i in seq_len(nrow(data_decision))) {
    group_key <- data_decision$duplicate_group_key[[i]]
    idx <- which(herbFuns$match_key(data_herb[[group_col]]) == group_key)

    if (length(idx) == 0L) {
      next
    }

    data_herb$duplicate_group_key[idx] <- group_key
    data_herb$duplicate_group_n[idx] <- length(idx)

    if (isTRUE(data_decision$merge_decision[[i]])) {
      merged_id <- herbFuns$make_source_herb_id(
        source,
        "merged",
        herbFuns$first_non_missing(data_herb[[group_col]][idx]),
        herbFuns$first_non_missing(data_herb$herb_en_name[idx]),
        herbFuns$first_non_missing(data_herb$herb_latin_name[idx])
      )
      data_herb$source_herb_id[idx] <- merged_id
      data_herb$duplicate_resolution[idx] <- "merged"
    } else {
      data_herb$duplicate_resolution[idx] <- "kept_separate"
    }
  }

  list(
    data_herb = data_herb,
    duplicate_herb_names = data_dup,
    duplicate_decisions = data_decision
  )
}

# --------------------------------------------------------------------------
# TCMSP normalizer
# --------------------------------------------------------------------------

herbFuns$normalize_tcmsp <- function(data_tcmsp, herbs, verbose = TRUE)
{
  query <- herbFuns$prepare_query(herbs)
  herbFuns$log_progress(verbose, "Normalizing TCMSP records for ", length(herbs), " query herbs.")

  data_rel <- data_tcmsp$herb_ingredients
  col_cn <- herbFuns$resolve_col(data_rel, c("herb_cn_name", "Chinese.Name"), required = TRUE)
  col_en <- herbFuns$resolve_col(data_rel, c("herb_en_name", "English.Name"), required = FALSE)
  col_pinyin <- herbFuns$resolve_col(data_rel, c("herb_pinyin", "Pinyin.Name"), required = FALSE)
  col_mol <- herbFuns$resolve_col(data_rel, c("MOL_ID", "Mol ID", "molecule_ID"), required = TRUE)
  col_name <- herbFuns$resolve_col(data_rel, c("molecule_name", "Name", "compound_name"), required = FALSE)

  data_rel_std <- tibble::tibble(
    herb_cn_name = herbFuns$clean_text(data_rel[[col_cn]]),
    herb_en_name = herbFuns$clean_text(herbFuns$get_col(data_rel, col_en)),
    herb_pinyin = herbFuns$clean_text(herbFuns$get_col(data_rel, col_pinyin)),
    herb_latin_name = NA_character_,
    compound_source_id = herbFuns$clean_text(data_rel[[col_mol]]),
    compound_name = herbFuns$clean_text(herbFuns$get_col(data_rel, col_name))
  )

  data_rel_std$source <- "tcmsp"
  data_rel_std$source_herb_id <- vapply(seq_len(nrow(data_rel_std)),
    function(i) {
      herbFuns$make_source_herb_id(
        "tcmsp",
        data_rel_std$herb_cn_name[[i]],
        data_rel_std$herb_pinyin[[i]],
        data_rel_std$herb_en_name[[i]]
      )
    }, character(1L))

  data_herb <- dplyr::distinct(
    data_rel_std,
    source_herb_id,
    herb_cn_name,
    herb_en_name,
    herb_pinyin,
    herb_latin_name
  )

  data_match <- herbFuns$match_herb_records(
    data_herb,
    query = query,
    source = "tcmsp",
    fields = c("herb_cn_name", "herb_en_name", "herb_pinyin"),
    field_labels = c("herb_cn_name", "herb_en_name", "herb_pinyin")
  )

  if (nrow(data_match) == 0L) {
    herbFuns$log_progress(verbose, "TCMSP matched 0 query herbs.")
    return(herbFuns$new_empty_source_result(query))
  }

  herbFuns$log_progress(
    verbose,
    "TCMSP matched ", herbFuns$format_n(dplyr::n_distinct(data_match$query_herb)),
    " query herbs and ", herbFuns$format_n(dplyr::n_distinct(data_match$source_herb_id)),
    " source herb records."
  )

  herb_map <- dplyr::select(
    data_match,
    query_id,
    query_herb,
    source,
    source_herb_id,
    herb_cn_name,
    herb_en_name,
    herb_pinyin,
    herb_latin_name,
    match_field,
    match_type,
    match_score,
    dplyr::any_of(c(
      "source_herb_id_original", "duplicate_group_key",
      "duplicate_group_n", "duplicate_resolution"
    ))
  )

  data_rel_hit <- dplyr::inner_join(
    data_rel_std,
    dplyr::select(
      herb_map,
      query_id,
      query_herb,
      source,
      source_herb_id,
      herb_cn_name,
      herb_en_name,
      herb_pinyin,
      herb_latin_name
    ),
    by = c(
      "source", "source_herb_id", "herb_cn_name", "herb_en_name",
      "herb_pinyin", "herb_latin_name"
    )
  )

  herb_compound <- dplyr::select(
    data_rel_hit,
    query_id,
    query_herb,
    source,
    source_herb_id,
    herb_cn_name,
    herb_en_name,
    herb_pinyin,
    herb_latin_name,
    compound_source_id,
    compound_name
  )

  herb_compound$compound_name_source <- "tcmsp"
  herb_compound$relationship_source <- "tcmsp"
  herb_compound$evidence_source <- "TCMSP"

  data_cid <- data_tcmsp$cid_info
  col_cid <- herbFuns$resolve_col(data_cid, c("Mol ID", "MOL_ID", "mol_id"), required = TRUE)
  col_smiles <- herbFuns$resolve_col(data_cid, c("Smiles", "SMILES", "smiles"), required = FALSE)
  col_inchikey <- herbFuns$resolve_col(data_cid, c("InChIKey", "inchikey"), required = FALSE)
  col_inchi <- herbFuns$resolve_col(data_cid, c("InChI", "inchi"), required = FALSE)
  col_formula <- herbFuns$resolve_col(data_cid, c("Molecular Formula", "molecular_formula", "Formula"), required = FALSE)
  col_weight <- herbFuns$resolve_col(data_cid, c("Molecular Weight", "Molecular_Weight", "Weight"), required = FALSE)

  compound <- tibble::tibble(
    source = "tcmsp",
    compound_source_id = herbFuns$clean_text(data_cid[[col_cid]]),
    compound_name = NA_character_,
    mol_id = herbFuns$clean_text(data_cid[[col_cid]]),
    pubchem_cid = NA_character_,
    smiles = herbFuns$clean_text(herbFuns$get_col(data_cid, col_smiles)),
    inchi = herbFuns$clean_text(herbFuns$get_col(data_cid, col_inchi)),
    inchikey = herbFuns$clean_text(herbFuns$get_col(data_cid, col_inchikey)),
    molecular_formula = herbFuns$clean_text(herbFuns$get_col(data_cid, col_formula)),
    molecular_weight = suppressWarnings(as.numeric(herbFuns$get_col(data_cid, col_weight)))
  )

  compound <- dplyr::filter(compound, compound_source_id %in% herb_compound$compound_source_id)

  compound_name_map <- dplyr::group_by(herb_compound, compound_source_id)
  compound_name_map <- dplyr::summarise(
    compound_name_map,
    compound_name_rel = herbFuns$first_non_missing(compound_name),
    .groups = "drop"
  )
  compound <- dplyr::left_join(compound, compound_name_map, by = "compound_source_id")
  compound$compound_name <- ifelse(
    is.na(compound$compound_name) | !nzchar(compound$compound_name),
    compound$compound_name_rel,
    compound$compound_name
  )
  compound$compound_name_rel <- NULL
  compound$compound_key <- herbFuns$make_compound_key(
    inchikey = compound$inchikey,
    smiles = compound$smiles,
    pubchem_cid = compound$pubchem_cid,
    source = compound$source,
    compound_source_id = compound$compound_source_id,
    compound_name = compound$compound_name
  )
  compound$source_priority <- 1L
  compound <- herbFuns$add_structure_status(compound)

  herb_compound <- dplyr::left_join(
    herb_compound,
    dplyr::select(compound, source, compound_source_id, compound_key),
    by = c("source", "compound_source_id")
  )

  herbFuns$log_progress(
    verbose,
    "TCMSP collected ", herbFuns$format_n(nrow(herb_compound)),
    " herb-compound relationships and ", herbFuns$format_n(nrow(compound)),
    " compound structure rows."
  )

  list(
    query = query,
    herb_map = herb_map,
    herb_compound = herb_compound,
    compound = compound,
    logs = list()
  )
}

# --------------------------------------------------------------------------
# BATMAN-TCM normalizer
# --------------------------------------------------------------------------

herbFuns$parse_batman_ingredient_item <- function(item)
{
  item <- herbFuns$clean_text(item)

  if (is.na(item)) {
    return(tibble::tibble())
  }

  item <- gsub("^[【\\[]", "", item)
  item <- gsub("[】\\]]$", "", item)
  item <- trimws(item)

  if (grepl("^NULL\\s*\\(\\s*NA\\s*\\)$", item, ignore.case = TRUE) ||
      grepl("^NA\\s*\\(\\s*NA\\s*\\)$", item, ignore.case = TRUE) ||
      grepl("^NULL$", item, ignore.case = TRUE) ||
      grepl("^NA$", item, ignore.case = TRUE)) {
    return(tibble::tibble())
  }

  has_cid <- grepl("\\(([0-9]+)\\)\\s*$", item)

  if (!has_cid) {
    return(tibble::tibble(
      ingredient_item = item,
      compound_name = item,
      compound_source_id = NA_character_,
      pubchem_cid = NA_character_
    ))
  }

  cid <- sub("^.*\\(([0-9]+)\\)\\s*$", "\\1", item)
  name <- sub("\\s*\\([0-9]+\\)\\s*$", "", item)

  tibble::tibble(
    ingredient_item = item,
    compound_name = herbFuns$clean_text(name),
    compound_source_id = herbFuns$clean_text(cid),
    pubchem_cid = herbFuns$clean_text(cid)
  )
}

herbFuns$parse_batman_ingredients <- function(data_herb)
{
  if (!"Ingredients" %in% colnames(data_herb)) {
    stop("BATMAN herb table must contain `Ingredients`.")
  }

  lst <- lapply(seq_len(nrow(data_herb)), function(i) {
    ingredients <- herbFuns$clean_text(data_herb$Ingredients[[i]])

    if (is.na(ingredients)) {
      return(tibble::tibble())
    }

    items <- unlist(strsplit(ingredients, "\\|"), use.names = FALSE)
    data_item <- dplyr::bind_rows(lapply(items, herbFuns$parse_batman_ingredient_item))

    if (nrow(data_item) == 0L) {
      return(tibble::tibble())
    }

    data_item$source_herb_id <- data_herb$source_herb_id[[i]]
    data_item
  })

  dplyr::bind_rows(lst)
}

herbFuns$normalize_batman <- function(data_batman, herbs,
  duplicate_policy = c("merge_if_consistent", "interactive", "keep_all", "error"),
  verbose = TRUE)
{
  duplicate_policy <- match.arg(duplicate_policy)
  query <- herbFuns$prepare_query(herbs)
  herbFuns$log_progress(verbose, "Normalizing BATMAN-TCM records for ", length(herbs), " query herbs.")

  data_herb_raw <- data_batman$herb_browse
  col_cn <- herbFuns$resolve_col(data_herb_raw, c("Chinese.Name", "herb_cn_name"), required = TRUE)
  col_en <- herbFuns$resolve_col(data_herb_raw, c("English.Name", "herb_en_name"), required = FALSE)
  col_pinyin <- herbFuns$resolve_col(data_herb_raw, c("Pinyin.Name", "herb_pinyin"), required = FALSE)
  col_latin <- herbFuns$resolve_col(data_herb_raw, c("Latin.Name", "herb_latin_name"), required = FALSE)

  data_herb <- data.frame(
    herb_cn_name = herbFuns$clean_text(data_herb_raw[[col_cn]]),
    herb_en_name = herbFuns$clean_text(herbFuns$get_col(data_herb_raw, col_en)),
    herb_pinyin = herbFuns$clean_text(herbFuns$get_col(data_herb_raw, col_pinyin)),
    herb_latin_name = herbFuns$clean_text(herbFuns$get_col(data_herb_raw, col_latin)),
    Ingredients = herbFuns$clean_text(data_herb_raw$Ingredients),
    stringsAsFactors = FALSE
  )

  data_herb$source <- "batman"
  data_herb$source_herb_id <- vapply(seq_len(nrow(data_herb)),
    function(i) {
      herbFuns$make_source_herb_id(
        "batman",
        data_herb$herb_cn_name[[i]],
        data_herb$herb_pinyin[[i]],
        data_herb$herb_en_name[[i]],
        data_herb$herb_latin_name[[i]]
      )
    }, character(1L))

  match_fields <- c("herb_cn_name", "herb_en_name", "herb_pinyin", "herb_latin_name")
  match_labels <- c("Chinese.Name", "English.Name", "Pinyin.Name", "Latin.Name")
  data_herb <- herbFuns$filter_herb_candidate_records(
    data_herb,
    query = query,
    fields = match_fields
  )

  herbFuns$log_progress(
    verbose,
    "BATMAN-TCM candidate herb records after query matching: ",
    herbFuns$format_n(nrow(data_herb)),
    "."
  )

  if (nrow(data_herb) == 0L) {
    herbFuns$log_progress(verbose, "BATMAN-TCM matched 0 query herbs.")
    return(herbFuns$new_empty_source_result(query))
  }

  duplicate_info <- herbFuns$resolve_duplicate_herb_names(
    data_herb,
    source = "batman",
    group_col = "herb_cn_name",
    duplicate_policy = duplicate_policy
  )
  data_herb <- duplicate_info$data_herb

  herbFuns$log_progress(
    verbose,
    "BATMAN-TCM duplicated Chinese-name groups among matched records: ",
    herbFuns$format_n(nrow(duplicate_info$duplicate_herb_names)),
    "; merged groups: ",
    herbFuns$format_n(sum(duplicate_info$duplicate_decisions$merge_decision, na.rm = TRUE)),
    "."
  )

  data_match <- herbFuns$match_herb_records(
    data_herb,
    query = query,
    source = "batman",
    fields = match_fields,
    field_labels = match_labels
  )

  if (nrow(data_match) == 0L) {
    herbFuns$log_progress(verbose, "BATMAN-TCM matched 0 query herbs.")
    return(herbFuns$new_empty_source_result(query))
  }

  herbFuns$log_progress(
    verbose,
    "BATMAN-TCM matched ", herbFuns$format_n(dplyr::n_distinct(data_match$query_herb)),
    " query herbs and ", herbFuns$format_n(dplyr::n_distinct(data_match$source_herb_id)),
    " source herb records."
  )

  herb_map <- dplyr::select(
    data_match,
    query_id,
    query_herb,
    source,
    source_herb_id,
    herb_cn_name,
    herb_en_name,
    herb_pinyin,
    herb_latin_name,
    match_field,
    match_type,
    match_score,
    dplyr::any_of(c(
      "source_herb_id_original", "duplicate_group_key",
      "duplicate_group_n", "duplicate_resolution"
    ))
  )

  data_herb_hit <- dplyr::semi_join(
    data_herb,
    dplyr::select(herb_map, source_herb_id),
    by = "source_herb_id"
  )

  herbFuns$log_progress(
    verbose,
    "Parsing BATMAN-TCM Ingredients for ", herbFuns$format_n(nrow(data_herb_hit)),
    " matched source herb rows only."
  )

  data_ingredient <- herbFuns$parse_batman_ingredients(data_herb_hit)

  herbFuns$log_progress(
    verbose,
    "BATMAN-TCM parsed ", herbFuns$format_n(nrow(data_ingredient)),
    " ingredient-CID rows after removing empty NULL/NA items."
  )

  data_rel_hit <- dplyr::inner_join(
    data_ingredient,
    dplyr::select(
      herb_map,
      query_id,
      query_herb,
      source,
      source_herb_id,
      herb_cn_name,
      herb_en_name,
      herb_pinyin,
      herb_latin_name
    ),
    by = "source_herb_id"
  )

  herb_compound <- dplyr::select(
    data_rel_hit,
    query_id,
    query_herb,
    source,
    source_herb_id,
    herb_cn_name,
    herb_en_name,
    herb_pinyin,
    herb_latin_name,
    compound_source_id,
    pubchem_cid,
    compound_name
  )

  herb_compound$compound_name_source <- "batman_ingredients"
  herb_compound$relationship_source <- "batman"
  herb_compound$evidence_source <- "BATMAN-TCM"

  data_cid <- data_batman$cid_info
  col_cid <- herbFuns$resolve_col(data_cid, c("Compound_CID", "CID", "cid"), required = TRUE)
  col_name <- herbFuns$resolve_col(data_cid, c("Name", "compound_name"), required = FALSE)
  col_smiles <- herbFuns$resolve_col(data_cid, c("SMILES", "Smiles", "smiles"), required = FALSE)
  col_inchikey <- herbFuns$resolve_col(data_cid, c("InChIKey", "inchikey"), required = FALSE)
  col_inchi <- herbFuns$resolve_col(data_cid, c("InChI", "inchi"), required = FALSE)
  col_formula <- herbFuns$resolve_col(data_cid, c("Molecular_Formula", "Molecular Formula", "Formula"), required = FALSE)
  col_weight <- herbFuns$resolve_col(data_cid, c("Molecular_Weight", "Molecular Weight", "Weight"), required = FALSE)

  compound <- tibble::tibble(
    source = "batman",
    compound_source_id = herbFuns$clean_text(data_cid[[col_cid]]),
    compound_name = herbFuns$clean_text(herbFuns$get_col(data_cid, col_name)),
    mol_id = NA_character_,
    pubchem_cid = herbFuns$clean_text(data_cid[[col_cid]]),
    smiles = herbFuns$clean_text(herbFuns$get_col(data_cid, col_smiles)),
    inchi = herbFuns$clean_text(herbFuns$get_col(data_cid, col_inchi)),
    inchikey = herbFuns$clean_text(herbFuns$get_col(data_cid, col_inchikey)),
    molecular_formula = herbFuns$clean_text(herbFuns$get_col(data_cid, col_formula)),
    molecular_weight = suppressWarnings(as.numeric(herbFuns$get_col(data_cid, col_weight)))
  )

  compound <- dplyr::filter(compound, compound_source_id %in% herb_compound$compound_source_id)
  compound <- dplyr::left_join(
    compound,
    dplyr::select(herb_compound, source, compound_source_id, compound_name_rel = compound_name),
    by = c("source", "compound_source_id")
  )
  compound$compound_name <- ifelse(
    is.na(compound$compound_name) | !nzchar(compound$compound_name),
    compound$compound_name_rel,
    compound$compound_name
  )
  compound$compound_name_rel <- NULL
  compound <- dplyr::distinct(compound)

  compound$compound_key <- herbFuns$make_compound_key(
    inchikey = compound$inchikey,
    smiles = compound$smiles,
    pubchem_cid = compound$pubchem_cid,
    source = compound$source,
    compound_source_id = compound$compound_source_id,
    compound_name = compound$compound_name
  )
  compound$source_priority <- 2L
  compound <- herbFuns$add_structure_status(compound)

  herb_compound <- dplyr::left_join(
    herb_compound,
    dplyr::select(compound, source, compound_source_id, compound_key),
    by = c("source", "compound_source_id")
  )

  herbFuns$log_progress(
    verbose,
    "BATMAN-TCM collected ", herbFuns$format_n(nrow(herb_compound)),
    " herb-compound relationships and ", herbFuns$format_n(nrow(compound)),
    " compound structure rows."
  )

  list(
    query = query,
    herb_map = herb_map,
    herb_compound = herb_compound,
    compound = compound,
    logs = list(
      duplicate_herb_names = duplicate_info$duplicate_herb_names,
      duplicate_decisions = duplicate_info$duplicate_decisions
    )
  )
}

# --------------------------------------------------------------------------
# Literature helpers and normalizer
# --------------------------------------------------------------------------

herbFuns$new_literature_template <- function(herbs = NULL, n = 0L)
{
  if (!is.null(herbs)) {
    herbs <- as.character(herbs)
    n <- length(herbs)
  }

  if (n < 0L) {
    stop("`n` must be non-negative.")
  }

  if (is.null(herbs)) {
    herbs <- rep(NA_character_, n)
  }

  tibble::tibble(
    herb = herbs,
    compound_name = rep(NA_character_, n),
    smiles = rep(NA_character_, n),
    inchikey = rep(NA_character_, n),
    inchi = rep(NA_character_, n),
    pubchem_cid = rep(NA_character_, n),
    molecular_formula = rep(NA_character_, n),
    molecular_weight = rep(NA_real_, n),
    pmid = rep(NA_character_, n),
    doi = rep(NA_character_, n),
    reference = rep(NA_character_, n),
    evidence_note = rep(NA_character_, n)
  )
}

herbFuns$new_literature_table <- function(herb, compound_name,
  smiles = NA_character_, inchikey = NA_character_, inchi = NA_character_,
  pubchem_cid = NA_character_, molecular_formula = NA_character_,
  molecular_weight = NA_real_, pmid = NA_character_, doi = NA_character_,
  reference = NA_character_, evidence_note = NA_character_)
{
  n <- max(length(herb), length(compound_name))

  data <- tibble::tibble(
    herb = rep(herb, length.out = n),
    compound_name = rep(compound_name, length.out = n),
    smiles = rep(smiles, length.out = n),
    inchikey = rep(inchikey, length.out = n),
    inchi = rep(inchi, length.out = n),
    pubchem_cid = rep(pubchem_cid, length.out = n),
    molecular_formula = rep(molecular_formula, length.out = n),
    molecular_weight = suppressWarnings(as.numeric(rep(molecular_weight, length.out = n))),
    pmid = rep(pmid, length.out = n),
    doi = rep(doi, length.out = n),
    reference = rep(reference, length.out = n),
    evidence_note = rep(evidence_note, length.out = n)
  )

  herbFuns$clean_literature_table(data)
}

herbFuns$read_literature_table <- function(file)
{
  herbFuns$clean_literature_table(herbFuns$read_table(file))
}

herbFuns$clean_literature_table <- function(data, col_map = NULL)
{
  data <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)

  if (!is.null(col_map)) {
    for (i in names(col_map)) {
      if (!col_map[[i]] %in% colnames(data)) {
        stop("Column specified by `col_map` was not found: ", col_map[[i]])
      }
      colnames(data)[colnames(data) == col_map[[i]]] <- i
    }
  }

  col_query <- herbFuns$resolve_col(
    data,
    c("herb", "Herb", "query_herb", "herb_name", "Herb.Name", "herb_cn_name", "Chinese.Name", "Chinese_Name", "中文名", "中药", "药物"),
    required = TRUE
  )
  col_name <- herbFuns$resolve_col(
    data,
    c("compound_name", "compound", "Compound", "Name", "name", "ingredient", "Ingredient", "化合物", "成分"),
    required = TRUE
  )
  col_smiles <- herbFuns$resolve_col(data, c("smiles", "SMILES", "Smiles"), required = FALSE)
  col_inchikey <- herbFuns$resolve_col(data, c("inchikey", "InChIKey", "Inchikey", "INCHIKEY"), required = FALSE)
  col_inchi <- herbFuns$resolve_col(data, c("inchi", "InChI"), required = FALSE)
  col_cid <- herbFuns$resolve_col(data, c("pubchem_cid", "PubChem_CID", "cid", "CID", "Compound_CID", "PubChem CID"), required = FALSE)
  col_formula <- herbFuns$resolve_col(data, c("molecular_formula", "Molecular_Formula", "Molecular Formula", "formula", "Formula"), required = FALSE)
  col_weight <- herbFuns$resolve_col(data, c("molecular_weight", "Molecular_Weight", "Molecular Weight", "weight", "Weight"), required = FALSE)
  col_pmid <- herbFuns$resolve_col(data, c("pmid", "PMID", "PubMed_ID", "PubMed ID"), required = FALSE)
  col_doi <- herbFuns$resolve_col(data, c("doi", "DOI"), required = FALSE)
  col_ref <- herbFuns$resolve_col(data, c("reference", "Reference", "ref", "Ref"), required = FALSE)
  col_note <- herbFuns$resolve_col(data, c("evidence_note", "note", "Note", "comment", "Comment"), required = FALSE)

  data_out <- tibble::tibble(
    herb = herbFuns$clean_text(data[[col_query]]),
    compound_name = herbFuns$clean_text(data[[col_name]]),
    smiles = herbFuns$clean_text(herbFuns$get_col(data, col_smiles)),
    inchikey = herbFuns$clean_text(herbFuns$get_col(data, col_inchikey)),
    inchi = herbFuns$clean_text(herbFuns$get_col(data, col_inchi)),
    pubchem_cid = herbFuns$clean_text(herbFuns$get_col(data, col_cid)),
    molecular_formula = herbFuns$clean_text(herbFuns$get_col(data, col_formula)),
    molecular_weight = suppressWarnings(as.numeric(herbFuns$get_col(data, col_weight))),
    pmid = herbFuns$clean_text(herbFuns$get_col(data, col_pmid)),
    doi = herbFuns$clean_text(herbFuns$get_col(data, col_doi)),
    reference = herbFuns$clean_text(herbFuns$get_col(data, col_ref)),
    evidence_note = herbFuns$clean_text(herbFuns$get_col(data, col_note))
  )

  data_out <- dplyr::filter(data_out, !is.na(herb), !is.na(compound_name))
  dplyr::distinct(data_out)
}

herbFuns$check_literature_table <- function(data)
{
  data <- herbFuns$clean_literature_table(data)

  list(
    n_rows = nrow(data),
    n_herbs = dplyr::n_distinct(data$herb),
    n_compounds = dplyr::n_distinct(data$compound_name),
    n_with_smiles = sum(!is.na(data$smiles)),
    n_with_inchikey = sum(!is.na(data$inchikey)),
    n_with_pubchem_cid = sum(!is.na(data$pubchem_cid)),
    table = data
  )
}

herbFuns$normalize_literature <- function(data_literature, herbs)
{
  query <- herbFuns$prepare_query(herbs)
  data_literature <- herbFuns$clean_literature_table(data_literature)

  data_literature$query_key <- herbFuns$match_key(data_literature$herb)
  data_hit <- dplyr::inner_join(
    data_literature,
    dplyr::select(query, query_id, query_herb, query_key),
    by = "query_key"
  )

  if (nrow(data_hit) == 0L) {
    return(herbFuns$new_empty_source_result(query))
  }

  data_hit$source <- "literature"
  data_hit$source_herb_id <- vapply(seq_len(nrow(data_hit)),
    function(i) {
      herbFuns$make_source_herb_id("literature", data_hit$query_herb[[i]])
    }, character(1L))

  herb_map <- dplyr::distinct(
    dplyr::select(
      data_hit,
      query_id,
      query_herb,
      source,
      source_herb_id
    )
  )
  herb_map$herb_cn_name <- herb_map$query_herb
  herb_map$herb_en_name <- NA_character_
  herb_map$herb_pinyin <- NA_character_
  herb_map$herb_latin_name <- NA_character_
  herb_map$match_field <- "literature_herb"
  herb_map$match_type <- "exact"
  herb_map$match_score <- 1

  data_hit$compound_source_id <- ifelse(
    !is.na(data_hit$pubchem_cid) & nzchar(data_hit$pubchem_cid),
    paste0("CID:", data_hit$pubchem_cid),
    paste0("LIT:", seq_len(nrow(data_hit)))
  )

  compound <- tibble::tibble(
    source = "literature",
    compound_source_id = data_hit$compound_source_id,
    compound_name = data_hit$compound_name,
    mol_id = NA_character_,
    pubchem_cid = data_hit$pubchem_cid,
    smiles = data_hit$smiles,
    inchi = data_hit$inchi,
    inchikey = data_hit$inchikey,
    molecular_formula = data_hit$molecular_formula,
    molecular_weight = data_hit$molecular_weight
  )
  compound$compound_key <- herbFuns$make_compound_key(
    inchikey = compound$inchikey,
    smiles = compound$smiles,
    pubchem_cid = compound$pubchem_cid,
    source = compound$source,
    compound_source_id = compound$compound_source_id,
    compound_name = compound$compound_name
  )
  compound$source_priority <- 3L
  compound <- herbFuns$add_structure_status(compound)
  compound <- dplyr::distinct(compound)

  herb_compound <- dplyr::select(
    data_hit,
    query_id,
    query_herb,
    source,
    source_herb_id,
    compound_source_id,
    compound_name,
    pmid,
    doi,
    reference,
    evidence_note
  )
  herb_compound$herb_cn_name <- herb_compound$query_herb
  herb_compound$herb_en_name <- NA_character_
  herb_compound$herb_pinyin <- NA_character_
  herb_compound$herb_latin_name <- NA_character_
  herb_compound$compound_name_source <- "literature"
  herb_compound$relationship_source <- "literature"
  herb_compound$evidence_source <- ifelse(
    !is.na(herb_compound$pmid) & nzchar(herb_compound$pmid),
    paste0("PMID:", herb_compound$pmid),
    "literature"
  )

  herb_compound <- dplyr::left_join(
    herb_compound,
    dplyr::select(compound, source, compound_source_id, compound_key),
    by = c("source", "compound_source_id")
  )

  list(
    query = query,
    herb_map = herb_map,
    herb_compound = herb_compound,
    compound = compound,
    logs = list()
  )
}

# --------------------------------------------------------------------------
# Collection assembly
# --------------------------------------------------------------------------

herbFuns$collect_herbs <- function(herbs, sources = c("tcmsp", "batman"),
  data_literature = NULL, data_tcmsp = NULL, data_batman = NULL,
  read_local = TRUE,
  batman_duplicate_policy = c("merge_if_consistent", "interactive", "keep_all", "error"),
  verbose = TRUE)
{
  batman_duplicate_policy <- match.arg(batman_duplicate_policy)
  sources <- unique(sources)
  query <- herbFuns$prepare_query(herbs)

  herbFuns$log_progress(
    verbose,
    "Start collecting herb compounds from sources: ", paste(sources, collapse = ", "),
    "; query herbs: ", paste(query$query_herb, collapse = ", "), "."
  )

  lst_piece <- list()

  if ("tcmsp" %in% sources) {
    if (is.null(data_tcmsp)) {
      if (!isTRUE(read_local)) {
        stop("`data_tcmsp` must be provided when `read_local = FALSE`.")
      }
      data_tcmsp <- herbFuns$read_tcmsp(verbose = verbose)
    }
    lst_piece$tcmsp <- herbFuns$normalize_tcmsp(data_tcmsp, herbs, verbose = verbose)
  }

  if ("batman" %in% sources) {
    if (is.null(data_batman)) {
      if (!isTRUE(read_local)) {
        stop("`data_batman` must be provided when `read_local = FALSE`.")
      }
      data_batman <- herbFuns$read_batman(verbose = verbose)
    }
    lst_piece$batman <- herbFuns$normalize_batman(
      data_batman,
      herbs,
      duplicate_policy = batman_duplicate_policy,
      verbose = verbose
    )
  }

  if (!is.null(data_literature)) {
    lst_piece$literature <- herbFuns$normalize_literature(data_literature, herbs)
  }

  collection <- herbFuns$merge_source_results(query, lst_piece)

  herbFuns$log_progress(
    verbose,
    "Finished collection: ", herbFuns$format_n(nrow(collection$herb_map)),
    " herb-map rows; ", herbFuns$format_n(nrow(collection$herb_compound)),
    " herb-compound relationships; ", herbFuns$format_n(nrow(collection$compound_unique)),
    " unique compounds."
  )

  collection
}

herbFuns$merge_source_results <- function(query, lst_piece)
{
  if (length(lst_piece) == 0L) {
    out <- herbFuns$new_empty_source_result(query)
    out$compound_unique <- tibble::tibble()
    out$source_tables <- lst_piece
    out$logs <- herbFuns$make_logs(out)
    class(out) <- "herbs_collection"
    return(out)
  }

  herb_map <- dplyr::bind_rows(lapply(lst_piece, function(x) x$herb_map))
  herb_compound <- dplyr::bind_rows(lapply(lst_piece, function(x) x$herb_compound))
  compound <- dplyr::bind_rows(lapply(lst_piece, function(x) x$compound))

  herb_compound <- herbFuns$ensure_relationship_columns(herb_compound)
  compound <- herbFuns$ensure_compound_columns(compound)

  herb_map <- dplyr::distinct(herb_map)
  herb_compound <- dplyr::distinct(herb_compound)
  compound <- dplyr::distinct(compound)
  compound <- herbFuns$add_structure_status(compound)

  compound_unique <- herbFuns$deduplicate_compounds(compound)
  compound_unique <- herbFuns$add_structure_status(compound_unique)

  out <- list(
    query = query,
    herb_map = herb_map,
    herb_compound = herb_compound,
    compound = compound,
    compound_unique = compound_unique,
    source_tables = lst_piece
  )
  out$logs <- herbFuns$make_logs(out)
  class(out) <- "herbs_collection"
  out
}

herbFuns$deduplicate_compounds <- function(compound)
{
  if (nrow(compound) == 0L) {
    return(compound)
  }

  compound <- herbFuns$ensure_columns(
    compound,
    c("compound_name_cleaned", "pubchem_name_query")
  )
  compound <- dplyr::arrange(compound, compound_key, source_priority)
  data_group <- dplyr::group_by(compound, compound_key)
  data_unique <- dplyr::summarise(
    data_group,
    compound_name = herbFuns$first_non_missing(compound_name),
    compound_name_cleaned = herbFuns$first_non_missing(compound_name_cleaned),
    pubchem_name_query = herbFuns$first_non_missing(pubchem_name_query),
    smiles = herbFuns$first_non_missing(smiles),
    inchikey = herbFuns$first_non_missing(inchikey),
    inchi = herbFuns$first_non_missing(inchi),
    pubchem_cid = herbFuns$first_non_missing(pubchem_cid),
    molecular_formula = herbFuns$first_non_missing(molecular_formula),
    molecular_weight = suppressWarnings(as.numeric(herbFuns$first_non_missing(molecular_weight))),
    source_list = herbFuns$collapse_unique(source),
    source_compound_ids = herbFuns$collapse_unique(compound_source_id),
    source_count = dplyr::n_distinct(source),
    .groups = "drop"
  )

  data_unique
}

herbFuns$make_logs <- function(collection)
{
  query <- collection$query
  herb_map <- collection$herb_map
  herb_compound <- collection$herb_compound
  compound <- collection$compound
  compound_unique <- collection$compound_unique

  matched <- unique(herb_map$query_herb)
  unmatched_herbs <- dplyr::filter(query, !query_herb %in% matched)

  compounds_without_smiles <- dplyr::filter(compound_unique, !has_smiles)
  compounds_without_inchikey <- dplyr::filter(compound_unique, !has_inchikey)

  if (nrow(herb_compound) > 0L) {
    source_summary <- dplyr::group_by(herb_compound, source)
    source_summary <- dplyr::summarise(
      source_summary,
      n_relationship = dplyr::n(),
      n_compound = dplyr::n_distinct(compound_key),
      .groups = "drop"
    )
  } else {
    source_summary <- tibble::tibble(
      source = character(0L),
      n_relationship = integer(0L),
      n_compound = integer(0L)
    )
  }

  source_logs <- lapply(collection$source_tables,
    function(z) {
      if (is.null(z$logs)) {
        return(list())
      }
      z$logs
    })

  duplicate_herb_names <- dplyr::bind_rows(lapply(source_logs,
    function(z) {
      if (is.null(z$duplicate_herb_names)) {
        return(tibble::tibble())
      }
      z$duplicate_herb_names
    }))

  duplicate_decisions <- dplyr::bind_rows(lapply(source_logs,
    function(z) {
      if (is.null(z$duplicate_decisions)) {
        return(tibble::tibble())
      }
      z$duplicate_decisions
    }))

  list(
    unmatched_herbs = unmatched_herbs,
    compounds_without_smiles = compounds_without_smiles,
    compounds_without_inchikey = compounds_without_inchikey,
    duplicate_herb_names = duplicate_herb_names,
    duplicate_decisions = duplicate_decisions,
    source_summary = source_summary,
    source_logs = source_logs,
    n_query = nrow(query),
    n_matched_query = length(matched),
    n_relationship = nrow(herb_compound),
    n_compound_source_rows = nrow(compound),
    n_compound_unique = nrow(compound_unique)
  )
}

herbFuns$inject_literature_table <- function(x, data_literature)
{
  herbFuns$validate_collection(x, stop_if_invalid = TRUE)

  herbs <- x$query$query_herb
  piece_lit <- herbFuns$normalize_literature(data_literature, herbs)
  lst_piece <- x$source_tables
  lst_piece$literature <- piece_lit

  herbFuns$merge_source_results(x$query, lst_piece)
}

# --------------------------------------------------------------------------
# Collection diagnostics and accessors
# --------------------------------------------------------------------------

herbFuns$validate_collection <- function(x, stop_if_invalid = TRUE)
{
  ok <- inherits(x, "herbs_collection") &&
    all(c("query", "herb_map", "herb_compound", "compound", "compound_unique", "logs") %in% names(x))

  if (!ok && isTRUE(stop_if_invalid)) {
    stop("Input is not a valid `herbs_collection` object.")
  }

  ok
}

herbFuns$as_compound_table <- function(x, unique = TRUE)
{
  herbFuns$validate_collection(x, stop_if_invalid = TRUE)

  if (isTRUE(unique)) {
    return(x$compound_unique)
  }

  x$compound
}

herbFuns$as_herb_compound_table <- function(x)
{
  herbFuns$validate_collection(x, stop_if_invalid = TRUE)
  x$herb_compound
}

herbFuns$get_smiles <- function(x, unique = TRUE)
{
  data <- herbFuns$as_compound_table(x, unique = unique)
  data <- dplyr::filter(data, has_smiles)
  data[, c("compound_key", "compound_name", "smiles"), drop = FALSE]
}

herbFuns$get_unmatched_herbs <- function(x)
{
  herbFuns$validate_collection(x, stop_if_invalid = TRUE)
  x$logs$unmatched_herbs
}

herbFuns$get_duplicate_herbs <- function(x)
{
  herbFuns$validate_collection(x, stop_if_invalid = TRUE)

  if (is.null(x$logs$duplicate_decisions)) {
    return(tibble::tibble())
  }

  x$logs$duplicate_decisions
}

herbFuns$stat_collection <- function(x)
{
  herbFuns$validate_collection(x, stop_if_invalid = TRUE)

  data.frame(
    n_query = x$logs$n_query,
    n_matched_query = x$logs$n_matched_query,
    n_relationship = x$logs$n_relationship,
    n_compound_source_rows = x$logs$n_compound_source_rows,
    n_compound_unique = x$logs$n_compound_unique,
    n_without_smiles = nrow(x$logs$compounds_without_smiles),
    n_without_inchikey = nrow(x$logs$compounds_without_inchikey),
    stringsAsFactors = FALSE
  )
}

herbFuns$diagnose_collection_counts <- function(x)
{
  herbFuns$validate_collection(x, stop_if_invalid = TRUE)

  data_hc <- x$herb_compound
  data_comp <- x$compound
  data_unique <- x$compound_unique

  data_overall <- tibble::tibble(
    n_herb_compound_relationships = nrow(data_hc),
    n_source_compound_rows = nrow(data_comp),
    n_unique_compounds = nrow(data_unique),
    n_relationships_minus_source_compounds = nrow(data_hc) - nrow(data_comp),
    n_source_compounds_minus_unique_compounds = nrow(data_comp) - nrow(data_unique)
  )

  if (nrow(data_hc) > 0L) {
    data_source <- dplyr::group_by(data_hc, source)
    data_source <- dplyr::summarize(
      data_source,
      n_relationships = dplyr::n(),
      n_relationship_compounds = dplyr::n_distinct(compound_key),
      .groups = "drop"
    )
  } else {
    data_source <- tibble::tibble()
  }

  if (nrow(data_comp) > 0L) {
    data_comp_source <- dplyr::group_by(data_comp, source)
    data_comp_source <- dplyr::summarize(
      data_comp_source,
      n_source_compound_rows = dplyr::n(),
      n_source_compounds = dplyr::n_distinct(compound_key),
      .groups = "drop"
    )
    data_source <- dplyr::left_join(data_source, data_comp_source, by = "source")
  }

  list(
    overall = data_overall,
    by_source = data_source
  )
}


herbFuns$count_manual_cid_input <- function(manual_cid)
{
  if (is.null(manual_cid)) {
    return(0L)
  }

  if (is.vector(manual_cid) && !is.null(names(manual_cid))) {
    data_manual <- data.frame(
      name_original = names(manual_cid),
      cid = as.character(unname(manual_cid)),
      stringsAsFactors = FALSE
    )
  } else {
    data_manual <- as.data.frame(manual_cid, stringsAsFactors = FALSE)
    if (!"cid" %in% colnames(data_manual)) {
      return(0L)
    }
    if (!"name_original" %in% colnames(data_manual)) {
      if ("compound_name" %in% colnames(data_manual)) {
        data_manual$name_original <- data_manual$compound_name
      } else {
        data_manual$name_original <- seq_len(nrow(data_manual))
      }
    }
  }

  data_manual$name_original <- herbFuns$clean_text(data_manual$name_original)
  data_manual$cid <- herbFuns$clean_text(data_manual$cid)
  data_manual <- data_manual[!is.na(data_manual$name_original) & !is.na(data_manual$cid), ]
  nrow(unique(data_manual[, c("name_original", "cid")]))
}

herbFuns$get_structure_status_label <- function(has_smiles, has_pubchem_cid, has_inchikey)
{
  ifelse(
    has_smiles,
    "SMILES available",
    ifelse(
      has_pubchem_cid | has_inchikey,
      "Identifier available",
      "Name record only"
    )
  )
}

herbFuns$make_final_compound_catalog <- function(x)
{
  herbFuns$validate_collection(x, stop_if_invalid = TRUE)

  data_hc <- herbFuns$ensure_relationship_columns(x$herb_compound)
  data_comp <- herbFuns$add_structure_status(x$compound)

  if (nrow(data_hc) == 0L) {
    return(tibble::tibble(
      Herb = character(0L),
      `Latin name` = character(0L),
      `Compound name` = character(0L),
      `Evidence sources` = character(0L),
      `PubChem CID` = character(0L),
      SMILES = character(0L),
      InChIKey = character(0L),
      InChI = character(0L),
      Formula = character(0L),
      `Molecular weight` = character(0L),
      `Structure status` = character(0L)
    ))
  }

  data_status <- herbFuns$ensure_compound_columns(data_comp)
  data_status <- dplyr::select(
    data_status,
    source,
    compound_source_id,
    compound_key,
    compound_name_structure = compound_name,
    pubchem_cid,
    smiles,
    inchi,
    inchikey,
    molecular_formula,
    molecular_weight,
    has_smiles,
    has_inchikey,
    has_pubchem_cid
  )
  data_status <- dplyr::distinct(data_status)

  data <- dplyr::left_join(
    data_hc,
    data_status,
    by = c("source", "compound_source_id", "compound_key")
  )

  data$compound_name_final <- herbFuns$fill_missing_text(
    data$compound_name,
    data$compound_name_structure
  )
  data$has_smiles <- ifelse(is.na(data$has_smiles), FALSE, data$has_smiles)
  data$has_inchikey <- ifelse(is.na(data$has_inchikey), FALSE, data$has_inchikey)
  data$has_pubchem_cid <- ifelse(is.na(data$has_pubchem_cid), FALSE, data$has_pubchem_cid)

  data$source_label <- herbFuns$format_source_label(data$source)

  data_group <- dplyr::group_by(data, query_herb, compound_key)
  data_out <- dplyr::summarise(
    data_group,
    `Latin name` = herbFuns$collapse_unique(herb_latin_name, sep = "; "),
    `Compound name` = herbFuns$first_non_missing(compound_name_final),
    `Evidence sources` = herbFuns$collapse_unique(source_label, sep = "; "),
    `PubChem CID` = herbFuns$collapse_unique(pubchem_cid, sep = "; "),
    SMILES = herbFuns$first_non_missing(smiles),
    InChIKey = herbFuns$first_non_missing(inchikey),
    InChI = herbFuns$first_non_missing(inchi),
    Formula = herbFuns$first_non_missing(molecular_formula),
    `Molecular weight` = herbFuns$first_non_missing(molecular_weight),
    has_smiles = any(has_smiles, na.rm = TRUE),
    has_inchikey = any(has_inchikey, na.rm = TRUE),
    has_pubchem_cid = any(has_pubchem_cid, na.rm = TRUE),
    .groups = "drop"
  )

  data_out$`Structure status` <- herbFuns$get_structure_status_label(
    data_out$has_smiles,
    data_out$has_pubchem_cid,
    data_out$has_inchikey
  )
  data_out$Herb <- data_out$query_herb
  data_out$`Latin name` <- ifelse(is.na(data_out$`Latin name`), "Not available", data_out$`Latin name`)

  data_out <- dplyr::select(
    data_out,
    Herb,
    `Latin name`,
    `Compound name`,
    `Evidence sources`,
    `PubChem CID`,
    SMILES,
    InChIKey,
    InChI,
    Formula,
    `Molecular weight`,
    `Structure status`
  )

  dplyr::arrange(data_out, Herb, `Compound name`)
}

herbFuns$stat_report_pubchem_final_summary <- function(x)
{
  data_catalog <- herbFuns$make_final_compound_catalog(x)

  if (nrow(data_catalog) == 0L) {
    return(tibble::tibble(
      Herb = character(0L),
      `Latin name` = character(0L),
      `Evidence sources` = character(0L),
      `Total compounds` = integer(0L),
      `SMILES coverage` = character(0L),
      `PubChem CID coverage` = character(0L),
      `Structure-ready compounds` = integer(0L)
    ))
  }

  data_catalog$has_smiles <- !is.na(data_catalog$SMILES) & nzchar(data_catalog$SMILES)
  data_catalog$has_pubchem_cid <- !is.na(data_catalog$`PubChem CID`) & nzchar(data_catalog$`PubChem CID`)

  lst_source <- lapply(seq_len(nrow(data_catalog)), function(i) {
    vec_source <- unlist(strsplit(data_catalog$`Evidence sources`[i], ";\\s*"))
    vec_source <- herbFuns$clean_text(vec_source)
    vec_source <- vec_source[!is.na(vec_source) & nzchar(vec_source)]
    if (length(vec_source) == 0L) {
      return(data.frame(Herb = data_catalog$Herb[i], source_label = NA_character_))
    }
    data.frame(
      Herb = rep(data_catalog$Herb[i], length(vec_source)),
      source_label = vec_source,
      stringsAsFactors = FALSE
    )
  })
  data_source <- dplyr::bind_rows(lst_source)
  data_source <- dplyr::group_by(data_source, Herb)
  data_source <- dplyr::summarise(
    data_source,
    `Evidence sources` = herbFuns$collapse_unique(source_label, sep = "; "),
    .groups = "drop"
  )

  data_group <- dplyr::group_by(data_catalog, Herb)
  data_out <- dplyr::summarise(
    data_group,
    `Latin name` = herbFuns$collapse_unique(`Latin name`, sep = "; "),
    `Total compounds` = dplyr::n(),
    n_with_smiles = sum(has_smiles, na.rm = TRUE),
    n_with_pubchem_cid = sum(has_pubchem_cid, na.rm = TRUE),
    .groups = "drop"
  )
  data_out <- dplyr::left_join(data_out, data_source, by = "Herb")

  data_out$`SMILES coverage` <- paste0(data_out$n_with_smiles, "/", data_out$`Total compounds`)
  data_out$`PubChem CID coverage` <- paste0(data_out$n_with_pubchem_cid, "/", data_out$`Total compounds`)
  data_out$`Structure-ready compounds` <- data_out$n_with_smiles

  data_out <- dplyr::select(
    data_out,
    Herb,
    `Latin name`,
    `Evidence sources`,
    `Total compounds`,
    `SMILES coverage`,
    `PubChem CID coverage`,
    `Structure-ready compounds`
  )

  dplyr::arrange(data_out, Herb)
}


# --------------------------------------------------------------------------
# Step3 target annotation helpers
# --------------------------------------------------------------------------

herbFuns$prepare_current_compound_target_index <- function(x)
{
  herbFuns$validate_collection(x, stop_if_invalid = TRUE)

  data_hc <- herbFuns$ensure_relationship_columns(x$herb_compound)
  data_comp <- herbFuns$add_structure_status(x$compound)

  if (nrow(data_hc) == 0L || nrow(data_comp) == 0L) {
    return(tibble::tibble())
  }

  data_status <- herbFuns$ensure_compound_columns(data_comp)
  data_status <- dplyr::select(
    data_status,
    source,
    compound_source_id,
    compound_key,
    compound_name_structure = compound_name,
    pubchem_cid,
    smiles,
    inchikey,
    inchi,
    molecular_formula,
    molecular_weight
  )
  data_status <- dplyr::distinct(data_status)

  data <- dplyr::left_join(
    data_hc,
    data_status,
    by = c("source", "compound_source_id", "compound_key")
  )

  data$compound_name_final <- herbFuns$fill_missing_text(
    data$compound_name,
    data$compound_name_structure
  )
  data$source_label <- herbFuns$format_source_label(data$source)

  data <- dplyr::select(
    data,
    query_id,
    query_herb,
    herb_cn_name,
    herb_latin_name,
    source,
    source_label,
    source_herb_id,
    compound_source_id,
    compound_key,
    compound_name = compound_name_final,
    pubchem_cid,
    smiles,
    inchikey,
    inchi,
    molecular_formula,
    molecular_weight
  )

  dplyr::distinct(data)
}

herbFuns$summarise_compound_target_index <- function(data_index)
{
  if (nrow(data_index) == 0L) {
    return(tibble::tibble())
  }

  data_group <- dplyr::group_by(data_index, query_herb, compound_key)
  data_out <- dplyr::summarise(
    data_group,
    herb_latin_name = herbFuns$collapse_unique(herb_latin_name, sep = "; "),
    compound_name = herbFuns$first_non_missing(compound_name),
    compound_evidence_sources = herbFuns$collapse_unique(source_label, sep = "; "),
    compound_source_ids = herbFuns$collapse_unique(
      paste0(source, ":", compound_source_id),
      sep = "; "
    ),
    pubchem_cid = herbFuns$collapse_unique(pubchem_cid, sep = "; "),
    smiles = herbFuns$first_non_missing(smiles),
    inchikey = herbFuns$first_non_missing(inchikey),
    .groups = "drop"
  )

  data_out
}

herbFuns$empty_compound_target_table <- function()
{
  tibble::tibble(
    query_herb = character(0L),
    herb_latin_name = character(0L),
    compound_key = character(0L),
    compound_name = character(0L),
    compound_evidence_sources = character(0L),
    compound_source_ids = character(0L),
    pubchem_cid = character(0L),
    smiles = character(0L),
    inchikey = character(0L),
    target_gene = character(0L),
    target_name = character(0L),
    target_id = character(0L),
    drugbank_id = character(0L),
    target_source = character(0L),
    target_evidence_type = character(0L),
    prediction_score = numeric(0L),
    svm_score = numeric(0L),
    rf_score = numeric(0L),
    validated = character(0L)
  )
}

herbFuns$collect_tcmsp_targets <- function(data_index, data_tcmsp_targets = NULL,
  verbose = TRUE)
{
  if (is.null(data_tcmsp_targets) || nrow(data_index) == 0L) {
    return(herbFuns$empty_compound_target_table())
  }

  data_idx <- data_index[data_index$source == "tcmsp", , drop = FALSE]
  data_idx <- data_idx[!is.na(data_idx$compound_source_id) & nzchar(data_idx$compound_source_id), , drop = FALSE]

  if (nrow(data_idx) == 0L) {
    return(herbFuns$empty_compound_target_table())
  }

  data_target <- as.data.frame(data_tcmsp_targets, stringsAsFactors = FALSE, check.names = FALSE)
  col_cn <- herbFuns$resolve_col(data_target, c("herb_cn_name", "Chinese.Name"), required = FALSE)
  col_mol <- herbFuns$resolve_col(data_target, c("MOL_ID", "Mol ID", "molecule_ID"), required = TRUE)
  col_name <- herbFuns$resolve_col(data_target, c("molecule_name", "compound_name", "Name"), required = FALSE)
  col_target <- herbFuns$resolve_col(data_target, c("target_name", "Target name"), required = FALSE)
  col_gene <- herbFuns$resolve_col(data_target, c("gene_name", "Gene", "gene"), required = TRUE)
  col_target_id <- herbFuns$resolve_col(data_target, c("target_ID", "target_id", "Target ID"), required = FALSE)
  col_drugbank <- herbFuns$resolve_col(data_target, c("drugbank_ID", "drugbank_id", "DrugBank ID"), required = FALSE)
  col_validated <- herbFuns$resolve_col(data_target, c("validated", "Validated"), required = FALSE)
  col_svm <- herbFuns$resolve_col(data_target, c("SVM_score", "svm_score"), required = FALSE)
  col_rf <- herbFuns$resolve_col(data_target, c("RF_score", "rf_score"), required = FALSE)

  data_target_std <- tibble::tibble(
    herb_match_key = herbFuns$match_key(herbFuns$get_col(data_target, col_cn)),
    compound_source_id = herbFuns$clean_text(data_target[[col_mol]]),
    compound_name_target = herbFuns$clean_text(herbFuns$get_col(data_target, col_name)),
    target_name = herbFuns$clean_text(herbFuns$get_col(data_target, col_target)),
    target_gene = herbFuns$clean_text(data_target[[col_gene]]),
    target_id = herbFuns$clean_text(herbFuns$get_col(data_target, col_target_id)),
    drugbank_id = herbFuns$clean_text(herbFuns$get_col(data_target, col_drugbank)),
    validated = herbFuns$clean_text(herbFuns$get_col(data_target, col_validated)),
    svm_score = suppressWarnings(as.numeric(herbFuns$get_col(data_target, col_svm))),
    rf_score = suppressWarnings(as.numeric(herbFuns$get_col(data_target, col_rf)))
  )

  data_target_std <- dplyr::filter(
    data_target_std,
    !is.na(compound_source_id),
    nzchar(compound_source_id),
    compound_source_id %in% data_idx$compound_source_id
  )

  if (nrow(data_target_std) == 0L) {
    return(herbFuns$empty_compound_target_table())
  }

  data_idx$herb_match_key <- herbFuns$match_key(data_idx$herb_cn_name)
  data_idx_strict <- dplyr::select(
    data_idx,
    query_herb,
    compound_key,
    compound_source_id,
    herb_match_key
  )
  data_link_key <- dplyr::inner_join(
    data_idx_strict,
    data_target_std,
    by = c("compound_source_id", "herb_match_key")
  )

  if (nrow(data_link_key) == 0L) {
    data_idx_loose <- dplyr::select(
      data_idx,
      query_herb,
      compound_key,
      compound_source_id
    )
    data_link_key <- dplyr::inner_join(
      data_idx_loose,
      data_target_std,
      by = "compound_source_id"
    )
  }

  if (nrow(data_link_key) == 0L) {
    return(herbFuns$empty_compound_target_table())
  }

  data_compound <- herbFuns$summarise_compound_target_index(data_index)
  data_out <- dplyr::left_join(
    data_link_key,
    data_compound,
    by = c("query_herb", "compound_key")
  )

  data_out$compound_name <- herbFuns$fill_missing_text(
    data_out$compound_name,
    data_out$compound_name_target
  )
  data_out$target_source <- "TCMSP"
  data_out$target_evidence_type <- "TCMSP annotation"
  data_out$prediction_score <- NA_real_

  data_out <- dplyr::select(
    data_out,
    query_herb,
    herb_latin_name,
    compound_key,
    compound_name,
    compound_evidence_sources,
    compound_source_ids,
    pubchem_cid,
    smiles,
    inchikey,
    target_gene,
    target_name,
    target_id,
    drugbank_id,
    target_source,
    target_evidence_type,
    prediction_score,
    svm_score,
    rf_score,
    validated
  )

  herbFuns$log_progress(
    verbose,
    "TCMSP target links collected: ", herbFuns$format_n(nrow(data_out)), "."
  )

  dplyr::distinct(data_out)
}

herbFuns$split_batman_known_targets <- function(data_known)
{
  if (is.null(data_known) || nrow(data_known) == 0L) {
    return(tibble::tibble(
      pubchem_cid = character(0L),
      target_gene = character(0L)
    ))
  }

  data_known <- as.data.frame(data_known, stringsAsFactors = FALSE, check.names = FALSE)
  col_cid <- herbFuns$resolve_col(data_known, c("PubChem_CID", "PubChem CID", "cid"), required = TRUE)
  col_target <- herbFuns$resolve_col(
    data_known,
    c("known_target_proteins", "known target proteins", "Targets"),
    required = TRUE
  )

  lst <- lapply(seq_len(nrow(data_known)), function(i) {
    cid_i <- herbFuns$clean_text(data_known[[col_cid]][[i]])
    target_i <- herbFuns$clean_text(data_known[[col_target]][[i]])

    if (is.na(cid_i) || !nzchar(cid_i) || is.na(target_i) || !nzchar(target_i)) {
      return(tibble::tibble())
    }

    vec_target <- unlist(strsplit(target_i, "\\|"), use.names = FALSE)
    vec_target <- herbFuns$clean_text(vec_target)
    vec_target <- vec_target[!is.na(vec_target) & nzchar(vec_target)]

    tibble::tibble(
      pubchem_cid = rep(cid_i, length(vec_target)),
      target_gene = vec_target
    )
  })

  dplyr::distinct(dplyr::bind_rows(lst))
}

herbFuns$split_batman_predicted_targets <- function(data_predicted)
{
  if (is.null(data_predicted) || nrow(data_predicted) == 0L) {
    return(tibble::tibble(
      pubchem_cid = character(0L),
      target_id = character(0L),
      prediction_score = numeric(0L)
    ))
  }

  data_predicted <- herbFuns$parse_batman_predicted_table(data_predicted)

  lst <- lapply(seq_len(nrow(data_predicted)), function(i) {
    cid_i <- herbFuns$clean_text(data_predicted$PubChem_CID[[i]])
    target_i <- herbFuns$clean_text(data_predicted$predicted_target_proteins[[i]])

    if (is.na(cid_i) || !nzchar(cid_i) || is.na(target_i) || !nzchar(target_i)) {
      return(tibble::tibble())
    }

    vec_item <- unlist(strsplit(target_i, "\\|"), use.names = FALSE)
    vec_item <- herbFuns$clean_text(vec_item)
    vec_item <- vec_item[!is.na(vec_item) & nzchar(vec_item)]

    target_id <- sub("\\(.*$", "", vec_item)
    score_text <- ifelse(
      grepl("\\(([^()]*)\\)", vec_item),
      sub("^.*\\(([^()]*)\\).*$", "\\1", vec_item),
      NA_character_
    )

    tibble::tibble(
      pubchem_cid = rep(cid_i, length(target_id)),
      target_id = herbFuns$clean_text(target_id),
      prediction_score = suppressWarnings(as.numeric(score_text))
    )
  })

  data_out <- dplyr::bind_rows(lst)
  data_out <- dplyr::filter(data_out, !is.na(target_id), nzchar(target_id))
  dplyr::distinct(data_out)
}

herbFuns$map_entrez_to_symbol <- function(entrez_id)
{
  entrez_id <- unique(herbFuns$clean_text(entrez_id))
  entrez_id <- entrez_id[!is.na(entrez_id) & nzchar(entrez_id)]

  if (length(entrez_id) == 0L) {
    return(tibble::tibble(
      target_id = character(0L),
      target_gene = character(0L)
    ))
  }

  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("Package `AnnotationDbi` is required for Entrez ID conversion.")
  }
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("Package `org.Hs.eg.db` is required for Entrez ID conversion.")
  }

  data_map <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = entrez_id,
    columns = c("SYMBOL"),
    keytype = "ENTREZID"
  )

  data_out <- tibble::tibble(
    target_id = herbFuns$clean_text(data_map$ENTREZID),
    target_gene = herbFuns$clean_text(data_map$SYMBOL)
  )
  data_out <- dplyr::filter(data_out, !is.na(target_id), nzchar(target_id))
  dplyr::distinct(data_out)
}

herbFuns$collect_batman_targets <- function(data_index, data_batman_targets = NULL,
  verbose = TRUE)
{
  if (is.null(data_batman_targets) || nrow(data_index) == 0L) {
    return(herbFuns$empty_compound_target_table())
  }

  data_idx <- data_index[!is.na(data_index$pubchem_cid) & nzchar(data_index$pubchem_cid), , drop = FALSE]
  if (nrow(data_idx) == 0L) {
    return(herbFuns$empty_compound_target_table())
  }

  vec_cid <- unique(data_idx$pubchem_cid)
  vec_cid <- vec_cid[!is.na(vec_cid) & nzchar(vec_cid)]

  data_known <- data_batman_targets$known
  data_predicted <- data_batman_targets$predicted

  data_known <- herbFuns$split_batman_known_targets(data_known)
  data_known <- data_known[data_known$pubchem_cid %in% vec_cid, , drop = FALSE]

  data_predicted <- herbFuns$parse_batman_predicted_table(data_predicted)
  data_predicted <- data_predicted[data_predicted$PubChem_CID %in% vec_cid, , drop = FALSE]
  data_predicted <- herbFuns$split_batman_predicted_targets(data_predicted)

  data_symbol <- herbFuns$map_entrez_to_symbol(data_predicted$target_id)
  if (nrow(data_predicted) > 0L) {
    data_predicted <- dplyr::left_join(data_predicted, data_symbol, by = "target_id")
    data_predicted$target_gene <- ifelse(
      is.na(data_predicted$target_gene) | !nzchar(data_predicted$target_gene),
      paste0("ENTREZ:", data_predicted$target_id),
      data_predicted$target_gene
    )
  } else {
    data_predicted$target_gene <- character(0L)
  }

  data_idx_cid <- dplyr::select(
    data_idx,
    query_herb,
    compound_key,
    pubchem_cid
  )
  data_idx_cid <- dplyr::distinct(data_idx_cid)

  data_compound <- herbFuns$summarise_compound_target_index(data_index)
  data_compound <- dplyr::select(data_compound, -dplyr::any_of("pubchem_cid"))

  if (nrow(data_known) > 0L) {
    data_known_link <- dplyr::inner_join(data_idx_cid, data_known, by = "pubchem_cid")
    data_known_link <- dplyr::left_join(
      data_known_link,
      data_compound,
      by = c("query_herb", "compound_key")
    )
    data_known_link$target_name <- NA_character_
    data_known_link$target_id <- NA_character_
    data_known_link$drugbank_id <- NA_character_
    data_known_link$target_source <- "BATMAN-TCM"
    data_known_link$target_evidence_type <- "known target"
    data_known_link$prediction_score <- NA_real_
    data_known_link$svm_score <- NA_real_
    data_known_link$rf_score <- NA_real_
    data_known_link$validated <- NA_character_
  } else {
    data_known_link <- herbFuns$empty_compound_target_table()
  }

  if (nrow(data_predicted) > 0L) {
    data_pred_link <- dplyr::inner_join(data_idx_cid, data_predicted, by = "pubchem_cid")
    data_pred_link <- dplyr::left_join(
      data_pred_link,
      data_compound,
      by = c("query_herb", "compound_key")
    )
    data_pred_link$target_name <- NA_character_
    data_pred_link$drugbank_id <- NA_character_
    data_pred_link$target_source <- "BATMAN-TCM"
    data_pred_link$target_evidence_type <- "predicted target"
    data_pred_link$svm_score <- NA_real_
    data_pred_link$rf_score <- NA_real_
    data_pred_link$validated <- NA_character_
  } else {
    data_pred_link <- herbFuns$empty_compound_target_table()
  }

  data_out <- dplyr::bind_rows(data_known_link, data_pred_link)
  if (nrow(data_out) == 0L) {
    return(herbFuns$empty_compound_target_table())
  }

  data_out <- dplyr::select(
    data_out,
    query_herb,
    herb_latin_name,
    compound_key,
    compound_name,
    compound_evidence_sources,
    compound_source_ids,
    pubchem_cid,
    smiles,
    inchikey,
    target_gene,
    target_name,
    target_id,
    drugbank_id,
    target_source,
    target_evidence_type,
    prediction_score,
    svm_score,
    rf_score,
    validated
  )

  herbFuns$log_progress(
    verbose,
    "BATMAN-TCM target links collected: ", herbFuns$format_n(nrow(data_out)), "."
  )

  dplyr::distinct(data_out)
}

herbFuns$make_compound_target_catalog <- function(data_target)
{
  if (is.null(data_target) || nrow(data_target) == 0L) {
    return(tibble::tibble(
      Herb = character(0L),
      `Latin name` = character(0L),
      `Compound name` = character(0L),
      `Compound evidence` = character(0L),
      `PubChem CID` = character(0L),
      `Target gene` = character(0L),
      `Target name` = character(0L),
      `Target source` = character(0L),
      `Target evidence` = character(0L),
      `Prediction score` = numeric(0L)
    ))
  }

  data_out <- tibble::tibble(
    Herb = data_target$query_herb,
    `Latin name` = data_target$herb_latin_name,
    `Compound name` = data_target$compound_name,
    `Compound evidence` = data_target$compound_evidence_sources,
    `PubChem CID` = data_target$pubchem_cid,
    `Target gene` = data_target$target_gene,
    `Target name` = data_target$target_name,
    `Target source` = data_target$target_source,
    `Target evidence` = data_target$target_evidence_type,
    `Prediction score` = data_target$prediction_score
  )

  data_out$`Latin name` <- ifelse(
    is.na(data_out$`Latin name`) | !nzchar(data_out$`Latin name`),
    "Not available",
    data_out$`Latin name`
  )
  data_out$`Target name` <- ifelse(
    is.na(data_out$`Target name`) | !nzchar(data_out$`Target name`),
    "Not available",
    data_out$`Target name`
  )
  data_out$`Prediction score` <- ifelse(
    is.na(data_out$`Prediction score`),
    NA_real_,
    round(data_out$`Prediction score`, 3L)
  )

  dplyr::arrange(
    dplyr::distinct(data_out),
    Herb,
    `Compound name`,
    `Target source`,
    `Target evidence`,
    `Target gene`
  )
}

herbFuns$stat_report_compound_target_summary <- function(data_target)
{
  if (is.null(data_target) || nrow(data_target) == 0L) {
    return(tibble::tibble(
      Herb = character(0L),
      `Latin name` = character(0L),
      `Compounds with targets` = integer(0L),
      `Compound-target links` = integer(0L),
      `Unique target genes` = integer(0L),
      `Target evidence sources` = character(0L),
      `TCMSP links` = integer(0L),
      `BATMAN known links` = integer(0L),
      `BATMAN predicted links` = integer(0L)
    ))
  }

  data_work <- data_target
  data_work$target_gene_count_key <- ifelse(
    is.na(data_work$target_gene) | !nzchar(data_work$target_gene),
    data_work$target_id,
    data_work$target_gene
  )
  data_work$target_evidence_label <- paste(
    data_work$target_source,
    data_work$target_evidence_type,
    sep = " "
  )

  data_group <- dplyr::group_by(data_work, query_herb)
  data_out <- dplyr::summarise(
    data_group,
    `Latin name` = herbFuns$collapse_unique(herb_latin_name, sep = "; "),
    `Compounds with targets` = dplyr::n_distinct(compound_key),
    `Compound-target links` = dplyr::n(),
    `Unique target genes` = dplyr::n_distinct(target_gene_count_key),
    `Target evidence sources` = herbFuns$collapse_unique(target_evidence_label, sep = "; "),
    `TCMSP links` = sum(target_source == "TCMSP", na.rm = TRUE),
    `BATMAN known links` = sum(
      target_source == "BATMAN-TCM" & target_evidence_type == "known target",
      na.rm = TRUE
    ),
    `BATMAN predicted links` = sum(
      target_source == "BATMAN-TCM" & target_evidence_type == "predicted target",
      na.rm = TRUE
    ),
    .groups = "drop"
  )

  data_out$Herb <- data_out$query_herb
  data_out$`Latin name` <- ifelse(
    is.na(data_out$`Latin name`) | !nzchar(data_out$`Latin name`),
    "Not available",
    data_out$`Latin name`
  )

  data_out <- dplyr::select(
    data_out,
    Herb,
    `Latin name`,
    `Compounds with targets`,
    `Compound-target links`,
    `Unique target genes`,
    `Target evidence sources`,
    `TCMSP links`,
    `BATMAN known links`,
    `BATMAN predicted links`
  )

  dplyr::arrange(data_out, Herb)
}

herbFuns$stat_target_source_summary <- function(data_target)
{
  if (is.null(data_target) || nrow(data_target) == 0L) {
    return(tibble::tibble(
      `Target source` = character(0L),
      `Target evidence` = character(0L),
      `Compounds with targets` = integer(0L),
      `Compound-target links` = integer(0L),
      `Unique target genes` = integer(0L)
    ))
  }

  data_work <- data_target
  data_work$target_gene_count_key <- ifelse(
    is.na(data_work$target_gene) | !nzchar(data_work$target_gene),
    data_work$target_id,
    data_work$target_gene
  )

  data_group <- dplyr::group_by(data_work, target_source, target_evidence_type)
  data_out <- dplyr::summarise(
    data_group,
    `Compounds with targets` = dplyr::n_distinct(compound_key),
    `Compound-target links` = dplyr::n(),
    `Unique target genes` = dplyr::n_distinct(target_gene_count_key),
    .groups = "drop"
  )
  data_out$`Target source` <- data_out$target_source
  data_out$`Target evidence` <- data_out$target_evidence_type

  data_out <- dplyr::select(
    data_out,
    `Target source`,
    `Target evidence`,
    `Compounds with targets`,
    `Compound-target links`,
    `Unique target genes`
  )

  dplyr::arrange(data_out, `Target source`, `Target evidence`)
}

# --------------------------------------------------------------------------
# Step1 statistics helpers
# --------------------------------------------------------------------------

herbFuns$stat_report_herb_summary <- function(x)
{
  herbFuns$validate_collection(x, stop_if_invalid = TRUE)

  data_hc <- herbFuns$ensure_relationship_columns(x$herb_compound)
  data_comp <- herbFuns$add_structure_status(x$compound)

  if (nrow(data_hc) == 0L) {
    return(tibble::tibble(
      Herb = character(0L),
      `Latin name` = character(0L),
      `Evidence sources` = character(0L),
      `Total compounds` = integer(0L),
      `Compounds by source` = character(0L),
      `Structure coverage` = character(0L)
    ))
  }

  data_status <- dplyr::select(
    data_comp,
    source,
    compound_source_id,
    compound_key,
    has_smiles,
    has_inchikey,
    has_pubchem_cid
  )
  data_status <- dplyr::distinct(data_status)

  data <- dplyr::left_join(
    data_hc,
    data_status,
    by = c("source", "compound_source_id", "compound_key")
  )
  data$has_smiles <- ifelse(is.na(data$has_smiles), FALSE, data$has_smiles)
  data$has_inchikey <- ifelse(is.na(data$has_inchikey), FALSE, data$has_inchikey)
  data$has_pubchem_cid <- ifelse(is.na(data$has_pubchem_cid), FALSE, data$has_pubchem_cid)

  data_source <- dplyr::group_by(data, query_herb, source)
  data_source <- dplyr::summarize(
    data_source,
    n_compound = dplyr::n_distinct(compound_key),
    .groups = "drop"
  )

  data_source$source_item <- paste0(
    herbFuns$format_source_label(data_source$source),
    ": ",
    data_source$n_compound
  )

  data_source_text <- dplyr::group_by(data_source, query_herb)
  data_source_text <- dplyr::summarize(
    data_source_text,
    `Evidence sources` = paste(herbFuns$format_source_label(source), collapse = "; "),
    `Compounds by source` = paste(source_item, collapse = "; "),
    .groups = "drop"
  )

  data_latin <- dplyr::group_by(data, query_herb)
  data_latin <- dplyr::summarize(
    data_latin,
    `Latin name` = herbFuns$collapse_unique(herb_latin_name, sep = "; "),
    .groups = "drop"
  )

  data_overall <- dplyr::group_by(data, query_herb)
  data_overall <- dplyr::summarize(
    data_overall,
    `Total compounds` = dplyr::n_distinct(compound_key),
    n_with_smiles = dplyr::n_distinct(compound_key[has_smiles]),
    n_with_inchikey = dplyr::n_distinct(compound_key[has_inchikey]),
    n_with_pubchem_cid = dplyr::n_distinct(compound_key[has_pubchem_cid]),
    .groups = "drop"
  )

  data_overall$`Structure coverage` <- paste0(
    "SMILES ", data_overall$n_with_smiles, "/", data_overall$`Total compounds`,
    "; InChIKey ", data_overall$n_with_inchikey, "/", data_overall$`Total compounds`,
    "; PubChem CID ", data_overall$n_with_pubchem_cid, "/", data_overall$`Total compounds`
  )

  data_out <- dplyr::left_join(data_overall, data_latin, by = "query_herb")
  data_out <- dplyr::left_join(data_out, data_source_text, by = "query_herb")
  data_out$`Latin name` <- ifelse(is.na(data_out$`Latin name`), "Not available", data_out$`Latin name`)
  data_out$Herb <- data_out$query_herb

  data_out <- dplyr::select(
    data_out,
    Herb,
    `Latin name`,
    `Evidence sources`,
    `Total compounds`,
    `Compounds by source`,
    `Structure coverage`
  )

  dplyr::arrange(data_out, Herb)
}

herbFuns$stat_source_compounds <- function(x)
{
  herbFuns$validate_collection(x, stop_if_invalid = TRUE)

  data_hc <- herbFuns$ensure_relationship_columns(x$herb_compound)
  data_comp <- herbFuns$add_structure_status(x$compound)

  if (nrow(data_hc) > 0L) {
    data_rel <- dplyr::group_by(data_hc, source)
    data_rel <- dplyr::summarize(
      data_rel,
      source_label = herbFuns$format_source_label(dplyr::first(source)),
      n_query_herb = dplyr::n_distinct(query_herb),
      n_source_herb = dplyr::n_distinct(source_herb_id),
      n_relationship = dplyr::n(),
      n_relationship_compound = dplyr::n_distinct(compound_key),
      .groups = "drop"
    )
  } else {
    data_rel <- tibble::tibble(
      source = character(0L),
      source_label = character(0L),
      n_query_herb = integer(0L),
      n_source_herb = integer(0L),
      n_relationship = integer(0L),
      n_relationship_compound = integer(0L)
    )
  }

  if (nrow(data_comp) > 0L) {
    data_struct <- dplyr::group_by(data_comp, source)
    data_struct <- dplyr::summarize(
      data_struct,
      n_compound_source_rows = dplyr::n(),
      n_unique_compound = dplyr::n_distinct(compound_key),
      n_with_smiles = sum(has_smiles, na.rm = TRUE),
      n_without_smiles = sum(!has_smiles, na.rm = TRUE),
      n_with_inchikey = sum(has_inchikey, na.rm = TRUE),
      n_without_inchikey = sum(!has_inchikey, na.rm = TRUE),
      n_with_pubchem_cid = sum(has_pubchem_cid, na.rm = TRUE),
      n_without_pubchem_cid = sum(!has_pubchem_cid, na.rm = TRUE),
      n_complete_structure = sum(structure_status == "complete_structure", na.rm = TRUE),
      n_partial_structure = sum(structure_status == "partial_structure", na.rm = TRUE),
      n_cid_only = sum(structure_status == "cid_only", na.rm = TRUE),
      n_name_only = sum(structure_status == "name_only", na.rm = TRUE),
      .groups = "drop"
    )
  } else {
    data_struct <- tibble::tibble(source = character(0L))
  }

  data_out <- dplyr::full_join(data_rel, data_struct, by = "source")
  data_out$source_label <- ifelse(
    is.na(data_out$source_label),
    herbFuns$format_source_label(data_out$source),
    data_out$source_label
  )
  data_out <- dplyr::arrange(data_out, source)
  data_out
}

herbFuns$stat_herb_compounds <- function(x)
{
  herbFuns$validate_collection(x, stop_if_invalid = TRUE)

  data_hc <- herbFuns$ensure_relationship_columns(x$herb_compound)
  data_comp <- herbFuns$add_structure_status(x$compound)

  if (nrow(data_hc) == 0L) {
    return(tibble::tibble())
  }

  data_status <- dplyr::select(
    data_comp,
    source,
    compound_source_id,
    compound_key,
    has_smiles,
    has_inchikey,
    has_pubchem_cid,
    structure_status
  )
  data_status <- dplyr::distinct(data_status)

  data <- dplyr::left_join(
    data_hc,
    data_status,
    by = c("source", "compound_source_id", "compound_key")
  )
  data$has_smiles <- ifelse(is.na(data$has_smiles), FALSE, data$has_smiles)
  data$has_inchikey <- ifelse(is.na(data$has_inchikey), FALSE, data$has_inchikey)
  data$has_pubchem_cid <- ifelse(is.na(data$has_pubchem_cid), FALSE, data$has_pubchem_cid)
  data$structure_status <- ifelse(
    is.na(data$structure_status),
    "missing_structure",
    data$structure_status
  )

  data_group <- dplyr::group_by(data, query_herb, source)
  data_out <- dplyr::summarize(
    data_group,
    source_label = herbFuns$format_source_label(dplyr::first(source)),
    n_relationship = dplyr::n(),
    n_unique_compound = dplyr::n_distinct(compound_key),
    n_with_smiles = sum(has_smiles, na.rm = TRUE),
    n_without_smiles = sum(!has_smiles, na.rm = TRUE),
    n_with_inchikey = sum(has_inchikey, na.rm = TRUE),
    n_without_inchikey = sum(!has_inchikey, na.rm = TRUE),
    n_with_pubchem_cid = sum(has_pubchem_cid, na.rm = TRUE),
    n_without_pubchem_cid = sum(!has_pubchem_cid, na.rm = TRUE),
    n_complete_structure = sum(structure_status == "complete_structure", na.rm = TRUE),
    n_partial_structure = sum(structure_status == "partial_structure", na.rm = TRUE),
    n_cid_only = sum(structure_status == "cid_only", na.rm = TRUE),
    n_name_only = sum(structure_status == "name_only", na.rm = TRUE),
    .groups = "drop"
  )

  dplyr::arrange(data_out, query_herb, source)
}

herbFuns$diagnose_missing_structures <- function(x,
  missing = c("smiles", "inchikey", "pubchem_cid", "any"))
{
  missing <- match.arg(missing)
  herbFuns$validate_collection(x, stop_if_invalid = TRUE)

  data_hc <- herbFuns$ensure_relationship_columns(x$herb_compound)
  data_comp <- herbFuns$add_structure_status(x$compound)

  if (nrow(data_hc) == 0L || nrow(data_comp) == 0L) {
    return(tibble::tibble())
  }

  data_status <- dplyr::select(
    data_comp,
    source,
    compound_source_id,
    compound_key,
    smiles,
    inchikey,
    pubchem_cid,
    has_smiles,
    has_inchikey,
    has_pubchem_cid,
    structure_status
  )
  data_status <- dplyr::distinct(data_status)

  data <- dplyr::left_join(
    data_hc,
    data_status,
    by = c("source", "compound_source_id", "compound_key")
  )
  data$has_smiles <- ifelse(is.na(data$has_smiles), FALSE, data$has_smiles)
  data$has_inchikey <- ifelse(is.na(data$has_inchikey), FALSE, data$has_inchikey)
  data$has_pubchem_cid <- ifelse(is.na(data$has_pubchem_cid), FALSE, data$has_pubchem_cid)
  data$structure_status <- ifelse(
    is.na(data$structure_status),
    "missing_structure",
    data$structure_status
  )

  if (identical(missing, "smiles")) {
    data <- dplyr::filter(data, !has_smiles)
  }

  if (identical(missing, "inchikey")) {
    data <- dplyr::filter(data, !has_inchikey)
  }

  if (identical(missing, "pubchem_cid")) {
    data <- dplyr::filter(data, !has_pubchem_cid)
  }

  if (identical(missing, "any")) {
    data <- dplyr::filter(data, !has_smiles | !has_inchikey | !has_pubchem_cid)
  }

  data$missing_type <- ifelse(
    !data$has_smiles & !data$has_inchikey & !data$has_pubchem_cid,
    "name_only",
    ifelse(
      !data$has_smiles,
      "missing_smiles",
      ifelse(!data$has_inchikey, "missing_inchikey", "missing_pubchem_cid")
    )
  )

  data <- dplyr::select(
    data,
    query_herb,
    source,
    source_herb_id,
    herb_cn_name,
    compound_source_id,
    compound_name,
    compound_key,
    pubchem_cid,
    smiles,
    inchikey,
    structure_status,
    missing_type,
    evidence_source
  )
  dplyr::distinct(data)
}

herbFuns$write_collection <- function(x, path)
{
  herbFuns$validate_collection(x, stop_if_invalid = TRUE)
  saveRDS(x, path)
  invisible(path)
}

herbFuns$read_collection <- function(path)
{
  x <- readRDS(path)
  herbFuns$validate_collection(x, stop_if_invalid = TRUE)
  x
}
