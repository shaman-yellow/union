# ==========================================================================
# workflow of mimedb
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_mimedb <- setClass("job_mimedb", 
  contains = c("job"),
  prototype = prototype(
    pg = "mimedb",
    info = c("https://mimedb.org/downloads"),
    cite = "",
    method = "",
    tag = "mimedb",
    analysis = "MiMeDB 代谢物数据挖掘"
    ))

job_mimedb <- function(patterns,
  ignore.case = TRUE,
  fixed = TRUE,
  prefer_name_contains = TRUE,
  auto_fuzzy = TRUE,
  fuzzy_max_dist = .25,
  fuzzy_method = "jw",
  use_broad = TRUE,
  unique_ids = TRUE,
  print_data = TRUE,
  print_stat = TRUE,
  stop_if_empty = TRUE
)
{
  .escape_regex <- function(x)
  {
    gsub(
      "([][{}()+*^$|\\\\?.])",
      "\\\\\\1",
      x,
      perl = TRUE
    )
  }

  .normalize_mimedb_string <- function(x,
    ignore.case = TRUE,
    sort_words = FALSE
  )
  {
    x <- as.character(x)
    x[is.na(x)] <- ""

    if (isTRUE(ignore.case)) {
      x <- tolower(x)
    }

    x <- gsub("-", " ", x, fixed = TRUE)
    x <- gsub("_", " ", x, fixed = TRUE)
    x <- gsub("[[:punct:]]+", " ", x, perl = TRUE)
    x <- gsub("\\s+", " ", x, perl = TRUE)
    x <- trimws(x)

    if (isTRUE(sort_words)) {
      x <- vapply(
        strsplit(x, "\\s+"),
        function(words) {
          words <- words[words != ""]
          paste(sort(words), collapse = " ")
        },
        character(1)
      )
    }

    x
  }

  .get_pattern_use <- function(pattern, fixed = TRUE)
  {
    if (isTRUE(fixed)) {
      return(.escape_regex(pattern))
    }

    pattern
  }

  .make_index <- function(pattern,
    row_id,
    match_type,
    match_field,
    fuzzy_name = NA_character_,
    fuzzy_distance = NA_real_
  )
  {
    if (!length(row_id)) {
      return(NULL)
    }

    data.frame(
      pattern = rep(pattern, length(row_id)),
      .row_id = row_id,
      match_type = rep(match_type, length(row_id)),
      match_field = match_field,
      fuzzy_name = fuzzy_name,
      fuzzy_distance = fuzzy_distance,
      stringsAsFactors = FALSE
    )
  }

  patterns <- unique(as.character(patterns))
  patterns <- patterns[!is.na(patterns) & patterns != ""]

  if (!length(patterns)) {
    stop("No valid pattern was supplied.")
  }

  data_raw <- ftibble(get_url_data(
    "mimedb_metabolites_v2.csv",
    "https://mimedb.org/system/downloads/2.0/mimedb_metabolites_v2.csv",
    "mimedb_metabolites_v2",
    dir = .prefix("mimedb", "db"),
    fun_decompress = NULL
  ))

  vec_need <- c("mime_id", "name", "description")
  vec_miss <- setdiff(vec_need, colnames(data_raw))

  if (length(vec_miss) > 0L) {
    stop(glue::glue(
      "Missing required column(s): {paste(vec_miss, collapse = ', ')}."
    ))
  }

  vec_name <- as.character(data_raw$name)
  vec_description <- as.character(data_raw$description)

  vec_name[is.na(vec_name)] <- ""
  vec_description[is.na(vec_description)] <- ""

  vec_name_norm <- .normalize_mimedb_string(
    vec_name,
    ignore.case = ignore.case,
    sort_words = FALSE
  )

  vec_description_norm <- .normalize_mimedb_string(
    vec_description,
    ignore.case = ignore.case,
    sort_words = FALSE
  )

  vec_pattern_norm <- .normalize_mimedb_string(
    patterns,
    ignore.case = ignore.case,
    sort_words = FALSE
  )

  stat_match <- data.frame(
    pattern = patterns,
    pattern_norm = vec_pattern_norm,
    matched = FALSE,
    final_match_type = NA_character_,
    n_name_exact = 0L,
    n_name_contains = 0L,
    n_name_fuzzy = 0L,
    n_broad = 0L,
    n_total = 0L,
    fuzzy_attempted = FALSE,
    fuzzy_name = NA_character_,
    fuzzy_distance = NA_real_,
    broad_attempted = FALSE,
    stringsAsFactors = FALSE
  )

  lst_index <- list()
  i_index <- 0L
  vec_pattern_done <- rep(FALSE, length(patterns))

  lst_name_index <- split(
    seq_along(vec_name_norm),
    vec_name_norm
  )

  for (i in seq_along(patterns)) {
    row_id <- lst_name_index[[vec_pattern_norm[i]]]

    if (is.null(row_id)) {
      row_id <- integer(0)
    }

    if (length(row_id) > 0L) {
      i_index <- i_index + 1L

      lst_index[[i_index]] <- .make_index(
        pattern = patterns[i],
        row_id = row_id,
        match_type = "name_exact",
        match_field = rep("name", length(row_id))
      )

      vec_pattern_done[i] <- TRUE
    }
  }

  id_need <- which(!vec_pattern_done)

  if (isTRUE(prefer_name_contains) && length(id_need) > 0L) {
    patterns_need <- vec_pattern_norm[id_need]

    patterns_need_use <- vapply(
      patterns_need,
      .get_pattern_use,
      character(1),
      fixed = fixed
    )

    pattern_all_name <- paste0(
      "(",
      paste(patterns_need_use, collapse = ")|("),
      ")"
    )

    vec_keep_name <- grepl(
      pattern = pattern_all_name,
      x = vec_name_norm,
      ignore.case = FALSE,
      perl = TRUE
    )

    vec_keep_name[is.na(vec_keep_name)] <- FALSE

    id_name_base <- which(vec_keep_name)
    vec_name_base <- vec_name_norm[id_name_base]

    if (length(id_name_base) > 0L) {
      for (i in seq_along(id_need)) {
        pattern_use <- .get_pattern_use(
          vec_pattern_norm[id_need[i]],
          fixed = fixed
        )

        vec_hit <- grepl(
          pattern = pattern_use,
          x = vec_name_base,
          ignore.case = FALSE,
          perl = TRUE
        )

        vec_hit[is.na(vec_hit)] <- FALSE
        row_id <- id_name_base[which(vec_hit)]

        if (length(row_id) > 0L) {
          i_index <- i_index + 1L

          lst_index[[i_index]] <- .make_index(
            pattern = patterns[id_need[i]],
            row_id = row_id,
            match_type = "name_contains",
            match_field = rep("name", length(row_id))
          )

          vec_pattern_done[id_need[i]] <- TRUE
        }
      }
    }
  }

  id_need <- which(!vec_pattern_done)

  if (isTRUE(auto_fuzzy) && length(id_need) > 0L) {
    if (!requireNamespace("stringdist", quietly = TRUE)) {
      stop('Package "stringdist" is required when `auto_fuzzy = TRUE`.')
    }

    stat_match$fuzzy_attempted[id_need] <- TRUE

    vec_name_unique <- unique(vec_name_norm)
    vec_name_unique <- vec_name_unique[
      !is.na(vec_name_unique) &
        vec_name_unique != ""
    ]

    if (length(vec_name_unique) > 0L) {
      vec_pattern_fuzzy <- .normalize_mimedb_string(
        patterns[id_need],
        ignore.case = ignore.case,
        sort_words = TRUE
      )

      vec_name_fuzzy <- .normalize_mimedb_string(
        vec_name_unique,
        ignore.case = ignore.case,
        sort_words = TRUE
      )

      idx <- stringdist::amatch(
        x = vec_pattern_fuzzy,
        table = vec_name_fuzzy,
        method = fuzzy_method,
        maxDist = fuzzy_max_dist
      )

      vec_distance <- rep(NA_real_, length(vec_pattern_fuzzy))
      valid <- !is.na(idx)

      if (any(valid)) {
        vec_distance[valid] <- stringdist::stringdist(
          a = vec_pattern_fuzzy[valid],
          b = vec_name_fuzzy[idx[valid]],
          method = fuzzy_method
        )
      }

      for (i in seq_along(id_need)) {
        if (is.na(idx[i])) {
          next
        }

        id_pattern <- id_need[i]
        matched_name_norm <- vec_name_unique[idx[i]]
        row_id <- which(vec_name_norm == matched_name_norm)

        if (!length(row_id)) {
          next
        }

        fuzzy_name <- vec_name[row_id]
        fuzzy_name[is.na(fuzzy_name) | fuzzy_name == ""] <- matched_name_norm

        i_index <- i_index + 1L

        lst_index[[i_index]] <- .make_index(
          pattern = patterns[id_pattern],
          row_id = row_id,
          match_type = "name_fuzzy",
          match_field = rep("name", length(row_id)),
          fuzzy_name = fuzzy_name,
          fuzzy_distance = rep(vec_distance[i], length(row_id))
        )

        stat_match$fuzzy_name[id_pattern] <- paste(
          unique(fuzzy_name),
          collapse = "; "
        )

        stat_match$fuzzy_distance[id_pattern] <- vec_distance[i]

        vec_pattern_done[id_pattern] <- TRUE
      }
    }
  }

  id_need <- which(!vec_pattern_done)

  if (isTRUE(use_broad) && length(id_need) > 0L) {
    stat_match$broad_attempted[id_need] <- TRUE

    patterns_need <- vec_pattern_norm[id_need]

    patterns_need_use <- vapply(
      patterns_need,
      .get_pattern_use,
      character(1),
      fixed = fixed
    )

    pattern_all <- paste0(
      "(",
      paste(patterns_need_use, collapse = ")|("),
      ")"
    )

    vec_text_norm <- paste(vec_name_norm, vec_description_norm, sep = " ")

    vec_keep <- grepl(
      pattern = pattern_all,
      x = vec_text_norm,
      ignore.case = FALSE,
      perl = TRUE
    )

    vec_keep[is.na(vec_keep)] <- FALSE

    id_base <- which(vec_keep)
    vec_name_base <- vec_name_norm[id_base]
    vec_description_base <- vec_description_norm[id_base]

    if (length(id_base) > 0L) {
      for (i in seq_along(id_need)) {
        pattern_use <- .get_pattern_use(
          vec_pattern_norm[id_need[i]],
          fixed = fixed
        )

        vec_hit_name <- grepl(
          pattern = pattern_use,
          x = vec_name_base,
          ignore.case = FALSE,
          perl = TRUE
        )

        vec_hit_description <- grepl(
          pattern = pattern_use,
          x = vec_description_base,
          ignore.case = FALSE,
          perl = TRUE
        )

        vec_hit_name[is.na(vec_hit_name)] <- FALSE
        vec_hit_description[is.na(vec_hit_description)] <- FALSE

        vec_hit <- vec_hit_name | vec_hit_description
        row_id <- id_base[which(vec_hit)]

        if (!length(row_id)) {
          next
        }

        match_field <- ifelse(
          vec_hit_name[vec_hit],
          "name",
          "description"
        )

        match_field <- ifelse(
          vec_hit_name[vec_hit] & vec_hit_description[vec_hit],
          "name;description",
          match_field
        )

        i_index <- i_index + 1L

        lst_index[[i_index]] <- .make_index(
          pattern = patterns[id_need[i]],
          row_id = row_id,
          match_type = "broad",
          match_field = match_field
        )

        vec_pattern_done[id_need[i]] <- TRUE
      }
    }
  }

  lst_index <- Filter(Negate(is.null), lst_index)

  if (!length(lst_index)) {
    data_index <- data.frame(
      pattern = character(0),
      .row_id = integer(0),
      match_type = character(0),
      match_field = character(0),
      fuzzy_name = character(0),
      fuzzy_distance = numeric(0),
      stringsAsFactors = FALSE
    )
  } else {
    data_index <- do.call(rbind, lst_index)
    rownames(data_index) <- NULL
  }

  if (nrow(data_index) > 0L) {
    tab_type <- table(data_index$pattern, data_index$match_type)

    for (pattern in rownames(tab_type)) {
      id_pattern <- match(pattern, stat_match$pattern)

      if (is.na(id_pattern)) {
        next
      }

      if ("name_exact" %in% colnames(tab_type)) {
        stat_match$n_name_exact[id_pattern] <- tab_type[
          pattern,
          "name_exact"
        ]
      }

      if ("name_contains" %in% colnames(tab_type)) {
        stat_match$n_name_contains[id_pattern] <- tab_type[
          pattern,
          "name_contains"
        ]
      }

      if ("name_fuzzy" %in% colnames(tab_type)) {
        stat_match$n_name_fuzzy[id_pattern] <- tab_type[
          pattern,
          "name_fuzzy"
        ]
      }

      if ("broad" %in% colnames(tab_type)) {
        stat_match$n_broad[id_pattern] <- tab_type[
          pattern,
          "broad"
        ]
      }
    }
  }

  stat_match$n_total <- stat_match$n_name_exact +
    stat_match$n_name_contains +
    stat_match$n_name_fuzzy +
    stat_match$n_broad

  stat_match$matched <- stat_match$n_total > 0L

  stat_match$final_match_type <- ifelse(
    stat_match$n_name_exact > 0L,
    "name_exact",
    ifelse(
      stat_match$n_name_contains > 0L,
      "name_contains",
      ifelse(
        stat_match$n_name_fuzzy > 0L,
        "name_fuzzy",
        ifelse(
          stat_match$n_broad > 0L,
          "broad",
          "unmatched"
        )
      )
    )
  )

  stat_match <- tibble::as_tibble(stat_match)

  if (nrow(data_index) == 0L) {
    data <- data_raw[0L, , drop = FALSE]

    data <- cbind(
      data.frame(
        pattern = character(0),
        match_type = character(0),
        match_field = character(0),
        fuzzy_name = character(0),
        fuzzy_distance = numeric(0),
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),
      data
    )

    data <- tibble::as_tibble(data)

    x <- .job_mimedb()
    x$ids <- stats::setNames(character(0), character(0))
    x$data <- data
    x$stat_match <- stat_match

    if (isTRUE(print_stat)) {
      print(stat_match)
    }

    if (isTRUE(stop_if_empty)) {
      stop("!nrow(data).")
    }

    return(x)
  }

  data <- data_raw[data_index$.row_id, , drop = FALSE]

  data <- cbind(
    data.frame(
      pattern = data_index$pattern,
      match_type = data_index$match_type,
      match_field = data_index$match_field,
      fuzzy_name = data_index$fuzzy_name,
      fuzzy_distance = data_index$fuzzy_distance,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    data
  )

  rownames(data) <- NULL
  data <- tibble::as_tibble(data)

  message(glue::glue("Got data ({nrow(data)})."))

  if (isTRUE(print_stat)) {
    message("Pattern matching summary:")
    print(stat_match)
  }

  if (isTRUE(print_data)) {
    print(data)
  }

  x <- .job_mimedb()

  if (isTRUE(unique_ids)) {
    data_ids <- data[!duplicated(data$mime_id), , drop = FALSE]
    x$ids <- stats::setNames(data_ids$mime_id, data_ids$name)
  } else {
    x$ids <- stats::setNames(data$mime_id, data$name)
  }

  x$data <- data
  x$stat_match <- stat_match

  x$dir_data <- create_job_cache_dir(x, class(x))
  x$file_data <- file.path(x$dir_data, "data_matches.csv")
  data.table::fwrite(x$data, x$file_data)
  gett_file(x$file_data)

  message(glue::glue("The file {x$file_data} is ready in clipboard."))
  x$patterns <- patterns

  return(x)
}

setMethod("step0", signature = c(x = "job_mimedb"),
  function(x){
    step_message("Prepare your data with function `job_mimedb`.")
  })

setMethod("step1", signature = c(x = "job_mimedb"),
  function(x, index = NULL)
  {
    step_message("Download data.")
    if (!is.null(index)) {
      mimeFuns$mimedb_open_urls(x$data$mime_id[index])
    }
    return(x)
  })

setMethod("step2", signature = c(x = "job_mimedb"),
  function(x,
    index = NULL,
    dir = "~/Downloads",
    add_detail_tables = TRUE
  )
  {
    step_message("Collate MiMeDB proteins and microbial sources.")

    if (is.null(x$data) || nrow(x$data) == 0L) {
      warning("No MiMeDB matched data was found in x$data.")
      return(x)
    }

    data_match <- as.data.frame(x$data, stringsAsFactors = FALSE)

    if (is.null(index)) {
      index <- seq_len(nrow(data_match))
    }

    index <- unique(as.integer(index))
    index <- index[!is.na(index)]
    index <- index[index >= 1L & index <= nrow(data_match)]

    if (!length(index)) {
      warning("No valid row index was supplied.")
      return(x)
    }

    data_selected <- data_match[index, , drop = FALSE]
    data_selected$.selected_index <- index

    if (!"pattern" %in% colnames(data_selected)) {
      data_selected$pattern <- data_selected$name
    }

    if (!"match_type" %in% colnames(data_selected)) {
      data_selected$match_type <- NA_character_
    }

    patterns_input <- if (!is.null(x$patterns)) {
      unique(as.character(x$patterns))
    } else {
      unique(as.character(data_match$pattern))
    }

    patterns_input <- patterns_input[
      !is.na(patterns_input) &
        patterns_input != ""
    ]

    files <- mimeFuns$get_html_files_by_id(
      dir = dir,
      id_pattern = "MMDBc[0-9]+"
    )

    data_file <- mimeFuns$prepare_file_status(
      ids = data_selected$mime_id,
      index = data_selected$.selected_index,
      files = files
    )

    if (any(!data_file$file_found)) {
      warning(glue::glue(
        "Missing local HTML files for {sum(!data_file$file_found)} selected metabolite(s)."
      ))
    }

    lst_tables <- lapply(
      seq_len(nrow(data_selected)),
      function(i) {
        mimeFuns$prepare_mimedb_detail_tables(
          file = data_file$file[i]
        )
      }
    )

    data_selected$.item_label <- ifelse(
      data_selected$pattern == data_selected$name,
      data_selected$pattern,
      paste0(data_selected$pattern, " -> ", data_selected$name)
    )

    names(lst_tables) <- make.unique(
      formal_name(data_selected$.item_label)
    )

    data_stat_metabolite <- mimeFuns$summarize_mimedb_detail_tables(
      data_items = data_selected,
      lst_tables = lst_tables,
      data_file = data_file
    )

    data_stat_pattern <- mimeFuns$summarize_mimedb_pattern_selection(
      patterns = patterns_input,
      data_candidates = data_match,
      data_selected = data_selected,
      data_item_stat = data_stat_metabolite,
      stat_match = x$stat_match
    )

    x$selected_index <- index
    x$selected_data <- tibble::as_tibble(data_selected)
    x$stat_selected_pattern <- data_stat_pattern
    x$stat_selected_metabolite <- data_stat_metabolite
    x$stat_selected_file <- tibble::as_tibble(data_file)
    x$lst_mimedb_tables <- lst_tables

    if (is.null(x$lst_refine)) {
      x$lst_refine <- list()
    }

    x$lst_refine$mimedb_detail <- list(
      selected_data = x$selected_data,
      stat_pattern = data_stat_pattern,
      stat_metabolite = data_stat_metabolite,
      stat_file = x$stat_selected_file,
      tables = lst_tables
    )

    metabolite <- as.character(data_selected$.item_label)
    metabolite[is.na(metabolite) | metabolite == ""] <- as.character(
      data_selected$name[is.na(metabolite) | metabolite == ""]
    )

    n_protein_entry <- data_stat_metabolite$n_protein_entry
    n_unique_protein <- data_stat_metabolite$n_unique_protein
    n_unique_uniprot <- data_stat_metabolite$n_unique_uniprot
    n_microbe_entry <- data_stat_metabolite$n_microbe_entry
    n_phylum <- data_stat_metabolite$n_phylum
    n_total_species <- data_stat_metabolite$n_total_species
    n_total_species <- ifelse(
      is.na(n_total_species),
      "NA",
      as.character(n_total_species)
    )

    ts.proteins <- lapply(
      lst_tables,
      function(data_item) {
        data_item$Human_Proteins_and_Enzymes
      }
    )

    names(ts.proteins) <- names(lst_tables)

    lab_proteins <- as.character(glue::glue(
      "{x@sig} {metabolite} MiMeDB human proteins and enzymes"
    ))

    labs_proteins <- as.character(glue::glue(
      "MiMeDB {metabolite} 人类相关蛋白与酶表|||",
      "该表展示经人工确认后的 MiMeDB 代谢物条目 {metabolite} 对应的人类相关蛋白与酶信息，",
      "共获取 {n_protein_entry} 条蛋白/酶记录，涉及 {n_unique_protein} 个唯一蛋白名称和 {n_unique_uniprot} 个 Uniprot ID，",
      "表中包括蛋白类型、蛋白名称及 Uniprot ID。"
    ))

    names(lab_proteins) <- names(ts.proteins)
    names(labs_proteins) <- names(ts.proteins)

    ts.proteins <- set_lab_legend(
      ts.proteins,
      lab_proteins,
      labs_proteins
    )

    ts.microbe <- lapply(
      lst_tables,
      function(data_item) {
        data_item$Microbial_Sources
      }
    )

    names(ts.microbe) <- names(lst_tables)

    lab_microbe <- as.character(glue::glue(
      "{x@sig} {metabolite} MiMeDB microbial sources"
    ))

    labs_microbe <- as.character(glue::glue(
      "MiMeDB {metabolite} 微生物来源表|||",
      "该表展示经人工确认后的 MiMeDB 代谢物条目 {metabolite} 对应的微生物来源信息，",
      "共获取 {n_microbe_entry} 条微生物来源记录，涉及 {n_phylum} 个门水平分类，相关物种数量合计为 {n_total_species}，",
      "表中包括微生物界/超界、门水平分类、物种数量及宿主/体内位点信息。"
    ))

    names(lab_microbe) <- names(ts.microbe)
    names(labs_microbe) <- names(ts.microbe)

    ts.microbe <- set_lab_legend(
      ts.microbe,
      lab_microbe,
      labs_microbe
    )

    t.pattern <- set_lab_legend(
      data_stat_pattern,
      glue::glue("{x@sig} MiMeDB input-level retrieval summary"),
      glue::glue(
        "MiMeDB 输入代谢物匹配与检索统计表|||",
        "该表以最初输入的代谢物名称为单位，汇总 MiMeDB 初步候选匹配数量、",
        "人工确认后的真实条目数量、MiMeDB 获取情况，以及对应获取到的人类蛋白/酶和微生物来源记录数。"
      )
    )

    t.metabolite <- set_lab_legend(
      data_stat_metabolite,
      glue::glue("{x@sig} MiMeDB metabolite-level retrieval summary"),
      glue::glue(
        "MiMeDB 代谢物条目层面检索统计表|||",
        "该表以人工确认后的 MiMeDB 代谢物条目为单位，统计每个代谢物对应的人类蛋白/酶记录数、",
        "唯一蛋白数、Uniprot ID 数、微生物来源记录数、微生物分类数量及物种数量。"
      )
    )

    x <- tablesAdd(
      x,
      t.pattern = t.pattern,
      t.metabolite = t.metabolite
    )

    p.stat <- mimeFuns$plot_mimedb_record_summary(
      data_stat = data_stat_metabolite
    )

    p.stat <- set_lab_legend(
      p.stat,
      glue::glue("{x@sig} MiMeDB retrieval statistics"),
      glue::glue(
        "MiMeDB 代谢物关联记录统计图|||",
        "该图展示经人工确认后的各 MiMeDB 代谢物条目对应的人类蛋白/酶记录数和微生物来源记录数。",
        "该图用于展示数据库关联信息的获取规模，不代表差异分析或富集分析结果。"
      )
    )

    x <- plotsAdd(x, p.stat)

    if (isTRUE(add_detail_tables)) {
      x <- tablesAdd(
        x,
        ts.microbe,
        ts.proteins
      )
    }

    n_input <- length(patterns_input)
    n_candidate <- nrow(data_match)
    n_selected <- nrow(data_selected)
    n_selected_pattern <- length(unique(data_selected$pattern))
    n_html_found <- sum(data_file$file_found, na.rm = TRUE)
    n_protein_entry <- sum(
      data_stat_metabolite$n_protein_entry,
      na.rm = TRUE
    )
    n_microbe_entry <- sum(
      data_stat_metabolite$n_microbe_entry,
      na.rm = TRUE
    )

    x <- methodAdd(
      x,
      glue::glue(
        "为从代谢物层面追溯其潜在蛋白关联与微生物来源，本研究基于 MiMeDB 数据库（<https://mimedb.org/>）",
        "对目标代谢物进行检索与人工确认。首先以输入代谢物名称为检索词，在 MiMeDB 代谢物总表中进行名称优先匹配；",
        "随后根据人工确认的 MiMeDB 条目编号获取对应代谢物页面，并整理其中的 Human Proteins and Enzymes ",
        "和 Microbial Sources 信息。对于每个经确认的代谢物条目，分别统计其关联的人类蛋白/酶记录数、",
        "唯一蛋白名称、Uniprot ID、微生物来源记录数、微生物分类信息及相关物种数量。",
        "该流程用于构建“代谢物–蛋白/酶–微生物来源”的数据库证据链。"
      )
    )

    x <- snapAdd(
      x,
      glue::glue(
        "MiMeDB 检索共输入 {n_input} 个代谢物名称，初步获得 {n_candidate} 条候选匹配记录；",
        "根据人工确认的 index 保留 {n_selected} 条真实 MiMeDB 代谢物记录，覆盖 {n_selected_pattern} 个输入代谢物，",
        "其中 {n_html_found} 条成功获取 MiMeDB 数据记录。"
      )
    )

    x <- snapAdd(
      x,
      glue::glue(
        "经整理后，本分析共获得 {n_protein_entry} 条人类蛋白/酶关联记录和 {n_microbe_entry} 条微生物来源记录，",
        "各代谢物对应记录数量见统计图{aref(p.stat)}，输入代谢物层面的匹配与检索情况见表格{aref(t.pattern)}。"
      )
    )

    return(x)
  })

setMethod("step3", signature = c(x = "job_mimedb"),
  function(x,
    layout = "anchor_orbit",
    directed = FALSE,
    include_proteins = TRUE,
    include_microbes = TRUE,
    microbe_level = c("phylum", "superkingdom_phylum"),
    label_node_types = c("Metabolite"),
    label_intersection_nodes = TRUE,
    label_intersection_types = c("Protein", "Microbe"),
    label_intersection_min_degree = 4L,
    label_wrap_width = 18L,
    label_max_nodes = 50L,
    add_network_tables = TRUE,
    seed = NULL,
    ...
  )
  {
    step_message("Build MiMeDB evidence network.")

    microbe_level <- match.arg(microbe_level)

    if (is.null(x$selected_data) || nrow(x$selected_data) == 0L) {
      stop("No selected MiMeDB data was found. Please run step2 first.")
    }

    if (is.null(x$lst_mimedb_tables) || !length(x$lst_mimedb_tables)) {
      stop("No MiMeDB detail table was found. Please run step2 first.")
    }

    if (is.null(seed)) {
      seed <- x$seed
    }

    if (is.null(seed) || length(seed) == 0L || is.na(seed)) {
      seed <- 1L
    }

    data_network <- mimeFuns$prepare_mimedb_evidence_network(
      data_items = x$selected_data,
      lst_tables = x$lst_mimedb_tables,
      include_proteins = include_proteins,
      include_microbes = include_microbes,
      microbe_level = microbe_level
    )

    data_nodes <- data_network$nodes
    data_edges <- data_network$edges

    if (nrow(data_nodes) == 0L || nrow(data_edges) == 0L) {
      warning("No valid network node or edge was obtained.")
      return(x)
    }

    graph <- mimeFuns$prepare_typed_network_graph(
      data_edges = data_edges,
      data_nodes = data_nodes,
      directed = directed
    )

    p.network <- mimeFuns$plot_typed_network(
      graph = graph,
      layout = layout,
      directed = directed,
      label_col = "label",
      type_col = "type",
      label_node_types = label_node_types,
      label_intersection_nodes = label_intersection_nodes,
      label_intersection_types = label_intersection_types,
      label_intersection_min_degree = label_intersection_min_degree,
      layout_intersection_min_degree = label_intersection_min_degree,
      layout_intersection_types = label_intersection_types,
      label_wrap_width = label_wrap_width,
      label_max_nodes = label_max_nodes,
      seed = seed,
      ...
    )

    tab_node <- table(data_nodes$type)
    tab_edge <- table(data_edges$relation)

    .get_count <- function(tab, name)
    {
      if (name %in% names(tab)) {
        return(as.integer(tab[[name]]))
      }

      0L
    }

    n_metabolite <- .get_count(tab_node, "Metabolite")
    n_protein <- .get_count(tab_node, "Protein")
    n_microbe <- .get_count(tab_node, "Microbe")
    n_edge_protein <- .get_count(tab_edge, "Metabolite-protein")
    n_edge_microbe <- .get_count(tab_edge, "Metabolite-microbe")

    data_stat <- data.frame(
      metric = c(
        "metabolite_node",
        "protein_node",
        "microbe_node",
        "metabolite_protein_edge",
        "metabolite_microbe_edge",
        "total_node",
        "total_edge"
      ),
      count = c(
        n_metabolite,
        n_protein,
        n_microbe,
        n_edge_protein,
        n_edge_microbe,
        nrow(data_nodes),
        nrow(data_edges)
      ),
      stringsAsFactors = FALSE
    )

    x$data_network_nodes <- tibble::as_tibble(data_nodes)
    x$data_network_edges <- tibble::as_tibble(data_edges)
    x$stat_network <- tibble::as_tibble(data_stat)
    x$graph_network <- graph

    if (is.null(x$lst_refine)) {
      x$lst_refine <- list()
    }

    x$lst_refine$mimedb_network <- list(
      nodes = x$data_network_nodes,
      edges = x$data_network_edges,
      stat = x$stat_network,
      graph = graph
    )

    p.network <- set_lab_legend(
      p.network,
      glue::glue("{x@sig} MiMeDB metabolite-protein-microbe evidence network"),
      glue::glue(
        "MiMeDB 代谢物-蛋白/酶-微生物来源证据网络|||",
        "该图基于经人工确认的 MiMeDB 代谢物条目构建证据链网络。",
        "图中节点表示代谢物、人类相关蛋白/酶或微生物来源；边表示 MiMeDB 页面中记录的代谢物与蛋白/酶或微生物来源之间的数据库关联证据，",
        "不代表调控方向或因果关系。为提高网络可读性，网络默认采用环绕式分层布局：代谢物节点位于外围，蛋白/酶节点位于中心区域，",
        "仅连接单一代谢物的微生物来源节点环绕对应代谢物，交叉证据节点则靠近中心区域；默认仅标注代谢物节点及连接多个代谢物的交叉证据节点。",
        "网络共包含 {nrow(data_nodes)} 个节点和 {nrow(data_edges)} 条边。"
      )
    )

    x <- plotsAdd(x, p.network)

    if (isTRUE(add_network_tables)) {
      t.network.stat <- set_lab_legend(
        x$stat_network,
        glue::glue("{x@sig} MiMeDB evidence network summary"),
        glue::glue(
          "MiMeDB 证据链网络统计表|||",
          "该表汇总 MiMeDB 代谢物-蛋白/酶-微生物来源证据网络中的节点和边数量，",
          "用于展示网络构建后的整体规模。"
        )
      )

      t.network.nodes <- set_lab_legend(
        x$data_network_nodes,
        glue::glue("{x@sig} MiMeDB evidence network nodes"),
        glue::glue(
          "MiMeDB 证据链网络节点表|||",
          "该表列出 MiMeDB 证据链网络中的全部节点，",
          "包括节点 ID、展示名称、节点类型及对应来源信息。"
        )
      )

      t.network.edges <- set_lab_legend(
        x$data_network_edges,
        glue::glue("{x@sig} MiMeDB evidence network edges"),
        glue::glue(
          "MiMeDB 证据链网络边表|||",
          "该表列出 MiMeDB 证据链网络中的全部边，",
          "包括代谢物与人类蛋白/酶或微生物来源之间的数据库关联关系。"
        )
      )

      x <- tablesAdd(
        x,
        t.network.stat = t.network.stat,
        t.network.nodes = t.network.nodes,
        t.network.edges = t.network.edges
      )
    }

    x <- methodAdd(
      x,
      glue::glue(
        "为整合 MiMeDB 数据库中的代谢物相关证据链，本研究将经人工确认的 MiMeDB 代谢物条目作为代谢物节点，",
        "将 Human Proteins and Enzymes 表中的蛋白或酶作为蛋白/酶节点，",
        "并将 Microbial Sources 表中的微生物来源信息作为微生物来源节点。",
        "若某代谢物条目中记录了相应蛋白/酶或微生物来源，则在代谢物节点与对应节点之间建立一条边，",
        "由此构建代谢物-蛋白/酶-微生物来源多类型证据网络。",
        "该网络通过 R 包 `igraph` ⟦pkgInfo('igraph')⟧ 构建，并采用 `ggraph` ⟦pkgInfo('ggraph')⟧ 和 `ggplot2` ⟦pkgInfo('ggplot2')⟧ 进行可视化。",
        "网络边表示数据库关联证据，不表示调控方向或因果效应。可视化时默认采用环绕式分层布局，将代谢物置于外围，将蛋白/酶节点置于中心区域，并将仅连接单一代谢物的微生物来源节点环绕在对应代谢物附近，以提高多类型证据链网络的可读性。"
      )
    )

    x <- snapAdd(
      x,
      glue::glue(
        "基于 MiMeDB 证据链构建代谢物-蛋白/酶-微生物来源网络{aref(p.network)}，",
        "共获得 {nrow(data_nodes)} 个节点和 {nrow(data_edges)} 条边；",
        "其中包括 {n_metabolite} 个代谢物节点、{n_protein} 个蛋白/酶节点和 {n_microbe} 个微生物来源节点。"
      )
    )

    x <- snapAdd(
      x,
      glue::glue(
        "网络中共包含 {n_edge_protein} 条代谢物-蛋白/酶关联边和 {n_edge_microbe} 条代谢物-微生物来源关联边，",
        "用于展示目标代谢物可追溯到的蛋白功能关联与微生物来源证据。"
      )
    )

    return(x)
  })





mimeFuns <- new.env(parent = emptyenv())

mimeFuns$mimedb_make_link_page <- function(ids,
  file = "material/mimedb_links.html",
  delay_ms = 800L
)
{
  ids <- unique(as.character(ids))
  ids <- ids[!is.na(ids) & ids != ""]

  if (!length(ids)) {
    stop("No valid MiMeDB ID was supplied.")
  }

  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)

  urls <- paste0("https://mimedb.org/metabolites/", ids)

  data_url <- data.frame(
    mime_id = ids,
    url = urls,
    stringsAsFactors = FALSE
  )

  html_links <- paste(
    sprintf(
      '<tr><td>%s</td><td><a href="%s" target="_blank">%s</a></td></tr>',
      data_url$mime_id,
      data_url$url,
      data_url$url
    ),
    collapse = "\n"
  )

  js_urls <- paste(
    sprintf('"%s"', data_url$url),
    collapse = ",\n"
  )

  html <- paste0(
    '<!doctype html>\n',
    '<html>\n',
    '<head>\n',
    '<meta charset="utf-8">\n',
    '<title>MiMeDB links</title>\n',
    '<style>\n',
    'body { font-family: sans-serif; margin: 24px; }\n',
    'table { border-collapse: collapse; width: 100%; }\n',
    'td, th { border: 1px solid #ddd; padding: 6px 8px; }\n',
    'button { padding: 8px 14px; margin: 0 8px 16px 0; }\n',
    '.note { color: #555; margin-bottom: 16px; }\n',
    '</style>\n',
    '</head>\n',
    '<body>\n',
    '<h2>MiMeDB metabolite pages</h2>\n',
    '<p class="note">If only one tab is opened, allow pop-ups for this local page and click again.</p>\n',
    '<button onclick="openAllTabs()">Open all pages</button>\n',
    '<button onclick="openNextTab()">Open next page</button>\n',
    '<span id="status"></span>\n',
    '<table>\n',
    '<tr><th>mime_id</th><th>url</th></tr>\n',
    html_links,
    '\n</table>\n',
    '<script>\n',
    'const urls = [\n',
    js_urls,
    '\n];\n',
    'let nextIndex = 0;\n',
    'function setStatus(text) {\n',
    '  document.getElementById("status").innerText = text;\n',
    '}\n',
    'function openAllTabs() {\n',
    '  const wins = [];\n',
    '  for (let i = 0; i < urls.length; i++) {\n',
    '    wins[i] = window.open("about:blank", "_blank");\n',
    '  }\n',
    '  let nOpened = wins.filter(w => w !== null).length;\n',
    '  setStatus("Opened blank tabs: " + nOpened + " / " + urls.length);\n',
    '  for (let i = 0; i < wins.length; i++) {\n',
    '    if (wins[i] !== null) {\n',
    '      setTimeout(() => {\n',
    '        wins[i].location.href = urls[i];\n',
    '      }, i * ',
    delay_ms,
    ');\n',
    '    }\n',
    '  }\n',
    '}\n',
    'function openNextTab() {\n',
    '  if (nextIndex >= urls.length) {\n',
    '    setStatus("All pages have been opened.");\n',
    '    return;\n',
    '  }\n',
    '  window.open(urls[nextIndex], "_blank");\n',
    '  nextIndex += 1;\n',
    '  setStatus("Opened " + nextIndex + " / " + urls.length);\n',
    '}\n',
    '</script>\n',
    '</body>\n',
    '</html>\n'
  )

  writeLines(html, con = file, useBytes = TRUE)
  utils::browseURL(normalizePath(file), "xdg-open")

  return(data_url)
}

mimeFuns$mimedb_open_urls <- function(ids,
  wait = 1.2,
  browser = c("xdg-open", "google-chrome-stable", "chromium", "firefox")
)
{
  browser <- match.arg(browser)

  ids <- unique(as.character(ids))
  ids <- ids[!is.na(ids) & ids != ""]

  if (!length(ids)) {
    stop("No valid MiMeDB ID was supplied.")
  }

  urls <- paste0("https://mimedb.org/metabolites/", ids)

  if (browser %in% c("google-chrome-stable", "chromium")) {
    for (url in urls) {
      system2(
        browser,
        args = c("--new-tab", url),
        wait = FALSE
      )
      Sys.sleep(wait)
    }
  } else if (browser == "firefox") {
    for (url in urls) {
      system2(
        browser,
        args = c("--new-tab", url),
        wait = FALSE
      )
      Sys.sleep(wait)
    }
  } else {
    for (url in urls) {
      system2(
        "xdg-open",
        args = url,
        wait = FALSE
      )
      Sys.sleep(wait)
    }
  }

  invisible(urls)
}

mimeFuns$.format_mimedb_html_data <- function(file) {
  lst <- get_table.html(readLines(file))
  lst <- lst[ !vapply(lst, is.null, logical(1)) ]
  lst <- lapply(lst, 
    function(data) {
      data <- as_tibble(data)
      if (colnames(data)[1] == "V1") {
        if (all(!is.na(names <- unlist(data[1, ])))) {
          data <- setNames(data[-1, ], names)
        }
      }
      data
    })
  colFirst <- vapply(lst, function(x) colnames(x)[1], character(1))
  names <- dplyr::recode(
    colFirst, "Superkingdom/Kingdom" = "Microbial Sources",
    "Health Outcome/Bioactivity" = "Disease",
    "Protein ID" = "Human Proteins and Enzymes"
  )
  setNames(lst, formal_name(names))
}



mimeFuns$select_existing_cols <- function(data,
  cols
)
{
  if (is.null(data)) {
    data <- data.frame(stringsAsFactors = FALSE)
  }

  data <- as.data.frame(data, stringsAsFactors = FALSE)

  for (col in cols) {
    if (!col %in% colnames(data)) {
      data[[col]] <- rep(NA_character_, nrow(data))
    }
  }

  data <- data[, cols, drop = FALSE]
  tibble::as_tibble(data)
}

mimeFuns$count_unique_non_empty <- function(x)
{
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]
  length(unique(x))
}

mimeFuns$parse_number_safe <- function(x)
{
  x <- as.character(x)
  x <- gsub(",", "", x, fixed = TRUE)
  x <- gsub("[^0-9.\\-]+", "", x, perl = TRUE)
  suppressWarnings(as.numeric(x))
}

mimeFuns$get_html_files_by_id <- function(dir = "~/Downloads",
  pattern = "\\.html?$",
  id_pattern = "MMDBc[0-9]+"
)
{
  files <- list.files(
    dir,
    pattern = pattern,
    full.names = TRUE
  )

  if (!length(files)) {
    return(stats::setNames(character(0), character(0)))
  }

  id_files <- strx(files, id_pattern)
  keep <- !is.na(id_files) & id_files != ""
  files <- files[keep]
  id_files <- id_files[keep]

  if (!length(files)) {
    return(stats::setNames(character(0), character(0)))
  }

  id_keep <- !duplicated(id_files)
  stats::setNames(files[id_keep], id_files[id_keep])
}

mimeFuns$prepare_file_status <- function(ids,
  index = NULL,
  files
)
{
  ids <- as.character(ids)

  if (is.null(index)) {
    index <- seq_along(ids)
  }

  data_file <- lapply(
    seq_along(ids),
    function(i) {
      file <- unname(files[ids[i]])

      if (!length(file)) {
        file <- NA_character_
      }

      data.frame(
        selected_index = index[i],
        mime_id = ids[i],
        file = file,
        file_found = !is.na(file) && file.exists(file),
        stringsAsFactors = FALSE
      )
    }
  )

  data_file <- do.call(rbind, data_file)
  rownames(data_file) <- NULL
  data_file
}

mimeFuns$prepare_mimedb_detail_tables <- function(file)
{
  cols_protein <- c(
    "Protein Type",
    "Protein Name",
    "Uniprot ID"
  )

  cols_microbe <- c(
    "Superkingdom/Kingdom",
    "Phylum",
    "Total species",
    "Hosts and Body Sites"
  )

  if (is.na(file) || !file.exists(file)) {
    data_protein <- mimeFuns$select_existing_cols(NULL, cols_protein)
    data_microbe <- mimeFuns$select_existing_cols(NULL, cols_microbe)

    return(list(
      Human_Proteins_and_Enzymes = data_protein,
      Microbial_Sources = data_microbe
    ))
  }

  lst <- mimeFuns$.format_mimedb_html_data(file)

  if ("Human_Proteins_and_Enzymes" %in% names(lst)) {
    data_protein <- lst$Human_Proteins_and_Enzymes
  } else {
    data_protein <- NULL
  }

  if ("Microbial_Sources" %in% names(lst)) {
    data_microbe <- lst$Microbial_Sources
  } else {
    data_microbe <- NULL
  }

  data_protein <- mimeFuns$select_existing_cols(
    data_protein,
    cols_protein
  )

  data_protein <- dplyr::arrange(
    data_protein,
    .data[["Protein Type"]]
  )

  data_protein <- dplyr::group_by(
    data_protein,
    .data[["Protein Type"]]
  )

  data_microbe <- mimeFuns$select_existing_cols(
    data_microbe,
    cols_microbe
  )

  data_microbe <- dplyr::arrange(
    data_microbe,
    .data[["Superkingdom/Kingdom"]],
    .data[["Phylum"]]
  )

  data_microbe <- dplyr::group_by(
    data_microbe,
    .data[["Superkingdom/Kingdom"]]
  )

  list(
    Human_Proteins_and_Enzymes = data_protein,
    Microbial_Sources = data_microbe
  )
}

mimeFuns$summarize_mimedb_detail_tables <- function(data_items,
  lst_tables,
  data_file
)
{
  if (nrow(data_items) == 0L) {
    return(tibble::tibble())
  }

  data_stat <- lapply(
    seq_len(nrow(data_items)),
    function(i) {
      data_protein <- lst_tables[[i]]$Human_Proteins_and_Enzymes
      data_microbe <- lst_tables[[i]]$Microbial_Sources

      vec_total_species <- mimeFuns$parse_number_safe(
        data_microbe[["Total species"]]
      )

      n_total_species <- if (length(vec_total_species) == 0L ||
          all(is.na(vec_total_species))) {
        NA_real_
      } else {
        sum(vec_total_species, na.rm = TRUE)
      }

      data.frame(
        selected_index = data_items$.selected_index[i],
        pattern = data_items$pattern[i],
        mime_id = data_items$mime_id[i],
        metabolite_name = data_items$name[i],
        match_type = data_items$match_type[i],
        file_found = data_file$file_found[i],
        n_protein_entry = nrow(data_protein),
        n_protein_type = mimeFuns$count_unique_non_empty(
          data_protein[["Protein Type"]]
        ),
        n_unique_protein = mimeFuns$count_unique_non_empty(
          data_protein[["Protein Name"]]
        ),
        n_unique_uniprot = mimeFuns$count_unique_non_empty(
          data_protein[["Uniprot ID"]]
        ),
        n_microbe_entry = nrow(data_microbe),
        n_superkingdom = mimeFuns$count_unique_non_empty(
          data_microbe[["Superkingdom/Kingdom"]]
        ),
        n_phylum = mimeFuns$count_unique_non_empty(
          data_microbe[["Phylum"]]
        ),
        n_total_species = n_total_species,
        stringsAsFactors = FALSE
      )
    }
  )

  data_stat <- do.call(rbind, data_stat)
  rownames(data_stat) <- NULL

  data_stat$metabolite_label <- ifelse(
    data_stat$pattern == data_stat$metabolite_name,
    data_stat$pattern,
    paste0(data_stat$pattern, " -> ", data_stat$metabolite_name)
  )

  tibble::as_tibble(data_stat)
}

mimeFuns$summarize_mimedb_pattern_selection <- function(patterns,
  data_candidates,
  data_selected,
  data_item_stat,
  stat_match = NULL
)
{
  patterns <- unique(as.character(patterns))
  patterns <- patterns[!is.na(patterns) & patterns != ""]

  data_candidates <- as.data.frame(
    data_candidates,
    stringsAsFactors = FALSE
  )

  data_selected <- as.data.frame(
    data_selected,
    stringsAsFactors = FALSE
  )

  data_item_stat <- as.data.frame(
    data_item_stat,
    stringsAsFactors = FALSE
  )

  data_stat <- lapply(
    patterns,
    function(pattern) {
      id_all <- which(data_candidates$pattern == pattern)
      id_selected <- which(data_selected$pattern == pattern)
      id_stat <- which(data_item_stat$pattern == pattern)

      selected_names <- unique(data_selected$name[id_selected])
      selected_names <- selected_names[
        !is.na(selected_names) &
          selected_names != ""
      ]

      selected_ids <- unique(data_selected$mime_id[id_selected])
      selected_ids <- selected_ids[
        !is.na(selected_ids) &
          selected_ids != ""
      ]

      data.frame(
        pattern = pattern,
        n_candidate_match = length(unique(data_candidates$mime_id[id_all])),
        n_selected_match = length(selected_ids),
        selected_mime_id = paste(selected_ids, collapse = "; "),
        selected_name = paste(selected_names, collapse = "; "),
        n_html_found = sum(data_item_stat$file_found[id_stat], na.rm = TRUE),
        n_protein_entry = sum(data_item_stat$n_protein_entry[id_stat], na.rm = TRUE),
        n_unique_protein_sum = sum(data_item_stat$n_unique_protein[id_stat], na.rm = TRUE),
        n_microbe_entry = sum(data_item_stat$n_microbe_entry[id_stat], na.rm = TRUE),
        n_phylum_sum = sum(data_item_stat$n_phylum[id_stat], na.rm = TRUE),
        matched_after_selection = length(selected_ids) > 0L,
        stringsAsFactors = FALSE
      )
    }
  )

  data_stat <- do.call(rbind, data_stat)
  rownames(data_stat) <- NULL

  if (!is.null(stat_match) && nrow(stat_match) > 0L) {
    stat_match <- as.data.frame(stat_match, stringsAsFactors = FALSE)
    id_match <- match(data_stat$pattern, stat_match$pattern)

    data_stat$initial_match_type <- stat_match$final_match_type[id_match]
    data_stat$initial_n_total_match <- stat_match$n_total[id_match]
    data_stat$initial_fuzzy_attempted <- stat_match$fuzzy_attempted[id_match]
    data_stat$initial_broad_attempted <- stat_match$broad_attempted[id_match]
  }

  tibble::as_tibble(data_stat)
}

mimeFuns$plot_mimedb_record_summary <- function(data_stat,
  label_size = 3.2
)
{
  data_plot <- rbind(
    data.frame(
      metabolite_label = data_stat$metabolite_label,
      metric = "Human proteins / enzymes",
      count = data_stat$n_protein_entry,
      stringsAsFactors = FALSE
    ),
    data.frame(
      metabolite_label = data_stat$metabolite_label,
      metric = "Microbial sources",
      count = data_stat$n_microbe_entry,
      stringsAsFactors = FALSE
    )
  )

  data_order <- stats::aggregate(
    count ~ metabolite_label,
    data = data_plot,
    FUN = max
  )

  data_order <- data_order[
    order(data_order$count, data_order$metabolite_label),
    ,
    drop = FALSE
  ]

  data_plot$metabolite_label <- factor(
    data_plot$metabolite_label,
    levels = data_order$metabolite_label
  )

  ggplot2::ggplot(
    data_plot,
    ggplot2::aes(
      x = metabolite_label,
      y = count,
      fill = metric
    )
  ) +
    ggplot2::geom_col(width = .72, show.legend = FALSE) +
    ggplot2::geom_text(
      ggplot2::aes(label = count),
      hjust = -.15,
      size = label_size
    ) +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(
      ~metric,
      scales = "free_y"
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, .18))
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Number of records"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey95"),
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 9)
    )
}
mimeFuns$collapse_unique_non_empty <- function(x,
  sep = "; "
)
{
  x <- as.character(x)
  x <- trimws(x)
  x <- x[!is.na(x) & x != ""]

  if (!length(x)) {
    return(NA_character_)
  }

  paste(unique(x), collapse = sep)
}

mimeFuns$make_node_table <- function(name,
  label,
  type,
  subtype = NA_character_,
  source_id = NA_character_,
  pattern = NA_character_,
  mime_id = NA_character_,
  external_id = NA_character_,
  node_detail = NA_character_
)
{
  .fit_character <- function(x, n)
  {
    x <- as.character(x)

    if (length(x) == n) {
      return(x)
    }

    if (length(x) == 1L) {
      return(rep(x, n))
    }

    if (n == 0L) {
      return(character(0))
    }

    stop("Input length is incompatible with node table length.")
  }

  n <- length(name)

  data.frame(
    name = .fit_character(name, n),
    label = .fit_character(label, n),
    type = .fit_character(type, n),
    subtype = .fit_character(subtype, n),
    source_id = .fit_character(source_id, n),
    pattern = .fit_character(pattern, n),
    mime_id = .fit_character(mime_id, n),
    external_id = .fit_character(external_id, n),
    node_detail = .fit_character(node_detail, n),
    stringsAsFactors = FALSE
  )
}

mimeFuns$make_edge_table <- function(from,
  to,
  relation,
  evidence_table,
  from_label,
  to_label,
  pattern = NA_character_,
  mime_id = NA_character_,
  edge_detail = NA_character_,
  edge_weight = NA_real_
)
{
  .fit_character <- function(x, n)
  {
    x <- as.character(x)

    if (length(x) == n) {
      return(x)
    }

    if (length(x) == 1L) {
      return(rep(x, n))
    }

    if (n == 0L) {
      return(character(0))
    }

    stop("Input length is incompatible with edge table length.")
  }

  .fit_numeric <- function(x, n)
  {
    x <- as.numeric(x)

    if (length(x) == n) {
      return(x)
    }

    if (length(x) == 1L) {
      return(rep(x, n))
    }

    if (n == 0L) {
      return(numeric(0))
    }

    stop("Input length is incompatible with edge table length.")
  }

  n <- length(from)

  data.frame(
    from = .fit_character(from, n),
    to = .fit_character(to, n),
    relation = .fit_character(relation, n),
    evidence_table = .fit_character(evidence_table, n),
    from_label = .fit_character(from_label, n),
    to_label = .fit_character(to_label, n),
    pattern = .fit_character(pattern, n),
    mime_id = .fit_character(mime_id, n),
    edge_detail = .fit_character(edge_detail, n),
    edge_weight = .fit_numeric(edge_weight, n),
    stringsAsFactors = FALSE
  )
}

mimeFuns$deduplicate_network_nodes <- function(data_nodes)
{
  data_nodes <- as.data.frame(data_nodes, stringsAsFactors = FALSE)

  if (nrow(data_nodes) == 0L) {
    return(data_nodes)
  }

  lst_nodes <- split(data_nodes, data_nodes$name)

  data_nodes <- lapply(
    lst_nodes,
    function(data_item) {
      data.frame(
        name = data_item$name[1L],
        label = data_item$label[1L],
        type = data_item$type[1L],
        subtype = mimeFuns$collapse_unique_non_empty(data_item$subtype),
        source_id = mimeFuns$collapse_unique_non_empty(data_item$source_id),
        pattern = mimeFuns$collapse_unique_non_empty(data_item$pattern),
        mime_id = mimeFuns$collapse_unique_non_empty(data_item$mime_id),
        external_id = mimeFuns$collapse_unique_non_empty(data_item$external_id),
        node_detail = mimeFuns$collapse_unique_non_empty(data_item$node_detail),
        stringsAsFactors = FALSE
      )
    }
  )

  data_nodes <- do.call(rbind, data_nodes)
  rownames(data_nodes) <- NULL
  data_nodes
}

mimeFuns$prepare_mimedb_evidence_network <- function(data_items,
  lst_tables,
  include_proteins = TRUE,
  include_microbes = TRUE,
  microbe_level = c("phylum", "superkingdom_phylum")
)
{
  microbe_level <- match.arg(microbe_level)

  data_items <- as.data.frame(data_items, stringsAsFactors = FALSE)

  if (nrow(data_items) == 0L) {
    return(list(
      nodes = mimeFuns$make_node_table(character(0), character(0), character(0)),
      edges = mimeFuns$make_edge_table(character(0), character(0), character(0), character(0), character(0), character(0))
    ))
  }

  if (!".item_label" %in% colnames(data_items)) {
    data_items$.item_label <- ifelse(
      data_items$pattern == data_items$name,
      data_items$pattern,
      paste0(data_items$pattern, " -> ", data_items$name)
    )
  }

  lst_nodes <- list()
  lst_edges <- list()
  i_node <- 0L
  i_edge <- 0L

  for (i in seq_len(nrow(data_items))) {
    item_label <- as.character(data_items$.item_label[i])

    if (is.na(item_label) || item_label == "") {
      item_label <- as.character(data_items$name[i])
    }

    mime_id <- as.character(data_items$mime_id[i])
    pattern <- as.character(data_items$pattern[i])
    match_type <- if ("match_type" %in% colnames(data_items)) {
      as.character(data_items$match_type[i])
    } else {
      NA_character_
    }

    metabolite_node <- paste0("metabolite:", mime_id)

    i_node <- i_node + 1L
    lst_nodes[[i_node]] <- mimeFuns$make_node_table(
      name = metabolite_node,
      label = item_label,
      type = "Metabolite",
      subtype = match_type,
      source_id = mime_id,
      pattern = pattern,
      mime_id = mime_id,
      external_id = mime_id,
      node_detail = as.character(data_items$name[i])
    )

    data_tables <- lst_tables[[i]]

    if (isTRUE(include_proteins) &&
        !is.null(data_tables$Human_Proteins_and_Enzymes)) {
      data_protein <- as.data.frame(
        data_tables$Human_Proteins_and_Enzymes,
        stringsAsFactors = FALSE
      )

      data_protein <- mimeFuns$select_existing_cols(
        data_protein,
        c("Protein Type", "Protein Name", "Uniprot ID")
      )

      data_protein <- as.data.frame(data_protein, stringsAsFactors = FALSE)

      protein_name <- as.character(data_protein[["Protein Name"]])
      uniprot_id <- as.character(data_protein[["Uniprot ID"]])
      protein_type <- as.character(data_protein[["Protein Type"]])

      protein_name[is.na(protein_name)] <- ""
      uniprot_id[is.na(uniprot_id)] <- ""
      protein_type[is.na(protein_type)] <- ""

      protein_label <- ifelse(
        protein_name != "",
        protein_name,
        uniprot_id
      )

      keep <- !is.na(protein_label) & protein_label != ""

      if (any(keep)) {
        protein_label <- protein_label[keep]
        protein_type <- protein_type[keep]
        uniprot_id <- uniprot_id[keep]

        protein_id <- ifelse(
          uniprot_id != "",
          paste0("protein:", uniprot_id),
          paste0("protein:", protein_label)
        )

        for (j in seq_along(protein_label)) {
          i_node <- i_node + 1L
          lst_nodes[[i_node]] <- mimeFuns$make_node_table(
            name = protein_id[j],
            label = protein_label[j],
            type = "Protein",
            subtype = protein_type[j],
            source_id = uniprot_id[j],
            pattern = pattern,
            mime_id = mime_id,
            external_id = uniprot_id[j],
            node_detail = protein_type[j]
          )

          i_edge <- i_edge + 1L
          lst_edges[[i_edge]] <- mimeFuns$make_edge_table(
            from = metabolite_node,
            to = protein_id[j],
            relation = "Metabolite-protein",
            evidence_table = "Human Proteins and Enzymes",
            from_label = item_label,
            to_label = protein_label[j],
            pattern = pattern,
            mime_id = mime_id,
            edge_detail = protein_type[j],
            edge_weight = NA_real_
          )
        }
      }
    }

    if (isTRUE(include_microbes) &&
        !is.null(data_tables$Microbial_Sources)) {
      data_microbe <- as.data.frame(
        data_tables$Microbial_Sources,
        stringsAsFactors = FALSE
      )

      data_microbe <- mimeFuns$select_existing_cols(
        data_microbe,
        c(
          "Superkingdom/Kingdom",
          "Phylum",
          "Total species",
          "Hosts and Body Sites"
        )
      )

      data_microbe <- as.data.frame(data_microbe, stringsAsFactors = FALSE)

      superkingdom <- as.character(data_microbe[["Superkingdom/Kingdom"]])
      phylum <- as.character(data_microbe[["Phylum"]])
      total_species <- mimeFuns$parse_number_safe(
        data_microbe[["Total species"]]
      )
      host_site <- as.character(data_microbe[["Hosts and Body Sites"]])

      superkingdom[is.na(superkingdom)] <- ""
      phylum[is.na(phylum)] <- ""
      host_site[is.na(host_site)] <- ""

      microbe_label <- if (microbe_level == "superkingdom_phylum") {
        ifelse(
          superkingdom != "" & phylum != "",
          paste0(superkingdom, " | ", phylum),
          ifelse(phylum != "", phylum, superkingdom)
        )
      } else {
        ifelse(phylum != "", phylum, superkingdom)
      }

      keep <- !is.na(microbe_label) & microbe_label != ""

      if (any(keep)) {
        microbe_label <- microbe_label[keep]
        superkingdom <- superkingdom[keep]
        phylum <- phylum[keep]
        total_species <- total_species[keep]
        host_site <- host_site[keep]
        microbe_id <- paste0("microbe:", microbe_label)

        for (j in seq_along(microbe_label)) {
          node_detail <- paste(
            c(superkingdom[j], phylum[j]),
            collapse = "; "
          )

          i_node <- i_node + 1L
          lst_nodes[[i_node]] <- mimeFuns$make_node_table(
            name = microbe_id[j],
            label = microbe_label[j],
            type = "Microbe",
            subtype = superkingdom[j],
            source_id = phylum[j],
            pattern = pattern,
            mime_id = mime_id,
            external_id = phylum[j],
            node_detail = node_detail
          )

          i_edge <- i_edge + 1L
          lst_edges[[i_edge]] <- mimeFuns$make_edge_table(
            from = metabolite_node,
            to = microbe_id[j],
            relation = "Metabolite-microbe",
            evidence_table = "Microbial Sources",
            from_label = item_label,
            to_label = microbe_label[j],
            pattern = pattern,
            mime_id = mime_id,
            edge_detail = host_site[j],
            edge_weight = total_species[j]
          )
        }
      }
    }
  }

  if (length(lst_nodes) > 0L) {
    data_nodes <- do.call(rbind, lst_nodes)
    data_nodes <- mimeFuns$deduplicate_network_nodes(data_nodes)
  } else {
    data_nodes <- mimeFuns$make_node_table(
      character(0),
      character(0),
      character(0)
    )
  }

  if (length(lst_edges) > 0L) {
    data_edges <- do.call(rbind, lst_edges)
    data_edges <- unique(data_edges)
    rownames(data_edges) <- NULL
  } else {
    data_edges <- mimeFuns$make_edge_table(
      character(0),
      character(0),
      character(0),
      character(0),
      character(0),
      character(0)
    )
  }

  list(
    nodes = tibble::as_tibble(data_nodes),
    edges = tibble::as_tibble(data_edges)
  )
}

mimeFuns$prepare_typed_network_graph <- function(data_edges,
  data_nodes,
  directed = FALSE
)
{
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop('Package "igraph" is required.')
  }

  data_edges <- as.data.frame(data_edges, stringsAsFactors = FALSE)
  data_nodes <- as.data.frame(data_nodes, stringsAsFactors = FALSE)

  vec_edge_need <- c("from", "to")
  vec_node_need <- c("name", "label", "type")

  vec_edge_miss <- setdiff(vec_edge_need, colnames(data_edges))
  vec_node_miss <- setdiff(vec_node_need, colnames(data_nodes))

  if (length(vec_edge_miss) > 0L) {
    stop(glue::glue(
      "Missing edge column(s): {paste(vec_edge_miss, collapse = ', ')}."
    ))
  }

  if (length(vec_node_miss) > 0L) {
    stop(glue::glue(
      "Missing node column(s): {paste(vec_node_miss, collapse = ', ')}."
    ))
  }

  igraph::graph_from_data_frame(
    d = data_edges,
    directed = directed,
    vertices = data_nodes
  )
}

mimeFuns$wrap_text_vector <- function(x,
  width = 28L
)
{
  x <- as.character(x)

  if (is.null(width) || length(width) == 0L || is.na(width) ||
      !is.finite(width) || width <= 0L) {
    return(x)
  }

  vapply(
    x,
    function(text) {
      if (is.na(text) || text == "") {
        return(text)
      }

      paste(strwrap(text, width = width), collapse = "\n")
    },
    character(1)
  )
}


mimeFuns$get_coords_spiral <- function(num = 100L,
  width_curve = 1,
  n_cir = NULL,
  min_rad = 1,
  max_rad = 3,
  rank_by = NULL,
  desc = TRUE,
  center_x = 0,
  center_y = 0,
  angle_offset = 0
)
{
  num <- as.integer(num)

  if (num <= 0L) {
    return(data.frame(x = numeric(0), y = numeric(0)))
  }

  if (is.null(n_cir) || length(n_cir) == 0L || is.na(n_cir)) {
    n_cir <- max(1L, num %/% 20L)
  }

  ang <- seq(
    angle_offset,
    by = 2 * pi * n_cir / max(num, 1L),
    length.out = num
  )

  t <- seq(0, 1, length.out = num)
  t_transformed <- t ^ width_curve
  rad <- min_rad + (max_rad - min_rad) * t_transformed

  data <- data.frame(
    x = sin(ang) * rad + center_x,
    y = cos(ang) * rad + center_y,
    stringsAsFactors = FALSE
  )

  if (!is.null(rank_by)) {
    rank_by <- as.numeric(rank_by)

    if (length(rank_by) != num) {
      stop("The length of rank_by should be equal to num.")
    }

    if (isTRUE(desc)) {
      rank_by <- -rank_by
    }

    rank <- rank(rank_by, ties.method = "first", na.last = "keep")
    data <- data[rank, , drop = FALSE]
  }

  rownames(data) <- NULL
  data
}

mimeFuns$get_named_layout_value <- function(x,
  name,
  default
)
{
  if (is.null(x)) {
    return(default)
  }

  if (is.list(x)) {
    if (name %in% names(x)) {
      return(x[[name]])
    }

    return(default)
  }

  if (is.atomic(x) && !is.null(names(x))) {
    if (name %in% names(x)) {
      return(x[[name]])
    }

    return(default)
  }

  default
}

mimeFuns$prepare_typed_spiral_layout <- function(graph,
  type_col = "type",
  layout_intersection_types = NULL,
  layout_intersection_min_degree = 3L,
  layout_type_centers = NULL,
  layout_type_min_rad = NULL,
  layout_type_max_rad = NULL,
  layout_type_width_curve = NULL,
  layout_type_n_cir = NULL,
  layout_seed = 1L
)
{
  if (!requireNamespace("ggraph", quietly = TRUE)) {
    stop('Package "ggraph" is required.')
  }

  set.seed(layout_seed)

  data_nodes <- igraph::as_data_frame(
    graph,
    what = "vertices"
  )

  if (!type_col %in% colnames(data_nodes)) {
    data_nodes[[type_col]] <- "Node"
  }

  if (!"name" %in% colnames(data_nodes)) {
    data_nodes$name <- igraph::V(graph)$name
  }

  vec_degree <- igraph::degree(graph)
  data_nodes$.degree <- as.numeric(
    vec_degree[match(data_nodes$name, names(vec_degree))]
  )

  data_nodes$.degree[is.na(data_nodes$.degree)] <- 0

  if (is.null(layout_intersection_types)) {
    layout_intersection_types <- unique(as.character(data_nodes[[type_col]]))
  }

  data_nodes$.layout_intersection <-
    data_nodes[[type_col]] %in% layout_intersection_types &
    data_nodes$.degree >= layout_intersection_min_degree

  data_nodes$.layout_group <- as.character(data_nodes[[type_col]])
  data_nodes$.layout_group[data_nodes$.layout_intersection] <- "Intersection"

  default_centers <- list(
    Metabolite = c(0, .75),
    Protein = c(1.75, -.25),
    Microbe = c(-1.2, -.45),
    Intersection = c(0, -.05),
    Node = c(0, 0)
  )

  if (is.null(layout_type_centers)) {
    layout_type_centers <- default_centers
  } else {
    for (name in names(default_centers)) {
      if (!name %in% names(layout_type_centers)) {
        layout_type_centers[[name]] <- default_centers[[name]]
      }
    }
  }

  default_min_rad <- c(
    Metabolite = .08,
    Protein = .12,
    Microbe = .08,
    Intersection = .75,
    Node = .08
  )

  default_max_rad <- c(
    Metabolite = .45,
    Protein = 1.35,
    Microbe = .65,
    Intersection = 1.35,
    Node = 1
  )

  default_width_curve <- c(
    Metabolite = .9,
    Protein = .75,
    Microbe = .8,
    Intersection = .9,
    Node = 1
  )

  data_coord <- data.frame(
    name = data_nodes$name,
    x = 0,
    y = 0,
    stringsAsFactors = FALSE
  )

  groups <- unique(data_nodes$.layout_group)
  groups <- groups[!is.na(groups) & groups != ""]

  for (group in groups) {
    id_group <- which(data_nodes$.layout_group == group)
    n_group <- length(id_group)

    center <- mimeFuns$get_named_layout_value(
      layout_type_centers,
      group,
      c(0, 0)
    )

    center <- as.numeric(center)

    if (length(center) < 2L || any(is.na(center[1:2]))) {
      center <- c(0, 0)
    }

    min_rad <- mimeFuns$get_named_layout_value(
      layout_type_min_rad,
      group,
      mimeFuns$get_named_layout_value(
        default_min_rad,
        group,
        .08
      )
    )

    max_rad <- mimeFuns$get_named_layout_value(
      layout_type_max_rad,
      group,
      mimeFuns$get_named_layout_value(
        default_max_rad,
        group,
        1
      )
    )

    width_curve <- mimeFuns$get_named_layout_value(
      layout_type_width_curve,
      group,
      mimeFuns$get_named_layout_value(
        default_width_curve,
        group,
        1
      )
    )

    n_cir <- mimeFuns$get_named_layout_value(
      layout_type_n_cir,
      group,
      max(1L, n_group %/% 18L)
    )

    coord <- mimeFuns$get_coords_spiral(
      num = n_group,
      width_curve = width_curve,
      n_cir = n_cir,
      min_rad = min_rad,
      max_rad = max_rad,
      rank_by = data_nodes$.degree[id_group],
      desc = TRUE,
      center_x = center[1L],
      center_y = center[2L],
      angle_offset = stats::runif(1L, 0, 2 * pi)
    )

    data_coord$x[id_group] <- coord$x
    data_coord$y[id_group] <- coord$y
  }

  ggraph::create_layout(
    graph,
    layout = "manual",
    x = data_coord$x,
    y = data_coord$y
  )
}


mimeFuns$get_coords_circle <- function(num = 100L,
  radius = 1,
  center_x = 0,
  center_y = 0,
  angle_offset = 0,
  rank_by = NULL,
  desc = TRUE
)
{
  num <- as.integer(num)

  if (num <= 0L) {
    return(data.frame(x = numeric(0), y = numeric(0)))
  }

  ang <- seq(
    angle_offset,
    angle_offset + 2 * pi,
    length.out = num + 1L
  )

  ang <- ang[seq_len(num)]

  data <- data.frame(
    x = cos(ang) * radius + center_x,
    y = sin(ang) * radius + center_y,
    stringsAsFactors = FALSE
  )

  if (!is.null(rank_by)) {
    rank_by <- as.numeric(rank_by)

    if (length(rank_by) != num) {
      stop("The length of rank_by should be equal to num.")
    }

    if (isTRUE(desc)) {
      rank_by <- -rank_by
    }

    rank <- rank(rank_by, ties.method = "first", na.last = "keep")
    data <- data[rank, , drop = FALSE]
  }

  rownames(data) <- NULL
  data
}

mimeFuns$get_coords_grid <- function(num = 100L,
  cell_width = .18,
  cell_height = .18,
  center_x = 0,
  center_y = 0,
  rank_by = NULL,
  desc = TRUE
)
{
  num <- as.integer(num)

  if (num <= 0L) {
    return(data.frame(x = numeric(0), y = numeric(0)))
  }

  n_col <- ceiling(sqrt(num))
  n_row <- ceiling(num / n_col)

  data <- expand.grid(
    col = seq_len(n_col),
    row = seq_len(n_row)
  )

  data <- data[seq_len(num), , drop = FALSE]

  data$x <- (data$col - mean(range(data$col))) * cell_width + center_x
  data$y <- (data$row - mean(range(data$row))) * cell_height + center_y

  data <- data[, c("x", "y"), drop = FALSE]

  if (!is.null(rank_by)) {
    rank_by <- as.numeric(rank_by)

    if (length(rank_by) != num) {
      stop("The length of rank_by should be equal to num.")
    }

    if (isTRUE(desc)) {
      rank_by <- -rank_by
    }

    rank <- rank(rank_by, ties.method = "first", na.last = "keep")
    data <- data[rank, , drop = FALSE]
  }

  rownames(data) <- NULL
  data
}

mimeFuns$get_neighbor_names <- function(graph,
  node
)
{
  id_node <- which(igraph::V(graph)$name == node)

  if (length(id_node) != 1L) {
    return(character(0))
  }

  vec_neighbor <- igraph::neighbors(
    graph,
    v = id_node,
    mode = "all"
  )

  igraph::as_ids(vec_neighbor)
}

mimeFuns$make_circle_path <- function(data_center,
  x_col = "x",
  y_col = "y",
  radius_col = "radius",
  n = 120L
)
{
  data_center <- as.data.frame(data_center, stringsAsFactors = FALSE)

  if (nrow(data_center) == 0L) {
    return(data.frame(x = numeric(0), y = numeric(0), group = character(0)))
  }

  theta <- seq(0, 2 * pi, length.out = n + 1L)

  data_path <- lapply(
    seq_len(nrow(data_center)),
    function(i) {
      radius <- as.numeric(data_center[[radius_col]][i])

      if (!is.finite(radius) || radius <= 0) {
        radius <- 1
      }

      data.frame(
        x = as.numeric(data_center[[x_col]][i]) + cos(theta) * radius,
        y = as.numeric(data_center[[y_col]][i]) + sin(theta) * radius,
        group = paste0("circle_", i),
        stringsAsFactors = FALSE
      )
    }
  )

  data_path <- do.call(rbind, data_path)
  rownames(data_path) <- NULL
  data_path
}

mimeFuns$prepare_anchor_orbit_layout <- function(graph,
  type_col = "type",
  anchor_type = "Metabolite",
  core_type = "Protein",
  satellite_type = "Microbe",
  shared_types = c("Protein", "Microbe"),
  shared_min_degree = 2L,
  anchor_radius = 4.2,
  anchor_angle_offset = pi / 2,
  core_layout = c("spiral", "grid"),
  core_min_rad = .12,
  core_max_rad = 1.55,
  core_grid_cell_width = .18,
  core_grid_cell_height = .18,
  shared_min_rad = 1.85,
  shared_max_rad = 2.55,
  satellite_min_rad = .24,
  satellite_max_rad = .78,
  satellite_width_curve = .82,
  other_min_rad = 2.7,
  other_max_rad = 3.1,
  layout_seed = 1L
)
{
  if (!requireNamespace("ggraph", quietly = TRUE)) {
    stop('Package "ggraph" is required.')
  }

  core_layout <- match.arg(core_layout)
  set.seed(layout_seed)

  data_nodes <- igraph::as_data_frame(
    graph,
    what = "vertices"
  )

  if (!type_col %in% colnames(data_nodes)) {
    data_nodes[[type_col]] <- "Node"
  }

  if (!"name" %in% colnames(data_nodes)) {
    data_nodes$name <- igraph::V(graph)$name
  }

  vec_degree <- igraph::degree(graph)
  data_nodes$.degree <- as.numeric(
    vec_degree[match(data_nodes$name, names(vec_degree))]
  )

  data_nodes$.degree[is.na(data_nodes$.degree)] <- 0

  data_coord <- data.frame(
    name = data_nodes$name,
    x = 0,
    y = 0,
    .layout_role = "other",
    .anchor_name = NA_character_,
    stringsAsFactors = FALSE
  )

  id_anchor <- which(data_nodes[[type_col]] == anchor_type)
  id_core <- which(data_nodes[[type_col]] == core_type)
  id_satellite <- which(data_nodes[[type_col]] == satellite_type)

  if (length(id_anchor) > 0L) {
    coord_anchor <- mimeFuns$get_coords_circle(
      num = length(id_anchor),
      radius = anchor_radius,
      center_x = 0,
      center_y = 0,
      angle_offset = anchor_angle_offset,
      rank_by = data_nodes$.degree[id_anchor],
      desc = TRUE
    )

    data_coord$x[id_anchor] <- coord_anchor$x
    data_coord$y[id_anchor] <- coord_anchor$y
    data_coord$.layout_role[id_anchor] <- "anchor"
  }

  if (length(id_core) > 0L) {
    if (core_layout == "grid") {
      coord_core <- mimeFuns$get_coords_grid(
        num = length(id_core),
        cell_width = core_grid_cell_width,
        cell_height = core_grid_cell_height,
        center_x = 0,
        center_y = 0,
        rank_by = data_nodes$.degree[id_core],
        desc = TRUE
      )
    } else {
      coord_core <- mimeFuns$get_coords_spiral(
        num = length(id_core),
        width_curve = .72,
        n_cir = max(2L, length(id_core) %/% 26L),
        min_rad = core_min_rad,
        max_rad = core_max_rad,
        rank_by = data_nodes$.degree[id_core],
        desc = TRUE,
        center_x = 0,
        center_y = 0,
        angle_offset = stats::runif(1L, 0, 2 * pi)
      )
    }

    data_coord$x[id_core] <- coord_core$x
    data_coord$y[id_core] <- coord_core$y
    data_coord$.layout_role[id_core] <- "core"
  }

  node_index <- stats::setNames(seq_len(nrow(data_nodes)), data_nodes$name)
  anchor_names <- data_nodes$name[id_anchor]
  lst_satellite_by_anchor <- stats::setNames(
    vector("list", length(anchor_names)),
    anchor_names
  )

  id_shared <- integer(0)
  id_other <- integer(0)

  for (id in id_satellite) {
    node_name <- data_nodes$name[id]
    neighbor_names <- mimeFuns$get_neighbor_names(
      graph = graph,
      node = node_name
    )

    neighbor_id_raw <- node_index[neighbor_names]
    valid_neighbor <- !is.na(neighbor_id_raw)
    neighbor_names_valid <- neighbor_names[valid_neighbor]
    neighbor_id <- node_index[neighbor_names_valid]
    neighbor_id <- neighbor_id[!is.na(neighbor_id)]

    neighbor_anchor <- neighbor_names_valid[
      data_nodes[[type_col]][neighbor_id] == anchor_type
    ]

    if (length(neighbor_anchor) == 1L && data_nodes$.degree[id] < shared_min_degree) {
      anchor_name <- neighbor_anchor[1L]
      lst_satellite_by_anchor[[anchor_name]] <- c(
        lst_satellite_by_anchor[[anchor_name]],
        id
      )
      data_coord$.anchor_name[id] <- anchor_name
    } else if (data_nodes[[type_col]][id] %in% shared_types ||
        data_nodes$.degree[id] >= shared_min_degree) {
      id_shared <- c(id_shared, id)
    } else {
      id_other <- c(id_other, id)
    }
  }

  for (anchor_name in names(lst_satellite_by_anchor)) {
    id_group <- lst_satellite_by_anchor[[anchor_name]]

    if (!length(id_group)) {
      next
    }

    id_anchor_one <- node_index[[anchor_name]]
    n_group <- length(id_group)
    extra_rad <- .05 * sqrt(max(n_group - 8L, 0L))

    coord_satellite <- mimeFuns$get_coords_spiral(
      num = n_group,
      width_curve = satellite_width_curve,
      n_cir = max(1L, n_group %/% 14L),
      min_rad = satellite_min_rad,
      max_rad = satellite_max_rad + extra_rad,
      rank_by = data_nodes$.degree[id_group],
      desc = TRUE,
      center_x = data_coord$x[id_anchor_one],
      center_y = data_coord$y[id_anchor_one],
      angle_offset = stats::runif(1L, 0, 2 * pi)
    )

    data_coord$x[id_group] <- coord_satellite$x
    data_coord$y[id_group] <- coord_satellite$y
    data_coord$.layout_role[id_group] <- "satellite"
  }

  id_shared <- setdiff(unique(id_shared), id_core)

  if (length(id_shared) > 0L) {
    coord_shared <- mimeFuns$get_coords_spiral(
      num = length(id_shared),
      width_curve = .95,
      n_cir = max(1L, length(id_shared) %/% 14L),
      min_rad = shared_min_rad,
      max_rad = shared_max_rad,
      rank_by = data_nodes$.degree[id_shared],
      desc = TRUE,
      center_x = 0,
      center_y = 0,
      angle_offset = stats::runif(1L, 0, 2 * pi)
    )

    data_coord$x[id_shared] <- coord_shared$x
    data_coord$y[id_shared] <- coord_shared$y
    data_coord$.layout_role[id_shared] <- "shared"
  }

  id_done <- which(data_coord$.layout_role != "other")
  id_remaining <- setdiff(seq_len(nrow(data_nodes)), id_done)
  id_remaining <- union(id_remaining, id_other)
  id_remaining <- unique(id_remaining)

  if (length(id_remaining) > 0L) {
    coord_other <- mimeFuns$get_coords_spiral(
      num = length(id_remaining),
      width_curve = 1,
      n_cir = max(1L, length(id_remaining) %/% 16L),
      min_rad = other_min_rad,
      max_rad = other_max_rad,
      rank_by = data_nodes$.degree[id_remaining],
      desc = TRUE,
      center_x = 0,
      center_y = 0,
      angle_offset = stats::runif(1L, 0, 2 * pi)
    )

    data_coord$x[id_remaining] <- coord_other$x
    data_coord$y[id_remaining] <- coord_other$y
    data_coord$.layout_role[id_remaining] <- "other"
  }

  ggraph::create_layout(
    graph,
    layout = "manual",
    x = data_coord$x,
    y = data_coord$y
  )
}

mimeFuns$plot_typed_network <- function(graph,
  layout = "anchor_orbit",
  directed = FALSE,
  label_col = "label",
  type_col = "type",
  label_node_types = c("Metabolite"),
  label_intersection_nodes = TRUE,
  label_intersection_types = NULL,
  label_intersection_min_degree = 4L,
  label_wrap_width = 18L,
  label_type_priority = c("Metabolite", "Protein", "Microbe"),
  label_max_nodes = 50L,
  layout_intersection_types = NULL,
  layout_intersection_min_degree = 3L,
  layout_type_centers = NULL,
  layout_type_min_rad = NULL,
  layout_type_max_rad = NULL,
  layout_type_width_curve = NULL,
  layout_type_n_cir = NULL,
  layout_anchor_type = "Metabolite",
  layout_core_type = "Protein",
  layout_satellite_type = "Microbe",
  layout_shared_types = c("Protein", "Microbe"),
  layout_shared_min_degree = 2L,
  layout_anchor_radius = 4.2,
  layout_core_layout = c("spiral", "grid"),
  layout_core_min_rad = .12,
  layout_core_max_rad = 1.55,
  layout_shared_min_rad = 1.85,
  layout_shared_max_rad = 2.55,
  layout_satellite_min_rad = .24,
  layout_satellite_max_rad = .78,
  show_anchor_orbits = TRUE,
  anchor_orbit_linetype = 3,
  anchor_orbit_color = "grey75",
  anchor_orbit_alpha = .5,
  node_size_values = c(
    Metabolite = 7,
    Protein = 2.4,
    Microbe = 4.2
  ),
  node_shape_values = c(
    Metabolite = 16,
    Protein = 17,
    Microbe = 15
  ),
  edge_color = "grey78",
  edge_width = .25,
  edge_alpha = .28,
  label_size = 2.8,
  label_force = .45,
  label_force_pull = 2.5,
  seed = 1L
)
{
  if (!requireNamespace("ggraph", quietly = TRUE)) {
    stop('Package "ggraph" is required.')
  }

  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop('Package "ggrepel" is required.')
  }

  layout_core_layout <- match.arg(layout_core_layout)

  set.seed(seed)

  if (identical(layout, "anchor_orbit") ||
      identical(layout, "metabolite_orbit")) {
    layout_graph <- mimeFuns$prepare_anchor_orbit_layout(
      graph = graph,
      type_col = type_col,
      anchor_type = layout_anchor_type,
      core_type = layout_core_type,
      satellite_type = layout_satellite_type,
      shared_types = layout_shared_types,
      shared_min_degree = layout_shared_min_degree,
      anchor_radius = layout_anchor_radius,
      core_layout = layout_core_layout,
      core_min_rad = layout_core_min_rad,
      core_max_rad = layout_core_max_rad,
      shared_min_rad = layout_shared_min_rad,
      shared_max_rad = layout_shared_max_rad,
      satellite_min_rad = layout_satellite_min_rad,
      satellite_max_rad = layout_satellite_max_rad,
      layout_seed = seed
    )
  } else if (identical(layout, "typed_spiral") ||
      identical(layout, "grouped_spiral")) {
    layout_graph <- mimeFuns$prepare_typed_spiral_layout(
      graph = graph,
      type_col = type_col,
      layout_intersection_types = layout_intersection_types,
      layout_intersection_min_degree = layout_intersection_min_degree,
      layout_type_centers = layout_type_centers,
      layout_type_min_rad = layout_type_min_rad,
      layout_type_max_rad = layout_type_max_rad,
      layout_type_width_curve = layout_type_width_curve,
      layout_type_n_cir = layout_type_n_cir,
      layout_seed = seed
    )
  } else {
    layout_graph <- ggraph::create_layout(
      graph,
      layout = layout
    )
  }

  vec_degree <- igraph::degree(graph)
  layout_graph$.degree <- vec_degree[match(layout_graph$name, names(vec_degree))]
  layout_graph$.degree[is.na(layout_graph$.degree)] <- 0

  if (isTRUE(directed)) {
    layer_edge <- ggraph::geom_edge_link(
      edge_width = edge_width,
      color = edge_color,
      alpha = edge_alpha,
      show.legend = FALSE,
      end_cap = ggraph::circle(5, "mm"),
      arrow = grid::arrow(length = grid::unit(.9, "mm"))
    )
  } else {
    layer_edge <- ggraph::geom_edge_link(
      edge_width = edge_width,
      color = edge_color,
      alpha = edge_alpha,
      show.legend = FALSE
    )
  }

  p <- ggraph::ggraph(layout_graph) +
    layer_edge

  if (isTRUE(show_anchor_orbits) &&
      (identical(layout, "anchor_orbit") || identical(layout, "metabolite_orbit")) &&
      type_col %in% colnames(layout_graph)) {
    data_anchor <- as.data.frame(layout_graph)
    data_anchor <- data_anchor[
      data_anchor[[type_col]] == layout_anchor_type,
      ,
      drop = FALSE
    ]

    if (nrow(data_anchor) > 0L) {
      data_anchor$radius <- layout_satellite_max_rad * 1.15
      data_circle <- mimeFuns$make_circle_path(
        data_center = data_anchor,
        x_col = "x",
        y_col = "y",
        radius_col = "radius"
      )

      p <- p +
        ggplot2::geom_path(
          data = data_circle,
          inherit.aes = FALSE,
          color = anchor_orbit_color,
          alpha = anchor_orbit_alpha,
          linetype = anchor_orbit_linetype,
          linewidth = .25,
          ggplot2::aes(
            x = x,
            y = y,
            group = group
          )
        )
    }
  }

  p <- p +
    ggraph::geom_node_point(
      ggplot2::aes(
        color = .data[[type_col]],
        size = .data[[type_col]],
        shape = .data[[type_col]]
      ),
      alpha = .9,
      stroke = .3
    ) +
    ggplot2::scale_size_manual(
      values = node_size_values
    ) +
    ggplot2::scale_shape_manual(
      values = node_shape_values
    ) +
    ggplot2::guides(
      size = "none",
      shape = "none",
      color = ggplot2::guide_legend(
        override.aes = list(size = 4)
      )
    ) +
    ggplot2::labs(color = "Type") +
    ggplot2::theme_void()

  data_label <- as.data.frame(layout_graph)

  if (!label_col %in% colnames(data_label)) {
    label_col <- "name"
  }

  if (!type_col %in% colnames(data_label)) {
    type_col <- "type"
  }

  if (is.null(label_intersection_types)) {
    label_intersection_types <- unique(as.character(data_label[[type_col]]))
  }

  data_label$.label_by_type <- data_label[[type_col]] %in% label_node_types
  data_label$.label_by_intersection <- isTRUE(label_intersection_nodes) &
    data_label[[type_col]] %in% label_intersection_types &
    data_label$.degree >= label_intersection_min_degree

  data_label <- data_label[
    data_label$.label_by_type |
      data_label$.label_by_intersection,
    ,
    drop = FALSE
  ]

  data_label$.label <- as.character(data_label[[label_col]])
  data_label <- data_label[
    !is.na(data_label$.label) &
      data_label$.label != "",
    ,
    drop = FALSE
  ]

  if (nrow(data_label) > 0L) {
    data_label$.label <- mimeFuns$wrap_text_vector(
      data_label$.label,
      width = label_wrap_width
    )
  }

  if (nrow(data_label) > 0L && !is.null(label_max_nodes) &&
      length(label_max_nodes) > 0L && is.finite(label_max_nodes) &&
      nrow(data_label) > label_max_nodes) {
    data_label$.type_priority <- match(
      data_label[[type_col]],
      label_type_priority
    )

    data_label$.type_priority[is.na(data_label$.type_priority)] <-
      length(label_type_priority) + 1L

    data_label$.label_reason_priority <- ifelse(
      data_label$.label_by_type,
      1L,
      2L
    )

    data_label <- data_label[
      order(
        data_label$.label_reason_priority,
        data_label$.type_priority,
        -data_label$.degree,
        data_label$.label
      ),
      ,
      drop = FALSE
    ]

    data_label <- data_label[
      seq_len(min(label_max_nodes, nrow(data_label))),
      ,
      drop = FALSE
    ]
  }

  if (nrow(data_label) > 0L) {
    p <- p +
      ggrepel::geom_label_repel(
        data = data_label,
        inherit.aes = FALSE,
        seed = seed,
        size = label_size,
        force = label_force,
        force_pull = label_force_pull,
        max.overlaps = Inf,
        box.padding = .12,
        point.padding = .08,
        label.size = .15,
        min.segment.length = 0,
        segment.alpha = .5,
        ggplot2::aes(
          x = x,
          y = y,
          label = .label
        )
      )
  }

  p
}

