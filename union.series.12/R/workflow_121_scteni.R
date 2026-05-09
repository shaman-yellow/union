# ==========================================================================
# workflow of scteni
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_scteni <- setClass("job_scteni", 
  contains = c("job"),
  prototype = prototype(
    pg = "scteni",
    info = c("https://github.com/cailab-tamu/scTenifoldKnk"),
    cite = "",
    method = "",
    tag = "scteni",
    analysis = "scTenifoldKnk 虚拟敲除"
    ))


setGeneric("asjob_scteni",
  function(x, ...) standardGeneric("asjob_scteni"))

setMethod("asjob_scteni", signature = c(x = "ANY"),
  function(x, ref){
    pr <- list()
    if (is(x, "job_seurat")) {
      pr <- params(x)
      if (any(!ref %in% rownames(object(x)))) {
        stop('any(!ref %in% rownames(object(x))).')
      }
      pr$dir_seurat <- create_job_cache_dir(x, "sctenifoldknk", path = ".")
      pr$file_seurat <- file.path(pr$dir_seurat, glue::glue("seurat_{sig(x)}.rds"))
      if (!file.exists(pr$file_seurat)) {
        message(glue::glue("Save file: {pr$file_seurat}"))
        saveRDS(object(x), pr$file_seurat)
      }
    }
    fea <- resolve_feature_snapAdd_onExit("x", ref)
    x <- .job_scteni(object = fea)
    x@params <- append(x@params, pr)
    x$features <- ref
    return(x)
  })

setMethod("step0", signature = c(x = "job_scteni"),
  function(x){
    step_message("Prepare your data with function `job_scteni`.")
  })

setMethod("step1", signature = c(x = "job_scteni"),
  function(x){
    step_message("Do nothing")
    x <- methodAdd(x, "**scTenifoldKnk** 是一种基于单细胞转录组数据构建基因调控网络并进行 *in silico* 虚拟敲除的计算方法，其核心目的是在无需实际基因编辑实验的情况下，评估特定基因在细胞系统中的功能重要性及其对全局转录调控网络的影响。具体而言，该方法通过对目标基因进行网络层面的“移除”，并比较敲除前后基因调控网络结构及表达模式的变化，从而识别潜在的下游调控通路及受影响的关键基因模块。该分析有助于优先筛选具有重要生物学意义的候选关键基因，为后续机制研究和实验验证提供理论依据。")
    x <- methodAdd(x, "以 R 包 `scTenifoldKnk` ⟦pkgInfo('scTenifoldKnk')⟧ 模拟敲除分析，以推测关键基因下游功能影响。构建好的基因调控网络（gene regulatory network，GRN）中进行无监督虚拟敲低，通过模拟关键基因节点的缺失来模拟其敲低状态，进而筛选出因关键基因敲低而调控关系发生显著改变的基因。")
    return(x)
  })

setMethod("step2", signature = c(x = "job_scteni"),
  function(x, use.p = c("p.adj", "p.value"), cut.p = .05, 
    cut.z = 2, recode = NULL, dir = x$dir_seurat,
    lst_diff = .read_scteni_results(dir, sig(x)))
  {
    step_message("Load results file and draw volcano plot.")
    use.p <- match.arg(use.p)
    if (!is.null(recode)) {
      names(lst_diff) <- dplyr::recode(names(lst_diff), !!!recode)
      lst_diff_no_filter <- lapply(lst_diff, 
        function(data) {
          dplyr::mutate(data, gene = dplyr::recode(gene, !!!recode))
        })
    } else {
      lst_diff_no_filter <- lst_diff
    }
    lst_diff <- sapply(names(lst_diff_no_filter), simplify = FALSE,
      function(name) {
        data <- lst_diff_no_filter[[name]]
        message(glue::glue("All data, nrow: {nrow(data)}"))
        data <- dplyr::filter(
          data, gene != !!name, !!rlang::sym(use.p) < !!cut.p
        )
        message(glue::glue("After filter by {use.p}, nrow: {nrow(data)}"))
        if (!is.null(cut.z)) {
          data <- dplyr::filter(data, abs(Z) > !!cut.z)
        }
        message(glue::glue("After filter by |Z|, nrow: {nrow(data)}"))
        data
      })
    lst_diff <- set_lab_legend(
      lst_diff,
      glue::glue("{x@sig} data of genes affected by {names(lst_diff)} knockout"),
      glue::glue("受 {names(lst_diff)} 敲除影响最显著的基因")
    )
    x <- tablesAdd(x, t.all_diff = lst_diff)
    feature(x) <- as_feature(
      lapply(lst_diff, function(x) x$gene), "受关键基因敲除而显著影响的基因"
    )
    ps.volcano <- sapply(names(lst_diff_no_filter), simplify = FALSE,
      function(name) {
        data <- lst_diff_no_filter[[name]]
        cut.fc <- if (is.null(cut.z)) 0 else cut.z
        p <- plot_volcano(
          data, "gene", use = use.p, use.fc = "Z",
          fc = cut.fc, label.fc = "Z-score", f.nudge = .5,
          mode_fc = 1, HLs = head(data$gene, n = 10), 
          show_legend = FALSE, use_break = FALSE
        )
        set_lab_legend(
          wrap(p, 5, 6),
          glue::glue("{x@sig} genes affected by {name} knockout"),
          glue::glue("受 {name} 敲除影响最显著的基因|||横坐标为标准化Z分数，Z 的绝对值越大表示该基因受 KO 扰动越显著；纵坐标为下游基因。")
        )
      })
    x <- plotsAdd(x, ps.volcano)
    snap <- .stat_table_by_pvalue(
      dplyr::bind_rows(lst_diff, .id = "ID"),
      n = 10, split = "ID", use.p = use.p, colName = "gene", 
      target = "基因", by = "被敲除后显著影响到", needSum = FALSE
    )
    x <- snapAdd(
      x, "关键基因敲除导致调控关系发生变化，{snap}"
    )
    if (!is.null(cut.z)) {
      snap.z <- glue::glue("& |Z| &gt; {cut.z}")
    } else {
      snap.z <- ""
    }
    x <- methodAdd(x, "筛选因关键基因敲低而调控关系发生显著改变的基因（阈值为{use.p} &lt; {cut.p} {snap.z}）。")
    return(x)
  })

.read_scteni_results <- function(dir, sig, type = "3000HVG",
  pattern = glue::glue("sctenifoldknk_{sig}_{type}_"), 
  ignore.case = TRUE
)
{
  files <- list.files(
    dir, pattern, full.names = TRUE, ignore.case = ignore.case
  )
  res <- lapply(files,
    function(file) {
      message(glue::glue("File exists {file.exists(file)} -> {file}"))
      readRDS(file)$diffRegulation
    })
  setNames(
    res, s(
      tools::file_path_sans_ext(basename(files)), 
      pattern, "", ignore.case = TRUE
    )
  )
}

.read_scteni_results.hb <- function(dir, type = "3000HVG",
  pattern = paste0("scTenifoldKnk_", type, "_"))
{
  .read_scteni_results(dir = dir, type = type, pattern = pattern)
}

run_remote_sctenifoldknk.huibang <- function(x,
  remote = "graphBan", remote_from = "remote")
{
  remote_dir <- glue::glue(
    "~/{s(guess_project(), '^[0-9]+_', '')}"
  )
  sig <- s(rlang::expr_text(substitute(x)), "^[^.]+\\.", "")
  dir_save <- dirname(x$file_seurat)
  dir.create(dir_save, FALSE)
  if (!file.exists(x$file_seurat)) {
    if (!is_sshfs_mount(remote_from)) {
      stop('!is_sshfs_mount(remote_from).')
    }
    message(glue::glue("Get `file_seurat` from: '{remote_from}'"))
    file_from <- file.path(remote_from, x$file_seurat)
    file.copy(file_from, x$file_seurat)
  }
  file_script <- file.path(dir_save, glue::glue("sctenifoldknk_{sig}.R"))
  x$prefix <- prefix <- glue::glue("sctenifoldknk_{sig}")
  if (TRUE) {
    lines_fun <- capture.output(dput(.run_scTenifoldKnk))
    lines_fun[1] <- paste(".run_scTenifoldKnk <- ", lines_fun[1])
    values <- bind(paste0("'", x$features, "'"))
    lines_run <- glue::glue(
      ".run_scTenifoldKnk('{basename(x$file_seurat)}', c({values}), prefix = '{prefix}')"
    )
    lines_script <- c(lines_fun, "", lines_run)
    writeLines(lines_script, file_script)
  }
  files <- paste(x$file_seurat, file_script)
  cmd_prepare <- glue::glue(
    "ssh {remote} 'mkdir {remote_dir}'"
  )
  cmd_send <- glue::glue("scp {files} {remote}:{remote_dir}")
  pg <- glue::glue(
    "/data/nas2/software/miniconda3/bin/conda run -n R4.3.3 Rscript"
  )
  cmd_run <- glue::glue("ssh {remote} 'cd {remote_dir} && nohup {pg} {basename(file_script)} > task.log 2>&1 &'")
  cdRun(cmd_prepare)
  cdRun(cmd_send)
  cdRun(cmd_run, wait = FALSE)
}

get_remote_sctenifoldknk.huibang <- function(x,
  sig = s(rlang::expr_text(substitute(x)), "^[^.]+\\.", ""),
  pattern = glue::glue("sctenifoldknk_{sig}*"),
  remote = "graphBan", remote_to = "remote", expect = length(x$features))
{
  dir_seurat <- x$dir_seurat
  dir.create(dir_seurat, FALSE)
  existFiles <- list.files(dir_seurat, pattern)
  if (!length(existFiles) || length(existFiles) < expect) {
    remote_dir <- glue::glue(
      "~/{s(guess_project(), '^[0-9]+_', '')}"
    )
    files <- pattern
    cmd_get <- glue::glue("scp {remote}:{remote_dir}/{files} {dir_seurat}")
    cdRun(cmd_get)
  }
  existFiles <- list.files(dir_seurat, pattern, full.names = TRUE)
  if (!is_sshfs_mount(remote_to)) {
    stop('!is_sshfs_mount(remote_to).')
  }
  toDir <- file.path(remote_to, x$dir_seurat)
  toDir_existFiles <- list.files(toDir, pattern)
  if (!length(toDir_existFiles) || length(toDir_existFiles) < expect) {
    message(glue::glue("Send to '{remote_to}'"))
    file.copy(existFiles, toDir)
  }
}

.run_scTenifoldKnk <- function(
  rds_file,
  genes_to_analyze,
  qc = FALSE,
  prefix = "scTenifoldKnk"
)
{
  # Check required packages
  required_pkgs <- c("Seurat", "SeuratObject", "scTenifoldKnk", "glue")

  invisible(
    lapply(required_pkgs, function(pkg) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Missing required package: ", pkg, call. = FALSE)
      }
    })
  )

  # Read Seurat object
  message(glue::glue("Reading Seurat object: {rds_file}"))
  seu <- readRDS(rds_file)

  # Extract HVG
  hvg <- Seurat::VariableFeatures(seu)

  if (length(hvg) == 0L) {
    stop("VariableFeatures(seu) is empty.")
  }

  # Extract count matrix
  counts <- SeuratObject::GetAssayData(
    object = seu,
    layer = "counts"
  )

  rn <- rownames(counts)

  # Pre-build gene pool and convert once
  genes_pool <- unique(c(genes_to_analyze, hvg))
  genes_pool <- intersect(rn, genes_pool)

  message("Convert as matrix...")
  mat0 <- as.matrix(
    counts[genes_pool, , drop = FALSE]
  )

  rn0 <- rownames(mat0)

  # Worker
  .one_gene <- function(gene)
  {
    message(glue::glue("Run with gene: {gene}"))
    if (match(gene, rn0, nomatch = 0L) == 0L) {
      message(gene, ": gene not found")
      return("Gene not found")
    }

    genes_use <- unique(c(gene, hvg))
    genes_use <- genes_use[genes_use %in% rn0]
    mat <- mat0[genes_use, , drop = FALSE]

    tryCatch(
      {
        res <- scTenifoldKnk::scTenifoldKnk(
          countMatrix = mat,
          gKO = gene,
          qc = qc
        )
        saveRDS(res, glue::glue("{prefix}_{length(hvg)}HVG_{gene}.rds"))
        message(gene, ": done")
        res
      },
      error = function(e) {
        message(gene, ": failed - ", e$message)
        e$message
      }
    )
  }
  invisible(stats::setNames(
    lapply(genes_to_analyze, .one_gene),
    genes_to_analyze
  ))
}

