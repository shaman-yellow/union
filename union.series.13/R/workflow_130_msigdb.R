# ==========================================================================
# workflow of msigdb
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_msigdb <- setClass("job_msigdb", 
  contains = c("job"),
  prototype = prototype(
    pg = "msigdb",
    info = c(""),
    cite = "",
    method = "",
    tag = "msigdb",
    analysis = "MSigDB 基因集获取"
    ))

job_msigdb <- function(mode)
{
  x <- .job_msigdb()
  x$mode <- mode
  x <- methodAdd(x, "MSigDB 整合了多种来源的基因集资源，包括 Hallmark、KEGG、Reactome、GO 及免疫相关特征基因集等，可根据研究目的筛选相应基因集并构建候选集合。所得基因集可进一步应用于 GSEA、GSVA、ssGSEA、AUCell 及通路富集分析，从而揭示样本间分子功能差异及潜在调控机制。")
  return(x)
}

setMethod("step0", signature = c(x = "job_msigdb"),
  function(x){
    step_message("Prepare your data with function `job_msigdb`.")
  })

setMethod("step1", signature = c(x = "job_msigdb"),
  function(x, pattern = NULL, name = pattern, join = TRUE, mode = x$mode,
    sub = NULL, species = "Homo sapiens")
  {
    step_message("Got data.")
    x <- .set_msig_db(x, mode, sub, species)
    sets <- as_feature(
      lapply(split(x$db_anno$gene_symbol, x$db_anno$gs_name), unique),
      glue::glue("MSigDB {name} 基因集")
    )
    if (!is.null(pattern)) {
      sets <- sets[ grp(names(sets), pattern, TRUE) ]
      methodAdd_onExit("x", "在基因集中获取与 {pattern} 相关的基因子集。")
      snap <- stat_features(sets, glue::glue("MSigDB {name}"), join, assign = "sets")
      methodAdd_onExit("x", "{snap}")
    }
    x$.feature <- sets
    return(x)
  })

.set_msig_db <- function(x, mode, sub = NULL, species = "Homo sapiens") {
  fun_data <- function(mode) {
    if (packageVersion("msigdbr") < "10.0.0") {
      db_anno <- e(msigdbr::msigdbr(species = species, category = mode))
    } else {
      db_anno <- e(msigdbr::msigdbr(species = species, collection = mode))
    }
  }
  if (length(mode) == 1L && mode != "all") {
    db_anno <- fun_data(mode)
    x <- methodAdd(
      x, "以 R 包 `msigdbr` ⟦pkgInfo('msigdbr')⟧ 获取 MSigDB 数据库 {mode} 基因集。"
    )
    if (!is.null(sub)) {
      select <- c("CP:REACTOME", "CP:KEGG", "CP:WIKIPATHWAYS")
      db_anno <- dplyr::filter(db_anno, gs_subcat %in% !!select)
      x <- methodAdd(x, "该基因集包含多个子集：{try_snap(db_anno, 'gs_subcat', 'gs_name')}。")
      x <- methodAdd(x, "选取 {bind(select)} 子集用于后续分析。")
    }
  } else {
    if (length(mode) == 1 && mode == "all") {
      mode <- c("H", paste0("C", 1:8))
    }
    db_anno <- lapply(mode, 
      function(type) {
        fun_data(type)
      })
    # db_anno <- dplyr::bind_rows(db_anno, .id = "collection")
    db_anno <- dplyr::bind_rows(db_anno)
    x <- methodAdd(
      x, "以 R 包 `msigdbr` ⟦pkgInfo('msigdbr')⟧ 获取 MSigDB 数据库 {bind(mode)} 基因集。"
    )
  }
  x$mode <- mode
  x$db_anno <- db_anno
  x$msig_db <- dplyr::select(db_anno, gs_id, symbol = gene_symbol)
  return(x)
}

setMethod("clear", signature = c(x = "job_msigdb"),
  function(x, save = FALSE, lite = TRUE, suffix = NULL, name = substitute(x, parent.frame(1)))
  {
    eval(name)
    if (save) {
      callNextMethod(
        x, save = save, lite = FALSE, suffix = suffix, name = name
      )
    }
    x$db_anno <- NULL
    x$msig_db <- NULL
    if (lite) {
      callNextMethod(
        x, save = FALSE, lite = TRUE, suffix = suffix, name = name
      )
    }
    return(x)
  })

