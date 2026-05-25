# ==========================================================================
# workflow of ssgsea
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_ssgsea <- setClass("job_ssgsea", 
  contains = c("job"),
  prototype = prototype(
    pg = "ssgsea",
    info = c("https://www.bioconductor.org/packages/release/bioc/html/GSVA.html"),
    cite = "[@Systematic_RNA_Barbie_2009]",
    method = "",
    tag = "ssgsea",
    analysis = "ssGSEA 单样本GSEA富集分析"
    ))

setGeneric("asjob_ssgsea",
   function(x, ...) standardGeneric("asjob_ssgsea"))

setMethod("asjob_ssgsea", signature = c(x = "job_limma"),
  function(x, use.filter = NULL, use = .guess_symbol(x), 
    use.format = TRUE, ...)
  {
    if (x@step < 1L) {
      stop('x@step < 1L.')
    }
    cli::cli_alert_info("extract_unique_genes.job_limma")
    if (FALSE) {
      object <- extract_unique_genes.job_limma(
        x, use.filter, use, use.format = use.format, ...
      )
    } else {
      object <- x$normed_data
      if (is(object, "DGEList")) {
        object <- new_from_package(
          "EList", "limma", list(E = object$counts, targets = object$samples, genes = object$genes)
        )
      }
    }
    mtx <- object$E
    genes <- object$genes
    if (use.format) {
      rownames(mtx) <- gname(genes[[use]])
    } else {
      rownames(mtx) <- genes[[use]]
    }
    metadata <- object$targets
    if (x@step >= 2L) {
      contrasts <- list(.guess_compare_limma(x, 1L))
    } else {
      contrasts <- NULL
      message("Can not match 'contrasts' in `x`.")
    }
    # mtx <- mtx[!is.na(rownames(mtx)), ]
    x <- job_ssgsea(mtx)
    x$metadata <- metadata
    x$contrasts <- contrasts
    return(x)
  })

job_ssgsea <- function(x)
{
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    BiocManager::install("GSVA")
  }
  if (!is(x, "GsvaExprData")) {
    stop('!is(x, "GsvaExprData").')
  }
  .job_ssgsea(object = x)
}

setMethod("step0", signature = c(x = "job_ssgsea"),
  function(x){
    step_message("Prepare your data with function `job_ssgsea`.")
  })

as_collection <- function(sets) {
  if (is.null(names(sets))) {
    stop('is.null(names(sets)).')
  }
  if (is(sets, "feature")) {
    sets <- setNames(sets@.Data, names(sets))
  }
  if (!is(sets, "list") && all(vapply(sets, is.character, logical(1L)))) {
    stop('!is(sets, "list") && all(vapply(sets, is.character, logical(1L))).')
  }
  sets <- mapply(names(sets), sets, SIMPLIFY = FALSE,
    FUN = function(name, genes) {
      GSEABase::GeneSet(genes, setName = name)
    })
  e(GSEABase::GeneSetCollection(sets))
}

setMethod("step1", signature = c(x = "job_ssgsea"),
  function(x, mode = c("matrisome"), org = c("human", "mouse"), sets)
  {
    step_message("Calculate ssGSEA enrichment score.")
    x <- methodAdd(x, "以 R 包 `GSVA` ⟦pkgInfo('GSVA')⟧ 用于 ssGSEA 分析。")
    if (missing(sets)) {
      mode <- match.arg(mode)
      if (mode == "matrisome") {
        db <- .job_matrisome(sig = x@sig)
        db <- step1(db, org)
        x <- snapAdd(x, db)
        # extract results.
        x$db <- db$db
        sets <- x$.feature <- feature(db)
        x <- snapAdd(x, "将{snap(sets)}用于 ssGSEA 富集分析。")
        sets <- list(matrisome = unlist(sets, use.names = FALSE))
      }
    } else {
      if (is(sets, "feature_char")) {
        x <- snapAdd(x, "将{snap(sets)}用于 ssGSEA 富集分析。")
        sets <- setNames(list(resolve_feature(sets)), sets@snap)
      }
      if (is(sets, "feature_list")) {
        x <- snapAdd(x, "将{snap(sets)}用于 ssGSEA 富集分析。")
        sets <- sets@.Data
        if (any(vapply(sets, class, "") != "character")) {
          stop("`feature_list` should not be nest list.")
        }
      }
    }
    if (is(sets, "list")) {
      sets <- as_collection(sets)
    }
    fun_compute <- function(...) {
      param <- e(GSVA::ssgseaParam(object(x), sets))
      res <- e(GSVA::gsva(param))
    }
    res <- expect_local_data(
      "tmp", "ssgsea", fun_compute,
      list(colnames(object(x)), rownames(object(x)), sets, object(x)[[1L]])
    )
    data <- dplyr::select(x$metadata, group, sample)
    if (!identical(data$sample, colnames(res))) {
      stop('!identical(data$sample, colnames(res)).')
    }
    data <- cbind(
      data, setNames(data.frame(t(res), check.names = FALSE), names(sets))
    )
    x$data <- tidyr::pivot_longer(data, dplyr::all_of(names(sets)), names_to = "type", values_to = "score")
    if (!is.null(x$contrasts)) {
      p.scores <- lapply(x$contrasts, 
        function(group) {
          data <- dplyr::filter(x$data, group %in% !!group)
          data <- dplyr::mutate(
            data, group = factor(group, levels = !!group)
          )
          p <- .map_boxplot2(
            x$data, TRUE, y = "score", ylab = "Enrichment score", ids = "type"
          )
          p <- set_lab_legend(
            wrap(p, 2.5, 3),
            glue::glue("{x@sig} {bind(group, co = ' ')} boxplot of enrichment score"),
            glue::glue("{bind(group, co = ' ')} 富集评分箱形图。")
          )
          snap(p) <- .stat_compare_by_pvalue(
            p, group, "富集", mode = "enrichment"
          )
          p
        })
      if (length(p.scores) == 1L) {
        x <- snapAdd(x, snap(p.scores[[1]]))
      }
      names(p.scores) <- vapply(
        x$contrasts, bind, character(1), co = " "
      )
      x <- plotsAdd(x, p.scores = p.scores)
    }
    return(x)
  })
