# ==========================================================================
# semantic extractor for workflow step
# ==========================================================================

setGeneric("description",
  function(x, ref, ...) standardGeneric("description")
)

setMethod("description", signature = c(x = "genericFunction", ref = "character"),
  function(x, ref, ...)
  {
    method <- methods::selectMethod(x, ref, ...)

    expr_body <- body(
      .get_local_fun(method)
    )

    message(glue::glue(
      "Extract semantic description from method: ",
      "{x@generic}.{ref}"
    ))

    lst_branch <- .extract_branch_semantic(
      expr_body
    )

    lst_res <- list(
      generic = x@generic,
      signature = ref,
      branches = lst_branch
    )

    return(lst_res)
  }
)

# ==========================================================================
# extract branch semantic
# ==========================================================================

.extract_branch_semantic <- function(
  expr_body,
  parent = NULL
)
{
  if (is.call(expr_body) && identical(expr_body[[1]], as.name("{"))) {

    lst_expr <- as.list(expr_body)[-1L]

  } else {

    lst_expr <- list(expr_body)
  }

  lst_branch <- list(
    branch0 = .new_semantic_branch(
      parent = parent,
      conditions = NULL,
      type = "root"
    )
  )

  obj_state <- Reduce(
    function(obj_state, expr_now)
    {
      .walk_semantic_expr(
        expr_now = expr_now,
        obj_state = obj_state
      )
    },
    lst_expr,
    init = list(
      lst_branch = lst_branch,
      n_branch = 0L,
      branch = "branch0"
    )
  )

  return(obj_state$lst_branch)
}

# ==========================================================================
# create empty branch
# ==========================================================================

.new_semantic_branch <- function(
  parent,
  conditions,
  type = "normal"
)
{
  return(list(
    parent = parent,
    type = type,
    methods = list(),
    snaps = list(),
    plots = list(),
    tables = list(),
    features = list(),
    messages = character(),
    conditions = conditions,
    terminal = FALSE
  ))
}

# ==========================================================================
# recursive semantic parser
# ==========================================================================

.walk_semantic_expr <- function(
  expr_now,
  obj_state
)
{
  if (!is.call(expr_now)) {
    return(obj_state)
  }

  lst_branch <- obj_state$lst_branch
  n_branch <- obj_state$n_branch
  branch <- obj_state$branch

  fn_name <- as.character(expr_now[[1]])

  # ------------------------------------------------------------------------
  # return semantic
  # ------------------------------------------------------------------------

  if (identical(fn_name, "return")) {

    lst_branch[[ branch ]]$terminal <- TRUE

    return(list(
      lst_branch = lst_branch,
      n_branch = n_branch,
      branch = branch
    ))
  }

  # ------------------------------------------------------------------------
  # if branch parser
  # ------------------------------------------------------------------------

  if (identical(fn_name, "if")) {

    str_condition <- paste(
      deparse(expr_now[[2]]),
      collapse = ""
    )

    # ----------------------------------------------------------------------
    # create TRUE branch
    # ----------------------------------------------------------------------

    n_branch <- n_branch + 1L

    str_branch_true <- glue::glue(
      "branch{n_branch}"
    )

    lst_branch[[ str_branch_true ]] <- .new_semantic_branch(
      parent = branch,
      conditions = str_condition,
      type = "if_true"
    )

    obj_true <- .walk_semantic_expr(
      expr_now = expr_now[[3]],
      obj_state = list(
        lst_branch = lst_branch,
        n_branch = n_branch,
        branch = str_branch_true
      )
    )

    lst_branch <- obj_true$lst_branch
    n_branch <- obj_true$n_branch

    str_branch_true_end <- obj_true$branch

    # ----------------------------------------------------------------------
    # create FALSE branch
    # ----------------------------------------------------------------------

    str_branch_false_end <- branch

    if (length(expr_now) >= 4L) {

      n_branch <- n_branch + 1L

      str_branch_false <- glue::glue(
        "branch{n_branch}"
      )

      lst_branch[[ str_branch_false ]] <- .new_semantic_branch(
        parent = branch,
        conditions = glue::glue(
          "!( {str_condition} )"
        ),
        type = "if_false"
      )

      obj_false <- .walk_semantic_expr(
        expr_now = expr_now[[4]],
        obj_state = list(
          lst_branch = lst_branch,
          n_branch = n_branch,
          branch = str_branch_false
        )
      )

      lst_branch <- obj_false$lst_branch
      n_branch <- obj_false$n_branch

      str_branch_false_end <- obj_false$branch
    }

    # ----------------------------------------------------------------------
    # create convergence branch
    # ----------------------------------------------------------------------

    n_branch <- n_branch + 1L

    str_branch_merge <- glue::glue(
      "branch{n_branch}"
    )

    lst_branch[[ str_branch_merge ]] <- .new_semantic_branch(
      parent = c(
        str_branch_true_end,
        str_branch_false_end
      ),
      conditions = NULL,
      type = "merge"
    )

    # ----------------------------------------------------------------------
    # terminal merge detection
    # ----------------------------------------------------------------------

    if (
      all(vapply(
        c(
          str_branch_true_end,
          str_branch_false_end
        ),
        function(x)
        {
          isTRUE(
            lst_branch[[ x ]]$terminal
          )
        },
        FUN.VALUE = logical(1)
      ))
    ) {

      lst_branch[[ str_branch_merge ]]$terminal <- TRUE
    }

    return(list(
      lst_branch = lst_branch,
      n_branch = n_branch,
      branch = str_branch_merge
    ))
  }

  # ------------------------------------------------------------------------
  # methodAdd semantic
  # ------------------------------------------------------------------------

  if (identical(fn_name, "methodAdd")) {

    lst_branch[[ branch ]]$methods <- c(
      lst_branch[[ branch ]]$methods,
      list(
        paste(
          deparse(expr_now),
          collapse = "\n"
        )
      )
    )
  }

  # ------------------------------------------------------------------------
  # snapAdd semantic
  # ------------------------------------------------------------------------

  if (identical(fn_name, "snapAdd")) {

    lst_branch[[ branch ]]$snaps <- c(
      lst_branch[[ branch ]]$snaps,
      list(
        paste(
          deparse(expr_now),
          collapse = "\n"
        )
      )
    )
  }

  # ------------------------------------------------------------------------
  # plotsAdd semantic
  # ------------------------------------------------------------------------

  if (identical(fn_name, "plotsAdd")) {

    vec_plot <- as.character(expr_now)[-c(1L, 2L)]

    lst_branch[[ branch ]]$plots <- unique(c(
      lst_branch[[ branch ]]$plots,
      vec_plot
    ))
  }

  # ------------------------------------------------------------------------
  # tablesAdd semantic
  # ------------------------------------------------------------------------

  if (identical(fn_name, "tablesAdd")) {

    vec_table <- as.character(expr_now)[-c(1L, 2L)]

    lst_branch[[ branch ]]$tables <- unique(c(
      lst_branch[[ branch ]]$tables,
      vec_table
    ))
  }

  # ------------------------------------------------------------------------
  # feature semantic
  # ------------------------------------------------------------------------

  if (identical(fn_name, "<-")) {

    expr_left <- expr_now[[2]]
    expr_right <- expr_now[[3]]

    str_left <- paste(
      deparse(expr_left),
      collapse = ""
    )

    if (
      grepl("^x\\$\\.feature", str_left) &&
      is.call(expr_right)
    ) {

      fn_right <- as.character(expr_right[[1]])

      if (identical(fn_right, "as_feature")) {

        str_semantic <- NULL
        str_nature <- NULL

        if (length(expr_right) >= 3L) {

          str_semantic <- paste(
            deparse(expr_right[[3]]),
            collapse = ""
          )
        }

        if ("nature" %in% names(expr_right)) {

          str_nature <- paste(
            deparse(expr_right[["nature"]]),
            collapse = ""
          )
        }

        lst_branch[[ branch ]]$features <- c(
          lst_branch[[ branch ]]$features,
          list(
            list(
              name = str_left,
              semantic = str_semantic,
              nature = str_nature
            )
          )
        )
      }
    }
  }

  # ------------------------------------------------------------------------
  # feature(x) <- semantic
  # ------------------------------------------------------------------------

  if (identical(fn_name, "feature<-")) {

    expr_right <- expr_now[[3]]

    if (
      is.call(expr_right) &&
      identical(
        as.character(expr_right[[1]]),
        "as_feature"
      )
    ) {

      str_semantic <- NULL
      str_nature <- NULL

      if (length(expr_right) >= 3L) {

        str_semantic <- paste(
          deparse(expr_right[[3]]),
          collapse = ""
        )
      }

      if ("nature" %in% names(expr_right)) {

        str_nature <- paste(
          deparse(expr_right[["nature"]]),
          collapse = ""
        )
      }

      lst_branch[[ branch ]]$features <- c(
        lst_branch[[ branch ]]$features,
        list(
          list(
            name = "feature(x)",
            semantic = str_semantic,
            nature = str_nature
          )
        )
      )
    }
  }

  # ------------------------------------------------------------------------
  # runtime message semantic
  # ------------------------------------------------------------------------

  if (identical(fn_name, "message")) {

    lst_branch[[ branch ]]$messages <- c(
      lst_branch[[ branch ]]$messages,
      paste(
        deparse(expr_now),
        collapse = ""
      )
    )
  }

  # ------------------------------------------------------------------------
  # recursive parse
  # ------------------------------------------------------------------------

  lst_child <- as.list(expr_now)[-1L]

  for (expr_child in lst_child) {

    obj_state <- .walk_semantic_expr(
      expr_now = expr_child,
      obj_state = list(
        lst_branch = lst_branch,
        n_branch = n_branch,
        branch = branch
      )
    )

    lst_branch <- obj_state$lst_branch
    n_branch <- obj_state$n_branch
    branch <- obj_state$branch

    if (
      isTRUE(
        lst_branch[[ branch ]]$terminal
      )
    ) {
      break
    }
  }

  return(list(
    lst_branch = lst_branch,
    n_branch = n_branch,
    branch = branch
  ))
}

# ==========================================================================
# helper for extracting local function
# ==========================================================================

.get_local_fun <- function(m)
{
  if (!is(m, "MethodDefinition")) {

    rlang::abort(
      "Input is not a MethodDefinition."
    )
  }

  expr_body <- body(m)

  if (
    is.call(expr_body) &&
    identical(expr_body[[1]], as.name("{")) &&
    length(expr_body) >= 2L &&
    is.call(expr_body[[2]]) &&
    identical(expr_body[[2]][[1]], as.name("<-")) &&
    identical(expr_body[[2]][[2]], as.name(".local"))
  ) {

    fun <- try(
      eval(expr_body[[2]][[3]]),
      TRUE
    )

  } else {

    fun <- try(
      m@.Data,
      TRUE
    )
  }

  if (inherits(fun, "try-error")) {

    rlang::abort(
      "Can not get local function from method."
    )
  }

  if (!is.function(fun)) {

    rlang::abort(
      "Resolved local object is not function."
    )
  }

  return(fun)
}

# ==========================================================================
# semantic tree table
# ==========================================================================

.get_semantic_tree <- function(lst_desc)
{
  data_branch <- tibble::tibble(
    branch = names(lst_desc$branches),

    parent = vapply(
      lst_desc$branches,
      function(x)
      {
        if (is.null(x$parent)) {
          return(NA_character_)
        }

        paste(
          x$parent,
          collapse = " | "
        )
      },
      FUN.VALUE = character(1)
    ),

    type = vapply(
      lst_desc$branches,
      function(x)
      {
        x$type
      },
      FUN.VALUE = character(1)
    ),

    terminal = vapply(
      lst_desc$branches,
      function(x)
      {
        isTRUE(x$terminal)
      },
      FUN.VALUE = logical(1)
    ),

    conditions = vapply(
      lst_desc$branches,
      function(x)
      {
        if (is.null(x$conditions)) {
          return(NA_character_)
        }

        x$conditions
      },
      FUN.VALUE = character(1)
    ),

    n_methods = vapply(
      lst_desc$branches,
      function(x)
      {
        length(x$methods)
      },
      FUN.VALUE = integer(1)
    ),

    n_snaps = vapply(
      lst_desc$branches,
      function(x)
      {
        length(x$snaps)
      },
      FUN.VALUE = integer(1)
    ),

    n_plots = vapply(
      lst_desc$branches,
      function(x)
      {
        length(x$plots)
      },
      FUN.VALUE = integer(1)
    ),

    n_tables = vapply(
      lst_desc$branches,
      function(x)
      {
        length(x$tables)
      },
      FUN.VALUE = integer(1)
    ),

    n_features = vapply(
      lst_desc$branches,
      function(x)
      {
        length(x$features)
      },
      FUN.VALUE = integer(1)
    ),

    n_messages = vapply(
      lst_desc$branches,
      function(x)
      {
        length(x$messages)
      },
      FUN.VALUE = integer(1)
    )
  )

  # ------------------------------------------------------------------------
  # semantic score
  # ------------------------------------------------------------------------

  data_branch$score_semantic <-
    data_branch$n_features * 10L +
    data_branch$n_methods * 4L +
    data_branch$n_plots * 2L +
    data_branch$n_tables * 2L +
    data_branch$n_snaps * 1L

  # ------------------------------------------------------------------------
  # runtime branch detection
  # ------------------------------------------------------------------------

  vec_runtime_pattern <- c(
    "workers",
    "future",
    "multicore",
    "remote",
    "parallel",
    "threads",
    "cores",
    "progress",
    "verbose"
  )

  data_branch$is_runtime <- vapply(
    data_branch$conditions,
    function(x)
    {
      if (is.na(x)) {
        return(FALSE)
      }

      any(grepl(
        paste(vec_runtime_pattern, collapse = "|"),
        x,
        ignore.case = TRUE
      ))
    },
    FUN.VALUE = logical(1)
  )

  # ------------------------------------------------------------------------
  # semantic signature
  # ------------------------------------------------------------------------

  data_branch$signature <- vapply(
    seq_len(nrow(data_branch)),
    function(i)
    {
      x <- lst_desc$branches[[ i ]]

      paste(
        c(
          paste0("M:", x$methods),
          paste0("S:", x$snaps),
          paste0("P:", x$plots),
          paste0("T:", x$tables),
          vapply(
            x$features,
            function(y)
            {
              paste0(
                "F:",
                y$name,
                "::",
                y$semantic
              )
            },
            FUN.VALUE = character(1)
          )
        ),
        collapse = " || "
      )
    },
    FUN.VALUE = character(1)
  )

  # ------------------------------------------------------------------------
  # semantic noise detection
  # ------------------------------------------------------------------------

  data_branch$is_semantic_noise <- FALSE

  idx_noise <- which(
    data_branch$score_semantic == 0L &
    !data_branch$terminal &
    data_branch$type != "merge"
  )

  data_branch$is_semantic_noise[idx_noise] <- TRUE

  idx_runtime <- which(
    data_branch$is_runtime &
    data_branch$score_semantic == 0L
  )

  data_branch$is_semantic_noise[idx_runtime] <- TRUE

  # ------------------------------------------------------------------------
  # importance
  # ------------------------------------------------------------------------

  data_branch$importance <- dplyr::case_when(
    data_branch$n_features > 0L ~ "critical",
    data_branch$n_methods > 0L &
      data_branch$n_plots > 0L ~ "major",
    data_branch$score_semantic > 0L ~ "minor",
    TRUE ~ "noise"
  )

  # ------------------------------------------------------------------------
  # semantic type
  # ------------------------------------------------------------------------

  data_branch$semantic_type <- dplyr::case_when(
    data_branch$is_runtime ~ "runtime",
    data_branch$terminal ~ "terminal",
    data_branch$n_features > 0L ~ "feature",
    data_branch$n_plots > 0L ~ "visualization",
    data_branch$n_methods > 0L ~ "analysis",
    data_branch$type == "merge" ~ "merge",
    TRUE ~ "other"
  )

  data_branch <- dplyr::relocate(
    data_branch,
    branch,
    parent,
    type,
    terminal,
    score_semantic,
    conditions
  )

  return(data_branch)
}
