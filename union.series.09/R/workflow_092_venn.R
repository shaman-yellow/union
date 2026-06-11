# ==========================================================================
# workflow of venn
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

.job_venn <- setClass("job_venn", 
  contains = c("job"),
  prototype = prototype(
    pg = "venn",
    info = c("..."),
    cite = "",
    method = "",
    tag = "venn",
    analysis = "Venn 交集"
    ))

.job_vennDEGs <- setClass("job_vennDEGs", contains = c("job_venn"))

setGeneric("asjob_venn",
  function(x, ...) standardGeneric("asjob_venn"))

job_vennDEGs <- function(pattern, exclude = NULL,
  mode = c("raw", "pro"), name = "guess", ...,
  pattern_dataset = ".*(GSE[0-9]+).*", gp = NULL)
{
  mode <- match.arg(mode)
  if (mode == "raw") {
    fun_extract <- function(x) x@tables$step2$tops[[1]][, 1:3]
    snapAdd_onExit("x", "数据集的差异表达基因的交集。")
  } else {
    stop(glue::glue("..."))
  }
  degs <- collate_dataset_DEGs(
    pattern, name = name, exclude = exclude, fun_extract = fun_extract, ...
  )
  degs_versus <- collate(
    pattern, function(x) names(x@tables$step2$tops)[1], exclude, ...
  )
  projects <- collate(
    pattern, function(x) x$project, exclude, ...
  )
  names(degs_versus) <- s(names(degs_versus), pattern_dataset, "\\1")
  degs <- split(degs$symbol, degs$Dataset)
  names(degs) <- s(names(degs), pattern_dataset, "\\1")
  x <- .job_vennDEGs(job_venn(lst = degs))
  x$degs_versus <- degs_versus
  if (!is.null(gp)) {
    metadata <- as_df.lst(x$degs_versus, "project", "versus")
    groups <- group_strings(unlist(x$degs_versus), gp, "versus")
    x$metadata <- map(
      metadata, "versus", groups, "versus", "group", col = "group"
    )
    x$metadata <- tibble::as_tibble(x$metadata)
    x$metadata <- set_lab_legend(
      x$metadata, "metadata of mutiple datasets", "数据集的元数据信息 (分组信息) "
    )
    x$groups <- nl(x$metadata$project, x$metadata$group)
  }
  x$projects <- projects
  return(x)
}

job_venn <- function(..., mode = c("key", "candidates", "ck", "other"),
  analysis = NULL, lst = NULL, name = NULL, fun_map = function(x) x)
{
  if (is.null(lst)) {
    object <- list(...)
  } else {
    object <- lst
  }
  nature <- "基因集"
  if (all(vapply(object, is, logical(1), "feature"))) {
    snaps <- paste0(
      "- ", vapply(
        object, snap, character(1), enumerate = FALSE, unlist = TRUE
      )
    )
    methodAdd_onExit("x", "数据集为：\n\n{bind(snaps, co = '\n')}\n\n\n\n")
    nature <- object[[1]]@nature
    message(glue::glue("Use nature as: {nature}"))
    object <- lapply(object, function(x) fun_map(unlist(x@.Data)))
    if (is.null(analysis)) {
      if (is.null(name)) {
        if (missing(mode)) {
          warning(crayon::red("`mode` is missing, use default: 'key'"))
        }
        mode <- match.arg(mode)
        name <- switch(
          mode, key = "关键", candidates = "候选", ck = "候选关键"
        )
      }
      analysis <- glue::glue("{name}{nature}")
    }
  } else {
    message(glue::glue("Some were not 'feature', be carefull! If you need automatic snap ..."))
  }
  x <- .job_venn(object = lapply(object, unlist))
  x <- methodAdd(x, "以 R 包 `ggVennDiagram` ⟦pkgInfo('ggVennDiagram')⟧ 对{nature}取交集。")
  x$nature <- nature
  x$analysis <- analysis
  return(x)
}

setMethod("step0", signature = c(x = "job_venn"),
  function(x){
    step_message("Prepare your data with function `job_venn`.")
  })

setMethod("step1", signature = c(x = "job_venn"),
  function(x, ...){
    step_message("Intersection.")
    p.venn <- new_venn(lst = object(x), force_upset = FALSE, ...)
    p.venn <- set_lab_legend(
      p.venn,
      glue::glue("{x@sig} intersection of {bind(names(object(x)), co = ' with ')}"),
      glue::glue("{bind(names(object(x)))} 交集维恩图|||不同颜色圆圈代表不同数据集，中间重叠部分表示同时存在多个集合中。图中 {length(p.venn$ins)} 交集为：{less(p.venn$ins, 20)}。")
    )
    x$.append_heading <- FALSE
    if (identical(parent.frame(1), .GlobalEnv)) {
      job_append_heading(
        x, heading = glue::glue(
          "汇总: ", bind(names(object(x)), co = " + ")
        )
      )
    }
    if (length(p.venn$ins) < 10) {
      iter <- glue::glue(" ({bind(p.venn$ins)}) ")
    } else {
      iter <- ""
    }
    x <- snapAdd(x, "对{bind(names(object(x)))} 取交集，得到{length(p.venn$ins)}个交集{iter}{aref(p.venn)}。")
    x <- plotsAdd(x, p.venn)
    x$.feature_sets <- as_feature(
      lapply(object(x), unique), "All sets for Venn"
    )
    if (!is.null(x$analysis)) {
      feature(x) <- as_feature(p.venn$ins, x$analysis, nature = x$nature, ...)
    } else {
      x$.feature <- as_feature(p.venn$ins, x, nature = x$nature, ...)
    }
    return(x)
  })

new_venn <- function(..., lst = NULL, wrap = TRUE,
  fun_pre = rm.no, force_upset = NULL, n = NULL,
  venn_palette = "pastel_pathway",
  venn_alpha = 0.68,
  venn_overlap_lighten = 0.28,
  venn_edge_color = "white",
  venn_edge_size = 0.35,
  venn_name_position = c("outside", "inside", "both", "none"),
  venn_name_outside_mode = c("top", "radial"),
  venn_name_push_x = 1,
  venn_name_push_y = 1,
  venn_name_outside_offset = 0.06,
  venn_name_wrap = 30L,
  venn_name_size = 4.5,
  venn_name_color = "#333333",
  venn_label_fill = NA,
  venn_label_border_size = 0
)
{
  venn_name_position <- match.arg(venn_name_position)
  venn_name_outside_mode <- match.arg(venn_name_outside_mode)

  .get_palette <- function(palette)
  {
    lst_palette <- list(
      pastel_pathway = c("#b8d8ad", "#f8f4eb", "#ead8df"),
      red_soft = c("#ce7e73", "#f5a596", "#f8d5ce"),
      blue_pink = c("#38546d", "#cba3b2", "#fbbe85"),
      wine_orange = c("#782c1f", "#df8d15", "#f2bc94"),
      blue_blush = c("#526fb3", "#f3d0cf", "#e3abac"),
      coral = c("#eb6d82", "#feb4b5", "#fed9be"),
      blue_mint = c("#8390ca", "#c0e0db", "#f7a7a6"),
      magenta = c("#7d1339", "#ee3862", "#fba8ad"),
      brown = c("#6a5853", "#e88e80", "#f1c5bf"),
      green_pink = c("#3ea88a", "#d4d2bd", "#eda4c8")
    )

    if (length(palette) > 1L) {
      return(as.character(palette))
    }

    if (!palette %in% names(lst_palette)) {
      message(glue::glue("Unknown Venn palette: {palette}. Fallback to pastel_pathway."))
      return(lst_palette[["pastel_pathway"]])
    }

    lst_palette[[palette]]
  }

  .mix_color <- function(vec_color)
  {
    if (length(vec_color) == 1L) {
      return(vec_color)
    }

    mat_rgb <- grDevices::col2rgb(vec_color) / 255
    vec_rgb <- rowMeans(mat_rgb)
    vec_rgb <- vec_rgb + (1 - vec_rgb) * venn_overlap_lighten

    grDevices::rgb(vec_rgb[1L], vec_rgb[2L], vec_rgb[3L])
  }

  .get_region_members <- function(id, vec_set_name)
  {
    id <- as.character(id)

    if (is.na(id) || !nzchar(id)) {
      return(character())
    }

    vec_token <- unlist(strsplit(id, "\\s*/\\s*|\\s*&\\s*|\\s*,\\s*"))
    vec_token <- trimws(vec_token)
    vec_token <- vec_token[nzchar(vec_token)]

    vec_member <- intersect(vec_token, vec_set_name)

    if (length(vec_member)) {
      return(vec_member)
    }

    if (grepl("^[01]+$", id) && nchar(id) == length(vec_set_name)) {
      vec_bit <- strsplit(id, "")[[1L]]
      return(vec_set_name[vec_bit == "1"])
    }

    vec_index <- suppressWarnings(as.integer(vec_token))
    vec_index <- vec_index[!is.na(vec_index)]
    vec_index <- vec_index[vec_index >= 1L & vec_index <= length(vec_set_name)]

    if (length(vec_index)) {
      return(vec_set_name[vec_index])
    }

    character()
  }

  .style_venn_fill <- function(p, lst)
  {
    vec_layer <- which(vapply(seq_along(p$layers), function(i) {
      data_layer <- p$layers[[i]]$data

      if (inherits(data_layer, "waiver") || is.null(data_layer)) {
        data_layer <- p$data
      }

      if (is.null(names(data_layer))) {
        return(FALSE)
      }

      has_region_data <- "count" %in% names(data_layer) &&
        any(c("id", "name") %in% names(data_layer))

      has_region_geom <- inherits(p$layers[[i]]$geom, "GeomSf") ||
        inherits(p$layers[[i]]$geom, "GeomPolygon")

      has_region_data && has_region_geom
    }, logical(1L)))

    if (!length(vec_layer)) {
      message("Skip Venn fill style: no valid region layer was found.")
      return(p)
    }

    n_layer <- vec_layer[1L]
    data_layer <- p$layers[[n_layer]]$data

    if (inherits(data_layer, "waiver") || is.null(data_layer)) {
      data_layer <- p$data
    }

    vec_palette <- .get_palette(venn_palette)

    if (length(vec_palette) < length(lst)) {
      vec_palette <- grDevices::colorRampPalette(vec_palette)(length(lst))
    }

    vec_palette <- vec_palette[seq_len(length(lst))]
    names(vec_palette) <- names(lst)

    vec_id <- if ("id" %in% names(data_layer)) {
      as.character(data_layer$id)
    } else {
      as.character(data_layer$name)
    }

    data_layer$venn_fill <- vapply(vec_id, function(id) {
      vec_member <- .get_region_members(id, names(lst))

      if (!length(vec_member)) {
        return("#f2f2f2")
      }

      .mix_color(vec_palette[vec_member])
    }, character(1L))

    p$layers[[n_layer]]$data <- data_layer
    p$layers[[n_layer]]$mapping$fill <- rlang::quo(venn_fill)
    p$layers[[n_layer]]$aes_params$alpha <- venn_alpha
    p$layers[[n_layer]]$aes_params$colour <- venn_edge_color
    p$layers[[n_layer]]$aes_params$linewidth <- venn_edge_size
    p$layers[[n_layer]]$aes_params$size <- venn_edge_size

    p + ggplot2::scale_fill_identity(guide = "none")
  }

  .style_venn_text <- function(p)
  {
    .clean_text <- function(x)
    {
      x <- as.character(x)
      x <- gsub("\\s+", " ", x)
      trimws(x)
    }

    .get_xy_data <- function(data_layer)
    {
      if (is.null(data_layer) || inherits(data_layer, "waiver")) {
        return(NULL)
      }

      if (all(c("x", "y") %in% names(data_layer))) {
        return(data.frame(
            x = as.numeric(data_layer$x),
            y = as.numeric(data_layer$y)
            ))
      }

      if (all(c("X", "Y") %in% names(data_layer))) {
        return(data.frame(
            x = as.numeric(data_layer$X),
            y = as.numeric(data_layer$Y)
            ))
      }

      if ("geometry" %in% names(data_layer) &&
        requireNamespace("sf", quietly = TRUE)) {
        mat_coord <- tryCatch(
          sf::st_coordinates(data_layer$geometry),
          error = function(e) NULL
        )

        if (!is.null(mat_coord) && all(c("X", "Y") %in% colnames(mat_coord))) {
          return(data.frame(
              x = as.numeric(mat_coord[, "X"]),
              y = as.numeric(mat_coord[, "Y"])
              ))
        }
      }

      NULL
    }

    .get_union_xy <- function(p_base)
    {
      lst_xy <- lapply(seq_along(p_base$layers), function(i) {
        if (is(p_base$layers[[i]]$geom, "GeomText") ||
          is(p_base$layers[[i]]$geom, "GeomLabel")) {
          return(NULL)
        }

        data_layer <- p_base$layers[[i]]$data

        if (is.null(data_layer) || inherits(data_layer, "waiver")) {
          data_layer <- p_base$data
        }

        data_xy <- .get_xy_data(data_layer)

        if (is.null(data_xy) || nrow(data_xy) < 6L) {
          return(NULL)
        }

        data_xy
        })

      lst_xy <- Filter(Negate(is.null), lst_xy)

      if (!length(lst_xy)) {
        data_build <- tryCatch(
          ggplot2::ggplot_build(p_base),
          error = function(e) NULL
        )

        if (is.null(data_build)) {
          return(NULL)
        }

        lst_xy <- lapply(data_build$data, function(data_layer) {
          if ("label" %in% names(data_layer)) {
            return(NULL)
          }

          data_xy <- .get_xy_data(data_layer)

          if (is.null(data_xy) || nrow(data_xy) < 6L) {
            return(NULL)
          }

          data_xy
        })

        lst_xy <- Filter(Negate(is.null), lst_xy)
      }

      if (!length(lst_xy)) {
        return(NULL)
      }

      data_xy <- dplyr::bind_rows(lst_xy)
      data_xy <- data_xy[is.finite(data_xy$x) & is.finite(data_xy$y), , drop = FALSE]

      if (!nrow(data_xy)) {
        return(NULL)
      }

      data_xy
    }

    .get_set_name_layers <- function(p_base)
    {
      vec_set_name <- .clean_text(names(lst))

      which(vapply(seq_along(p_base$layers), function(i) {
          if (!is(p_base$layers[[i]]$geom, "GeomText")) {
            return(FALSE)
          }

          data_layer <- tryCatch(
            ggplot2::layer_data(p_base, i),
            error = function(e) NULL
          )

          if (is.null(data_layer)) {
            return(FALSE)
          }

          if ("label" %in% names(data_layer)) {
            vec_label <- .clean_text(data_layer$label)
          } else if ("name" %in% names(data_layer)) {
            vec_label <- .clean_text(data_layer$name)
          } else {
            return(FALSE)
          }

          any(vec_label %in% vec_set_name)
          }, logical(1L)))
    }

    .get_single_set_anchor <- function(data_label)
    {
      data_xy <- .get_xy_data(data_label)

      if (is.null(data_xy)) {
        return(NULL)
      }

      vec_id <- if ("id" %in% names(data_label)) {
        as.character(data_label$id)
      } else if ("name" %in% names(data_label)) {
        as.character(data_label$name)
      } else {
        rep(NA_character_, nrow(data_label))
      }

      lst_member <- lapply(vec_id, .get_region_members, vec_set_name = names(lst))
      vec_n_member <- vapply(lst_member, length, integer(1L))
      keep <- vec_n_member == 1L

      if (!any(keep) && "name" %in% names(data_label)) {
        vec_name <- .clean_text(data_label$name)
        keep <- vec_name %in% .clean_text(names(lst))
        lst_member <- lapply(vec_name, function(x) {
          n_match <- match(x, .clean_text(names(lst)))

          if (is.na(n_match)) {
            return(character())
          }

          names(lst)[n_match]
        })
      }

      if (!any(keep)) {
        return(NULL)
      }

      data_anchor <- data.frame(
        name = vapply(lst_member[keep], function(x) x[1L], character(1L)),
        x = data_xy$x[keep],
        y = data_xy$y[keep]
      )

      data_anchor <- data_anchor[!duplicated(data_anchor$name), , drop = FALSE]
      data_anchor
    }

    .get_fallback_anchor <- function(data_union)
    {
      n_set <- length(lst)
      x_range <- range(data_union$x, na.rm = TRUE)
      y_range <- range(data_union$y, na.rm = TRUE)
      x_mid <- mean(x_range)
      y_mid <- mean(y_range)
      x_span <- diff(x_range)
      y_span <- diff(y_range)

      if (n_set == 2L) {
        vec_angle <- c(pi / 2, -pi / 2)
      } else if (n_set == 3L) {
        vec_angle <- c(5 * pi / 6, pi / 6, -pi / 2)
      } else if (n_set == 4L) {
        vec_angle <- c(3 * pi / 4, pi / 4, -3 * pi / 4, -pi / 4)
      } else {
        vec_angle <- seq(
          from = pi / 2,
          to = pi / 2 - 2 * pi * (n_set - 1L) / n_set,
          length.out = n_set
        )
      }

      data.frame(
        name = names(lst),
        x = x_mid + cos(vec_angle) * x_span * 0.35,
        y = y_mid + sin(vec_angle) * y_span * 0.35
      )
    }

    .complete_anchor <- function(data_anchor, data_union)
    {
      data_fallback <- .get_fallback_anchor(data_union)

      if (is.null(data_anchor) || !nrow(data_anchor)) {
        return(data_fallback)
      }

      data_out <- data_fallback
      n_match <- match(data_out$name, data_anchor$name)
      keep <- !is.na(n_match)

      data_out$x[keep] <- data_anchor$x[n_match[keep]]
      data_out$y[keep] <- data_anchor$y[n_match[keep]]

      data_out
    }

    .make_outside_name_data <- function(data_anchor, data_union)
    {
      data_anchor <- .complete_anchor(data_anchor, data_union)

      x_range <- range(data_union$x, na.rm = TRUE)
      y_range <- range(data_union$y, na.rm = TRUE)
      x_mid <- mean(x_range)
      y_mid <- mean(y_range)
      max_span <- max(diff(x_range), diff(y_range))

      data_text <- do.call(rbind, lapply(seq_len(nrow(data_anchor)), function(i) {
          set_name <- data_anchor$name[i]
          dx <- data_anchor$x[i] - x_mid
          dy <- data_anchor$y[i] - y_mid
          dist <- sqrt(dx ^ 2 + dy ^ 2)

          if (!is.finite(dist) || dist < .Machine$double.eps) {
            n_index <- match(set_name, names(lst))
            angle <- pi / 2 - 2 * pi * (n_index - 1L) / length(lst)
            dx <- cos(angle)
            dy <- sin(angle)
            dist <- 1
          }

          ux <- dx / dist
          uy <- dy / dist

          vec_proj <- (data_union$x - x_mid) * ux + (data_union$y - y_mid) * uy
          boundary_proj <- max(vec_proj, na.rm = TRUE)

          if (venn_name_outside_mode == "top") {
            if (uy > 0.2) {
              x_label <- data_anchor$x[i] + dx * 0.05
              y_label <- y_range[2L] + max_span * venn_name_outside_offset * venn_name_push_y

              hjust <- if (ux > 0.15) {
                0
              } else if (ux < -0.15) {
                1
              } else {
                0.5
              }

              vjust <- 0
            } else if (uy < -0.2) {
              x_label <- x_mid
              y_label <- y_range[1L] - max_span * venn_name_outside_offset * venn_name_push_y
              hjust <- 0.5
              vjust <- 1
            } else {
              x_label <- x_mid + ux * boundary_proj +
                ux * max_span * venn_name_outside_offset * venn_name_push_x

              y_label <- y_mid + uy * boundary_proj +
                uy * max_span * venn_name_outside_offset * venn_name_push_y

              hjust <- if (ux > 0.15) {
                0
              } else if (ux < -0.15) {
                1
              } else {
                0.5
              }

              vjust <- 0.5
            }
          } else {
            x_label <- x_mid + ux * boundary_proj +
              ux * max_span * venn_name_outside_offset * venn_name_push_x

            y_label <- y_mid + uy * boundary_proj +
              uy * max_span * venn_name_outside_offset * venn_name_push_y

            hjust <- if (ux > 0.15) {
              0
            } else if (ux < -0.15) {
              1
            } else {
              0.5
            }

            vjust <- if (uy > 0.15) {
              0
            } else if (uy < -0.15) {
              1
            } else {
              0.5
            }
          }
          data.frame(
            name = set_name,
            label = stringr::str_wrap(set_name, venn_name_wrap),
            x = x_label,
            y = y_label,
            hjust = hjust,
            vjust = vjust
          )
          }))

      data_text
    }

    which <- NA_integer_

    for (i in rev(seq_along(p$layers))) {
      if (is(p$layers[[i]]$geom, "GeomLabel")) {
        which <- i
        break
      }
    }

    data_anchor <- NULL

    if (!is.na(which)) {
      data_label <- p$layers[[which]]$data
      data_anchor <- .get_single_set_anchor(data_label)

      vec_id <- if ("id" %in% names(data_label)) {
        as.character(data_label$id)
      } else {
        as.character(data_label$name)
      }

      lst_member <- lapply(vec_id, .get_region_members, vec_set_name = names(lst))
      is_single_set <- vapply(lst_member, length, integer(1L)) == 1L

      if (!any(is_single_set)) {
        is_single_set <- !grepl("/", vec_id, fixed = TRUE)
      }

      if (venn_name_position %in% c("inside", "both")) {
        data_label[is_single_set, ] <- dplyr::mutate(
          data_label[is_single_set, ],
          both = paste0(
            stringr::str_wrap(name, venn_name_wrap), "\n",
            count, " (", percent, ")"
          )
        )

        data_label[!is_single_set, ] <- dplyr::mutate(
          data_label[!is_single_set, ],
          both = paste0(count, " (", percent, ")")
        )
      } else {
        data_label <- dplyr::mutate(
          data_label,
          both = paste0(count, " (", percent, ")")
        )
      }

      p$layers[[which]]$data <- data_label
      p$layers[[which]]$aes_params$fill <- venn_label_fill
      p$layers[[which]]$aes_params$alpha <- 1
      p$layers[[which]]$geom_params$label.size <- venn_label_border_size
    }

    vec_name_layer <- .get_set_name_layers(p)

    if (length(vec_name_layer)) {
      p$layers[vec_name_layer] <- NULL
    }

    if (venn_name_position %in% c("outside", "both")) {
      data_union <- .get_union_xy(p)

      if (is.null(data_union)) {
        message("Skip outside Venn names: failed to extract Venn geometry.")
        return(p)
      }

      data_text <- .make_outside_name_data(data_anchor, data_union)

      if (!is.null(data_text) && nrow(data_text)) {
        x_range <- range(c(data_union$x, data_text$x), na.rm = TRUE)
        y_range <- range(c(data_union$y, data_text$y), na.rm = TRUE)
        x_pad <- diff(x_range) * 0.10
        y_pad <- diff(y_range) * 0.10

        p <- p +
          ggplot2::coord_equal(
            xlim = c(x_range[1L] - x_pad, x_range[2L] + x_pad),
            ylim = c(y_range[1L] - y_pad, y_range[2L] + y_pad),
            clip = "off"
            ) +
          ggplot2::geom_text(
            data = data_text,
            mapping = ggplot2::aes(
              x = x,
              y = y,
              label = label,
              hjust = hjust,
              vjust = vjust
              ),
            inherit.aes = FALSE,
            size = venn_name_size,
            colour = venn_name_color,
            show.legend = FALSE
            ) +
          ggplot2::theme(
            plot.margin = ggplot2::margin(8, 16, 8, 16)
          )
      }
    }

    p
  }

  if (!is.null(lst) && length(list(...))) {
    lst <- c(lst, list(...))
  }

  if (is.null(lst)) {
    lst <- list(...)
  }

  lst <- lst_clear0(lst)
  lst <- lapply(lst, function(x) as.character(fun_pre(x)))

  if (is.null(force_upset)) {
    if (length(lst) > 3L) {
      force_upset <- TRUE
    } else {
      force_upset <- FALSE
    }
  }

  if (force_upset) {
    p <- ggVennDiagram::ggVennDiagram(
      lst, force_upset = TRUE, nintersects = n
    )

    p$plotlist[[2L]]$layers[[1L]]$geom_params$width <- 0.4
    p$plotlist[[3L]]$layers[[1L]]$geom_params$width <- 0.7
  } else {
    p <- ggVennDiagram::ggVennDiagram(
      lst,
      label_percent_digit = 1L,
      set_color = venn_edge_color,
      edge_size = 0.3
    )

    p <- .style_venn_fill(p, lst)

    p <- p +
      ggplot2::theme_void() +
      ggplot2::guides(fill = "none") +
      ggplot2::theme(
        axis.text = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank()
      )

    p <- .style_venn_text(p)
  }

  if (wrap) {
    p <- wrap(p, 5, 6)
  }

  attr(p, "ins") <- ins <- ins(lst = lst)
  attr(p, "lich") <- new_lich(list(All_intersection = ins))
  lab(p) <- paste0("Intersection of ", paste0(names(lst), collapse = " with "))
  p <- setLegend(p, "为 {bind(names(lst))} 各自交集。")
  p
}

setMethod("meta", signature = c(x = "job_vennDEGs"),
  function(x, group = NULL, bind = TRUE, arrange = TRUE, get = "project")
  {
    if (is.null(x$metadata)) {
      stop('is.null(x$metadata).')
    }
    if (arrange) {
      x$metadata <- dplyr::arrange(x$metadata, group)
    }
    if (!is.null(group)) {
      res <- dplyr::filter(x$metadata, group %in% !!group)
      if (!nrow(res)) {
        stop('!nrow(res). No any results.')
      }
      res <- res[[ get ]]
    } else {
      res <- x$metadata[[ get ]]
    }
    if (bind) {
      bind(res)
    } else {
      res
    }
  })

setMethod("feature", signature = c(x = "job_vennDEGs"),
  function(x, group, intersect = TRUE){
    if (missing(group)) {
      callNextMethod(x)
    } else {
      projects <- meta(x, group, bind = FALSE)
      sets <- object(x)[ names(object(x)) %in% projects ]
      if (intersect) {
        sets <- ins(lst = sets)
      }
      as_feature(
        sets, x, analysis = " Venn 交集 ({bind(group)}: {bind(projects)})"
      )
    }
  })

alias_intersect_multi <- function(..., mode = c("main", "all", "index"), main = 1L, sep = "///")
{
  require(data.table)

  mode <- match.arg(mode)

  inputs <- list(...)
  k <- length(inputs)

  # Expand all inputs into long alias tables
  dt_list <- lapply(seq_along(inputs), function(i) {
    dt <- .expand_alias_dt(inputs[[i]], paste0("s", i, "_"), sep)
    dt[, set_id := i]
    dt
  })

  # Combine all tables
  dt_all <- rbindlist(dt_list)

  # Count how many distinct sets each alias appears in
  alias_sets <- dt_all[, .(set_count = uniqueN(set_id)), by = alias]

  # Keep only aliases present in all sets
  valid_alias <- alias_sets[set_count == k, alias]

  # If no overlap, return empty structure
  if (length(valid_alias) == 0) {
    if (mode == "main") {
      return(inputs[[main]][ integer(0) ])
    } else {
      return(vector("list", k))
    }
  }

  # Filter hits
  hit <- dt_all[alias %in% valid_alias]

  # Extract unique indices per set
  res_idx <- hit[, .(idx = unique(idx)), by = set_id]

  # Return based on mode
  if (mode == "main") {
    idx <- res_idx[set_id == main, idx]
    inputs[[main]][idx]
  } else if (mode == "index") {
    split(res_idx$idx, res_idx$set_id)
  } else if (mode == "all") {
    lapply(seq_along(inputs), function(i) {
      idx <- res_idx[set_id == i, idx]
      inputs[[i]][idx]
    })
  }
}


# Main function: compute intersection based on alias overlap
alias_intersect <- function(x, y, mode = c("x", "y", "both", "index"), sep = "///")
{
  require(data.table)

  mode <- match.arg(mode)

  # Expand both inputs into long alias tables
  x_dt <- .expand_alias_dt(x, "x_", sep)
  y_dt <- .expand_alias_dt(y, "y_", sep)

  # Set keys for fast join on alias
  setkey(x_dt, alias)
  setkey(y_dt, alias)

  # Perform join to find overlapping aliases
  hit <- x_dt[y_dt, nomatch = 0]

  # Extract matched indices directly (faster than gid parsing)
  x_idx <- unique(hit$idx)
  y_idx <- unique(hit$i.idx)

  # Return results based on mode
  if (mode == "x") {
    x[x_idx]
  } else if (mode == "y") {
    y[y_idx]
  } else if (mode == "both") {
    list(x = x[x_idx], y = y[y_idx])
  } else {
    list(x_index = x_idx, y_index = y_idx)
  }
}


# Expand input into a long-format data.table (robust to vector/list, names optional)
.expand_alias_dt <- function(lst, prefix, sep = "\\s*///\\s*") {

  n <- length(lst)

  # Handle names safely (vector without names → auto-generate)
  nm <- names(lst)
  if (is.null(nm)) {
    nm <- paste0("V", seq_len(n))
  } else {
    empty <- nm == "" | is.na(nm)
    nm[empty] <- paste0("V", which(empty))
  }

  # Build base table
  dt <- data.table(
    name = nm,
    idx  = seq_len(n)
  )

  # Split aliases (vectorized, avoid lapply overhead)
  alias_list <- strsplit(as.character(lst), sep)

  # Flatten once
  alias_vec <- trimws(unlist(alias_list, use.names = FALSE))

  # Map aliases back to original indices
  lens <- lengths(alias_list)
  dt <- dt[rep.int(seq_len(n), lens)]

  dt[, alias := alias_vec]

  # Remove invalid aliases
  dt <- dt[alias != "" & !is.na(alias)]

  # Deduplicate within each group
  dt <- unique(dt, by = c("idx", "alias"))

  # Create group ID (kept for compatibility, though no longer used downstream)
  dt[, gid := paste0(prefix, idx)]

  dt[, .(gid, name, idx, alias)]
}

