
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(hiseqr))
# plot_profile, plot_profile_ss


#------------------------------------------------------------------------------#
# make_metaplot (profile)
# save to "${out_dir}/5.matrix2profile/${prefix}_plotProfile.pdf
#' @param x character path to the YAML file
#' @param reverse_strand bool default strandness is dUTP, if is not, use `TRUE`
#' to reverse the strand
#'
#' @export
make_metaplot <- function(x, ...) {
  args <- .setup_profile_params(config = x, ...)
  # 4. output
  out_dir <- ifelse(rlang::has_name(args, "out_dir"), args$out_dir, "./")
  prefix  <- ifelse(rlang::has_name(args, "prefix"), args$prefix, "metaplot")
  mat <- file.path(out_dir, "2.bw2matrix", glue::glue("{prefix}.mat.gz"))
  mat_sens <- gsub(".mat.gz$", "_sens.mat.gz", mat)
  mat_anti <- gsub(".mat.gz$", "_anti.mat.gz", mat)
  # 5. make plot file
  sub_dir <- file.path(out_dir, "5.matrix2profile_R")
  sub_pdf <- file.path(sub_dir, glue::glue("{prefix}_plotProfile.pdf"))
  if(file.exists(sub_pdf) & ! isTRUE(args$overwrite)) {
    message(glue::glue("File exists: {sub_pdf}"))
    return(sub_pdf)
  }
  if(! dir.exists(sub_dir)) {
    dir.create(sub_dir, recursive = TRUE)
  }
  if(file.exists(mat_sens) & file.exists(mat_anti)) {
    if(isTRUE(args$reverse_strand)) {
      m1 <- mat_anti
      m2 <- mat_sens
    } else {
      m1 <- mat_sens
      m2 <- mat_anti
    }
    p <- plot_profile_ss(m1, m2, sub_pdf, !!!args)
  } else if(file.exists(mat)) {
    p <- plot_profile(mat, sub_pdf, !!!args)
    # x <- 1
  } else {
    warning(glue::glue("Could not found matrix file: {mat}"))
    return(NULL)
  }
  return(NULL)
}


#------------------------------------------------------------------------------#
# read yaml file
# expect output files
# 2.bw2matrix
# 3.matrxi2profile
# 4.matrix2heatmap
# 5.matrix2profile_R
.is_valid_avg_func <- function(x) {
  avg_list <- c("mean", "median", "min", "max", "sum", "std")
  if(inherits(x, "character")) {
    if(length(x) > 1) {
      message(glue::glue("Multiple avg_func found, choose the first: [{x[1]}]"))
      x <- x[1]
    }
    if(x %in% avg_list) {
      return(TRUE)
    } else {
      aa <- paste(avg_list, collapse = ", ")
      message(glue::glue("bin_avg_type expect [{aa}], got: {x}"))
    }
  } else {
    message(glue::glue("Unknown avg_func: {x}, expect character"))
  }
  return(FALSE)
}


.plot_profile_args <- function(...) {
  # default arguments
  args <- list(
    # averageType   = "mean", # from header, avg_func
    # colors        = NULL, # auto
    # start_label   = "TSS",
    # end_label     = "TES",
    # point_label   = "TSS",  # same as start label
    # group_labels  = NULL, # from header, regionsLabel
    # linewidth     = 0.6,
    # return_data   = FALSE, # return data.frame
    # sample_labels = NULL, # from header, samplesLabel
    # sample_list   = NULL,
    # x_title       = NULL, # from header
    # plot_theme    = NULL, # default: theme_bw()
    # y_max         = NULL, # auto
    # y_min         = NULL, # auto
    # y_title       = NULL  # from header
    # yMax       = NULL, # from deeptool.plotProfile
    # yMin       = NULL  # from deeptool.plotProfile
    reverse_strand = FALSE,
    out_dir    = "./",
    prefix     = "metaplot",
    width      = 5,
    height     = 2.5,
    units      = "in",
    dpi        = 300,
    n_per_row  = 3,    # from config.yaml
    per_group  = TRUE, # from config.yaml
    plot_title = "metaplot",
    overwrite  = FALSE,
    add_x_ticks_extra = FALSE,
    add_genebody_bar  = TRUE
  )
  purrr::list_modify(args, ...) # items of list to environment
}


#' Count digits of float number
#'
#' @param x float numbers
#'
#' @return
#' @export
#'
#' @examples
#' count_digits(c(0.1, 0.02))
count_digits <- function(x) {
  sapply(x, function(i) {
    if(i %% 1 != 0) {
      s <- strsplit(gsub("0+$", "", as.character(i)), ".", fixed = TRUE)
      nchar(unlist(s)[2])
    } else {
      0
    }
  })
}


#' Fix y-axis ticks
#'
#' Force to show min/max ticks
#'
#' @param x ggplot A ggplot
#' @param axis character x or y axis, default: [y],
#'   currently, only support y-axis.#'
#' @return
#' @export
#'
#' @examples
fix_axis_text <- function(x, y_min = NULL, y_max = NULL) {
  # check arguments
  if(! inherits(x, "ggplot")) {
    warning("Not a ggplot input, skipped ...")
    return(x)
  }
  # check ticks/breaks on axis
  axis <- get_ggplot_axis(x)
  y_breaks <- purrr::discard(axis$y_breaks, is.na)
  y_limits <- axis$y_limits
  y_limits <- c(
    min(y_limits, y_breaks, y_min),
    max(y_limits, y_breaks, y_max)
  )
  y_breaks <- scales::breaks_extended(only.loose = TRUE)(y_limits)
  y_limits <- c(
    min(y_limits, y_breaks, y_min),
    max(y_limits, y_breaks, y_max)
  )
  # fix digits
  n_digits <- max(count_digits(y_breaks))
  accuracy <- 10^(- n_digits)
  # update y-axis
  suppressMessages(
    x +
      scale_y_continuous(
        limits = y_limits,
        breaks = y_breaks,
        labels = scales::number_format(accuracy = accuracy),
        expand = expansion(mult = c(0, 0.06))
      )
  )
}


#' Read arguments from profile YAML
#'
#' @param x character path to the `.yaml` file
#' @param ... extra options for arguments
#'
#' @return
#' @export
#'
#' @examples
#' # read yaml
#' args <- read_profile_yaml(x)
read_profile_yaml <- function(x, ...) {
  # read config.yaml for profile/metaplot
  # see deeptools.plotProfile arguments
  args_default <- list(
    width      = 5,
    height     = 2.5,
    units      = "in",
    dpi        = 300,
    linewidth  = 0.6,
    overwrite  = FALSE,
    n_per_row  = 3   # from config.yaml
  )
  args <- .plot_profile_args(!!!args_default)
  # read yaml
  args1 <- tryCatch(
    {
      yaml::read_yaml(x)
    },
    error = function(cond) {
      message(glue::glue("Failed to YAML file: {x}"))
      message(cond)
      # Choose a return value in case of error
      return(NULL)
    }
  )
  args <- purrr::list_modify(args, !!!args1) # overwrite by yaml
  # args1 <- purrr::list_modify(args1, ...) # customized
  ## renmae argument names for downstream plots
  args2 <- rlang::list2(
    # add_x_ticks_extra = FALSE,
    extra       = args$add_x_ticks_extra,
    bin_size    = args$binSize,
    plot_title  = args$plotTitle,
    start_label = args$startLabel,
    point_label = args$refPointLabel,
    end_label   = args$endLabel,
    x_breaks    = args$x_ticks,
    x_title     = args$xAxisLabel,
    y_title     = args$yAxisLabel,
    y_min       = args$yMin,
    y_max       = args$yMax,
    per_group   = args$perGroup,
    n_per_row   = args$numPlotsPerRow,
    group_labels  = args$regionsLabel,
    sample_labels = args$samplesLabel,
    bin_avg_type  = args$averageType
  )
  args <- purrr::list_modify(args, !!!args2) # rename arguments
  args <- purrr::list_modify(args, ...) # overwrite by ...
  args <- purrr::discard(args, is.null) # remove null arguments
  # fix colors
  if(inherits(args$colors, "character")) {
    args$colors <- strsplit(args$colors, "\\s+", perl = TRUE) %>% unlist
  }
  args
}


#' Initate params for matrix2plot
#' combine parameters from matrix and config.yaml file
#'
#' @param m str Path to matrix file or config file, default NULL
#' @param ... optinoal arguments from config.yaml file
#'
#' @return
#' @export
#'
#' @examples
#'  .setup_param()
.setup_profile_params <- function(m = NULL, config = NULL, ...) {
  # load defaults
  args <- .plot_profile_args() # default
  # read matrix header
  # args <- list() # empty
  if(inherits(m, "character")) {
    if(endsWith(m, ".mat.gz")) {
      h1 <- .read_matrix_header(m) # extract from matrix header
      h2 <- rlang::list2(
        sample_labels = h1$smp_labels,
        group_labels  = h1$group_labels,
        averageType   = h1$bin_avg_type,
        bin_size      = h1$bs,
        x_breaks      = h1$x_ticks
        # per_group     = TRUE # default
      )
      h <- modifyList(h1, h2)
      args <- modifyList(args, h)
      # args <- purrr::list_modify(args, !!!h1)
      # args <- purrr::list_modify(args, !!!h2) # args1
    }
  }
  # read config
  if(inherits(config, "character")) {
    if(endsWith(config, ".yaml")) {
      args2 <- read_profile_yaml(config)
      args  <- purrr::list_modify(args, !!!args2) # updated
    }
  }
  # read from ...
  purrr::list_modify(args, ...)
  # # Extract params from matrix file
  # if(inherits(x, "character")) {
  #   if(endsWith(x, ".mat.gz")) {
  #     h1 <- .read_matrix_header(x) # extract from matrix header
  #     h2 <- rlang::list2(
  #       sample_labels = h1$smp_labels,
  #       group_labels  = h1$group_labels,
  #       averageType   = h1$bin_avg_type,
  #       bin_size      = h1$bs,
  #       per_group     = TRUE # default
  #     )
  #     h <- purrr::list_modify(h2, !!!h1)
  #   } else if(endsWith(x, "yaml")) {
  #     h <- tryCatch(
  #       {
  #         yaml::read_yaml(x)
  #       },
  #       error = function(cond) {
  #         message(glue::glue("Failed to YAML file: {x}"))
  #         message(cond)
  #         # Choose a return value in case of error
  #         return(NULL)
  #       }
  #     )
  #   } else {
  #     h <- list()
  #   }
  # } else {
  #   h <- list() # blank
  # }
  # # optional from config.yaml
  # # see deeptools.plotProfile arguments
  # ## default
  # args_default <- list(
  #   width         = 5,
  #   height        = 2.5,
  #   units         = "in",
  #   dpi           = 300,
  #   line_size     = 0.6,
  #   numPlotsPerRow = 3,   # from config.yaml
  #   overwrite     = FALSE
  # )
  # args1 <- purrr::list_modify(args_default, ...)
  # # args1 <- .plot_profile_args(...) # default + custome
  # ## renmae argument names for downstream plots
  # args2 <- rlang::list2(
  #   add_x_ticks_extra = FALSE,
  #   bin_size    = args1$binSize,
  #   plot_title  = args1$plotTitle,
  #   start_label = args1$startLabel,
  #   point_label = args1$refPointLabel,
  #   end_label   = args1$endLabel,
  #   x_title     = args1$xAxisLabel,
  #   y_title     = args1$yAxisLabel,
  #   y_min       = args1$yMin,
  #   y_max       = args1$yMax,
  #   per_group   = args1$perGroup,
  #   group_labels  = args1$regionsLabel,
  #   sample_labels = args1$samplesLabel,
  #   bin_avg_type  = args1$averageType
  # )
  # args <- purrr::list_modify(args2, !!!args1) # modify
  # args <- purrr::discard(args, is.null) # remove null arguments
  # # fix colors
  # if(inherits(args$colors, "character")) {
  #   args$colors <- strsplit(args$colors, "\\s+", perl = TRUE) %>% unlist
  # }
  # # updated by optional arguments
  # purrr::list_modify(h, !!!args) # overwrite header
  # # purrr::list_modify(args, ...)
}


#' Update axis range
#'
#' @param x numeric values to plot on y-axis
#' @param y_min numeric minimum on y-axis, default: NULL
#' @param y_max numeric maximum on y-axis, default: NULL
#' @param ...
#'
#' @return numeric
#' @export
#'
#' @examples
#' update_axis_range_y(x, y_min = 0, y_max = 5)
update_axis_range_y <- function(x, y_min = NULL, y_max = NULL) {
  if(! inherits(x, "numeric")) {
    message(glue::glue(
      "`update_axis_range_y()` skipped, ",
      "expect numbers, got {class(x)[1]}"
    ))
    return(NULL)
  }
  # breaks <- scales::breaks_extended()(c(min(x), max(x)))
  # # fix y_min and y_max
  # y_min <- min(c(y_min, head(breaks, 1), min(x)))
  # y_max <- max(c(y_max, tail(breaks, 1), max(x)))
  y_min <- ifelse(inherits(y_min, "numeric"), y_min, min(x))
  y_max <- ifelse(inherits(y_max, "numeric"), y_max, max(x))
  # output
  list(y_min = y_min, y_max = y_max)
}


# read ggplot data
get_ggplot_axis <- function(p) {
  if(! inherits(p, "ggplot")) {
    message(glue::glue("Expect ggplot, got {class(p)[1]}"))
    return(p)
  }
  p_build <- ggplot2::ggplot_build(p)
  # output
  list(
    ## x axis
    x_discrete = p_build$layout$panel_params[[1]]$x.sec$is_discrete(),
    x_name   = p_build$layout$panel_params[[1]]$x.sec$name,      # name
    x_limits = p_build$layout$panel_params[[1]]$x.sec$get_limits(), # limits
    x_breaks = p_build$layout$panel_params[[1]]$x.sec$get_breaks(), # breaks
    x_labels = p_build$layout$panel_params[[1]]$x.sec$get_labels(), # breaks
    x_range  = p_build$layout$panel_params[[1]]$x.range, # range
    ## y axis
    y_discrete = p_build$layout$panel_params[[1]]$y.sec$is_discrete(),
    y_name   = p_build$layout$panel_params[[1]]$y.sec$name,         # name
    y_limits = p_build$layout$panel_params[[1]]$y.sec$get_limits(), # limits
    y_breaks = p_build$layout$panel_params[[1]]$y.sec$get_breaks(), # breaks
    y_labels = p_build$layout$panel_params[[1]]$y.sec$get_labels(), # breaks
    y_range  = p_build$layout$panel_params[[1]]$y.range # range
  )
}


#' Return the refpoints on x axis
#'
#' @param x_breaks numeric breaks on x axis
#' @param x_labels character labels on x axis
#'
#' @return
#' @export
#'
#' @examples
get_x_refpoint <- function(x_breaks, x_labels) {
  if(inherits(x_breaks, "numeric") & inherits(x_labels, "character")) {
    if(length(x_breaks) == length(x_labels)) {
      i <- grep("(TSS|TES|PAS)", x_labels, ignore.case = TRUE)
      if(length(i) %in% 1:2) {
        return(
          list(
            x_breaks = x_breaks[i],
            x_labels = x_labels[i]
          )
        )
      }
    }
  }
}


#' Breaks and labels matched
#'
#' @param x_breaks numeric breaks on axis
#' @param x_labels character labels on axis
#'
#' @return
#' @export
#'
#' @examples
is_valid_labels <- function(x_breaks = NULL, x_labels = NULL) {
  if(inherits(x_breaks, "numeric") & inherits(x_labels, "character")) {
    out <- length(x_breaks) == length(x_labels)
  } else {
    out <- FALSE
  }
  return(out)
}


#' Basic profile plot
#'
#' Add main lines and titles
#'
#' This is the sub-function of `plot_profile`, the following arguments
#' @param df data.frame is the matrix produced by `deeptools.computeMatrix`
#'   and processed by function `read_matrix()`, contain the following
#'   columns: x, score, sample_labels, group_labels
#' @param x_breaks int The position of ticks on x-axis,
#' @param x_labels character The labels on x-axis, should be identical
#'   in length with `x_ticks`
#' @param y_min int Manually select the minimum value of y-axis, default [NULL]
#'   automatically generated
#' @param y_max int As `y_min`, default [NULL]
#' @param dashed_lines bool Whether plot dashed line on genebody/TSS/TES ticks
#'   default [TRUE]
#' @param color_by character Color the lines by samples or groups, choices:
#'   sample_labels, group_labels, default [sample_lables]
#' @param ... optional arguments
#'
#' @return
#' @export
#'
#' @examples
#' # generate a basic plot with config.yaml and matrix file
#' x <- "pol2.genebody.mat.gz"
#' df <- .read_matrix(x, bin_avg_type = "mean")
#' plot_profile_basic(
#'   df, x_ticks = c(0, 40, 140, 340),
#'   x_labels = c("-2", "TSS", "TES", "+10"),
#'   y_title = "Mean of log2(IP/input)",
#'   x_tilte = "Genomic regions (kb)"
#' )
#'
plot_profile_basic <- function(df, x_breaks = NULL, x_labels = NULL,
                               y_min = NULL, y_max = NULL,
                               add_dashed_vlines = TRUE,
                               add_dashed_hlines = FALSE,
                               yintercept = 0,
                               color_by = "sample_labels", ...) {
  # required data.frame
  if(! inherits(df, "data.frame")) {
    message(glue::glue("Expect data.frame, got {class(df)[1]}"))
    return(NULL)
  }
  req_cols <- c("x", "score", "sample_labels", "group_labels")
  if(! all(req_cols %in% colnames(df))) {
    req_str <- paste(req_cols, collapse = ", ")
    df_str  <- paste(colnames(df)[1:4], collapse = ", ")
    message(glue::glue(
      "Expect columns: {req_str}, got: [{df_str}]"
    ))
    return(NULL)
  }
  #----------------------------------------------------------------------------#
  # 0. Arguments
  #  for dashed lines
  args <- list(
    color     = "grey50",
    linewidth = 0.3,
    linetype  = "dashed"
  )
  # dots <- rlang::list2(...)
  # dots <- dots[names(args)] %>% purrr::discard(is.null)
  # args <- purrr::list_modify(args, !!!dots) # update
  #----------------------------------------------------------------------------#
  # 1. Basic plot
  p <- ggplot(df)
  y_range <- update_axis_range_y(df$score, y_min, y_max)
  gb <- get_x_refpoint(x_breaks, x_labels)
  if(inherits(gb, "list") & isTRUE(add_dashed_vlines)) {
    # segment: dashed lines
    seg_v = list(
      geom = "segment",
      x    = gb$x_breaks,
      xend = gb$x_breaks,
      y    = -Inf,
      yend = y_range$y_max
    )
    # geom_vline
    vline = list(
      xintercept = gb$x_breaks
    )
    # if(length(gb$x_breaks) == 0) {
    #   p <- p + do.call(annotate, purrr::list_modify(seg_v, !!!args))
    # } else {
    #   p <- p + do.call(geom_vline, purrr::list_modify(vline, !!!args))
    # }
    p <- p +
      do.call(annotate, purrr::list_modify(seg_v, !!!args)) # +
      # scale_x_continuous(breaks = x_breaks, labels = x_labels)
  }
  if(inherits(gb, "list") & isTRUE(add_dashed_hlines)) {
    # segment: dashed line
    seg_h = list(
      geom = "segment",
      x    = -Inf,
      xend = max(x_breaks),
      y    = yintercept,
      yend = yintercept
    )
    p <- p +
      do.call(annotate, purrr::list_modify(seg_h, !!!args))
  }
  #----------------------------------------------------------------------------#
  # 2. Add main lines
  if(! color_by %in% c("group_labels", "sample_labels")) {
    message(glue::glue(
      "unknown `color_by=` {color_by}, ",
      "expect [sample_labels, group_labels]"
    ))
    color_by <- "sample_labels"
  }
  p <- p +
    geom_line(aes(x, score, color = .data[[color_by]]))
  #----------------------------------------------------------------------------#
  # 2. Optimize x axis ticks/labels
  dots  <- rlang::list2(...)
  if(is_valid_labels(x_breaks, x_labels)) {
    extra <- ifelse(rlang::has_name(dots, "extra"), dots$extra, TRUE)
    p <- update_axis_ticks_x(p, x_breaks, x_labels, extra)
  }
  #----------------------------------------------------------------------------#
  # 3. Optimize y axis range
  y_limits <- c(y_range$y_min, y_range$y_max)
  y_breaks <- scales::breaks_extended(only.loose = TRUE)(y_limits)
  y_limits2 <- c(min(y_limits, y_breaks), max(y_limits, y_breaks)) # updated
  n_digits <- count_digits(y_breaks)
  accuracy <- 10^(-max(n_digits))
  p <- p +
    scale_y_continuous(
      limits = y_limits2,
      breaks = y_breaks,
      labels = scales::number_format(accuracy = accuracy),
      expand = expansion(mult = c(0, .05))
    )
  #----------------------------------------------------------------------------#
  # 4. Add titles
  x_title <- ifelse(rlang::has_name(dots, "x_title"), dots$x_title,
                    "Genomic region (kb)")
  y_title <- ifelse(rlang::has_name(dots, "y_title"), dots$y_title,
                    "Score")
  plot_title <- ifelse(rlang::has_name(dots, "plot_title"), dots$plot_title,
                       "Profile")
  p + labs(x = x_title, y = y_title,title = plot_title)
}


#' Update ticks and labels on X-axis
#'
#' update x_ticks, x_labels
#' original:
#' 0 40 240;
#' "-2" "TES" "+10"
#' updated:
#' 0 40 80 120 160 200 240;
#' "-2"  "TES" "+2"  "+4"  "+6"  "+8"  "+10"
#'
#' @param x_breaks numeric breaks positions on x-axis
#' @param x_labels string labels on x-axis
#' @param extra bool Add extra x-ticks, for 2-kb breaks, default [TRUE]
#' @param ...
#'
#' @return
#' @examples
#' # update
#' ticks <- update_axis_ticks_x(x_ticks, x_labels)
update_axis_ticks_x <- function(x, x_breaks = NULL, x_labels = NULL,
                                extra = TRUE) {
  # Check arguments
  if(! inherits(x, "ggplot")) {
    warning("Not a ggplot input, skipped ...")
    return(x)
  }
  axis <- get_ggplot_axis(x) # read original layout
  if(is.null(x_breaks)) {
    x_breaks <- axis$x_breaks
  }
  if(is.null(x_labels)) {
    x_labels <- axis$x_labels
  }
  if(inherits(x_breaks, "numeric") & inherits(x_labels, "character")) {
    if(length(x_breaks) == length(x_labels)) {
      axis$x_breaks <- x_breaks
      axis$x_labels <- x_labels
      axis$x_limits <- c(
        min(x_breaks, axis$x_limits),
        max(x_breaks, axis$x_limits)
      )
      # Check the end of x ticks
      x_end <- tail(x_labels, 2) # last two items
      x_end <- gsub("[^0-9]+", "", x_end) %>% as.numeric()
      # Update ticks and labels
      if(is.na(x_end[1]) & x_end[2] > 2 & isTRUE(extra)) {
        n_breaks <- floor(x_end[2] / 2) + 1
        breaks   <- tail(x_breaks, 2) # ticks
        breaks_extra <- seq(breaks[1], breaks[2], length.out = n_breaks)
        labels_extra <- seq(0, x_end[2], length.out = n_breaks)
        labels_extra <- paste0('+', labels_extra)
        # remove first one
        axis$x_breaks <- c(head(x_breaks, -1), tail(breaks_extra, -1))
        axis$x_labels <- c(head(x_labels, -1), tail(labels_extra, -1))
      }
    }
  }
  # Update x-axis
  x_limits2 <- c(min(axis$x_breaks), max(axis$x_breaks))
  suppressMessages(
    x +
      scale_x_continuous(
        limits = x_limits2,
        breaks = axis$x_breaks,
        labels = axis$x_labels,
        expand = expansion(mult = c(0, 0))
      )
  )
}


#' Update axis lines by segment
#'
#' keep axis lines between axis ticks
#'
#' @param x ggplot
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' add_axis_line(x)
update_axis_line <- function(x, ...) {
  # Check arguments
  args <- list(
    color     = "grey20",
    linewidth = .3,
    linetype  = "solid"
  )
  dots <- rlang::list2(...)
  dots <- dots[names(args)] %>% purrr::discard(is.null)
  args <- purrr::list_modify(args, !!!dots)
  if(! inherits(x, "ggplot")) {
    warning("Not a ggplot input, skipped ...")
    return(x)
  }
  #----------------------------------------------------------------------------#
  # Replace axis-line by segments
  # Fix error: panel.border was clipped
  # Segment on axis position was clipped (smaller than linewidth)
  # use: coord_cartesian(clip = "off")
  axis <- get_ggplot_axis(x)
  y_breaks <- scales::breaks_extended(only.loose = TRUE)(axis$y_limits)
  # # axis line
  seg_hv <- list(
    geom = "segment",
    x    = c(-Inf, -Inf),
    xend = c(tail(axis$x_breaks, 1), -Inf),
    y    = c(-Inf, -Inf),
    yend = c(-Inf, tail(axis$y_breaks, 1))
  )
  x +
    do.call(annotate, c(seg_hv, args)) +
    coord_cartesian(clip = "off") +
    theme(axis.line.x = element_blank(),
          axis.ticks = element_line(linewidth = args$linewidth))
}


#' Add genebody bar on top of metaplot
#'
#' @param x
#' @param expand
#' @param margin
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
update_genebody_bar <- function(x, expand = 0.1, margin = 0, ...) {
  #----------------------------------------------------------------------------#
  # Check arguments
  if(! inherits(x, "ggplot")) {
    warning("Not a ggplot input, skipped ...")
    return(x)
  }
  args <- list(
    fill      = "grey50",
    color     = "grey50",
    linewidth = .2,
    linetype  = "solid"
  )
  dots <- rlang::list2(...)
  dots <- dots[names(args)] %>% purrr::discard(is.null)
  args <- purrr::list_modify(args, !!!dots)
  #----------------------------------------------------------------------------#
  # 1. Fetch axis info from plot
  axis <- get_ggplot_axis(x)
  # Checkpoint-2: x axis
  if(isTRUE(axis$x_discrete) || isTRUE(axis$y_discrete)) {
    warning(glue::glue(
      "Expect numeric axises: x-axis: ",
      "{! args$x_discrete}, y-axis: {! args$y_discrete}, ",
      "update_genebody_bar() skipped..."
    ))
    return(x)
  }
  #----------------------------------------------------------------------------#
  # 2. Fetch TSS/TES/PAS ticks, for genebody bar
  gb <- get_x_refpoint(axis$x_breaks, axis$x_labels)
  x_str <- paste(axis$x_labels, collapse = ", ")
  if(inherits(gb, "list")) {
    if(length(gb$x_breaks) == 2) {
      gb_breaks <- gb$x_breaks
    } else {
      message(glue::glue("Genebody not found: {x_str}"))
      return(x)
    }
  } else {
    message(glue::glue("Genebody not found: {x_str}"))
  }
  #----------------------------------------------------------------------------#
  # 3. Update y-axis range/limits/labels/ticks
  # Expand y-axis by `expand=`, `margin=`
  y_margin  <- abs(diff(axis$y_limits)) * margin
  y_exp     <- abs(diff(axis$y_limits)) * expand
  y_seg     <- axis$y_limits[2] + y_margin + y_exp * c(0.1, 0.4, 0.5, 0.9)
  y_limits  <- c(axis$y_limits[1], axis$y_limits[2] + y_margin + y_exp)
  # y_breaks  <- scales::breaks_extended(only.loose = TRUE)(y_limits)
  # update limits
  # y_limits2 <- c(min(y_limits, y_breaks), max(y_limits, y_breaks))
  n_digits  <- max(count_digits(axis$y_breaks)) # fix digits width
  accuracy  <- 10^(- n_digits)
  #--------------------------#
  # 3.1 update y-axis ticks
  suppressMessages(
    x2 <- x +
      scale_y_continuous(
        limits = y_limits,
        breaks = axis$y_breaks,
        labels = scales::number_format(accuracy = accuracy),
        expand = expansion(mult = c(0, .05))
      )
  )
  #--------------------------#
  # 3.2 Add genebody-bar, arrows, lines
  # bar in horizontal
  rect_h <- list(
    "geom" = "rect",
    "xmin" = gb_breaks[1],
    "xmax" = gb_breaks[2],
    "ymin" = y_seg[1],
    "ymax" = y_seg[2]
  )
  suppressWarnings(
    p2 <- x2 +
      do.call(annotate, modifyList(rect_h, args))
  )
  # vertical line, with gap, + horizontal line + arrow
  seg_v <- list(
    "geom" = "segment",
    "x"    = gb_breaks[1],
    "xend" = gb_breaks[1],
    "y"    = y_seg[3],
    "yend" = y_seg[4]
  )
  seg_h <- list(
    "geom" = "segment",
    "x"    = gb_breaks[1],
    "xend" = gb_breaks[1] + abs(diff(gb_breaks)) * 0.15, # 15% of genebody
    "y"    = y_seg[4],
    "yend" = y_seg[4],
    arrow = arrow(length = unit(c(0.04, 0.04), "npc"))
  )
  suppressWarnings(
    p2 <- p2 +
      do.call(annotate, modifyList(seg_v, args)) +
      do.call(annotate, modifyList(seg_h, args))
  )
  #----------------------------------------------------------------------------#
  # output
  # update_axis_line(p2)
  p2
}


#------------------#
# theme global
theme_metaplot <- function(x) {
  x +
    theme_classic() +
    theme(
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.grid   = element_blank(),
      legend.title = element_blank(),
      rect         = element_blank(),
      plot.title   = element_text(hjust = .5, size = 10),
      axis.title   = element_text(color = "black", size = 9),
      axis.text    = element_text(color = "black", size = 8),
      axis.line    = element_blank(),
      axis.ticks   = element_line(color = "grey30", linewidth = .3),
      text         = element_text(family = "Helvetica")
    )
}


save_as_pdf <- function(x, file, ...) {
  #----------------------------------------------------------------------------#
  # 1. Check arguments
  if(! inherits(x, "ggplot")) {
    warning(glue::glue(
      "`save_as_pdf(x=)` is {class(x)[1]} but not ggplot, failed ..."
    ))
    return(NULL)
  }
  if(! inherits(file, "character")) {
    warning(glue::glue(
      "`save_as_pdf(file=)`, expect character, ",
      "got {class(file)}, failed ..."
    ))
    return(NULL)
  }
  #----------------------------------------------------------------------------#
  # 2. Save plot
  ## args for export::graph2pdf()
  args <- list(
    width  = 5,
    height = 2.5,
    font   = "Arial",
    bg     = "transparent"
  )
  dots  <- rlang::list2(...)
  dots2 <- dots[names(args)] %>% purrr::discard(is.null)
  args  <- purrr::list_modify(args, !!!dots2)
  overwrite <- ifelse(rlang::has_name(dots, "overwrite"), dots$overwrite, FALSE)
  if(file.exists(file) && ! isTRUE(overwrite)) {
    message(glue::glue("file exists: {file}"))
  } else {
    pdf(NULL) # prevent generating empty file: "Rplot.pdf"
    message(glue::glue(
      "export pdf size, width={args$width}, height={args$height}"
    ))
    do.call(export::graph2pdf, c(args, list(x = x, file = file)))
    # export::graph2pdf(
    #   x = p, file = filename, width = width, height = height,
    #   font = "Arial", bg = "transparent"
    # )
    dev.off() # close empty pdf
    rds <- gsub("\\.[a-z]+$", ".rds", file, perl = TRUE)
    saveRDS(x, file = rds)
  }
}


# reference-point
#' see `computeMatrix` of deeptools at
#' https://deeptools.readthedocs.io/en/develop/content/tools/plotProfile.html
#'
#' @param x Character path to the matrix file
#' @param header_only bool return the head only
#'
#' sl # sample_labels
#' ss # list of sample_labels
#' x_axis   = x_axis,
#' x_ticks  = x_ticks,
#' x_labels = x_labels,
#' x_title  = x_title,
#' y_title  = y_title,
#' x_sect   = x_sect
#'
#' @export
.read_matrix_header <- function(x) {
  h <- readLines(x, n = 1)
  j <- jsonlite::parse_json(gsub("^@", "", h))
  # 0. average type
  bin_avg_type <- j$`bin avg type` # !!!
  # 1. sample labels (bw files)
  smp_labels <- unlist(j$sample_labels)
  smp_bound  <- unlist(j$sample_boundaries)
  # 2. group labels (region files)
  group_labels <- unlist(j$group_labels)
  group_bound  <- unlist(j$group_boundaries)
  # 3. axis
  y_title <- "Mean of score"
  ## up, TSS, TES, down
  us <- ifelse(rlang::has_name(j, "upstream"), j$upstream[[1]], 0)
  gb <- ifelse(rlang::has_name(j, "body"), j$body[[1]], 1000)
  ds <- ifelse(rlang::has_name(j, "downstream"), j$downstream[[1]], 0)
  bs <- ifelse(rlang::has_name(j, "bin size"), j$`bin size`[[1]], 10)
  u5 <- ifelse(rlang::has_name(j, "unscaled 5 prime"), j$`unscaled 5 prime`[1], 0)
  u3 <- ifelse(rlang::has_name(j, "unscaled 3 prime"), j$`unscaled 3 prime`[1], 0)
  # x-axis ticks labels
  # warnings: us or ds < 1000 ?
  usl <- ifelse(any(c(us, ds) > 1000), paste0("-", round(us / 1000, 1)),
                paste0("-", round(us, 0)))
  dsl <- ifelse(any(c(us, ds) > 1000), paste0("+", round(ds / 1000, 1)),
                paste0("+", round(ds, 0)))
  x_title  <- ifelse(any(c(us, ds) > 1000), "Genomic region (kb)",
                     "Genome region (bp)")
  ## x ticks
  ref <- j$`ref point`[[1]]
  if(is.null(ref)) {
    matrix_type <- "scale-regions" # see deeptools.computeMatrix
    x_labels <- c(usl, "TSS", "TES", dsl)
    x_ticks  <- cumsum(round(c(0, us, gb, ds) / bs, 0))
  } else {
    matrix_type <- "reference-point" # see deeptools.computeMatrix
    x_labels <- c(usl, ref, dsl)
    x_ticks  <- cumsum(round(c(0, us, ds) / bs, 0))
  }
  # 4. output
  # result <- list(bin_avg_type, smp_labels, smp_bound, group_labels, group_bound,
  #                x_title, x_ticks, x_labels, y_title, matrix_type)
  # return(result)
  return(as.list(environment()))
}


#----------------------------------------------------------------------------#
# 2. load matrix
.read_matrix <- function(x, bin_avg_type = NULL) {
  # 1. loading header
  header <- .read_matrix_header(x)
  if(is.null(bin_avg_type)) bin_avg_type <- header$bin_avg_type #
  if(! .is_valid_avg_func(bin_avg_type)) return(NULL) # exception
  # 2. loading matrix
  message(glue::glue("Loading matrix from file: {basename(x)}"))
  df <- read.delim(x, header = FALSE, sep = "\t", comment.char = "@")
  df1 <- subset(df, select = -c(1:6)) # remove bed6 record
  hd_col <- max(header$smp_bound)
  hd_row <- max(header$group_bound)
  # check matrix size: col-row
  if(! all(dim(df1) == c(hd_row, hd_col))) {
    msg <- glue::glue(
      "Error, expect matrix size [{dim(df)[1]}, {dim(df)[2]}]",
      ", actual matrix size is [{hd_row}, {hd_col}]"
    )
    warning(msg)
    return(NULL)
  }
  # 3. split data.frame by sample-group
  # to-do: split data.frame more elegantly/efficiently using base R functions?
  ## group list
  gl <- mapply(rep, x = header$group_labels, times = diff(header$group_bound))
  gl <- unlist(gl)
  df2 <- lapply(seq_len(length(header$smp_labels)), function(i) {
    smp_bd <- header$smp_bound
    s1  <- smp_bd[i] + 1
    s2  <- smp_bd[i + 1]
    dfs <- subset(df1, select = c(s1:s2))
    dfs$group_label <- gl
    # split by groups
    dfs_list <- split(dfs, dfs$group_label) # split by group_labels
    dx <- lapply(dfs_list, function(j) {
      dj <- subset(j, select = -group_label)
      data.frame(
        row.names = NULL,
        x     = seq_len(smp_bd[2]),
        score = apply(dj, 2, match.fun(bin_avg_type)), # avg values by group
        sample_labels = header$smp_labels[i],
        group_labels  = j$group_label[1]
      )
    }) %>%
      Reduce(rbind, .)
    # dplyr::bind_rows()
  }) %>%
    Reduce(rbind, .)
  # dplyr::bind_rows()
  # # to-do: split data.frame more elegantly/efficiently using base R functions?
  #
  # df2 <- lapply(seq_len(length(header$smp_labels)), function(i) {
  #   # choose the boundaries
  #   smp_bd <- header$smp_bound
  #   s1  <- smp_bd[i] + 1
  #   s2  <- smp_bd[i + 1]
  #   dfs <- subset(df1, select = c(s1:s2))
  #   # return
  #   lapply(seq_len(length(header$group_labels)), function(j) {
  #     group_bd <- header$group_bound
  #     g1  <- group_bd[j] + 1
  #     g2  <- group_bd[j + 1]
  #     dfg <- dfs[g1:g2, ]
  #     ss <- data.frame(
  #       row.names = NULL,
  #       x     = seq_len(smp_bd[2]),
  #       score = apply(dfg, 2, match.fun(bin_avg_type)),
  #       sample_label = header$smp_labels[i],
  #       group_label  = header$group_labels[j]
  #     )
  #   }) %>%
  #     dplyr::bind_rows()
  # }) %>%
  #   dplyr::bind_rows()
  # 4. fix order
  df2$sample_labels <- factor(df2$sample_labels, levels = header$smp_labels)
  df2$group_labels  <- factor(df2$group_labels, levels = header$group_labels)
  df2
}


#' Generate metaplot for regions/reference-point
#'
#' Reference
#' 1. Kamieniarz-Gdula, K. et al. Selective Roles of Vertebrate PCF11 in
#' Premature and Full-Length Transcript Termination. Molecular Cell 74,
#' 158-172.e9 (2019). doi: 10.1016/j.molcel.2019.01.027
#' Figure 1B
#' url: https://doi.org/10.1016/j.molcel.2019.01.027
#'
#' Structures of metaplot
#' Add the layers as the following order
#' 1. `plot_profile_basic()` add main lines, x, y title and plot title
#'
#'
#'
#' see `plotProfile` of deeptools at
#' https://deeptools.readthedocs.io/en/develop/content/tools/plotProfile.html
#'
#' @param m character path to the matrix file, output of `computeMatrix`
#'   in `deeptools`
#' @param filename path to the plot file
#'
#' @importFrom rlang list2 is_empty
#' @importFrom purrr list_modify
#' @importFrom jsonlite parse_json
#' @importFrom ggthemes scale_color_few theme_few
#' @importFrom ggplot2 ggplot geom_vline geom_line scale_x_continuous ggtitle
#'   theme scale_color_manual scale_y_continuous theme_bw theme
#'
#' @export
plot_profile <- function(m, filename = NULL, ...) {
  # parse arugments from matrix_header + ...
  args <- .setup_profile_params(m = m, config = NULL, ...) # ... from config?
  # loading header from matrix
  df <- .read_matrix(x = m, bin_avg_type = args$bin_avg_type) #
  #----------------------------------------------------------------------------#
  # Checkpoint-1: update group labels
  group_labels <- unique(df$group_labels)
  if(inherits(args$group_labels, "character")) {
    if(length(args$group_labels) == length(group_labels)) {
      gl <- setNames(object = args$group_labels, nm = group_labels)
      df$group_labels <- gl[df$group_label]
    }
  }
  color_by <- ifelse(isTRUE(args$per_group), "sample_labels", "group_labels")
  group_by <- ifelse(isTRUE(args$per_group), "group_labels", "sample_labels")
  # 3.0 determine fig numbers
  n_plots  <- length(unique(df[[group_by]]))
  n_width  <- ifelse(n_plots > args$n_per_row, args$n_per_row, n_plots)
  n_height <- ceiling(n_plots / args$n_per_row)
  args$width    <- n_width * args$width
  args$height   <- n_height * args$height
  #----------------------------------------------------------------------------#
  # Checkpoint-1. Basic plot
  # p <- plot_profile_basic(df, !!!args)
  args$df <- df
  args$color_by <- color_by
  p <- do.call(plot_profile_basic, args)
  #----------------------------------------------------------------------------#
  # Checkpoint-2. Update axis lines
  axis_linewidth <- ifelse(
    rlang::has_name(args, "axis_linewidth"), args$axis_linewidth, 0.3
  )
  axis_color <- ifelse(
    rlang::has_name(args, "axis_color"), args$axis_color, "grey20"
  )
  p <- update_axis_line(p, linewidth = axis_linewidth, color = axis_color) #
  #----------------------------------------------------------------------------#
  # Checkpoint-3. Update genebody-bar
  if(args$matrix_type == "scale-regions" & isTRUE(args$add_genebody_bar)) {
    p <- update_genebody_bar(p, expand = .1, margin = 0)
  }
  #----------------------------------------------------------------------------#
  # Checkpoint-4. Update colors for main lines
  n_colors <- length(unique(df[[color_by]]))
  if(inherits(args$colors, "character")) {
    colors <- rep(args$colors, 100)[1:n_colors] #
    p <- p + scale_color_manual(values = colors)
  }
  #----------------------------------------------------------------------------#
  # Checkpoint-5. Update themes
  p <- theme_metaplot(p)
  #----------------------------------------------------------------------------#
  # Checkpoint-6. Facet for multiple groups
  if(length(group_labels) > 1) {
    p <- p +
      facet_wrap(as.formula(paste("~", group_by)), ncol = n_width)
  }
  #----------------------------------------------------------------------------#
  # Checkpoint-7. Save to PDF
  save_as_pdf(p, filename, !!!args)
  # output
  p
}


#' Generate strand-specific metaplot
#'
#' see `plotProfile` of deeptools at
#' https://deeptools.readthedocs.io/en/develop/content/tools/plotProfile.html
#'
#' @param m1 matrix file for sense-strand, see `computeMatrix`
#' @param m2 matrix file for anti-sense-strand,
#' @param ... additional parameters, see arguments in `plot_profile` function
#'
#' @export
plot_profile_ss <- function(m1, m2,  filename = NULL, ...) {
  # parse arugments from matrix_header + ...
  args <- .setup_profile_params(m = m1, config = NULL, ...) # ... from config?
  #----------------------------------------------------------------------------#
  # 0. Loading matrix from file
  # sample_labels were saved in matrix header, group_labels not !!!
  # update group_labels by values from args
  # sens
  df1 <- .read_matrix(x = m1, bin_avg_type = args$bin_avg_type) # sene
  h1  <- .read_matrix_header(m1)
  if(length(unique(df1$group_labels)) == length(h1$group_labels)) {
    gl  <- setNames(args$group_labels, nm = unique(df1$group_labels))
    df1 <- dplyr::mutate(df1, group_labels = gl[group_labels])
  }
  # anti
  df2 <- .read_matrix(x = m2, bin_avg_type = args$bin_avg_type) # anti
  h2  <- .read_matrix_header(m2)
  if(length(unique(df2$group_labels)) == length(h2$group_labels)) {
    gl  <- setNames(args$group_labels, nm = unique(df2$group_labels))
    df2 <- dplyr::mutate(df2, group_labels = gl[group_labels])
  }
  df2 <- dplyr::mutate(df2, score = -score)
  #----------------------------------------------------------------------------#
  # 1. Update `color_by` column; add `_fwd` and `_rev` suffix
  color_by <- ifelse(isTRUE(args$per_group), "sample_labels", "group_labels")
  df1x <- dplyr::mutate(df1, !!color_by := paste0(.data[[color_by]], "_fwd"))
  df2x <- dplyr::mutate(df2, !!color_by := paste0(.data[[color_by]], "_rev"))
  # factor
  lv <- paste0(
    rep(args[[color_by]], times = 2),
    rep(c("_fwd", "_rev"), each = length(args[[color_by]]))
  )
  df <- dplyr::bind_rows(df1x, df2x) %>%
    dplyr::mutate(
      !!color_by := factor(.data[[color_by]], levels = lv)
    )
  # #----------------------------------------------------------------------------#
  # # 2. Update `group_by` column
  # group_labels <- unique(df[[group_by]])
  # if(inherits(args$group_labels, "character")) {
  #   if(length(args$group_labels) == length(group_labels)) {
  #     gl <- setNames(object = args$group_labels, nm = group_labels)
  #     df[[group_by]] <- gl[df[[group_by]]]
  #   }
  # }
  #----------------------------------------------------------------------------#
  # 3.Determine fig arguments
  group_by <- ifelse(isTRUE(args$per_group), "group_labels", "sample_labels")
  n_plots  <- length(unique(df[[group_by]]))
  n_width  <- ifelse(n_plots > args$n_per_row, args$n_per_row, n_plots)
  n_height <- ceiling(n_plots / args$n_per_row)
  args$width  <- n_width * args$width
  args$height <- n_height * args$height
  #----------------------------------------------------------------------------#
  # Checkpoint-1. Basic plot
  args$df <- df
  args$color_by <- color_by
  args$add_dashed_hlines <- TRUE
  p <- do.call(plot_profile_basic, args)
  #----------------------------------------------------------------------------#
  # Checkpoint-2. Update axis lines
  axis_linewidth <- ifelse(
    rlang::has_name(args, "axis_linewidth"), args$axis_linewidth, 0.3
  )
  axis_color <- ifelse(
    rlang::has_name(args, "axis_color"), args$axis_color, "grey30"
  )
  p <- update_axis_line(p, linewidth = axis_linewidth, color = axis_color) #
  #----------------------------------------------------------------------------#
  # Checkpoint-3. Update genebody-bar
  if(args$matrix_type == "scale-regions" & isTRUE(args$add_genebody_bar)) {
    p <- update_genebody_bar(p, expand = .1, margin = 0)
  }
  #----------------------------------------------------------------------------#
  # Checkpoint-4. Update colors for main lines
  n_colors <- length(unique(df[[color_by]]))
  if(inherits(args$colors, "character")) {
    if(length(args$colors) >= n_colors) {
      colors <- args$colors
    } else if(length(args$colors) >= n_colors / 2) {
      colors <- rep(args$colors[seq_len(n_colors / 2)], times = 2)
    } else {
      colors <- rep(args$colors, 100)[1:n_colors]
    }
    # clip colors by length
    colors <- colors[seq_len(n_colors)]
    p <- p + scale_color_manual(values = colors)
  }
  #----------------------------------------------------------------------------#
  # Checkpoint-5. Update themes
  p <- theme_metaplot(p)
  #----------------------------------------------------------------------------#
  # Checkpoint-6. Facet for multiple groups
  if(n_plots > 1) {
    p <- p +
      facet_wrap(as.formula(paste("~", group_by)), ncol = n_width,
                 scales = "free")
  }
  #----------------------------------------------------------------------------#
  # Checkpoint-7. Save to PDF
  save_as_pdf(p, filename, !!!args)
  # output
  p
}


# scale-regions
# x = '/data/yulab/wangming/work/yu_2022/projects/20221229_dlj_ChrRNA_yy218/results/flanking_genes/results/fig1.gs_6k/2.bw2matrix/fig1.ChrRNA_YY218.gs_6k_sens.mat.gz'
# m1 = '/data/yulab/wangming/work/yu_2022/projects/20221229_dlj_ChrRNA_yy218/results/flanking_genes/results/fig1.gs_6k/2.bw2matrix/fig1.ChrRNA_YY218.gs_6k_sens.mat.gz'
# m2 = '/data/yulab/wangming/work/yu_2022/projects/20221229_dlj_ChrRNA_yy218/results/flanking_genes/results/fig1.gs_6k/2.bw2matrix/fig1.ChrRNA_YY218.gs_6k_anti.mat.gz'
# x = "data/config/metaplot/paused_all/fig1A/fig1.A.ChIP_8WG16_60m.metaplot.tss.yaml"

# am = "results/metaplot/paused_all/fig1A/2.bw2matrix/fig1.A.ChIP_8WG16_60m.tes.mat.gz"
# (p <- plot_profile(am, filename = "tmp.cnt.pdf", colors = c("black", "red"), overwrite = T, add_x_ticks_extra = T))

# ay = "data/config/metaplot/paused_all/fig1A/fig1.A.ChIP_8WG16_60m.metaplot.tes.yaml"
# make_metaplot(ay, overwrite = T, add_x_ticks_extra = T)
#
# am = "results/metaplot/paused_all/fig1E/2.bw2matrix/fig1.E.CnT_Ser2P_vs_total.genebody.mat.gz"
# ay = "data/config/metaplot/paused_all/fig1E/fig1.E.CnT_Ser2P_vs_total.metaplot.genebody.yaml"

# m1 = "results/metaplot/paused_all/fig1G/2.bw2matrix/fig1.G.ChrRNA_rmPcf11.genebody_sens.mat.gz"
# m2 = "results/metaplot/paused_all/fig1G/2.bw2matrix/fig1.G.ChrRNA_rmPcf11.genebody_anti.mat.gz"
# ay = "data/config/metaplot/paused_all/fig1G/fig1.G.ChrRNA_rmPcf11.metaplot.genebody.yaml"
#
