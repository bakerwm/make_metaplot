
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
  args <- .setup_profile_params(m = NULL, config = x, ...)
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
    add_x_ticks_extra = FALSE
  )
  purrr::list_modify(args, ...) # items of list to environment
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
    linewidth = 0.6,
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
    bin_size    = args$binSize,
    plot_title  = args$plotTitle,
    start_label = args$startLabel,
    point_label = args$refPointLabel,
    end_label   = args$endLabel,
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
        bin_size      = h1$bs
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
#' @param y_max numeric  maximum on y-axis, default: NULL
#' @param ...
#'
#' @return numeric
#' @export
#'
#' @examples
#' update_axis_range_y(x, y_min = 0, y_max = 5)
update_axis_range_y <- function(x, y_min = NULL, y_max = NULL) {
  if(! inherits(x, "numeric")) {
    message(glue::glue("Expect numbers, got {class(x)[1]}"))
    return(NULL)
  }
  # s  <- summary(x) #
  # scales::breaks_pretty(n=5)(c(min(x), max(x)))
  breaks <- scales::breaks_extended()(c(min(x), max(x)))
  # fix y_min and y_max
  y_min <- min(c(y_min, head(breaks, 1), min(x)))
  y_max <- max(c(y_max, tail(breaks, 1), max(x)))
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
    x_name   = p_build$layout$panel_params[[1]]$x.sec$name,      # name
    x_limits = p_build$layout$panel_params[[1]]$x.sec$get_limits(), # limits
    x_breaks = p_build$layout$panel_params[[1]]$x.sec$get_breaks(), # breaks
    x_labels = p_build$layout$panel_params[[1]]$x.sec$get_labels(), # breaks
    x_range  = p_build$layout$panel_params[[1]]$x.range, # range
    ## y axis
    y_name   = p_build$layout$panel_params[[1]]$y.sec$name,         # name
    y_limits = p_build$layout$panel_params[[1]]$y.sec$get_limits(), # limits
    y_breaks = p_build$layout$panel_params[[1]]$y.sec$get_breaks(), # breaks
    y_labels = p_build$layout$panel_params[[1]]$y.sec$get_labels(), # breaks
    y_range  = p_build$layout$panel_params[[1]]$y.range # range
  )
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
fix_axis_ticks <- function(x, y_min = NULL, y_max = NULL) {
  # check arguments
  if(! inherits(x, "ggplot")) {
    warning("Not a ggplot input, skipped ...")
    return(x)
  }
  #------------------------------------------#
  # count the number of digits
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
  #------------------------------------------#
  # check ticks/breaks on axis
  axis <- get_ggplot_axis(x)
  if(any(is.na(axis$y_breaks))) {
    y_limits <- axis$y_limits
    y_breaks <- scales::breaks_extended()(axis$y_limits)
    y_limits <- c(
      min(y_limits, y_breaks, y_min),
      max(y_limits, y_breaks, y_max)
    )
    ## update all
    y_breaks <- scales::breaks_extended(only.loose = TRUE)(y_limits) # updated
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
  } else {
    y_limits <- axis$y_limits
    y_breaks <- axis$y_breaks
    y_limits <- c(min(y_limits, y_breaks), max(y_limits, y_breaks))
    # fix digits
    n_digits <- max(count_digits(y_breaks))
    accuracy <- 10^(- n_digits)
    suppressMessages(
      x +
        scale_y_continuous(
          limits = y_limits,
          breaks = y_breaks,
          labels = scales::number_format(accuracy = accuracy) ,
          expand = expansion(mult = c(0, 0.06))
        )
    )
  }
}


#' Basic profile plot
#'
#' This is the sub-function of `plot_profile`, the following arguments
#' @param df data.frame is the matrix produced by `deeptools.computeMatrix`
#'   and processed by function `read_matrix()`, contain the following
#'   columns: x, score, sample_labels, group_labels
#' @param x_ticks int The position of ticks on x-axis,
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
plot_profile_basic <- function(df, x_ticks, x_labels, y_min = NULL,
                               y_max = NULL, add_dashed_lines = TRUE,
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
  # determine y-axis
  dots <- rlang::list2(...)
  gar  <- update_axis_range_y(df$score, y_min, y_max)
  # x-axis ticks, labels
  if(! all(
    inherits(x_ticks, "numeric"),
    inherits(x_labels, "character"),
    length(x_ticks) == length(x_labels)
  )) {
    message(glue::glue(
      "Expect x_ticks (numeric) and x_labels (character), ",
      "x_labels is {class(x_labels)}, ",
      "x_ticks is {class(x_ticks)}"
    ))
    return(NULL)
  }
  # x_sect <- head(tail(x_ticks, -1), -1) # remove first, last element !!!
  x_sect <- x_ticks[grep("[A-Z]", x_labels)]
  # segment, x_sect from outer
  p <- ggplot(df) # basic
  #----------------------------------------------------------------------------#
  # 1. Add dashed lines
  if(isTRUE(add_dashed_lines)) {
    if(length(x_sect) == 2) {
      # genebody: scale-regions
      p <- p +
        annotate(
          "segment",
          x     = x_sect,
          xend  = x_sect,
          y     = -Inf,
          yend  = gar$y_max,
          color = "grey50",
          linewidth = 0.3,
          linetype  = "dashed"
        )
    } else {
      # refpoint: reference-point
      p <- p +
        geom_vline(
          xintercept = x_sect,
          color      = "grey50",
          linewidth  = 0.3,
          linetype   = "dashed"
        )
    }
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
  # fix y_axis
  p <- fix_axis_ticks(p, y_min, y_max)
  #----------------------------------------------------------------------------#
  # 3. Add titles
  dots <- rlang::list2(...)
  x_title <- ifelse(rlang::has_name(dots, "x_title"), dots$x_title,
                    "Genomic region (kb)")
  y_title <- ifelse(rlang::has_name(dots, "y_title"), dots$y_title,
                    "Score")
  plot_title <- ifelse(rlang::has_name(dots, "plot_title"), dots$plot_title,
                       "Profile")
  p + labs(x = x_title, y = y_title,title = plot_title)
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
add_genebody_bar <- function(x, expand = 0.1, margin = 0, ...) {
  # check arguments
  if(! inherits(x, "ggplot")) {
    warning("Not a ggplot input, skipped ...")
    return(x)
  }
  # default arguments
  default_args <- list(
    fill  = "grey50",
    color = "grey50",
    linewidth = 0.3
  )
  args <- purrr::list_modify(default_args, ...)
  # Extract arguments from plot
  p_build <- ggplot_build(x)
  # x axis
  x_discrete <- p_build$layout$panel_params[[1]]$x.sec$is_discrete()
  y_discrete <- p_build$layout$panel_params[[1]]$y.sec$is_discrete()
  # Checkpoint-2: x axis
  if(isTRUE(x_discrete) || isTRUE(y_discrete)) {
    warning(glue::glue(
      "Expect numeric axises: x-axis: {! x_discrete}, y-axis: {! y_discrete}, ",
      "add_genebody_bar() skipped..."
    ))
    return(x)
  }
  # guess genebody range (TSS to TES/PAS)
  x_breaks <- p_build$layout$panel_params[[1]]$x.sec$breaks
  x_labels <- p_build$layout$panel_params[[1]]$x.sec$get_labels()
  gb   <- FALSE # init
  gb_i <- 2
  for(i in seq_len(length(x_labels))) {
    if(grepl("^TSS$", x_labels[i], ignore.case = TRUE)) {
      gb <- grepl("^(TES|PAS)$", x_labels[i+1], ignore.case = TRUE)
      if(gb) gb_i <- i
    }
  }
  if(! isTRUE(gb)) {
    message(glue::glue(
      "Could not find genebody region on x-axis, ",
      "expect TSS-TES, TSS-PAS",
    ))
    return(x)
  }
  gb_range <- x_breaks[gb_i:(gb_i + 1)]
  # y_range <- p_build$layout$panel_params[[1]]$y.range
  y_name   <- p_build$layout$panel_params[[1]]$y.sec$name
  y_limits <- p_build$layout$panel_params[[1]]$y.sec$limits ## skipped
  y_breaks <- p_build$layout$panel_params[[1]]$y.sec$get_breaks() ## skipped
  y_labels <- p_build$layout$panel_params[[1]]$y.sec$get_labels() ## skipped
  # Expand y-axis at y-max by 10%
  #--------------------------#
  # 1. update margin
  y_limits2 <- y_limits * c(1, 1 + margin) * (1 + expand)
  y_mg <- y_limits * c(1, 1 + margin)
  y_gb <- y_mg[2] + diff(y_mg) * expand * c(0.1, 0.4, 0.5, 0.9)
  #--------------------------#
  # 1. update y-axis
  suppressMessages(
    # overwrite old y-axis settings # !!! side-effect ?
    x2 <- x +
      scale_y_continuous(
        name = y_name,
        breaks = y_breaks,
        labels = y_labels,
        limits = y_limits2
      )
  )
  #--------------------------#
  # Locate the box, lines
  rect_h <- list(
    "geom" = "rect",
    "xmin" = gb_range[1],
    "xmax" = gb_range[2],
    "ymin" = y_gb[1],
    "ymax" = y_gb[2]
  )
  seg_v <- list(
    "geom" = "segment",
    "x"    = gb_range[1],
    "xend" = gb_range[1],
    "y"    = y_gb[3],
    "yend" = y_gb[4]
  )
  seg_h <- list(
    "geom" = "segment",
    "x"    = gb_range[1],
    "xend" = gb_range[1] + diff(gb_range) * 0.15,
    "y"    = y_gb[4],
    "yend" = y_gb[4],
    arrow = arrow(length = unit(c(0.04, 0.04), "npc"))
  )
  #------------------------#
  y_name     <- p_build$layout$panel_params[[1]]$y.sec$name
  y_range    <- p_build$layout$panel_params[[1]]$y.range
  y_limits   <- p_build$layout$panel_params[[1]]$y.sec$limits ## skipped
  y_breaks   <- p_build$layout$panel_params[[1]]$y.sec$get_breaks()
  y_labels   <- p_build$layout$panel_params[[1]]$y.sec$get_labels()
  #-------------------------#
  # 1. add genebody grey bar
  p2 <- x2 +
    do.call(annotate, modifyList(rect_h, args))
  # 2. add arrow (horizontal)
  suppressWarnings(
    # suppress, fill for segment()
    p2 <- p2 +
      do.call(annotate, modifyList(seg_v, args)) +
      do.call(annotate, modifyList(seg_h, args))
  )
  p2
}


#' update X and Y axis line style
#'
#' 1. remove original x-axis, y-axis
#' 2. add segment on x-axis, y-axis
#' 3. add extra x_ticks
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
  # defaults arguments
  default_args <- list(
    linewidth = .5,
    linetype  = "solid",
    color     = "grey20"
  )
  args <- purrr::list_modify(default_args, ...)
  # check arguments
  if(! inherits(x, "ggplot")) {
    warning("Not a ggplot input, skipped ...")
    return(x)
  }
  # Extract arguments from plot
  p_build <- ggplot_build(x)
  # Extract x-axis
  x_breaks <- p_build$layout$panel_params[[1]]$x.sec$breaks
  x_breaks <- x_breaks[! is.na(x_breaks)]
  y_breaks <- p_build$layout$panel_params[[1]]$y.sec$breaks
  y_breaks <- y_breaks[! is.na(y_breaks)]
  #-------------------------------------#
  # axis line
  seg_h <- list(
    "geom" = "segment",
    "x"    = c(-Inf, -Inf),
    "xend" = c(tail(x_breaks, 1), -Inf),
    "y"    = c(-Inf, -Inf),
    "yend" = c(-Inf, tail(y_breaks, 1))
  )
  # add lines
  # warning(glue::glue(
  #   "!A-1, linewidth: {args$linewidth}"
  # ))
  x +
    do.call(annotate, modifyList(seg_h, args)) +
    theme(
      axis.line  = element_blank(),
      axis.ticks = element_line(
        linewidth = args$linewidth * .5, # 60% of axis line
        color     = args$color
      )
    )
}


#' update x-axis ticks, labels
#'
#' update x_ticks, x_labels
#' original:
#' 0 40 240;
#' "-2" "TES" "+10"
#' updated:
#' 0 40 80 120 160 200 240;
#' "-2"  "TES" "+2"  "+4"  "+6"  "+8"  "+10"
#'
#' @param x_labels string labels
#' @param x_ticks numeric positions on x-axis
#' @param ...
#'
#' @return
#' x_ticks, x_labels
update_axis_ticks_x <- function(x_labels, x_ticks) {
  # default
  x_end <- tail(x_labels, 1)
  x_end <- gsub("[^0-9]+", "", x_end) %>% as.numeric()
  if(x_end > 2) {
    n_ticks <- floor(x_end / 2) + 1
    tes     <- tail(x_ticks, 2)
    ticks_extra  <- seq(tes[1], tes[2], length.out = n_ticks)
    labels_extra <- seq(0, x_end, length.out = n_ticks)
    labels_extra <- paste0('+', labels_extra)
    # remove first one
    list(
      x_ticks = c(head(x_ticks, -1), tail(ticks_extra, -1)),
      x_labels = c(head(x_labels, -1), tail(labels_extra, -1))
    )
  } else {
    list(
      x_ticks  = x_ticks,
      x_labels = x_labels
    )
  }
}


# reference-point
# x = '/data/yulab/wangming/work/yu_2022/projects/20221229_dlj_ChrRNA_yy218/results/flanking_genes/results/fig2.tss.gs_6k/2.bw2matrix/fig2.ChrRNA_YY218.gs_6k_anti.mat.gz'
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
  # Checkpoint-1: group labels
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
  width    <- n_width * args$width
  height   <- n_height * args$height
  x_sect   <- head(args$x_ticks[-1], -1) # remove first, last element
  #----------------------------------------------------------------------------#
  # Checkpoint-1. basic plot
  # p <- plot_profile_basic(df, !!!args)
  args$df <- df
  p <- do.call(plot_profile_basic, args)
  #----------------------------------------------------------------------------#
  # Checkpoint-2. update x_axis, x_ticks
  if(isTRUE(args$add_x_ticks_extra)) {
    tix <- update_axis_ticks_x(args$x_labels, args$x_ticks)
  } else {
    tix <- list(x_ticks = args$x_ticks, x_labels = args$x_labels)
  }
  p <- p +
    scale_x_continuous(
      name   = args$x_title,
      limits = c(min(tix$x_ticks), max(tix$x_ticks)),
      breaks = tix$x_ticks,
      labels = tix$x_labels,
      expand = expansion(mult = c(0, 0))
    )
  #----------------------------------------------------------------------------#
  # Checkpoint-3. update axis lines and genebody-bar
  p <- update_axis_line(p, linewidth = args$linewidth) #
  if(args$matrix_type == "scale-regions") {
    p <- add_genebody_bar(p)
  }
  # #----------------------------------------------------------------------------#
  # # # Checkpoint-3. update y_axis, y_ticks
  # gar <- update_axis_range_y(df$score, args$y_min, args$y_max)
  # y_limits <- c(gar$y_min, gar$y_max * 1.1) # expand
  # suppressMessages(
  #   p <- p +
  #     scale_y_continuous(
  #       name   = args$y_title,
  #       limits = y_limits,
  #       expand = expansion(mult = c(0.01, 0)) # already expanded
  #     )
  # )
  #----------------------------------------------------------------------------#
  # Checkpoint-4. colors + themes
  n_colors <- length(unique(df[[color_by]]))
  if(inherits(args$colors, "character")) {
    colors <- rep(args$colors, 100)[1:n_colors] #
    p <- p + scale_color_manual(values = colors)
  }
  # add themes
  p <- p +
    theme_classic() +
    theme(
      # panel.border = element_blank(),
      # panel.grid = element_blank(),
      # rect       = element_blank(),
      axis.text  = element_text(color = "black"),
      axis.line  = element_blank(),
      axis.ticks = element_line(
        linewidth = args$linewidth * .5, # 60% of axis line
        color     = args$color
      )
      # axis.ticks = element_line(linewidth = .4, color = "grey30"),
    )
  # message(glue::glue(
  #   "!B-1, axis linewidth: {args$linewidth}"
  # ))
  # 4. Facet by group_labels
  if(length(group_labels) > 1) {
    p <- p +
      facet_wrap(as.formula(paste("~", group_by)), ncol = n_width)
  }
  #----------------------------------------------------------------------------#
  # Checkpoint-6. save plot to file
  if(inherits(filename, "character")) {
    if(file.exists(filename) && ! isTRUE(args$overwrite)) {
      message(glue::glue("file exists: {filename}"))
    } else {
      pdf(NULL) # prevent generating empty file: "Rplot.pdf"
      message(glue::glue(
        "export pdf size, width={width}, height={height}"
      ))
      export::graph2pdf(
        x = p, file = filename, width = width, height = height,
        font = "Arial", bg = "transparent"
      )
      dev.off() # close empty pdf
      rds <- gsub("\\.[a-z]+$", ".rds", filename, perl = T)
      saveRDS(p, file = rds)
    }
  }
  # output
  p
}


#' Generate strand-specific metaplot for regions/reference-point
#'
#' see `plotProfile` of deeptools at
#' https://deeptools.readthedocs.io/en/develop/content/tools/plotProfile.html
#'
#' @param x1 matrix file for sense-strand, see `computeMatrix`
#' @param x2 matrix file for anti-sense-strand,
#' @param ... additional parameters, see arguments in `plot_profile` function
#'
#' @export
plot_profile_ss <- function(x1, x2,  filename = NULL, ...) {
  #----------------------------------------------------------------------------#
  args <- .plot_profile_args(...) # default + ... arguments
  # to variable
  for(name in names(args)) {
    if(rlang::is_empty(name)) next # skip NULL variables
    assign(name, args[[name]]) # items to environment
  }
  #----------------------------------------------------------------------------#
  # loading header from matrix
  ## sens
  hd1 <- .read_matrix_header(x1)
  df1 <- .read_matrix(x1, averageType) # averageType checked
  ## anti
  hd2 <- .read_matrix_header(x2)
  df2 <- .read_matrix(x2, averageType) # averageType checked
  df2$score <- -df2$score
  # check groups
  if(isTRUE(perGroup)) {
    color_by <- "sample_label"
    group_by <- "group_label"
  } else {
    color_by <- "group_label"
    group_by <- "sample_label"
  }
  # 2.2 check color_by, group_by
  cb <- rep(levels(df1[[color_by]]), times = 2)
  cs <- rep(c("_fwd", "_rev"), each = length(cb) / 2)
  cv <- paste0(cb, cs)
  df1[[color_by]] <- paste0(df1[[color_by]], "_fwd")
  df2[[color_by]] <- paste0(df2[[color_by]], "_rev")
  df <- rbind(df1, df2)
  df[[color_by]] <- factor(df[[color_by]], levels = cv)
  # required:
  hs <- c("x_ticks", "x_labels", "x_title", "y_title",
          "smp_labels", "bin_avg_type")
  for(name in hs) {
    if(rlang::is_empty(name)) next
    assign(name, hd1[[name]])
  }
  # update group labels
  if(inherits(group_labels, "character")) {
    if(length(group_labels) == length(hd1$group_labels)) {
      # update
      gl <- setNames(object = group_labels, nm = hd1$group_labels)
      df$group_label <- gl[df$group_label]
    }
  }
  # 3.0 determine fig numbers
  n_plots <- length(unique(df[[group_by]]))
  n_width <- ifelse(n_plots > numPlotsPerRow, numPlotsPerRow, n_plots)
  n_height <- ceiling(n_plots / numPlotsPerRow)
  width <- n_width * width
  height <- n_height * height
  ## 3.1 :basic
  x_sect <- head(hd1$x_ticks[-1], -1) # remove first, last element
  p <- ggplot(df, aes(x, score, color = .data[[color_by]])) +
    geom_vline(xintercept = x_sect, linewidth = .5,
               color = "grey50", linetype = 2) +
    geom_hline(yintercept = 0, linewidth = .5, color = "grey30") +
    geom_line(linewidth = line_size) +
    scale_x_continuous(
      name   = hd1$x_title,
      breaks = hd1$x_ticks,
      labels = hd1$x_labels
    ) +
    facet_wrap(as.formula(paste("~", group_by)), ncol = numPlotsPerRow,
               scales = "free") +
    ggtitle(plot_title)
  # 3.2 :colors
  if(inherits(colors, "character")) {
    color_len <- unique(df1[[color_by]])
    if(length(colors) == length(color_len)) {
      colors_ss <- rep(colors, times = 2)
      p <- p + scale_color_manual(values = colors_ss)
    }
  }
  ## 3.3 :yaxis
  if(inherits(c(y_min, y_max), "numeric")) {
    p <- p + scale_y_continuous(limits = c(-y_max, y_max))
  }
  ## 3.4 :theme
  if(inherits(plot_theme, "character")) {
    if(plot_theme %in% c("few")) {
      p <- p +
        ggthemes::scale_color_few() +
        ggthemes::theme_few()
    }
  }else {
    p <- p +
      ggplot2::theme_bw() +
      theme(panel.grid = element_blank())
  }
  # 4. save to files
  message(glue::glue("save to file: {filename}"))
  if(inherits(filename, "character")) {
    ggsave(filename, p, width = width, height = height, units = units, dpi = dpi)
    rds <- gsub("\\.[a-z]+$", ".rds", filename, perl = T)
    saveRDS(p, file = rds)
  }
  p # return
}


# scale-regions
# x = '/data/yulab/wangming/work/yu_2022/projects/20221229_dlj_ChrRNA_yy218/results/flanking_genes/results/fig1.gs_6k/2.bw2matrix/fig1.ChrRNA_YY218.gs_6k_sens.mat.gz'
# x1 = '/data/yulab/wangming/work/yu_2022/projects/20221229_dlj_ChrRNA_yy218/results/flanking_genes/results/fig1.gs_6k/2.bw2matrix/fig1.ChrRNA_YY218.gs_6k_sens.mat.gz'
# x2 = '/data/yulab/wangming/work/yu_2022/projects/20221229_dlj_ChrRNA_yy218/results/flanking_genes/results/fig1.gs_6k/2.bw2matrix/fig1.ChrRNA_YY218.gs_6k_anti.mat.gz'
# x = "data/config/metaplot/paused_all/fig1A/fig1.A.ChIP_8WG16_60m.metaplot.tss.yaml"

# am = "results/metaplot/paused_all/fig1A/2.bw2matrix/fig1.A.ChIP_8WG16_60m.tes.mat.gz"
# (p <- plot_profile(am, filename = "tmp.cnt.pdf", colors = c("black", "red"), overwrite = T, add_x_ticks_extra = T))

# ay = "data/config/metaplot/paused_all/fig1A/fig1.A.ChIP_8WG16_60m.metaplot.tes.yaml"
# make_metaplot(ay, overwrite = T, add_x_ticks_extra = T)

