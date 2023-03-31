
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
  # 1. load config
  j <- tryCatch(
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
  # 2. check-point
  if(! inherits(j, "list")) {
    return(NULL)
  }
  # 3. check args
  j$reverse_strand <- FALSE # default
  args <- .setup_metaplot_params(x = NULL, !!!j)
  args <- purrr::list_modify(args, ...)
  # 4. output
  out_dir <- ifelse(rlang::has_name(j, "out_dir"), j$out_dir, "./")
  prefix  <- ifelse(rlang::has_name(j, "prefix"), j$prefix, "metaplot")
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
    # p <- matrix2profile2(m1, m2, sub_pdf, !!!args)
    p <- plot_profile_ss(m1, m2, sub_pdf, !!!args)
  } else if(file.exists(mat)) {
    # p <- matrix2profile(mat, sub_pdf, !!!args)
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
.plot_profile_args <- function(...) {
  # default arguments
  args <- list(
    averageType   = NULL, # from header, avg_func
    colors        = NULL, # auto
    width         = 6,
    height        = 3,
    units         = "in",
    dpi           = 300,
    # start_label   = "TSS",
    # end_label     = "TES",
    # point_label   = "TSS",  # same as start label
    # group_labels  = NULL, # from header, regionsLabel
    line_size     = 0.6,
    numPlotsPerRow = 3,   # from config.yaml
    overwrite     = FALSE,
    perGroup      = TRUE, # from config.yaml
    plot_theme    = NULL, # default: theme_bw()
    plot_title    = "metaplot",
    # return_data   = FALSE, # return data.frame
    # sample_labels = NULL, # from header, samplesLabel
    sample_list   = NULL,
    x_title       = NULL, # from header
    yMax          = NULL, # from deeptool.plotProfile
    yMin          = NULL, # from deeptool.plotProfile
    # y_max         = NULL, # auto
    # y_min         = NULL, # auto
    y_title       = NULL  # from header
  )
  # update y_min, y_max
  if(is.null(args[["y_min"]]) & inherits(args[['yMin']], "numeric")) {
    args[["y_min"]] <- args[["yMin"]]
  }
  if(is.null(args[["y_max"]]) & inherits(args[["yMax"]], "numeric")) {
    args[["y_max"]] <- args[["yMax"]]
  }
  # dots <- rlang::list2(...)
  purrr::list_modify(args, ...) # items of list to environment
}


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


#' Initate params for matrix2plot
#' combine parameters from matrix and config.yaml file
#'
#' @param x str Path to matrix file
#' @param ... optinoal arguments from config.yaml file
#'
#' @return
#' @export
#'
#' @examples
#'  .setup_param()
.setup_metaplot_params <- function(x = NULL, ...) {
  # Extract params from matrix file
  if(inherits(x, "character")) {
    h1 <- .read_matrix_header(x) # extract from matrix header
    h2 <- rlang::list2(
      sample_labels = h1$smp_labels,
      group_labels  = h1$group_labels,
      averageType   = h1$bin_avg_type,
      bin_size      = h1$bs,
      per_group     = TRUE # default
    )
    h <- purrr::list_modify(h2, !!!h1)
  } else {
    h <- list() # blank
  }
  # optional from config.yaml
  # see deeptools.plotProfile arguments
  ## default
  args_default <- list(
    width         = 6,
    height        = 2.5,
    units         = "in",
    dpi           = 300,
    line_size     = 0.6,
    numPlotsPerRow = 3,   # from config.yaml
    overwrite     = FALSE
  )
  args1 <- purrr::list_modify(args_default, ...)
  # args1 <- .plot_profile_args(...) # default + custome
  ## renmae argument names for downstream plots
  args2 <- rlang::list2(
    add_x_ticks_extra = FALSE,
    bin_size    = args1$binSize,
    plot_title  = args1$plotTitle,
    start_label = args1$startLabel,
    point_label = args1$refPointLabel,
    end_label   = args1$endLabel,
    x_title     = args1$xAxisLabel,
    y_title     = args1$yAxisLabel,
    y_min       = args1$yMin,
    y_max       = args1$yMax,
    per_group   = args1$perGroup,
    group_labels  = args1$regionsLabel,
    sample_labels = args1$samplesLabel,
    bin_avg_type  = args1$averageType
  )
  args <- purrr::list_modify(args2, !!!args1) # modify
  args <- purrr::discard(args, is.null) # remove null arguments
  # fix colors
  if(inherits(args$colors, "character")) {
    args$colors <- strsplit(args$colors, "\\s+", perl = TRUE) %>% unlist
  }
  # updated by optional arguments
  purrr::list_modify(h, !!!args) # overwrite header
  # purrr::list_modify(args, ...)
}


#' Add genebody bar on top of metaplot
#'
#' @param x ggplot2 plot
#' @param ...
#'
#' @return ggplot
#' @export
#'
#' @examples
#' add_genebody_bar(x)
add_genebody_bar <- function(x, ...) {
  # check arguments
  if(! inherits(x, "ggplot")) {
    warning("Not a ggplot input, skipped ...")
    return(x)
  }
  # Extract arguments from plot
  p_build <- ggplot_build(x)
  # x axis
  x_discrete <- p_build$layout$panel_params[[1]]$x.sec$is_discrete()
  x_breaks <- p_build$layout$panel_params[[1]]$x.sec$breaks
  x_labels <- p_build$layout$panel_params[[1]]$x.sec$get_labels()
  # y axis
  y_discrete <- p_build$layout$panel_params[[1]]$y.sec$is_discrete()
  # y_range  <- p_build$layout$panel_params[[1]]$y.range
  y_limits <- p_build$layout$panel_params[[1]]$y.sec$limits ## skipped
  # Checkpoint-2: x axis
  if(isTRUE(x_discrete) || isTRUE(y_discrete)) {
    warning(glue::glue(
      "Expect numeric axises: x-axis: {! x_discrete}, y-axis: {! y_discrete}, ",
      "add_genebody_bar() skipped..."
    ))
    return(x)
  }
  if(! x_labels[2] == "TSS" || ! x_labels[3] %in% c("TES", "PAS")) {
    warning("Could not found TSS--TES/PAS labels on axis")
    return(x)
  }
  #-------------------------#
  # 1. add genebody grey bar
  x_sect <- x_breaks[2:3] # TSS - TES
  y_sect <- y_limits[2] * c(.92, .95)
  p <- x +
    annotate(
      geom  = "rect",
      xmin  = x_sect[1],
      xmax  = x_sect[2],
      ymin  = y_sect[1],
      ymax  = y_sect[2],
      fill  = "grey50"
    )
  # 2. add arrow (horizontal)
  x_seg1 <- c(x_sect[1], x_sect[1] + diff(x_sect) * 0.15)
  y_seg1 <- y_limits[2] * c(.99, .99)
  p <- p +
    annotate(
      "segment",
      x     = x_seg1[1],
      xend  = x_seg1[2],
      y     = y_seg1[1],
      yend  = y_seg1[2],
      color = "grey50",
      linewidth = 0.3,
      arrow = arrow(length = unit(0.02, "npc"))
    )
  # 3. add vertical line
  x_seg2 <- c(x_sect[1], x_sect[1])
  y_seg2 <- y_limits[2] * c(.95, .99)
  p <- p +
    annotate(
      "segment",
      x     = x_seg2[1],
      xend  = x_seg2[2],
      y     = y_seg2[1],
      yend  = y_seg2[2],
      color = "grey50",
      linewidth = 0.3
    )
  # output
  p
}


#' Replace X and Y axises style
#'
#' 1. remove original x-axis, y-axis
#' 2. add segment on x-axis, y-axis
#' 3. add extra x_ticks
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' add_axis_line(x)
add_axis_line <- function(x, ...) {
  dots <- rlang::list2(...)
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
  x_labels <- p_build$layout$panel_params[[1]]$x.sec$get_labels()
  x_limits <- p_build$layout$panel_params[[1]]$x.sec$get_limits()
  x_name   <- p_build$layout$panel_params[[1]]$x.sec$name
  # x_limits <- c(min(x_breaks), max(x_breaks))
  # Extract y-axis
  y_breaks <- p_build$layout$panel_params[[1]]$y.sec$breaks
  y_breaks <- y_breaks[! is.na(y_breaks)]
  # y_labels <- p_build$layout$panel_params[[1]]$y.sec$get_labels()
  # y_limits <- p_build$layout$panel_params[[1]]$y.sec$get_limits()
  # Fix axis-line
  p <- x +
    annotate(
      "segment",
      x     = -Inf,
      xend  = tail(x_breaks, 1),
      y     = -Inf,
      yend  = -Inf,
      linewidth = .5,
      linetype  = 1,
      color = "grey20"
    ) +
    annotate(
      "segment",
      x     = -Inf,
      xend  = -Inf,
      y     = -Inf,
      yend  = tail(y_breaks, 1),
      linewidth = 0.5,
      linetype  = 1,
      color = "grey20"
    ) +
    theme(
      axis.line  = element_blank(),
      axis.ticks = element_line(linewidth = .5, color = "grey20")
    )
  # update x-axis ticks
  if(isTRUE(dots$add_x_ticks_extra)) {
    x_end <- tail(x_labels, 1)
    x_end <- gsub("[^0-9]+", "", x_end) %>% as.numeric()
    if(x_end > 2) {
      n_ticks_extra <- floor(x_end / 2) + 1
      t1 <- seq(rev(x_breaks)[2], rev(x_breaks)[1], length.out = n_ticks_extra)
      t2 <- seq(0, x_end, length.out = n_ticks_extra)
      t2 <- paste0('+', t2)
      # remove first one
      t1 <- t1[-1]
      t2 <- t2[-1]
      x_breaks_extra <- c(rev(rev(x_breaks)[-1]), t1)
      x_labels_extra <- c(rev(rev(x_labels)[-1]), t2)
      # update x-axis
      message("update x-axis")
      suppressMessages(
        p <- p +
          scale_x_continuous(
            name   = x_name,
            limits = x_limits,
            breaks = x_breaks_extra,
            labels = x_labels_extra,
            expand = expansion(mult = c(0, 0))
          )
      )
    }
  }
  p
  # x +
  #   geom_segment(
  #     aes(x = -Inf, xend = tail(x_breaks, 1), y = -Inf, yend = -Inf),
  #     linewidth = .3,
  #     color     = "grey30"
  #   ) +
  #   geom_segment(
  #     aes(x = -Inf, xend = -Inf, y = -Inf, yend = tail(y_breaks, 1)),
  #     linewidth = .3,
  #     color     = "grey30"
  #   ) + theme(
  #     axis.line  = element_blank(),
  #     axis.ticks = element_line(linewidth = .5, color = "black")
  #   )
}

#----------------------------------------------------------------------------#


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
#' see `plotProfile` of deeptools at
#' https://deeptools.readthedocs.io/en/develop/content/tools/plotProfile.html
#'
#' @param x path to the matrix file, output of `computeMatrix` in `deeptools`
#' @param filename path to the plot file
#' @param colors character List of colors to use for the plotted lines, [auto]
#' @param plot_title character the title of the plot, default: [metaplot]
#' @param start_label character label on start point, default [TSS]
#' @param end_label character label on end point, default [TES]
#' @param point_label character label on reference point, default [TSS]
#' @param sample_labels character Labels for the samples plotted [NULL]
#' @param group_labels character labels for regions [NULL]
#' @param averageType character The type of statistic used for profile,
#' options ["mean", "median", "min", "max", "sum", "std"], default ["mean"]
#' @param x_title character Title on x axis, default ["Genomic region (kb)"]
#' @param y_title character Title on y axis, default ["score"]
#' @param y_min numeric Minimum value on y axis, default [auto]
#' @param y_max numeric Maximum value on y axis, default [auto]
#' @param overwrite bool Overwrite the exists plot file
#' @param plot_theme character Add theme to the plot, options ["few"], default [NULL]
#' @param width numeric Set the width of plot, default [7] inches
#' @param height numeric Set the height of plot, default [3] inches
#' @param units character Set the unit for plot, ["cm", "mm", "in", "px"], default [in]
#' @param dpi numeric Plot resolution, see `ggsave()`
#' @param perGroup bool value default TRUE
#' @param numPlotPerRow numeric how many plots each row, default [3]
#' @param return_data bool return the data.frame for plot
#'
#' @importFrom rlang list2 is_empty
#' @importFrom purrr list_modify
#' @importFrom jsonlite parse_json
#' @importFrom ggthemes scale_color_few theme_few
#' @importFrom ggplot2 ggplot geom_vline geom_line scale_x_continuous ggtitle theme scale_color_manual scale_y_continuous theme_bw theme
#'
#' @export
plot_profile <- function(x, filename = NULL, ...) {
  # parse arugments from matrix_header + ...
  args <- .setup_metaplot_params(x, ...) # default + ... arguments
  # loading header from matrix
  df <- .read_matrix(x, bin_avg_type = args$bin_avg_type) # averageType checked
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
  n_width  <- ifelse(n_plots > args$numPlotsPerRow, args$numPlotsPerRow, n_plots)
  n_height <- ceiling(n_plots / args$numPlotsPerRow)
  width    <- n_width * args$width
  height   <- n_height * args$height
  x_sect   <- head(args$x_ticks[-1], -1) # remove first, last element
  #----------------------------------------------------------------------------#
  # 1. Proudfoot lab metaplot
  ## guess y_limits
  tmp <- ggplot(df) +
    geom_line(aes(x, score, color = .data[[color_by]])) +
    scale_y_continuous(expand = expansion(mult = c(0.01, 0.1))) # expand 10%
  tmp_build <- ggplot_build(tmp)
  # tmp_x_breaks <- tmp_build$layout$panel_params[[1]]$x.sec$get_breaks()
  tmp_y_breaks <- tmp_build$layout$panel_params[[1]]$y.sec$get_breaks()
  ## add vlines
  p <- ggplot(df) +
    annotate(
      "segment",
      x     = x_sect,
      xend  = x_sect,
      y     = c(-Inf, -Inf),
      yend  = rep(tail(tmp_y_breaks, 1), 2),
      linewidth = 0.3,
      linetype  = 2,
      color = "grey50"
    ) +
    geom_line(aes(x, score, color = .data[[color_by]]))
  ## fix x-axis
  if(isTRUE(args$add_x_ticks_extra)) {
    # update args$x_ticks
    message("Add extra ticks on x-axis")
  }
  p <- p +
    scale_x_continuous(
      name   = args$x_title,
      limits = c(args$x_ticks[1], tail(args$x_ticks, 1)),
      breaks = args$x_ticks,
      labels = args$x_labels,
      expand = expansion(mult = c(0, 0))
    ) +
    ggtitle(args$plot_title)
  # 2. colors
  n_colors <- length(unique(df[[color_by]]))
  if(inherits(args$colors, "character")) {
    colors <- rep(args$colors, 100)[1:n_colors] #
    p <- p + scale_color_manual(values = colors)
  }
  # 3. update y-axis limits
  ## update: y-axis limits
  # extract y axis range from plot
  p_build  <- ggplot2::ggplot_build(p)
  # y_limits <- p_build$layout$panel_params[[1]]$y.sec$get_limits() # limits
  y_limits  <- p_build$layout$panel_params[[1]]$y.range # range
  y_limits[2] <- y_limits[2] * 1.1 # expand y-axis upper 10%
  if(inherits(args$y_min, "numeric")) {
    y_limits[1] <- args$y_min
  }
  if(inherits(args$y_max, "numeric")) {
    y_limits[2] <- args$y_max
  }
  p <- p +
    scale_y_continuous(
      name   = args$y_title,
      limits = y_limits,
      expand = expansion(mult = c(0.01, 0)) # already expanded
    )
  # 4. Facet by group_labels
  if(length(group_labels) > 1) {
    p <- p +
      facet_wrap(as.formula(paste("~", group_by)), ncol = n_width)
  }
  # 5. theme
  p <- p +
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      rect       = element_blank(),
      axis.line  = element_line(linewidth = .8, color = "black"),
      axis.ticks = element_line(linewidth = .8, color = "black"),
      axis.text  = element_text(color = "black")
    )
  # 6. add genebody bar: for scale-regions
  p <- add_axis_line(p, add_x_ticks_extra = args$add_x_ticks_extra) #
  if(args$matrix_type == "scale-regions") {
    p <- add_genebody_bar(p)
  }
  # 7. save to files
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
  p # return
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

# x = "results/metaplot/paused_all/fig1A/2.bw2matrix/fig1.A.ChIP_8WG16_60m.genebody.mat.gz"
# (p <- plot_profile(x, filename = "tmp.cnt.pdf", colors = c("black", "red"), add_x_ticks_extra = T))

x = "data/config/metaplot/paused_all/fig1A/fig1.A.ChIP_8WG16_60m.metaplot.genebody.yaml"
make_metaplot(x, overwrite = T, add_x_ticks_extra = T)

