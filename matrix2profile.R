

# setwd("~/work/yu_2021/pcf11_lxj/results/20220110_pas_metaplot/deeptools/results/TTseq_YY122.anti")
suppressPackageStartupMessages(library(dplyr))
# library(readr)
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(hiseqr))



#------------------------------------------------------------------------------#
# make_metaplot (profile)
# save to "${out_dir}/5.matrix2profile/${prefix}_plotProfile.pdf
#' @param x character path to the YAML file
#' @param reverse_strand bool default strandness is dUTP, if is not, use `TRUE` 
#' to reverse the strand
#'
#' @export
make_metaplot <- function(x, reverse_strand = FALSE, ...) {
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
  dots <- rlang::list2(...)
  args <- purrr:::update_list(j, !!!dots)
  # args <- purrr::discard(args, is.null) # why?
  # fix colors
  args$colors <- strsplit(args$colors, "\\s+", perl = TRUE) %>% unlist
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
    dir.create(sub_dir)
  }
  if(file.exists(mat_sens) & file.exists(mat_anti)) {
    if(isTRUE(reverse_strand)) {
      m1 <- mat_anti
      m2 <- mat_sens
    } else {
      m1 <- mat_sens
      m2 <- mat_anti
    }
    p <- matrix2profile2(m1, m2, sub_pdf, !!!args)
  } else if(file.exists(mat)) {
    p <- matrix2profile(mat, sub_pdf, !!!args)
    # x <- 1
  } else {
    warning(glue::glue("Could not found matrix file: {mat}"))
    return(NULL)
  }
  return(NULL)
}


#------------------------------------------------------------------------------#


#' Generate metaplot for regions/reference-point
#'
#' see `plotProfile` of deeptools at
#' https://deeptools.readthedocs.io/en/develop/content/tools/plotProfile.html
#'
#' @param x path to the matrix file, output of `computeMatrix` in `deeptools`
#' @param filename path to the plot file
#'
#' @param start_label character label on start point, default [TSS]
#' @param end_label character label on end point, default [TES]
#' @param point_label character label on reference point, default [TSS]
#' @param sample_labels character Labels for the samples plotted
#' @param plotTitle character the title of the plot, default: [metaplot]
#' @param colors character List of colors to use for the plotted lines, [auto]
#' @param x_title character Title on x axis, default ["Genomic region (kb)"]
#' @param y_title character Title on y axis, default ["score"]
#' @param y_min numeric Minimum value on y axis, default [auto]
#' @param y_max numeric Maximum value on y axis, default [auto]
#' @param return_data bool return the data.frame for plot
#' @param overwrite bool Overwrite the exists plot file
#' @param avg_func character The type of statistic used for profile,
#' options ["mean", "median", "min", "max", "sum", "std"], default ["mean"]
#' @param plot_theme character Add theme to the plot, options ["few"], default [NULL]
#' @param width numeric Set the width of plot, default [7] inches
#' @param height numeric Set the height of plot, default [3] inches
#' @param units character Set the unit for plot, ["cm", "mm", "in", "px"], default [in]
#' @param dpi numeric Plot resolution, see `ggsave()`
#'
#' @importFrom rlang list2 is_empty
#' @importFrom purrr update_list
#' @importFrom jsonlite parse_json
#' @importFrom ggthemes scale_color_few theme_few
#' @importFrom ggplot2 ggplot geom_vline geom_line scale_x_continuous ggtitle theme scale_color_manual scale_y_continuous theme_bw theme
#'
#' @export
matrix2profile <- function(x, filename = NULL, ...) {
  # default arguments
  args   <- .metaplot_args(...)
  header <- load_matrix_header(x)
  args   <- purrr::update_list(header, !!!args)
  # to variable
  for(name in names(args)) {
    if(rlang::is_empty(name)) next
    assign(name, args[[name]])
  }
  # load matrix data.frame
  df <- load_matrix(x, args$avg_func) # avg_func: default args
  # check arguments
  if(isTRUE(args$return_data)) { # return_data: default args
    return(c(
      header,
      list(df = df, plotTitle = args$plotTitle, colors = args$colors)
    ))
  }
  #----------------------------------------------------------------------------#
  # 3. plot
  ## 3.1 :basic
  ## fix y-max:
  if(is.numeric(args$yMin)) {
  # if(inherits(args$yMin, "numeric")) {
    y_min <- args$yMin
  } else if(min(df$score) > 0) {
    y_min <- 0
  } else {
    y_min <- min(df$score) * 1.1
  }
  if(is.numeric(args$yMax)) {
  # if(inherits(args$yMax, "numeric")) {
    y_max <- args$yMax
  } else if(max(df$score) < 0) {
    y_max <- 0
  } else {
    y_max <- max(df$score) * 1.1
  }
  # :main
  p <- ggplot(df, aes(x, score, color = label)) +
    geom_vline(xintercept = args$x_sect, size = .3, color = "grey50", linetype = 2) +
    geom_line(size = .5) +
    scale_x_continuous(
      name   = args$xAxisLabel,
      breaks = args$x_ticks,
      labels = args$x_labels
    ) +
    scale_y_continuous(
      name = args$yAxisLabel,
      limits = c(y_min, y_max)
    ) +
    # ylab(args$yAxisLabel) +
    ggtitle(args$plotTitle)
  ## 3.2 :colors
  if(inherits(colors, "character")) {
    colors <- rep(colors, 100)[1:length(sl)] #
    p <- p + scale_color_manual(values = args$colors)
  }
  ## 3.3 :yaxis
  if(inherits(c(args$yMin, args$yMax), "numeric")) {
    p <- p + 
      scale_y_continuous(
        name   = args$yAxisLabel,
        limits = c(args$yMin, args$yMax)
      )
  }
  ## 3.4 :theme
  if(inherits(args$plot_theme, "character")) {
    if(plot_theme %in% c("few")) {
      p <- p +
        ggthemes::scale_color_few() +
        ggthemes::theme_few() + 
        theme(
          text = element_text(family = "Arial")
        )
    }
  }else {
    p <- p +
      ggplot2::theme_bw() +
      theme(
        panel.grid = element_blank(),
        axis.title = element_text(size = 10, color = "black"),
        axis.text  = element_text(size = 10, color = "black"))
  }
  #----------------------------------------------------------------------------#
  # 4. save to files
  if(inherits(filename, "character")) {
    # ggsave(filename, plot = p, width = width, height = height, units = units, dpi = dpi)
    pdf(NULL) # prevent generating empty file: "Rplot.pdf"
    export::graph2pdf(x = p, file = filename, width = width, height = height,
                      font = "Arial", bg = "transparent")
    rds <- gsub("\\.[a-z]+$", ".rds", filename, perl = T)
    saveRDS(p, file = rds)
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
matrix2profile2 <- function(x1, x2, filename = NULL, ...) {
  dots <- rlang::list2(...)
  args <- list(
    colors    = NULL,
    ss_colors = FALSE
  )
  args <- purrr::update_list(args, !!!dots)
  # 0. :fix colors for sens/anti
  if(isTRUE(args$ss_colors) & inherits(args$colors, "character")) {
    cc <- split(args$colors, f = c("sens", "anti"))
    args_sens <- purrr::update_list(args, colors = cc$sens)
    args_anti <- purrr::update_list(args, colors = cc$anti)
  } else {
    args_sens <- args_anti <- args # copy
  }
  # 1. :basic plot
  p1 <- matrix2profile(x1, filename = NULL, !!!args_sens) # sens
  p2 <- matrix2profile(x2, filename = NULL, !!!args_anti) # anti
  #----------------------------------------------------------------------------#
  # 1. :auto-define y-axis, range
  ymax <- max(p1$data$score, p2$data$score)
  ymin <- min(p1$data$score, p2$data$score)
  if(inherits(args[["ymax"]], "numeric")) {
    ymax <- args[["ymax"]]
  }
  # 2. :fix sense (on top)
  p1 <- p1 +
    scale_y_continuous(
      name   = "Sense",
      limits = c(ymin, ymax),
      # breaks = scales::pretty_breaks()(ymax, 5)
    ) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      panel.border = element_rect(color = "black", size = .5),
      plot.margin  = margin(b = 0, unit = "pt"))
  # 3. :fix anti-sense (on bottom)
  p2 <- p2 +
    scale_y_continuous(
      name   = "Antisense",
      limits = c(ymax, ymin),
      # breaks = scales::pretty_breaks()(ymax, 5),
      trans  = scales::reverse_trans()) +
    theme(
      # legend.position = "none",
      plot.title = element_blank(),
      panel.border = element_rect(color = "black", size = .5),
      plot.margin = margin(t = 0, unit = "pt"))
  # # 4. :fix color
  # # add colors
  # if(inherits(fish_color, "character")) {
  #   if(fish_color %in% fishualize::fishcolors[, "option"]) {
  #     p1 <- p1 + fishualize::scale_color_fish_d(option = fish_color)
  #     p2 <- p2 + fishualize::scale_color_fish_d(option = fish_color)
  #   }
  # }
  # 5. :merge sense+anti
  p12 <- patchwork::wrap_plots(p1, p2, ncol = 1, guides = "collect")
  #----------------------------------------------------------------------------#
  # 4. save to files
  if(inherits(filename, "character")) {
    # ggsave(filename, plot = p12, width = width, height = height, units = units, dpi = dpi)
    pdf(NULL) # prevent generating empty file: "Rplot.pdf"\
    px <- export::graph2pdf(x = p12, file = filename, width = args$width,
                            height = args$height * 1.2,
                            font = "Arial", bg = "transparent")
    rds <- gsub("\\.[a-z]+$", ".rds", filename, perl = T)
    saveRDS(p12, file = rds)
  }
  p12 # return
}


#------------------------------------------------------------------------------#
# read yaml file
# expect output files
# 2.bw2matrix
# 3.matrxi2profile
# 4.matrix2heatmap
# 5.matrix2profile_R

#' Generate default arguments for metaplot
#' 
#' default arguments for all profile
.metaplot_args <- function(...) {
  # default arguments
  args <- list(
    avg_func      = "mean",
    colors        = NULL, # auto
    dpi           = 300,
    end_label     = "TES",
    height        = 3,
    linesAtTickMarks = FALSE,
    out_dir       = "results/fig1_metaplot",
    overwrite     = FALSE,
    perGroup      = TRUE,
    plotHeight    = 8,
    ploWidth      = 12,
    plot_theme    = NULL, # default: theme_bw()
    plotTitle     = "metaplot",
    plotType      = "lines",
    point_label   = "TSS",  # same as start label
    prefix        = NULL,
    refPointLabel = "TSS",
    referencePoint = "TSS",
    return_data   = FALSE, # return data.frame
    sample_labels = NULL, # default
    sample_list   = NULL,
    startLabel    = "TSS",
    units         = "in",
    width         = 5,
    xAxisLabel    = "Genomic region (kb)",
    yAxisLabel    = "Score",
    yMax          = NULL, # auto
    yMin          = NULL, # auto
    y_title       = "score"
  )
  dots <- rlang::list2(...)
  purrr::update_list(args, !!!dots)
}


#' see `computeMatrix` of deeptools at
#' https://deeptools.readthedocs.io/en/develop/content/tools/plotProfile.html
#'
#' @param x Character path to the matrix file
#' @param header_only bool return the head only
#'
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
#'
#' @export
load_matrix_header <- function(x) {
  h <- readLines(x, n = 1)
  j <- jsonlite::parse_json(gsub("^@", "", h))
  ## 1.1 sample labels (bw files)
  sr1 <- unlist(j$sample_boundaries)
  sr  <- tail(sr1, -1) - head(sr1, -1)
  sl  <- unlist(j$sample_labels) # !!! save order
  ss  <- rep(sl, sr) # labels, global variable
  #-------------------------#
  # # update sample_labeles, from global_env
  # if(inherits(sample_labels, "character")) {
  #   if(length(s) == length(sample_labels)) {
  #     s <- sample_labels
  #   }
  # }
  #-------------------------#
  ## 1.2 x-axis, labels
  ## to-do: unscaled 5 prime: !!!
  ## up, TSS, TES, down
  us <- ifelse(rlang::has_name(j, "upstream"), j$upstream[[1]], 0)
  gb <- ifelse(rlang::has_name(j, "body"), j$body[[1]], 1000)
  ds <- ifelse(rlang::has_name(j, "downstream"), j$downstream[[1]], 0)
  bs <- ifelse(rlang::has_name(j, "bin size"), j$`bin size`[[1]], 10)
  u5 <- ifelse(rlang::has_name(j, "unscaled 5 prime"), j$`unscaled 5 prime`[1], 0)
  u3 <- ifelse(rlang::has_name(j, "unscaled 3 prime"), j$`unscaled 3 prime`[1], 0)
  ## lables on x axis
  usl <- paste0("-", round(us / 1000, 1))
  dsl <- paste0("+", round(ds / 1000, 1))
  ## x-tick labels
  if(rlang::has_name(j, "ref point")) {
    ref <- j$`ref point`[[1]]
  } else {
    ref <- NULL
  }
  if(is.null(ref)) {
    start_label <- "TSS"
    end_label   <- "TES"
    x_labels <- c(usl, start_label, end_label, dsl)
    x_list      <- c(us, gb, ds) / bs
    # x_ticks  <- Reduce(f = "+", x = c(us, gb, ds) / bs, accumulate = TRUE)
  } else {
    start_label <- NULL
    end_label   <- NULL
    x_labels <- c(usl, ref, dsl)
    x_list      <- c(us, ds) / bs
    # x_ticks  <- Reduce(f = "+", x = c(us, ds) / bs, accumulate = TRUE)
  }
  # ref <- switch(is.null(a)+1,"notNullHihi",NULL)
  x_ticks  <- Reduce(f = "+", x = x_list, accumulate = TRUE) # add 0
  x_ticks  <- c(0, x_ticks)
  ## 1.3 Axis, labels
  x_axis   <- unlist(lapply(sr, seq))
  # x_ticks  <- Reduce(f = "+", x = c(us, gb, ds) / bs, accumulate = TRUE) # add 0
  # x_ticks  <- c(0, x_ticks)
  # x_labels <- c(usl, start_label, end_label, dsl)
  x_title  <- "Genomic region (kb)"
  y_title  <- "Mean of score"
  x_sect   <- head(x_ticks, -1) %>% tail(-1) # TSS, TES
  # print(paste0("!AAAA-1", f, x_ticks))
  # return values
  list(
    sl = sl,   # sample_labels
    ss = ss,  # list of sample_labels
    x_axis   = x_axis,
    x_ticks  = x_ticks,
    x_labels = x_labels,
    xAxisLabel  = x_title,
    yAxisLabel  = y_title,
    x_sect   = x_sect
  )
  # x_sect, x_title, x_ticks, x_labels, plotTitle, s, colors, y_min, y_max, plot_theme
}


#----------------------------------------------------------------------------#
# 2. load matrix
# load_matrix <- function(x, avg_func = "mean") {
#   ## 2.1 load file
#   message("Reading matrix...")
#   df1 <- read.delim(x, header = FALSE, sep = "\t", comment.char = "@")
#   ma  <- df1 %>%
#     dplyr::select(-c(1:6)) %>%
#     as.matrix
#   ## 2.2 load header
#   header <- load_matrix_header(x)
#   ## 2.2 meta data
#   # "mean", "median", "min", "max", "sum" and "std"; default: [mean]
#   score <- apply(ma, 2, match.fun(avg_func))
#   tibble::tibble(score = score, label = header$ss, x = header$x_axis) %>%
#     dplyr::mutate(label = factor(label, levels = header$sl))
# }


#' read matrix
#' 
#' @export
load_matrix <- function(x, avg_func = "mean") {
  ## 2.1 load file
  message("Reading matrix...")
  df1 <- read.delim(x, header = FALSE, sep = "\t", comment.char = "@")
  ma  <- df1 %>%
    dplyr::select(-c(1:6)) %>%
    as.matrix
  ## 2.2 load header
  header <- load_matrix_header(x)
  ## 2.2 meta data
  # "mean", "median", "min", "max", "sum" and "std"; default: [mean]
  func <- c("mean", "median", "min", "max", "sum", "std")
  if(avg_func %in% func) {
    score <- apply(ma, 2, match.fun(avg_func))
    tibble::tibble(score = score, label = header$ss, x = header$x_axis) %>%
      dplyr::mutate(label = factor(label, levels = header$sl))
  } else {
    # matrix: 1-6
    dx <- sapply(header$sl, function(i) {
      col <- grep(i, header$ss)
      col <- c(1:6, col + 6)
      dplyr::select(df1, all_of(col))
    }, USE.NAMES = TRUE, simplify = FALSE)
  }
}






#------------------------------------------------------------------------------#












#' 
#' #' Generate metaplot for regions/reference-point
#' #'
#' #' see `plotProfile` of deeptools at
#' #' https://deeptools.readthedocs.io/en/develop/content/tools/plotProfile.html
#' #'
#' #' @param x path to the matrix file, output of `computeMatrix` in `deeptools`
#' #' @param filename path to the plot file
#' #'
#' #' @param start_label character label on start point, default [TSS]
#' #' @param end_label character label on end point, default [TES]
#' #' @param point_label character label on reference point, default [TSS]
#' #' @param sample_labels character Labels for the samples plotted
#' #' @param plotTitle character the title of the plot, default: [metaplot]
#' #' @param colors character List of colors to use for the plotted lines, [auto]
#' #' @param x_title character Title on x axis, default ["Genomic region (kb)"]
#' #' @param y_title character Title on y axis, default ["score"]
#' #' @param y_min numeric Minimum value on y axis, default [auto]
#' #' @param y_max numeric Maximum value on y axis, default [auto]
#' #' @param return_data bool return the data.frame for plot
#' #' @param overwrite bool Overwrite the exists plot file
#' #' @param avg_func character The type of statistic used for profile,
#' #' options ["mean", "median", "min", "max", "sum", "std"], default ["mean"]
#' #' @param plot_theme character Add theme to the plot, options ["few"], default [NULL]
#' #' @param width numeric Set the width of plot, default [7] inches
#' #' @param height numeric Set the height of plot, default [3] inches
#' #' @param units character Set the unit for plot, ["cm", "mm", "in", "px"], default [in]
#' #' @param dpi numeric Plot resolution, see `ggsave()`
#' #'
#' #' @importFrom rlang list2 is_empty
#' #' @importFrom purrr update_list
#' #' @importFrom jsonlite parse_json
#' #' @importFrom ggthemes scale_color_few theme_few
#' #' @importFrom ggplot2 ggplot geom_vline geom_line scale_x_continuous ggtitle theme scale_color_manual scale_y_continuous theme_bw theme
#' #'
#' #' @export
#' plot_profile <- function(x, filename = NULL, ...) {
#'   # default arguments
#'   args <- list(
#'     start_label   = "TSS",
#'     end_label     = "TES",
#'     point_label   = "TSS",  # same as start label
#'     sample_labels = NULL, # default
#'     plotTitle    = "metaplot",
#'     colors        = NULL, # auto
#'     x_title       = "Genomic region (kb)",
#'     y_title       = "score",
#'     y_min         = NULL, # auto
#'     y_max         = NULL, # auto
#'     return_data     = FALSE, # return data.frame
#'     overwrite     = FALSE,
#'     avg_func      = "mean",
#'     plot_theme    = NULL, # default: theme_bw()
#'     width         = 7,
#'     height        = 3,
#'     units         = "in",
#'     dpi           = 300,
#'     sample_list   = NULL
#'   )
#'   dots <- rlang::list2(...)
#'   args <- purrr::update_list(args, !!!dots)
#'   # to variable
#'   for(name in names(args)) {
#'     if(rlang::is_empty(name)) next
#'     assign(name, args[[name]])
#'   }
#'   #----------------------------------------------------------------------------#
#'   # load matrix header
#'   header <- .read_matrix_header(x)
#'   for(name in names(header)) {
#'     if(rlang::is_empty(name)) next
#'     assign(name, header[[name]])
#'   }
#'   # load matrix data.frame
#'   df <- .read_matrix(x)
#'   # check arguments
#'   if(isTRUE(return_data)) {
#'     return(c(
#'       header,
#'       list(df = df, plotTitle = plotTitle, colors = colors)
#'     ))
#'   }
#'   #----------------------------------------------------------------------------#
#'   # 3. plot
#'   ## 3.1 :basic
#'   p <- ggplot(df, aes(x, score, color = label)) +
#'     geom_vline(xintercept = x_sect, size = .5, color = "grey50", linetype = 2) +
#'     geom_line(size = .7) +
#'     scale_x_continuous(
#'       name   = x_title,
#'       breaks = x_ticks,
#'       labels = x_labels
#'     ) +
#'     ggtitle(plotTitle)
#'   ## 3.2 :colors
#'   if(inherits(colors, "character")) {
#'     colors <- rep(colors, 100)[1:length(sl)] #
#'     p <- p + scale_color_manual(values = colors)
#'   }
#'   ## 3.3 :yaxis
#'   if(inherits(c(y_min, y_max), "numeric")) {
#'     p <- p + scale_y_continuous(limits = c(y_min, y_max))
#'   }
#'   ## 3.4 :theme
#'   if(inherits(plot_theme, "character")) {
#'     if(plot_theme %in% c("few")) {
#'       p <- p +
#'         ggthemes::scale_color_few() +
#'         ggthemes::theme_few()
#'     }
#'   }else {
#'     p <- p +
#'       ggplot2::theme_bw() +
#'       theme(panel.grid = element_blank())
#'   }
#'   #----------------------------------------------------------------------------#
#'   # 4. save to files
#'   if(inherits(filename, "character")) {
#'     ggsave(filename, width = width, height = height, units = units, dpi = dpi)
#'   }
#'   p
#' }
#' 
#' 
#' #' Generate strand-specific metaplot for regions/reference-point
#' #'
#' #' see `plotProfile` of deeptools at
#' #' https://deeptools.readthedocs.io/en/develop/content/tools/plotProfile.html
#' #'
#' #' @param x1 matrix file for sense-strand, see `computeMatrix`
#' #' @param x2 matrix file for anti-sense-strand,
#' #' @param ... additional parameters, see arguments in `plot_profile` function
#' #'
#' #' @export
#' plot_profile_ss <- function(x1, x2, ...) {
#'   dots <- rlang::list2(...)
#'   args <- list(
#'     colors    = NULL,
#'     ss_colors = FALSE
#'   )
#'   args <- purrr::update_list(args, !!!dots)
#'   # 0. :fix colors for sens/anti
#'   if(isTRUE(args$ss_colors) & inherits(args$colors, "character")) {
#'     cc <- split(args$colors, f = c("sens", "anti"))
#'     args_sens <- purrr::update_list(args, colors = cc$sens)
#'     args_anti <- purrr::update_list(args, colors = cc$anti)
#'   } else {
#'     args_sens <- args_anti <- args # copy
#'   }
#'   # 1. :basic plot
#'   p1 <- plot_profile(x1, filename = NULL, !!!args_sens) # sens
#'   p2 <- plot_profile(x2, filename = NULL, !!!args_anti) # anti
#'   #----------------------------------------------------------------------------#
#'   # 1. :auto-define y-axis, range
#'   ymax <- max(p1$data$score, p2$data$score)
#'   ymin <- min(p1$data$score, p2$data$score)
#'   if(inherits(args[["ymax"]], "numeric")) {
#'     ymax <- args[["ymax"]]
#'   }
#'   # 2. :fix sense (on top)
#'   p1 <- p1 +
#'     scale_y_continuous(
#'       name   = "Sense",
#'       limits = c(ymin, ymax),
#'       # breaks = scales::pretty_breaks()(ymax, 5)
#'     ) +
#'     theme(
#'       axis.title.x = element_blank(),
#'       axis.text.x  = element_blank(),
#'       axis.ticks.x = element_blank(),
#'       panel.border = element_rect(color = "black", size = .5),
#'       plot.margin  = margin(b = 0, unit = "pt"))
#'   # 3. :fix anti-sense (on bottom)
#'   p2 <- p2 +
#'     scale_y_continuous(
#'       name   = "Antisense",
#'       limits = c(ymax, ymin),
#'       # breaks = scales::pretty_breaks()(ymax, 5),
#'       trans  = scales::reverse_trans()) +
#'     theme(
#'       # legend.position = "none",
#'       plot.title = element_blank(),
#'       panel.border = element_rect(color = "black", size = .5),
#'       plot.margin = margin(t = 0, unit = "pt"))
#'   # # 4. :fix color
#'   # # add colors
#'   # if(inherits(fish_color, "character")) {
#'   #   if(fish_color %in% fishualize::fishcolors[, "option"]) {
#'   #     p1 <- p1 + fishualize::scale_color_fish_d(option = fish_color)
#'   #     p2 <- p2 + fishualize::scale_color_fish_d(option = fish_color)
#'   #   }
#'   # }
#'   # 5. :merge sense+anti
#'   patchwork::wrap_plots(p1, p2, ncol = 1, guides = "collect")
#' }
#' 
#' 
#' #' Generate metaplot for regions/reference-point
#' #'
#' #' see `plotProfile` of deeptools at
#' #' https://deeptools.readthedocs.io/en/develop/content/tools/plotProfile.html
#' #'
#' #' @param m path to the matrix file, output of `computeMatrix` in `deeptools`
#' #' @param filename path to the plot file
#' #'
#' #' @param start_label character label on start point, default [TSS]
#' #' @param end_label character label on end point, default [TES]
#' #' @param point_label character label on reference point, default [TSS]
#' #' @param sample_labels character Labels for the samples plotted
#' #' @param plotTitle character the title of the plot, default: [metaplot]
#' #' @param colors character List of colors to use for the plotted lines, [auto]
#' #' @param x_title character Title on x axis, default ["Genomic region (kb)"]
#' #' @param y_title character Title on y axis, default ["score"]
#' #' @param y_min numeric Minimum value on y axis, default [auto]
#' #' @param y_max numeric Maximum value on y axis, default [auto]
#' #' @param return_data bool return the data.frame for plot
#' #' @param overwrite bool Overwrite the exists plot file
#' #' @param avg_func character The type of statistic used for profile,
#' #' options ["mean", "median", "min", "max", "sum", "std"], default ["mean"]
#' #' @param plot_theme character Add theme to the plot, options ["few"], default [NULL]
#' #' @param width numeric Set the width of plot, default [7] inches
#' #' @param height numeric Set the height of plot, default [3] inches
#' #' @param units character Set the unit for plot, ["cm", "mm", "in", "px"], default [in]
#' #' @param dpi numeric Plot resolution, see `ggsave()`
#' #'
#' #' @importFrom rlang list2 is_empty
#' #' @importFrom purrr update_list
#' #' @importFrom jsonlite parse_json
#' #' @importFrom ggthemes scale_color_few theme_few
#' #' @importFrom ggplot2 ggplot geom_vline geom_line scale_x_continuous ggtitle theme scale_color_manual scale_y_continuous theme_bw theme
#' #'
#' #' @export
#' plot_profile2 <- function(m, filename = NULL, ...) {
#'   # default arguments
#'   args <- list(
#'     start_label   = "TSS",
#'     end_label     = "TES",
#'     point_label   = "TSS",  # same as start label
#'     sample_labels = NULL, # default
#'     plotTitle    = "metaplot",
#'     colors        = NULL, # auto
#'     x_title       = "Genomic region (kb)",
#'     y_title       = "score",
#'     y_min         = NULL, # auto
#'     y_max         = NULL, # auto
#'     return_data     = FALSE, # return data.frame
#'     overwrite     = FALSE,
#'     avg_func      = "mean",
#'     plot_theme    = NULL, # default: theme_bw()
#'     width         = 7,
#'     height        = 3,
#'     units         = "in",
#'     dpi           = 300,
#'     sample_list   = NULL
#'   )
#'   dots <- rlang::list2(...)
#'   args <- purrr::update_list(args, !!!dots)
#'   # to variable
#'   for(name in names(args)) {
#'     if(rlang::is_empty(name)) next
#'     assign(name, args[[name]])
#'   }
#'   #----------------------------------------------------------------------------#
#'   # 1. check arguments
#'   ## :m
#'   if(missing(m)) {
#'     warning("m=, missing")
#'     return(NULL)
#'   }
#'   if(! inherits(m, "character")) {
#'     warning(glue::glue("illegal args, m, expect character, got {class(f)}"))
#'     return(NULL)
#'   }
#'   if(length(m) > 1) {
#'     warning("multiple files detected, use the first one only, m = ")
#'     m <- m[1]
#'   }
#'   if(! file.exists(m)) {
#'     warning(glue::glue("file not exists, m = {m}"))
#'     return(NULL)
#'   }
#'   ## :filename (pdf/png) + return_data
#'   if(isTRUE(return_data)) {
#'     # filename <- tempfile("plot", fileext = ".png")
#'     filename <- NULL
#'   } else {
#'     # if(missing(filename)) {
#'     #   warning("filename=, missing")
#'     #   return(NULL)
#'     # }
#'     # if(! inherits(filename, "character")) {
#'     #   warning(glue::glue(
#'     #     "illegal args, filename, expect character, got {class(filename)}"
#'     #   ))
#'     #   return(NULL)
#'     # }
#'     if(inherits(filename, "character")) {
#'       if(length(filename) > 1) {
#'         warning("multiple files detected, use the first one only, filename = ")
#'         filename <- filename[1]
#'       }
#'       if(! dir.exists(dirname(filename))) {
#'         dir.create(dirname(filename), recursive = TRUE)
#'       }
#'       if(file.exists(filename) & ! isTRUE(overwrite)) {
#'         message("output file exists, set overwrite=TRUE to re-create the file")
#'         return(NULL)
#'       }
#'     }
#'   }
#'   ## :avg_func
#'   func_list <- c("mean", "median", "min", "max", "sum", "std")
#'   if(! inherits(avg_func, "character")) {
#'     warning(glue::glue(
#'       "unknown avg_func=, expect character, got {class(avg_func)}"
#'     ))
#'     return(NULL)
#'   }
#'   if(! avg_func %in% func_list) {
#'     func_str <- paste(func_list, collapse = ", ")
#'     warning(glue::glue(
#'       "unknown avg_func={avg_func}, ",
#'       "choose from [{func_str}]"
#'     ))
#'     return(NULL)
#'   }
#'   
#'   #----------------------------------------------------------------------------#
#'   # 2. parse header
#'   #' parse the header from matrix file
#'   #' return list
#'   .parse_header <- function(f) {
#'     h <- readLines(f, n = 1)
#'     j <- jsonlite::parse_json(gsub("^@", "", h))
#'     ## 1.1 labels (bw files)
#'     x  <- unlist(j$sample_boundaries)
#'     x1 <- tail(x, -1) - head(x, -1)
#'     s  <- unlist(j$sample_labels) # !!! save order
#'     #-------------------------#
#'     # update sample_labeles, from global_env
#'     if(inherits(sample_labels, "character")) {
#'       if(length(s) == length(sample_labels)) {
#'         s <- sample_labels
#'       }
#'     }
#'     #-------------------------#
#'     ss <- rep(s, x1) # labels, global variable
#'     ## 1.2 x-axis, labels
#'     ## to-do: unscaled 5 prime: !!!
#'     ## up, TSS, TES, down
#'     us <- ifelse(rlang::has_name(j, "upstream"), j$upstream[[1]], 0)
#'     gb <- ifelse(rlang::has_name(j, "body"), j$body[[1]], 1000)
#'     ds <- ifelse(rlang::has_name(j, "downstream"), j$downstream[[1]], 0)
#'     bs <- ifelse(rlang::has_name(j, "bin size"), j$`bin size`[[1]], 10)
#'     u5 <- ifelse(rlang::has_name(j, "unscaled 5 prime"), j$`unscaled 5 prime`[[1]], 0)
#'     u3 <- ifelse(rlang::has_name(j, "unscaled 3 prime"), j$`unscaled 3 prime`[[1]], 0)
#'     ## lables on x axis
#'     usl <- paste0("-", round(us / 1000, 1))
#'     dsl <- paste0("+", round(ds / 1000, 1))
#'     ## 1.3 Axis, labels
#'     x_axis   <- unlist(lapply(x1, seq))
#'     x_ticks  <- Reduce(f = "+", x = c(us, gb, ds) / bs, accumulate = TRUE) # add 0
#'     x_ticks  <- c(0, x_ticks)
#'     x_labels <- c(usl, start_label, end_label, dsl)
#'     x_title  <- "Genomic region (kb)"
#'     y_title  <- "Mean of score"
#'     x_sect   <- head(x_ticks, -1) %>% tail(-1) # TSS, TES
#'     # print(paste0("!AAAA-1", f, x_ticks))
#'     # return values
#'     list(
#'       s  = s,   # sample_labels
#'       ss = ss,  # list of sample_labels
#'       x_axis   = x_axis,
#'       x_ticks  = x_ticks,
#'       x_labels = x_labels,
#'       x_title  = x_title,
#'       y_title  = y_title,
#'       x_sect   = x_sect
#'     )
#'     # x_sect, x_title, x_ticks, x_labels, plotTitle, s, colors, y_min, y_max, plot_theme
#'   }
#'   header <- tryCatch(
#'     {
#'       .parse_header(m) # global variables
#'     },
#'     error = function(cond) {
#'       message("Failed to parse header")
#'       message(cond)
#'       # Choose a return value in case of error
#'       return(NULL)
#'     }
#'   )
#'   # set to local variable
#'   for(name in names(header)) {
#'     if(rlang::is_empty(name)) next
#'     assign(name, header[[name]])
#'   }
#'   
#'   #----------------------------------------------------------------------------#
#'   # 2. load matrix
#'   .load_matrix <- function(f) {
#'     ## 2.1 load file
#'     df1 <- read.delim(m, header = FALSE, sep = "\t", comment.char = "@")
#'     ma  <- df1 %>%
#'       dplyr::select(-c(1:6)) %>%
#'       as.matrix
#'     ## 2.2 meta data
#'     # "mean", "median", "min", "max", "sum" and "std"; default: [mean]
#'     score <- apply(ma, 2, match.fun(avg_func))
#'     tibble::tibble(score = score, label = ss, x = x_axis) %>%
#'       dplyr::mutate(label = factor(label, levels = s))
#'   }
#'   df <- tryCatch(
#'     {
#'       .load_matrix(m) # skip 1st line
#'     },
#'     error = function(cond) {
#'       message("Failed to parse matrix")
#'       message(cond)
#'       # Choose a return value in case of error
#'       return(NULL)
#'     }
#'   )
#'   #--#
#'   if(isTRUE(return_data)) {
#'     return(c(
#'       header,
#'       list(df = df, plotTitle = plotTitle, colors = colors)
#'     ))
#'   }
#'   
#'   #----------------------------------------------------------------------------#
#'   # 3. plot
#'   ## 3.1 basic plot
#'   if(inherits(sample_list, "character")) {
#'     df <- dplyr::filter(df, label %in% sample_list)
#'   }
#'   if(nrow(df) == 0) {
#'     warning("no samples found")
#'     return(NULL)
#'   }
#'   p <- ggplot(df, aes(x, score, color = label)) +
#'     geom_vline(xintercept = x_sect, size = .5, color = "grey50", linetype = 2) +
#'     geom_line(size = .7) +
#'     scale_x_continuous(
#'       name   = x_title,
#'       breaks = x_ticks,
#'       labels = x_labels
#'     ) +
#'     ggtitle(plotTitle)
#'   
#'   ## 3.2 :colors
#'   if(inherits(colors, "character")) {
#'     # valid colors !!!
#'     colors <- rep(colors, 100)[1:length(s)] #
#'     p <- p +
#'       scale_color_manual(values = colors)
#'   }
#'   
#'   ## 3.3 :yaxis
#'   if(inherits(c(y_min,y_max), "numeric")) {
#'     p +
#'       scale_y_continuous(limits = c(y_min, y_max))
#'   }
#'   
#'   ## 3.4 :theme
#'   if(inherits(plot_theme, "character")) {
#'     if(plot_theme %in% c("few")) {
#'       p <- p +
#'         ggthemes::scale_color_few() +
#'         ggthemes::theme_few()
#'     }
#'   }else {
#'     p <- p +
#'       ggplot2::theme_bw() +
#'       theme(panel.grid = element_blank())
#'   }
#'   
#'   #----------------------------------------------------------------------------#
#'   # 4. save to files
#'   if(inherits(filename, "character")) {
#'     ggsave(filename, width = width, height = height, units = units, dpi = dpi)
#'   }
#'   p
#' }
#' 
#' 
#' combine_profile <- function(sens, anti, fish_color = NULL, ...) {
#'   p1 <- plot_profile(sens, ...)  # sens
#'   p2 <- plot_profile(anti, ...)  # anti
#'   ymax <- max(p1$data$score, p2$data$score)
#'   ymin <- min(p1$data$score, p2$data$score)
#'   # add colors
#'   if(inherits(fish_color, "character")) {
#'     if(fish_color %in% fishualize::fishcolors[, "option"]) {
#'       p1 <- p1 +
#'         fishualize::scale_color_fish_d(option = fish_color)
#'       p2 <- p2 +
#'         fishualize::scale_color_fish_d(option = fish_color)
#'     }
#'   }
#'   p1 <- p1 +
#'     scale_y_continuous(
#'       name   = "Sense",
#'       limits = c(ymin, ymax),
#'       breaks = scales::pretty_breaks()(ymax, 5)
#'     ) +
#'     theme(
#'       axis.title.x = element_blank(),
#'       axis.text.x  = element_blank(),
#'       axis.ticks.x = element_blank(),
#'       panel.border = element_rect(color = "black", size = .5),
#'       plot.margin = margin(b = 0, unit = "pt"))
#'   p2 <- p2 +
#'     scale_y_continuous(
#'       name   = "Antisense",
#'       limits = c(ymax, ymin),
#'       breaks = scales::pretty_breaks()(ymax, 5),
#'       trans  = scales::reverse_trans()) +
#'     theme(
#'       legend.position = "none",
#'       plot.title = element_blank(),
#'       panel.border = element_rect(size = .5),
#'       plot.margin = margin(t = 0, unit = "pt"))
#'   patchwork::wrap_plots(p1, p2, ncol = 1, guides = "collect")
#' }
#' 
#' 
#' combine_profile2 <- function(f1, f2, colors = NULL, ...) {
#'   # load data
#'   df1 <- plot_profile(f1, return_data = T, ...)
#'   df2 <- plot_profile(f2, return_data = T, ...)
#'   
#'   # fix anti/sens
#'   df <- dplyr::bind_rows(
#'     dplyr::mutate(df1$df, label = paste0(label, ":sens")),
#'     dplyr::mutate(df2$df, label = paste0(label, ":anti"), score = -score)
#'   )
#'   
#'   # fix levels
#'   v <- levels(df1$df$label)
#'   v <- c(paste0(v, ":sens"), paste0(v, ":anti"))
#'   df <- df %>%
#'     dplyr::mutate(label = factor(label, levels = v))
#'   
#'   # common variables
#'   p <- ggplot(df, aes(x, score, color = label)) +
#'     geom_vline(xintercept = df1$x_sect, size = .5, color = "grey50", linetype = 2) +
#'     geom_hline(yintercept = 0, size = .5, color = "grey30") +
#'     geom_line(size = .7) +
#'     scale_x_continuous(
#'       name   = df1$x_title,
#'       breaks = df1$x_ticks,
#'       labels = df1$x_labels
#'     ) +
#'     ggtitle(df1$plotTitle) +
#'     ggplot2::theme_bw() +
#'     theme(panel.grid = element_blank()) +
#'     annotate("text",
#'              x = max(df1$x_ticks),
#'              y = c(-0.03, 0.03),
#'              label = c("antisense", "sense"),
#'              hjust = 1)
#'   
#'   # fix colors
#'   if(inherits(colors, "character")) {
#'     colors2 <- rep(colors, 2)
#'     p <- p + scale_color_manual(values = colors2)
#'   }
#'   p
#' }
#' 
#' 





#------------------------------------------------------------------------------#
# tmp func
#
#' 
#' #------------------------------------------------------------------------------#
#' #' see `computeMatrix` of deeptools at
#' #' https://deeptools.readthedocs.io/en/develop/content/tools/plotProfile.html
#' #'
#' #' @param x path to the matrix file
#' #'
#' .read_matrix_header <- function(x) {
#'   h <- readLines(x, n = 1)
#'   j <- jsonlite::parse_json(gsub("^@", "", h))
#'   ## 1.1 sample labels (bw files)
#'   sr1 <- unlist(j$sample_boundaries)
#'   sr  <- tail(sr1, -1) - head(sr1, -1)
#'   sl  <- unlist(j$sample_labels) # !!! save order
#'   ss  <- rep(sl, sr) # labels, global variable
#'   #-------------------------#
#'   # # update sample_labeles, from global_env
#'   # if(inherits(sample_labels, "character")) {
#'   #   if(length(s) == length(sample_labels)) {
#'   #     s <- sample_labels
#'   #   }
#'   # }
#'   #-------------------------#
#'   ## 1.2 x-axis, labels
#'   ## to-do: unscaled 5 prime: !!!
#'   ## up, TSS, TES, down
#'   us <- ifelse(rlang::has_name(j, "upstream"), j$upstream[[1]], 0)
#'   gb <- ifelse(rlang::has_name(j, "body"), j$body[[1]], 1000)
#'   ds <- ifelse(rlang::has_name(j, "downstream"), j$downstream[[1]], 0)
#'   bs <- ifelse(rlang::has_name(j, "bin size"), j$`bin size`[[1]], 10)
#'   u5 <- ifelse(rlang::has_name(j, "unscaled 5 prime"), j$`unscaled 5 prime`[1], 0)
#'   u3 <- ifelse(rlang::has_name(j, "unscaled 3 prime"), j$`unscaled 3 prime`[1], 0)
#'   ## lables on x axis
#'   usl <- paste0("-", round(us / 1000, 1))
#'   dsl <- paste0("+", round(ds / 1000, 1))
#'   ## x-tick labels
#'   if(rlang::has_name(j, "ref point")) {
#'     ref <- j$`ref point`[[1]]
#'   } else {
#'     ref <- NULL
#'   }
#'   if(is.null(ref)) {
#'     start_label <- "TSS"
#'     end_label   <- "TES"
#'     x_labels <- c(usl, start_label, end_label, dsl)
#'     x_list      <- c(us, gb, ds) / bs
#'     # x_ticks  <- Reduce(f = "+", x = c(us, gb, ds) / bs, accumulate = TRUE)
#'   } else {
#'     start_label <- NULL
#'     end_label   <- NULL
#'     x_labels <- c(usl, ref, dsl)
#'     x_list      <- c(us, ds) / bs
#'     # x_ticks  <- Reduce(f = "+", x = c(us, ds) / bs, accumulate = TRUE)
#'   }
#'   # ref <- switch(is.null(a)+1,"notNullHihi",NULL)
#'   x_ticks  <- Reduce(f = "+", x = x_list, accumulate = TRUE) # add 0
#'   x_ticks  <- c(0, x_ticks)
#'   ## 1.3 Axis, labels
#'   x_axis   <- unlist(lapply(sr, seq))
#'   # x_ticks  <- Reduce(f = "+", x = c(us, gb, ds) / bs, accumulate = TRUE) # add 0
#'   # x_ticks  <- c(0, x_ticks)
#'   # x_labels <- c(usl, start_label, end_label, dsl)
#'   x_title  <- "Genomic region (kb)"
#'   y_title  <- "Mean of score"
#'   x_sect   <- head(x_ticks, -1) %>% tail(-1) # TSS, TES
#'   # print(paste0("!AAAA-1", f, x_ticks))
#'   # return values
#'   list(
#'     sl = sl,   # sample_labels
#'     ss = ss,  # list of sample_labels
#'     x_axis   = x_axis,
#'     x_ticks  = x_ticks,
#'     x_labels = x_labels,
#'     x_title  = x_title,
#'     y_title  = y_title,
#'     x_sect   = x_sect
#'   )
#'   # x_sect, x_title, x_ticks, x_labels, plotTitle, s, colors, y_min, y_max, plot_theme
#' }
#' 
#' 
#' #----------------------------------------------------------------------------#
#' # 2. load matrix
#' .read_matrix <- function(x, avg_func="mean") {
#'   ## 2.1 load file
#'   df1 <- read.delim(x, header = FALSE, sep = "\t", comment.char = "@")
#'   ma  <- df1 %>%
#'     dplyr::select(-c(1:6)) %>%
#'     as.matrix
#'   ## 2.2 load header
#'   header <- .read_matrix_header(x)
#'   ## 2.2 meta data
#'   # "mean", "median", "min", "max", "sum" and "std"; default: [mean]
#'   score <- apply(ma, 2, match.fun(avg_func))
#'   tibble::tibble(score = score, label = header$ss, x = header$x_axis) %>%
#'     dplyr::mutate(label = factor(label, levels = header$sl))
#' }
#' 



























#------------------------------------------------------------------------------#
# deprecated

# f <- "TTseq_YY122.anti.mat.gz"
# a <- plot_profile(f, "abc.png", overwrite = T, return_data = T)

# plot_profile_fun <- function(df, ...) {
#   dots <- rlang::list2(...)
#   # required
#   # x_sect, x_title, x_ticks, x_labels, plotTitle, s, colors, y_min, y_max, plot_theme
#   if(! inherits(df, "data.frame")) {
#     warning("require data.frame, df=, failed.")
#     return(NULL)
#   }
#   r_cols <- c("score", "label", "x")
#   if(! all(r_cols %in% names(df))) {
#     r_str <- paste(r_cols, collapse = ", ")
#     warning(glue::glue("missing required columns, see [{r_str}]"))
#     return(NULL)
#   }
#
#   # to variable
#   for(name in names(args)) {
#     if(rlang::is_empty(name)) next
#     assign(name, args[[name]])
#   }
#   # 3.1 basic plot
#   p <- ggplot(df, aes(x, score, color = label)) +
#     geom_vline(xintercept = x_sect, size = .5, color = "grey50", linetype = 2) +
#     geom_line(size = .7) +
#     scale_x_continuous(
#       name   = x_title,
#       breaks = x_ticks,
#       labels = x_labels
#     ) +
#     ggtitle(plotTitle)
#
#   ## 3.2 :colors
#   if(inherits(colors, "character")) {
#     # valid colors !!!
#     colors <- rep(colors, 100)[1:length(s)] #
#     p <- p +
#       scale_color_manual(values = colors)
#   }
#
#   ## 3.3 :yaxis
#   if(inherits(c(y_min,y_max), "numeric")) {
#     p +
#       scale_y_continuous(limits = c(y_min, y_max))
#   }
#
#   ## 3.4 :theme
#   if(inherits(plot_theme, "character")) {
#     if(plot_theme %in% c("few")) {
#       p <- p +
#         ggthemes::scale_color_few() +
#         ggthemes::theme_few()
#     }
#   }else {
#     p <- p +
#       ggplot2::theme_bw() +
#       theme(panel.grid = element_blank())
#   }
#
#   p
# }









# #
# f <- "TTseq_YY122.anti.mat.gz"
# a <- plot_profile(f, "abc.png", overwrite = T)
# a <- plot_profile(f, "abc.png", colors = c("black", "red"), width = 10, overwrite = T, return_data = TRUE)

























# #------------------------------------------------------------------------------#
# # 0. arguments
# sl <- "TSS"
# el <- "TES"

# #------------------------------------------------------------------------------#
# # 1. load header
# h <- readLines(f, n = 1)
# j <- jsonlite::parse_json(gsub("^@", "", h))
# j$upstream[[1]]
# ## 1.1 labels (bw files)
# x  <- unlist(j$sample_boundaries)
# x1 <- tail(x, -1) - head(x, -1)
# s  <- unlist(j$sample_labels)
# ss <- rep(s, x1) # labels
#
# ## 1.2 x-axis, labels
# ## to-do: unscaled 5 prime: !!!
# ##
# ## up, TSS, TES, down
# us <- ifelse(rlang::has_name(j, "upstream"), j$upstream[[1]], 0)
# gb <- ifelse(rlang::has_name(j, "body"), j$body[[1]], 1000)
# ds <- ifelse(rlang::has_name(j, "downstream"), j$downstream[[1]], 0)
# bs <- ifelse(rlang::has_name(j, "bin size"), j$`bin size`[[1]], 10)
# u5 <- ifelse(rlang::has_name(j, "unscaled 5 prime"), j$`unscaled 5 prime`[[1]], 0)
# u3 <- ifelse(rlang::has_name(j, "unscaled 3 prime"), j$`unscaled 3 prime`[[1]], 0)
# ## lables on x axis
# usl <- paste0("-", round(us / 1000, 1))
# dsl <- paste0("+", round(ds / 1000, 1))
#
# ## 1.3 Axis, labels
# x_ticks  <- Reduce(f = "+", x = c(us, gb, ds) / bs, accumulate = TRUE) # add 0
# x_ticks  <- c(0, x_ticks)
# x_labels <- c(usl, sl, el, dsl)
# x_title  <- "Genomic region (kb)"
# y_title  <- "Mean of score"
# x_sect   <- head(x_ticks, -1) %>% tail(-1) # TSS, TES

# #------------------------------------------------------------------------------#
# # 2. load matrix
# df <- read.delim(f, header = FALSE, sep = "\t", comment.char = "@")
# ma <- df %>%
#   dplyr::select(-c(1:6)) %>%
#   as.matrix
#
# ## 2.1 average matrix
# # "mean", "median", "min", "max", "sum" and "std"; default: [mean]
# score <- apply(ma, 2, mean)
# df1 <- tibble::tibble(score = score, label = ss, x = x_axis)
#
# #------------------------------------------------------------------------------#
# # 3. plot
# ## 3.1 basic plot
# p <- ggplot(df1, aes(x, score, color = label)) +
#   geom_vline(xintercept = x_sect, size = .5, color = "grey50", linetype = 2) +
#   geom_line(size = .7) +
#   scale_x_continuous(
#     name   = x_title,
#     breaks = x_ticks,
#     labels = x_labels
#   )
#
# ## 3.2 add theme
# p <- p +
#   ggthemes::scale_color_few() +
#   ggthemes::theme_few()


# f <- "TTseq_YY122.anti.mat.gz"
# df <- readr::read_delim(f, "\t", col_names = F, comment = "@")
#
# df1 <- df %>%
#   dplyr::select(-c(1:6)) %>%
#   colMeans() %>%
#   as.data.frame() %>%
#   dplyr::mutate(sample = rep(c("s1", "s2"), each = 170),
#                 pos    = rep(c(1:170), 2))
#
#
# ggplot(df1, aes(pos, ., color = sample)) +
#   geom_line()
#
#
# s <- readr::read_lines(f, n_max = 1)
# s <- gsub("^@", "", s)
# d <- jsonlite::parse_json(s)




#genes:4413
#downstream:10000       upstream:2000   body:5000       bin size:100    unscaled 5 prime:0      unscaled 3 prime:0


# ## load matrix, gz file
# f2 <- "TTseq_YY122.anti.mat.tab"
# df2 <- readr::read_delim(f2, delim = "\t", col_names = T, comment = "#")
