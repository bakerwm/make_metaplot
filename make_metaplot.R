
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
  dots <- rlang::list2(...)
  args <- purrr:::list_modify(j, !!!dots)
  args$group_labels  <- args$regionsLabel
  args$sample_labesl <- args$samplesLabel
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
    # filename      = NULL,
    colors        = NULL, # auto
    plot_title    = "metaplot",
    start_label   = "TSS",
    end_label     = "TES",
    point_label   = "TSS",  # same as start label
    sample_labels = NULL, # from header, samplesLabel
    group_labels  = NULL, # from header, regionsLabel
    averageType   = NULL, # from header, avg_func
    x_title       = NULL, # from header
    y_title       = NULL, # from header
    y_min         = NULL, # auto
    y_max         = NULL, # auto
    yMin          = NULL, # from deeptool.plotProfile
    yMax          = NULL, # from deeptool.plotProfile
    overwrite     = FALSE,
    plot_theme    = NULL, # default: theme_bw()
    line_size     = 0.6,
    width         = 5,
    height        = 3,
    units         = "in",
    dpi           = 300,
    sample_list   = NULL,
    perGroup      = TRUE, # from config.yaml
    numPlotsPerRow = 3,   # from config.yaml
    return_data   = FALSE # return data.frame
  )
  # update y_min, y_max
  if(is.null(args[["y_min"]]) & inherits(args[['yMin']], "numeric")) {
    args[["y_min"]] <- args[["yMin"]]
  }
  if(is.null(args[["y_max"]]) & inherits(args[["yMax"]], "numeric")) {
    args[["y_max"]] <- args[["yMax"]]
  }
  dots <- rlang::list2(...)
  purrr::list_modify(args, !!!dots) # items of list to environment
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
      message(glue::glue("Unknown avg_func: {x}, expect one of [{aa}]"))
    }
  } else {
    message(glue::glue("Unknown avg_func: {x}, expect character"))
  }
  return(FALSE)
}


# scale-regions
# x = '/data/yulab/wangming/work/yu_2022/projects/20221229_dlj_ChrRNA_yy218/results/flanking_genes/results/fig1.gs_6k/2.bw2matrix/fig1.ChrRNA_YY218.gs_6k_sens.mat.gz'
# x1 = '/data/yulab/wangming/work/yu_2022/projects/20221229_dlj_ChrRNA_yy218/results/flanking_genes/results/fig1.gs_6k/2.bw2matrix/fig1.ChrRNA_YY218.gs_6k_sens.mat.gz'
# x2 = '/data/yulab/wangming/work/yu_2022/projects/20221229_dlj_ChrRNA_yy218/results/flanking_genes/results/fig1.gs_6k/2.bw2matrix/fig1.ChrRNA_YY218.gs_6k_anti.mat.gz'

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
.read_matrix <- function(x, avg_func = NULL) {
  # 1. loading header
  header <- .read_matrix_header(x)
  if(is.null(avg_func)) avg_func <- header$bin_avg_type #
  if(! .is_valid_avg_func(avg_func)) return(NULL) # exception
  # 2. loading matrix
  message(glue::glue("Loading matrix from file: {basename(x)}"))
  df <- read.delim(x, header = FALSE, sep = "\t", comment.char = "@")
  df1 <- subset(df, select = -c(1:6)) # remove bed6 record
  hd_col <- max(header$smp_bound)
  hd_row <- max(header$group_bound)
  # check matrix size: col-row
  if(! all(dim(df1) == c(hd_row, hd_col))) {
    msg <- glue::glue("Error, expect matrix size [{dim(df)[1]}, {dim(df)[2]}]",
                      ", actual matrix size is [{hd_row}, {hd_col}]")
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
        score = apply(dj, 2, match.fun(avg_func)), # avg values by group
        sample_label = header$smp_labels[i],
        group_label  = j$group_label[1]
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
  #       score = apply(dfg, 2, match.fun(avg_func)),
  #       sample_label = header$smp_labels[i],
  #       group_label  = header$group_labels[j]
  #     )
  #   }) %>%
  #     dplyr::bind_rows()
  # }) %>%
  #   dplyr::bind_rows()
  # 4. fix order
  df2$sample_label <- factor(df2$sample_label, levels = header$smp_labels)
  df2$group_label  <- factor(df2$group_label, levels = header$group_labels)
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
  #----------------------------------------------------------------------------#
  args <- .plot_profile_args(...) # default + ... arguments
  # to variable
  for(name in names(args)) {
    if(rlang::is_empty(name)) next # skip NULL variables
    assign(name, args[[name]]) # items to environment
  }
  # message(glue::glue("1. {args[['line_size']]}, 2. {line_size}"))
  #----------------------------------------------------------------------------#
  # loading header from matrix
  header <- .read_matrix_header(x)
  df <- .read_matrix(x, avg_func = averageType) # averageType checked
  # required:
  hs <- c("x_ticks", "x_labels", "x_title", "y_title", "smp_labels",
          "bin_avg_type")
  for(name in hs) {
    if(rlang::is_empty(name)) next
    assign(name, header[[name]])
  }
  # update group labels
  if(inherits(group_labels, "character")) {
    if(length(group_labels) == length(header$group_labels)) {
      # update
      gl <- setNames(object = group_labels, nm = header$group_labels)
      df$group_label <- gl[df$group_label]
    }
  }
  # check groups
  if(isTRUE(perGroup)) {
    color_by <- "sample_label"
    group_by <- "group_label"
  } else {
    color_by <- "group_label"
    group_by <- "sample_label"
  }
  # columns: x, score, sample_label, group_label
  if(inherits(group_labels, "character")) {
    if(length(group_labels) == length(levels(df[[group_by]]))) {
      group_df <- setNames(group_labels, nm = levels(df[[group_by]]))
    }
  }
  # fix y-axis
  ## !!! to-do
  # 3.0 determine fig numbers
  n_plots <- length(unique(df[[group_by]]))
  n_width <- ifelse(n_plots > numPlotsPerRow, numPlotsPerRow, n_plots)
  n_height <- ceiling(n_plots / numPlotsPerRow)
  width <- n_width * width
  height <- n_height * height
  # 3.1 :basic
  x_sect <- head(header$x_ticks[-1], -1) # remove first, last element
  p <- ggplot(df, aes(x, score, color = .data[[color_by]])) +
    geom_vline(xintercept = x_sect, linewidth = .5,
               color = "grey50", linetype = 2) +
    geom_line(linewidth = line_size) +
    scale_x_continuous(
      name   = header$x_title,
      breaks = header$x_ticks,
      labels = header$x_labels
    ) +
    facet_wrap(as.formula(paste("~", group_by)), ncol = numPlotsPerRow) +
    ggtitle(plot_title)
  ## 3.2 :colors
  n_colors <- length(unique(df[[color_by]]))
  if(inherits(colors, "character")) {
    colors <- rep(colors, 100)[1:n_colors] #
    p <- p + scale_color_manual(values = colors)
  }
  ## 3.3 :yaxis
  if(inherits(c(y_min, y_max), "numeric")) {
    p <- p + scale_y_continuous(limits = c(y_min, y_max))
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
  if(inherits(filename, "character")) {
    if(file.exists(filename) && ! isTRUE(overwrite)) {
      message(glue::glue("file exists: {filename}"))
    } else {
      pdf(NULL) # prevent generating empty file: "Rplot.pdf"
      export::graph2pdf(x = p, file = filename, width = width, height = height,
                        font = "Arial", bg = "transparent")
      rds <- gsub("\\.[a-z]+$", ".rds", filename, perl = T)
      saveRDS(p, file = rds)
    # ggsave(filename, p, width = width, height = height, units = units, dpi = dpi)
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
  df1 <- .read_matrix(x1, avg_func = averageType) # averageType checked
  ## anti
  hd2 <- .read_matrix_header(x2)
  df2 <- .read_matrix(x2, avg_func = averageType) # averageType checked
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

