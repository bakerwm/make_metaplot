

# Description
#
# Plot pausing index data


# ECDF plot
read_pi <- function(f, gene = NULL) {
  name <- gsub("CnT_mHaploid_0_1F_|AM_", "", basename(dirname(f)))
  df1 <- readr::read_delim(f, "\t", col_names = F, show_col_types = F)
  df2 <- df1 %>%
    # dplyr::filter(X1 %in% gene) %>%
    dplyr::select(X1, X16) %>%
    dplyr::arrange(X16) %>%
    dplyr::mutate(
      pidx = X16,
      name = name
    ) %>%
    dplyr::filter(pidx >= 1) %>%
    dplyr::mutate(log2pidx = log2(pidx))
  if(inherits(gene, "character")) {
    df2 <- dplyr::filter(df2, X1 %in% gene)
  }
  # fraction
  df2 <- df2 %>%
    dplyr::mutate(
      rank = row_number(),
      frac = rank / length(df2)
    )
  df2
}

plot_pi <- function(df, title = "Pausing Index") {
  ggplot(df, aes(log2(pidx), frac, color = name)) +
    geom_line() +
    xlim(0, 6) +
    scale_color_manual(values = cc) +
    xlab("Log2(Pausing index)") +
    ylab("Fraction") +
    ggtitle(title) +
    ggthemes::theme_clean() +
    theme(
      legend.position = c(0.7, 0.2),
      legend.title = element_blank(),
      legend.box.background = element_blank(),
      legend.background = element_blank(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
}

prj_dir <- "/data/yulab/wangming/work/yu_2023/projects/20230615_dlj_PRO_yy337/demo2/results/pause_index"
pi_dir  <- file.path(prj_dir, "results/3.pausing_index/count_files")
pattern <- "pausing_index.rpk.txt"
f_list  <- list.files(pi_dir, pattern, recursive = T, full.names = T)
df <- lapply(f_list, read_pi) %>%
  dplyr::bind_rows()

## subset: m13
df1 <- dplyr::filter(df, grepl("m13", name))
df1$name <- hiseqr::deseq_sanitize_str(df1$name, n_max = 10)
cc <- c("black", "red")
p  <- plot_pi(df1, "Pausing Index (m13)")
pdf_out <- "pausing_index.m13.pdf"
ggsave(pdf_out, p, width = 4, height = 3)

## subset: m24
df2 <- dplyr::filter(df, grepl("m24", name))
df2$name <- hiseqr::deseq_sanitize_str(df2$name, n_max = 10)
cc <- c("black", "red")
p  <- plot_pi(df2, "Pausing Index (m24)")
pdf_out <- "pausing_index.m24.pdf"
ggsave(pdf_out, p, width = 4, height = 3)


# dev.off()
# ggsave("pcf11_ChIP_pausing_index.CnT_YY136.pdf", p, width = 4.5, height = 3.5)







