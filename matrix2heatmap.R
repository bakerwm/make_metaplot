

# setwd("~/work/yu_2021/pcf11_lxj/results/20220110_pas_metaplot/deeptools/results/TTseq_YY122.anti")
suppressPackageStartupMessages(library(dplyr))
# library(readr)
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(hiseqr))
source("/data/biosoft/make_metaplot/matrix2profile.R") # main functions


f <- "results/paused_all/fig1_metaplot/2.bw2matrix/fig1.A.CnT_4H8.genebody.mat.gz"
# 
# df <- load_matrix(f, avg_func = FALSE)
# 
# # sortUsingSample: 
# # sortUsing: 
# # sortUsingSamples: 
# 
# df1 <- df[[1]] %>% dplyr::select(-c(1:3, 5:6)) # save V4 as id
# ma1 <- df1 %>%
#   tibble::column_to_rownames("V4") %>%
#   as.matrix()
# df1x <- df1 %>%
#   dplyr::rename(y = "V4") %>%
#   tidyr::pivot_longer(-1, names_to = "x", values_to = "value") 

# 
# 
# ggplot(df1x, aes(x = x, y = y, fill = value)) +
#   geom_tile() +
#   coord_fixed() +
#   theme(legend.position = "none")
# 
# # pheatmap
# library(RColorBrewer)
# cc <- colorRampPalette(
#   rev(brewer.pal(n = 7, name = "RdYlBu")))
# 
# breaks <- seq(min(unlist(c(ma1, ma1))), 
#              max(unlist(c(ma1, ma1))), length.out=100)
# 
# library(pheatmap)
# p1 = pheatmap(
#   ma1, color=cc(100), breaks=breaks, 
#   show_rownames = F, 
#   show_colnames = F,
#   cluster_cols = F,
#   cluster_rows = F,
#   silent = T
#   )
# 
# pdf('aaa.pdf', width = 4, height = 15)
# print(p1)
# dev.off()








