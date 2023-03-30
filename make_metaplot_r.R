#!/usr/bin/env Rscripts
# Generate metaplot using R function, from matrix
#
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  print("Usage: Rscript make_metaplot_r.R config.yaml <reverse_trand>")
  print("")
  print("Options:")
  print("    config.yaml   The config file for metaplot, in YAML format")
  print("reverse_strand    Reverse strand if required, 1=yes, 0=no (optional)")
  stop("arguments failed")
}

cfg <- args[1] # config
reverse_strand <- args[2] == "1"
suppressPackageStartupMessages(library(hiseqr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(fishualize))
suppressPackageStartupMessages(library(ggthemes))
# source("/data/biosoft/make_metaplot/matrix2profile.R") # main functions
source("/data/biosoft/make_metaplot/make_metaplot.R") # make_metaplot()

p <- make_metaplot(cfg)
# p <- make_metaplot(cfg, overwrite = T, choose_samples = "t2_rep2", colors = c("blue"))

## END
