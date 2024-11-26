#### Jinxin Meng, 20241027, 20241028 ####

setwd("F:/proj/20240731_BC_IBD_Zhangyn/submit/code_availability/GX_IPA_biosynthesis/")
pacman::p_load(dplyr, tidyr, tibble, purrr, stringr, ggplot2, ggpubr)
library(Biostrings)
library(rlang)

#### plot ####
source("/code/R_func/plot_msa.R")

plot_msa("acdA.msa.trim", width = 60, font_size = 2)
ggsave("acdA.msa.trim.pdf", width = 8, height = 8)

plot_msa("fldC.msa.trim", width = 60, font_size = 2)
ggsave("fldC.msa.trim.pdf", width = 8, height = 8)

plot_msa("fldH.msa.trim", width = 60, font_size = 2)
ggsave("fldH.msa.trim.pdf", width = 8, height = 8)