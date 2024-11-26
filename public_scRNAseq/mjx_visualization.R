#### Jinxin Meng, 20240711, 20240817 ####

setwd("F:/proj/20240731_BC_IBD_Zhangyn/submit/code_availability/public_scRNAseq/")
pacman::p_load(dplyr, tidyr, tibble, purrr, stringr, ggplot2, ggpubr)
source("/code/R_func/difference_analysis.R")
library(Seurat)
library(harmony)
library(cowplot)
library(data.table)

#### 细胞分群 ####
seurat_obj <- readRDS("seurat.rds")

DimPlot(seurat_obj, reduction = "tsne") +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = c(1, 1),
        legend.justification = c(1, 1))

allmarkers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top_markers <- allmarkers %>%
  group_by(cluster) %>%
  top_n(n = 8, wt = avg_log2FC) %>% 
  pull()
top_markers[grepl("Hist1h1b", top_markers)] <- "Lgr5"

DotPlot(seurat_obj, features = unique(top_markers), cols = c("#abd9e9", "#fdae61"), dot.scale = 2.5) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 8, hjust = 1)) +
  NoLegend()

#### co-expression Lgr5 Hmgcs2 ####

data <- FetchData(seurat_obj, vars = c("Lgr5", "Hmgcs2")) %>% 
  add_column(cells = seurat_obj@active.ident) %>% 
  rownames_to_column("name") %>% 
  left_join(rownames_to_column(data.frame(seurat_obj@reductions$tsne@cell.embeddings), "name"), by = "name") %>% 
  mutate(label = ifelse(cells == "CBCs" & Hmgcs2 > 0, "in_CBC",
                        ifelse(cells != "CBCs" & Hmgcs2 > 0, "in_other", "not_exp")))

ggplot(data, aes(x = tSNE_1, y = tSNE_2, color = label)) +
  geom_point(size = .3) +
  scale_color_manual(values = c("#f8766d", "#b2e0d2", "grey85")) +
  theme_classic() +
  theme(legend.position = c(1, 1),
        legend.justification = c(1, 1))

data.frame(x = filter(data, label != "not_exp") %>% pull(cells) %>% table,
           y = data %>% pull(cells) %>% table %>% as.numeric()) %>% 
  rename(cells = 1, x = 2) %>% 
  mutate(prec = x / y * 100) %>% 
  ggbarplot("cells", "prec", fill = "cells", legend = "none", xlab = "", width = .6,
            ylab = "Precentage (%)",  sort.val = "asc", sort.by.groups = F, rotate = T,
            title = "Expression pattern of Hmgcs2 in ISCs and other cells") +
  scale_fill_brewer(palette = "Set3") +
  theme(aspect.ratio = 1)
