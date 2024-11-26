#### Jinxin Meng, 20240711, 20240805 ####

setwd("F:/proj/20240731_BC_IBD_Zhangyn/submit/code_availability/scRNAseq/")
pacman::p_load(dplyr, tidyr, tibble, purrr, stringr, ggplot2, ggpubr)
library(Seurat)
library(clusterProfiler)
library(harmony)
library(data.table)
library(cowplot)

#### 细胞组成 ####
seurat_obj <- readRDS("seurats.rds")
seurat_obj$orig.ident <- factor(seurat_obj$orig.ident, c("PBS", "BC"))

data <- seurat_obj@meta.data %>% 
  select(group = orig.ident) %>% 
  add_column(cells = Idents(seurat_obj)) %>% 
  group_by(group, cells) %>% 
  summarise(value = n()) %>% 
  group_modify(~.x %>% mutate(prec = value / sum(value) * 100)) %>% 
  ungroup

ggbarplot(data, "group", "prec", fill = "cells", xlab = "", ylab = "Cell precentage (%)",
          legend = "right") +
  scale_fill_brewer(palette = "Set3") +
  theme(aspect.ratio = 1.8)
# ggsave("cell_prec.pdf", width = 4, height = 4)

#### Hopx+ cell 分析 ####
data <- seurat_obj@meta.data %>% 
  select(group = orig.ident) %>% 
  add_column(Hopx = FetchData(seurat_obj, vars = "Hopx") %>% unlist(use.names = F)) %>% 
  mutate(Hopx = ifelse(Hopx == 0, "neg", "pos"))

chisq.test(table(data$Hopx, seurat_obj$orig.ident))

plot_data <- data %>% 
  group_by(group, Hopx) %>% 
  summarise(value = n()) %>% 
  group_modify(~.x %>% mutate(prec = value / sum(value) * 100)) %>% 
  ungroup %>% 
  filter(Hopx == "pos")

ggbarplot(plot_data, "group", "prec", fill = "group", xlab = "", ylab = "Cell precentage (%)",
          legend = "right") +
  geom_signif(annotations = "0.0692", y_position = 45, xmin = 1, xmax = 2, tip_length = 0.15, size = .5) +
  scale_fill_brewer(palette = "Set2")

#### scRNA BC vs PBS 上调基因富集分析 ####
# pathway
data <- read.csv("CBCs BC vs PBS eKEGG.csv") %>% 
  dplyr::select(all_of(c("Description","GeneRatio","pvalue"))) %>% 
  mutate(GeneRatio = as.numeric(stringr::str_split_i(GeneRatio, "/", 1))/
           as.numeric(stringr::str_split_i(GeneRatio, "/", 2))) %>% 
  arrange(desc(GeneRatio)) %>% 
  mutate(Description = factor(Description, (Description)))

ggbarplot(data, "Description", "GeneRatio", fill = "pvalue", rotate = T, width = .6, 
          ylab = "Gene ratio", xlab = "", legend = "right",
          title = "KEGG pathways enrichment for up-regulated\ngenes in stem cells after BC-colonization") +
  scale_fill_gradient(low = "#abd9e9", high = "#fdae61")

# module
diff <- read.csv("CBCs BC vs PBS deg.csv") %>% 
  filter(change == "up")

bg_data <- read.delim("mmu_gene2mod_enrichment.tsv")

eKEGG <- enricher(gene = diff$ENTREZID, TERM2GENE = bg_data, minGSSize = 1, pvalueCutoff = 1, qvalueCutoff = 1)
eKEGG_out <- data.frame(eKEGG@result)

data <- read.delim("scRNA.CBCs.up_gene.ekegg.module.for_plot.tsv") %>% 
  mutate(GeneRatio = as.numeric(stringr::str_split_i(GeneRatio, "/", 1))/
           as.numeric(stringr::str_split_i(GeneRatio, "/", 2))) %>% 
  arrange((GeneRatio))

ggbarplot(data, "Description", "GeneRatio", fill = "pvalue", rotate = T, width = .6, 
          ylab = "Gene ratio", xlab = "", legend = "right",
          title = "KEGG pathways enrichment for up-regulated\ngenes in stem cells after BC-colonization") +
  scale_fill_gradient(low = "#abd9e9", high = "#fdae61")

#### CBCs ppara ####
# source from CBCs BC vs PBS deg.csv
# pct.1 0.087, pct.2 0.092 pvalue 0.7022
plot_data <- data.frame(group = c("BC", "PBS"), value = c(0.087, 0.092))

ggdotchart(plot_data, "group", "value", color = "group", add = "segments", legend = "none",
           palette = c("#55d5d5", "#a0a0a4"), dot.size = 6, ylab = "Expression of Ppara",
           size = 1) +
  geom_signif(annotations = "0.7022", xmin = 1, xmax = 2, y_position = 0.1, tip_length = 0.02)
# ggsave("scRNA.CBCs.ppara_exp.pdf", width = 3, height = 3)

#### lgr5 and hmgcs2 共表达 ####
seurat_obj <- RunTSNE(seurat_obj, dims = 1:5, reduction = "harmony")

# 分群
DimPlot(seurat_obj, reduction = "tsne") +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = c(1, 1),
        legend.justification = c(1, 1))

allmarkers <- FindAllMarkers(seurat_obj, only.pos = TRUE)

top_markers <- allmarkers %>%
  group_by(cluster) %>%
  top_n(n = 8, wt = avg_log2FC) %>% 
  pull
top_markers[grepl("Jph1", top_markers)] <- "Lgr5"

DotPlot(seurat_obj, features = unique(top_markers), cols = c("#abd9e9", "#fdae61"), dot.scale = 3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 8, hjust = 1)) +
  NoLegend()

# 共表达
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
