#### Jinxin Meng, 20240801, 20240806 ####

setwd("F:/proj/20240731_BC_IBD_Zhangyn/submit/code_availability/RNAseq/")
pacman::p_load(tidyr, dplyr, tibble, purrr, ggpubr)
library(DESeq2) %>% suppressMessages()
library(ComplexHeatmap) %>% suppressMessages()
library(clusterProfiler) %>% suppressMessages()
source("/Code/R_func/profile_process.R")

#### day 7 ####
# gene ab
fpkm <- read.delim("fpkm_profile.txt", row.names = 1, check.names = F) %>% 
  select(contains("7d"))

group <- read.delim("mapping.txt") %>% 
  filter(sample %in% colnames(fpkm))

group_levels <- c("Veh_7d","BC_7d")

markers <- list(
  aISC = c("Lgr5","Sox9","Ascl2"),
  rISC = c("Hopx","Clu","Bmi1"),
  PC = c("Ang4","Lyz1"),
  GC = c("Muc2","Ccl9"),
  TC = c("Dclk1","Cd24a","Trpm5"),
  EEC = c("Chga", "Neurog3")
)

marker_levels <- c("aISC", "rISC", "PC", "GC", "TC", "EEC")

fpkm <- filter(fpkm, rownames(fpkm) %in% unlist(markers, use.names = F)) 

annotation_row <- markers %>% 
  map2_dfr(., names(.), \(x, y) data.frame(name = x, class = y)) %>% 
  mutate(name = factor(name, rownames(fpkm))) %>% 
  arrange(name) %>% 
  column_to_rownames("name") %>% 
  mutate(class = factor(class, marker_levels))

annotation_col <- group %>% 
  select(sample, group) %>% 
  mutate(sample = factor(sample, colnames(fpkm))) %>% 
  arrange(sample) %>% 
  column_to_rownames("sample") %>% 
  mutate(group = factor(group, group_levels))

split_row <- annotation_row$class
split_col <- annotation_col$group

pheatmap(fpkm, 
         scale = "row",
         color = colorRampPalette(c("#307cc0", "white", "#e43589"))(100),
         split = split_row, column_split = split_col,
         cluster_row_slices = F, cluster_column_slices = F, 
         row_title_rot = 0,
         row_title_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8),
         row_gap = unit(1, "mm"), column_gap = unit(1, "mm"),
         treeheight_col = 15, treeheight_row = 15,
         cellheight = 12, cellwidth = 12,
         border_color = "white",
         border_gp = gpar(col = "black")
)

# gene diff
rc <- read.delim("rc_profile.txt", row.names = 1, check.names = F) %>% 
  select(all_of(group$sample))
rc <- rc[rowSums(rc > 5) >= 2, ]

group_levels <- c("BC_7d", "Veh_7d")

metadata <- group %>% 
  select(sample, group) %>% 
  mutate(sample = factor(sample, colnames(rc))) %>% 
  arrange(sample) %>% 
  column_to_rownames(var = "sample") %>% 
  mutate(group = factor(group, group_levels))

dds <- DESeqDataSetFromMatrix(countData = rc, colData = metadata, design = ~ group)
DEseq_obj <- DESeq(dds)

diff <- results(DEseq_obj, contrast = c("group", group_levels)) %>% 
  data.frame(.) %>% 
  rownames_to_column(var = "name") %>% 
  mutate(enriched = ifelse(log2FoldChange > 1 & pvalue < 0.05, "up", 
                           ifelse(log2FoldChange < -1 & pvalue < 0.05, "down", "none"))) %>% 
  filter(enriched == "up")

gene_info <- read.delim("gene_info.txt")

genes <- gene_info %>%
  filter(Gene%in%diff$name & EntrezGene!="-") %>% 
  dplyr::select(1:2)

# kegg
eKEGG <- enrichKEGG(gene = na.omit(genes$EntrezGene), organism = "mmu", pvalueCutoff = 1, qvalueCutoff = 1)
eKEGG_out <- data.frame(eKEGG) %>% 
  mutate(geneName = lapply((strsplit(x = geneID, "/")), 
                           function(x) genes$Gene[match(x, genes$EntrezGene)]) %>% 
           sapply(., function(x) paste(x, collapse = "/")), .after = geneID, 
         name = stringr::str_replace(Description, " - Mus.*$", ""))

# GO 
library(org.Mm.eg.db)
eGO <- enrichGO(gene = na.omit(genes$EntrezGene), org.Mm.eg.db, pvalueCutoff = 1, qvalueCutoff = 1)
eGO_out <- data.frame(eGO) %>%
  mutate(geneName = lapply((strsplit(x = geneID, "/")),
                           function(x) genes$Gene[match(x, genes$EntrezGene)]) %>%
           sapply(., function(x) paste(x, collapse = "/")), .after = geneID)

# custom
custom_path <- read.delim("custom_path_library.tsv")

enricher <- enricher(genes$Gene, TERM2GENE = custom_path)
enricher_out <- data.frame(enricher@result) 


plot_data <- rbind(eKEGG_out %>% 
        filter(ID %in% c("mmu04151","mmu04310","mmu04530","mmu04010","mmu04152","mmu00512")) %>% 
        mutate(Description = gsub(" -.*", "", Description)) %>% 
        dplyr::select(Description, GeneRatio, pvalue) %>%
        data.frame(row.names = NULL),
      eGO_out %>% 
        filter(ID %in% c("GO:0005200","GO:0008083","GO:0050431")) %>% 
        dplyr::select(Description, GeneRatio, pvalue) %>% 
        data.frame(row.names = NULL), 
      enricher_out %>% 
        filter(ID == "regulation_of_BMP_signaling_pathway") %>% 
        dplyr::select(Description, GeneRatio, pvalue) %>% 
        data.frame(row.names = NULL)) %>% 
  mutate(GeneRatio = as.numeric(stringr::str_split_i(GeneRatio, "/", 1))/
           as.numeric(stringr::str_split_i(GeneRatio, "/", 2))) %>% 
  arrange(desc(GeneRatio)) %>% 
  filter(pvalue < 0.05)

ggscatter(plot_data, "Description", "GeneRatio", fill = "pvalue", rotate = T, shape = 21, size = 4, 
          x.text.angle = 90, legend = "right") +
  scale_fill_viridis_c(direction = -1, begin = .3, end = 1)

#### day 3 ####
fpkm <- read.delim("fpkm_profile.txt", row.names = 1, check.names = F) %>% 
  dplyr::select(contains("3d"))

group <- read.delim("mapping.txt") %>% 
  filter(sample %in% colnames(fpkm))

group_levels <- c("Veh_3d","BC_3d")

markers <- list(
  aISC = c("Lgr5","Sox9","Ascl2"),
  rISC = c("Hopx","Clu","Bmi1"),
  PC = c("Ang4","Lyz1"),
  GC = c("Muc2","Ccl9"),
  TC = c("Dclk1","Cd24a","Trpm5"),
  EEC = c("Chga", "Neurog3")
)

marker_levels <- c("aISC", "rISC", "PC", "GC", "TC", "EEC")

fpkm <- filter(fpkm, rownames(fpkm) %in% unlist(markers, use.names = F)) 

annotation_row <- markers %>% 
  map2_dfr(., names(.), \(x, y) data.frame(name = x, class = y)) %>% 
  mutate(name = factor(name, rownames(fpkm))) %>% 
  arrange(name) %>% 
  column_to_rownames("name") %>% 
  mutate(class = factor(class, marker_levels))

annotation_col <- group %>% 
  dplyr::select(sample, group) %>% 
  mutate(sample = factor(sample, colnames(fpkm))) %>% 
  arrange(sample) %>% 
  column_to_rownames("sample") %>% 
  mutate(group = factor(group, group_levels))

split_row <- annotation_row$class
split_col <- annotation_col$group

pheatmap(fpkm, 
         scale = "row",
         color = colorRampPalette(c("#307cc0", "white", "#e43589"))(100),
         split = split_row, column_split = split_col,
         cluster_row_slices = F, cluster_column_slices = F, 
         row_title_rot = 0,
         row_title_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8),
         row_gap = unit(1, "mm"), column_gap = unit(1, "mm"),
         treeheight_col = 15, treeheight_row = 15,
         cellheight = 12, cellwidth = 12,
         border_color = "white",
         border_gp = gpar(col = "black")
)

# gene diff
rc <- read.delim("rc_profile.txt", row.names = 1, check.names = F) %>% 
  select(all_of(group$sample))
rc <- rc[rowSums(rc > 5) >= 2, ]

group_levels <- c("BC_3d", "Veh_3d")

metadata <- group %>% 
  select(sample, group) %>% 
  mutate(sample = factor(sample, colnames(rc))) %>% 
  arrange(sample) %>% 
  column_to_rownames(var = "sample") %>% 
  mutate(group = factor(group, group_levels))

dds <- DESeqDataSetFromMatrix(countData = rc, colData = metadata, design = ~ group)
DEseq_obj <- DESeq(dds)

diff <- results(DEseq_obj, contrast = c("group", group_levels)) %>% 
  data.frame(.) %>% 
  rownames_to_column(var = "name") %>% 
  mutate(enriched = ifelse(log2FoldChange > 1 & pvalue < 0.05, "up", 
                           ifelse(log2FoldChange < -1 & pvalue < 0.05, "down", "none")))
write.table(diff, "RNAseq.D3.diff.tsv", sep = "\t", quote = F, row.names = F)

#### day 0 ####
fpkm <- read.delim("fpkm_profile.txt", row.names = 1, check.names = F) %>% 
  dplyr::select(contains("0d"))

group <- read.delim("mapping.txt") %>% 
  filter(sample %in% colnames(fpkm))

group_levels <- c("Veh_0d","BC_0d")

markers <- list(
  aISC = c("Lgr5","Sox9","Ascl2"),
  rISC = c("Hopx","Bmi1","Clu"),
  PC = c("Ang4","Lyz1"),
  GC = c("Muc2","Ccl9"),
  TC = c("Dclk1","Cd24a","Trpm5"),
  EEC = c("Chga", "Neurog3")
)

marker_levels <- c("aISC", "rISC", "PC", "GC", "TC", "EEC")

fpkm <- filter(fpkm, rownames(fpkm) %in% unlist(markers, use.names = F)) 

annotation_row <- markers %>% 
  map2_dfr(., names(.), \(x, y) data.frame(name = x, class = y)) %>% 
  mutate(name = factor(name, rownames(fpkm))) %>% 
  arrange(name) %>% 
  column_to_rownames("name") %>% 
  mutate(class = factor(class, marker_levels))

annotation_col <- group %>% 
  dplyr::select(sample, group) %>% 
  mutate(sample = factor(sample, colnames(fpkm))) %>% 
  arrange(sample) %>% 
  column_to_rownames("sample") %>% 
  mutate(group = factor(group, group_levels))

split_row <- annotation_row$class
split_col <- annotation_col$group

pheatmap(fpkm, 
         scale = "row",
         color = colorRampPalette(c("#307cc0", "white", "#e43589"))(100),
         split = split_row, column_split = split_col,
         cluster_row_slices = F, cluster_column_slices = F, 
         row_title_rot = 0,
         row_title_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8),
         row_gap = unit(1, "mm"), column_gap = unit(1, "mm"),
         treeheight_col = 15, treeheight_row = 15,
         cellheight = 12, cellwidth = 12,
         border_color = "white",
         border_gp = gpar(col = "black")
)
