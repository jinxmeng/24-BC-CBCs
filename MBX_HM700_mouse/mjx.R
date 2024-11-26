#### Jinxin Meng, 20240801, 20241113 ####

setwd("F:/proj/20240731_BC_IBD_Zhangyn/submit/code_availability/MBX_HM700_mouse/")
pacman::p_load(tidyr, dplyr, tibble, purrr, ggpubr)
library(ropls)
library(clusterProfiler)
source("/code/R_func/difference_analysis.R")

#### 差异分析 ####

metabolite_info <- read.delim("metabolite_info.txt")
profile <- read.delim("profile.txt", row.names = 1)
group <- read.delim("group.txt")

gps <- list(
  c("BC_0d", "Veh_0d"),
  c("BC_3d", "Veh_3d"))

gp <- gps[[2]]

samples <- group %>% 
  filter(group %in% gp) %>% 
  pull(sample)
profile_x <- profile %>% 
  select(all_of(samples))
group_x <- group %>% 
  filter(sample %in% samples)

# 监督分组以疾病-健康为例，orthoI = NA 时执行 OPLS-DA, 自动计算正交分量的数量；默认设置为0，执行的是PLS分析
oplsda <- opls(x = data.frame(t(profile_x), check.names = F), pull(group_x, group), orthoI = NA, predI = 1)

diff <- difference_analysis(profile_x, group_x, gp = gp, diff_method = "wilcox")

out <- oplsda@vipVn %>%
  data.frame(vip = .) %>% 
  rownames_to_column(var = "cpd_id") %>% 
  left_join(diff, ., by = c("name" = "cpd_id")) %>% 
  mutate(enriched = ifelse(vip > 1 & log2FC > 1, gp[1], 
                           ifelse(vip > 1 & log2FC < -1, gp[2], "none")))

table(out$enriched)

plot_data <- out %>% 
  select(name, log2FC, vip, enriched) %>% 
  mutate(enriched = factor(enriched, c("BC_3d", "none", "Veh_3d")),
         log2FC = ifelse(log2FC > 3, 3, log2FC),
         log2FC = ifelse(log2FC < -3, -3, log2FC)) %>% 
  left_join(metabolite_info %>% select(1:2, cpd_class), by = c("name" = "cpd_id")) %>% 
  mutate(label = ifelse(enriched != "none", cpd_name, ""))

ggscatter(plot_data, "log2FC", "vip", color = "enriched", legend = "right",
          palette = c("#eeacec", "grey", "#21c1dc"), size = 2) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_text(aes(label = label), size = 1, fontface = "italic")

#### 富集分析 通路 ####
bg_data <- read.delim("cpd2path_enrichment.tsv")

cpds <- out %>% 
  filter(enriched == gp[1]) %>% 
  select(cpd_id = name) %>% 
  left_join(metabolite_info, by = "cpd_id") %>% 
  filter(cpd_KEGG != "") %>% 
  pull(cpd_KEGG)

eKEGG <- enricher(gene = cpds, TERM2GENE = bg_data, minGSSize = 1, pvalueCutoff = 1, qvalueCutoff = 1)
eKEGG_out <- data.frame(eKEGG@result)

plot_data <- read.delim("eKEGG.BC_3d_vs_Veh_3d.up.for_plot.tsv") %>% 
  dplyr::select(all_of(c("Description","GeneRatio","pvalue"))) %>% 
  mutate(GeneRatio = as.numeric(stringr::str_split_i(GeneRatio, "/", 1))/
           as.numeric(stringr::str_split_i(GeneRatio, "/", 2))) %>% 
  arrange(desc(GeneRatio)) %>% 
  filter(pvalue < 0.05) %>% 
  mutate(Description = factor(Description, rev(Description)))

ggscatter(plot_data, "GeneRatio", "Description", fill = "pvalue", shape = 21, size = 4, 
          x.text.angle = 90, legend = "right") +
  scale_fill_viridis_c(direction = -1, begin = .3, end = 1)

plot_data <- read.delim("eKEGG.bc_3d_vs_pbs_3d.up.for_plot.tsv") %>% 
  dplyr::select(all_of(c("Description","GeneRatio","pvalue"))) %>% 
  mutate(GeneRatio = as.numeric(stringr::str_split_i(GeneRatio, "/", 1))/
           as.numeric(stringr::str_split_i(GeneRatio, "/", 2))) %>% 
  arrange(desc(GeneRatio)) %>% 
  filter(pvalue < 0.05) %>% 
  mutate(Description = factor(Description, (Description)))

#### 富集分析 模块 ####
bg_data <- read.delim("cpd2mod_enrichment.tsv")

cpds <- out %>% 
  filter(enriched == gp[1]) %>% 
  select(cpd_id = name) %>% 
  left_join(metabolite_info, by = "cpd_id") %>% 
  filter(cpd_KEGG != "") %>% 
  pull(cpd_KEGG)

eKEGG <- enricher(gene = cpds, TERM2GENE = bg_data, minGSSize = 1, pvalueCutoff = 1, qvalueCutoff = 1)
eKEGG_out <- data.frame(eKEGG@result)

plot_data <- read.delim("eKEGG.BC_3d_vs_Veh_3d.up.mod.for_plot.tsv") %>% 
  dplyr::select(all_of(c("Description","GeneRatio","pvalue"))) %>% 
  mutate(GeneRatio = as.numeric(stringr::str_split_i(GeneRatio, "/", 1))/
           as.numeric(stringr::str_split_i(GeneRatio, "/", 2))) %>% 
  arrange(desc(GeneRatio)) %>% 
  filter(pvalue < 0.05) %>% 
  mutate(Description = factor(Description, rev(Description)))

ggbarplot(plot_data, "Description", "pvalue", fill = "pvalue", rotate = T, width = .6) +
  scale_fill_gradient(low = "#abd9e9", high = "#fdae61")

#### 相关性分析 ####
# PPAR genes & fatty acids
diff <- read.delim("../RNAseq/RNAseq.D3.diff.tsv") %>% 
  filter(pvalue < 0.05 & abs(log2FoldChange) > 1)

PPAR_genes <- read.delim("PPAR_path_genes.txt")
fpkm_genes <- read.delim("../RNAseq/fpkm_profile.txt", row.names = 1) %>% 
  filter(rownames(.) %in% PPAR_genes$symbol) %>% 
  select(contains("3d")) %>% 
  filter(colSums(apply(., 1, \(x) x > 10)) == 6)

fatty_acids <- read.delim("metabolite_info.txt") %>% 
  filter(grepl("Fatty acid", cpd_class)) %>% 
  pull(1)

ab_metabolites <- read.delim("profile.txt", row.names = 1) %>% 
  filter(rownames(.) %in% fatty_acids) %>% 
  select(contains("3d"), -BC_3d_4, -Veh_3d_4)

source("corr_process.R")
cor_obj <- psych::corr.test(t(fpkm_genes), t(ab_metabolites), method = "spearman")
cor_out <- corr_merge(cor_obj$r, cor_obj$p, rho = .5, pval = .05)

metabolite_info <- read.delim("metabolite_info.txt")

edges <- cor_out %>% 
  rownames_to_column("gene") %>% 
  gather(key = "cpd_id", value = "rval", -gene) %>% 
  filter(rval != 0) %>% 
  mutate(metabolite = metabolite_metabolite_info$cpd_name[match(cpd_id, metabolite_metabolite_info$cpd_id)]) %>% 
  select(from = gene, to = metabolite, rval)
write.table(edges, "nwk.edges.tsv", sep = "\t", quote = F, row.names = F)

edges <- read.delim("nwk.edges.tsv")
nodes <- read.delim("nwk.nodes.tsv")

pacman::p_load(igraph, ggraph, tidygraph)

graph <- tbl_graph(nodes = nodes, edges = edges)

ggraph(graph, layout = "linear", circular = T) + 
  geom_edge_arc(aes(color = factor(dir, levels = c("pos", "neg"))), width = .4, ) +
  scale_edge_color_manual(values = c("#f1a340", "#998ec3")) +
  geom_node_point(aes(shape = type, fill = factor(enriched, c("bc", "none", "pbs"))), 
                  color = "#000000", size = 1.7, show.legend = F, stroke = .2) +
  geom_node_text(aes(x = x*1.03, y = y*1.03, label = node, angle = atan(y/x)*360/(2*pi), 
                     hjust = "outward", color = type), size = 1.2) +
  scale_color_manual(values = c("#f1a340", "#998ec3")) +
  scale_fill_manual(values = c("#f1a340", "grey", "#998ec3")) +
  scale_shape_manual(values = c(21, 24)) +
  labs(color = "enriched", edge_color = "enriched") +
  coord_fixed() +
  theme_graph()

#### BHB MS 可视化 ####
source("ms_utilites.R")
library(xcms)
library(MsExperiment)

# 公司给的数据，应该数据SRM（选择性反应监测，Selected Reaction Monitoring）的数据
# MSnbase包 MChromatograms对象
# https://bioconductor.org/packages/release/bioc/vignettes/MSnbase/inst/doc/v01-MSnbase-demo.html
# data <- readSRMData("../data/submit/submit/BC_0d_1.mzML") # 读入色谱数据 
# polarity(data) # 计算离子的极性
# precursorMz(data) # 前体分子的质荷比
# productMz(data) # 产物分子的质荷比
# plot(data[4]) # 画图 rt和强度
# rtime(data[4]) # 查看某一个光谱的停留时间
# intensity(data[4]) # 查看某一个光谱的强度

# Indole-3-propionic acid （HMDBHMDB0002302）的RT=7.16 min， 
# 2-Hydroxybutyric acid （HMDB0000008）的RT=3.8 min，
# 3-Hydroxybutyric acid （HMDB0000011）的RT=3.48 min； 
rt = 3.48

files <- list.files("../submit/submit/", pattern = "3d", full.names = T)
files <- files[grepl("mzML", files)]
names <- list.files("../submit/submit/", pattern = "3d")
names <- names[grepl("mzML", names)]
names <- sub(".mzML", "", names)

group_col <- structure(c(rep('#f8766d', 4), rep('#00bfc4', 4)),
                       names = c("BC_3d_1","BC_3d_2","BC_3d_3","BC_3d_4",
                                 "Veh_3d_1","Veh_3d_2","Veh_3d_3","Veh_3d_4"))

ms_data <-  map(files, \(x) readSRMData(x)) %>% setNames(names)

ms_parser <- map(ms_data, \(x) {
  parser_x <- parse_MChromatograms(x)
  filter(parser_x, rt_min < rt & rt_max > rt) } )

plot_list <- map(ms_parser$BC_3d_1$name, \(x) {
  plot_data <- map2_df(ms_data, names(ms_data), \(m, n) 
                       data.frame(rtime = rtime(m[x]), 
                                  intensity = intensity(m[x])) %>% 
                         add_column(class = n) ) %>% 
    mutate(rtime = as.numeric(rtime))
  
  ggplot(plot_data, aes(rtime, intensity, group = class)) +
    geom_line(aes(color = class), show.legend = F) +
    geom_vline(xintercept = rt, linetype = "dashed") +
    scale_color_manual(values = group_col) +
    labs(title = paste0("ms: ", x)) +
    theme_pubr() +
    theme(aspect.ratio = .7) } )
cowplot::plot_grid(plotlist = plot_list, align = "v", nrow = 10)

# MS 43 为 BHB
plot_data <- map2_df(ms_data, names(ms_data), \(m, n) 
                     data.frame(rtime = rtime(m["43"]), 
                                intensity = intensity(m["43"])) %>% 
                       add_column(class = n) ) %>% 
  mutate(rtime = as.numeric(rtime))

options(scipen = 10)

plot_data %>% 
  filter(rtime < 4.2) %>%
  filter(rtime > 3.53 | rtime < 3.49) %>% 
  filter(!intensity %in% c(263547, 228244, 203055, 234293)) %>%
  ggplot(aes(rtime, intensity, group = class)) +
  ggalt::geom_xspline(aes(color = class), size = 1, spline_shape = 1.1) +
  scale_color_manual(values = group_col) +
  labs(x = 'Retention time (min)', y = 'Intensity') +
  theme_pubr() +
  theme(aspect.ratio = .9) +
  guides(color = guide_legend(position = "right"))
MS kes# ggsave("BHB.pdf", width = 6, height = 4)
