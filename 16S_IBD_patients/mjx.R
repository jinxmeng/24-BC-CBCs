#### Jinxin Meng, 20240801, 20241112 ####

setwd("F:/proj/20240731_BC_IBD_Zhangyn/submit/code_availability/Figure 1/")
pacman::p_load(tidyr, dplyr, tibble, purrr, ggpubr)
source("/code/R_func/difference_analysis.R")
source("/code/R_func/calcu_diff.R")
source("/code/R_func/taxa.R")

#### violin plot ####

profile <- read.delim("profile_genus.txt", row.names = 1, check.names = F) %>% 
  profile_transRA()
group <- read.delim("sample_group.txt")

plot_data <- profile %>% 
  filter(rownames(.) %in% "g__Blautia") %>% 
  t() %>% 
  data.frame() %>% 
  dplyr::rename(value = 1) %>% 
  rownames_to_column("sample") %>% 
  left_join(group, by = "sample") %>%
  mutate(group = factor(group, c("Control", "IBD")))

ggviolin(plot_data, "group", "value", color = "group", legend = "none", xlab = "", size = 1, width = .6,
         ylab = "Relative abundance(%)", palette = c("#69a7bc","#E69F00"),
         outlier.shape = NA) +
  geom_boxplot(aes(color = group), width = .1, size = 1, outlier.shape = NA) + 
  lims(y = c(NA, 30)) +
  geom_jitter(aes(color = group), width = .2) +
  stat_compare_means(comparisons = list(c("IBD", "Control")), label = "p.signif", label.y = 27,
                     tip.length = .01, step.increase = .02, vjust = 1, size = 3) +
  theme(aspect.ratio = 1.5)
# ggsave("g__Blautia.rela_ab.pdf", width = 4, height = 4)

#### rf + AUC ####
source("/code/R_func/model_randomforest.R")
source("/code/R_func/plot_roc.R")

profile <- read.delim("profile_genus.txt", row.names = 1, check.names = F) %>% 
  profile_transRA()
group <- read.delim("sample_group.txt")

pred <- rf_loom(profile, group, seed = 2024, ntree = 1000)
pred <- left_join(pred, group, by = "sample")
roc <- roc(pred$group, pred$IBD)
plot_roc(roc, plot_se = T)
# ggsave("roc.plot.pdf", width = 4, height = 4)

diff <- difference_analysis(profile, group) %>% 
  mutate(enriched = ifelse(Control_ab > IBD_ab, "Control", "IBD"))

data <- rf_variable_rank(profile, group, seed = 2024, ntree = 1000)

left_join(data, diff, by = c("variable" = "name")) %>% 
  head(10) %>% 
  ggbarplot("variable", "MeanDecreaseAccuracy", fill = "enriched", rotate = T, sort.val = "asc",
            sort.by.groups = F, width = .6, palette = c("#69a7bc","#E69F00"), xlab = "") +
  theme(aspect.ratio = 1.4)
# ggsave("variable_rank.pdf", width = 5.5, height = 5)

#### lefse ####
# microeco是基于R6class开发的
library(ggtree)
library(microeco)
library(magrittr)
set.seed(2024)

profile <- read.delim("table_rarefied.tsv", row.names = 1, check.names = F) %>%
  profile_transRA()
group <- read.delim("sample_group.txt", row.names = 1)
tax <- read.delim("taxonomy.tsv", row.names = 1)
tr <- read.tree("tree.nwk")

#通常，分类表中有一些杂乱的信息，例如NA，未识别和未知需要清理
tax %<>% tidy_taxonomy

#构建microeco dataset
dataset <- microtable$new(sample_table = group, otu_table = profile, tax_table = tax, phylo_tree = tr)

#lefse分析
lefse <- trans_diff$new(dataset = dataset, method = "lefse", group = "group", alpha = 0.05, lefse_subgroup = NULL)

#lefse柱状图
lefse$plot_diff_bar(use_number = 1:30, width = 0.6, group_order = c("Control", "IBD"), 
                    color_values = c("#69a7bc","#E69F00"))
# ggsave("lefse.bar.pdf", width = 8, height = 10)

#lefse进化树
#此数据集中的分类单元太多。例如，我们只选择树中丰度最高的200个分类群和50个差异特征。
lefse$plot_diff_cladogram(use_taxa_num = 100, use_feature_num = 30, clade_label_level = 4, 
                          group_order = c("Control", "IBD"), color = c("#69a7bc","#E69F00"))
# ggsave("lefse.cladogram.pdf", width = 10, height = 8)

