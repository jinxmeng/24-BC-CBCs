#### Jinxin Meng, 20240807, 20240807 ####

setwd("F:/proj/20240731_BC_IBD_Zhangyn/submit/code_availability/Figure 1/public_MGX/")
pacman::p_load(tidyr, dplyr, tibble, purrr, ggpubr)
source("/code/R_func/calcu_diff.R")

#### SchirmerM_2018.PRJNA389280 ####

group <- read.delim("SchirmerM_2018.PRJNA389280.sample_group") %>% 
  select(sample = run, group = group2)
data <- read.delim("SchirmerM_2018.PRJNA389280.profile") %>% 
  mutate(rela_ab = rc / libs * 100) %>% 
  left_join(group, by = "sample") %>% 
  mutate(group = factor(group, c("CT", "CD", "UC")))

comparisons <- calcu_diff(select(data, sample, value = rela_ab), group) %>% 
  filter(pval < 0.05) %>% 
  pull(gp) %>% 
  strsplit(split = "_vs_")

ggviolin(data, "group", "rela_ab", color = "group", legend = "none", xlab = "", size = 1, width = .5,
         ylab = "Relative abundance(%)",  palette = c("#69a7bc","#E69F00","#f08178"),
         outlier.shape = NA, title = "SchirmerM_2018") +
  geom_boxplot(aes(color = group), width = .2, size = 1, outlier.shape = NA) + 
  geom_jitter(aes(color = group), width = .2) +
  lims(y = c(NA, .85)) +
  stat_compare_means(comparisons = comparisons, label = "p.signif", label.y = .75,
                     tip.length = .01, step.increase = .042, vjust = .9, size = 3)
# ggsave("SchirmerM_2018.pdf", width = 3.5, height = 4)

#### WengY_2019.PRJNA429990 ####

group <- read.delim("WengY_2019.PRJNA429990.sample_group")
data <- read.delim("WengY_2019.PRJNA429990.profile") %>% 
  mutate(rela_ab = rc / libs * 100) %>% 
  left_join(group, by = "sample") %>% 
  mutate(group = factor(group, c("Control", "CD", "UC")))

ggviolin(data, "group", "rela_ab", color = "group", legend = "none", xlab = "", size = 1, width = .5,
         ylab = "Relative abundance(%)",  palette = c("#69a7bc","#E69F00","#f08178"),
         outlier.shape = NA, title = "WengY_2019") +
  geom_boxplot(aes(color = group), width = .2, size = 1, outlier.shape = NA) + 
  geom_jitter(aes(color = group), width = .2) +
  lims(y = c(NA, 1)) +
  stat_compare_means(comparisons = list(c("CD", "Control"), c("UC", "Control"), c("CD", "UC")),
                     label.y = .38, tip.length = .003, step.increase = .01, vjust = .9, size = 3) +
  lims(y = c(NA, .5))
# ggsave("WengY_2019.pdf", width = 3.5, height = 4)

#### HeQ_2017.PRJEB15371 ####
group <- read.delim("HeQ_2017.PRJEB15371.sample_group")
data <- read.delim("HeQ_2017.PRJEB15371.profile") %>% 
  mutate(rela_ab = rc / libs * 100) %>% 
  left_join(group, by = "sample") %>% 
  mutate(group = factor(group, c("Control", "CD")))

comparisons <- calcu_diff(select(data, sample, value = rela_ab), group) %>% 
  filter(pval < 0.05) %>% 
  pull(gp) %>% 
  strsplit(split = "_vs_")

ggviolin(data, "group", "rela_ab", color = "group", legend = "none", xlab = "", size = 1, width = .5,
         ylab = "Relative abundance(%)",  palette = c("#69a7bc","#E69F00","#f08178"),
          outlier.shape = NA, title = "HeQ_2017") +
  geom_boxplot(aes(color = group), width = .2, size = 1, outlier.shape = NA) + 
  geom_jitter(aes(color = group), width = .2) +
  lims(y = c(NA, .8)) +
  stat_compare_means(comparisons = comparisons, label = "p.signif", label.y = .65,
                     tip.length = .004, step.increase = .042, vjust = .9, size = 3)
# ggsave("HeQ_2017.boxplot.pdf", width = 2.5, height = 4)

#### YanQ_2023c.PRJEB67456 ####
group <- read.delim("YanQ_2023c.PRJEB67456.sample_group") %>% 
  select(sample = run, group)

data <- read.delim("YanQ_2023c.PRJEB67456.profile") %>% 
  mutate(rela_ab = rc / libs * 100) %>% 
  left_join(group, by = "sample") %>% 
  mutate(group = factor(group, c("healthy", "CD", "UC")))

ggviolin(data, "group", "rela_ab", color = "group", legend = "none", xlab = "", size = 1, width = .5,
         ylab = "Relative abundance(%)",  palette = c("#69a7bc","#E69F00","#f08178"),
         outlier.shape = NA, title = "YanQ_2023") +
  geom_boxplot(aes(color = group), width = .2, size = 1, outlier.shape = NA) + 
  geom_jitter(aes(color = group), width = .2) +
  lims(y = c(NA, 1.6)) +
  stat_compare_means(comparisons = list(c("healthy", "UC"), c("healthy", "CD")), label.y = 1.4,
                     tip.length = .01, step.increase = .042, vjust = .9, size = 3)
# ggsave("YanQ_2023c.pdf", width = 3.5, height = 4)

