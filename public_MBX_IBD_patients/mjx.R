#### Jinxin Meng, 20240801, 20240814 ####

setwd("F:/proj/20240731_BC_IBD_Zhangyn/submit/code_availability/public_MBX_IBD_patients/")
pacman::p_load(tidyr, dplyr, tibble, purrr, ggpubr, ggpmisc)
library(ropls)
source("/code/R_func/difference_analysis.R")
source("/code/R_func/calcu_diff.R")

#### difference analysis ####

metadata <- read.delim("metadata.txt")
profile <- read.delim("profile.txt", row.names = 1, check.names = F) %>% 
  profile_transRA()
group <- read.delim("sample_group.txt")

plot_data <- profile_smp2grp_2(profile, group, method = "mean") %>% 
  rowwise() %>% 
  mutate(avg_ab = (CD + UC + nonIBD) / 3) %>% 
  left_join(metadata, by = c("name" = "cpd_id")) %>% 
  mutate(label = ifelse(cpd_name == "phenyllactate", cpd_name, ""))

library(ggtern)
ggtern(plot_data, aes(CD, UC, nonIBD)) + 
  geom_point(aes(size = avg_ab), color = "#fc8d62") +
  geom_text(aes(label = label), size = 1) +
  scale_size_continuous(range = c(1, 6)) +
  theme_bw() +  
  theme_bvbw()

plot_data <- data.frame(t(profile["cpd_242",])) %>% 
  rename(value = 1) %>% 
  rownames_to_column("sample") %>% 
  left_join(group, by = "sample") %>% 
  mutate(group = factor(group, c("nonIBD", "CD", "UC")))

comparisons <- calcu_diff(dplyr::select(plot_data, sample, value), group) %>% 
  filter(pval < 0.05) %>% 
  pull(gp) %>% 
  strsplit(split = "_vs_")

ggviolin(plot_data, "group", "value", color = "group", legend = "none", xlab = "", size = 1, width = .5,
         ylab = "Relative abundance(%)", palette = c("#69a7bc","#E69F00","#f08178"),
         outlier.shape = NA, title = "LloydPriceJ_2019") +
  geom_boxplot(aes(color = group), width = .2, size = 1, outlier.shape = NA) + 
  geom_jitter(aes(color = group), width = .2) +
  lims(y = c(NA, .01)) +
  stat_compare_means(comparisons = comparisons, label = "p.signif", label.y = .008,
                     tip.length = .004, step.increase = .042, vjust = .9, size = 3)
