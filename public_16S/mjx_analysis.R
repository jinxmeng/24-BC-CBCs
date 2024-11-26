#### Jinxin Meng, 20240801, 20240815 ####

setwd("F:/proj/20240731_BC_IBD_Zhangyn/submit/code_availability/Figure 1/public_16S/")
pacman::p_load(tidyr, dplyr, tibble, purrr, ggpubr)
source("/code/R_func/profile_process.R")

#### code ####
group <- read.delim("sample_group.txt") %>% 
  dplyr::select(sample = external, group) %>% 
  mutate(sample = as.character(sample))

profile <- read.delim("profile.txt", row.names = 1, check.names = F)

taxonomy <- data.frame(name = rownames(profile), taxonomy = profile$taxonomy) %>% 
  taxa_split(sep = ";") %>% 
  dplyr::select(-species) %>% 
  map_df(\(x) ifelse(grepl("__\\w$", x), "Unknown", gsub("^\\s+[_]+", "", x)))

profile <- profile[,-ncol(profile)]

data <- profile_filter(profile, dplyr::select(group, sample, group), by_group = T, n_group = 1, min_n = 3) %>% 
  profile_transRA()

features <- taxonomy %>% 
  filter(genus == "Blautia") %>% 
  pull(name)

plot_data <- data %>% 
  filter(rownames(.) %in% features) %>% 
  colSums() %>% 
  data.frame(value = .) %>% 
  rownames_to_column("sample") %>% 
  left_join(group, by = "sample") %>% 
  mutate(group = factor(group, c("nonIBD", "CD", "UC")))

ggviolin(plot_data, "group", "value", color = "group", legend = "none", xlab = "", size = 1, width = .5,
         ylab = "Relative abundance(%)", palette = c("#69a7bc","#E69F00","#f08178"),
         outlier.shape = NA, title = "Lloyd-Price et al., 2019") +
  geom_boxplot(aes(color = group), width = .18, size = 1, outlier.shape = NA) + 
  geom_jitter(aes(color = group), width = .2) +
  lims(y = c(NA, 10)) +
  stat_compare_means(comparisons = list(c("nonIBD", "CD"), c("nonIBD", "UC")), label.y = 8.5,
                     tip.length = .01, step.increase = .022, vjust = .9, size = 3) +
  theme(aspect.ratio = 1.4)
# ggsave("g__Blautia.rela_ab.pdf", width = 3.5, height = 4)
