#### Jinxin Meng, 20240801, 20241113 ####

setwd("F:/proj/20240731_BC_IBD_Zhangyn/submit/code_availability/MBX_HM350_BC/")
pacman::p_load(tidyr, dplyr, tibble, purrr, ggpubr)

#### IPA ####
source("/code/R_func/ms_utilites.R")
library(xcms) %>% suppressMessages()
library(MsExperiment) %>% suppressMessages()

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

# 22P17750007和22P17750008这2个样本在合同F22FTSECKF11912中对应的产品下，
# indole-3-propionic acid (HMDB0002302)的RT也是7.16 min

rt = 7.16

# BC group
BC_ms <- readSRMData("../rawdata/BC.mzML")
BC_ms_parser <- parse_MChromatograms(BC_ms)

BC_ms_data <- BC_ms_parser %>% 
  filter(rt_min < rt & rt_max > rt)

# BC control group
Con_ms <- readSRMData("../rawdata/BC_control.mzML")
Con_ms_parser <- parse_MChromatograms(Con_ms)
Con_ms_data <- Con_ms_parser %>% 
  filter(rt_min < rt & rt_max > rt)

plot_list <- map(Con_ms_data$name, \(x) {
  plot_data <- rbind(data.frame(rtime = rtime(BC_ms[x]), 
                                intensity = intensity(BC_ms[x])) %>% 
                       add_column(class = "BC"),
                     data.frame(rtime = rtime(Con_ms[x]), 
                                intensity = intensity(Con_ms[x])) %>% 
                       add_column(class = "Con")) %>% 
    mutate(rtime = as.numeric(rtime))
  
  ggplot(plot_data, aes(rtime, intensity, group = class)) +
    geom_line(aes(color = class), show.legend = F) +
    geom_vline(xintercept = 7.16, linetype = "dashed") +
    scale_color_manual(values = c("BC" = "#f8766d", "Con" = "#00bfc4")) +
    labs(title = paste0("ms: ", x)) +
    theme_pubr() +
    theme(aspect.ratio = .7) } )
cowplot::plot_grid(plotlist = plot_list, align = "v", nrow = 6)
ggsave("ms_plot.pdf", width = 30, height = 24)

# MS 388
plot_data <- rbind(data.frame(rtime = rtime(BC_ms["388"]), 
                              intensity = intensity(BC_ms["388"])) %>% 
                     add_column(class = "BC"),
                   data.frame(rtime = rtime(Con_ms["388"]), 
                              intensity = intensity(Con_ms["388"])) %>% 
                     add_column(class = "Con")) %>% 
  mutate(rtime = as.numeric(rtime))

options(scipen = 10)

ggplot(plot_data, aes(rtime, intensity, group = class)) +
  ggalt::geom_xspline(aes(color = class), size = 3, spline_shape = 1) +
  scale_color_manual(values = c("BC" = "#f8766d", "Con" = "#00bfc4")) +
  labs(x = 'Retention time (min)', y = 'Intensity') +
  theme_pubr() +
  theme(aspect.ratio = .9) +
  guides(color = guide_legend(position = "right"))
