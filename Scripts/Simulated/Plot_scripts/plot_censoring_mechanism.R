# set working directory to directory of this script
# setwd(directory_name)
setwd("/home/reto/polybox/ETH/Master_Thesis/Code/Framework/Simulation/Scripts/Simulated/Run_scripts")


################################################################################
## load packages and data
library(tidyverse)
library(ggpubr)
library(censcyt)

plot_dir <- "../../../Plots/Simulated/"
dir_save <- "../../../Data/Simulated/"

################################################################################
## plots

dfplt <- readRDS(paste0(dir_save,"censoring_mechanism_df.rds"))

plt <- ggplot(dfplt) + 
  geom_hline(aes(yintercept=mean_censoring))+
  geom_line(aes(levels_begin,density,color="density"),size=1) +
  geom_line(aes(levels_begin,density*(1-censrate),color="density*(1-censoring rate)"),size=1) +
  geom_smooth(aes(levels_begin,censrate,color="censoring rate"),span=0.2,level=0) +
  geom_point(aes(levels_begin,censrate,color="censoring rate")) +
  labs(x="x",y="censoring rate / density", colour = "") + 
  facet_wrap(transform_fn~cov_dep_cens,scales = "free_x",nrow = 3) + 
  scale_y_continuous(breaks = seq(0,1,0.1)) +
  coord_cartesian(ylim = c(0,1)) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_line(size=0.2)) 

plt
ggsave(paste0(plot_dir,"censoring_mechanism_comparison.png"),plt, width = 12,height=12)
