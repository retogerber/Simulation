# set working directory to directory of this script
# setwd(directory_name)
setwd("/home/reto/polybox/ETH/Master_Thesis/Code/Framework/Simulation/Scripts/Simulated/Plot_scripts")


################################################################################
## load packages and data
library(tidyverse)
library(ggpubr)

dir_save <- "../../../Data/Simulated/"
plot_dir <- "../../../Plots/Simulated/"

test_results <- readRDS(paste0(dir_save,"test_results_normal_null.rds"))
plot_width <- 20
plot_height <- 15


test_results_proc <- test_results[[2]] %>% 
  mutate(impType = factor(impType, levels = c("cc","pmm","km", "rs","mrl","ppd")),
         simType = factor(simType, levels = c("lm","glm","glmer")),
         cov_dep_cens = factor(cov_dep_cens, levels = c(0,1),labels=c("MCAR","MAR")))

################################################################################
## p-value plot
pval_plot <- function(tib2plt_n_dat,facetting_var,legend_label){
  facetting_var <- enquo(facetting_var)
  ggplot(data=tib2plt_n_dat) + 
    geom_rect(aes(fill=impType,xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf)) +
    guides(fill=FALSE)+
    scale_fill_manual(values=c(rep(c("lightgrey","grey"),3))) +
    # geom_boxplot(aes(1, pval), size=0.2)+
    # geom_violin(aes(1, pval), size=0.2)+
    geom_histogram(aes(!!facetting_var),fill=RColorBrewer::brewer.pal(5,"Blues")[5],bins=10)+
    geom_hline(yintercept = 100)+
    # geom_density(aes(pval))+
    # geom_point(aes(1, pval, color=factor(!!facetting_var)), size=1) + 
    labs(color = "legend_label", y = "", x = "") +
    facet_wrap(n_dat ~ impType,nrow = 1) +
    theme(
      # axis.text.x = element_blank(),
          # panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(size=0.3),
          panel.spacing = unit(0.1, "lines"),
          # axis.ticks.x = element_blank(),
          plot.margin = margin(t=15,r=5,b=-10,l=5))+ 
    scale_x_continuous(breaks = c(0,0.5,1))+
    # scale_y_continuous(limits = c(0,300)) +
    # coord_cartesian(ylim=c(0,300))+
    scale_color_manual(values=RColorBrewer::brewer.pal(5,"Blues")[3:5]  )
}

################################################################################
#### regression coeff example
cens_rate <- 0.7
# n_samples <- 20
tib2plt_reg_coef <- test_results_proc %>%   
  filter(error_variance==0,
         variance_fixeff == 0.5,
         variance_raneff == 0.1,
         simType == "glmer",
         # b1_True == -1.5,
         # n_dat==n_samples,
         censoring_rate==cens_rate,
         cov_dep_cens=="MCAR",
         transform_fn =="log_positive"
  ) 
plt_reg_coef_pval_n10 <-  pval_plot(tib2plt_reg_coef %>% filter(n_dat==10),pval,"") + 
  theme(plot.margin = margin(t=15,r=5,b=-10,l=0),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  scale_fill_manual(values=c(rep(c("grey","lightgrey"),3)))
plt_reg_coef_pval_n20 <-  pval_plot(tib2plt_reg_coef %>% filter(n_dat==20),pval,"")  + 
  theme(plot.margin = margin(t=0,r=5,b=-10,l=5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
plt_reg_coef_pval_n50 <-  pval_plot(tib2plt_reg_coef %>% filter(n_dat==50),pval,"")  + 
  theme(plot.margin = margin(t=0,r=5,b=-10,l=5))
plt_empty <- ggplot(tib2plt_reg_coef) + geom_blank() + theme_bw() + 
  theme(line = element_blank(),rect = element_blank(),plot.margin = margin(t=15,r=0,b=-10,l=0))

plt_reg_coef <- ggarrange(ggarrange(plt_empty,plt_reg_coef_pval_n10,widths = c(1,5)),plt_reg_coef_pval_n20,plt_reg_coef_pval_n50,
                          common.legend = TRUE,ncol = 1,align = "v",heights = c(1,1,1))
plt_reg_coef <- annotate_figure(plt_reg_coef,
                                fig.lab = "P-value distribution under the Null model",fig.lab.size = 14,fig.lab.face = "bold")
plt_reg_coef
ggsave(paste0(plot_dir,"simulation_normal_effect_reg_coef_null.png"), plt_reg_coef, width = plot_width, height = plot_height, units = "cm")

################################################################################
