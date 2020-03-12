# set working directory to directory of this script
# setwd(directory_name)
setwd("/home/reto/polybox/ETH/Master_Thesis/Code/Framework/Simulation/Scripts")


################################################################################
## load packages and data
library(tidyverse)
library(ggpubr)
dir_save <- "../../../Data/Simulated/"
plot_dir <- "../../../Plots/Simulated/"

test_results <- readRDS(paste0(dir_save,"test_results_normal_transformation.rds"))
tibplt <- test_results[[2]] %>% 
  mutate(impType = factor(impType, levels = c("cc","pmm","km", "rs","mrl","ppd")),
         simType = factor(simType, levels = c("lm","glm","glmer")),
         cov_dep_cens = factor(cov_dep_cens, levels = c(0,1),labels=c("MCAR","MAR"))) %>%
  filter(error_variance == 0,
         simType == "glmer",
         b1_True == -1,
  )

################################################################################
## do plots
for(cens_rate in c(0.3,0.5,0.7)){
  conds <- expand.grid( transform_fn_tmp=c("identity","sqrt","boxcox_positive","log_positive"),variance_fixeff_tmp = c(0.5),variance_raneff_tmp = c(0.5),n_dat_tmp=c(20))
  plts <- purrr::pmap(conds, function(transform_fn_tmp,variance_fixeff_tmp,variance_raneff_tmp,n_dat_tmp){
    ggplot(tibplt %>% 
             mutate(sim_id = factor(sim_id)) %>% 
             filter(n_dat==n_dat_tmp, 
                    variance_fixeff == variance_fixeff_tmp, 
                    variance_raneff == variance_raneff_tmp, 
                    censoring_rate==cens_rate, 
                    transform_fn==transform_fn_tmp)) + 
      geom_errorbar(aes(x=sim_id,ymin=ci_lower,ymax=ci_upper),size=0.2) + 
      geom_hline(yintercept = -1,size=0.2) +
      geom_point(aes(sim_id,b1),size=0.2) +
      coord_cartesian(ylim=c(-5, 3))  +
      scale_y_continuous(breaks = c(seq(-5,3,by=2),0)) +
      facet_grid(cov_dep_cens ~ impType ) + 
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.grid.major.x = element_line(size=0.2)) +
      labs(title=transform_fn_tmp, x="",y="estimated b1")
  })
  plt_comb <- ggarrange(plotlist = plts,common.legend = TRUE)
  plt_comb <- annotate_figure(plt_comb,top=text_grob(paste0("Sample size = 20, censoring rate = ",cens_rate*100,"%"), color = "black", face = "bold", size = 14))
  ggsave(paste0(plot_dir,"ci_mult_conditions_c",100*cens_rate,"_transformations.png"), plt_comb, width = 40, height = 30, units = "cm")
}


for(cens_rate in c(0.3,0.5,0.7)){
  conds <- expand.grid( transform_fn_tmp=c("identity","sqrt","boxcox_positive","log_positive"),variance_fixeff_tmp = c(0.5),variance_raneff_tmp = c(0.5),n_dat_tmp=c(20))
  plts <- purrr::pmap(conds, function(transform_fn_tmp,variance_fixeff_tmp,variance_raneff_tmp,n_dat_tmp){
    ggplot(tibplt %>% 
             mutate(sim_id = factor(sim_id),
                    cens_eff = n_cens/n_dat) %>% 
             filter(n_dat==n_dat_tmp, variance_fixeff == variance_fixeff_tmp, 
                    variance_raneff == variance_raneff_tmp, 
                    censoring_rate==cens_rate, transform_fn==transform_fn_tmp)) + 
      geom_hline(yintercept = 0.05,size=0.2) +
      geom_point(aes(sim_id,pval),size=1) +
      scale_y_continuous(breaks = c(0,0.05,seq(0.2,1,0.2))) +
      coord_cartesian(ylim=c(0, 1)) +
      facet_grid(cov_dep_cens ~ impType ) + 
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.grid.major.x = element_line(size=0.2)) +
      labs(title=transform_fn_tmp, x="",y="p-value")
  })
  plt_comb <- ggarrange(plotlist = plts,common.legend = TRUE)
  plt_comb <- annotate_figure(plt_comb,top=text_grob(paste0("Sample size = 20, censoring rate = ",cens_rate*100,"%"), color = "black", face = "bold", size = 14))
  ggsave(paste0(plot_dir,"pval_mult_conditions_c",100*cens_rate,"_transformations.png"), plt_comb, width = 40, height = 30, units = "cm")
}

