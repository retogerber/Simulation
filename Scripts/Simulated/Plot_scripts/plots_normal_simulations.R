# set working directory to directory of this script
# setwd(directory_name)
setwd("/home/reto/polybox/ETH/Master_Thesis/Code/Framework/Simulation/Scripts/Simulated/Plot_scripts")


################################################################################
## load packages and data
library(tidyverse)
library(ggpubr)

dir_save <- "../../../Data/Simulated/"
plot_dir <- "../../../Plots/Simulated/"

test_results <- readRDS(paste0(dir_save,"test_results_normal.rds"))
plot_width <- 20
plot_height <- 15


test_results_proc <- test_results[[2]] %>% 
  mutate(impType = factor(impType, levels = c("cc","pmm","km", "rs","mrl","ppd")),
         simType = factor(simType, levels = c("lm","glm","glmer")),
         cov_dep_cens = factor(cov_dep_cens, levels = c(0,1),labels=c("MCAR","MAR")))

################################################################################
## b1 plot
est_plot <- function(tib2plt_n_dat, facetting_var,legend_label){
  facetting_var <- enquo(facetting_var)
  ggplot(data=tib2plt_n_dat) + 
    geom_rect(aes(fill=impType,xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf)) +
    guides(fill=FALSE)+
    scale_fill_manual(values=c(rep(c("lightgrey","grey"),3))) +
    geom_boxplot(aes(1, b1), size=0.2)+
    geom_point(aes(1, b1, color=factor(!!facetting_var)), size=1) + 
    geom_hline(yintercept = -0.5) +
    labs(color = legend_label, y = expression(hat(beta)["1"]), x = "") +
    facet_wrap(vars(!!facetting_var,impType),nrow = 1) +
    theme(axis.text.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(size=0.3),
          panel.spacing = unit(0.1, "lines"),
          axis.ticks.x = element_blank(),
          plot.margin = margin(r=7,b=-10,l=5))+ 
    scale_color_manual(values=RColorBrewer::brewer.pal(5,"Blues")[3:5]  )
}

################################################################################
## mse plot
mse_plot <- function(tib2plt, facetting_var,legend_label,b1_True_val=-0.5){
  mse_fn <- function(x,b1_True_val=b1_True_val){
    ceiling((1/length(x))*sum((x-b1_True_val)^2)*1000)/1000
  }
  facetting_var <- enquo(facetting_var)
  mse_b1_tib <- tib2plt %>% 
    group_by(b1_True,!!facetting_var,impType) %>% 
    select(b1_True,!!facetting_var,impType,b1) %>% 
    summarise(mse=mse_fn(b1,b1_True)) %>% 
    ungroup()
  mse_b1 <- mse_b1_tib %>% 
    select(mse) %>% 
    as_vector()
  mse_limits <- c(min(mse_b1)-0.1*max(mse_b1),max(mse_b1)+0.1*max(mse_b1))
  
  ggplot(data=mse_b1_tib) + 
  geom_rect(aes(fill=impType,xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf)) +
  guides(fill=FALSE)+
  geom_label(aes(x=1,y=mse,label=..y..),size=2,
             label.padding = unit(0.1, "lines"),
             label.r = unit(0.05, "lines")) +
  scale_fill_manual(values=c(rep(c("lightgrey","grey"),3))) +
  labs(color = legend_label, y = expression(italic(MSE)), x = "") +
  facet_wrap(vars(!!facetting_var,impType),nrow = 1) +
  theme(axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size=0.3),
        panel.spacing = unit(0.1, "lines"),
        axis.ticks.x = element_blank(),
        plot.margin = margin(r=5,b=-10,l=5),
        strip.background = element_blank(),
        strip.text = element_blank())+ 
  scale_y_continuous(breaks = c(0,ceiling(mse_limits[2]*10)/10),limits = c(mse_limits[1],ceiling(mse_limits[2]*10)/10)) +
  scale_color_manual(values=RColorBrewer::brewer.pal(5,"Blues")[3:5]  )
}

################################################################################
## p-value plot
pval_plot <- function(tib2plt_n_dat, facetting_var,legend_label){
  facetting_var <- enquo(facetting_var)
  ggplot(data=tib2plt_n_dat) + 
  geom_rect(aes(fill=impType,xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf)) +
  guides(fill=FALSE)+
  scale_fill_manual(values=c(rep(c("lightgrey","grey"),3))) +
  geom_boxplot(aes(1, pval), size=0.2)+
  geom_point(aes(1, pval, color=factor(!!facetting_var)), size=1) + 
  labs(color = "legend_label", y = "p-value", x = "") +
  facet_wrap(vars(!!facetting_var,impType),nrow = 1) +
  theme(axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size=0.3),
        panel.spacing = unit(0.1, "lines"),
        axis.ticks.x = element_blank(),
        plot.margin = margin(r=5,b=-10,l=5),
        strip.background = element_blank(),
        strip.text = element_blank())+ 
  scale_y_continuous(breaks = c(0,0.05,0.2,0.4,0.6,0.8,1),limits = c(0,1)) +
  scale_color_manual(values=RColorBrewer::brewer.pal(5,"Blues")[3:5]  )
}

################################################################################
#### sample size example
cens_rate <- 0.7
tib2plt_n_dat <- test_results_proc %>%   
  filter(error_variance==0,
         variance_fixeff == 0.5, 
         variance_raneff == 0.1,
         simType == "glmer",
         b1_True == -0.5,
         censoring_rate==cens_rate,
         cov_dep_cens=="MCAR", 
         transform_fn =="log_positive"
         ) 

plt_n_dat_b1 <- est_plot(tib2plt_n_dat,n_dat,"") +
  scale_y_continuous(breaks = c(seq(-90,90,by=0.25),-0.5))
plt_n_dat_mse <- mse_plot(tib2plt_n_dat,n_dat,"")
plt_n_dat_pval <-  pval_plot(tib2plt_n_dat,n_dat,"")

plt_n_dat <- ggarrange(plt_n_dat_b1,plt_n_dat_mse,plt_n_dat_pval,
                       common.legend = TRUE,ncol = 1,align = "v",heights = c(5,1,4))
plt_n_dat <- annotate_figure(plt_n_dat,
                             # top=text_grob(paste0("censoring rate = ",cens_rate*100,"%"), color = "black", face = "bold", size = 8,hjust = 5),
                             fig.lab = "Sample size",fig.lab.size = 14,fig.lab.face = "bold")
ggsave(paste0(plot_dir,"simulation_normal_effect_size.png"), plt_n_dat, width = plot_width, height = plot_height, units = "cm")




################################################################################
#### censoring rate example
n_samples <- 20
tib2plt_cens <- test_results_proc %>%   
  filter(error_variance==0,
         variance_fixeff == 0.5, 
         variance_raneff == 0.1,
         simType == "glmer",
         b1_True == -0.5,
         n_dat==n_samples,
         cov_dep_cens=="MCAR", 
         transform_fn =="log_positive"
  )
plt_cens_b1 <- est_plot(tib2plt_cens,censoring_rate,"") +
  scale_y_continuous(breaks = c(seq(-95,90,by=0.25),-0.5))
plt_cens_mse <- mse_plot(tib2plt_cens,censoring_rate,"")
plt_cens_pval <-  pval_plot(tib2plt_cens,censoring_rate,"")

plt_cens <- ggarrange(plt_cens_b1,plt_cens_mse,plt_cens_pval,
                      common.legend = TRUE,ncol = 1,align = "v",heights = c(5,1,4))
plt_cens <- annotate_figure(plt_cens,
                            fig.lab = "Censoring rate",fig.lab.size = 14,fig.lab.face = "bold")
ggsave(paste0(plot_dir,"simulation_normal_effect_cens.png"), plt_cens, width = plot_width, height = plot_height, units = "cm")





################################################################################
#### cov_dep_cens example
cens_rate <- 0.7
n_samples <- 50
tib2plt_cov_dep <- test_results_proc %>%   
  filter(error_variance==0,
         variance_fixeff == 0.5, 
         variance_raneff == 0.1,
         simType == "glmer",
         b1_True == -0.5,
         n_dat==n_samples,
         censoring_rate==cens_rate,
         # cov_dep_cens=="MCAR", 
         transform_fn =="log_positive"
  )
plt_cov_dep_b1 <- est_plot(tib2plt_cov_dep,cov_dep_cens,"") +
  scale_y_continuous(breaks = c(seq(-125,120,by=0.25),-0.5)) 
plt_cov_dep_mse <- mse_plot(tib2plt_cov_dep,cov_dep_cens,"")
plt_cov_dep_pval <-  pval_plot(tib2plt_cov_dep,cov_dep_cens,"")

plt_cov_dep <- ggarrange(plt_cov_dep_b1,plt_cov_dep_mse,plt_cov_dep_pval,
                         common.legend = TRUE,ncol = 1,align = "v",heights = c(5,1,4))
plt_cov_dep <- annotate_figure(plt_cov_dep, 
                               fig.lab = "Censoring type",fig.lab.size = 14,fig.lab.face = "bold")
ggsave(paste0(plot_dir,"simulation_normal_effect_cov_dep.png"), plt_cov_dep, width = plot_width, height = plot_height, units = "cm")



################################################################################
#### transformations example
cens_rate <- 0.7
n_samples <- 20
tib2plt_transform <- test_results_proc %>%   
  filter(error_variance==0,
         variance_fixeff == 0.5, 
         variance_raneff == 0.1,
         simType == "glmer",
         b1_True == -0.5,
         n_dat==n_samples,
         censoring_rate==cens_rate,
         cov_dep_cens=="MCAR"
         # transform_fn =="identity"
  ) %>% 
  mutate(transform_fn=factor(transform_fn,levels = c("identity","boxcox_positive","log_positive"),labels = c("none","boxcox","log")))

plt_transform_b1 <- est_plot(tib2plt_transform,transform_fn,"") +
  scale_y_continuous(breaks = c(seq(-122,120,by=3),-0.5)) 
plt_transform_mse <- mse_plot(tib2plt_transform,transform_fn,"")
plt_transform_pval <-  pval_plot(tib2plt_transform,transform_fn,"")

plt_transform <- ggarrange(plt_transform_b1,plt_transform_mse,plt_transform_pval,
                         common.legend = TRUE,ncol = 1,align = "v",heights = c(5,1,4))
plt_transform <- annotate_figure(plt_transform, 
                               fig.lab = "Data transformation",fig.lab.size = 14,fig.lab.face = "bold")
ggsave(paste0(plot_dir,"simulation_normal_effect_transform.png"), plt_transform, 
       width = plot_width, height = plot_height, units = "cm",dpi=600)



################################################################################
#### variance random effects example
cens_rate <- 0.7
n_samples <- 20
tib2plt_ranvar <- test_results_proc %>%   
  filter(error_variance==0,
         variance_fixeff == 0.5, 
         # variance_raneff == 0.1,
         simType == "glmer",
         b1_True == -0.5,
         n_dat==n_samples,
         censoring_rate==cens_rate,
         cov_dep_cens=="MCAR",
         transform_fn =="log_positive"
  ) 
plt_ranvar_b1 <- est_plot(tib2plt_ranvar,variance_raneff,"") +
  scale_y_continuous(breaks = c(seq(-120,120,by=0.5),-0.5)) 
plt_ranvar_mse <- mse_plot(tib2plt_ranvar,variance_raneff,"")
plt_ranvar_pval <-  pval_plot(tib2plt_ranvar,variance_raneff,"")

plt_ranvar <- ggarrange(plt_ranvar_b1,plt_ranvar_mse,plt_ranvar_pval,
                           common.legend = TRUE,ncol = 1,align = "v",heights = c(5,1,4))
plt_ranvar <- annotate_figure(plt_ranvar, 
                                 fig.lab = "Variance of random effect",fig.lab.size = 14,fig.lab.face = "bold")
ggsave(paste0(plot_dir,"simulation_normal_effect_ranvar.png"), plt_ranvar, width = plot_width, height = plot_height, units = "cm")

################################################################################
#### variance fixed effects example
cens_rate <- 0.7
n_samples <- 20
tib2plt_fixeff <- test_results_proc %>%   
  filter(error_variance==0,
         # variance_fixeff == 0.5, 
         variance_raneff == 0.1,
         simType == "glmer",
         b1_True == -0.5,
         n_dat==n_samples,
         censoring_rate==cens_rate,
         cov_dep_cens=="MCAR",
         transform_fn =="log_positive"
  ) 
plt_fixeff_b1 <- est_plot(tib2plt_fixeff,variance_fixeff,"") +
  scale_y_continuous(breaks = c(seq(-120,120,by=0.25),-0.5)) 
plt_fixeff_mse <- mse_plot(tib2plt_fixeff,variance_fixeff,"")
plt_fixeff_pval <-  pval_plot(tib2plt_fixeff,variance_fixeff,"")

plt_fixeff <- ggarrange(plt_fixeff_b1,plt_fixeff_mse,plt_fixeff_pval,
                        common.legend = TRUE,ncol = 1,align = "v",heights = c(5,1,4))
plt_fixeff <- annotate_figure(plt_fixeff, 
                              fig.lab = "Variance of fixed effect",fig.lab.size = 14,fig.lab.face = "bold")
ggsave(paste0(plot_dir,"simulation_normal_effect_fixeff.png"), plt_fixeff, width = plot_width, height = plot_height, units = "cm")

################################################################################
#### regression coeff example
cens_rate <- 0.7
n_samples <- 20
tib2plt_reg_coef <- test_results_proc %>%   
  filter(error_variance==0,
         variance_fixeff == 0.5,
         variance_raneff == 0.1,
         simType == "glmer",
         # b1_True == -1.5,
         n_dat==n_samples,
         censoring_rate==cens_rate,
         cov_dep_cens=="MCAR",
         transform_fn =="log_positive"
  ) 
plt_reg_coef_b1 <- est_plot(tib2plt_reg_coef,b1_True,"") +
  geom_hline(yintercept = -1.5) +
  geom_hline(yintercept = -1) +
  scale_y_continuous(breaks = c(seq(-120,120,by=3),-1.5,-0.5,-1)) 
plt_reg_coef_mse <- mse_plot(tib2plt_reg_coef,b1_True,"")
plt_reg_coef_pval <-  pval_plot(tib2plt_reg_coef,b1_True,"")

plt_reg_coef <- ggarrange(plt_reg_coef_b1,plt_reg_coef_mse,plt_reg_coef_pval,
                        common.legend = TRUE,ncol = 1,align = "v",heights = c(5,1,4))
plt_reg_coef <- annotate_figure(plt_reg_coef, 
                              fig.lab = "Regression coefficient",fig.lab.size = 14,fig.lab.face = "bold")
ggsave(paste0(plot_dir,"simulation_normal_effect_reg_coef.png"), plt_reg_coef, width = plot_width, height = plot_height, units = "cm")


################################################################################
#### simtype example
cens_rate <- 0.7
n_samples <- 20
tib2plt_simtype <- test_results_proc %>%   
  filter(error_variance==0,
         variance_fixeff == 0.5, 
         variance_raneff == 0.1,
         # simType == "glmer",
         b1_True == -0.5,
         n_dat==n_samples,
         censoring_rate==cens_rate,
         cov_dep_cens=="MCAR",
         transform_fn =="log_positive"
  )  %>% 
  mutate(simType=factor(simType,levels = c("lm","glm","glmer"),labels = c("LM","GLM","GLMM")))

plt_simtype_b1 <- est_plot(tib2plt_simtype,simType,"") +
  scale_y_continuous(breaks = c(seq(-125,120,by=0.25),-0.5))
plt_simtype_mse <- mse_plot(tib2plt_simtype,simType,"")
plt_simtype_pval <-  pval_plot(tib2plt_simtype,simType,"")

plt_simtype <- ggarrange(plt_simtype_b1,plt_simtype_mse,plt_simtype_pval,
                         common.legend = TRUE,ncol = 1,align = "v",heights = c(5,1,4))
plt_simtype <- annotate_figure(plt_simtype, 
                               fig.lab = "Regression type",fig.lab.size = 14,fig.lab.face = "bold")
ggsave(paste0(plot_dir,"simulation_normal_effect_simtype.png"), plt_simtype, 
       width = plot_width, height = plot_height, units = "cm",dpi=600)




