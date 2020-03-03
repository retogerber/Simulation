library(tidyverse)
# library(scales)
library(cowplot)
library(ggpubr)
# dir_save <- "/home/reto/polybox/ETH/Master_Thesis/Portmacquarie/tmp/"
dir_save <- "/home/reto/polybox/ETH/Master_Thesis/Code/Framework/Simulation/"
plot_dir <- "/home/reto/polybox/ETH/Master_Thesis/plots/Simulations/normal/"
# test_results <- readRDS(paste0(dir_save,"test_results_normal_new_new.rds"))
test_results <- readRDS(paste0(dir_save,"test_results_normal_glmer_only.rds"))
# test_results_params_ind <- map(test_results_params, ~ is.data.frame(.x))
# errors <- (test_results_params[!unlist(test_results_params_ind)])
# errors_uniq <- unique(errors)
# errors_ind <- map(errors_uniq, function(y) unlist(map(errors, ~identical(.x,y))))
# map(1:16,~table(errors_ind[[.x]]))
# errors_uniq[[1]]
# 
# 
# params <- bind_rows(test_results_params[unlist(test_results_params_ind)])
# test_results <- list(simulations_data_processor(params), params)
# 
# grouping_quo <- dplyr::quos(simType, formula, impType, n_dat, censoring_rate, 
#                             b0_True, b1_True,b2_True,mi_rep, cov_dep_cens, error_variance)
# 
# tmp <- expand.grid(grouping_quo,grouping_quo)
# size <- test_results[[2]] %>% tally() %>% as_vector()
# combs <- pmap(tmp, function(Var1,Var2){
#   tmptib <- test_results[[2]] %>%
#     dplyr::group_by(!!Var1,!!Var2) %>%
#     dplyr::tally()
#   tmptibsize <- dim(tmptib)[1]
#   tmptib <- tmptib %>% 
#     mutate(n_rel=n/size, n_lower = n < 0.9*size/tmptibsize,  n_upper = n>1.1*size/tmptibsize)
#   return(tmptib)
# })
# test_results <- out_ls

tib2plt <- test_results[[1]] %>% 
  filter(parameter == "b1") %>% 
  mutate(impType = factor(impType, levels = c("cc","pmm","km", "rs","mrl","ppd")),
         simType = factor(simType, levels = c("lm","glm","glmer")),
         absolute_difference = abs(mean-b1_True), 
         b1_True = as.factor(b1_True), 
         n_dat= as.factor(n_dat),
         cov_dep_cens = factor(cov_dep_cens, levels = c(0,1),labels=c("MCAR","MAR"))) %>% 
  filter(error_variance==0,
         # b2_True==0.5, 
         variance_fixeff == 0.5, 
         variance_raneff == 0.1
         )

# # censoring rate
# ggplot(data=tib2plt) + geom_point(aes(censoring_rate,proportion_censored, color=n_dat))

# mse vs. true parameter, color by censoring rate
# impTypes <- c("mrl","km","rs","cc")
# simulation_conditions <- tib2plt %>% select(simType,formula) %>% unique()
# parameters_to_test <- expand.grid(list(impTypes=impTypes,rownames=rownames(as.data.frame(simulation_conditions))))
# parameters_to_test <- as_tibble(parameters_to_test) %>% 
#   inner_join(simulation_conditions %>% mutate(rownames=rownames(as.data.frame(simulation_conditions)))) %>% 
#   select(-rownames)

# mse vs. true parameter, color by censoring rate
plt1 <- ggplot(data=tib2plt %>% filter(n_dat==20)) + 
  # geom_jitter(aes(b1_True, mse, color=proportion_censored,shape=n_dat), size=2, width = 0.3) + 
  geom_point(aes(b1_True, mse, color=proportion_censored,shape=n_dat), size=2) + 
  ylim(c(0,10)) +
  labs(color = "censoring\nrate", shape = "sample\nsize", y = "MSE", x = "b1") +#, title = "MSE vs. b1 for large error variance") +
  facet_grid(simType ~ impType + cov_dep_cens) + 
  # scale_colour_gradientn(colours=rainbow(3))
  scale_colour_gradient(low = "green", high = "black")
plt1
ggsave(paste0(plot_dir,"mse_vs_parameter_colby_censoring_size20.png"), plt1, width = 20, height = 20, units = "cm")


plt2 <- ggplot(data=tib2plt %>% 
                 filter( n_dat==50)) + 
  # geom_jitter(aes(b1_True, mse, color=proportion_censored,shape=n_dat), size=2, width = 0.3) + 
  geom_point(aes(b1_True, mse, color=proportion_censored,shape=n_dat), size=2) + 
  ylim(c(0,10)) +
  labs(color = "censoring\nrate", shape = "sample\nsize", y = "MSE", x = "b1") +#, title = "MSE vs. b1 for large error variance") +
  facet_grid(simType ~ impType + cov_dep_cens) + 
  # scale_colour_gradientn(colours=rainbow(3))
  scale_colour_gradient(low = "green", high = "black")
plt2
ggsave(paste0(plot_dir,"mse_vs_parameter_colby_censoring_size100.png"), plt2, width = 20, height = 20, units = "cm")


# mse vs. true parameter, color by log snr
ggplot(data=tib2plt %>% 
         filter(n_dat==50, simType!="glmer") %>% 
         mutate(absolute_difference = abs(mean-b1_True), b1_True = as.factor(b1_True), 
                snr_log = log(snr), n_dat= as.factor(n_dat))) + 
  geom_jitter(aes(b1_True, mse, color=snr_log,shape=n_dat), size=2) + 
  ylim(c(0,10)) + 
  labs(color = "log(snr)", shape = "sample\nsize", y = "MSE", x = "b1", title = "MSE") +
  facet_grid(vars(simType), vars(impType)) + 
  scale_colour_gradient(low = "green", high = "black")

# mse vs. true parameter, color mi reps
ggplot(data=tib2plt %>% 
         filter(n_dat==50) %>% 
         mutate(mi_rep = as.factor(mi_rep))) + 
  geom_jitter(aes(b1_True, mse, color=mi_rep,shape=n_dat), size=2) + 
  ylim(c(0,10)) +
  labs(color = "mi rep", shape = "sample\nsize", y = "MSE", x = "b1", title = "MSE") +
  facet_grid(vars(simType), vars(impType)) 


# data size vs. mse
ggplot(data=tib2plt %>% 
         # filter(proportion_censored>0.5) %>% 
         mutate(absolute_difference = abs(mean-b1_True), b1_True = as.factor(b1_True), 
                mi_rep = as.factor(mi_rep), n_dat=as.factor(n_dat))) + 
  geom_jitter(aes(n_dat, mse, color=proportion_censored), size=2) + 
  ylim(c(0,50)) +
  labs(color = "censoring\nrate", y = "MSE", x = "sample size", title = "MSE") +
  facet_grid(cols = vars(impType)) 



# proportion censored vs. mse
ggplot(data=tib2plt %>% 
         # filter(proportion_censored>0.5) %>% 
         mutate(absolute_difference = abs(mean-b1_True), b1_True = as.factor(b1_True), 
                mi_rep = as.factor(mi_rep), n_dat=as.factor(n_dat))) + 
  geom_jitter(aes(proportion_censored, mse, color=n_dat), size=2) + 
  ylim(c(0,100)) +
  labs(color = "sample\nsize", y = "MSE", x = "proportion censored", title = "MSE") +
  facet_grid(cols = vars(impType), rows = vars(simType)) 

# p value vs. true parameter
plt3 <- ggplot(data=tib2plt %>% 
         filter( n_dat==50)) + 
  geom_point(aes(b1_True, pval, color=proportion_censored,shape=n_dat), size=2, width = 0.3) + 
  ylim(c(0,1)) +
  labs(color = "censoring\nrate", shape = "sample\nsize", y = "pval", x = "b1", title = "P-value") +
  facet_grid(simType~impType+ cov_dep_cens) + 
  # scale_colour_gradientn(colours=rainbow(3))
  scale_colour_gradient(low = "green", high = "black") +
  geom_hline(yintercept = 0.05)
plt3
ggsave(paste0(plot_dir,"pval_vs_parameter_colby_censoring_size100.png"), plt3, width = 20, height = 20, units = "cm")

# p value vs. true parameter
plt4 <- ggplot(data=tib2plt %>% 
                 filter( n_dat==20)) + 
  geom_point(aes(b1_True, pval, color=proportion_censored,shape=n_dat), size=2, width = 0.3) + 
  ylim(c(0,1)) +
  labs(color = "censoring\nrate", shape = "sample\nsize", y = "pval", x = "b1", title = "P-value") +
  facet_grid(simType~impType+ cov_dep_cens) + 
  # scale_colour_gradientn(colours=rainbow(3))
  scale_colour_gradient(low = "green", high = "black") +
  geom_hline(yintercept = 0.05)
plt4
ggsave(paste0(plot_dir,"pval_vs_parameter_colby_censoring_size20.png"), plt4, width = 20, height = 20, units = "cm")

# number of NA in b1 depends on n_dat
# params_mean <- test_results[[2]]  %>% 
#   # select(-iter) %>% 
#   group_by(simType, impType, n_dat,censoring_rate) %>% 
#   summarise_all(mean) %>%  
#   mutate(proportion_censored = n_cens/n_dat) %>% 
#   select(-starts_with("Var"),-n_cens)
# params_mean %>% mutate(b1_is_NA = ifelse(is.na(b1),TRUE,FALSE)) %>%  group_by(n_dat,b1_is_NA) %>% tally()



tibplt <- test_results[[2]] %>% 
  # filter(parameter == "b1") %>% 
  mutate(impType = factor(impType, levels = c("cc","pmm","km", "rs","mrl","ppd")),
         simType = factor(simType, levels = c("lm","glm","glmer")),
  #        absolute_difference = abs(mean-b1_True), 
  #        b1_True = as.factor(b1_True), 
  #        n_dat= as.factor(n_dat),
         cov_dep_cens = factor(cov_dep_cens, levels = c(0,1),labels=c("MCAR","MAR"))) %>%
  filter(error_variance == 0,
         # variance_fixeff == 1, 
         # variance_raneff == 0.5, 
         # mi_rep == 100,
         simType == "glmer",
         # b0_True == -2,
         b1_True == -1,
         # b2_True == 1,
         # n_dat == 20,
         # censoring_rate == 0.3
         # cov_dep_cens == "TRUE"
  )

# tibplt %>% filter(sim_id==1)
for(cens_rate in c(0.3,0.5,0.7)){
  # unique(tibplt$transform_fn)
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
      # scale_x_continuous(breaks = c(0.01,0.05,0.1)) +
      facet_grid(cov_dep_cens ~ impType ) + 
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.grid.major.x = element_line(size=0.2)) +
      labs(title=transform_fn_tmp, x="",y="estimated b1")
      # scale_x_continuous(breaks = c(0.01,0.05,0.1)) +
      # facet_grid(cov_dep_cens ~ impType ) +
      # labs(title = paste0("n",n_dat_tmp,"_vf",variance_fixeff_tmp,"_vr",variance_raneff_tmp,"_",transform_fn_tmp))
  })
  plt_comb <- ggarrange(plotlist = plts,common.legend = TRUE)
  plt_comb <- annotate_figure(plt_comb,top=text_grob(paste0("Sample size = 20, censoring rate = ",cens_rate*100,"%"), color = "black", face = "bold", size = 14))
  # plt_comb <- ggarrange(plotlist = plts)
  ggsave(paste0(plot_dir,"ci_mult_conditions_c",100*cens_rate,".png"), plt_comb, width = 40, height = 30, units = "cm")
}


# tibplt %>% filter(sim_id==1)
for(cens_rate in c(0.3,0.5,0.7)){
  conds <- expand.grid( transform_fn_tmp=c("identity","sqrt","boxcox_positive","log_positive"),variance_fixeff_tmp = c(0.5),variance_raneff_tmp = c(0.5),n_dat_tmp=c(20))
  plts <- purrr::pmap(conds, function(transform_fn_tmp,variance_fixeff_tmp,variance_raneff_tmp,n_dat_tmp){
    ggplot(tibplt %>% 
             mutate(sim_id = factor(sim_id),
                    cens_eff = n_cens/n_dat) %>% 
             filter(n_dat==n_dat_tmp, variance_fixeff == variance_fixeff_tmp, 
                    variance_raneff == variance_raneff_tmp, 
                    censoring_rate==cens_rate, transform_fn==transform_fn_tmp)) + 
      # geom_errorbar(aes(x=sim_id,ymin=ci_lower,ymax=ci_upper),size=0.2) + 
      # geom_hline(yintercept = -1,size=0.2) +
      # geom_point(aes(sim_id,b1),size=0.2) +
      # coord_cartesian(ylim=c(-5, 3)) +
      # geom_boxplot(aes(b1_True,pval),size=0.2) +
      geom_hline(yintercept = 0.05,size=0.2) +
      geom_point(aes(sim_id,pval),size=1) +
      scale_y_continuous(breaks = c(0,0.05,seq(0.2,1,0.2))) +
      coord_cartesian(ylim=c(0, 1)) +
      # scale_x_continuous(breaks = c(0.01,0.05,0.1)) +
      facet_grid(cov_dep_cens ~ impType ) + 
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.grid.major.x = element_line(size=0.2)) +
      labs(title=transform_fn_tmp, x="",y="p-value")
      # labs(title = paste0("n",n_dat_tmp,"_vf",variance_fixeff_tmp,"_vr",variance_raneff_tmp,"_",transform_fn_tmp))
  })
  plt_comb <- ggarrange(plotlist = plts,common.legend = TRUE)
  plt_comb <- annotate_figure(plt_comb,top=text_grob(paste0("Sample size = 20, censoring rate = ",cens_rate*100,"%"), color = "black", face = "bold", size = 14))
  ggsave(paste0(plot_dir,"pval_mult_conditions_c",100*cens_rate,".png"), plt_comb, width = 40, height = 30, units = "cm")
}



# plt5 <- ggplot(tibplt %>% 
#                  mutate(rep_id = factor(ci_lower,labels = seq_along(ci_lower)) ,sim_id = factor(sim_id)) %>% 
#                  filter(n_dat==20, variance_fixeff == 1, variance_raneff == 0.5)) + 
#   geom_errorbar(aes(x=sim_id,ymin=ci_lower,ymax=ci_upper)) + 
#   geom_hline(yintercept = -1) +
#   geom_point(aes(sim_id,b1)) +
#   coord_cartesian(ylim=c(-5, 3)) +
#   # scale_x_continuous(breaks = c(0.01,0.05,0.1)) +
#   facet_grid(cov_dep_cens ~ impType ) 
# plt5
# ggsave(paste0(plot_dir,"ci_one_condition_size100.png"), plt5, width = 20, height = 20, units = "cm")


plt6 <- ggplot(tibplt %>% mutate(b1_True = as.factor(b1_True) ,sim_id = factor(sim_id))) + 
  # geom_errorbar(aes(x=sim_id,ymin=ci_lower,ymax=ci_upper)) + 
  geom_boxplot(aes(b1_True,pval)) +
  geom_hline(yintercept = 0.05) +
  geom_jitter(aes(b1_True,pval)) +
  coord_cartesian(ylim=c(0, 1)) +
  # scale_x_continuous(breaks = c(0.01,0.05,0.1)) +
  facet_grid(cov_dep_cens ~ impType )
plt6
ggsave(paste0(plot_dir,"pval_one_condition_size100.png"), plt6, width = 20, height = 20, units = "cm")








plt7 <- ggplot(tibplt %>% mutate(b1_True = as.factor(b1_True) ,sim_id = factor(sim_id))) + 
  # geom_errorbar(aes(x=sim_id,ymin=ci_lower,ymax=ci_upper)) + 
  geom_boxplot(aes(b1_True,mean_mse_fits)) +
  # geom_hline(yintercept = 0.05) +
  geom_jitter(aes(b1_True,mean_mse_fits)) +
  # coord_cartesian(ylim=c(0, 0.1)) +
  # scale_x_continuous(breaks = c(0.01,0.05,0.1)) +
  facet_grid(cov_dep_cens ~ impType )
plt7
ggsave(paste0(plot_dir,"mean_mse_one_condition_size100.png"), plt7, width = 20, height = 20, units = "cm")
