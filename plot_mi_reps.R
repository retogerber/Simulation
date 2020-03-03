# devtools::load_all("/home/reto/polybox/ETH/Master_Thesis/Code/censcyt")
# # simulations_data_processor(out_ls[[2]])
# out_ls <- run_simulations_wrapper(reps = 10,
#                                   nr_cores = 1,
#                                   n_obs_c = c(20),
#                                   simDataTypes = c("glmer"),
#                                   formulas = list(glmer = list(formula(Y~Surv(X,I) + Z + (1|R)))),
#                                   n_levels_fixeff = 2,
#                                   n_levels_raneff = NULL,
#                                   betas = list(b0=c(-2),b1=c(1),b2=c(1)),
#                                   # betas = list(b0=c(0),b1=c(0.5),b2=c(0)),
#                                   censoring_rates = c(0.3),
#                                   # censoring_rates = c(0.3),
#                                   imputation_types = c("rs"),
#                                   run_with_covariate_dependent_censoring = FALSE,
#                                   mi_reps = c(2,3,4,5,10),
#                                   error_variance = c(0),
#                                   variance_fixeff = c(0.5),
#                                   variance_raneff = c(0.5),
#                                   verbose = TRUE, seed = 1234)
# number_of_clusters = 5)#,
# size_of_random_subset_for_testing = 10)

# out_ls[[2]] %>% 
#   ggplot() +
#   geom_point(aes(mi_rep,pval))

devtools::load_all("/home/retger/FlowCap/scripts/censcyt")
# devtools::load_all("/home/reto/polybox/ETH/Master_Thesis/Code/censcyt")
library(ggplot2)
library(tidyverse)
set.seed(123)
formula_cens_glmm <- formula(y~Surv(X,I)+z+(1|r))
censoring_params <- censcyt::scales_for_censoring(c(0.3,0.5,0.7),log_ratio_val = c(0,0.4))
# censoring_params <- censoring_params %>% filter(censoring==0.7, log_ratio==0)
censoring_list <- list(censoring = censoring_params$censoring, 
                       log_ratio = censoring_params$log_ratio, 
                       C1 = censoring_params$C1, 
                       C2 = censoring_params$C2)
mi_reps <- c(2,3,4,5,7,10,20,30,40,50,75,100,125,150,200,300,400,500)
# mi_reps <- c(2,5,10)
methods_ests <- c("rs","km","mrl","pmm","ppd")
# methods_ests <- c("ppd")
res <- par_map(1:6,function(x){
  data_sim <- simulate_data(n = 50, formula = formula_cens_glmm, type = "glmer",
                            error_variance = 0,
                            b = c(-2,1,1),
                            weibull_params = list(X = list(shape = 0.5, scale = 0.25),
                                                  C = list(shape = 1, scale = censoring_list$C1[x])),
                            variance_fixeff = 0.5,
                            variance_raneff = 0.5,
                            n_levels_raneff = NULL,
                            n_levels_fixeff = 2,
                            censoring_dependent_on_covariate = ifelse(censoring_list$log_ratio[x]==0,FALSE,TRUE),
                            weibull_params_covariate_dependent_censoring = list(shape = 0.5, scale = censoring_list$C2[x]))
  outs <- purrr::map(methods_ests,function(method_est){
    outs <- purrr::map(mi_reps, function(mi_rep){
      set.seed(123)
      tmptime <- system.time({
      tmp <- conditional_multiple_imputation(data_sim,formula = formula_cens_glmm,
                                             regression_type = "glmer",
                                             repetitions = mi_rep,
                                             method_est = method_est,
                                             weights = data_sim$size_tot,
                                             contrasts = c(0,1,0),
                                             family = "binomial")
      tmp_out <- summary(mice::pool(tmp[["fits"]]))
      })
      return(data.frame(estimate = tmp_out$estimate[2],
                        std_err = tmp_out$std.error[2], 
                        df = tmp_out$df[2],
                        pval=tmp_out$p.value[2],
                        mi_rep = mi_rep,
                        sim_id = x,
                        method_est=method_est,
                        user_cpu = tmptime[1],
                        system_cpu = tmptime[2],
                        elapsed_time = tmptime[3]))
    # },mc.cores=length(mi_reps))
  })
    return(bind_rows(outs))
  # },mc.cores=length(methods_ests))
  })
  return(bind_rows(outs))
},mc.cores=6)
# })
# res
saveRDS(res,"/home/retger/FlowCap/results/Simulations/res_plot_mi_reps.rds")



################################################################################
##locale
#################################
res <- readRDS("/home/reto/polybox/ETH/Master_Thesis/Portmacquarie/FlowCap/results/Simulations/res_plot_mi_reps.rds")

res_df <- bind_rows(res) %>% 
  rename("p-value":=pval)
res_long <- pivot_longer(data=res_df,cols=c("estimate","std_err","df","p-value","user_cpu","system_cpu","elapsed_time"))  %>% 
  arrange(method_est,mi_rep) %>% 
  mutate(sim_cond = factor(sim_id,levels = 1:6,labels=c("c0.3-MCAR","c0.5-MCAR","c0.7-MCAR","c0.3-MAR","c0.5-MAR","c0.7-MAR")))
# res_df %>% arrange(mi_rep)


plt0 <- ggplot(expand.grid(censoring=c(0.3,0.5,0.7),mdm = c("MCAR","MAR"))) +
  geom_tile(aes(x = factor(censoring), y = factor(mdm),fill = factor(1:6))) +
  # scale_fill_gradientn(colours=c(“blue”,“white”,“red”)) +
  scale_fill_manual(values = c(RColorBrewer::brewer.pal(4,"Greens")[-1],RColorBrewer::brewer.pal(4,"Blues")[-1])) +
  scale_x_discrete(position = "top") +
  theme_bw() +
  theme(line = element_blank(),rect = element_blank(),legend.position = "none",axis.text = element_text(size=10)) +
  labs(x=element_blank(),y=element_blank())

plt1 <- ggplot(res_long %>% filter(name %in% c("p-value"))) +
  geom_point(aes(mi_rep,value,color=sim_cond))+
  geom_line(aes(mi_rep,value,color=sim_cond)) +
  scale_color_manual(values = c(RColorBrewer::brewer.pal(4,"Greens")[-1],RColorBrewer::brewer.pal(4,"Blues")[-1])) +
  scale_x_continuous(breaks = c(0,mi_reps[-c(1,2,3,4,5,6,7,8,9,11,13)])) +
  scale_y_continuous(breaks = c(0,0.05,0.25,0.5,0.75,1)) +
  facet_wrap(name~method_est,ncol = 5) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(size=0.2),
        panel.grid.major.x = element_line(size=0.3),
        legend.position = "none") +
  guides(color=guide_legend(nrow=2,byrow = TRUE)) +
  labs(x="",y="",color="Simulation\nCondition")
# ggsave("/home/reto/polybox/ETH/Master_Thesis/plots/Simulations/pval_mi_reps.png",plt1, width = 12)

plt2 <- ggplot(res_long %>% filter(name %in% c("estimate"))) +
  geom_point(aes(mi_rep,value,color=sim_cond))+
  geom_line(aes(mi_rep,value,color=sim_cond)) +
  scale_color_manual(values = c(RColorBrewer::brewer.pal(4,"Greens")[-1],RColorBrewer::brewer.pal(4,"Blues")[-1])) +
  scale_x_continuous(breaks = c(0,mi_reps[-c(1,2,3,4,5,6,7,8,9,11,13)])) +
  scale_y_continuous(breaks = c(0,2,4,6)) +
  facet_wrap(name~method_est,ncol = 5) +
  geom_hline(yintercept = 1) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(size=0.2),
        panel.grid.major.x = element_line(size=0.3),
        legend.position = "none") +
  # guides(color=guide_legend(nrow=2,byrow = TRUE)) +
  labs(x="Multiple Imputation repetitions",y="",color="Simulation\nCondition")
library(ggpubr)
plt_empty <- ggplot(res_long) + geom_blank() + theme_bw() + theme(line = element_blank(),rect = element_blank())
# plt_comb <- ggarrange(plt1,plt2,plt_empty,plt0,plt_empty,nrow=3,ncol = 3, widths = c(5,5,2,1,2),heights = c(1,1,1,1,1))
plt_comb <- ggarrange(ggarrange(plt_empty,plt0,plt_empty,ncol = 3),plt1,plt2,nrow=3, heights = c(1,5,5))
plt_comb
ggsave("/home/reto/polybox/ETH/Master_Thesis/plots/Simulations/pval_estimate_mi_reps.png",plt_comb, width = 12,height = 8)




ggplot(res_long %>% filter(name %in% c("elapsed_time"))) +
  geom_point(aes(mi_rep,value,color=sim_cond))+
  geom_line(aes(mi_rep,value,color=sim_cond)) +
  scale_color_manual(values = c(RColorBrewer::brewer.pal(4,"Greens")[-1],RColorBrewer::brewer.pal(4,"Blues")[-1])) +
  scale_x_continuous(breaks = c(0,mi_reps[-c(1,2,3,4,5)])) +
  facet_wrap(name~method_est,ncol = 5) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(size=0.2),
        panel.grid.major.x = element_line(size=0.3),
        legend.position = "top") +
  guides(color=guide_legend(nrow=2,byrow = TRUE)) +
  labs(x="",y="",color="Simulation\nCondition")
