# set working directory to directory of this script
# setwd(directory_name)
setwd("/home/reto/polybox/ETH/Master_Thesis/Code/Framework/Simulation/Scripts/Simulated/Run_scripts")


################################################################################
## load packages 
library(censcyt)
library(tidyverse)
library(ggpubr)
plot_dir <- "../../../Plots/Simulated/"


################################################################################
## data simulations and plots
set.seed(123)
formula_cens_glmm <- formula(y~Surv(X,I)+z+(1|r))
censoring_params <- censcyt::scales_for_censoring(c(0.3,0.5,0.7),log_ratio_val = c(0,0.4))
censoring_params <- censoring_params %>% filter(censoring==0.7, log_ratio==0)
censoring_list <- list(censoring = censoring_params$censoring, 
                       log_ratio = censoring_params$log_ratio, 
                       C1 = censoring_params$C1, 
                       C2 = censoring_params$C2)
res_all <- purrr::map(c("identity","log_positive","boxcox_positive"), function(transform_fn){
    
  # transform_fn <- "identity"
  res <- purrr::map(1:100,function(x){
    data_sim <- simulate_data(n = 50, formula = formula_cens_glmm, type = "glmer",
                              error_variance = 0,
                              b = c(-2,1,1),
                              weibull_params = list(X = list(shape = 0.5, scale = 0.25),
                                                    C = list(shape = 1, scale = censoring_list$C1)),
                              variance_fixeff = 0.5,
                              variance_raneff = 0.5,
                              n_levels_raneff = NULL,
                              n_levels_fixeff = 2,
                              censoring_dependent_on_covariate = ifelse(censoring_list$log_ratio==0,FALSE,TRUE),
                              weibull_params_covariate_dependent_censoring = list(shape = 0.5, scale = censoring_list$C2),
                              transform_fn = transform_fn)
    
    data_sim$I_ran <- sample(data_sim$I)
    cc <- complete_case(data = data_sim,censored_variable = "X",censoring_indicator = "I",
                        formula = y~X+z+(1|r),regression_type = "glmer",weights = "size_tot",family = "binomial")
    cc_su <- multcomp::glht(cc[[4]], matrix(c(0,1,0),nrow=1))
    cc_ran <- complete_case(data = data_sim,censored_variable = "TrVal",censoring_indicator = "I_ran",
                            formula = y~TrVal+z+(1|r),regression_type = "glmer",weights = "size_tot",family = "binomial")
    cc_ran_su <- multcomp::glht(cc_ran[[4]], matrix(c(0,1,0),nrow=1))
    return(tibble::tibble(pvalue = c(summary(cc_su)$test$pval,summary(cc_ran_su)$test$pval),
                   estimate = c(summary(cc_su)$test$coefficients,summary(cc_ran_su)$test$coefficients),
                   nr_censored = rep(dim(data_sim)[1]-sum(data_sim$I),2),
                   type = c("normal","random"),
                   transform_fn = transform_fn))
  })
  res <- dplyr::bind_rows(res)
  res_long <- pivot_longer(data=res,cols=c("pvalue","estimate")) %>% 
    mutate(nr_censored_fac = as.factor(nr_censored),
           nr_censored_rel = nr_censored/50,
           nr_censored_rel_fac = as.factor(nr_censored_rel))
})
res_all_comb <- dplyr::bind_rows(res_all) %>% 
  mutate(transform_fn=factor(transform_fn,levels = c("identity","log_positive","boxcox_positive")))

cc_plot <- function(res_long,facetting_var,y_breaks=c(seq(-10,0,by=2),0,0.05,0.1,seq(0.1,1,by=0.2),1,seq(1,20,by=2))){
  cens_levels <- as.numeric(levels(factor(res_long$nr_censored_rel_fac)))
  cens_levels_red <- cens_levels[c(rep(c(TRUE,FALSE),length(cens_levels)/2))]
  facetting_var <- enquo(facetting_var)
  ggplot(res_long) + 
  geom_jitter(aes(type,value,size=nr_censored_rel, color=nr_censored_rel),width = 0.3) + 
  scale_color_continuous(breaks=cens_levels_red,
                         low="lightgrey", high="black"
                         ) +
  scale_size_continuous(breaks=cens_levels_red,
                        limits = c(min(cens_levels_red),max(cens_levels_red)),
                        range=c(0,2)) +
  scale_y_continuous(breaks = y_breaks) +
  facet_wrap(vars(!!facetting_var,transform_fn)) +
  labs(x="",y="",color = "censoring rate", size="censoring rate") +
  guides(color= guide_legend(), size=guide_legend()) +
  theme(panel.grid.major.x = element_line(size=0.2),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(size=0.3))
}
plt_est <- cc_plot(res_all_comb %>% filter(name=="estimate"),name,y_breaks = c(seq(-30,30,by=2)))  
plt_pval <- cc_plot(res_all_comb %>% filter(name=="pvalue"),name)#,y_breaks = c(seq(-30,30,by=2)))  
plt_comb <- ggpubr::ggarrange(plt_est,plt_pval,nrow = 2, legend = "top",common.legend = TRUE)  

ggsave(paste0(plot_dir,"cc_fits_normal_and_random.png"),plt_comb, width = 12,height=12)
