# set working directory to directory of this script
# setwd(directory_name)
setwd("/home/reto/polybox/ETH/Master_Thesis/Code/Framework/Simulation/Scripts")


################################################################################
## load packages and data
library(tidyverse)
library(ggpubr)
devtools::load_all("/home/retger/FlowCap/scripts/censcyt")

# number of datapoints to simulate, high number to get smoother curve
nn <- 10000000


################################################################################
## plots

set.seed(123)
transform_fns <- c("identity","log_positive","boxcox_positive")
formula_cens_glmm <- formula(y~Surv(X,I)+z+(1|r))
censoring_params <- censcyt::scales_for_censoring(c(0.3,0.5,0.7),log_ratio_val = c(0,0.4))
censoring_params <- censoring_params %>% filter(censoring==0.7)
censoring_list <- list(censoring = rep(censoring_params$censoring,length(transform_fns)), 
                       log_ratio = rep(censoring_params$log_ratio,length(transform_fns)), 
                       C1 = rep(censoring_params$C1,length(transform_fns)), 
                       C2 = rep(censoring_params$C2,length(transform_fns)),
                       transform_fn = rep(transform_fns,each=length(censoring_params$censoring)))
df_ls <- purrr::pmap(censoring_list, function(censoring, log_ratio, C1, C2, transform_fn){
  data_sim <- simulate_data(n = nn, formula = formula_cens_glmm, type = "glmer",
                            error_variance = 0,
                            b = c(-2,1,1),
                            weibull_params = list(X = list(shape = 0.5, scale = 0.25),
                                                  C = list(shape = 1, scale = C1)),
                            variance_fixeff = 0.5,
                            variance_raneff = 0.5,
                            n_levels_raneff = NULL,
                            n_levels_fixeff = 2,
                            censoring_dependent_on_covariate = TRUE,
                            weibull_params_covariate_dependent_censoring = list(shape = 0.5, scale = C2),
                            transform_fn = transform_fn)
  minmax <- list(minx=ifelse(min(data_sim$X)==0,-0.0001,min(data_sim$X)*0.9999),maxx=max(data_sim$X)*1.0001)
  minmax$diff <- minmax$maxx-minmax$minx
  seq_labs <-  seq(minmax$minx,minmax$maxx,by=minmax$diff/100)
  data_sim$X_br <- cut(data_sim$X,breaks = seq_labs,labels=seq_labs[-101])
  censrate <- unlist(purrr::map(levels(data_sim$X_br), function(x){
    I <- data_sim$I[data_sim$X_br==x]
    return(1-sum(I)/length(I))
  }))
  nr_x <- unlist(purrr::map(levels(data_sim$X_br), function(x){
    return(sum(data_sim$X_br==x)/length(data_sim$I))
  }))
  levels_begin <- as.numeric(levels(data_sim$X_br))
  return(data.frame(censrate = censrate, 
                    nr_x = nr_x, 
                    levels_begin=levels_begin, 
                    log_ratio=log_ratio,
                    censoring=censoring, 
                    transform_fn=transform_fn,
                    mean_censoring = 1-sum(data_sim$I)/dim(data_sim)[1]))
  
})
df_ls <- purrr::map(df_ls,function(x){
  x <- arrange(x,levels_begin)
  mutate(x,
         cut_length = rep(diff(x$levels_begin[1:2]),100),
         density=x$nr_x/cut_length)
})
dfplt <- bind_rows(df_ls) %>%
  mutate(cov_dep_cens = factor(log_ratio,levels=c(0,0.4),labels=c("MCAR","MAR")),
         transform_fn = factor(transform_fn,levels = transform_fns,labels = transform_fns))
dir_save <- paste0(head(strsplit(tempdir(), "/")[[1]], -1), collapse = "/")
saveRDS(dfplt,paste0(dir_save,"/censoring_mechanism_df.rds"))