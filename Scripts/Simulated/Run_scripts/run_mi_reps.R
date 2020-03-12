### run on cluster with 6 cores
devtools::load_all("/home/retger/FlowCap/scripts/censcyt")
library(tidyverse)

set.seed(123)
formula_cens_glmm <- formula(y~Surv(X,I)+z+(1|r))
censoring_params <- censcyt::scales_for_censoring(c(0.3,0.5,0.7),log_ratio_val = c(0,0.4))
censoring_list <- list(censoring = censoring_params$censoring, 
                       log_ratio = censoring_params$log_ratio, 
                       C1 = censoring_params$C1, 
                       C2 = censoring_params$C2)
mi_reps <- c(2,3,4,5,7,10,20,30,40,50,75,100,125,150,200,300,400,500)
methods_ests <- c("rs","km","mrl","pmm","ppd")
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
    })
    return(bind_rows(outs))
  })
  return(bind_rows(outs))
},mc.cores=6)

# save to temporary directory
dir_save <- paste0(head(strsplit(tempdir(), "/")[[1]], -1), collapse = "/")
saveRDS(res,paste0(dir_save,"/res_plot_mi_reps.rds"))
