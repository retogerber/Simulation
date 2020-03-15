### run on cluster with many cores
devtools::load_all("/home/retger/FlowCap/scripts/censcyt")
out_ls <- suppressMessages(run_simulations_wrapper(reps = 1000,
                                                   nr_cores = 40,
                                                   n_obs_c = c(10,20,50),
                                                   simDataTypes = c("glmer"),
                                                   formulas = list(glmer = list(formula(Y~Surv(X,I) + Z + (1|R)))
                                                   ),
                                                   n_levels_fixeff = 2,
                                                   n_levels_raneff = NULL,
                                                   betas = list(b0=c(-5),b1=c(0),b2=c(0.5)),
                                                   censoring_rates = c(0.7),
                                                   imputation_types = c("mrl","km","rs","cc","pmm","ppd"),
                                                   run_with_covariate_dependent_censoring = FALSE,
                                                   mi_reps = c(50),
                                                   error_variance = c(0),
                                                   variance_fixeff = c(0.5),
                                                   variance_raneff = c(0.1),
                                                   verbose = TRUE, seed = 123,
                                                   transform_fn=c("log_positive"),
))
# save to temporary directory
dir_save <- paste0(head(strsplit(tempdir(), "/")[[1]], -1), collapse = "/")
saveRDS(out_ls, paste0(dir_save,"/test_results_normal_null.rds"))
