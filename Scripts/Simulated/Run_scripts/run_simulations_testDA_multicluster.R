### run on cluster with many cores
library(censcyt)
out_ls <- run_simulations_wrapper(reps = 10,
                                  nr_cores = 60,
                                  n_obs_c = c(10,20,50),
                                  simDataTypes = c("glmer"),
                                  formulas = list(glmer = list(formula(Y~Surv(X,I) + Z + (1|R)))),
                                  n_levels_fixeff = 2,
                                  n_levels_raneff = NULL,
                                  betas = list(b0=c(-4),b1=c(-1),b2=c(0.5)),
                                  censoring_rates = c(0.3,0.5,0.7),
                                  imputation_types = c("mrl","cc","pmm","rs","km","ppd"),
                                  run_with_covariate_dependent_censoring = TRUE,
                                  mi_reps = c(50),
                                  error_variance = c(0),
                                  variance_fixeff = c(0.5),
                                  variance_raneff = c(0.5),
                                  number_of_clusters = 20, 
                                  number_of_differential_clusters=9,
                                  multiple_betas_in_list_for_diff_clusters = list(
                                    b1=c(-4,-0.5,0.4),b2=c(-4,-0.6,0.4),b3=c(-4,-0.7,0.4),
                                    b4=c(-4,-0.5,0.5),b5=c(-4,-0.6,0.5),b6=c(-4,-0.7,0.5),
                                    b7=c(-4,-0.5,0.6),b8=c(-4,-0.6,0.6),b9=c(-4,-0.7,0.6)
                                  ),
                                  verbose = TRUE, seed = 123,
                                  transform_fn=c("identity","log_positive","boxcox_positive"))
# save to temporary directory
dir_save <- paste0(head(strsplit(tempdir(), "/")[[1]], -1), collapse = "/")
saveRDS(out_ls, paste0(dir_save,"/test_results_testDA_multicluster_10reps.rds"))
