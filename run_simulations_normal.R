# library(censcyt)
devtools::load_all("/home/retger/FlowCap/scripts/censcyt")
out_ls <- run_simulations_wrapper(reps = 10,
                             nr_cores = 60,
                             n_obs_c = c(20),
                             simDataTypes = c("glmer"),
                             formulas = list(glmer = list(formula(Y~Surv(X,I) + Z + (1|R)))),
                             n_levels_fixeff = 2,
                             n_levels_raneff = NULL,
                             betas = list(b0=c(-5),b1=c(-1),b2=c(0.5)),
                             # betas = list(b0=c(0),b1=c(0.5),b2=c(0)),
                             censoring_rates = c(0.3,0.5,0.7),
                             # censoring_rates = c(0.3),
                             imputation_types = c("mrl","km","rs","cc","pmm","ppd"),
                             # imputation_types = c("ppd"),
                             run_with_covariate_dependent_censoring = TRUE,
                             mi_reps = c(50),
                             error_variance = c(0),
                             variance_fixeff = c(0.5),
                             variance_raneff = c(0.5),
                             verbose = FALSE, seed = 123,
                             transform_fn=c("identity","boxcox","sqrt","log","log_positive","boxcox_positive"))
                             # size_of_random_subset_for_testing = 120)
dir_save <- paste0(head(strsplit(tempdir(), "/")[[1]], -1), collapse = "/")
saveRDS(out_ls, paste0(dir_save,"/test_results_normal_glmer_only.rds"))


# library(censcyt)
# devtools::load_all("/home/reto/polybox/ETH/Master_Thesis/Code/censcyt")
# # simulations_data_processor(out_ls[[2]])
# out_ls <- run_simulations_wrapper(reps = 2,
#                                   nr_cores = 1,
#                                   n_obs_c = c(20),
#                                   simDataTypes = c("glmer"),
#                                   formulas = list(glmer = list(formula(Y~Surv(X,I) + Z + (1|R)))),
#                                   n_levels_fixeff = 2,
#                                   n_levels_raneff = NULL,
#                                   betas = list(b0=c(-2),b1=c(1),b2=c(1)),
#                                   # betas = list(b0=c(0),b1=c(0.5),b2=c(0)),
#                                   censoring_rates = c(0.7),
#                                   # censoring_rates = c(0.3),
#                                   imputation_types = c("cc","rs","km"),
#                                   run_with_covariate_dependent_censoring = TRUE,
#                                   mi_reps = c(10),
#                                   error_variance = c(0.1),
#                                   variance_fixeff = c(0.1),
#                                   variance_raneff = c(0.1),
#                                   verbose = TRUE, seed = 1234)
# number_of_clusters = 5)#,
                                  # size_of_random_subset_for_testing = 10)

# devtools::load_all("/home/reto/polybox/ETH/Master_Thesis/Code/censcyt")
# out_ls <- run_simulations_wrapper(reps = 1,
#                                   n_obs_c = c(50),
#                                   simDataTypes = c("glmer"),
#                                   formulas = list(glmer = list(formula(Y~Surv(X,I) + Z + (1|R)))),
#                                   n_levels_fixeff = 2,
#                                   n_levels_raneff = NULL,
#                                   betas = list(b0=c(1),b1=c(1),b2=c(1)),
#                                   censoring_rates = c(0.3),
#                                   # censoring_rates = c(0.3),
#                                   imputation_types = c("rs"),
#                                   mi_reps = c(5),
#                                   error_variance = c(0.1),
#                                   number_of_cluster = 4,
#                                   verbose = TRUE, seed = 123) 

# around 5.5 min for cens 0.5, n_obs 500, glmer, mi_reps 200

