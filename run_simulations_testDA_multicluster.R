# library(censcyt)
devtools::load_all("/home/retger/FlowCap/scripts/censcyt")
out_ls <- run_simulations_wrapper(reps = 1,
                                  nr_cores = 36,
                                  n_obs_c = c(20,50,200),
                                  simDataTypes = c("glmer"),
                                  formulas = list(glmer = list(formula(Y~Surv(X,I) + Z + (1|R)))),
                                  n_levels_fixeff = 2,
                                  n_levels_raneff = NULL,
                                  betas = list(b0=c(-5),b1=c(-1),b2=c(0.2)),
                                  # betas = list(b0=c(0),b1=c(0.5),b2=c(0)),
                                  censoring_rates = c(0.3,0.5,0.7),
                                  # censoring_rates = c(0.3),
                                  imputation_types = c("mrl","cc","pmm","rs","km","ppd"),
                                  # imputation_types = c("ppd"),
                                  run_with_covariate_dependent_censoring = TRUE,
                                  mi_reps = c(50),
                                  error_variance = c(0),
                                  variance_fixeff = c(0.5),
                                  variance_raneff = c(0.5),
                                  number_of_clusters = 20, 
                                  number_of_differential_clusters=9,
                                  multiple_betas_in_list_for_diff_clusters = list(
                                    b0=c(-5,-0.5,0.5),b2=c(-5,-1,0.5),b3=c(-5,-1.5,0.5),#b5=c(-5,-2,0.2),
                                    b0=c(-5,-0.5,0.5),b2=c(-5,-1,0.5),b3=c(-5,-1.5,0.5),#b5=c(-5,-2,0.2),
                                    b0=c(-5,-0.5,0.5),b2=c(-5,-1,0.5),b3=c(-5,-1.5,0.5)#b5=c(-5,-2,0.2)),
                                  ),
                                  # size_of_random_subset_for_testing = 1,
                                  verbose = TRUE, seed = 123,
                                  transform_fn=c("identity","log_positive"))
# size_of_random_subset_for_testing=100)
dir_save <- paste0(head(strsplit(tempdir(), "/")[[1]], -1), collapse = "/")
saveRDS(out_ls, paste0(dir_save,"/test_results_testDA_multicluster.rds"))

# devtools::load_all("/home/reto/polybox/ETH/Master_Thesis/Code/censcyt")
